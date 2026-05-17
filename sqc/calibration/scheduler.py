"""sqc.calibration.scheduler — CalibrationScheduler (control room).

Lightweight dispatcher that maintains a registry of calibration tasks with
dependency tracking.  Manual control today; DAG automation (Kelly 2018
``check_state → maintain → calibrate``) is reserved as interface stubs
for future implementation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from sqc.calibration.base import Calibration, CalibrationTable


@dataclass
class CalibrationScheduler:
    """Qubit calibration control room.

    Maintains a registry of named calibration tasks with optional
    dependencies.  Call :meth:`run` to execute a task by name, or
    :meth:`run_next` to execute the next ready task.

    Pre-registered tasks (auto-populated on init)::

        flux_response_ramsey       → FluxResponseCalibration(method="ramsey")
        frequency_measurement      → FrequencyMeasurement (method="ramsey" by default)
        frequency_closed_loop      → SinglePointFrequencyCalibration(method="closed_loop")
        waveform_transfer_function → WaveformCalibration(method="transfer_function")
        waveform_predistortion     → WaveformCalibration(method="predistortion")

    Dependencies
    ------------
    - ``frequency_closed_loop`` depends on ``flux_response_ramsey``
      (provides bracketing [V_a, V_b]).
    - ``waveform_predistortion`` depends on ``waveform_transfer_function``
      (provides the measured transfer function).

    DAG interfaces (stubs for future Kelly 2018 automation)
    --------------------------------------------------------
    - :meth:`check_state` — inspect qubit state, return diagnostics dict.
    - :meth:`maintain` — decide whether re-calibration is needed.
    - :meth:`auto_calibrate` — run the full maintain→calibrate loop.
    """

    # registry: task_name → Calibration *class* (not instance)
    tasks: dict[str, type[Calibration]] = field(default_factory=dict)
    # dependency graph: task_name → set of prerequisite task names
    dependencies: dict[str, set[str]] = field(default_factory=dict)
    # completed results: task_name → CalibrationTable
    completed: dict[str, CalibrationTable] = field(default_factory=dict)
    # currently running task
    current_task: str | None = None

    # ------------------------------------------------------------------
    # Registration
    # ------------------------------------------------------------------

    def register(
        self,
        name: str,
        cal_class: type[Calibration],
        depends_on: list[str] | None = None,
    ) -> None:
        """Register a calibration task.

        Parameters
        ----------
        name : str
            Unique task name.
        cal_class : type[Calibration]
            Calibration subclass (uninstantiated).
        depends_on : list[str] or None
            Names of prerequisite tasks.
        """
        self.tasks[name] = cal_class
        if depends_on:
            self.dependencies[name] = set(depends_on)

    def register_defaults(self) -> None:
        """Populate the scheduler with the standard calibration registry."""
        from sqc.calibration.frequency import (
            FluxResponseCalibration,
            FrequencyMeasurement,
            SinglePointFrequencyCalibration,
        )
        from sqc.calibration.waveform import WaveformCalibration

        self.register("flux_response_ramsey", FluxResponseCalibration)
        self.register("frequency_measurement", FrequencyMeasurement)
        self.register(
            "frequency_closed_loop",
            SinglePointFrequencyCalibration,
            depends_on=["flux_response_ramsey"],
        )
        self.register("waveform_transfer_function", WaveformCalibration)
        self.register(
            "waveform_predistortion",
            WaveformCalibration,
            depends_on=["waveform_transfer_function"],
        )

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def run(self, name: str, **override_kwargs: Any) -> CalibrationTable:
        """Run a calibration task by name.

        Parameters
        ----------
        name : str
            Registered task name.
        **override_kwargs
            Passed to the Calibration constructor, overriding defaults.

        Returns
        -------
        CalibrationTable

        Raises
        ------
        KeyError
            If *name* is not registered.
        RuntimeError
            If dependencies are not satisfied.
        """
        if name not in self.tasks:
            raise KeyError(
                f"Unknown task '{name}'. Registered: {list(self.tasks)}"
            )

        # Check dependencies
        deps = self.dependencies.get(name, set())
        missing = [d for d in deps if d not in self.completed]
        if missing:
            raise RuntimeError(
                f"Task '{name}' depends on {missing}, which have not been "
                f"completed. Run them first or call run_next()."
            )

        cal_class = self.tasks[name]
        self.current_task = name

        # Build instance — merge dependency results into kwargs when relevant
        cal_kwargs: dict[str, Any] = dict(override_kwargs)
        cal_kwargs = self._inject_dependencies(name, cal_kwargs)

        instance = cal_class(**cal_kwargs)
        result = instance.calibrate()
        self.completed[name] = result
        self.current_task = None
        return result

    def run_next(self) -> CalibrationTable | None:
        """Run the first ready task that has not yet been completed.

        A task is *ready* when all its dependencies are satisfied.

        Returns
        -------
        CalibrationTable or None
            Result of the executed task, or None if no tasks are ready.
        """
        for name in self.tasks:
            if name in self.completed:
                continue
            deps = self.dependencies.get(name, set())
            if deps.issubset(self.completed.keys()):
                return self.run(name)
        return None

    def _inject_dependencies(
        self, name: str, kw: dict[str, Any],
    ) -> dict[str, Any]:
        """Inject results from completed dependencies into kwargs.

        Avoids the user having to manually thread calibration tables.
        """
        deps = self.dependencies.get(name, set())

        if "flux_response_ramsey" in deps:
            table = self.completed["flux_response_ramsey"]
            if "V_a" not in kw and "V_b" not in kw and table.inputs is not None:
                # Use the flux scan range as voltage bracket
                kw.setdefault("V_a", float(np.min(table.inputs)))
                kw.setdefault("V_b", float(np.max(table.inputs)))
            if "f_target" not in kw:
                # Default: use the qubit bare frequency (caller should
                # override with a specific target)
                pass

        if "waveform_transfer_function" in deps:
            tf_table = self.completed["waveform_transfer_function"]
            if "transfer_model" not in kw:
                fit = tf_table.fit_params
                if "amplitudes" in fit:
                    from sqc.hardware.distortion import MultiExponentialDistortion
                    kw["transfer_model"] = MultiExponentialDistortion(
                        amplitudes=fit["amplitudes"],
                        taus=fit["taus"],
                    )
                elif "amplitude" in fit:
                    from sqc.hardware.distortion import SingleExponentialDistortion
                    kw["transfer_model"] = SingleExponentialDistortion(
                        amplitude=fit["amplitude"],
                        tau=fit["tau"],
                    )

        return kw

    # ------------------------------------------------------------------
    # Query
    # ------------------------------------------------------------------

    def get_result(self, name: str) -> CalibrationTable:
        """Retrieve a completed calibration result.

        Raises
        ------
        KeyError
            If *name* has not been completed.
        """
        if name not in self.completed:
            raise KeyError(
                f"Task '{name}' has not been completed. "
                f"Completed: {list(self.completed)}"
            )
        return self.completed[name]

    def status(self) -> dict:
        """Return current scheduler status.

        Returns
        -------
        dict
            Keys: "current_task", "completed", "pending", "ready".
        """
        ready = [
            n for n in self.tasks
            if n not in self.completed
            and self.dependencies.get(n, set()).issubset(self.completed.keys())
        ]
        pending = [
            n for n in self.tasks
            if n not in self.completed
            and n not in ready
        ]
        return {
            "current_task": self.current_task,
            "completed": list(self.completed),
            "ready": ready,
            "pending": pending,
        }

    # ------------------------------------------------------------------
    # DAG stubs (Kelly 2018 — future automation)
    # ------------------------------------------------------------------

    def check_state(self) -> dict:
        """Inspect qubit state, return diagnostics.

        Placeholder for Kelly 2018 DAG ``check_state`` node.  In the
        future this will read out qubit parameters (frequency, T1, T2,
        etc.) and return a diagnostics dict that :meth:`maintain` uses
        to decide whether re-calibration is needed.

        Returns
        -------
        dict
            Diagnostics snapshot (empty for now).
        """
        return {}

    def maintain(self) -> str | None:
        """Decide which calibration task to run next.

        Placeholder for Kelly 2018 DAG ``maintain`` node.  Consumes the
        dict from :meth:`check_state` and returns the name of the
        highest-priority task that needs re-running, or None if
        everything is within spec.

        Returns
        -------
        str or None
            Task name to run, or None if no action needed.
        """
        return None

    def auto_calibrate(self, target_state: dict | None = None) -> list[CalibrationTable]:
        """Run the full maintain → calibrate loop.

        Placeholder.  Future implementation will:
        1. Call :meth:`check_state`
        2. Call :meth:`maintain` to get the next task
        3. Call :meth:`run` on that task
        4. Repeat until maintain returns None

        Parameters
        ----------
        target_state : dict or None
            Desired qubit state (frequency, etc.).  Reserved.

        Returns
        -------
        list[CalibrationTable]
            Results of all calibration tasks that were run.
        """
        results: list[CalibrationTable] = []
        while True:
            task_name = self.maintain()
            if task_name is None:
                break
            results.append(self.run(task_name))
        return results
