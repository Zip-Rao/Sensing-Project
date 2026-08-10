"""sqc.workflows.frequency_runtime — synchronous runtime for the frequency
calibration state machine.

Orchestrates the command→event loop, manages the :class:`DampedSecantTracker`
for the Track state, tracks budgets, and provides checkpoint / journal
support.

Usage::

    runtime = FrequencyCalibrationRuntime(qubit, f_target=..., config=...)
    result = runtime.run()
    print(result["state"], result["run_status"])
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import json as _json
import time as _time

import numpy as np

from sqc.calibration.frequency_control import (
    DampedSecantTracker,
    FrequencyEstimate,
)
from sqc.workflows.frequency_backends import SQCExecutor
from sqc.workflows.frequency_state_machine import (
    AcquireFrequency,
    Command,
    FrequencyCalibrationConfig,
    FrequencyState,
    FrequencyStateMachine,
    MonitorFrequency,
    ReasonCode,
    RunStatus,
    SafeHold,
    StateMachineSnapshot,
    TrackFrequency,
    VerifyFrequency,
)


# ===================================================================
# Runtime
# ===================================================================


@dataclass
class FrequencyCalibrationRuntime:
    """Synchronous runtime for the frequency calibration state machine.

    Parameters
    ----------
    qubit : TransmonQubit
        The qubit to calibrate.
    f_target : float
        Target qubit frequency (angular, rad·GHz).
    config : FrequencyCalibrationConfig or None
        Protocol thresholds and budgets.
    executor : SQCExecutor or None
        Optional custom executor (default creates an :class:`SQCExecutor`).
    """

    qubit: object
    f_target: float
    config: FrequencyCalibrationConfig | None = None
    executor: SQCExecutor | None = None

    # ---- internal state --------------------------------------------------
    _machine: FrequencyStateMachine | None = field(default=None, repr=False)
    _tracker: DampedSecantTracker | None = field(default=None, repr=False)
    _journal: list[dict] = field(default_factory=list, repr=False)
    _checkpoints: list[StateMachineSnapshot] = field(default_factory=list, repr=False)

    def __post_init__(self):
        if self.config is None:
            self.config = FrequencyCalibrationConfig()
        if self.executor is None:
            self.executor = SQCExecutor(
                qubit=self.qubit,
                config=self.config,
                f_target=self.f_target,
            )

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self) -> dict:
        """Run the calibration event loop to completion.

        Returns
        -------
        dict
            ``state``, ``run_status``, ``f_final``, ``candidate_bias``,
            ``transition_log``, ``journal``, ``n_commands``, ``elapsed``.
        """
        t_start = _time.time()
        self._machine = FrequencyStateMachine(
            config=self.config, f_target=self.f_target,
        )
        self._machine.start()
        self._journal = []

        while self._machine.run_status == RunStatus.RUNNING:
            cmd = self._machine.next_command()

            # Terminal command → exit loop
            if isinstance(cmd, SafeHold):
                self._record_journal(cmd, None)
                break

            # Pre-process: for Track commands, use the tracker to refine the proposal
            cmd = self._enrich_command(cmd)

            # Execute
            event = self.executor.execute(cmd)
            self._record_journal(cmd, event)

            # Handle event → state transition
            self._machine.handle(event)

            # Post-process: update tracker after Track measurements
            self._post_process(cmd, event)

            # Checkpoint periodically
            if len(self._journal) % 10 == 0:
                self._checkpoints.append(self._machine.snapshot())

        elapsed = _time.time() - t_start
        return {
            "state": self._machine.state.value,
            "run_status": self._machine.run_status.value,
            "f_final": self._machine._f_hat,
            "candidate_bias": self._machine._candidate_bias,
            "transition_log": self._machine._transition_log,
            "journal": self._journal,
            "n_commands": len(self._journal),
            "elapsed": elapsed,
            "checkpoints": len(self._checkpoints),
        }

    # ------------------------------------------------------------------
    # Command enrichment (tracker integration)
    # ------------------------------------------------------------------

    def _enrich_command(self, cmd: Command) -> Command:
        """Refine Track commands with the DampedSecantTracker proposal.

        For Track state, the state machine produces a baseline TrackFrequency
        command.  We use the tracker to propose a refined bias/drive.
        """
        if not isinstance(cmd, TrackFrequency):
            return cmd

        # Lazy-init tracker on first Track command
        if self._tracker is None:
            self._tracker = DampedSecantTracker(
                f_target=self.f_target,
                damping=self.config.damping,
                first_bias_step=self.config.first_bias_step,
                max_bias_step=self.config.max_bias_step,
                converge_streak=1,  # tracker checks within-step; SM checks overall
                max_iter=999,        # managed by SM
                epsilon_f=self.config.epsilon_enter,  # candidate threshold
            )

        # If we have a prior snapshot, propose next step
        if self._machine._tracker_snapshot is not None:
            from sqc.calibration.frequency_control import TrackSnapshot
            snap_dict = self._machine._tracker_snapshot
            snap = TrackSnapshot(**snap_dict)
            proposal = self._tracker.propose(snap)
            if not proposal.converged:
                # Predict drive using track policy
                predicted_drive = DampedSecantTracker.resolve_track_drive(
                    V=proposal.V_next,
                    V_prev=snap.V,
                    f_prev=snap.f,
                    s_hat=proposal.s_hat,
                )
                return TrackFrequency(
                    command_id=cmd.command_id,
                    candidate_bias=proposal.V_next,
                    predicted_drive=predicted_drive,
                    tracker_spec={
                        "s_hat": proposal.s_hat,
                        "step": proposal.step,
                    },
                )

        # First Track step or no tracker state: use the SM's default
        return cmd

    # ------------------------------------------------------------------
    # Post-processing (tracker update)
    # ------------------------------------------------------------------

    def _post_process(self, cmd: Command, event):
        """Update tracker state after a Track measurement."""
        if not isinstance(cmd, TrackFrequency):
            return

        from sqc.workflows.frequency_state_machine import MeasurementSucceeded

        if not isinstance(event, MeasurementSucceeded):
            return

        # Initialise tracker state from first Track measurement
        if self._machine._tracker_snapshot is None and self._tracker is not None:
            seed = FrequencyEstimate(
                frequency=event.frequency,
                uncertainty=event.uncertainty,
                valid=event.valid,
                ambiguous=event.ambiguous,
                method=event.method,
                shots=event.shots,
                elapsed_time=event.elapsed_time,
                diagnostics=event.diagnostics,
            )
            snap = self._tracker.initialize(seed, event.applied_bias)
            self._machine._tracker_snapshot = {
                "V": snap.V,
                "V_prev": snap.V_prev,
                "f": snap.f,
                "f_prev": snap.f_prev,
                "e": snap.e,
                "e_prev": snap.e_prev,
                "f_target": snap.f_target,
                "best_V": snap.best_V,
                "best_abs_e": snap.best_abs_e,
                "streak": snap.streak,
                "n_iter": snap.n_iter,
                "s_hat": snap.s_hat,
            }
            return

        # Update existing tracker state
        if self._machine._tracker_snapshot is not None and self._tracker is not None:
            from sqc.calibration.frequency_control import (
                TrackProposal,
                TrackSnapshot,
            )
            snap_dict = self._machine._tracker_snapshot
            old_snap = TrackSnapshot(**snap_dict)

            # Reconstruct the proposal from the command
            proposal = TrackProposal(
                V_next=cmd.candidate_bias,
                step=cmd.tracker_spec.get("step", 0.0),
                s_hat=cmd.tracker_spec.get("s_hat"),
                converged=False,
            )

            estimate = FrequencyEstimate(
                frequency=event.frequency,
                uncertainty=event.uncertainty,
            )
            result = self._tracker.accept(old_snap, estimate, proposal)
            new_snap = result.snapshot
            self._machine._tracker_snapshot = {
                "V": new_snap.V,
                "V_prev": new_snap.V_prev,
                "f": new_snap.f,
                "f_prev": new_snap.f_prev,
                "e": new_snap.e,
                "e_prev": new_snap.e_prev,
                "f_target": new_snap.f_target,
                "best_V": new_snap.best_V,
                "best_abs_e": new_snap.best_abs_e,
                "streak": new_snap.streak,
                "n_iter": new_snap.n_iter,
                "s_hat": new_snap.s_hat,
            }

    # ------------------------------------------------------------------
    # Journal
    # ------------------------------------------------------------------

    def _record_journal(self, cmd: Command, event=None):
        """Record a journal entry for this command→event cycle."""
        entry = {
            "command": type(cmd).__name__,
            "command_id": cmd.command_id,
            "event": type(event).__name__ if event is not None else None,
        }
        if hasattr(cmd, "candidate_bias"):
            entry["bias"] = getattr(cmd, "candidate_bias", None)
        elif hasattr(cmd, "frozen_bias"):
            entry["bias"] = getattr(cmd, "frozen_bias", None)
        elif hasattr(cmd, "locked_bias"):
            entry["bias"] = getattr(cmd, "locked_bias", None)

        if event is not None and hasattr(event, "frequency"):
            entry["frequency"] = event.frequency

        self._journal.append(entry)

    # ------------------------------------------------------------------
    # Checkpoint / save
    # ------------------------------------------------------------------

    def save_journal(self, path: str):
        """Save the journal to a JSONL file."""
        with open(path, "w") as f:
            for entry in self._journal:
                f.write(_json.dumps(entry, default=str) + "\n")

    def latest_checkpoint(self) -> StateMachineSnapshot | None:
        """Return the most recent checkpoint, or None."""
        return self._checkpoints[-1] if self._checkpoints else None
