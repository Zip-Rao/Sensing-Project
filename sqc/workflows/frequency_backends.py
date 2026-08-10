"""sqc.workflows.frequency_backends — SQC measurement executor for the
frequency calibration state machine.

Translates abstract :class:`Command` objects into real measurements on a
:class:`TransmonQubit` via the existing :mod:`sqc.calibration.frequency`
infrastructure.

Provides:
- :class:`SQCExecutor` — the main synchronous executor for simulation.
- :class:`FaultInjectionExecutor` — wraps an executor with configurable
  fault injection for testing recovery paths.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import time as _time

import numpy as np

from sqc.calibration.frequency import FrequencyMeasurement
from sqc.calibration.frequency_control import (
    DampedSecantTracker,
    FrequencyEstimate,
)
from sqc.workflows.frequency_state_machine import (
    AcquireFrequency,
    BudgetExhausted,
    CancelRequested,
    Command,
    ControlApplied,
    Event,
    FrequencyCalibrationConfig,
    InterlockTriggered,
    MeasurementRejected,
    MeasurementSucceeded,
    MeasurementTechnicalFailure,
    MonitorFrequency,
    ReasonCode,
    SafeHold,
    SafeHoldApplied,
    TimerElapsed,
    TrackFrequency,
    VerifyFrequency,
)


# ===================================================================
# SQC Executor
# ===================================================================


@dataclass
class SQCExecutor:
    """Synchronous executor that runs commands on a real TransmonQubit.

    Parameters
    ----------
    qubit : TransmonQubit
        The qubit to measure.
    config : FrequencyCalibrationConfig
        Protocol thresholds and budget parameters.
    f_target : float
        Target qubit frequency (angular, rad·GHz).
    """

    qubit: object
    config: FrequencyCalibrationConfig = field(default_factory=FrequencyCalibrationConfig)
    f_target: float = 0.0

    # ---- internal measurement objects (created lazily) --------------------
    _ramsey_meas: FrequencyMeasurement | None = field(default=None, repr=False)
    _transient_meas: FrequencyMeasurement | None = field(default=None, repr=False)
    _verify_meas: FrequencyMeasurement | None = field(default=None, repr=False)

    def execute(self, cmd: Command) -> Event:
        """Execute a single command and return the corresponding event.

        Parameters
        ----------
        cmd : Command
            The command to execute (from :meth:`FrequencyStateMachine.next_command`).

        Returns
        -------
        Event
            A :class:`MeasurementSucceeded`, :class:`MeasurementTechnicalFailure`,
            :class:`SafeHoldApplied`, etc.
        """
        match cmd:
            case AcquireFrequency():
                return self._execute_acquire(cmd)
            case TrackFrequency():
                return self._execute_track(cmd)
            case VerifyFrequency():
                return self._execute_verify(cmd)
            case MonitorFrequency():
                return self._execute_monitor(cmd)
            case SafeHold():
                return self._execute_safe_hold(cmd)
            case _:
                return MeasurementTechnicalFailure(
                    command_id=cmd.command_id,
                    reason=f"unknown command type: {type(cmd).__name__}",
                )

    # ------------------------------------------------------------------
    # Per-command executors
    # ------------------------------------------------------------------

    def _execute_acquire(self, cmd: AcquireFrequency) -> Event:
        """Wide-range Ramsey acquisition (absolute frequency measurement)."""
        t0 = _time.time()
        try:
            meas = self._get_ramsey()
            # Acquire at zero flux (sweet spot) — wide sweep for absolute f01
            f = meas.measure(flux=0.0)
            elapsed = _time.time() - t0
            n_shots = len(meas.tau_list) * (1 if meas.f_artificial is not None else 2)
            return MeasurementSucceeded(
                command_id=cmd.command_id,
                frequency=f,
                uncertainty=0.0,
                valid=True,
                ambiguous=False,
                method="ramsey",
                shots=n_shots,
                elapsed_time=elapsed,
                applied_bias=0.0,
                diagnostics={
                    "solver_calls": n_shots,
                    "role": cmd.role,
                },
            )
        except Exception as exc:
            return MeasurementTechnicalFailure(
                command_id=cmd.command_id,
                reason=str(exc),
                elapsed_time=_time.time() - t0,
            )

    def _execute_track(self, cmd: TrackFrequency) -> Event:
        """Transient single-point measurement at the proposed bias/drive."""
        t0 = _time.time()
        try:
            meas = self._get_transient()
            f = meas.measure(
                flux=cmd.candidate_bias,
                omega_d=cmd.predicted_drive,
            )
            elapsed = _time.time() - t0
            return MeasurementSucceeded(
                command_id=cmd.command_id,
                frequency=f,
                uncertainty=0.0,
                valid=True,
                ambiguous=False,
                method="transient",
                shots=2,  # 2 orthogonal readout solves
                elapsed_time=elapsed,
                applied_bias=cmd.candidate_bias,
                applied_drive=cmd.predicted_drive,
                diagnostics={
                    "solver_calls": 2,
                    "s_hat": cmd.tracker_spec.get("s_hat"),
                },
            )
        except Exception as exc:
            return MeasurementTechnicalFailure(
                command_id=cmd.command_id,
                reason=str(exc),
                elapsed_time=_time.time() - t0,
            )

    def _execute_verify(self, cmd: VerifyFrequency) -> Event:
        """Independent Ramsey verification at frozen bias."""
        t0 = _time.time()
        try:
            meas = self._get_verify()
            # Use double-sweep Ramsey (f_artificial=None) for unbiased verification
            f = meas.measure(flux=cmd.frozen_bias, omega_d=cmd.drive)
            elapsed = _time.time() - t0
            n_shots = 2 * len(meas.tau_list)  # double-sweep
            return MeasurementSucceeded(
                command_id=cmd.command_id,
                frequency=f,
                uncertainty=0.0,
                valid=True,
                ambiguous=False,
                method="ramsey",
                shots=n_shots,
                elapsed_time=elapsed,
                applied_bias=cmd.frozen_bias,
                applied_drive=cmd.drive,
                diagnostics={
                    "solver_calls": n_shots,
                    "verifier_spec": cmd.verifier_spec,
                },
            )
        except Exception as exc:
            return MeasurementTechnicalFailure(
                command_id=cmd.command_id,
                reason=str(exc),
                elapsed_time=_time.time() - t0,
            )

    def _execute_monitor(self, cmd: MonitorFrequency) -> Event:
        """Low-cost Lock monitor.

        First version: if no verified low-cost monitor is available, falls
        back to a single transient measurement at the locked bias.  The
        monitor is *not* a full Ramsey — it's a quick check with the
        transient discriminator.
        """
        t0 = _time.time()
        try:
            meas = self._get_transient()
            f = meas.measure(flux=cmd.locked_bias, omega_d=self.f_target)
            elapsed = _time.time() - t0
            return MeasurementSucceeded(
                command_id=cmd.command_id,
                frequency=f,
                uncertainty=0.0,
                valid=True,
                ambiguous=False,
                method="transient",
                shots=2,
                elapsed_time=elapsed,
                applied_bias=cmd.locked_bias,
                applied_drive=self.f_target,
                diagnostics={
                    "solver_calls": 2,
                    "monitor_spec": cmd.monitor_spec,
                },
            )
        except Exception as exc:
            return MeasurementTechnicalFailure(
                command_id=cmd.command_id,
                reason=str(exc),
                elapsed_time=_time.time() - t0,
            )

    def _execute_safe_hold(self, cmd: SafeHold) -> Event:
        """Apply safe hold bias."""
        # In simulation, this is a no-op (we just record the event).
        # In hardware, this would set the bias DAC to the safe value.
        return SafeHoldApplied(
            command_id=cmd.command_id,
            applied_bias=cmd.bias,
        )

    # ------------------------------------------------------------------
    # Lazy measurement objects
    # ------------------------------------------------------------------

    def _get_ramsey(self) -> FrequencyMeasurement:
        if self._ramsey_meas is None:
            self._ramsey_meas = FrequencyMeasurement(
                qubit=self.qubit,
                method="ramsey",
                f_artificial=0.1,  # single-sweep for acquisition
            )
        return self._ramsey_meas

    def _get_transient(self) -> FrequencyMeasurement:
        if self._transient_meas is None:
            self._transient_meas = FrequencyMeasurement(
                qubit=self.qubit,
                method="transient",
                order=1,  # linear — fast; set order=3 + g3_source for accuracy
            )
        return self._transient_meas

    def _get_verify(self) -> FrequencyMeasurement:
        if self._verify_meas is None:
            self._verify_meas = FrequencyMeasurement(
                qubit=self.qubit,
                method="ramsey",
                f_artificial=None,  # double-sweep — unbiased, signed
            )
        return self._verify_meas


# ===================================================================
# Fault-injection executor (for testing recovery paths)
# ===================================================================


class FaultInjectionExecutor:
    """Wraps an :class:`SQCExecutor` with configurable fault injection.

    Parameters
    ----------
    inner : SQCExecutor
        The real executor to delegate to.
    fail_on_state : FrequencyState or None
        Inject a technical failure when the machine is in this state.
    reject_on_state : FrequencyState or None
        Inject a measurement rejection (e.g. low confidence).
    ambiguous_on_state : FrequencyState or None
        Return an ambiguous measurement result.
    fail_count : int
        Number of times to inject the fault before falling through.
    """

    def __init__(
        self,
        inner: SQCExecutor,
        fail_on_state=None,
        reject_on_state=None,
        ambiguous_on_state=None,
        fail_count: int = 1,
    ):
        self._inner = inner
        self._fail_on_state = fail_on_state
        self._reject_on_state = reject_on_state
        self._ambiguous_on_state = ambiguous_on_state
        self._fail_count = fail_count
        self._fail_remaining = fail_count
        self._machine_state = None

    def set_machine_state(self, state):
        """Inform the executor of the current machine state (called by runtime)."""
        self._machine_state = state

    def execute(self, cmd: Command) -> Event:
        """Execute with fault injection."""
        # Technical failure injection
        if self._fail_on_state is not None and self._machine_state == self._fail_on_state:
            if self._fail_remaining > 0:
                self._fail_remaining -= 1
                return MeasurementTechnicalFailure(
                    command_id=cmd.command_id,
                    reason="injected fault",
                )

        # Rejection injection
        if self._reject_on_state is not None and self._machine_state == self._reject_on_state:
            if self._fail_remaining > 0:
                self._fail_remaining -= 1
                return MeasurementRejected(
                    command_id=cmd.command_id,
                    reason_code=ReasonCode.CONFIDENCE_INSUFFICIENT,
                )

        # Ambiguity injection
        if self._ambiguous_on_state is not None and self._machine_state == self._ambiguous_on_state:
            if self._fail_remaining > 0:
                self._fail_remaining -= 1
                # Run the real measurement but mark it ambiguous
                event = self._inner.execute(cmd)
                if isinstance(event, MeasurementSucceeded):
                    event.ambiguous = True
                return event

        return self._inner.execute(cmd)
