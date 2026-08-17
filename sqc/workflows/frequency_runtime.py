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

from dataclasses import asdict, dataclass, field, fields
import json as _json
import threading as _threading
import time as _time
import uuid as _uuid

import numpy as np

from sqc.calibration.frequency_control import (
    DampedSecantTracker,
    FrequencyEstimate,
)
from sqc.workflows.frequency_backends import SQCExecutor
from sqc.workflows.frequency_state_machine import (
    Budget,
    CancelRequested,
    Command,
    FrequencyCalibrationConfig,
    FrequencyState,
    FrequencyStateMachine,
    MeasurementRejected,
    MeasurementTechnicalFailure,
    MonitorFrequency,
    MeasurementSucceeded,
    RunStatus,
    SafeHold,
    SafeHoldApplied,
    StateMachineSnapshot,
    TrackFrequency,
)


# ===================================================================
# Cooperative interruption
# ===================================================================


class CancellationToken:
    """Thread-safe, first-request-wins cancellation signal.

    Signal handlers and UI/API threads should only call :meth:`request`.
    The runtime thread converts the request into a ``CancelRequested`` event
    at a command boundary or after the current synchronous measurement returns.
    """

    def __init__(self):
        self._event = _threading.Event()
        self._lock = _threading.Lock()
        self._request: dict | None = None

    def request(self, reason: str = "user_requested", source: str = "api") -> bool:
        """Set the cancellation signal; return ``True`` for the first request."""
        with self._lock:
            if self._event.is_set():
                return False
            self._request = {
                "reason": str(reason),
                "source": str(source),
                "requested_at": _time.time(),
            }
            self._event.set()
            return True

    @property
    def requested(self) -> bool:
        return self._event.is_set()

    def wait(self, timeout: float) -> bool:
        """Wait up to *timeout* seconds; return whether cancellation was set."""
        return self._event.wait(max(0.0, float(timeout)))

    def snapshot(self) -> dict | None:
        """Return a serialisable copy of the first cancellation request."""
        with self._lock:
            return dict(self._request) if self._request is not None else None


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
    cancellation_token: CancellationToken | None = None
    checkpoint_directory: str | None = None
    handle_keyboard_interrupt: bool = True
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))

    # ---- internal state --------------------------------------------------
    _machine: FrequencyStateMachine | None = field(default=None, repr=False)
    _tracker: DampedSecantTracker | None = field(default=None, repr=False)
    _journal: list[dict] = field(default_factory=list, repr=False)
    _cost_ledger: list[dict] = field(default_factory=list, repr=False)
    _checkpoints: list[StateMachineSnapshot] = field(default_factory=list, repr=False)
    _last_monitor_time: float = field(default=0.0, repr=False)
    _interrupt_handled: bool = field(default=False, repr=False)
    _interrupt_request: dict | None = field(default=None, repr=False)
    _interrupt_checkpoint_saved: bool = field(default=False, repr=False)
    _interrupt_checkpoint_error: str | None = field(default=None, repr=False)
    _safe_hold_confirmed: bool = field(default=False, repr=False)

    def __post_init__(self):
        if self.config is None:
            self.config = FrequencyCalibrationConfig()
        if self.executor is None:
            self.executor = SQCExecutor(
                qubit=self.qubit,
                config=self.config,
                f_target=self.f_target,
            )
        if self.cancellation_token is None:
            self.cancellation_token = CancellationToken()

    def request_cancel(
        self,
        reason: str = "user_requested",
        source: str = "api",
    ) -> bool:
        """Request cooperative cancellation from any thread.

        The request prevents the next science command from starting.  A
        synchronous command already inside ``executor.execute`` must return
        before the runtime can enter SafeStop.
        """
        return self.cancellation_token.request(reason=reason, source=source)

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
        if self._machine is None:
            self._machine = FrequencyStateMachine(
                config=self.config,
                f_target=self.f_target,
            )
            self._machine.start()
            self._journal = []
            self._cost_ledger = []
            self._last_monitor_time = 0.0
        elif self._machine.run_status == RunStatus.READY:
            self._machine.start()

        try:
            self._run_event_loop()
        except KeyboardInterrupt:
            if not self.handle_keyboard_interrupt:
                raise
            self.request_cancel(
                reason="keyboard_interrupt",
                source="keyboard",
            )
            self._process_cancel_request()
            self._run_event_loop()

        elapsed = _time.time() - t_start
        completion = next(
            (
                row for row in self._machine._transition_log
                if row.get("from") == FrequencyState.VERIFY.value
                and row.get("to") == FrequencyState.LOCK.value
            ),
            None,
        )
        return {
            "state": self._machine.state.value,
            "run_status": self._machine.run_status.value,
            "f_final": self._machine._f_hat,
            "candidate_bias": self._machine._candidate_bias,
            "transition_log": self._machine._transition_log,
            "journal": self._journal,
            "cost_ledger": self._cost_ledger,
            "cost_summary": self._cost_summary(),
            "n_commands": self._command_count(),
            "elapsed": elapsed,
            "checkpoints": len(self._checkpoints),
            "safe_hold_confirmed": self._safe_hold_confirmed,
            "hold_target_met": self._machine._hold_target_met,
            "hold_samples": self._machine._hold_samples,
            "hold_passes": self._machine._hold_passes,
            "interrupted": self._interrupt_handled,
            "interrupt": self._interrupt_request,
            "interrupt_checkpoint_saved": self._interrupt_checkpoint_saved,
            "interrupt_checkpoint_error": self._interrupt_checkpoint_error,
            "run_id": self.run_id,
            "backend_provenance": self._backend_provenance(),
            "completion_event": "verify_to_lock" if completion is not None else None,
            "completed_verified_lock": completion is not None,
        }

    def _run_event_loop(self):
        """Drive command/event cycles until the lifecycle is terminal."""
        safe_hold_attempts = 0
        while self._machine.run_status in (RunStatus.RUNNING, RunStatus.CALIBRATED) or (
            self._machine.state == FrequencyState.SAFE_STOP
            and not self._safe_hold_confirmed
            and safe_hold_attempts <= self.config.max_technical_retries
        ):
            if self._process_cancel_request():
                continue

            cmd = self._machine.next_command()

            # Pre-process: for Track commands, use the tracker to refine the proposal
            cmd = self._enrich_command(cmd)
            if self._machine._pending_command_id is not None:
                self._machine.replace_pending_command(cmd)

            preflight_failure = self._machine.validate_command(cmd)
            if preflight_failure is not None:
                state_before = self._machine.state.value
                version_before = self._machine.state_version
                event = MeasurementRejected(
                    command_id=cmd.command_id,
                    state_version=version_before,
                    reason_code=preflight_failure,
                    diagnostics={
                        "stage": "premeasurement",
                        "predicted_detuning": getattr(
                            cmd,
                            "predicted_detuning",
                            None,
                        ),
                        "predicted_detuning_uncertainty": getattr(
                            cmd,
                            "predicted_detuning_uncertainty",
                            None,
                        ),
                        "prediction_guard_complete": getattr(
                            cmd,
                            "prediction_guard_complete",
                            False,
                        ),
                    },
                )
                self._machine.handle(event)
                self._record_journal(
                    cmd,
                    event,
                    state_before=state_before,
                    version_before=version_before,
                )
                continue

            if not isinstance(cmd, SafeHold):
                estimate_cost = getattr(self.executor, "estimate_cost", None)
                cost = estimate_cost(cmd) if callable(estimate_cost) else {}
                if not self._machine.reserve_execution_budget(
                    estimated_shots=cost.get("shots", 0),
                    estimated_solver_calls=cost.get("solver_calls", 0),
                ):
                    cmd = self._machine.next_command()

            if self._process_cancel_request():
                continue

            if (
                isinstance(cmd, MonitorFrequency)
                and self.config.monitor_interval > 0
                and self._last_monitor_time > 0
            ):
                delay = self.config.monitor_interval - (
                    _time.time() - self._last_monitor_time
                )
                if delay > 0:
                    if self.cancellation_token.wait(delay):
                        self._process_cancel_request()
                        continue
                if not self._machine.reserve_execution_budget():
                    cmd = self._machine.next_command()

            if self._process_cancel_request():
                continue

            # Execute
            state_before = self._machine.state.value
            version_before = self._machine.state_version
            event = self.executor.execute(cmd)
            if hasattr(event, "state_version") and event.state_version == 0:
                event.state_version = version_before
            if isinstance(cmd, MonitorFrequency):
                self._last_monitor_time = _time.time()

            # Handle event → state transition
            self._machine.handle(event)
            self._record_journal(
                cmd,
                event,
                state_before=state_before,
                version_before=version_before,
            )

            # A synchronous executor cannot be pre-empted by the runtime.  If
            # cancellation arrived during execute(), retain the returned event
            # for cost/audit purposes, then stop before any next science command.
            if not isinstance(cmd, SafeHold) and self._process_cancel_request():
                continue

            # Terminal status is not enough: send the safe-bias command to the
            # executor and record its acknowledgement before leaving the loop.
            if isinstance(cmd, SafeHold):
                safe_hold_attempts += 1
                self._safe_hold_confirmed = isinstance(event, SafeHoldApplied)
                if self._safe_hold_confirmed:
                    self._persist_interrupt_checkpoint()
                    break
                continue

            # Post-process: update tracker after Track measurements
            self._post_process(cmd, event)

            # Checkpoint periodically
            if len(self._journal) % 10 == 0:
                self._checkpoints.append(self._machine.snapshot())

    def _process_cancel_request(self) -> bool:
        """Convert the shared cancellation flag into one state-machine event."""
        if self._interrupt_handled or not self.cancellation_token.requested:
            return False

        request = self.cancellation_token.snapshot() or {
            "reason": "user_requested",
            "source": "api",
            "requested_at": _time.time(),
        }
        state_before = self._machine.state.value
        version_before = self._machine.state_version
        pending_command = self._machine._pending_command
        self._machine.handle(CancelRequested(**request))
        self._interrupt_handled = True
        self._interrupt_request = request
        self._journal.append(
            self._json_safe(
                {
                    "command": None,
                    "command_id": None,
                    "event": "CancelRequested",
                    "state_before": state_before,
                    "state_after": self._machine.state.value,
                    "state_version": self._machine.state_version,
                    "run_status": self._machine.run_status.value,
                    "interrupt_reason": request["reason"],
                    "interrupt_source": request["source"],
                    "requested_at": request["requested_at"],
                    "interrupted_command": (
                        type(pending_command).__name__
                        if pending_command is not None
                        else None
                    ),
                    "interrupted_command_id": (
                        pending_command.command_id
                        if pending_command is not None
                        else None
                    ),
                    "transition_reason": (
                        self._machine._transition_log[-1]["reason"]
                        if self._machine.state_version != version_before
                        else None
                    ),
                }
            )
        )
        self._checkpoints.append(self._machine.snapshot())
        self._persist_interrupt_checkpoint()
        return True

    def _command_count(self) -> int:
        """Count executed command journal rows, excluding spontaneous events."""
        return sum(entry.get("command") is not None for entry in self._journal)

    def _persist_interrupt_checkpoint(self):
        """Best-effort persistence for an acknowledged cooperative interrupt."""
        if not self._interrupt_handled or self.checkpoint_directory is None:
            return
        try:
            self.save_run(self.checkpoint_directory)
        except Exception as exc:  # SafeHold must still be attempted.
            self._interrupt_checkpoint_saved = False
            self._interrupt_checkpoint_error = str(exc)
        else:
            self._interrupt_checkpoint_saved = True
            self._interrupt_checkpoint_error = None

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
                expected_sensitivity_sign=self.config.expected_sensitivity_sign,
                max_bias_step=self.config.max_bias_step,
                converge_streak=1,  # tracker checks within-step; SM checks overall
                max_iter=999,  # managed by SM
                # The FSM owns epsilon_enter and routes candidates to Verify.
                # The tracker must still propose a correction after a reliable
                # Verify miss, so its internal zero-step threshold is final.
                epsilon_f=self.config.epsilon_final,
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
                predicted_detuning = None
                predicted_detuning_uncertainty = None
                prediction_guard_complete = False
                prediction_source = "unavailable"

                if proposal.s_hat is None:
                    if self.config.blind_step_detuning_bound is not None:
                        predicted_detuning = self.config.blind_step_detuning_bound
                        predicted_detuning_uncertainty = 0.0
                        prediction_guard_complete = True
                        prediction_source = "configured_blind_step_bound"
                else:
                    observed_step = snap.V - snap.V_prev
                    if abs(observed_step) > 1e-15 and snap.uncertainty_prev is not None:
                        sensitivity_uncertainty = np.hypot(
                            snap.uncertainty,
                            snap.uncertainty_prev,
                        ) / abs(observed_step)
                        prediction_step = proposal.V_next - snap.V
                        predicted_detuning = 0.0
                        predicted_detuning_uncertainty = float(
                            np.hypot(
                                snap.uncertainty,
                                prediction_step * sensitivity_uncertainty,
                            )
                        )
                        prediction_guard_complete = True
                        prediction_source = "secant_propagation"

                return TrackFrequency(
                    command_id=cmd.command_id,
                    candidate_bias=proposal.V_next,
                    predicted_drive=predicted_drive,
                    predicted_detuning=predicted_detuning,
                    predicted_detuning_uncertainty=predicted_detuning_uncertainty,
                    prediction_guard_complete=prediction_guard_complete,
                    tracker_spec={
                        "s_hat": proposal.s_hat,
                        "step": proposal.step,
                        "prediction_source": prediction_source,
                        "sensitivity_direction_valid": proposal.diagnostics.get(
                            "reason"
                        )
                        != "sensitivity_direction_mismatch",
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

        if not isinstance(event, MeasurementSucceeded):
            return
        if not event.valid or event.ambiguous:
            return
        if self._machine.state == FrequencyState.REACQUIRE:
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
                "uncertainty": snap.uncertainty,
                "uncertainty_prev": snap.uncertainty_prev,
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
                "uncertainty": new_snap.uncertainty,
                "uncertainty_prev": new_snap.uncertainty_prev,
            }

    # ------------------------------------------------------------------
    # Journal
    # ------------------------------------------------------------------

    def _record_journal(
        self,
        cmd: Command,
        event=None,
        *,
        state_before: str | None = None,
        version_before: int | None = None,
    ):
        """Record a journal entry for this command→event cycle."""
        state_after = self._machine.state.value if self._machine is not None else None
        entry = {
            "command": type(cmd).__name__,
            "command_id": cmd.command_id,
            "event": type(event).__name__ if event is not None else None,
            "state_before": state_before,
            "state_after": state_after,
            "state_version": (
                self._machine.state_version if self._machine is not None else None
            ),
            "run_status": (
                self._machine.run_status.value if self._machine is not None else None
            ),
        }
        if hasattr(cmd, "candidate_bias"):
            entry["commanded_bias"] = getattr(cmd, "candidate_bias", None)
        elif hasattr(cmd, "frozen_bias"):
            entry["commanded_bias"] = getattr(cmd, "frozen_bias", None)
        elif hasattr(cmd, "locked_bias"):
            entry["commanded_bias"] = getattr(cmd, "locked_bias", None)

        for attribute in ("predicted_drive", "drive"):
            if hasattr(cmd, attribute):
                entry["commanded_drive"] = getattr(cmd, attribute)
                break

        if event is not None and hasattr(event, "frequency"):
            entry["frequency"] = event.frequency
            entry["uncertainty"] = event.uncertainty
            entry["residual"] = event.frequency - self.f_target
            entry["valid"] = event.valid
            entry["ambiguous"] = event.ambiguous
            entry["method"] = event.method
            entry["shots"] = event.shots
            entry["elapsed_time"] = event.elapsed_time
            entry["applied_bias"] = event.applied_bias
            entry["applied_drive"] = event.applied_drive
            entry["probe_detuning"] = event.probe_detuning
            entry["probe_detuning_uncertainty"] = event.probe_detuning_uncertainty
            entry["confidence_multiplier"] = self.config.confidence_multiplier
            if isinstance(cmd, MonitorFrequency):
                entry["hold_target_met"] = self._machine._hold_target_met
            entry["diagnostics"] = event.diagnostics
        elif event is not None:
            if hasattr(event, "elapsed_time"):
                entry["elapsed_time"] = event.elapsed_time
            if hasattr(event, "reason"):
                entry["failure_reason"] = event.reason
            if hasattr(event, "reason_code"):
                entry["failure_reason"] = event.reason_code

        if (
            self._machine is not None
            and version_before is not None
            and self._machine.state_version != version_before
            and self._machine._transition_log
        ):
            entry["transition_reason"] = self._machine._transition_log[-1]["reason"]

        self._journal.append(self._json_safe(entry))
        self._record_cost(cmd, event, state_before, state_after)

    def _record_cost(self, cmd, event, state_before: str, state_after: str):
        if cmd is None:
            return
        diagnostics = getattr(event, "diagnostics", {}) or {}
        breakdown = diagnostics.get("solver_call_breakdown", {}) or {}
        solver_calls = diagnostics.get("solver_calls", 0)
        if isinstance(event, MeasurementTechnicalFailure):
            solver_calls = event.solver_calls
        self._cost_ledger.append(
            self._json_safe(
                {
                    "command_id": cmd.command_id,
                    "command": type(cmd).__name__,
                    "state_before": state_before,
                    "state_after": state_after,
                    "event": type(event).__name__ if event is not None else None,
                    "success": isinstance(event, MeasurementSucceeded),
                    "shots": int(getattr(event, "shots", 0) or 0),
                    "circuits": int(diagnostics.get("circuits", 0) or 0),
                    "mesolve_calls": int(breakdown.get("mesolve", 0) or 0),
                    "sesolve_calls": int(breakdown.get("sesolve", 0) or 0),
                    "solver_calls": int(solver_calls or 0),
                    "elapsed_time": float(getattr(event, "elapsed_time", 0.0) or 0.0),
                }
            )
        )

    def _cost_summary(self) -> dict:
        keys = ("shots", "circuits", "mesolve_calls", "sesolve_calls", "solver_calls")
        return {
            key: sum(int(row.get(key, 0) or 0) for row in self._cost_ledger)
            for key in keys
        }

    def _backend_provenance(self) -> dict:
        provider = getattr(self.executor, "provenance", None)
        return provider() if callable(provider) else {"backend": type(self.executor).__name__}

    @classmethod
    def _json_safe(cls, value):
        if isinstance(value, dict):
            return {str(key): cls._json_safe(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [cls._json_safe(item) for item in value]
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if hasattr(value, "value") and isinstance(value.value, str):
            return value.value
        return value

    # ------------------------------------------------------------------
    # Persistence: save / load / resume
    # ------------------------------------------------------------------

    def save_journal(self, path: str):
        """Save the journal to a JSONL file."""
        with open(path, "w") as f:
            for entry in self._journal:
                f.write(_json.dumps(entry, default=str) + "\n")

    def save_cost_ledger(self, path: str):
        with open(path, "w") as f:
            for entry in self._cost_ledger:
                f.write(_json.dumps(entry, default=str) + "\n")

    def latest_checkpoint(self) -> StateMachineSnapshot | None:
        """Return the most recent checkpoint, or None."""
        return self._checkpoints[-1] if self._checkpoints else None

    def save_run(self, directory: str):
        """Persist the complete run state to *directory*.

        Writes:
        - ``config.json`` — protocol configuration
        - ``commands.jsonl`` — journal of commands issued
        - ``transitions.jsonl`` — state transition log
        - ``checkpoint.json`` — latest machine snapshot
        - ``result.json`` — final result summary
        """
        import os as _os

        _os.makedirs(directory, exist_ok=True)

        # config
        config_d = asdict(self.config)
        config_d["f_target"] = self.f_target
        config_d["run_id"] = self.run_id
        config_d["backend_provenance"] = self._backend_provenance()
        rng_state = getattr(self.executor, "rng_state", None)
        if callable(rng_state):
            config_d["rng_state"] = self._json_safe(rng_state())
        with open(_os.path.join(directory, "config.json"), "w") as f:
            _json.dump(config_d, f, indent=2, default=str)

        # commands journal
        self.save_journal(_os.path.join(directory, "commands.jsonl"))
        self.save_cost_ledger(_os.path.join(directory, "cost-ledger.jsonl"))

        # transitions
        if self._machine is not None:
            with open(_os.path.join(directory, "transitions.jsonl"), "w") as f:
                for t in self._machine._transition_log:
                    f.write(_json.dumps(t, default=str) + "\n")

        # checkpoint
        if self._machine is not None:
            snap = self._machine.snapshot()
            snap_d = {
                "state": snap.state.value,
                "run_status": snap.run_status.value,
                "state_version": snap.state_version,
                "candidate_bias": snap.candidate_bias,
                "candidate_bias_version": snap.candidate_bias_version,
                "candidate_source_event": snap.candidate_source_event,
                "f_hat": snap.f_hat,
                "uncertainty": snap.uncertainty,
                "tracker_snapshot": snap.tracker_snapshot,
                "verify_streak": snap.verify_streak,
                "monitor_streak": snap.monitor_streak,
                "lock_entry_time": snap.lock_entry_time,
                "lock_accumulated_time": snap.lock_accumulated_time,
                "last_audit_time": snap.last_audit_time,
                "hold_target_met": snap.hold_target_met,
                "hold_samples": snap.hold_samples,
                "hold_passes": snap.hold_passes,
                "budget": asdict(snap.budget),
                "transition_log": snap.transition_log,
                "last_event": snap.last_event,
                "pending_command": snap.pending_command,
                "pending_command_id": snap.pending_command_id,
                "technical_retries": snap.technical_retries,
            }
            with open(_os.path.join(directory, "checkpoint.json"), "w") as f:
                _json.dump(self._json_safe(snap_d), f, indent=2)

        # result
        result_d = {
            "state": self._machine.state.value if self._machine else "unknown",
            "run_status": self._machine.run_status.value
            if self._machine
            else "unknown",
            "f_final": self._machine._f_hat if self._machine else None,
            "candidate_bias": self._machine._candidate_bias if self._machine else None,
            "n_commands": self._command_count(),
            "safe_hold_confirmed": self._safe_hold_confirmed,
            "hold_target_met": (
                self._machine._hold_target_met if self._machine else None
            ),
            "hold_samples": self._machine._hold_samples if self._machine else 0,
            "hold_passes": self._machine._hold_passes if self._machine else 0,
            "interrupted": self._interrupt_handled,
            "interrupt": self._interrupt_request,
            "interrupt_checkpoint_saved": self._interrupt_checkpoint_saved,
            "interrupt_checkpoint_error": self._interrupt_checkpoint_error,
            "cost_summary": self._cost_summary(),
        }
        with open(_os.path.join(directory, "result.json"), "w") as f:
            _json.dump(result_d, f, indent=2, default=str)

    @classmethod
    def load_run(
        cls,
        directory: str,
        qubit: object,
        executor: SQCExecutor | None = None,
    ) -> "FrequencyCalibrationRuntime":
        """Restore a runtime from a previously saved run directory.

        The restored runtime can be resumed by calling :meth:`run` — the
        state machine will continue from the checkpoint.
        """
        import os as _os

        # config
        with open(_os.path.join(directory, "config.json")) as f:
            config_d = _json.load(f)
        f_target = config_d.pop("f_target")
        run_id = config_d.pop("run_id", str(_uuid.uuid4()))
        config_d.pop("backend_provenance", None)
        saved_rng_state = config_d.pop("rng_state", None)
        config_fields = {item.name for item in fields(FrequencyCalibrationConfig)}
        config = FrequencyCalibrationConfig(
            **{key: value for key, value in config_d.items() if key in config_fields}
        )

        runtime = cls(
            qubit=qubit,
            f_target=f_target,
            config=config,
            executor=executor,
            checkpoint_directory=directory,
            run_id=run_id,
        )
        restore_rng = getattr(runtime.executor, "restore_rng_state", None)
        if saved_rng_state is not None and callable(restore_rng):
            restore_rng(saved_rng_state)

        journal_path = _os.path.join(directory, "commands.jsonl")
        if _os.path.exists(journal_path):
            with open(journal_path) as f:
                runtime._journal = [_json.loads(line) for line in f if line.strip()]

        ledger_path = _os.path.join(directory, "cost-ledger.jsonl")
        if _os.path.exists(ledger_path):
            with open(ledger_path) as f:
                runtime._cost_ledger = [
                    _json.loads(line) for line in f if line.strip()
                ]

        result_path = _os.path.join(directory, "result.json")
        if _os.path.exists(result_path):
            with open(result_path) as f:
                result_d = _json.load(f)
            runtime._safe_hold_confirmed = bool(
                result_d.get("safe_hold_confirmed", False)
            )
            runtime._interrupt_handled = bool(result_d.get("interrupted", False))
            runtime._interrupt_request = result_d.get("interrupt")
            runtime._interrupt_checkpoint_saved = bool(
                result_d.get("interrupt_checkpoint_saved", False)
            )
            runtime._interrupt_checkpoint_error = result_d.get(
                "interrupt_checkpoint_error"
            )

        # checkpoint
        ckpt_path = _os.path.join(directory, "checkpoint.json")
        if _os.path.exists(ckpt_path):
            with open(ckpt_path) as f:
                ckpt_d = _json.load(f)
            budget_d = ckpt_d.get("budget", {})
            budget_fields = {item.name for item in fields(Budget)}
            budget = Budget(
                **{
                    key: value
                    for key, value in budget_d.items()
                    if key in budget_fields
                }
            )
            snap = StateMachineSnapshot(
                state=FrequencyState(ckpt_d["state"]),
                run_status=RunStatus(ckpt_d["run_status"]),
                state_version=ckpt_d.get("state_version", 0),
                candidate_bias=ckpt_d.get("candidate_bias"),
                candidate_bias_version=ckpt_d.get("candidate_bias_version", 0),
                candidate_source_event=ckpt_d.get("candidate_source_event", ""),
                f_hat=ckpt_d.get("f_hat"),
                uncertainty=ckpt_d.get("uncertainty", 0.0),
                tracker_snapshot=ckpt_d.get("tracker_snapshot"),
                verify_streak=ckpt_d.get("verify_streak", 0),
                monitor_streak=ckpt_d.get("monitor_streak", 0),
                lock_entry_time=ckpt_d.get("lock_entry_time", 0.0),
                lock_accumulated_time=ckpt_d.get("lock_accumulated_time", 0.0),
                last_audit_time=ckpt_d.get("last_audit_time", 0.0),
                hold_target_met=ckpt_d.get("hold_target_met"),
                hold_samples=ckpt_d.get("hold_samples", 0),
                hold_passes=ckpt_d.get("hold_passes", 0),
                budget=budget,
                transition_log=ckpt_d.get("transition_log", []),
                last_event=ckpt_d.get("last_event"),
                pending_command=ckpt_d.get("pending_command"),
                pending_command_id=ckpt_d.get("pending_command_id"),
                technical_retries=ckpt_d.get("technical_retries", 0),
            )
            runtime._machine = FrequencyStateMachine(
                config=runtime.config,
                f_target=runtime.f_target,
            )
            runtime._machine.restore(snap)
            runtime._checkpoints.append(snap)

        return runtime
