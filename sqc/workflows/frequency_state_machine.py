"""sqc.workflows.frequency_state_machine — event-driven frequency calibration
state machine (V2).

Implements the ``Acquire → Track → Verify → Lock`` state machine with
``Reacquire`` recovery and ``SafeStop`` safety hold, per the
*Frequency Calibration State Machine V2* specification.

The machine is **pure** — no I/O, no QuTiP, no measurement hardware.  Usage::

    machine = FrequencyStateMachine(config)
    while machine.run_status == RunStatus.RUNNING:
        cmd = machine.next_command()
        event = executor.execute(cmd)      # external I/O
        machine.handle(event)

Design rules (hard constraints)
-------------------------------
1. Only **Track** may propose a new working bias.
2. **Verify** freezes the candidate bias from entry to exit.
3. **Lock** monitors drift; it never adjusts bias directly.
   Small drift → Verify; large jump → Reacquire.
4. Periodic Ramsey audits in Lock go through ``Lock → Verify``, never
   duplicated inside Lock.
5. Analytic ``f(Φ)`` is a simulation oracle only — it never feeds the
   transition reducer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional
import time as _time
import uuid as _uuid


# ===================================================================
# Enums
# ===================================================================


class FrequencyState(str, Enum):
    """Calibration protocol state."""
    ACQUIRE = "acquire"
    TRACK = "track"
    VERIFY = "verify"
    LOCK = "lock"
    REACQUIRE = "reacquire"
    SAFE_STOP = "safe_stop"


class RunStatus(str, Enum):
    """Run-lifecycle status — orthogonal to the protocol state."""
    READY = "ready"
    RUNNING = "running"
    CALIBRATED = "calibrated"       # first Verify→Lock transition
    COMPLETED = "completed"         # bounded sim stopped by policy
    SAFE_STOPPED = "safe_stopped"
    FAILED = "failed"


class ReasonCode(str, Enum):
    """Stable reason codes for every transition and event."""
    TARGET_CANDIDATE = "TARGET_CANDIDATE"
    TARGET_VERIFIED = "TARGET_VERIFIED"
    TARGET_NOT_VERIFIED = "TARGET_NOT_VERIFIED"
    LOCAL_RANGE_LOST = "LOCAL_RANGE_LOST"
    INVERSE_FAILED = "INVERSE_FAILED"
    BRANCH_DISCONTINUITY = "BRANCH_DISCONTINUITY"
    SENSITIVITY_INVALID = "SENSITIVITY_INVALID"
    CONFIDENCE_INSUFFICIENT = "CONFIDENCE_INSUFFICIENT"
    DRIFT_SUSPECTED = "DRIFT_SUSPECTED"
    AUDIT_DUE = "AUDIT_DUE"
    LARGE_FREQUENCY_JUMP = "LARGE_FREQUENCY_JUMP"
    REFERENCE_LOST = "REFERENCE_LOST"
    TECHNICAL_TIMEOUT = "TECHNICAL_TIMEOUT"
    INTERLOCK = "INTERLOCK"
    BUDGET_EXHAUSTED = "BUDGET_EXHAUSTED"
    REACQUIRE_LIMIT = "REACQUIRE_LIMIT"
    CANCELLED = "CANCELLED"
    # Additional codes for internal use
    CANDIDATE_REACHED = "CANDIDATE_REACHED"
    SEED_RELIABLE = "SEED_RELIABLE"
    CORRECTION_REQUIRED = "CORRECTION_REQUIRED"
    MONITOR_CLEAR = "MONITOR_CLEAR"
    MONITOR_SUSPECT = "MONITOR_SUSPECT"
    VERIFY_PASSED = "VERIFY_PASSED"


# ===================================================================
# Protocol configuration
# ===================================================================


@dataclass
class FrequencyCalibrationConfig:
    """All thresholds, budgets, and policy parameters for one calibration run.

    Thresholds follow the hierarchy:

        epsilon_hold < epsilon_final < epsilon_enter < Delta_val

    where the physical *target* is ``epsilon_hold`` (the long-term stability
    bandwidth), but state transitions use the calibrated monitor thresholds.
    """

    # ---- physical target thresholds (angular, rad·GHz) ------------------
    epsilon_enter: float = 2.0 * 3.141592653589793 * 5e-3   # 5 MHz — enter Verify
    epsilon_final: float = 2.0 * 3.141592653589793 * 1e-4   # 100 kHz — pass Verify
    epsilon_hold: float = 2.0 * 3.141592653589793 * 1e-5    # 10 kHz — long-term goal
    Delta_val: float = 2.0 * 3.141592653589793 * 2e-2       # 20 MHz — validity window

    # ---- acquisition ---------------------------------------------------
    sigma_acquire_max: float = 2.0 * 3.141592653589793 * 1e-3   # 1 MHz max acq uncert
    max_acq_retries: int = 3

    # ---- track ---------------------------------------------------------
    damping: float = 0.8
    first_bias_step: float = 0.01
    max_bias_step: float = 0.02
    guard_margin: float = 2.0 * 3.141592653589793 * 1e-3  # 1 MHz pred guard
    S_min: float = 1.0     # min |sensitivity| (rad·GHz / Φ₀)
    S_max: float = 1e4     # max |sensitivity|

    # ---- verify ---------------------------------------------------------
    N_verify: int = 2               # consecutive passes required
    max_verify_attempts_per_episode: int = 5
    max_verify_shots: int = 10_000
    bias_freeze_tolerance: float = 1e-12   # Φ₀ — bit-level freeze

    # ---- lock monitor ---------------------------------------------------
    epsilon_mon_clear: float = 2.0 * 3.141592653589793 * 2e-4    # 200 kHz
    epsilon_mon_suspect: float = 2.0 * 3.141592653589793 * 1e-3  # 1 MHz
    Delta_mon_reacquire: float = 2.0 * 3.141592653589793 * 1e-2  # 10 MHz
    N_mon_suspect: int = 3            # grey-zone streak → Verify
    audit_interval: float = 0.0       # seconds between Ramsey audits (0→disabled)
    monitor_interval: float = 0.0     # seconds between monitor checks (0→every loop)

    # ---- budgets -------------------------------------------------------
    max_commands: int = 200
    max_wall_time: float = 3600.0     # seconds
    max_shots: int = 1_000_000
    max_solver_calls: int = 100_000
    max_reacquire_attempts: int = 5
    stop_after_lock_cycles: int = 0   # 0 → run indefinitely (long-run mode)

    # ---- technical retries ----------------------------------------------
    max_technical_retries: int = 3


# ===================================================================
# Budget
# ===================================================================


@dataclass
class Budget:
    """Mutable budget tracker — consumed by the runtime, checked by guards."""

    commands_issued: int = 0
    wall_time_start: float = 0.0
    total_shots: int = 0
    solver_calls: int = 0
    reacquire_attempts: int = 0
    verify_attempts_this_episode: int = 0
    lock_cycles_completed: int = 0

    def consume_command(self):
        self.commands_issued += 1

    def consume_shots(self, n: int):
        self.total_shots += n

    def consume_solver_calls(self, n: int):
        self.solver_calls += n

    def record_reacquire(self):
        self.reacquire_attempts += 1

    def record_verify_attempt(self):
        self.verify_attempts_this_episode += 1

    def reset_verify_episode(self):
        self.verify_attempts_this_episode = 0

    def record_lock_cycle(self):
        self.lock_cycles_completed += 1

    def check(self, config: FrequencyCalibrationConfig) -> ReasonCode | None:
        """Return the first budget-exhausted reason, or None."""
        if self.commands_issued >= config.max_commands:
            return ReasonCode.BUDGET_EXHAUSTED
        if self.total_shots >= config.max_shots:
            return ReasonCode.BUDGET_EXHAUSTED
        if self.solver_calls >= config.max_solver_calls:
            return ReasonCode.BUDGET_EXHAUSTED
        if self.reacquire_attempts > config.max_reacquire_attempts:
            return ReasonCode.REACQUIRE_LIMIT
        if self.verify_attempts_this_episode > config.max_verify_attempts_per_episode:
            return ReasonCode.BUDGET_EXHAUSTED
        elapsed = _time.time() - self.wall_time_start if self.wall_time_start > 0 else 0.0
        if elapsed >= config.max_wall_time:
            return ReasonCode.BUDGET_EXHAUSTED
        if 0 < config.stop_after_lock_cycles <= self.lock_cycles_completed:
            return ReasonCode.BUDGET_EXHAUSTED
        return None


# ===================================================================
# Commands
# ===================================================================


@dataclass
class AcquireFrequency:
    """Issue a wide-range Ramsey acquisition.

    Parameters
    ----------
    role : str
        ``"acquire"`` (initial) or ``"reacquire"`` (recovery).
    """
    command_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    role: str = "acquire"  # "acquire" | "reacquire"


@dataclass
class TrackFrequency:
    """Issue a single transient tracking step.

    Parameters
    ----------
    candidate_bias : float
        Bias at which to measure (Φ₀).
    predicted_drive : float
        Drive frequency (angular, rad·GHz) predicted by the tracker.
    """
    command_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    candidate_bias: float = 0.0
    predicted_drive: float = 0.0
    tracker_spec: dict = field(default_factory=dict)


@dataclass
class VerifyFrequency:
    """Issue an independent Ramsey verification at the frozen candidate bias.

    Parameters
    ----------
    frozen_bias : float
        Bias frozen at Verify entry (Φ₀).
    drive : float
        Target drive frequency (angular, rad·GHz) — typically ``f_target``.
    """
    command_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    frozen_bias: float = 0.0
    drive: float = 0.0
    verifier_spec: dict = field(default_factory=dict)


@dataclass
class MonitorFrequency:
    """Issue a low-cost Lock monitor check.

    Parameters
    ----------
    locked_bias : float
        The locked bias (Φ₀).
    """
    command_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    locked_bias: float = 0.0
    monitor_spec: dict = field(default_factory=dict)


@dataclass
class SafeHold:
    """Issue a safe-hold command — set bias to a known safe value and stop."""
    command_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    bias: float = 0.0
    reason: ReasonCode = ReasonCode.CANCELLED


# Command union type
Command = AcquireFrequency | TrackFrequency | VerifyFrequency | MonitorFrequency | SafeHold


# ===================================================================
# Events
# ===================================================================


@dataclass
class MeasurementSucceeded:
    """A science measurement completed successfully.

    Carries the frequency estimate, its uncertainty, and the raw data
    reference so the state machine can evaluate guards.
    """
    command_id: str
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    state_version: int = 0
    frequency: float = 0.0           # angular (rad·GHz)
    uncertainty: float = 0.0
    valid: bool = True
    ambiguous: bool = False
    method: str = ""
    shots: int = 0
    elapsed_time: float = 0.0
    applied_bias: float = 0.0
    applied_drive: float = 0.0
    diagnostics: dict = field(default_factory=dict)


@dataclass
class MeasurementTechnicalFailure:
    """A measurement failed for technical reasons (timeout, hardware, etc.)."""
    command_id: str
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    state_version: int = 0
    reason: str = "timeout"
    elapsed_time: float = 0.0


@dataclass
class MeasurementRejected:
    """A measurement succeeded technically but was rejected by post-hoc
    quality guards (e.g. fit failed, contrast too low)."""
    command_id: str
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    state_version: int = 0
    reason_code: ReasonCode = ReasonCode.CONFIDENCE_INSUFFICIENT
    diagnostics: dict = field(default_factory=dict)


@dataclass
class TimerElapsed:
    """A scheduled timer fired (e.g. audit interval)."""
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    state_version: int = 0
    timer_type: str = "AUDIT_DUE"   # "AUDIT_DUE" | "MONITOR"


@dataclass
class InterlockTriggered:
    """An external interlock was triggered."""
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    state_version: int = 0
    description: str = ""


@dataclass
class BudgetExhausted:
    """A budget was exhausted."""
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    state_version: int = 0
    reason_code: ReasonCode = ReasonCode.BUDGET_EXHAUSTED


@dataclass
class CancelRequested:
    """External cancel requested."""
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    state_version: int = 0


@dataclass
class ControlApplied:
    """Confirmation that a bias/drive command was applied."""
    command_id: str
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    state_version: int = 0
    applied_bias: float = 0.0
    applied_drive: float = 0.0


@dataclass
class SafeHoldApplied:
    """Confirmation that safe hold is active."""
    command_id: str
    run_id: str = field(default_factory=lambda: str(_uuid.uuid4()))
    state_version: int = 0
    applied_bias: float = 0.0


# Event union type
Event = (
    MeasurementSucceeded | MeasurementTechnicalFailure | MeasurementRejected
    | TimerElapsed | InterlockTriggered | BudgetExhausted
    | CancelRequested | ControlApplied | SafeHoldApplied
)


# ===================================================================
# State machine snapshot
# ===================================================================


@dataclass
class StateMachineSnapshot:
    """Full serialisable state for checkpoint / replay."""
    state: FrequencyState = FrequencyState.ACQUIRE
    run_status: RunStatus = RunStatus.READY
    state_version: int = 0

    # candidate tracking (frozen on Verify entry)
    candidate_bias: float | None = None
    candidate_bias_version: int = 0
    candidate_source_event: str = ""

    # frequency estimates
    f_hat: float | None = None          # most recent frequency estimate (rad·GHz)
    uncertainty: float = 0.0

    # tracking state
    tracker_snapshot: dict | None = None  # serialised TrackSnapshot

    # verify state
    verify_streak: int = 0

    # lock state
    lock_entry_time: float = 0.0
    lock_accumulated_time: float = 0.0
    monitor_streak: int = 0              # consecutive suspect counts
    last_audit_time: float = 0.0

    # budget
    budget: Budget = field(default_factory=Budget)

    # history
    transition_log: list[dict] = field(default_factory=list)
    last_event: dict | None = None


# ===================================================================
# Guard helpers (pure functions)
# ===================================================================


def _check_candidate(
    uncertainty: float, error: float, config: FrequencyCalibrationConfig,
) -> bool:
    """Candidate condition: U <= epsilon_enter."""
    U = abs(error) + uncertainty
    return U <= config.epsilon_enter


def _check_verified(
    uncertainty: float, error: float, config: FrequencyCalibrationConfig,
) -> bool:
    """Final verification condition: U <= epsilon_final."""
    U = abs(error) + uncertainty
    return U <= config.epsilon_final


def _check_local_validity(
    valid: bool, ambiguous: bool, uncertainty: float,
    error: float, config: FrequencyCalibrationConfig,
) -> bool:
    """Check whether the local estimate is still trustworthy."""
    if not valid or ambiguous:
        return False
    U = abs(error) + uncertainty
    return U <= config.Delta_val


def _check_monitor_clear(
    uncertainty: float, error: float, config: FrequencyCalibrationConfig,
) -> bool:
    """Monitor says 'clear' — keep Lock."""
    U = abs(error) + uncertainty
    return U <= config.epsilon_mon_clear


def _check_monitor_suspect(
    uncertainty: float, error: float, config: FrequencyCalibrationConfig,
) -> bool:
    """Monitor says 'suspect' — grey zone."""
    U = abs(error) + uncertainty
    return config.epsilon_mon_clear < U <= config.epsilon_mon_suspect


def _check_monitor_reacquire(
    uncertainty: float, error: float, config: FrequencyCalibrationConfig,
) -> bool:
    """Monitor says 'large jump' — go to Reacquire."""
    U = abs(error) + uncertainty
    return U >= config.Delta_mon_reacquire


# ===================================================================
# Frequency State Machine
# ===================================================================


class FrequencyStateMachine:
    """Event-driven frequency calibration state machine.

    Parameters
    ----------
    config : FrequencyCalibrationConfig
        Thresholds, budgets, and policy parameters.
    f_target : float
        Target qubit frequency (angular, rad·GHz).
    """

    def __init__(
        self,
        config: FrequencyCalibrationConfig | None = None,
        f_target: float = 0.0,
    ):
        self.config = config or FrequencyCalibrationConfig()
        self.f_target = float(f_target)

        # --- mutable state ---
        self._state = FrequencyState.ACQUIRE
        self._run_status = RunStatus.READY
        self._state_version = 0
        self._candidate_bias: float | None = None
        self._candidate_bias_version = 0
        self._candidate_source_event = ""
        self._f_hat: float | None = None
        self._uncertainty: float = 0.0
        self._verify_streak: int = 0
        self._lock_entry_time: float = 0.0
        self._lock_accumulated_time: float = 0.0
        self._monitor_streak: int = 0
        self._last_audit_time: float = 0.0
        self._budget = Budget()
        self._transition_log: list[dict] = []
        self._last_event: dict | None = None

        # --- pending command tracking ---
        self._pending_command: Command | None = None
        self._pending_command_id: str | None = None

        # --- tracker snapshot (serialised for replay) ---
        self._tracker_snapshot: dict | None = None

        # --- technical retry counter ---
        self._technical_retries: int = 0

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def state(self) -> FrequencyState:
        return self._state

    @property
    def run_status(self) -> RunStatus:
        return self._run_status

    @property
    def state_version(self) -> int:
        return self._state_version

    @property
    def candidate_bias(self) -> float | None:
        return self._candidate_bias

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def start(self):
        """Transition from READY to RUNNING and start the first command."""
        if self._run_status != RunStatus.READY:
            raise RuntimeError(f"cannot start: run_status={self._run_status}")
        self._run_status = RunStatus.RUNNING
        self._budget.wall_time_start = _time.time()
        self._log_transition(
            from_state=FrequencyState.ACQUIRE,
            to_state=FrequencyState.ACQUIRE,
            reason=ReasonCode.SEED_RELIABLE,
            note="run started",
        )

    def snapshot(self) -> StateMachineSnapshot:
        """Return a serialisable snapshot for checkpoint / replay."""
        return StateMachineSnapshot(
            state=self._state,
            run_status=self._run_status,
            state_version=self._state_version,
            candidate_bias=self._candidate_bias,
            candidate_bias_version=self._candidate_bias_version,
            candidate_source_event=self._candidate_source_event,
            f_hat=self._f_hat,
            uncertainty=self._uncertainty,
            tracker_snapshot=self._tracker_snapshot,
            verify_streak=self._verify_streak,
            lock_entry_time=self._lock_entry_time,
            lock_accumulated_time=self._lock_accumulated_time,
            monitor_streak=self._monitor_streak,
            last_audit_time=self._last_audit_time,
            budget=Budget(
                commands_issued=self._budget.commands_issued,
                wall_time_start=self._budget.wall_time_start,
                total_shots=self._budget.total_shots,
                solver_calls=self._budget.solver_calls,
                reacquire_attempts=self._budget.reacquire_attempts,
                verify_attempts_this_episode=self._budget.verify_attempts_this_episode,
                lock_cycles_completed=self._budget.lock_cycles_completed,
            ),
            transition_log=list(self._transition_log),
            last_event=self._last_event,
        )

    def restore(self, snap: StateMachineSnapshot):
        """Restore state from a previously saved snapshot."""
        self._state = snap.state
        self._run_status = snap.run_status
        self._state_version = snap.state_version
        self._candidate_bias = snap.candidate_bias
        self._candidate_bias_version = snap.candidate_bias_version
        self._candidate_source_event = snap.candidate_source_event
        self._f_hat = snap.f_hat
        self._uncertainty = snap.uncertainty
        self._tracker_snapshot = snap.tracker_snapshot
        self._verify_streak = snap.verify_streak
        self._lock_entry_time = snap.lock_entry_time
        self._lock_accumulated_time = snap.lock_accumulated_time
        self._monitor_streak = snap.monitor_streak
        self._last_audit_time = snap.last_audit_time
        self._budget = snap.budget
        self._transition_log = list(snap.transition_log)
        self._last_event = snap.last_event
        self._pending_command = None
        self._pending_command_id = None
        self._technical_retries = 0

    # ------------------------------------------------------------------
    # next_command — the "output" side of the state machine
    # ------------------------------------------------------------------

    def next_command(self) -> Command:
        """Return the next command the executor should run.

        Must be called after :meth:`start`.  Idempotent: returns the same
        command until :meth:`handle` is called with an event that carries
        the matching ``command_id``.
        """
        if self._run_status == RunStatus.READY:
            raise RuntimeError("call start() before next_command()")
        if self._run_status in (RunStatus.COMPLETED, RunStatus.SAFE_STOPPED,
                                RunStatus.FAILED):
            # No more commands — return a no-op SafeHold
            return SafeHold(bias=0.0, reason=ReasonCode.CANCELLED)

        # If a command is already pending, return it again (idempotent)
        if self._pending_command is not None:
            return self._pending_command

        cmd = self._dispatch_command()
        self._pending_command = cmd
        self._pending_command_id = cmd.command_id
        self._budget.consume_command()
        return cmd

    def _dispatch_command(self) -> Command:
        """Produce the command appropriate for the current state."""
        match self._state:
            case FrequencyState.ACQUIRE:
                return AcquireFrequency(role="acquire")
            case FrequencyState.REACQUIRE:
                return AcquireFrequency(role="reacquire")
            case FrequencyState.TRACK:
                # Build TrackFrequency from tracker snapshot
                bias = self._candidate_bias if self._candidate_bias is not None else 0.0
                drive = self._f_hat if self._f_hat is not None else self.f_target
                return TrackFrequency(
                    candidate_bias=bias,
                    predicted_drive=drive,
                    tracker_spec=self._tracker_snapshot or {},
                )
            case FrequencyState.VERIFY:
                bias = self._candidate_bias if self._candidate_bias is not None else 0.0
                return VerifyFrequency(
                    frozen_bias=bias,
                    drive=self.f_target,
                )
            case FrequencyState.LOCK:
                # Check if audit is due
                now = _time.time()
                if (
                    self.config.audit_interval > 0
                    and self._last_audit_time > 0
                    and (now - self._last_audit_time) >= self.config.audit_interval
                ):
                    return VerifyFrequency(
                        frozen_bias=self._candidate_bias or 0.0,
                        drive=self.f_target,
                    )
                # Otherwise, issue monitor
                bias = self._candidate_bias if self._candidate_bias is not None else 0.0
                return MonitorFrequency(locked_bias=bias)
            case FrequencyState.SAFE_STOP:
                return SafeHold(bias=0.0, reason=ReasonCode.CANCELLED)

    # ------------------------------------------------------------------
    # handle — the "input" side of the state machine
    # ------------------------------------------------------------------

    def handle(self, event: Event):
        """Process an event and transition state.

        The event's ``command_id`` must match the currently pending command
        (or be a spontaneous event like ``TimerElapsed``, ``InterlockTriggered``,
        ``CancelRequested``).
        """
        # --- check budget first (global guard) ---
        budget_reason = self._budget.check(self.config)
        if budget_reason is not None:
            was_in_lock = self._state == FrequencyState.LOCK
            self._transition_to(
                FrequencyState.SAFE_STOP, budget_reason,
            )
            if budget_reason == ReasonCode.BUDGET_EXHAUSTED:
                # If we were in Lock, budget expiry is COMPLETED; else FAILED
                self._run_status = (
                    RunStatus.COMPLETED if was_in_lock else RunStatus.FAILED
                )
            return

        # --- spontaneous events (no command_id match needed) ---
        match event:
            case InterlockTriggered():
                self._transition_to(
                    FrequencyState.SAFE_STOP, ReasonCode.INTERLOCK,
                )
                self._run_status = RunStatus.SAFE_STOPPED
                self._clear_pending()
                return
            case CancelRequested():
                self._transition_to(
                    FrequencyState.SAFE_STOP, ReasonCode.CANCELLED,
                )
                self._run_status = RunStatus.SAFE_STOPPED
                self._clear_pending()
                return
            case TimerElapsed(timer_type="AUDIT_DUE"):
                if self._state == FrequencyState.LOCK:
                    self._transition_to(
                        FrequencyState.VERIFY, ReasonCode.AUDIT_DUE,
                    )
                    self._last_audit_time = _time.time()
                self._clear_pending()
                return
            case _:
                pass

        # --- command-matched events ---
        cmd_id = getattr(event, "command_id", None)
        if cmd_id is not None and cmd_id != self._pending_command_id:
            # Stale or mismatched event — ignore (idempotent)
            return

        # --- budget tracking ---
        if isinstance(event, MeasurementSucceeded):
            self._budget.consume_shots(event.shots)
            self._budget.consume_solver_calls(
                event.diagnostics.get("solver_calls", 0)
            )

        # --- dispatch by current state ---
        match self._state:
            case FrequencyState.ACQUIRE:
                self._handle_acquire(event)
            case FrequencyState.TRACK:
                self._handle_track(event)
            case FrequencyState.VERIFY:
                self._handle_verify(event)
            case FrequencyState.LOCK:
                self._handle_lock(event)
            case FrequencyState.REACQUIRE:
                self._handle_reacquire(event)
            case FrequencyState.SAFE_STOP:
                # Already stopped — ignore further events
                pass

        self._record_event(event)

    # ------------------------------------------------------------------
    # Per-state handlers
    # ------------------------------------------------------------------

    def _handle_acquire(self, event: Event):
        match event:
            case MeasurementSucceeded(
                valid=True, ambiguous=False,
                frequency=f, uncertainty=u,
            ) if u <= self.config.sigma_acquire_max:
                self._f_hat = f
                self._uncertainty = u
                self._candidate_bias = event.applied_bias
                self._candidate_bias_version += 1
                self._candidate_source_event = "acquire"
                self._tracker_snapshot = None   # fresh start
                self._technical_retries = 0

                if _check_candidate(u, f - self.f_target, self.config):
                    self._transition_to(FrequencyState.VERIFY, ReasonCode.CANDIDATE_REACHED)
                else:
                    self._transition_to(FrequencyState.TRACK, ReasonCode.SEED_RELIABLE)

            case MeasurementSucceeded(valid=False) | MeasurementSucceeded(ambiguous=True):
                self._handle_technical_failure(event, "acquire")

            case MeasurementTechnicalFailure() | MeasurementRejected():
                self._handle_technical_failure(event, "acquire")

            case _:
                pass  # ignore unexpected events

        self._clear_pending()

    def _handle_track(self, event: Event):
        match event:
            case MeasurementSucceeded(
                valid=True, ambiguous=False,
                frequency=f, uncertainty=u,
            ):
                error = f - self.f_target
                self._f_hat = f
                self._uncertainty = u
                self._candidate_bias = event.applied_bias
                self._technical_retries = 0

                # Guard evaluation (fixed order per §5.2)
                if not _check_local_validity(
                    True, False, u, error, self.config,
                ):
                    # Check specific failure modes for diagnostics
                    U = abs(error) + u
                    if not (self.config.S_min <= abs(event.diagnostics.get("s_hat", 1.0)) <= self.config.S_max):
                        reason = ReasonCode.SENSITIVITY_INVALID
                    elif U > self.config.Delta_val:
                        reason = ReasonCode.LOCAL_RANGE_LOST
                    else:
                        reason = ReasonCode.CONFIDENCE_INSUFFICIENT
                    self._transition_to(FrequencyState.REACQUIRE, reason)
                elif _check_candidate(u, error, self.config):
                    # Candidate reached → freeze and verify
                    self._candidate_bias_version += 1
                    self._candidate_source_event = "track"
                    self._transition_to(FrequencyState.VERIFY, ReasonCode.CANDIDATE_REACHED)
                else:
                    # Continue tracking
                    self._stay(FrequencyState.TRACK, ReasonCode.CORRECTION_REQUIRED)

            case MeasurementSucceeded(valid=False) | MeasurementSucceeded(ambiguous=True):
                self._transition_to(FrequencyState.REACQUIRE, ReasonCode.CONFIDENCE_INSUFFICIENT)

            case MeasurementTechnicalFailure() | MeasurementRejected():
                self._handle_technical_failure(event, "track")

            case _:
                pass

        self._clear_pending()

    def _handle_verify(self, event: Event):
        match event:
            case MeasurementSucceeded(
                valid=True, ambiguous=False,
                frequency=f, uncertainty=u,
            ):
                error = f - self.f_target
                self._f_hat = f
                self._uncertainty = u
                self._budget.record_verify_attempt()
                self._technical_retries = 0

                if _check_verified(u, error, self.config):
                    self._verify_streak += 1
                    if self._verify_streak >= self.config.N_verify:
                        self._transition_to(FrequencyState.LOCK, ReasonCode.VERIFY_PASSED)
                        self._run_status = RunStatus.CALIBRATED
                        self._verify_streak = 0
                        self._budget.reset_verify_episode()
                        self._lock_entry_time = _time.time()
                        self._monitor_streak = 0
                        self._last_audit_time = _time.time()
                    else:
                        self._stay(FrequencyState.VERIFY, ReasonCode.VERIFY_PASSED)
                elif _check_local_validity(True, False, u, error, self.config):
                    # Reliable but not yet at final tolerance
                    self._verify_streak = 0
                    self._transition_to(
                        FrequencyState.TRACK, ReasonCode.CORRECTION_REQUIRED,
                    )
                    # Re-seed tracker with the verifier's absolute frequency
                    self._tracker_snapshot = None
                else:
                    # Ambiguous or outside local recovery range
                    self._verify_streak = 0
                    self._transition_to(
                        FrequencyState.REACQUIRE, ReasonCode.CONFIDENCE_INSUFFICIENT,
                    )

            case MeasurementSucceeded(valid=False) | MeasurementSucceeded(ambiguous=True):
                self._transition_to(FrequencyState.REACQUIRE, ReasonCode.CONFIDENCE_INSUFFICIENT)

            case MeasurementTechnicalFailure() | MeasurementRejected():
                self._handle_technical_failure(event, "verify")

            case _:
                pass

        self._clear_pending()

    def _handle_lock(self, event: Event):
        match event:
            case MeasurementSucceeded(
                valid=True, ambiguous=False,
                frequency=f, uncertainty=u,
            ):
                error = f - self.f_target
                self._f_hat = f
                self._uncertainty = u
                self._technical_retries = 0

                if _check_monitor_clear(u, error, self.config):
                    self._monitor_streak = 0
                    self._stay(FrequencyState.LOCK, ReasonCode.MONITOR_CLEAR)
                elif _check_monitor_suspect(u, error, self.config):
                    self._monitor_streak += 1
                    if self._monitor_streak >= self.config.N_mon_suspect:
                        self._transition_to(FrequencyState.VERIFY, ReasonCode.DRIFT_SUSPECTED)
                        self._monitor_streak = 0
                    else:
                        self._stay(FrequencyState.LOCK, ReasonCode.MONITOR_SUSPECT)
                elif _check_monitor_reacquire(u, error, self.config):
                    self._transition_to(
                        FrequencyState.REACQUIRE, ReasonCode.LARGE_FREQUENCY_JUMP,
                    )
                    self._monitor_streak = 0
                else:
                    # Between suspect and reacquire — go to Verify
                    self._transition_to(FrequencyState.VERIFY, ReasonCode.DRIFT_SUSPECTED)
                    self._monitor_streak = 0

            case MeasurementSucceeded(valid=False) | MeasurementSucceeded(ambiguous=True):
                self._transition_to(
                    FrequencyState.REACQUIRE, ReasonCode.REFERENCE_LOST,
                )

            case MeasurementTechnicalFailure() | MeasurementRejected():
                self._handle_technical_failure(event, "lock")

            case _:
                pass

        # After handling a Lock event, check if it was actually a Verify
        # triggered by the audit timer
        if isinstance(event, MeasurementSucceeded) and self._state == FrequencyState.VERIFY:
            self._last_audit_time = _time.time()

        self._clear_pending()

    def _handle_reacquire(self, event: Event):
        self._budget.record_reacquire()

        match event:
            case MeasurementSucceeded(
                valid=True, ambiguous=False,
                frequency=f, uncertainty=u,
            ) if u <= self.config.sigma_acquire_max:
                self._f_hat = f
                self._uncertainty = u
                self._candidate_bias = event.applied_bias
                self._candidate_bias_version += 1
                self._candidate_source_event = "reacquire"
                self._tracker_snapshot = None
                self._verify_streak = 0
                self._monitor_streak = 0
                self._technical_retries = 0

                if _check_candidate(u, f - self.f_target, self.config):
                    self._transition_to(FrequencyState.VERIFY, ReasonCode.CANDIDATE_REACHED)
                else:
                    self._transition_to(FrequencyState.TRACK, ReasonCode.SEED_RELIABLE)

            case MeasurementSucceeded(valid=False) | MeasurementSucceeded(ambiguous=True):
                self._handle_technical_failure(event, "reacquire")

            case MeasurementTechnicalFailure() | MeasurementRejected():
                self._handle_technical_failure(event, "reacquire")

            case _:
                pass

        self._clear_pending()

    # ------------------------------------------------------------------
    # Technical failure handling
    # ------------------------------------------------------------------

    def _handle_technical_failure(self, event: Event, state_name: str):
        """Retry on technical failures; escalate to SafeStop on exhaustion."""
        self._technical_retries += 1
        max_retries = self.config.max_technical_retries

        # Reacquire has its own budget
        if state_name == "reacquire":
            if self._budget.reacquire_attempts > self.config.max_reacquire_attempts:
                self._transition_to(FrequencyState.SAFE_STOP, ReasonCode.REACQUIRE_LIMIT)
                self._run_status = RunStatus.FAILED
                return

        if self._technical_retries > max_retries:
            self._transition_to(FrequencyState.SAFE_STOP, ReasonCode.TECHNICAL_TIMEOUT)
            self._run_status = RunStatus.FAILED
        else:
            # Stay in current state and retry
            self._stay(self._state, ReasonCode.TECHNICAL_TIMEOUT)

    # ------------------------------------------------------------------
    # Internal: transition helpers
    # ------------------------------------------------------------------

    def _transition_to(self, target: FrequencyState, reason: ReasonCode):
        """Record a state transition."""
        old = self._state
        self._state = target
        self._state_version += 1
        self._log_transition(old, target, reason)

    def _stay(self, state: FrequencyState, reason: ReasonCode):
        """Record a self-loop (state stays the same)."""
        self._state_version += 1
        self._log_transition(state, state, reason)

    def _clear_pending(self):
        """Clear the pending command after handling an event."""
        self._pending_command = None
        self._pending_command_id = None

    def _log_transition(
        self, from_state: FrequencyState, to_state: FrequencyState,
        reason: ReasonCode, note: str = "",
    ):
        """Append a transition record to the log."""
        entry: dict = {
            "version": self._state_version,
            "from": from_state.value,
            "to": to_state.value,
            "reason": reason.value,
            "note": note,
            "candidate_bias": self._candidate_bias,
            "f_hat": self._f_hat,
            "uncertainty": self._uncertainty,
        }
        self._transition_log.append(entry)

    def _record_event(self, event: Event):
        """Record the last event for replay."""
        self._last_event = {
            "type": type(event).__name__,
            "command_id": getattr(event, "command_id", None),
            "state_version": self._state_version,
        }
