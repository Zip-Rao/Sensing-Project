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
import time as _time
from contextlib import AbstractContextManager

import numpy as np
from scipy.optimize import least_squares

from sqc.calibration.frequency import (
    FrequencyMeasurement,
    _fft_peak,
    _solve_cubic_detuning,
)
from sqc.workflows.frequency_state_machine import (
    AcquireFrequency,
    Command,
    Event,
    FrequencyCalibrationConfig,
    MeasurementRejected,
    MeasurementSucceeded,
    MeasurementTechnicalFailure,
    MonitorFrequency,
    ReasonCode,
    SafeHold,
    SafeHoldApplied,
    TrackFrequency,
    VerifyFrequency,
)


def _ramsey_peak_with_uncertainty(
    populations: np.ndarray,
    tau_list: np.ndarray,
    shots_per_circuit: int,
) -> tuple[float | None, float]:
    """Fit a Ramsey fringe and return peak frequency and 1-sigma error in GHz."""
    values = np.asarray(populations, dtype=float)
    tau = np.asarray(tau_list, dtype=float)
    span = float(np.ptp(tau))
    fallback = 1.0 / max(span, 1.0)
    if values.size != tau.size or values.size < 5 or shots_per_circuit <= 0:
        return None, fallback
    dt = float(tau[1] - tau[0])
    peak = _fft_peak(values, dt)
    if peak is None:
        return None, fallback

    variance = np.maximum(values * (1.0 - values), 0.25 / shots_per_circuit)
    sigma = np.sqrt(variance / shots_per_circuit)
    omega0 = 2.0 * np.pi * peak
    design = np.column_stack(
        (np.ones_like(tau), np.cos(omega0 * tau), np.sin(omega0 * tau))
    )
    offset, cosine, sine = np.linalg.lstsq(design, values, rcond=None)[0]

    def residual(params):
        base, c_term, s_term, frequency = params
        phase = 2.0 * np.pi * frequency * tau
        model = base + c_term * np.cos(phase) + s_term * np.sin(phase)
        return (model - values) / sigma

    nyquist = 0.5 / dt
    fit = least_squares(
        residual,
        x0=np.array([offset, cosine, sine, peak]),
        bounds=([-0.5, -1.5, -1.5, 0.0], [1.5, 1.5, 1.5, nyquist]),
        method="trf",
    )
    if not fit.success or fit.jac.shape[0] <= fit.jac.shape[1]:
        return float(peak), fallback
    information = fit.jac.T @ fit.jac
    try:
        covariance = np.linalg.inv(information)
    except np.linalg.LinAlgError:
        return float(peak), fallback
    frequency_sigma = float(np.sqrt(max(covariance[3, 3], 0.0)))
    if not np.isfinite(frequency_sigma) or frequency_sigma <= 0:
        frequency_sigma = fallback
    return float(fit.x[3]), min(frequency_sigma, fallback)


# ===================================================================
# SQC Executor
# ===================================================================


class _SolverCallCounter(AbstractContextManager):
    """Instrument solver entry points used by one synchronous SQC command."""

    def __init__(self):
        self.counts = {"mesolve": 0, "sesolve": 0}
        self._targets = []

    def __enter__(self):
        import qutip
        import sqc.calibration.frequency as frequency_module
        import sqc.reconstruction.kernel as kernel_module

        self._wrap(frequency_module, "mesolve", "mesolve")
        self._wrap(kernel_module, "mesolve", "mesolve")
        self._wrap(qutip, "sesolve", "sesolve")
        return self

    def _wrap(self, module, name: str, kind: str):
        original = getattr(module, name)
        self._targets.append((module, name, original))

        def counted(*args, **kwargs):
            self.counts[kind] += 1
            return original(*args, **kwargs)

        setattr(module, name, counted)

    @property
    def total(self) -> int:
        return sum(self.counts.values())

    def __exit__(self, *exc):
        for module, name, original in reversed(self._targets):
            setattr(module, name, original)
        return False


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
    config: FrequencyCalibrationConfig = field(
        default_factory=FrequencyCalibrationConfig
    )
    f_target: float = 0.0
    pulse_protocol: dict = field(default_factory=dict)
    track_measurement: dict = field(
        default_factory=lambda: {"order": 3, "g3_source": "kernel_full"}
    )
    track_method: str = "transient"

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
        with _SolverCallCounter() as counter:
            match cmd:
                case AcquireFrequency():
                    event = self._execute_acquire(cmd)
                case TrackFrequency():
                    event = self._execute_track(cmd)
                case VerifyFrequency():
                    event = self._execute_verify(cmd)
                case MonitorFrequency():
                    event = self._execute_monitor(cmd)
                case SafeHold():
                    event = self._execute_safe_hold(cmd)
                case _:
                    event = MeasurementTechnicalFailure(
                        command_id=cmd.command_id,
                        reason=f"unknown command type: {type(cmd).__name__}",
                    )
        if isinstance(event, MeasurementSucceeded):
            event.diagnostics["solver_calls"] = counter.total
            event.diagnostics["solver_call_breakdown"] = dict(counter.counts)
        elif isinstance(event, MeasurementTechnicalFailure):
            event.solver_calls = counter.total
            event.diagnostics["solver_call_breakdown"] = dict(counter.counts)
        return event

    def estimate_cost(self, cmd: Command) -> dict[str, int]:
        """Return conservative preflight reservations, not measured cost."""
        if isinstance(cmd, AcquireFrequency):
            meas = self._get_ramsey()
            n = len(meas.tau_list) * (1 if meas.f_artificial is not None else 2)
        elif isinstance(cmd, VerifyFrequency):
            n = 2 * len(self._get_verify().tau_list)
        elif isinstance(cmd, (TrackFrequency, MonitorFrequency)):
            n = (
                2 * len(self._get_verify().tau_list)
                if isinstance(cmd, TrackFrequency) and self.track_method == "ramsey"
                else 2
            )
        else:
            n = 0
        solver_calls = n
        if isinstance(cmd, (TrackFrequency, MonitorFrequency)):
            solver_calls = (
                n if isinstance(cmd, TrackFrequency) and self.track_method == "ramsey"
                else int(self.track_measurement.get("estimated_solver_calls", 200))
            )
        return {"shots": n, "solver_calls": solver_calls}

    # ------------------------------------------------------------------
    # Per-command executors
    # ------------------------------------------------------------------

    def _execute_acquire(self, cmd: AcquireFrequency) -> Event:
        """Wide-range Ramsey acquisition (absolute frequency measurement)."""
        t0 = _time.time()
        try:
            meas = self._get_ramsey()
            # Acquire at zero flux (sweet spot) — wide sweep for absolute f01
            f = meas.measure(flux=cmd.bias)
            elapsed = _time.time() - t0
            n_shots = len(meas.tau_list) * (1 if meas.f_artificial is not None else 2)
            return MeasurementSucceeded(
                command_id=cmd.command_id,
                frequency=f,
                uncertainty=0.0,
                valid=True,
                ambiguous=False,
                method="ramsey",
                shots=0,
                elapsed_time=elapsed,
                applied_bias=cmd.bias,
                diagnostics={
                    "solver_calls": n_shots,
                    "role": cmd.role,
                    "circuits": n_shots,
                    "uncertainty_source": "deterministic_zero",
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
            meas = self._get_track_measurement()
            f = meas.measure(
                flux=cmd.candidate_bias,
                omega_d=cmd.predicted_drive,
            )
            elapsed = _time.time() - t0
            detuning = f - cmd.predicted_drive
            out_of_range = (
                self.config.linear_range is not None
                and abs(detuning) + self.config.guard_margin > self.config.linear_range
            )
            return MeasurementSucceeded(
                command_id=cmd.command_id,
                frequency=f,
                uncertainty=0.0,
                valid=True,
                ambiguous=False,
                method=self.track_method,
                shots=0,
                elapsed_time=elapsed,
                applied_bias=cmd.candidate_bias,
                applied_drive=cmd.predicted_drive,
                probe_detuning=detuning,
                probe_detuning_uncertainty=0.0,
                diagnostics={
                    "solver_calls": 2,
                    "s_hat": cmd.tracker_spec.get("s_hat"),
                    "circuits": (
                        2 * len(meas.tau_list) if self.track_method == "ramsey" else 2
                    ),
                    "detuning": detuning,
                    "out_of_range": out_of_range,
                    "confidence": None,
                    "uncertainty_source": "deterministic_zero",
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
                shots=0,
                elapsed_time=elapsed,
                applied_bias=cmd.frozen_bias,
                applied_drive=cmd.drive,
                diagnostics={
                    "solver_calls": n_shots,
                    "verifier_spec": cmd.verifier_spec,
                    "circuits": n_shots,
                    "uncertainty_source": "deterministic_zero",
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
                shots=0,
                elapsed_time=elapsed,
                applied_bias=cmd.locked_bias,
                applied_drive=self.f_target,
                diagnostics={
                    "solver_calls": 2,
                    "monitor_spec": cmd.monitor_spec,
                    "circuits": 2,
                    "uncertainty_source": "deterministic_zero",
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
                f_artificial=0.1,
                **self.pulse_protocol,
            )
        return self._ramsey_meas

    def _get_transient(self) -> FrequencyMeasurement:
        if self._transient_meas is None:
            measurement_config = {
                key: value
                for key, value in self.track_measurement.items()
                if key != "estimated_solver_calls"
            }
            self._transient_meas = FrequencyMeasurement(
                qubit=self.qubit,
                method="transient",
                **measurement_config,
                **self.pulse_protocol,
            )
        return self._transient_meas

    def _get_track_measurement(self) -> FrequencyMeasurement:
        if self.track_method == "transient":
            return self._get_transient()
        if self.track_method == "ramsey":
            return self._get_verify()
        raise ValueError("track_method must be 'transient' or 'ramsey'")

    def _get_verify(self) -> FrequencyMeasurement:
        if self._verify_meas is None:
            self._verify_meas = FrequencyMeasurement(
                qubit=self.qubit,
                method="ramsey",
                f_artificial=None,
                **self.pulse_protocol,
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

    def estimate_cost(self, cmd: Command) -> dict[str, int]:
        """Delegate preflight cost estimation to the wrapped executor."""
        return self._inner.estimate_cost(cmd)

    def execute(self, cmd: Command) -> Event:
        """Execute with fault injection."""
        # Technical failure injection
        if (
            self._fail_on_state is not None
            and self._machine_state == self._fail_on_state
        ):
            if self._fail_remaining > 0:
                self._fail_remaining -= 1
                return MeasurementTechnicalFailure(
                    command_id=cmd.command_id,
                    reason="injected fault",
                )

        # Rejection injection
        if (
            self._reject_on_state is not None
            and self._machine_state == self._reject_on_state
        ):
            if self._fail_remaining > 0:
                self._fail_remaining -= 1
                return MeasurementRejected(
                    command_id=cmd.command_id,
                    reason_code=ReasonCode.CONFIDENCE_INSUFFICIENT,
                )

        # Ambiguity injection
        if (
            self._ambiguous_on_state is not None
            and self._machine_state == self._ambiguous_on_state
        ):
            if self._fail_remaining > 0:
                self._fail_remaining -= 1
                # Run the real measurement but mark it ambiguous
                event = self._inner.execute(cmd)
                if isinstance(event, MeasurementSucceeded):
                    event.ambiguous = True
                return event

        return self._inner.execute(cmd)


@dataclass
class FiniteShotSQCExecutor(SQCExecutor):
    """SQC executor with branch-wise binomial population sampling."""

    master_seed: int = 0
    shots_per_circuit: int = 1024
    shots_by_role: dict[str, int] = field(default_factory=dict)
    assignment_p01: float = 0.0
    assignment_p10: float = 0.0
    _rngs: dict = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self):
        if self.shots_per_circuit <= 0:
            raise ValueError("shots_per_circuit must be positive")
        valid_roles = {"acquire", "track", "verify", "monitor"}
        unknown = set(self.shots_by_role) - valid_roles
        if unknown:
            raise ValueError(f"unknown finite-shot roles: {sorted(unknown)}")
        if any(value <= 0 for value in self.shots_by_role.values()):
            raise ValueError("role-specific shots per circuit must be positive")
        if not (0 <= self.assignment_p01 < 1 and 0 <= self.assignment_p10 < 1):
            raise ValueError("assignment errors must lie in [0, 1)")
        sequence = np.random.SeedSequence(self.master_seed)
        names = ("acquire", "track", "verify", "monitor", "process")
        self._rngs = {
            name: np.random.default_rng(child)
            for name, child in zip(names, sequence.spawn(len(names)))
        }

    def _shots(self, role: str) -> int:
        return int(self.shots_by_role.get(role, self.shots_per_circuit))

    def estimate_cost(self, cmd: Command) -> dict[str, int]:
        cost = super().estimate_cost(cmd)
        if isinstance(cmd, AcquireFrequency):
            role = "acquire"
        elif isinstance(cmd, TrackFrequency):
            role = "track"
        elif isinstance(cmd, VerifyFrequency):
            role = "verify"
        elif isinstance(cmd, MonitorFrequency):
            role = "monitor"
        else:
            return cost
        return {**cost, "shots": cost["shots"] * self._shots(role)}

    def _sample_population(self, value, role: str):
        p = np.clip(np.asarray(value, dtype=float), 0.0, 1.0)
        observed = self.assignment_p01 + (1 - self.assignment_p01 - self.assignment_p10) * p
        shots = self._shots(role)
        return self._rngs[role].binomial(shots, observed) / shots

    def _sample_transient(self, cmd, role: str) -> MeasurementSucceeded:
        meas = self._get_transient()
        bias = cmd.candidate_bias if isinstance(cmd, TrackFrequency) else cmd.locked_bias
        drive = cmd.predicted_drive if isinstance(cmd, TrackFrequency) else self.f_target
        details = meas.measure_details(flux=bias, omega_d=drive)
        plus = float(self._sample_population(details.populations["plus_x"], role))
        minus = float(self._sample_population(details.populations["minus_x"], role))
        p_diff = (plus - minus) / 2.0
        G1 = float(details.estimator["G1"])
        G3 = float(details.estimator["G3"])
        delta = _solve_cubic_detuning(p_diff, G1, G3) if G3 else p_diff / G1
        frequency = drive - delta
        shots = self._shots(role)
        sigma_p = np.sqrt(
            (plus * (1 - plus) + minus * (1 - minus))
            / (4 * shots)
        )
        uncertainty = abs(sigma_p / G1)
        return MeasurementSucceeded(
            command_id=cmd.command_id, frequency=frequency,
            uncertainty=uncertainty, method="transient_finite_shot",
            shots=2 * shots, applied_bias=bias,
            applied_drive=drive, probe_detuning=frequency - drive,
            probe_detuning_uncertainty=uncertainty,
            diagnostics={"circuits": 2, "sampling_role": role,
                         "uncertainty_source": "binomial_delta_method"},
        )

    def _sample_ramsey(self, cmd, role: str) -> MeasurementSucceeded:
        meas = self._get_ramsey() if isinstance(cmd, AcquireFrequency) else self._get_verify()
        bias = cmd.bias if isinstance(cmd, AcquireFrequency) else cmd.frozen_bias
        if isinstance(cmd, AcquireFrequency):
            drive = self.qubit.frequency
        elif role == "track":
            drive = cmd.drive
        else:
            drive = self.f_target
        details = meas.measure_details(flux=bias, omega_d=drive)
        shots = self._shots(role)
        plus = self._sample_population(details.populations["plus"], role)
        fp, sigma_fp = _ramsey_peak_with_uncertainty(
            plus, meas.tau_list, shots,
        )
        if "minus" in details.populations:
            minus = self._sample_population(details.populations["minus"], role)
            fn, sigma_fn = _ramsey_peak_with_uncertainty(
                minus, meas.tau_list, shots,
            )
            fa = float(details.estimator["f_artificial_ghz"])
            detuning = (fn**2 - fp**2) / (4 * fa) if fp is not None and fn is not None else 0.0
            sigma_detuning = (
                np.hypot(fn * sigma_fn, fp * sigma_fp) / (2 * abs(fa))
                if fp is not None and fn is not None else max(sigma_fp, sigma_fn)
            )
        else:
            sigma_fn = None
            fa = float(details.estimator["f_artificial_ghz"])
            detuning = fa - fp if fp is not None else 0.0
            sigma_detuning = sigma_fp
        frequency = drive + 2 * np.pi * detuning
        return MeasurementSucceeded(
            command_id=cmd.command_id, frequency=frequency,
            uncertainty=2 * np.pi * sigma_detuning,
            method="ramsey_finite_shot", shots=details.circuits * shots,
            applied_bias=bias, applied_drive=drive,
            diagnostics={"circuits": details.circuits, "sampling_role": role,
                         "uncertainty_source": "weighted_ramsey_fit_covariance",
                         "peak_sigma_plus_ghz": sigma_fp,
                         "peak_sigma_minus_ghz": sigma_fn},
        )

    def _execute_track(self, cmd):
        if self.track_method == "ramsey":
            proxy = VerifyFrequency(
                command_id=cmd.command_id,
                frozen_bias=cmd.candidate_bias,
                drive=cmd.predicted_drive,
            )
            event = self._sample_ramsey(proxy, "track")
            event.probe_detuning = event.frequency - cmd.predicted_drive
            event.probe_detuning_uncertainty = event.uncertainty
            return event
        return self._sample_transient(cmd, "track")

    def _execute_monitor(self, cmd):
        return self._sample_transient(cmd, "monitor")

    def _execute_acquire(self, cmd):
        return self._sample_ramsey(cmd, "acquire")

    def _execute_verify(self, cmd):
        return self._sample_ramsey(cmd, "verify")

    def rng_state(self) -> dict:
        return {name: rng.bit_generator.state for name, rng in self._rngs.items()}

    def restore_rng_state(self, state: dict):
        for name, value in state.items():
            if name in self._rngs:
                self._rngs[name].bit_generator.state = value

    def provenance(self) -> dict:
        return {
            "backend": type(self).__name__, "master_seed": self.master_seed,
            "shots_per_circuit": self.shots_per_circuit,
            "shots_by_role": {
                role: self._shots(role)
                for role in ("acquire", "track", "verify", "monitor")
            },
            "assignment_p01": self.assignment_p01,
            "assignment_p10": self.assignment_p10,
            "rng_roles": sorted(self._rngs),
        }


@dataclass
class ProcessSimulationExecutor:
    """Wrapper for simulated drift, controlled jumps, and reference loss."""

    inner: SQCExecutor
    drift_per_command: float = 0.0
    jump_schedule: dict[int, float] = field(default_factory=dict)
    reference_loss_commands: set[int] = field(default_factory=set)
    _command_index: int = 0
    _process_offset: float = 0.0

    def estimate_cost(self, cmd):
        return self.inner.estimate_cost(cmd)

    def set_machine_state(self, state):
        setter = getattr(self.inner, "set_machine_state", None)
        if callable(setter):
            setter(state)

    def execute(self, cmd):
        index = self._command_index
        self._command_index += 1
        if index in self.reference_loss_commands:
            return MeasurementRejected(
                command_id=cmd.command_id,
                reason_code=ReasonCode.CONFIDENCE_INSUFFICIENT,
                diagnostics={"process_event": "reference_loss", "command_index": index},
            )
        self._process_offset += self.drift_per_command
        self._process_offset += float(self.jump_schedule.get(index, 0.0))
        event = self.inner.execute(cmd)
        if isinstance(event, MeasurementSucceeded):
            event.frequency += self._process_offset
            if event.probe_detuning is not None:
                event.probe_detuning += self._process_offset
            event.diagnostics["process_offset"] = self._process_offset
            event.diagnostics["process_event"] = (
                "controlled_jump" if index in self.jump_schedule else "slow_drift"
            )
        return event

    def provenance(self) -> dict:
        return {
            "backend": type(self).__name__,
            "inner_backend": type(self.inner).__name__,
            "drift_per_command": self.drift_per_command,
            "jump_schedule": dict(self.jump_schedule),
            "reference_loss_commands": sorted(self.reference_loss_commands),
        }
