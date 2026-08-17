"""sqc.calibration.frequency_control — shared frequency-tracking controller.

Provides the stepwise damped-secant tracker that underpins both the legacy
:meth:`SinglePointFrequencyCalibration._closed_loop_gradient` and the new
event-driven :class:`FrequencyStateMachine` (Track state).

The controller is **pure**: no I/O, no qubit model, no mesolve.  The caller
is responsible for resolving drive frequencies, performing measurements, and
enriching diagnostic state (linear-range guards, cost accounting).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class FrequencyEstimate:
    """Single-point frequency measurement result with metadata.

    Parameters
    ----------
    frequency : float
        Measured qubit frequency (angular, rad·GHz), signed.
    uncertainty : float
        Estimated 1-σ uncertainty (rad·GHz).  0.0 for deterministic simulations.
    valid : bool
        Whether the measurement passed technical quality checks.
    ambiguous : bool
        Whether the measurement could not resolve sign / branch / aliasing.
    method : str
        Measurement method (``"ramsey"``, ``"transient"``).
    shots : int
        Number of repeated measurements (or solver calls).
    elapsed_time : float
        Wall-clock time for this measurement (s).
    diagnostics : dict
        Raw residuals, fit quality, kernel coefficients, etc.
    """

    frequency: float
    uncertainty: float = 0.0
    valid: bool = True
    ambiguous: bool = False
    method: str = ""
    shots: int = 0
    elapsed_time: float = 0.0
    diagnostics: dict | None = None

    def __post_init__(self):
        if self.diagnostics is None:
            self.diagnostics = {}


@dataclass
class TrackSnapshot:
    """Full tracker state at a point in time (immutable record).

    Carries the current and previous measurement, residual, best point,
    convergence streak, and the most recent sensitivity estimate.  The
    caller can serialise this for checkpoint / replay.
    """

    V: float  # current bias (Φ₀)
    V_prev: float  # previous bias (Φ₀)
    f: float | None  # current measured frequency (rad·GHz)
    f_prev: float | None  # previous measured frequency
    e: float | None  # residual = f − f_target (rad·GHz)
    e_prev: float | None  # previous residual
    f_target: float  # target frequency (rad·GHz)
    best_V: float  # bias of the best point seen so far
    best_abs_e: float  # |residual| at the best point
    streak: int  # consecutive in-tolerance iterations
    n_iter: int  # number of accepted measurements (≥ 1)
    s_hat: float | None = None  # local sensitivity ∂f/∂V (secant est.)
    uncertainty: float = 0.0  # current 1-sigma frequency uncertainty
    uncertainty_prev: float | None = None  # previous 1-sigma uncertainty


@dataclass
class TrackProposal:
    """Proposed next bias point from the controller.

    The caller **must** resolve the drive frequency (``omega_d``) for this
    point, perform the measurement, and pass the result to :meth:`accept`.
    """

    V_next: float  # proposed bias for next measurement (Φ₀)
    step: float  # bias step that will be applied (ΔΦ₀)
    s_hat: float | None  # secant sensitivity estimate at this step
    converged: bool  # True → no further steps needed
    diagnostics: dict = field(default_factory=dict)


@dataclass
class TrackStepResult:
    """Result returned by :meth:`DampedSecantTracker.accept`."""

    snapshot: TrackSnapshot  # updated tracker state
    is_best: bool  # this measurement is the best so far
    stopped_by: str  # "tolerance" | "max_iter" | "predicate" | "running"


# ---------------------------------------------------------------------------
# Controller
# ---------------------------------------------------------------------------


class DampedSecantTracker:
    """Stepwise damped-secant controller for single-point frequency tracking.

    Pure controller — no I/O, no qubit model.  Usage::

        tracker = DampedSecantTracker(f_target=..., ...)
        snapshot = tracker.initialize(seed_estimate, V_seed)

        while True:
            proposal = tracker.propose(snapshot)
            if proposal.converged:
                break
            # --- caller responsibilities ---
            omega_d = resolve_drive(proposal)
            estimate = measure(proposal.V_next, omega_d)
            # --------------------------------
            result = tracker.accept(snapshot, estimate, proposal)
            snapshot = result.snapshot
            if result.stopped_by != "running":
                break

    The same controller is used by:
    - ``SinglePointFrequencyCalibration._closed_loop_gradient()`` (legacy batch)
    - ``FrequencyStateMachine`` Track state (stepwise, one propose/accept
      per event loop iteration)

    Parameters
    ----------
    f_target : float
        Target qubit frequency (angular, rad·GHz).
    damping : float
        Damping factor ∈ (0, 1] for the secant step.  Default 0.8.
    first_bias_step : float
        Probe step size (Φ₀) on the first iteration.  Default 0.01.
    expected_sensitivity_sign : {-1, 1} or None
        Known sign of the local sensitivity ``∂f/∂V``.  When configured,
        blind/fallback steps follow this direction and secant estimates with
        the opposite sign are rejected before another measurement is issued.
        ``None`` preserves the legacy direction convention.
    max_bias_step : float
        Per-step clamp (Φ₀).  Default 0.02.
    V_lo, V_hi : float or None
        Optional hard flux bounds.  None disables clamping.
    converge_streak : int
        Number of consecutive in-tolerance iterations before declaring
        convergence.  Default 1 (stop on first hit).
    max_iter : int
        Hard iteration cap.  Default 20.
    epsilon_f : float
        Convergence tolerance on |residual| (rad·GHz).  Default 1e-4.
    stop_predicate : callable or None
        Optional early-stop hook called after each accept.  Receives a
        state dict; return True to stop with ``stopped_by="predicate"``.
    """

    def __init__(
        self,
        f_target: float,
        damping: float = 0.8,
        first_bias_step: float = 0.01,
        max_bias_step: float = 0.02,
        V_lo: float | None = None,
        V_hi: float | None = None,
        converge_streak: int = 1,
        max_iter: int = 20,
        epsilon_f: float = 1e-4,
        stop_predicate: Callable[[dict], bool] | None = None,
        expected_sensitivity_sign: int | None = None,
    ):
        self.f_target = float(f_target)
        self.damping = float(damping)
        self.first_bias_step = float(first_bias_step)
        if expected_sensitivity_sign not in (None, -1, 1):
            raise ValueError("expected_sensitivity_sign must be -1, 1, or None")
        self.expected_sensitivity_sign = expected_sensitivity_sign
        self.max_bias_step = float(max_bias_step)
        self.V_lo = float(V_lo) if V_lo is not None else None
        self.V_hi = float(V_hi) if V_hi is not None else None
        if self.V_lo is not None and self.V_hi is not None and self.V_lo > self.V_hi:
            raise ValueError("V_lo must not exceed V_hi")
        self.converge_streak = max(1, int(converge_streak))
        self.max_iter = int(max_iter)
        self.epsilon_f = float(epsilon_f)
        self.stop_predicate = stop_predicate

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def initialize(
        self,
        seed: FrequencyEstimate,
        V_seed: float,
    ) -> TrackSnapshot:
        """Create the initial snapshot from the first measurement.

        The caller measures at **V_seed** with the appropriate drive
        frequency, wraps the result as a :class:`FrequencyEstimate`, and
        passes it here.
        """
        e = seed.frequency - self.f_target
        abs_e = abs(e)
        streak = 1 if abs_e <= self.epsilon_f else 0
        return TrackSnapshot(
            V=float(V_seed),
            V_prev=float(V_seed),
            f=seed.frequency,
            f_prev=None,
            e=e,
            e_prev=None,
            f_target=self.f_target,
            best_V=float(V_seed),
            best_abs_e=abs_e,
            streak=streak,
            n_iter=0,
            s_hat=None,
            uncertainty=seed.uncertainty,
            uncertainty_prev=None,
        )

    def propose(self, snapshot: TrackSnapshot) -> TrackProposal:
        """Compute the next bias step from the current tracker state.

        Returns a :class:`TrackProposal` with ``converged=True`` when the
        convergence streak or iteration budget has been met — the caller
        should stop the loop without measuring.
        """
        n_streak = self.converge_streak

        # --- already converged? -----------------------------------------
        if snapshot.streak >= n_streak:
            return TrackProposal(
                V_next=snapshot.V,
                step=0.0,
                s_hat=snapshot.s_hat,
                converged=True,
                diagnostics={"reason": "already_converged"},
            )

        # --- budget exhausted? ------------------------------------------
        if snapshot.n_iter >= self.max_iter:
            return TrackProposal(
                V_next=snapshot.V,
                step=0.0,
                s_hat=snapshot.s_hat,
                converged=True,
                diagnostics={"reason": "max_iter"},
            )

        e = snapshot.e
        s_secant: float | None = None

        def probe_step() -> float:
            residual_sign = np.sign(e) if e is not None else 1.0
            sensitivity_sign = self.expected_sensitivity_sign or 1
            return float(self.first_bias_step * sensitivity_sign * residual_sign)

        # --- compute step -----------------------------------------------
        if e is not None and abs(e) <= self.epsilon_f:
            # In tolerance: hold position (only reached when
            # converge_streak > 1 and re-measuring to confirm).
            step = 0.0
        elif snapshot.n_iter == 0:
            # No gradient yet → fixed probe step.
            step = probe_step()
        else:
            de = (
                e - snapshot.e_prev
                if (e is not None and snapshot.e_prev is not None)
                else 0.0
            )
            dV = snapshot.V - snapshot.V_prev
            if abs(de) > 1e-15:
                grad_inv = dV / de  # dB/de
                step = self.damping * e * grad_inv  # damped secant
                if abs(dV) > 1e-15:
                    s_secant = de / dV  # local ∂f/∂V
            else:
                # Gradient undefined → reuse probe step.
                step = probe_step()

        direction_valid = (
            s_secant is None
            or self.expected_sensitivity_sign is None
            or np.sign(s_secant) == self.expected_sensitivity_sign
        )
        if not direction_valid:
            return TrackProposal(
                V_next=snapshot.V,
                step=0.0,
                s_hat=s_secant,
                converged=False,
                diagnostics={
                    "reason": "sensitivity_direction_mismatch",
                    "expected_sensitivity_sign": self.expected_sensitivity_sign,
                },
            )

        # --- clamp step -------------------------------------------------
        step = float(np.clip(step, -self.max_bias_step, self.max_bias_step))

        # --- apply step → V_next ----------------------------------------
        V_next = snapshot.V - step

        # --- hard bounds ------------------------------------------------
        if self.V_lo is not None:
            V_next = max(self.V_lo, V_next)
        if self.V_hi is not None:
            V_next = min(self.V_hi, V_next)

        return TrackProposal(
            V_next=V_next,
            step=step,
            s_hat=s_secant,
            converged=False,
            diagnostics={},
        )

    def accept(
        self,
        snapshot: TrackSnapshot,
        estimate: FrequencyEstimate,
        proposal: TrackProposal,
        *,
        extra_predicate_state: dict | None = None,
    ) -> TrackStepResult:
        """Accept a measurement taken at ``proposal.V_next``.

        Parameters
        ----------
        snapshot : TrackSnapshot
            The tracker state *before* the measurement.
        estimate : FrequencyEstimate
            The measurement result at ``proposal.V_next``.
        proposal : TrackProposal
            The proposal that prompted this measurement.
        extra_predicate_state : dict or None
            Optional fields merged into the predicate state dict (e.g.
            ``delta``, ``out_of_range``, ``history`` set by the caller).

        Returns
        -------
        TrackStepResult
        """
        V_next = proposal.V_next
        f_next = estimate.frequency
        e_next = f_next - self.f_target
        abs_e_next = abs(e_next)
        n_iter = snapshot.n_iter + 1

        # --- convergence de-bounce --------------------------------------
        streak = snapshot.streak + 1 if abs_e_next <= self.epsilon_f else 0

        # --- best-point tracking ----------------------------------------
        is_best = abs_e_next < snapshot.best_abs_e
        best_V = V_next if is_best else snapshot.best_V
        best_abs_e = abs_e_next if is_best else snapshot.best_abs_e

        new_snapshot = TrackSnapshot(
            V=V_next,
            V_prev=snapshot.V,
            f=f_next,
            f_prev=snapshot.f,
            e=e_next,
            e_prev=snapshot.e,
            f_target=self.f_target,
            best_V=best_V,
            best_abs_e=best_abs_e,
            streak=streak,
            n_iter=n_iter,
            s_hat=proposal.s_hat,
            uncertainty=estimate.uncertainty,
            uncertainty_prev=snapshot.uncertainty,
        )

        # --- determine stop reason --------------------------------------
        n_streak = self.converge_streak
        if streak >= n_streak:
            stopped_by = "tolerance"
        elif n_iter >= self.max_iter:
            stopped_by = "max_iter"
        elif self.stop_predicate is not None:
            state = self._build_predicate_state(
                new_snapshot,
                proposal,
                is_best,
                extra_predicate_state,
            )
            stopped_by = "predicate" if self.stop_predicate(state) else "running"
        else:
            stopped_by = "running"

        return TrackStepResult(
            snapshot=new_snapshot,
            is_best=is_best,
            stopped_by=stopped_by,
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _build_predicate_state(
        self,
        snapshot: TrackSnapshot,
        proposal: TrackProposal,
        is_best: bool,
        extra: dict | None,
    ) -> dict:
        """Build the state dict passed to ``stop_predicate``."""
        e_val = snapshot.e
        e_prev_val = snapshot.e_prev
        state: dict = {
            "iter": snapshot.n_iter,
            "V": float(snapshot.V),
            "residual": float(e_val) if e_val is not None else float("inf"),
            "abs_residual": float(abs(e_val)) if e_val is not None else float("inf"),
            "step": float(proposal.step),
            "de": float(e_val - e_prev_val)
            if (e_val is not None and e_prev_val is not None)
            else 0.0,
            "dV": float(snapshot.V - snapshot.V_prev),
            "delta": 0.0,  # caller overrides via extra
            "out_of_range": False,  # caller overrides via extra
            "best_abs_residual": float(snapshot.best_abs_e),
        }
        if extra:
            state.update(extra)
        return state

    @staticmethod
    def resolve_track_drive(
        V: float,
        V_prev: float,
        f_prev: float | None,
        s_hat: float | None,
    ) -> float:
        """Predict f_q at bias *V* using the last measurement and sensitivity.

        ``omega_d = f_prev + s_hat·(V − V_prev)`` — the standard "track"
        drive policy that keeps the measurement in its linear window.
        """
        if f_prev is None:
            return 0.0
        if s_hat is None:
            return float(f_prev)
        return float(f_prev + s_hat * (V - V_prev))
