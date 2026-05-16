"""sqc.calibration.frequency — Qubit frequency calibration.

Two major categories:
  - FluxResponseCalibration: f(Φ) curve via flux scan
  - SinglePointFrequencyCalibration: f01 at a single flux working point

See _sensing theory.md §闭环反馈控制 for closed-loop algorithm reference
(Vepsalainen 2022).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.calibration.base import Calibration, CalibrationTable
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse


# ---------------------------------------------------------------------------
# Shared helper: Ramsey FFT frequency fitting
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Internal helpers for Ramsey FFT fitting
# ---------------------------------------------------------------------------

def _fft_peak(p_e_vals: np.ndarray, dt: float) -> float | None:
    """Extract dominant oscillation frequency (GHz) via FFT + quadratic sub-bin
    interpolation.  Returns None when contrast is too low to fit reliably.
    """
    C = float(np.ptp(p_e_vals))
    if C < 0.01:
        return None

    from numpy.fft import rfft, rfftfreq

    p_centered = p_e_vals - np.mean(p_e_vals)
    n_fft = 2048
    fft_vals = np.abs(rfft(p_centered, n=n_fft))
    freqs = rfftfreq(n_fft, d=dt)
    peak_idx = int(np.argmax(fft_vals))

    df_bin = freqs[1] - freqs[0]  # GHz (1 / dt_ns)
    if 0 < peak_idx < len(fft_vals) - 1:
        y1, y2, y3 = fft_vals[peak_idx - 1], fft_vals[peak_idx], fft_vals[peak_idx + 1]
        denom = 2 * (y1 + y3 - 2 * y2)
        if abs(denom) > 1e-15:
            delta_bin = (y1 - y3) / denom
            return float((peak_idx + delta_bin) * df_bin)
        return float(freqs[peak_idx])
    return float(freqs[peak_idx])


def _run_ramsey_sweep(
    qubit: object,
    omega_d: float,
    tau_list: np.ndarray,
    t_rabi: np.ndarray,
    t_global: np.ndarray,
    f_art: float,
) -> np.ndarray:
    """Run one Ramsey τ-sweep with artificial detuning *f_art* (GHz).

    Artificial detuning is applied via phase ramping of the second π/2 pulse:
    phase2 = 2π·f_art·τ, equivalent to an extra detuning *f_art* in the
    rotating frame.
    """
    n_levels = qubit.n_levels
    psi_e = basis(n_levels, 1)
    p_e_vals = np.zeros(len(tau_list), dtype=float)

    for i, tau in enumerate(tau_list):
        phase2 = 2.0 * np.pi * f_art * float(tau)
        ctrl = create_ramsey_pulse(
            t_rabi, tau,
            omega_d=omega_d,
            phase1=0.0, phase2=phase2,
            qubit=qubit,
        )
        H = (
            QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1)
            + QobjEvo(ctrl.hamiltonian_on(t_global), tlist=t_global, order=1)
        )
        result = mesolve(
            H, qubit.state, t_global, [],
            e_ops=[psi_e * psi_e.dag()],
        )
        p_e_vals[i] = float(result.expect[0][-1])

    return p_e_vals


# ---------------------------------------------------------------------------
# Shared helper: Ramsey FFT frequency fitting (signed)
# ---------------------------------------------------------------------------

def _fit_ramsey_frequency(
    qubit: object,
    omega_d: float,
    tau_list: np.ndarray,
    t_rabi: np.ndarray,
    t_global: np.ndarray,
    flux: float = 0.0,
    f_artificial: float | None = 0.1,
) -> float:
    """Run Ramsey at a given DC flux, return fitted f01 (signed) via FFT peak.

    Two operating modes, controlled by *f_artificial*:

    **Single-sweep** (``f_artificial`` is a positive float, default 0.1 GHz):
      Applies a known artificial detuning *f_a* large enough to guarantee
      Δ + f_a > 0 (i.e. *f_a* > |Δ|_max).  The FFT-measured frequency is::

          f_meas = |Δ + f_a| = Δ + f_a   →   Δ = f_meas − f_a

      Fast (1×τ-sweep), but the caller must ensure *f_a* exceeds the
      worst-case |f_q − ω_d|.  Default 100 MHz suffices for measurements
      near the sweet spot.

    **Double-sweep** (``f_artificial=None``):
      Two sweeps with ±50 MHz artificial detuning.  Signed detuning
      follows from the identity::

          Δ = (f_p² − f_n²) / (4·f_a)    (f_a = 0.05 GHz)

      Robust for arbitrary |Δ|, at 2× the time cost.

    Parameters
    ----------
    qubit : TransmonQubit
    omega_d : float
        Drive frequency (angular, GHz·2π).
    tau_list : np.ndarray
        Free evolution times (ns).
    t_rabi : np.ndarray
        Rabi pulse time axis (ns).
    t_global : np.ndarray
        Global evolution time axis (ns).
    flux : float
        DC flux offset (Φ₀). Default 0.
    f_artificial : float or None
        Single-sweep artificial detuning (GHz).  None selects double-sweep.

    Returns
    -------
    float
        Fitted qubit frequency (angular, GHz·2π), signed.
    """
    # -- set qubit flux bias ----------------------------------------------
    t_sig = CONFIG.pulse.make_time(0, 300)
    Phi = FluxSignal(
        type=1 if flux != 0.0 else 0,
        t_list=t_sig,
        amplitude=float(flux),
        offset=0.0,
    )
    qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

    dt_val = float(tau_list[1] - tau_list[0])

    # -- single-sweep mode ------------------------------------------------
    if f_artificial is not None:
        p_e = _run_ramsey_sweep(
            qubit, omega_d, tau_list, t_rabi, t_global, f_artificial,
        )
        f_meas = _fft_peak(p_e, dt_val)
        if f_meas is None:
            return float(omega_d)
        detuning_ghz = f_meas - f_artificial  # Δ + f_a > 0 guaranteed by caller
        return float(omega_d + 2.0 * np.pi * detuning_ghz)

    # -- double-sweep mode (f_artificial is None) -------------------------
    _fa = 0.05  # internal f_a for double-sweep (GHz)
    p_plus = _run_ramsey_sweep(
        qubit, omega_d, tau_list, t_rabi, t_global, +_fa,
    )
    p_minus = _run_ramsey_sweep(
        qubit, omega_d, tau_list, t_rabi, t_global, -_fa,
    )

    f_p = _fft_peak(p_plus, dt_val)
    f_n = _fft_peak(p_minus, dt_val)

    if f_p is None and f_n is None:
        return float(omega_d)
    if f_p is None:
        detuning_ghz = -_fa  # |Δ + f_a| = 0  ⇒  Δ = −f_a
    elif f_n is None:
        detuning_ghz = +_fa  # |Δ − f_a| = 0  ⇒  Δ = +f_a
    else:
        detuning_ghz = (f_p ** 2 - f_n ** 2) / (4.0 * _fa)

    return float(omega_d + 2.0 * np.pi * detuning_ghz)


# ===================================================================
# FluxResponseCalibration — f(Φ) curve
# ===================================================================

@dataclass
class FluxResponseCalibration(Calibration):
    """Magnetic flux response calibration — measure f(Φ) lookup table.

    Replaces src/protocal.py:Calibration.calibrate case 1 (ramsey scan).

    Methods
    -------
    - ``"ramsey"``: scan DC flux steps, measure frequency via Ramsey at each
      point → f(Φ) CalibrationTable (kind="f_phi").
    - ``"transient"``: unknown transient signal → polynomial fit of Δω(Φ).
      **Requires Track B 1.2.**

    Parameters
    ----------
    qubit : TransmonQubit
    method : str
        "ramsey" (complete) or "transient" (Track B 1.2 stub).
    h_list : np.ndarray or None
        Flux values to scan (Φ₀). Default linspace(-0.03, 0.03, 51).
    tau : float
        Free precession time per flux point (ns). Default 100.0.
    t_rabi : np.ndarray
        Rabi pulse time axis (ns).
    """

    qubit: object
    method: Literal["ramsey", "transient"] = "ramsey"
    h_list: np.ndarray | None = None
    tau: float = 100.0
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )

    def __post_init__(self):
        if self.h_list is None:
            self.h_list = np.linspace(-0.03, 0.03, 51)

    def calibrate(self) -> CalibrationTable:
        match self.method:
            case "ramsey":
                return self._calibrate_ramsey()
            case "transient":
                return self._calibrate_transient()

    # ------------------------------------------------------------------
    def _calibrate_ramsey(self) -> CalibrationTable:
        """Scan flux Φ, fit Ramsey detuning at each step → f(Φ) table."""
        omega_d = self.qubit.frequency
        tau_list = CONFIG.pulse.make_time(0, 200)
        t_global = CONFIG.pulse.t_global.copy()
        frequency_list: list[float] = []

        for h in self.h_list:
            f01 = _fit_ramsey_frequency(
                self.qubit, omega_d, tau_list,
                self.t_rabi, t_global, flux=float(h),
            )
            frequency_list.append(f01)

        return CalibrationTable(
            name="flux_response_ramsey",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f_phi",
            inputs=np.asarray(self.h_list, dtype=float),
            outputs=np.asarray(frequency_list, dtype=float),
            fit_params={"method": "ramsey", "omega_d": float(omega_d)},
            metadata={},
        )

    # ------------------------------------------------------------------
    def _calibrate_transient(self) -> CalibrationTable:
        raise NotImplementedError(
            "FluxResponseCalibration(method='transient'): "
            "requires Track B 1.2 (transient frequency calibration)."
        )


# ===================================================================
# SinglePointFrequencyCalibration — single-point f01
# ===================================================================

@dataclass
class SinglePointFrequencyCalibration(Calibration):
    """Single-point qubit frequency calibration.

    Measures f01 at a single flux working point. Replaces
    src/protocal.py:Calibration.calibrate case 0 (ramsey).

    Methods
    -------
    - ``"ramsey"``: single Ramsey FFT → f01.
    - ``"closed_loop"``: secant-method iterative convergence to a target
      frequency (Vepsalainen 2022 closed-loop feedback). The per-iteration
      frequency measurement is controlled by ``measure_method``.
    - ``"transient"``: single-shot transient measurement → f01.
      **Requires Track B 1.2.**

    Parameters
    ----------
    qubit : TransmonQubit
    method : str
        "ramsey", "closed_loop", or "transient".
    tau_list : np.ndarray or None
        Free evolution times for Ramsey (ns). Default make_time(0, 200).
    t_rabi : np.ndarray
        Rabi pulse time axis (ns).
    t_global : np.ndarray or None
        Global evolution time axis (ns).

    Closed-loop parameters
    ----------------------
    f_target : float or None
        Target qubit frequency (angular, GHz·2π). Required for closed_loop.
    epsilon_f : float
        Convergence tolerance (GHz·2π). Default 1e-4.
    V_a : float or None
        Lower voltage bound bracketing target. Required for closed_loop.
    V_b : float or None
        Upper voltage bound bracketing target. Required for closed_loop.
    max_iter : int
        Maximum secant iterations. Default 20.
    measure_method : str
        Per-iteration frequency measurement technique:
        - ``"ramsey"``: Ramsey FFT (default).
        - ``"transient"``: transient measurement (Track B 1.2).
    bracket_tightening : bool
        If True (default), tighten the bracket [V_a, V_b] on each iteration
        (regula falsi) for faster convergence.  Set False to observe pure
        secant iteration trends at the cost of more iterations.
    step_method : str
        Root-finding algorithm for closed-loop iteration:
        - ``"secant"`` (default): secant method with bracket tightening.
          Fast convergence (1–3 iterations) when the measurement is smooth.
        - ``"bisection"``: simple bisection.  Converges in O(log₂(range/ε))
          iterations (~10–15), giving a clear exponential convergence trend
          for visualisation and diagnostics.
    """

    qubit: object
    method: Literal["ramsey", "closed_loop", "transient"] = "ramsey"

    # Ramsey params
    tau_list: np.ndarray | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    t_global: np.ndarray | None = None

    # Closed-loop params
    f_target: float | None = None
    epsilon_f: float = 1e-4
    V_a: float | None = None
    V_b: float | None = None
    max_iter: int = 20
    measure_method: Literal["ramsey", "transient"] = "ramsey"
    bracket_tightening: bool = True
    step_method: Literal["secant", "bisection"] = "secant"

    def __post_init__(self):
        if self.tau_list is None:
            self.tau_list = CONFIG.pulse.make_time(0, 200)
        if self.t_global is None:
            self.t_global = CONFIG.pulse.t_global.copy()

    # ------------------------------------------------------------------
    def calibrate(self) -> CalibrationTable:
        match self.method:
            case "ramsey":
                return self._calibrate_ramsey()
            case "closed_loop":
                return self._calibrate_closed_loop()
            case "transient":
                return self._calibrate_transient()

    # ------------------------------------------------------------------
    # Single-shot Ramsey
    # ------------------------------------------------------------------

    def _calibrate_ramsey(self) -> CalibrationTable:
        """Single Ramsey FFT → f01 at zero flux."""
        omega_d = self.qubit.frequency
        fitted_freq = _fit_ramsey_frequency(
            self.qubit, omega_d, self.tau_list,
            self.t_rabi, self.t_global, flux=0.0,
        )
        return CalibrationTable(
            name="frequency_ramsey",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f01",
            inputs=np.array([0.0]),
            outputs=np.array([fitted_freq]),
            fit_params={
                "method": "ramsey",
                "omega_d": float(omega_d),
            },
            metadata={},
        )

    # ------------------------------------------------------------------
    # Closed-loop feedback (Vepsalainen 2022)
    # ------------------------------------------------------------------

    def _calibrate_closed_loop(self) -> CalibrationTable:
        """Closed-loop frequency calibration.

        Iteratively adjusts flux bias voltage to drive qubit frequency to
        f_target.  Uses the secant method (default) or bisection for
        root-finding on r(V) = f_Q(V) - f_target = 0.

        Requires f_target, V_a, V_b to be set.  V_a and V_b must bracket
        the target: r(V_a) · r(V_b) < 0.
        """
        if self.f_target is None:
            raise ValueError("f_target is required for closed_loop method.")
        if self.V_a is None or self.V_b is None:
            raise ValueError(
                "V_a and V_b are required for closed_loop method. "
                "Run FluxResponseCalibration first to bracket the target."
            )

        f_target = self.f_target
        V_lo, V_hi = min(self.V_a, self.V_b), max(self.V_a, self.V_b)
        epsilon = self.epsilon_f

        if self.step_method == "bisection":
            return self._closed_loop_bisection(f_target, V_lo, V_hi, epsilon)
        else:
            return self._closed_loop_secant(f_target, V_lo, V_hi, epsilon)

    # ------------------------------------------------------------------
    def _closed_loop_secant(
        self, f_target: float, V_lo: float, V_hi: float, epsilon: float,
    ) -> CalibrationTable:
        """Secant-method closed-loop iteration."""

        # Initial midpoint
        V_n = (V_lo + V_hi) / 2.0
        r_n = self._measure_frequency(V_n) - f_target
        n_iter = 0

        V_prev = V_lo
        r_prev = self._measure_frequency(V_lo) - f_target

        history: list[dict] = []

        while abs(r_n) > epsilon and n_iter < self.max_iter:
            denom = r_n - r_prev
            if abs(denom) < 1e-15:
                V_next = (V_lo + V_hi) / 2.0
            else:
                V_next = V_n - r_n * (V_n - V_prev) / denom

            # Bound check
            if V_next < V_lo or V_next > V_hi:
                V_n = (V_lo + V_hi) / 2.0
            else:
                V_prev, r_prev = V_n, r_n
                V_n = V_next

            f_n = self._measure_frequency(V_n)
            r_n = f_n - f_target
            n_iter += 1

            history.append({
                "iter": n_iter, "V": float(V_n),
                "f": float(f_n), "residual": float(r_n),
            })

            # Bracket tightening (regula falsi)
            if self.bracket_tightening:
                if r_n > 0:
                    V_hi = V_n
                else:
                    V_lo = V_n

        return self._build_result(f_target, V_n, r_n, n_iter, epsilon, history)

    # ------------------------------------------------------------------
    def _closed_loop_bisection(
        self, f_target: float, V_lo: float, V_hi: float, epsilon: float,
    ) -> CalibrationTable:
        """Bisection-method closed-loop iteration.

        Each iteration halves the bracket width, giving a clean
        exponential convergence trend — ideal for visualisation.

        Handles even f(Φ) by auto-splitting the bracket at Φ=0 when
        f(V_lo) and f(V_hi) have the same sign (the root lies on one
        side of the sweet spot).
        """
        r_lo = self._measure_frequency(V_lo) - f_target
        r_hi = self._measure_frequency(V_hi) - f_target

        # Auto-split bracket if both ends have the same sign (even f(Φ))
        if r_lo * r_hi > 0:
            V_mid0 = (V_lo + V_hi) / 2.0
            r_mid0 = self._measure_frequency(V_mid0) - f_target
            if r_lo * r_mid0 < 0:
                V_hi, r_hi = V_mid0, r_mid0
            elif r_mid0 * r_hi < 0:
                V_lo, r_lo = V_mid0, r_mid0
            else:
                # f is flat or all same sign — can't bracket
                raise RuntimeError(
                    f"Bisection cannot bracket root: r(lo)={r_lo/(2*np.pi)*1e3:.1f}, "
                    f"r(mid)={r_mid0/(2*np.pi)*1e3:.1f}, "
                    f"r(hi)={r_hi/(2*np.pi)*1e3:.1f} MHz"
                )

        n_iter = 0
        history: list[dict] = []

        while (V_hi - V_lo) / 2.0 > 1e-15 and n_iter < self.max_iter:
            V_mid = (V_lo + V_hi) / 2.0
            r_mid = self._measure_frequency(V_mid) - f_target
            n_iter += 1

            history.append({
                "iter": n_iter, "V": float(V_mid),
                "f": float(r_mid + f_target),
                "residual": float(r_mid),
                "bracket_width": float(V_hi - V_lo),
            })

            if abs(r_mid) <= epsilon:
                return self._build_result(f_target, V_mid, r_mid, n_iter, epsilon, history)

            if r_mid * r_lo < 0:
                V_hi, r_hi = V_mid, r_mid
            else:
                V_lo, r_lo = V_mid, r_mid

        V_mid = (V_lo + V_hi) / 2.0
        r_mid = self._measure_frequency(V_mid) - f_target
        return self._build_result(f_target, V_mid, r_mid, n_iter, epsilon, history)

    # ------------------------------------------------------------------
    def _build_result(
        self, f_target: float, V_opt: float, residual: float,
        n_iter: int, epsilon: float, history: list[dict],
    ) -> CalibrationTable:
        """Package closed-loop result into a CalibrationTable."""
        return CalibrationTable(
            name="frequency_closed_loop",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f01",
            inputs=np.array([V_opt]),
            outputs=np.array([f_target + residual]),
            fit_params={
                "method": "closed_loop",
                "step_method": self.step_method,
                "measure_method": self.measure_method,
                "f_target": float(f_target),
                "V_opt": float(V_opt),
                "n_iter": n_iter,
                "residual": float(residual),
                "converged": abs(residual) <= epsilon,
                "history": history,
            },
            metadata={},
        )

    def _measure_frequency(self, flux: float) -> float:
        """Measure f01 at a given flux offset.

        Dispatches to Ramsey or transient measurement based on
        self.measure_method.
        """
        if self.measure_method == "ramsey":
            omega_d = self.qubit.frequency
            return _fit_ramsey_frequency(
                self.qubit, omega_d, self.tau_list,
                self.t_rabi, self.t_global, flux=flux,
                f_artificial=None,  # double-sweep: flux may be far from sweet spot
            )
        elif self.measure_method == "transient":
            return self._measure_frequency_transient(flux)
        else:
            raise ValueError(f"Unknown measure_method: {self.measure_method}")

    def _measure_frequency_transient(self, flux: float) -> float:
        """Transient-based single-point frequency measurement.

        Uses orthogonal Ramsey readout (R_y–R_x and R_y–R_{-x}) with tau=0
        to measure detuning via differential p_e and the control-pulse
        kernel sensitivity G_α = ∫ k(t) dt.

        Theory ref: _sensing theory.md §瞬态磁场协议与核函数策略
        """
        omega_d = self.qubit.frequency
        n_levels = self.qubit.n_levels

        # -- set qubit to target flux ------------------------------------
        t_sig = CONFIG.pulse.make_time(0, 300)
        Phi = FluxSignal(
            type=1 if flux != 0.0 else 0,
            t_list=t_sig,
            amplitude=float(flux),
            offset=0.0,
        )
        self.qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

        # -- build orthogonal readout pulses (tau=0) ---------------------
        tau = 0.0
        ctrl_x = create_ramsey_pulse(
            self.t_rabi, tau,
            omega_d=omega_d,
            phase1=np.pi / 2, phase2=0.0,
            qubit=self.qubit,
        )
        ctrl_mx = create_ramsey_pulse(
            self.t_rabi, tau,
            omega_d=omega_d,
            phase1=np.pi / 2, phase2=np.pi,
            qubit=self.qubit,
        )
        # -- run both measurements (P7: use hamiltonian_on on t_global) ---
        psi_e = basis(n_levels, 1)
        H_base = QobjEvo(
            self.qubit.H_list, tlist=self.qubit.mag_signal.t_list, order=1,
        )

        H_x = H_base + QobjEvo(
            ctrl_x.hamiltonian_on(self.t_global), tlist=self.t_global, order=1,
        )
        res_x = mesolve(
            H_x, self.qubit.state, self.t_global, [],
            e_ops=[psi_e * psi_e.dag()],
        )
        p_x = float(res_x.expect[0][-1])

        H_mx = H_base + QobjEvo(
            ctrl_mx.hamiltonian_on(self.t_global), tlist=self.t_global, order=1,
        )
        res_mx = mesolve(
            H_mx, self.qubit.state, self.t_global, [],
            e_ops=[psi_e * psi_e.dag()],
        )
        p_mx = float(res_mx.expect[0][-1])

        p_diff = (p_x - p_mx) / 2.0

        # -- kernel sensitivity G_α = ∫ k(t) dt --------------------------
        ctrl_x.get_kernel(self.qubit)
        kernel = np.asarray(ctrl_x.kernel, dtype=float)
        t_kernel = np.asarray(ctrl_x.t_samples, dtype=float)
        G_alpha = float(np.trapezoid(kernel, t_kernel))

        # restore qubit state after get_kernel side effects
        self.qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

        if abs(G_alpha) < 1e-15:
            return float(omega_d)

        delta_omega = p_diff / G_alpha
        return float(omega_d + delta_omega)

    # ------------------------------------------------------------------
    # Transient (stub)
    # ------------------------------------------------------------------

    def _calibrate_transient(self) -> CalibrationTable:
        raise NotImplementedError(
            "SinglePointFrequencyCalibration(method='transient'): "
            "requires Track B 1.2."
        )
