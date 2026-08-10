"""sqc.calibration.frequency — Qubit frequency calibration.
Three calibration classes:
  - FluxResponseCalibration: f(Φ) curve via flux scan
  - FrequencyMeasurement: single-point f01 measurement (ramsey or transient)
  - SinglePointFrequencyCalibration: single-point frequency *tuning*
        (currently closed-loop feedback; extensible for future methods)

See _sensing theory.md §闭环反馈控制 for closed-loop algorithm reference
(Vepsalainen 2022).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Literal, Optional

import numpy as np
from qutip import Qobj, QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.calibration.base import Calibration, CalibrationTable
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse


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
        # max_step prevents adaptive stepper from skipping over narrow
        # pulse windows on the long t_global axis.
        result = mesolve(
            H, qubit.state, t_global, [],
            e_ops=[psi_e * psi_e.dag()],
            options={"max_step": float(CONFIG.awg.dt)},
        )
        p_e_vals[i] = float(result.expect[0][-1])

    return p_e_vals


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
      f_a − Δ > 0 (i.e. *f_a* > |Δ|_max).  Due to the n+I/2 convention
      used by ``qubit_in_mag``, the FFT-measured frequency is::

          f_meas = |f_a − Δ| = f_a − Δ   →   Δ = f_a − f_meas

      Fast (1×τ-sweep), but the caller must ensure *f_a* exceeds the
      worst-case |f_q − ω_d|.  Default 100 MHz suffices for measurements
      near the sweet spot.

    **Double-sweep** (``f_artificial=None``):
      Two sweeps with ±50 MHz artificial detuning.  Signed detuning
      follows from the identity::

          Δ = (f_n² − f_p²) / (4·f_a)    (f_a = 0.05 GHz)

      where f_p = |f_a − Δ| and f_n = |−f_a − Δ|.  Robust for arbitrary
      |Δ|, at 2× the time cost.

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
        detuning_ghz = f_artificial - f_meas  # f_meas = |f_a - Δ| (n+I/2 convention)
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
        detuning_ghz = +_fa  # |f_a - Δ| = 0  ⇒  Δ = +f_a
    elif f_n is None:
        detuning_ghz = -_fa  # |-f_a - Δ| = 0  ⇒  Δ = -f_a
    else:
        detuning_ghz = (f_n ** 2 - f_p ** 2) / (4.0 * _fa)

    return float(omega_d + 2.0 * np.pi * detuning_ghz)


# ---------------------------------------------------------------------------
# Module-level G₃ Taylor coefficient cache
# ---------------------------------------------------------------------------
# Key: (hash of t_rabi.tobytes(), round(omega_d, 6))
# Value: (G1_fit, G3_taylor)
# The G₁, G₃ coefficients depend only on the pulse sequence shape (t_rabi
# duration + dt) and drive frequency — not on the qubit state.  Caching
# avoids re-running the ~42-mesolve calibration scan on every measurement.
_g3_cache: dict[tuple, tuple[float, float]] = {}


def _make_cache_key(
    t_rabi: np.ndarray, omega_d: float, tag: object = "adaptive",
) -> tuple:
    """Return a hashable cache key from the pulse time axis and drive frequency.

    ``tag`` distinguishes calibration variants (e.g. a manual ``delta_max``
    override vs the adaptive default) so they do not collide in the cache.
    """
    return (hash(t_rabi.tobytes()), round(float(omega_d), 6), tag)


# ---------------------------------------------------------------------------
# Internal helper: G₃ Taylor calibration via p_diff(Δ) polynomial fit
# ---------------------------------------------------------------------------

def _fit_g3_at_delta_max(
    qubit: object,
    omega_d: float,
    t_rabi: np.ndarray,
    t_global: np.ndarray,
    n_scan: int,
    delta_max_ghz: float,
) -> tuple[float, float]:
    """Single odd-polynomial fit of ``p_diff(Δ)`` over ``±delta_max_ghz``.

    Adds a constant ``-(Δ/2)·σ_z`` detuning at the sweet spot (known Δ of
    both signs), measures ``p_diff`` via orthogonal Ramsey readout, and
    fits ``p_diff = c1·Δ + c3·Δ³ + c5·Δ⁵`` (in the **detuning** convention,
    Δ = ω_q − ω_d).  Returns ``(G1_fit, G3_taylor) = (c1, 6·c3)``.

    No caching — the caller (:func:`_calibrate_g3_taylor`) handles that.
    Cost: ``2·(n_scan−1)`` mesolve calls.
    """
    n_levels = qubit.n_levels
    psi_e = basis(n_levels, 1)
    t_sig = CONFIG.pulse.make_time(0, 300)
    sigma_z = Qobj(np.diag([1.0, -1.0] + [0.0] * (n_levels - 2)))

    Phi_zero = FluxSignal(type=0, t_list=t_sig)
    qubit.qubit_in_mag(Phi_zero, frame=1, omega_d=omega_d)

    # Stretched grid: finer spacing near Δ=0 for derivative accuracy.
    x = np.linspace(-1.0, 1.0, n_scan)
    delta_scan_ghz = delta_max_ghz * np.sign(x) * (np.abs(x) ** 1.5)
    delta_scan_ghz = delta_scan_ghz[delta_scan_ghz != 0.0]

    delta_vals: list[float] = []
    p_diff_vals: list[float] = []

    for delta_ghz in delta_scan_ghz:
        delta_rad = 2.0 * np.pi * float(delta_ghz)
        delta_vals.append(delta_rad)

        detuning_coeff = -0.5 * delta_rad * np.ones_like(t_global)
        H_detuning = QobjEvo([[sigma_z, detuning_coeff]], tlist=t_global, order=1)

        tau = 0.0
        ctrl_x = create_ramsey_pulse(
            t_rabi, tau, omega_d=omega_d,
            phase1=np.pi / 2, phase2=0.0, qubit=qubit,
        )
        ctrl_mx = create_ramsey_pulse(
            t_rabi, tau, omega_d=omega_d,
            phase1=np.pi / 2, phase2=np.pi, qubit=qubit,
        )

        H_base = QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1)
        H_pulse_x = QobjEvo(ctrl_x.hamiltonian_on(t_global), tlist=t_global, order=1)
        H_pulse_mx = QobjEvo(ctrl_mx.hamiltonian_on(t_global), tlist=t_global, order=1)

        res_x = mesolve(
            H_base + H_detuning + H_pulse_x, qubit.state, t_global, [],
            e_ops=[psi_e * psi_e.dag()],
            options={"max_step": float(CONFIG.awg.dt)},
        )
        res_mx = mesolve(
            H_base + H_detuning + H_pulse_mx, qubit.state, t_global, [],
            e_ops=[psi_e * psi_e.dag()],
            options={"max_step": float(CONFIG.awg.dt)},
        )
        p_diff_vals.append((float(res_x.expect[0][-1]) - float(res_mx.expect[0][-1])) / 2.0)

    delta_arr = np.asarray(delta_vals, dtype=float)
    p_diff_arr = np.asarray(p_diff_vals, dtype=float)
    A = np.column_stack([delta_arr, delta_arr ** 3, delta_arr ** 5])
    coeffs, *_ = np.linalg.lstsq(A, p_diff_arr, rcond=None)
    return float(coeffs[0]), float(6.0 * coeffs[1])


def _calibrate_g3_taylor(
    qubit: object,
    omega_d: float,
    t_rabi: np.ndarray,
    t_global: np.ndarray,
    n_scan: int = 21,
    delta_max_ghz: float | None = None,
) -> tuple[float, float]:
    """Fit the odd-polynomial Taylor coefficients of ``p_diff(Δ)``.

    Returns ``(G1_fit, G3_taylor)`` in the **detuning** convention
    (``p_diff = G1_fit·Δ + (G3_taylor/6)·Δ³ + …``, Δ = ω_q − ω_d).

    ``delta_max_ghz``
        Scan half-range (GHz).  **The fit is only valid when this stays
        inside the cubic regime** (well below the p_diff turnover
        Δ_fold ≈ π/(2·T_eff)); a too-wide range lets fringe saturation
        corrupt the fit (e.g. G1 biased low by ~0.6× at 0.08 GHz for a
        10 ns π/2 pulse).

        - ``None`` (default): **adaptive** — shrink the range until G1
          converges, so the result auto-fits the qubit's linear regime.
        - ``float``: manual override (skips the adaptive search).

    Results are cached in ``_g3_cache`` keyed by
    ``(t_rabi_hash, omega_d, delta_max_tag)``.
    """
    import time

    tag = "adaptive" if delta_max_ghz is None else round(float(delta_max_ghz), 6)
    cache_key = _make_cache_key(t_rabi, omega_d, tag)
    if cache_key in _g3_cache:
        return _g3_cache[cache_key]

    t_start = time.time()

    if delta_max_ghz is not None:
        result = _fit_g3_at_delta_max(
            qubit, omega_d, t_rabi, t_global, n_scan, float(delta_max_ghz),
        )
    else:
        # Adaptive: shrink the scan range until G1 stabilizes (i.e. the
        # range has entered the cubic regime).  Descending candidates;
        # accept once consecutive G1 agree within 1.5%.
        candidates = (0.03, 0.02, 0.013, 0.008)
        G1_prev: float | None = None
        result = None
        for dmax in candidates:
            G1, G3 = _fit_g3_at_delta_max(
                qubit, omega_d, t_rabi, t_global, n_scan, dmax,
            )
            result = (G1, G3)
            if G1_prev is not None and abs(G1 - G1_prev) <= 0.015 * abs(G1):
                break
            G1_prev = G1

    _g3_cache[cache_key] = result

    elapsed = time.time() - t_start
    print(
        f"G3 Taylor calibration ({tag}): {elapsed:.1f}s, "
        f"G1={result[0]:.4f}, G3={result[1]:.1f}, "
        f"G3/G1={result[1] / result[0]:.1f} ns^2"
    )
    return result


# ---------------------------------------------------------------------------
# Internal helper: Route B — G1/G3 from the full off-diagonal sim kernel
# ---------------------------------------------------------------------------

def _calibrate_g3_kernel_full(
    qubit: object,
    omega_d: float,
    t_rabi: np.ndarray,
) -> tuple[float, float]:
    """Compute ``(G1, G3)`` from the full off-diagonal Heisenberg kernel.

    The constant-Δ cubic response needs the **triple** integral of the
    off-diagonal kernel, not the diagonal slice::

        G1      = ∫ k1(t) dt
        G3_full = ∭ k3(t1, t2, t3) dt1 dt2 dt3

    Both are taken on the orthogonal-Ramsey difference kernel
    ``(k_x − k_{-x})/2`` via ``KernelEstimator(method='sim',
    extract_off_diagonal=True)`` — one sesolve per pulse, exact.

    **Sign convention.**  The sim VZ kernel models the probe as
    ``+(φ/2)·σ_z`` while a physical detuning enters as ``−(Δ/2)·σ_z``;
    hence the raw kernel integrals are already in the ``δω = +φ = −Δ``
    convention used by the linear estimate ``G_freq`` (so they pair
    directly with ``delta_omega = −Δ`` in the Newton solve — no flip).

    Returns ``(G1, G3)`` ready to use as ``p_diff = G1·δω + (G3/6)·δω³``
    with ``δω = −Δ`` (the same convention as ``G_freq``).
    """
    from sqc.reconstruction.kernel import KernelEstimator

    ctrl_x = create_ramsey_pulse(
        t_rabi, 0.0, omega_d=omega_d, phase1=np.pi / 2, phase2=0.0, qubit=qubit,
    )
    ctrl_mx = create_ramsey_pulse(
        t_rabi, 0.0, omega_d=omega_d, phase1=np.pi / 2, phase2=np.pi, qubit=qubit,
    )
    est = KernelEstimator(
        mode='omega', method='sim', order=3, extract_off_diagonal=True,
    )
    rx = est.estimate_full(ctrl_x, qubit)
    rmx = est.estimate_full(ctrl_mx, qubit)
    tk = np.asarray(rx.t_samples, dtype=float)

    k1_diff = (np.asarray(rx.kernels[0]) - np.asarray(rmx.kernels[0])) / 2.0
    k3_diff = (np.asarray(rx.kernels[2]) - np.asarray(rmx.kernels[2])) / 2.0

    G1 = float(np.trapezoid(k1_diff, tk))
    G3 = float(np.trapezoid(np.trapezoid(np.trapezoid(k3_diff, tk, axis=0), tk, axis=0), tk, axis=0))
    return G1, G3


# ---------------------------------------------------------------------------
# Internal helper: transient single-point frequency measurement
# ---------------------------------------------------------------------------

def _solve_cubic_detuning(p_diff: float, G_lin: float, G3: float) -> float:
    """Invert ``p_diff = G_lin·δω + (G₃/6)·δω³`` for δω via Newton, with a
    fold guard.

    When ``G_lin`` and ``G₃`` have opposite signs, ``p_diff(δω)`` is
    non-monotone and turns over at ``|δω_fold| = sqrt(-2·G_lin/G₃)``. Beyond
    the fold the cubic is not invertible and Newton (seeded from the linear
    estimate) converges to the spurious far root of opposite sign — a large,
    sign-flipped error. If the solved δω lands past the fold, this returns the
    graceful order-1 linear estimate ``p_diff/G_lin`` instead.

    Falls back to the linear estimate when the cubic is degenerate
    (``G₃≈0`` or ``G_lin≈0``).
    """
    dw_lin = p_diff / G_lin if abs(G_lin) > 1e-12 else 0.0
    if abs(G3) <= 1e-10 or abs(G_lin) <= 1e-12:
        return dw_lin

    dw = dw_lin
    for _ in range(40):
        dw2 = dw * dw
        f_val = G_lin * dw + (G3 / 6.0) * dw * dw2 - p_diff
        f_prime = G_lin + (G3 / 2.0) * dw2
        if abs(f_prime) < 1e-15:
            break
        dw_new = dw - f_val / f_prime
        if abs(dw_new - dw) < 1e-12 * max(abs(dw), 1e-12):
            dw = dw_new
            break
        dw = dw_new

    fold_sq = -2.0 * G_lin / G3   # > 0 only when signs oppose
    if fold_sq > 0 and dw * dw > fold_sq:
        return dw_lin             # order-1 fallback (monotone, no overshoot)
    return dw


def _measure_frequency_transient(
    qubit: object,
    omega_d: float,
    t_rabi: np.ndarray,
    t_global: np.ndarray,
    flux: float = 0.0,
    order: int = 1,
    g3_source: Literal["fit", "kernel_full"] = "fit",
    g3_delta_max: float | None = None,
) -> float:
    """Transient-based single-point frequency measurement.

    Uses orthogonal Ramsey readout (R_y–R_x and R_y–R_{-x}) with tau=0
    to measure detuning via differential p_e and the control-pulse
    kernel sensitivity ``G_freq = ∫ k₁(t) dt``.

    Linear estimate: ``Δω = p_diff / G_freq`` → ``f = ω_d + Δ``.

    When ``order >= 3``, applies a cubic Newton correction solving::

        p_diff = G_freq·δω + (1/6)·G₃·δω³      (δω = −Δ convention)

    The cubic coefficient G₃ comes from one of two **equivalent** sources
    (both the constant-Δ triple-integral object, NOT the diagonal slice):

    - ``g3_source="fit"`` (default): odd-polynomial fit of ``p_diff(Δ)``
      (:func:`_calibrate_g3_taylor`), with adaptive scan range.
    - ``g3_source="kernel_full"``: full off-diagonal Heisenberg kernel
      ``∭k₃ dt³`` (:func:`_calibrate_g3_kernel_full`); no Δ scan needed.

    Both are converted into the δω = −Δ convention so the Newton result
    pairs consistently with the linear ``G_freq`` (``f = ω_d − δω``).

    ``g3_delta_max``
        Manual scan half-range (GHz) for ``g3_source="fit"``; ``None``
        uses the adaptive search.  Ignored for ``"kernel_full"``.

    Theory ref: idea/refactor/transient_frequency_theory.md §6;
                kernel/verify_transient_highorder.py (validated reference).
    """
    if order >= 3 and g3_source not in ("fit", "kernel_full"):
        raise ValueError(
            f"unknown g3_source={g3_source!r}; expected 'fit' or "
            f"'kernel_full' (the diagonal-kernel option was removed — "
            f"it is physically incorrect for constant-Δ inversion)."
        )

    n_levels = qubit.n_levels

    # -- set qubit to target flux ----------------------------------------
    t_sig = CONFIG.pulse.make_time(0, 300)
    Phi = FluxSignal(
        type=1 if flux != 0.0 else 0,
        t_list=t_sig,
        amplitude=float(flux),
        offset=0.0,
    )
    qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

    # -- build orthogonal readout pulses (tau=0) -------------------------
    tau = 0.0
    ctrl_x = create_ramsey_pulse(
        t_rabi, tau,
        omega_d=omega_d,
        phase1=np.pi / 2, phase2=0.0,
        qubit=qubit,
    )
    ctrl_mx = create_ramsey_pulse(
        t_rabi, tau,
        omega_d=omega_d,
        phase1=np.pi / 2, phase2=np.pi,
        qubit=qubit,
    )

    # -- run both measurements (use hamiltonian_on on t_global) ---------
    psi_e = basis(n_levels, 1)
    H_base = QobjEvo(
        qubit.H_list, tlist=qubit.mag_signal.t_list, order=1,
    )

    H_x = H_base + QobjEvo(
        ctrl_x.hamiltonian_on(t_global), tlist=t_global, order=1,
    )
    res_x = mesolve(
        H_x, qubit.state, t_global, [],
        e_ops=[psi_e * psi_e.dag()],
        options={"max_step": float(CONFIG.awg.dt)},
    )
    p_x = float(res_x.expect[0][-1])

    H_mx = H_base + QobjEvo(
        ctrl_mx.hamiltonian_on(t_global), tlist=t_global, order=1,
    )
    res_mx = mesolve(
        H_mx, qubit.state, t_global, [],
        e_ops=[psi_e * psi_e.dag()],
        options={"max_step": float(CONFIG.awg.dt)},
    )
    p_mx = float(res_mx.expect[0][-1])

    p_diff = (p_x - p_mx) / 2.0

    # -- linear sensitivity G_freq via order-1 omega kernel --------------
    # G_freq = dp_diff/d(δω) in the δω = +φ = −Δ convention (Virtual Z
    # probe +φ/2·σ_z vs detuning −Δ/2·σ_z).  Then δω = p_diff/G_freq = −Δ
    # and f = ω_d − δω = ω_d + Δ.  No κ unit conversion needed.
    from sqc.reconstruction.kernel import KernelEstimator

    estimator = KernelEstimator(
        mode='omega', method='exp', order=1, virtual_z_impl='math',
    )
    result_x = estimator.estimate_full(ctrl_x, qubit)
    k_omega_x = np.asarray(result_x.kernels[0], dtype=float)
    t_kernel = np.asarray(result_x.t_samples, dtype=float)
    result_mx = estimator.estimate_full(ctrl_mx, qubit)
    k_omega_mx = np.asarray(result_mx.kernels[0], dtype=float)

    # restore qubit state after kernel estimation side effects
    qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

    k_omega_diff = (k_omega_x - k_omega_mx) / 2.0
    G_freq = float(np.trapezoid(k_omega_diff, t_kernel))

    if abs(G_freq) < 1e-5:
        return float(omega_d)

    # -- linear estimate (always computed) ---------------------------------
    delta_omega = p_diff / G_freq

    # -- cubic Newton correction (order >= 3) -----------------------------
    # Solves  p_diff = G_lin·δω + (1/6)·G₃·δω³  in the **δω = −Δ** convention
    # (same as the linear estimate), so the final f = ω_d − δω is sign-
    # consistent.  k₂ ≈ 0 for the orthogonal Y-X sequence, so only the
    # cubic term enters.
    if order >= 3:
        if g3_source == "fit":
            # Route A: odd-poly fit gives (G1_fit, G3_taylor) in the
            # DETUNING (Δ) convention.  Convert to δω = −Δ by negating both.
            G1_fit, G3_taylor = _calibrate_g3_taylor(
                qubit, omega_d, t_rabi, t_global, delta_max_ghz=g3_delta_max,
            )
            G_lin = -G1_fit
            G3 = -G3_taylor
        elif g3_source == "kernel_full":
            # Route B: full off-diagonal triple integral, already in the
            # δω convention (raw sim-kernel integrals) — no flip.
            G_lin, G3 = _calibrate_g3_kernel_full(qubit, omega_d, t_rabi)
            # re-seat qubit after kernel side effects
            qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)
        else:
            raise ValueError(
                f"unknown g3_source={g3_source!r}; expected 'fit' or "
                f"'kernel_full' (the diagonal-kernel option was removed — "
                f"it is physically incorrect for constant-Δ inversion)."
            )

        # Solve the cubic p_diff = G_lin·δω + (G₃/6)·δω³ for δω, with a fold
        # guard that rejects the spurious far root past the turnover.
        if abs(G3) > 1e-10 and abs(G_lin) > 1e-12:
            delta_omega = _solve_cubic_detuning(p_diff, G_lin, G3)

    return float(omega_d - delta_omega)


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
        tau_list = CONFIG.pulse.make_time(0, 100)
        t_global = CONFIG.pulse.t_global.copy()
        frequency_list: list[float] = []

        for h in self.h_list:
            print(f"Calibrating f(Φ) at Φ={h:.4f} ...")
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
# FrequencyMeasurement — single-point f01 measurement (read-only)
# ===================================================================

@dataclass
class FrequencyMeasurement(Calibration):
    """Single-point qubit frequency measurement.

    Measures f01 at a single flux working point.  Returns the measured
    angular frequency (rad·GHz) and can also be packaged into a
    CalibrationTable for inclusion in calibration workflows.

    Methods
    -------
    - ``"ramsey"``: Ramsey τ-sweep + FFT peak → f01.  Default mode uses
      single-sweep (f_artificial=0.1 GHz) which assumes |Δ| < 0.1 GHz; for
      arbitrary |Δ| (e.g. inside a closed-loop search) set ``f_artificial``
      to None to use the double-sweep mode (slower but signed and unbounded).
    - ``"transient"``: τ=0 orthogonal Ramsey (R_y–R_x and R_y–R_{-x})
      differential readout + control-kernel sensitivity G_α = ∫k(t)dt to
      extract Δω directly.  Cheaper than Ramsey τ-sweep but relies on the
      weak-signal linear approximation; intended for |Δω| close to zero.

    Parameters
    ----------
    qubit : TransmonQubit
    method : str
        "ramsey" (default) or "transient".
    flux : float
        DC flux offset at which to measure (Φ₀).  Default 0.0 (sweet spot).
    tau_list : np.ndarray or None
        Ramsey free evolution times (ns).  Default make_time(0, 200).
        Ignored when method="transient".
    t_rabi : np.ndarray
        Rabi pulse time axis (ns).
    t_global : np.ndarray or None
        Global evolution time axis for mesolve (ns).
    f_artificial : float or None
        Single-sweep artificial detuning (GHz) for the Ramsey method.
        None selects double-sweep.  Default 0.1 GHz; closed-loop callers
        should pass None.  Ignored when method="transient".
    order : int
        Kernel order for transient method.  Default 1 (linear).  Set to
        3 for cubic Newton correction.  Values > 3 are accepted but use
        the same 5-point FD stencil.  Ignored when method="ramsey".
    g3_source : str
        Source of the cubic coefficient G₃ for the Newton correction when
        ``order >= 3``.  Both options are the **constant-Δ triple-integral**
        object (the correct one); the old diagonal-kernel shortcut was
        removed (physically wrong — ~170× too small):

        - ``"fit"`` (default): odd-polynomial fit of p_diff(Δ) with an
          adaptive scan range (:func:`_calibrate_g3_taylor`).
        - ``"kernel_full"``: full off-diagonal Heisenberg kernel
          ``∭k₃ dt³`` (:func:`_calibrate_g3_kernel_full`); no Δ scan,
          exact, useful as a cross-check of the fit.
    g3_delta_max : float or None
        Manual scan half-range (GHz) for ``g3_source="fit"``.  ``None``
        (default) uses the adaptive search that auto-fits the qubit's
        linear regime.  Ignored for ``"kernel_full"`` and ``order < 3``.
    """

    qubit: object
    method: Literal["ramsey", "transient"] = "ramsey"
    flux: float = 0.0
    order: int = 1
    g3_source: Literal["fit", "kernel_full"] = "fit"
    g3_delta_max: float | None = None

    tau_list: np.ndarray | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    t_global: np.ndarray | None = None
    f_artificial: float | None = 0.1
    omega_d: float | None = None

    def __post_init__(self):
        if self.tau_list is None:
            self.tau_list = CONFIG.pulse.make_time(0, 200)
        if self.t_global is None:
            self.t_global = CONFIG.pulse.t_global.copy()

    # ------------------------------------------------------------------
    def measure(self, flux: float | None = None,
                omega_d: float | None = None) -> float:
        """Return f01 (angular, rad·GHz) at the given flux.

        Parameters
        ----------
        flux : float or None
            Flux offset (Φ₀).  Defaults to self.flux when None.
        omega_d : float or None
            Drive/reference frequency (angular, rad·GHz).  The measurement
            returns ``omega_d + measured detuning``, so this sets the point the
            detuning is measured *relative to*.  Resolution order: this
            argument → ``self.omega_d`` → ``self.qubit.frequency`` (the
            sweet-spot construction frequency).

            Passing an explicit ``omega_d`` is what real experiments do: the
            true f01 is unknown, so the drive is set at a spectroscopy estimate
            and the measurement reports the residual detuning.  It also enables
            dynamic re-centering (feedback): feed the previous estimate back as
            ``omega_d`` to keep the transient measurement inside its linear
            regime.  The default (``qubit.frequency``) measures detuning
            relative to the sweet spot — convenient in simulation but assumes
            the answer is already known.

        Returns
        -------
        float
            Measured qubit frequency, signed (angular, rad·GHz).
        """
        flux_val = self.flux if flux is None else float(flux)
        omega_d = (
            omega_d if omega_d is not None
            else self.omega_d if self.omega_d is not None
            else self.qubit.frequency
        )
        match self.method:
            case "ramsey":
                return _fit_ramsey_frequency(
                    self.qubit, omega_d, self.tau_list,
                    self.t_rabi, self.t_global,
                    flux=flux_val, f_artificial=self.f_artificial,
                )
            case "transient":
                return _measure_frequency_transient(
                    self.qubit, omega_d,
                    self.t_rabi, self.t_global,
                    flux=flux_val, order=self.order,
                    g3_source=self.g3_source,
                    g3_delta_max=self.g3_delta_max,
                )

    # ------------------------------------------------------------------
    def calibrate(self) -> CalibrationTable:
        """Standalone single-point measurement → CalibrationTable."""
        f_meas = self.measure(self.flux)
        return CalibrationTable(
            name=f"frequency_measurement_{self.method}",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f01",
            inputs=np.array([self.flux]),
            outputs=np.array([f_meas]),
            fit_params={
                "method": self.method,
                "omega_d": float(self.qubit.frequency),
            },
            metadata={},
        )


# ===================================================================
# SinglePointFrequencyCalibration — single-point frequency tuning
# ===================================================================

@dataclass
class SinglePointFrequencyCalibration(Calibration):
    """Single-point qubit frequency tuning.

    Drives the qubit frequency f_q(V) to a target f_target by adjusting
    the flux bias V.  Currently supports closed-loop feedback
    (Vepsalainen 2022).  Designed to admit future single-point control
    strategies (e.g. gradient descent, feedforward) via the ``method``
    field — new methods should add a ``case`` branch in :meth:`calibrate`.

    Methods
    -------
    - ``"closed_loop"`` (default): secant, bisection, or gradient-based
      root-finding on r(V) = f_q(V) − f_target = 0.  Per-iteration
      measurement is delegated to an internal :class:`FrequencyMeasurement`
      instance, which can be Ramsey-based (robust, default) or transient-based.

    Parameters
    ----------
    qubit : TransmonQubit
    method : str
        Tuning strategy.  Currently only "closed_loop".

    Closed-loop parameters
    ----------------------
    f_target : float or None
        Target qubit frequency (angular, GHz·2π).  Required.
    V_a, V_b : float or None
        Voltage bounds bracketing the target.  Required; typically
        provided by a preceding :class:`FluxResponseCalibration`.
    epsilon_f : float
        Convergence tolerance on |r| (GHz·2π).  Default 1e-4.
    max_iter : int
        Maximum iterations.  Default 20.
    measure_method : str
        Per-iteration frequency measurement strategy passed to the inner
        :class:`FrequencyMeasurement`:

        - ``"ramsey"`` (default, robust): Ramsey FFT with double-sweep,
          handles arbitrary |Δω| during the search.
        - ``"transient"``: τ=0 orthogonal Ramsey + kernel sensitivity.
          Faster but relies on weak-signal linear approximation; risky
          when the search may probe far from the root.
    bracket_tightening : bool
        If True (default), tighten [V_a, V_b] each iteration based on the
        sign of r (regula falsi).  Improves robustness; disable to observe
        pure secant behaviour.
    step_method : str
        - ``"secant"`` (default): secant method, superlinear convergence
          (1–3 iterations typical) when paired with bracket_tightening.
        - ``"bisection"``: classical bisection, O(log₂(range/ε)) iterations
          (~10–15) with a clean exponential bracket-width trend, ideal
          for visualisation and diagnostics.  Auto-splits the bracket
          when both endpoints share the same residual sign (handles even
          f(Φ) crossing the sweet spot).
        - ``"gradient"``: damped secant (numerical-gradient Newton step)
          with clamping and best-point tracking.  Does **not** require
          pre-bracketing V_a/V_b — starts from V_seed with a probe step,
          then uses a damping factor to suppress overshoot from noisy
          gradient estimates.  See gradient-specific parameters below.

    Gradient-method parameters (step_method="gradient")
    ---------------------------------------------------
    V_seed : float or None
        Starting flux bias (Φ₀).  If None, falls back to (V_a+V_b)/2
        when V_a, V_b are both set.  Required when V_a, V_b are absent.
    damping : float
        Damping factor ∈ (0, 1] to suppress overshoot.  Default 0.8.
    first_bias_step : float
        Probe step size (Φ₀) used on the first iteration or when the
        gradient estimate is undefined (Δe=0).  Default 0.01.
    max_bias_step : float
        Per-step clamp (Φ₀).  Prevents a near-zero denominator from
        flinging the bias far away.  Default 0.05.

    Inner-measurement parameters (forwarded to FrequencyMeasurement)
    ----------------------------------------------------------------
    tau_list, t_rabi, t_global
        See :class:`FrequencyMeasurement`.  The inner measurement always
        uses ``f_artificial=None`` (double-sweep) because the search may
        probe flux far from the sweet spot.
    """

    qubit: object
    method: Literal["closed_loop"] = "closed_loop"

    # Closed-loop parameters
    f_target: float | None = None
    V_a: float | None = None
    V_b: float | None = None
    epsilon_f: float = 1e-4
    max_iter: int = 20
    measure_method: Literal["ramsey", "transient"] = "ramsey"
    bracket_tightening: bool = True
    step_method: Literal["secant", "bisection", "gradient"] = "secant"

    # Transient-measurement tuning (forwarded to the inner FrequencyMeasurement
    # when measure_method="transient"; ignored for "ramsey").
    order: int = 1
    g3_source: Literal["fit", "kernel_full"] = "fit"

    # Gradient-method parameters (step_method="gradient")
    V_seed: float | None = None
    damping: float = 0.8
    first_bias_step: float = 0.01
    # Per-step flux clamp. Near the sweet spot ∂f_q/∂V→0, so the damped-secant
    # step blows up; too large a clamp flings V across the sweet spot (into the
    # region where the transient discriminator is non-linear and mis-reads),
    # stalling the coarse search. 0.02 keeps the search local enough to descend.
    max_bias_step: float = 0.02

    # Optional early-stop predicate (step_method="gradient" only). Called at
    # the end of each iteration with a state dict; return True to stop early
    # (e.g. gradient flattened, stall, cost budget). Ignored by secant /
    # bisection, whose loops are not step-instrumented. Default None keeps the
    # loop byte-for-byte unchanged.
    stop_predicate: Optional[Callable[[dict], bool]] = None

    # Drive-frequency policy (gradient step only). Sets the omega_d the inner
    # measurement uses each iteration — the *reference* the detuning is measured
    # relative to. The control error is ALWAYS f_q − f_target (absolute), so the
    # policy only affects measurement linearity, never the convergence target.
    #   "sweet"  : omega_d = qubit.frequency (sweet spot). Default; back-compat.
    #   "target" : omega_d = f_target (LOCKED phase — measured detuning IS the
    #              control error, simplest closed form).
    #   "track"  : omega_d predicts f_q at the point about to be measured,
    #              omega_d = f_hat_prev + s_hat·(V − V_prev) (TRACKING phase —
    #              keeps |detuning| inside the protocol's linear window as the
    #              search moves far from the sweet spot).
    #   callable : user function(state) -> omega_d for arbitrary feedback laws.
    drive_policy: "str | Callable[[dict], float]" = "sweet"
    # Sensitivity source for "track" prediction and diagnostics.
    #   "secant": s_hat = de/dV from the last two measurements (model-free,
    #             matches a real experiment). "model": qubit.frequency_sensitivity.
    sensitivity_source: Literal["secant", "model"] = "secant"
    # Initial frequency estimate f̂_{q,0} for the "track" policy's first drive
    # (f_{d,0} = omega_d_seed). Without it a fresh track stage falls back to the
    # sweet spot on iter 0 — fatal when the target is far (the transient sees the
    # full offset). A staged pipeline passes the previous stage's estimate here.
    omega_d_seed: float | None = None
    # Linear-range guard: if set, each iteration flags out_of_range when
    # |detuning| > rho·linear_range, surfaced in the stop_predicate state so a
    # caller can trigger re-acquisition. None disables the guard.
    linear_range: float | None = None
    rho: float = 0.6
    # Convergence de-bounce (gradient step): require |e| ≤ epsilon on this many
    # consecutive iterations before declaring convergence. Default 1 = current
    # behaviour (stop on first hit). >1 re-measures in place to confirm, guarding
    # against a noisy measurement that momentarily dips under tolerance.
    converge_streak: int = 1

    # Inner-measurement parameters (forwarded to FrequencyMeasurement)
    tau_list: np.ndarray | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    t_global: np.ndarray | None = None

    def __post_init__(self):
        if self.tau_list is None:
            self.tau_list = CONFIG.pulse.make_time(0, 200)
        if self.t_global is None:
            self.t_global = CONFIG.pulse.t_global.copy()

        # Inner measurement: dual-sweep Ramsey so |Δω| is unconstrained
        # during the search.  For transient mode, forward the kernel order /
        # G₃ source so callers can widen the accurate range (order=3).
        self._meas = FrequencyMeasurement(
            qubit=self.qubit,
            method=self.measure_method,
            tau_list=self.tau_list,
            t_rabi=self.t_rabi,
            t_global=self.t_global,
            f_artificial=None,
            order=self.order,
            g3_source=self.g3_source,
        )

        # Non-default drive policies steer omega_d per iteration, which only the
        # gradient loop threads (secant/bisection are not step-instrumented).
        if (not (isinstance(self.drive_policy, str)
                 and self.drive_policy == "sweet")
                and self.step_method != "gradient"):
            raise ValueError(
                f"drive_policy={self.drive_policy!r} requires "
                f"step_method='gradient'; got {self.step_method!r}."
            )

    # ------------------------------------------------------------------
    def _resolve_omega_d(self, state: dict) -> float:
        """Return the drive frequency for the iteration described by *state*.

        The control error is always ``f_q − f_target`` (absolute); this only
        chooses the reference the *detuning* is measured against, i.e. it keeps
        the measurement inside its linear window. See the ``drive_policy`` field.
        """
        policy = self.drive_policy
        if callable(policy):
            return float(policy(state))
        if policy == "sweet":
            return float(self.qubit.frequency)
        if policy == "target":
            return float(self.f_target)
        if policy == "track":
            # First iteration (no prior estimate): use the seed estimate if
            # given (f_{d,0} = f̂_{q,0}), else fall back to the sweet spot.
            f_prev = state.get("f_hat_prev")
            if f_prev is None:
                return float(self.omega_d_seed
                             if self.omega_d_seed is not None
                             else self.qubit.frequency)
            s_hat = state.get("s_hat")
            if s_hat is None:
                return float(f_prev)          # no slope yet: reuse last estimate
            dV = state["V"] - state["V_prev"]
            return float(f_prev + s_hat * dV)  # predict f_q at the new point
        raise ValueError(f"unknown drive_policy={policy!r}.")

    # ------------------------------------------------------------------
    def calibrate(self) -> CalibrationTable:
        match self.method:
            case "closed_loop":
                return self._calibrate_closed_loop()

    # ------------------------------------------------------------------
    # Closed-loop feedback (Vepsalainen 2022)
    # ------------------------------------------------------------------

    def _calibrate_closed_loop(self) -> CalibrationTable:
        """Closed-loop frequency tuning.

        Iteratively adjusts flux bias voltage to drive qubit frequency
        to f_target via root-finding on r(V) = f_Q(V) - f_target = 0.

        For ``step_method`` secant/bisection: requires V_a, V_b bracketing
        the target.  For ``step_method="gradient"``: requires only V_seed
        (or V_a/V_b as fallback); no bracketing needed.
        """
        if self.f_target is None:
            raise ValueError("f_target is required for closed_loop method.")

        f_target = self.f_target
        epsilon = self.epsilon_f

        if self.step_method == "gradient":
            return self._closed_loop_gradient(f_target, epsilon)

        if self.V_a is None or self.V_b is None:
            raise ValueError(
                "V_a and V_b are required for closed_loop method. "
                "Run FluxResponseCalibration first to bracket the target."
            )

        V_lo, V_hi = min(self.V_a, self.V_b), max(self.V_a, self.V_b)

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
        r_n = self._meas.measure(V_n) - f_target
        n_iter = 0

        V_prev = V_lo
        r_prev = self._meas.measure(V_lo) - f_target

        # Sign of the residual at the lower bracket endpoint. Regula-falsi
        # must replace whichever endpoint shares r_n's sign — hardcoding
        # "r_n > 0 -> shrink V_hi" only holds when f increases with V, and
        # breaks for monotonically-decreasing f(Phi) (qubit freq vs |flux|).
        s_lo = 1.0 if r_prev > 0 else -1.0

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

            f_n = self._meas.measure(V_n)
            r_n = f_n - f_target
            n_iter += 1

            history.append({
                "iter": n_iter, "V": float(V_n),
                "f": float(f_n), "residual": float(r_n),
            })

            # Bracket tightening (regula falsi): replace the endpoint whose
            # residual shares r_n's sign, so the root stays bracketed
            # regardless of whether f increases or decreases with V.
            if self.bracket_tightening:
                if (r_n > 0) == (s_lo > 0):
                    V_lo = V_n
                else:
                    V_hi = V_n

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
        r_lo = self._meas.measure(V_lo) - f_target
        r_hi = self._meas.measure(V_hi) - f_target

        # Auto-split bracket if both ends have the same sign (even f(Φ))
        if r_lo * r_hi > 0:
            V_mid0 = (V_lo + V_hi) / 2.0
            r_mid0 = self._meas.measure(V_mid0) - f_target
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
            r_mid = self._meas.measure(V_mid) - f_target
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
        r_mid = self._meas.measure(V_mid) - f_target
        return self._build_result(f_target, V_mid, r_mid, n_iter, epsilon, history)

    # ------------------------------------------------------------------
    def _closed_loop_gradient(
        self, f_target: float, epsilon: float,
    ) -> CalibrationTable:
        """Gradient-based (damped-secant) closed-loop iteration.

        Each iteration measures f_n → residual e_n = f_n − f_target,
        then updates the bias via a damped numerical-gradient Newton step
        with per-step clamping and best-point tracking.

        Does **not** require pre-bracketing V_a/V_b: starts from V_seed
        (or the midpoint of V_a/V_b if V_seed is None).  The first step
        is a fixed probe step in the direction of −sign(e_n); subsequent
        steps use ΔV/Δe with damping.  When Δe == 0 the gradient is
        undefined and the probe step is reused.

        Internally delegates to :class:`DampedSecantTracker` — the shared
        stepwise controller also used by the event-driven state machine.
        """
        from sqc.calibration.frequency_control import (
            DampedSecantTracker,
            FrequencyEstimate,
        )

        # ---- seed --------------------------------------------------------
        if self.V_seed is not None:
            V_n = float(self.V_seed)
        elif self.V_a is not None and self.V_b is not None:
            V_n = (self.V_a + self.V_b) / 2.0
        else:
            raise ValueError(
                "V_seed is required for gradient method when V_a/V_b "
                "are not both set."
            )

        # Optional clamp bounds from V_a/V_b (if provided)
        V_lo = (
            float(min(self.V_a, self.V_b))
            if (self.V_a is not None and self.V_b is not None)
            else None
        )
        V_hi = (
            float(max(self.V_a, self.V_b))
            if (self.V_a is not None and self.V_b is not None)
            else None
        )

        # --- seed measurement: resolve drive, measure, record diagnostics ---
        omega_d = self._resolve_omega_d(
            {"iter": 0, "V": V_n, "V_prev": V_n,
             "f_hat_prev": None, "s_hat": None, "f_target": f_target}
        )
        f_n = self._meas.measure(V_n, omega_d=omega_d)
        e_n = f_n - f_target
        delta_n = f_n - omega_d                    # local detuning (metrology)

        # Track best point in case later iterations diverge
        best_V = V_n
        best_residual = e_n
        best_abs_e = abs(e_n)

        history: list[dict] = []
        history.append({
            "iter": 0, "V": float(V_n),
            "f": float(f_n), "residual": float(e_n),
            "omega_d": float(omega_d), "delta": float(delta_n),
            "is_best": True,
        })

        # ---- shared stepwise tracker ------------------------------------
        tracker = DampedSecantTracker(
            f_target=f_target,
            damping=self.damping,
            first_bias_step=self.first_bias_step,
            max_bias_step=self.max_bias_step,
            V_lo=V_lo, V_hi=V_hi,
            converge_streak=self.converge_streak,
            max_iter=self.max_iter,
            epsilon_f=epsilon,
            stop_predicate=None,   # checked manually for exact equivalence
        )
        seed_est = FrequencyEstimate(frequency=f_n)
        snapshot = tracker.initialize(seed_est, V_n)

        n_streak = max(1, int(self.converge_streak))

        # ---- main loop ---------------------------------------------------
        stopped_by: str | None = None
        while True:
            proposal = tracker.propose(snapshot)
            if proposal.converged:
                stopped_by = proposal.diagnostics.get("reason", "max_iter")
                if stopped_by == "already_converged":
                    stopped_by = "tolerance"
                break

            # Save pre-step state for omega_d resolution
            V_prev = snapshot.V
            f_prev_hat = snapshot.f
            n_iter = snapshot.n_iter

            # --- resolve drive frequency for this iteration ---
            if self.sensitivity_source == "model":
                s_hat = float(self.qubit.frequency_sensitivity(V_prev))
            else:
                s_hat = proposal.s_hat
            omega_d = self._resolve_omega_d({
                "iter": n_iter + 1, "V": proposal.V_next,
                "V_prev": V_prev,
                "f_hat_prev": f_prev_hat, "s_hat": s_hat,
                "f_target": f_target,
            })

            # Measure at the proposed point
            f_n = self._meas.measure(proposal.V_next, omega_d=omega_d)
            delta_n = f_n - omega_d

            # Linear-range guard
            out_of_range = (
                self.linear_range is not None
                and abs(delta_n) > self.rho * self.linear_range
            )

            # Accept measurement into tracker
            estimate = FrequencyEstimate(frequency=f_n)
            result = tracker.accept(
                snapshot, estimate, proposal,
                extra_predicate_state={
                    "delta": float(delta_n),
                    "out_of_range": bool(out_of_range),
                    "history": history,
                },
            )
            snapshot = result.snapshot
            e_n = snapshot.e
            n_iter = snapshot.n_iter

            # Track best point
            is_best = result.is_best
            if is_best:
                best_abs_e = abs(e_n) if e_n is not None else float("inf")
                best_residual = e_n if e_n is not None else float("inf")
                best_V = snapshot.V

            history.append({
                "iter": n_iter, "V": float(snapshot.V),
                "f": float(f_n), "residual": float(e_n) if e_n is not None else float("inf"),
                "omega_d": float(omega_d), "delta": float(delta_n),
                "out_of_range": bool(out_of_range),
                "is_best": is_best,
            })

            # Check tracker stop reason
            if result.stopped_by != "running":
                stopped_by = result.stopped_by
                break

            # Optional early stop: hand the caller a rich per-step state so it
            # can switch phase on gradient/stall/divergence/budget signals that
            # a plain residual threshold cannot express.
            if self.stop_predicate is not None:
                state = {
                    "iter": n_iter,
                    "V": float(snapshot.V),
                    "residual": float(e_n) if e_n is not None else float("inf"),
                    "abs_residual": float(abs(e_n)) if e_n is not None else float("inf"),
                    "step": float(proposal.step),
                    "de": float(e_n - snapshot.e_prev)
                    if (e_n is not None and snapshot.e_prev is not None)
                    else 0.0,
                    "dV": float(snapshot.V - snapshot.V_prev),
                    "delta": float(delta_n),
                    "out_of_range": bool(out_of_range),
                    "best_abs_residual": float(snapshot.best_abs_e),
                    "history": history,
                }
                if self.stop_predicate(state):
                    stopped_by = "predicate"
                    break

        if stopped_by is None:
            e_final = snapshot.e
            stopped_by = (
                "tolerance"
                if (e_final is not None and abs(e_final) <= epsilon)
                else "max_iter"
            )

        return self._build_result(
            f_target, best_V, best_residual,
            n_iter, epsilon, history,
            extra={
                "damping": self.damping,
                "first_bias_step": self.first_bias_step,
                "max_bias_step": self.max_bias_step,
                "stopped_by": stopped_by,
            },
        )

    # ------------------------------------------------------------------
    def _build_result(
        self, f_target: float, V_opt: float, residual: float,
        n_iter: int, epsilon: float, history: list[dict],
        extra: dict | None = None,
    ) -> CalibrationTable:
        """Package closed-loop result into a CalibrationTable."""
        fit_params = {
            "method": self.method,
            "step_method": self.step_method,
            "measure_method": self.measure_method,
            "f_target": float(f_target),
            "V_opt": float(V_opt),
            "n_iter": n_iter,
            "residual": float(residual),
            "converged": abs(residual) <= epsilon,
            "history": history,
        }
        if extra is not None:
            fit_params.update(extra)
        return CalibrationTable(
            name="frequency_closed_loop",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f01",
            inputs=np.array([V_opt]),
            outputs=np.array([f_target + residual]),
            fit_params=fit_params,
            metadata={},
        )
