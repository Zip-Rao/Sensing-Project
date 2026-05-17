"""sqc.reconstruction.dispersion — shared Transmon dispersion inversion.

Used by cryoscope and delay_ramsey reconstruction to convert
frequency shift Δω → flux offset h analytically, and to perform
model-guided phase unwrapping anchored to the analytical dispersion.
"""
from __future__ import annotations

import numpy as np


def _qubit_bias(qubit) -> float:
    """Return the qubit's DC flux bias (Φ₀)."""
    bias = getattr(qubit, "flux_bias", None)
    if bias is None:
        bias = getattr(qubit, "flux", 0.0)
    return float(bias)


def omega_q_at_flux(qubit, flux_total) -> np.ndarray:
    """Analytical Transmon frequency at absolute flux Φ_total.

    f_Q(Φ) = sqrt(8·EJ_0·|cos(π·Φ)|·EC) - EC    (ħ=1, all angular units)

    Parameters
    ----------
    qubit : TransmonQubit or QubitSpec
        Must provide EC, EJ (or EJ_0).
    flux_total : float or np.ndarray
        Absolute flux (Φ₀), not offset from bias.
    """
    EC = qubit.EC
    EJ0 = getattr(qubit, "EJ_0", qubit.EJ)
    return np.sqrt(8.0 * EJ0 * np.abs(np.cos(np.pi * np.asarray(flux_total, dtype=float))) * EC) - EC


def cryoscope_phase_theory(
    qubit, flux_at_each_point: np.ndarray, tau: float,
    omega_d: float | None = None,
) -> np.ndarray:
    """Theoretical accumulated phase for a constant-flux Ramsey readout.

    φ_theory[i] = (ω_q(Φ_bias + flux_at_each_point[i]) - ω_d) · τ

    Used by CryoscopeCalibration where each scan point applies a
    square pulse of constant height h for duration τ.
    """
    if omega_d is None:
        omega_d = float(qubit.frequency)
    flux_bias = _qubit_bias(qubit)
    omega_q = omega_q_at_flux(
        qubit, flux_bias + np.asarray(flux_at_each_point, dtype=float)
    )
    return (omega_q - float(omega_d)) * float(tau)


def cumulative_phase_theory(
    qubit, flux_signal_t: np.ndarray, flux_signal_h: np.ndarray,
    query_times: np.ndarray, omega_d: float | None = None,
    window: tuple[float, float] | None = None,
) -> np.ndarray:
    """Cumulative theoretical phase for an arbitrary flux waveform.

    Default (window=None):
        φ_theory(t_d) = ∫₀^{t_d} (ω_q(Φ_bias + h(t)) - ω_d) dt
        Used by CryoscopeExperiment.

    With window=(t0, t1):
        φ_theory(t_d) = ∫_{t0+t_d}^{t1+t_d} (ω_q(Φ_bias + h(t)) - ω_d) dt
        Used by DelayRamseyExperiment, where each t_d slides the
        free-evolution window across the flux waveform.

    Parameters
    ----------
    qubit : TransmonQubit
    flux_signal_t : np.ndarray
        Fine-grained time grid of the flux signal (ns).
    flux_signal_h : np.ndarray
        Flux offset h(t) (relative to qubit.flux_bias, Φ₀).
    query_times : np.ndarray
        Times at which to evaluate the integral.
    omega_d : float, optional
        Defaults to qubit.frequency.
    window : (t0, t1), optional
        Sliding-window mode: integrate over [t0+t_d, t1+t_d].
        If None, integrate over [0, t_d].
    """
    if omega_d is None:
        omega_d = float(qubit.frequency)
    flux_bias = _qubit_bias(qubit)
    flux_t = np.asarray(flux_signal_t, dtype=float)
    flux_h = np.asarray(flux_signal_h, dtype=float)
    omega_q = omega_q_at_flux(qubit, flux_bias + flux_h)
    integrand = omega_q - float(omega_d)
    cum_phase = np.concatenate((
        [0.0],
        np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(flux_t)),
    ))
    q = np.asarray(query_times, dtype=float)
    if window is None:
        return np.interp(q, flux_t, cum_phase)
    t0, t1 = window
    upper = np.interp(t1 + q, flux_t, cum_phase,
                      left=cum_phase[0], right=cum_phase[-1])
    lower = np.interp(t0 + q, flux_t, cum_phase,
                      left=cum_phase[0], right=cum_phase[-1])
    return upper - lower


def unwrap_phase_with_model(
    varphi_raw: np.ndarray, varphi_theory: np.ndarray,
) -> np.ndarray:
    """Model-guided phase unwrap.

    Single source of truth for the 2π branch-picking step used by
    cryoscope and delay_ramsey, ensuring the absolute-phase
    convention is consistent across experiment and calibration paths.

    The 2π branch is chosen to minimise |φ_theory - φ_raw| per point.

    Parameters
    ----------
    varphi_raw : np.ndarray
        Raw wrapped phases from arctan2 (rad).
    varphi_theory : np.ndarray
        Theoretical phases of the same shape, from any physics model.
        Use `cryoscope_phase_theory` for constant-flux calibration
        sweeps, `cumulative_phase_theory` for cryoscope truncation
        experiments.

    Returns
    -------
    np.ndarray
        Unwrapped phases (rad).
    """
    varphi_raw = np.asarray(varphi_raw, dtype=float)
    varphi_theory = np.asarray(varphi_theory, dtype=float)
    if varphi_theory.shape != varphi_raw.shape:
        raise ValueError(
            f"varphi_theory shape {varphi_theory.shape} does not match "
            f"varphi_raw shape {varphi_raw.shape}"
        )
    n_wraps = np.round((varphi_theory - varphi_raw) / (2.0 * np.pi))
    return varphi_raw + 2.0 * np.pi * n_wraps


def qubit_inverse_frequency(dphi_dt: np.ndarray, qubit) -> np.ndarray:
    """Map angular frequency shift Δω → flux h via analytical Transmon dispersion.

    f_Q(Φ) = sqrt(8·EJ·|cos(π·Φ)|·EC) - EC    (ħ=1, all angular units)

    Inverting:  cos(π·Φ) = (f_Q + EC)² / (8·EJ·EC)
                Φ = arccos(clip(ratio, 0, 1)) / π
                h = Φ - flux_bias

    Parameters
    ----------
    dphi_dt : np.ndarray
        Phase derivative (rad/ns), equals the angular frequency shift
        Δω = ω_Q(Φ_bias + h) - ω_Q(Φ_bias) in natural units.
    qubit : TransmonQubit or QubitSpec
        Must provide frequency, EC, EJ (or EJ_0), and flux (or flux_bias).

    Returns
    -------
    np.ndarray
        Flux offset h (Φ₀) corresponding to dphi_dt.
    """
    f_q = qubit.frequency
    f_target = f_q + dphi_dt
    EC = qubit.EC
    EJ0 = getattr(qubit, "EJ_0", qubit.EJ)
    ratio = (f_target + EC) ** 2 / (8.0 * EC * EJ0)
    ratio = np.clip(ratio, 0.0, 1.0)
    total_flux = np.arccos(ratio) / np.pi
    return total_flux - _qubit_bias(qubit)
