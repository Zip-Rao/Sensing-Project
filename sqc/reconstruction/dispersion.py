"""sqc.reconstruction.dispersion — shared Transmon dispersion inversion.

Used by cryoscope and delay_ramsey reconstruction to convert
frequency shift Δω → flux offset h analytically.
"""
from __future__ import annotations

import numpy as np


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
    bias = getattr(qubit, "flux_bias", getattr(qubit, "flux", 0.0))
    return total_flux - bias
