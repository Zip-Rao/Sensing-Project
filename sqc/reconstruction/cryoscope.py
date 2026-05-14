"""sqc.reconstruction.cryoscope — CryoscopeReconstruction.

Cryoscope waveform reconstruction from phase-vs-truncation data.

Replaces Analysis.get_signal_from_cryoscope + get_h_from_phi.

See _cryoscope_implementation.md §3.5 for algorithm details.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from sqc.calibration.base import CalibrationTable
from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.reconstruction.base import Reconstruction


def _qubit_inverse_frequency(
    dphi_dt: np.ndarray,
    qubit,
) -> np.ndarray:
    """Map angular frequency shift Δω → flux h via analytical Transmon dispersion.

    f_Q(Φ) = sqrt(8·EJ·|cos(π·Φ)|·EC) - EC    (all in angular units, ħ=1)

    Inverting:  cos(π·Φ) = (f_Q + EC)² / (8·EJ·EC)
               Φ = arccos(clip(ratio, 0, 1)) / π
               h = Φ - flux_bias

    Parameters
    ----------
    dphi_dt : np.ndarray
        Phase derivative (rad/ns), which equals the angular frequency
        shift Δω = ω_Q(Φ_bias+h) - ω_Q(Φ_bias) in natural units (ħ=1).
    qubit
        TransmonQubit or QubitSpec providing frequency, EC, EJ, and flux.

    Returns
    -------
    np.ndarray
        Flux offset h (Φ₀) corresponding to dphi_dt.
    """
    # Current qubit frequency at its bias point (angular, GHz·2π)
    f_q = qubit.frequency
    # Frequency after shift: ω_Q(Φ+h) = ω_Q(Φ) + Δω = f_q + dφ/dt
    f_target = f_q + dphi_dt

    EC = qubit.EC
    # EJ at zero flux (maximum)
    if hasattr(qubit, "EJ_0"):
        EJ0 = qubit.EJ_0
    elif hasattr(qubit, "EJ"):
        EJ0 = qubit.EJ
    else:
        raise TypeError("qubit must have EJ or EJ_0 attribute")

    ratio = (f_target + EC) ** 2 / (8.0 * EC * EJ0)
    ratio = np.clip(ratio, 0.0, 1.0)
    total_flux = np.arccos(ratio) / np.pi

    # Subtract qubit's DC bias flux
    if hasattr(qubit, "flux_bias"):
        bias = qubit.flux_bias
    elif hasattr(qubit, "flux"):
        bias = qubit.flux
    else:
        bias = 0.0

    return total_flux - bias


@dataclass
class CryoscopeReconstruction(Reconstruction):
    """Cryoscope waveform reconstruction.

    Core physics:  dφ/dt = 2π·Δf(h(t))

    Two inversion strategies:
    - ``"calibration"``: use pre-measured φ(h) lookup table via
      ``CalibrationTable.inverse()``. Requires ``calibration``.
    - ``"response"``: use analytical Transmon frequency-flux relation
      directly. Requires ``qubit`` (needs EC, EJ, flux).

    Optional Savitzky-Golay pre-smoothing on the phase derivative
    (``use_sg_filter=True``) for noisy data.

    Parameters
    ----------
    tau : float
        Calibration square-pulse length (ns).
    inversion : str
        - ``"calibration"`` — map φ → h via calibration table.
        - ``"response"`` — map Δf → h via qubit dispersion relation.
    calibration : CalibrationTable or None
        Required for inversion="calibration".
    qubit : TransmonQubit or None
        Required for inversion="response".
    use_sg_filter : bool
        If True, apply Savitzky-Golay smoothing before differentiation.
    sg_window : int
        SG window length (odd). Only used when use_sg_filter=True.
    sg_poly : int
        SG polynomial order. Only used when use_sg_filter=True.
    """

    tau: float = field(
        default_factory=lambda: CONFIG.reconstruction.cryoscope_tau
    )
    inversion: Literal["calibration", "response"] = "calibration"
    calibration: CalibrationTable | None = None
    qubit: object | None = None  # TransmonQubit or QubitSpec
    use_sg_filter: bool = False
    sg_window: int = 7
    sg_poly: int = 2

    def reconstruct(
        self,
        measurement,
        kernel=None,
        calibration=None,
        dt: float | None = None,
    ) -> FluxSignal:
        """Reconstruct flux waveform h(t) from cryoscope phase data.

        Physics:  dφ/dt = 2π·Δf(h(t))
        →  h(t) = f_Q⁻¹(Δf)   or   h(t) = φ_cal⁻¹(dφ/dt · tau)

        Parameters
        ----------
        measurement : ExperimentResult
            Must contain data["varphi"] and axes["trunc"].
        kernel : optional
            Not used.
        calibration : CalibrationTable, optional
            Override instance calibration table.
        dt : float, optional
            Time step (ns). Computed from trunc axis if not given.

        Returns
        -------
        FluxSignal
            Reconstructed flux signal h(t), type=8, t_list = trunc times.
        """
        varphi = np.asarray(measurement.data["varphi"], dtype=float)
        trunc = np.asarray(measurement.axes["trunc"], dtype=float)

        # CryoscopeExperiment stores trunc in decreasing order, but varphi
        # has already been reversed to increasing order via [::-1] at
        # experiments/cryoscope.py L117.  Only reverse the time axis.
        if len(trunc) > 1 and trunc[0] > trunc[-1]:
            trunc = trunc[::-1]

        if dt is None:
            dt = float(trunc[1] - trunc[0])

        # --- step 1: dφ/dt (rad/ns) with optional SG smoothing ---
        if self.use_sg_filter:
            from scipy.signal import savgol_filter

            window = self.sg_window
            if window >= len(varphi):
                window = max(3, len(varphi) // 2 * 2 - 1)
            dphi_dt = savgol_filter(
                varphi, window_length=window, polyorder=self.sg_poly,
                deriv=1, delta=dt,
            )
        else:
            dphi_dt = np.gradient(varphi, trunc)

        # --- step 2: φ → h ---
        if self.inversion == "calibration":
            cal = calibration if calibration is not None else self.calibration
            if cal is None:
                raise ValueError(
                    "inversion='calibration' requires a CalibrationTable. "
                    "Pass calibration=... or set use_sg_filter=False."
                )
            phi_equiv = dphi_dt * self.tau  # rad/ns * ns → rad
            h_recon = cal.inverse(phi_equiv)

        elif self.inversion == "response":
            if self.qubit is None:
                raise ValueError(
                    "inversion='response' requires qubit=... "
                    "(TransmonQubit with EC, EJ, flux)."
                )
            # dφ/dt is the angular frequency shift Δω (natural units, ħ=1)
            h_recon = _qubit_inverse_frequency(dphi_dt, self.qubit)

        else:
            raise ValueError(f"Unknown inversion: {self.inversion}")

        return FluxSignal(type=8, t_list=trunc, signal=h_recon)