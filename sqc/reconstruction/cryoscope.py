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
from sqc.control.flux_signal import FluxSignal
from sqc.reconstruction.base import Reconstruction


@dataclass
class CryoscopeReconstruction(Reconstruction):
    """Cryoscope waveform reconstruction.

    Replaces Analysis.get_signal_from_cryoscope + get_h_from_phi.

    Three reconstruction strategies:

    - ``"calib_inverse"``: discrete phase differences, scaled by tau/dt,
      mapped through calibration inverse. Best for precise calibration.
    - ``"SG_diff"``: Savitzky-Golay derivative → dφ/dt → phase_rate·tau
      → calibration inverse. Noise-robust.
    - ``"diff"``: np.gradient derivative → same pipeline. Simplest.

    Parameters
    ----------
    calibration : CalibrationTable
        Phase-to-flux calibration table (kind="phi_h") with inverse().
    tau : float
        Calibration square-pulse length (ns).
    method : str
        Reconstruction method. Default "calib_inverse".
    sg_window : int
        Savitzky-Golay window length (odd). Only for method="SG_diff".
    sg_poly : int
        Savitzky-Golay polynomial order. Only for method="SG_diff".
    """

    calibration: CalibrationTable
    tau: float
    method: Literal["calib_inverse", "SG_diff", "diff"] = "calib_inverse"
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

        Physics: dφ/dt = 2π·Δf(h(t)).  Calibration gives φ_cal(h) for a
        fixed tau, so h(t) = φ_cal⁻¹(dφ/dt · tau).

        Parameters
        ----------
        measurement : ExperimentResult
            Must contain data["varphi"] and axes["trunc"].
        kernel : optional
            Not used in cryoscope reconstruction.
        calibration : CalibrationTable, optional
            Override the instance calibration table.
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

        cal = calibration if calibration is not None else self.calibration

        # Compute dφ/dt (rad/ns), then scale to calibration-equivalent phase
        if self.method == "calib_inverse":
            # Discrete phase difference between consecutive truncations
            delta_phi = np.diff(varphi, prepend=0.0)
            phi_equiv = delta_phi * self.tau / dt  # → rad (as if held for tau)
        elif self.method == "SG_diff":
            from scipy.signal import savgol_filter

            window = self.sg_window
            if window >= len(varphi):
                window = max(3, len(varphi) // 2 * 2 - 1)
            dphi_dt = savgol_filter(
                varphi, window_length=window, polyorder=self.sg_poly,
                deriv=1, delta=dt,
            )
            phi_equiv = dphi_dt * self.tau  # rad/ns * ns → rad
        elif self.method == "diff":
            dphi_dt = np.gradient(varphi, trunc)
            phi_equiv = dphi_dt * self.tau  # rad/ns * ns → rad
        else:
            raise ValueError(f"Unknown method: {self.method}")

        h_recon = cal.inverse(phi_equiv)

        return FluxSignal(type=8, t_list=trunc, signal=h_recon)
