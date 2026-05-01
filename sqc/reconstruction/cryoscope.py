"""Cryoscope reconstruction helpers."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from sqc.reconstruction.base import Reconstruction


@dataclass
class CryoscopeReconstruction(Reconstruction):
    """Reconstruct a waveform from cryoscope phase measurements."""

    calibration: object
    tau: float
    method: Literal["calib_inverse", "SG_diff", "diff"] = "calib_inverse"

    def reconstruct(self, measurement, kernel=None, calibration=None, dt=None) -> np.ndarray:
        """Return reconstructed flux samples."""
        varphi = np.asarray(measurement.data["varphi"], dtype=float)
        trunc_list = np.asarray(measurement.axes["trunc"], dtype=float)
        if dt is None:
            dt = abs(trunc_list[0] - trunc_list[1])

        table = calibration or self.calibration
        if self.method == "calib_inverse":
            delta_phi = np.diff(varphi, prepend=0)
            delta_phi_norm = delta_phi * self.tau / dt
            return table.inverse(delta_phi_norm)

        if self.method == "SG_diff":
            from scipy.signal import savgol_filter

            d_varphi = savgol_filter(
                varphi, window_length=11, polyorder=3, deriv=1, delta=dt
            )
            return table.inverse(d_varphi)

        if self.method == "diff":
            d_varphi = np.gradient(varphi, dt)
            return table.inverse(d_varphi)

        raise ValueError(f"unknown method: {self.method}")


def build_calibration_curve(h_list, phi_list):
    """Return ``(phi_of_h, h_of_phi)`` interpolation callables."""
    from scipy.interpolate import interp1d

    phi_of_h = interp1d(h_list, phi_list, kind="cubic", fill_value="extrapolate")
    idx = np.argsort(phi_list)
    h_sorted = np.asarray(h_list)[idx]
    phi_sorted = np.asarray(phi_list)[idx]
    mask = np.diff(phi_sorted, prepend=-np.inf) > 1e-12
    h_of_phi = interp1d(
        phi_sorted[mask], h_sorted[mask], kind="cubic", fill_value="extrapolate"
    )
    return phi_of_h, h_of_phi
