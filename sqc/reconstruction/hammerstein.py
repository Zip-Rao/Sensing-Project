"""Hammerstein-Wiener reconstruction."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sqc.reconstruction.base import Reconstruction
from sqc.reconstruction.wiener import WienerReconstruction


@dataclass
class HammersteinWienerReconstruction(Reconstruction):
    """Recover flux after linear deconvolution of the frequency response."""

    qubit: object
    lambda_reg: float = 1.0

    def reconstruct(self, measurement, kernel, calibration=None, dt=None):
        """Return ``(t_list, B)`` using the legacy transmon inverse formula."""
        omega_signal = WienerReconstruction(self.lambda_reg).reconstruct(
            measurement, kernel, dt=dt
        )
        omega = omega_signal.signal
        ec = self.qubit.EC
        ej = self.qubit.EJ
        frequency = self.qubit.frequency() if callable(getattr(self.qubit, "frequency", None)) else self.qubit.frequency
        ratio = (omega + frequency + ec) ** 2 / (8 * ec * ej)
        ratio = np.clip(ratio, -1.0, 1.0)
        b_field = (1 / np.pi) * np.arccos(ratio)
        return omega_signal.t_list, b_field
