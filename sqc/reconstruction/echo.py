"""sqc.reconstruction.echo — differential echo reconstruction.

    φ(τ) = arcsin(2·p_e - 1)
    B(τ) = -φ / (2·k·κ·t_int)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sqc.reconstruction.base import Reconstruction


def _get_sensitivity(qubit) -> float:
    """Extract frequency sensitivity kappa = dω/dΦ."""
    if hasattr(qubit, "flux_bias"):
        return qubit.sensitivity()
    if hasattr(qubit, "frequency_sensitivity") and hasattr(qubit, "flux"):
        return qubit.frequency_sensitivity(qubit.flux)
    raise TypeError(
        f"qubit {type(qubit).__name__} has neither flux_bias (QubitSpec) "
        f"nor frequency_sensitivity + flux (legacy TransmonQubit)"
    )


@dataclass
class EchoReconstruction(Reconstruction):
    """Differential echo waveform reconstruction.

    Parameters
    ----------
    qubit : QubitSpec or TransmonQubit
    t_int : float
        Integration time per echo block (ns).
    k : int
        Number of π-pulse pairs.
    """

    qubit: object
    t_int: float
    k: int

    def reconstruct(self, measurement, **kwargs) -> np.ndarray:
        """Reconstruct B from differential echo p_e data.

        Parameters
        ----------
        measurement : ExperimentResult
            Must have data["p_e"].

        Returns
        -------
        np.ndarray
            Reconstructed magnetic field B.
        """
        p_e = np.asarray(measurement.data["p_e"], dtype=float)
        varphi = np.arcsin(2 * p_e - 1)
        kappa = _get_sensitivity(self.qubit)
        return -varphi / (2 * self.k * kappa * self.t_int)
