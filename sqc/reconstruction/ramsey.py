"""sqc.reconstruction.ramsey — Ramsey interferometry reconstruction.

    method="iq":      IQ dual-channel arctan2 phase extraction (needs p_e_I, p_e_Q)
    method="unwrap":  single-channel arccos + k-span phase unwrapping (needs p_e)
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

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
class RamseyReconstruction(Reconstruction):
    """Ramsey interferometry reconstruction.

    Parameters
    ----------
    qubit : QubitSpec or TransmonQubit
    method : str
        - ``"iq"``: IQ dual-channel arctan2 phase extraction.
        - ``"unwrap"``: single-channel arccos + k-span phase unwrapping.
    k_span : int
        Branch-search window for unwrap method. Default 3.
    """

    qubit: object
    method: Literal["iq", "unwrap"] = "unwrap"
    k_span: int = 3

    def reconstruct(self, measurement, **kwargs) -> np.ndarray:
        """Reconstruct B(τ) from Ramsey experiment data."""
        match self.method:
            case "iq":
                return self._reconstruct_iq(measurement)
            case "unwrap":
                return self._reconstruct_unwrap(measurement)
            case _:
                raise ValueError(f"Unknown method: {self.method}")

    # -- IQ ----------------------------------------------------------------

    def _reconstruct_iq(self, measurement) -> np.ndarray:
        p_e_I = np.asarray(measurement.data["p_e_I"], dtype=float)
        p_e_Q = np.asarray(measurement.data["p_e_Q"], dtype=float)
        tau_list = np.asarray(measurement.axes["tau"], dtype=float)

        C_I = np.max(p_e_I) - np.min(p_e_I)
        C_Q = np.max(p_e_Q) - np.min(p_e_Q)
        cosphi = (1 - 2 * p_e_I) / C_I
        sinphi = (1 - 2 * p_e_Q) / C_Q

        phi = np.unwrap(np.arctan2(sinphi, cosphi))
        B = np.gradient(phi, tau_list) / _get_sensitivity(self.qubit)
        return B

    # -- unwrap ------------------------------------------------------------

    def _reconstruct_unwrap(self, measurement) -> np.ndarray:
        p_e = np.asarray(measurement.data["p_e"], dtype=float)
        tau = np.asarray(measurement.axes["tau"], dtype=float)
        N = len(p_e)

        cos_phi = 1 - 2 * p_e
        theta = np.arccos(np.clip(cos_phi, -1, 1))
        ks = np.arange(-self.k_span, self.k_span + 1)

        def _candidates(i):
            return np.array(
                [2 * k * np.pi + s * theta[i] for k in ks for s in [1, -1]]
            )

        phi = np.zeros(N)
        cands = _candidates(0)
        phi[0] = cands[np.argmin(np.abs(cands))]
        cands = _candidates(1)
        phi[1] = cands[np.argmin(np.abs(cands - phi[0]))]

        for i in range(2, N):
            cands = _candidates(i)
            expected = 2 * phi[i - 1] - phi[i - 2]
            phi[i] = cands[np.argmin(np.abs(cands - expected))]

        B = np.gradient(phi, tau) / _get_sensitivity(self.qubit)
        return B
