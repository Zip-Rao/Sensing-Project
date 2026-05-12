"""sqc.reconstruction.hammerstein — Hammerstein-Wiener nonlinear reconstruction.

Port from src/analysis.py: Analysis.hammerstein_wiener_deconvolution.

Block structure: B(t) → ω(B) → φ → Δp

Step 1: linear Wiener deconvolution to obtain ω(t) = Δω(t).
Step 2: inverse Transmon dispersion relation to obtain B(t).

The inverse dispersion relation is:
    B = (1/π) arccos( (ω + ω₀ + EC)² / (8 EC EJ) )

derived from the Transmon frequency formula (Gao 2021 Eq. 18):
    ω_T = sqrt(8 EJ(Φ) EC) - EC
    with EJ(Φ) = EJ_0 |cos(π B)|
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.config import CONFIG
from sqc.reconstruction.base import Reconstruction
from sqc.reconstruction.wiener import WienerReconstruction


def _get_qubit_params(qubit) -> tuple[float, float, float]:
    """Extract (EC, EJ, frequency) from QubitSpec or legacy TransmonQubit."""
    if hasattr(qubit, "EC"):
        EC = qubit.EC
    else:
        raise TypeError(f"qubit {type(qubit).__name__} has no EC attribute")

    if hasattr(qubit, "EJ"):
        EJ = qubit.EJ
    elif hasattr(qubit, "EJ_0"):
        EJ = qubit.EJ_0
    else:
        raise TypeError(f"qubit {type(qubit).__name__} has no EJ/EJ_0 attribute")

    if hasattr(qubit, "frequency"):
        freq = qubit.frequency() if callable(qubit.frequency) else qubit.frequency
    else:
        freq = 0.0

    return EC, EJ, freq


@dataclass
class HammersteinWienerReconstruction(Reconstruction):
    """Hammerstein-Wiener nonlinear reconstruction.

    Block-structured B-field reconstruction that accounts for the
    nonlinear Transmon dispersion relation:
        B(t) → ω(B) → φ → Δp

    Replaces Analysis.hammerstein_wiener_deconvolution.

    Parameters
    ----------
    qubit
        QubitSpec or legacy TransmonQubit providing EC, EJ, frequency.
    lambda_reg : float
        Regularisation parameter for the inner Wiener filter.
        Default: ``CONFIG.reconstruction.lambda_reg``.
    """

    qubit: object
    lambda_reg: float = field(
        default_factory=lambda: CONFIG.reconstruction.lambda_reg
    )

    def reconstruct(
        self,
        measurement,
        kernel,
        calibration=None,
        dt: float | None = None,
    ) -> np.ndarray:
        """Block-structured B(t) reconstruction.

        Parameters
        ----------
        measurement : ExperimentResult or np.ndarray
            If ExperimentResult: extracts "delta_p" from measurement.data.
            If ndarray: treated directly as delta_p.
        kernel : np.ndarray
            Control kernel.
        calibration : ignored
        dt : float or None
            Time step. If None, inferred from measurement.

        Returns
        -------
        np.ndarray
            Reconstructed magnetic field B.
        """
        # Step 1: linear Wiener gives ω(t)
        wiener = WienerReconstruction(lambda_reg=self.lambda_reg)
        omega_signal = wiener.reconstruct(measurement, kernel, dt=dt)
        omega = np.asarray(omega_signal.signal, dtype=float)

        # Step 2: inverse Transmon dispersion
        # ω = sqrt(8 EC EJ(B)) - EC  ⇒  EJ(B) = (ω + EC)^2 / (8 EC)
        # EJ(B) = EJ_0 |cos(π B)|  ⇒  cos(π B) = EJ(B) / EJ_0
        # B = (1/π) arccos( (ω + EC)^2 / (8 EC EJ_0) )
        EC, EJ, freq = _get_qubit_params(self.qubit)

        B = (1.0 / np.pi) * np.arccos(
            (omega + freq + EC) ** 2 / (8 * EC * EJ)
        )
        return B
