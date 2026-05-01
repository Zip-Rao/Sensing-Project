"""sqc.hardware.transfer_matrix — TransferMatrix dataclass.

Models crosstalk between multiple control lines.
Full implementation in Phase 5.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class TransferMatrix:
    """Linear crosstalk matrix between N control lines.

    Phi_j(omega) = sum_i H_{ji}(omega) * V_i(omega)

    Attributes
    ----------
    n_lines : int
        Number of control lines.
    elements : dict[tuple[int, int], float]
        (i, j) -> H_{ji} DC value.
    frequency : np.ndarray or None
        Frequency axis (GHz). If None, DC matrix only.
    """

    n_lines: int
    elements: dict[tuple[int, int], float] = field(default_factory=dict)
    frequency: np.ndarray | None = None

    def apply(self, voltages: np.ndarray) -> np.ndarray:
        """Apply the transfer matrix to a voltage vector.

        Parameters
        ----------
        voltages : np.ndarray
            Shape (n_lines,) input voltages.

        Returns
        -------
        np.ndarray
            Shape (n_lines,) output fluxes.
        """
        raise NotImplementedError("implemented in Phase 5")
