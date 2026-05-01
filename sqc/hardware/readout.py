"""sqc.hardware.readout — ReadoutModel ABC.

Models IQ-demodulated dispersive readout.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class ReadoutModel(ABC):
    """Abstract readout model.

    Subclasses implement IQ-demodulated measurement of qubit state
    via dispersive interaction with a readout resonator.

    Full implementation in Phase 2.
    """

    @abstractmethod
    def measure(
        self, state, t_list: np.ndarray, **kwargs
    ) -> tuple[np.ndarray, np.ndarray]:
        """Perform IQ readout and return (I, Q) projections.

        Parameters
        ----------
        state : Qobj
            Quantum state to measure.
        t_list : np.ndarray
            Readout time axis.

        Returns
        -------
        I : np.ndarray
            In-phase component.
        Q : np.ndarray
            Quadrature component.
        """
        raise NotImplementedError("implemented in Phase 2")
