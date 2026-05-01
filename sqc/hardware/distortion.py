"""sqc.hardware.distortion — DistortionModel ABC.

Models the transfer function from AWG voltage to on-chip flux.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class DistortionModel(ABC):
    """Abstract model of control-line distortion.

    Subclasses implement step_response, frequency_response, etc.
    to model linear or nonlinear signal distortion.

    Full implementation in Phase 4.
    """

    @abstractmethod
    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        """Apply distortion to a waveform.

        Parameters
        ----------
        waveform : np.ndarray
            Input samples.
        dt : float
            Sample spacing (ns).

        Returns
        -------
        np.ndarray
            Distorted waveform.
        """
        raise NotImplementedError("implemented in Phase 4")

    @abstractmethod
    def step_response(self, t: np.ndarray) -> np.ndarray:
        """Step response of the distortion channel."""
        raise NotImplementedError("implemented in Phase 4")

    @abstractmethod
    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        """Impulse response of the distortion channel."""
        raise NotImplementedError("implemented in Phase 4")

    @abstractmethod
    def frequency_response(self, f: np.ndarray) -> np.ndarray:
        """Frequency response (complex) of the distortion channel."""
        raise NotImplementedError("implemented in Phase 4")
