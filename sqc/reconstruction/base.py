"""Waveform reconstruction abstractions."""
from __future__ import annotations

from abc import ABC, abstractmethod


class Reconstruction(ABC):
    """Base class for reconstruction algorithms."""

    @abstractmethod
    def reconstruct(self, *args, **kwargs):
        """Reconstruct a signal from measurement data."""
