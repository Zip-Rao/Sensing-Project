"""sqc.reconstruction.base — Reconstruction ABC.

Abstract interface for waveform/field reconstruction from
experimental measurement data.
"""
from __future__ import annotations

from abc import ABC, abstractmethod


class Reconstruction(ABC):
    """Abstract waveform reconstruction.

    Subclasses implement reconstruct() which takes experimental
    measurements and returns the estimated physical waveform.

    Invariants
    ----------
    - reconstruct() must be pure: it reads measurement data but does
      NOT run new simulations. If forward simulation is needed (e.g., LM),
      inject the qubit spec and control pulse as constructor args.
    - Output is always a FluxSignal (or Waveform) in physical units.

    Full implementation in Phase 2-3.
    """

    @abstractmethod
    def reconstruct(self, *args, **kwargs):
        """Reconstruct the physical signal from measurement data."""
        raise NotImplementedError("implemented in Phase 2/3")
