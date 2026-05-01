"""sqc.experiments.base — Experiment ABC.

Experiment = device + flux + sequence + readout + runner.
"""
from __future__ import annotations

from abc import ABC, abstractmethod


class Experiment(ABC):
    """Abstract sensing experiment.

    Subclasses implement build_sequence() and run().
    An experiment composes a device, flux signal, control
    sequence, readout model, and simulation runner.

    Full implementation in Phase 2.
    """

    @abstractmethod
    def build_sequence(self) -> None:
        """Build the control sequence for this experiment."""
        raise NotImplementedError("implemented in Phase 2")

    @abstractmethod
    def run(self, *args, **kwargs):
        """Execute the experiment and return results."""
        raise NotImplementedError("implemented in Phase 2")
