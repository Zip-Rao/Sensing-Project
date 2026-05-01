"""Experiment abstractions."""
from __future__ import annotations

from abc import ABC, abstractmethod


class Experiment(ABC):
    """Base class for experimental protocols."""

    @abstractmethod
    def build_sequence(self):
        """Build the control sequence."""

    @abstractmethod
    def run(self):
        """Execute the experiment."""
