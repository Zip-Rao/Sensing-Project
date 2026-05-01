"""Workflow abstractions."""
from __future__ import annotations

from abc import ABC, abstractmethod


class Workflow(ABC):
    """Base class for high-level workflows."""

    @abstractmethod
    def run(self):
        """Run the workflow."""
