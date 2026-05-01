"""sqc.workflows.base — Workflow ABC.

High-level workflows that compose experiments, calibrations,
and reconstructions.
"""
from __future__ import annotations

from abc import ABC, abstractmethod


class Workflow(ABC):
    """Abstract high-level experimental workflow.

    A workflow orchestrates multiple experiments, calibrations,
    and reconstructions to achieve a research goal.

    Full implementation in Phase 5.
    """

    @abstractmethod
    def run(self) -> dict:
        """Execute the workflow and return named results."""
        raise NotImplementedError("implemented in Phase 5")
