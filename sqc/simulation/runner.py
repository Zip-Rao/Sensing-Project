"""sqc.simulation.runner — Simulation runner ABCs.

RunnerBase: abstract interface for time evolution.
MesolveRunner: QuTiP mesolve-based solver.
SlidingMeasurementRunner: sliding measurement protocol (Phase 2).

See _refactor_plan.md §6.5.
"""
from __future__ import annotations

from abc import ABC, abstractmethod


class RunnerBase(ABC):
    """Abstract base for all simulation runners.

    Subclasses implement run() to perform time evolution
    and return an ExperimentResult.
    """

    @abstractmethod
    def run(self, *args, **kwargs):
        """Execute the simulation and return results."""
        raise NotImplementedError("implemented in Phase 2")


class MesolveRunner(RunnerBase):
    """QuTiP mesolve-based simulation runner.

    Full implementation in Phase 2.
    """

    def run(self, *args, **kwargs):
        raise NotImplementedError("implemented in Phase 2")


class SlidingMeasurementRunner(RunnerBase):
    """Sliding measurement protocol runner.

    Full implementation in Phase 2.
    """

    def run(self, *args, **kwargs):
        raise NotImplementedError("implemented in Phase 2")
