"""Simulation result data structures."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import pickle

import numpy as np


@dataclass
class ExperimentResult:
    """Container for experiment outputs and metadata."""

    data: dict[str, np.ndarray] = field(default_factory=dict)
    axes: dict[str, np.ndarray] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)
    config: dict = field(default_factory=dict)

    def save(self, path: str | Path) -> None:
        """Save this result as a pickle."""
        with open(path, "wb") as f:
            pickle.dump(self, f, protocol=5)

    @classmethod
    def load(cls, path: str | Path) -> "ExperimentResult":
        """Load an ExperimentResult pickle."""
        with open(path, "rb") as f:
            obj = pickle.load(f)
        if not isinstance(obj, cls):
            raise TypeError(f"Loaded {type(obj).__name__}, expected ExperimentResult")
        return obj


@dataclass(frozen=True)
class MeasurementTrace:
    """A one-dimensional measurement trace."""

    axis: np.ndarray
    p_e: np.ndarray
    p_e_iq: tuple[np.ndarray, np.ndarray] | None = None
    metadata: dict = field(default_factory=dict)


def extract_expectation(result, e_ops_index: int = 0) -> np.ndarray:
    """Extract expectation values from a QuTiP result object."""
    if hasattr(result, "expect"):
        if isinstance(result.expect, list) and len(result.expect) > e_ops_index:
            return np.asarray(result.expect[e_ops_index])
        return np.asarray(result.expect)
    raise ValueError("Result has no .expect attribute")


def extract_population(result, level: int) -> np.ndarray:
    """Extract population of a Fock level from result.states."""
    from qutip import basis, expect

    if not hasattr(result, "states"):
        raise ValueError("Result has no .states attribute")
    populations = []
    for state in result.states:
        n_levels = state.dims[0][0]
        proj = basis(n_levels, level) * basis(n_levels, level).dag()
        populations.append(expect(proj, state))
    return np.asarray(populations)
