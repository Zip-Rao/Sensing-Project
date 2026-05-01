"""sqc.simulation.result — Experiment result data structures.

ExperimentResult: serialisable container for simulation outputs.
MeasurementTrace: single measurement axis + populations.
Helper functions: extract_expectation, extract_population.

See _refactor_plan.md §5.5.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import pickle
from typing import Optional

import numpy as np


@dataclass
class ExperimentResult:
    """Container for experiment/simulation results.

    Attributes
    ----------
    data : dict[str, np.ndarray]
        Named data arrays (p_e, kernel, delta_p, etc.).
    axes : dict[str, np.ndarray]
        Named coordinate axes (tau_list, t_list, etc.).
    metadata : dict
        Free-form metadata (protocol type, qubit params, timestamp).
    config : dict
        Simulation configuration snapshot (EC, EJ, n_levels, etc.).
    """

    data: dict[str, np.ndarray] = field(default_factory=dict)
    axes: dict[str, np.ndarray] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)
    config: dict = field(default_factory=dict)

    def save(self, path: str | Path) -> None:
        """Pickle this result to a file.

        Parameters
        ----------
        path : str or Path
            Output file path.
        """
        with open(path, "wb") as f:
            pickle.dump(self, f, protocol=5)

    @classmethod
    def load(cls, path: str | Path) -> "ExperimentResult":
        """Load a pickled ExperimentResult.

        Parameters
        ----------
        path : str or Path
            Input file path.

        Returns
        -------
        ExperimentResult

        Raises
        ------
        TypeError
            If loaded object is not an ExperimentResult.
        """
        with open(path, "rb") as f:
            obj = pickle.load(f)
        if not isinstance(obj, cls):
            raise TypeError(
                f"Loaded object is {type(obj).__name__}, not ExperimentResult"
            )
        return obj


@dataclass(frozen=True)
class MeasurementTrace:
    """Single measurement trace (1-D sweep).

    Attributes
    ----------
    axis : np.ndarray
        Sweep axis (e.g., tau values).
    p_e : np.ndarray
        Excited-state population.
    p_e_iq : tuple[np.ndarray, np.ndarray] or None
        I/Q projections if IQ readout was used.
    metadata : dict
        Trace-specific metadata.
    """

    axis: np.ndarray
    p_e: np.ndarray
    p_e_iq: tuple[np.ndarray, np.ndarray] | None = None
    metadata: dict = field(default_factory=dict)


def extract_expectation(result, e_ops_index: int = 0) -> np.ndarray:
    """Extract expectation values from a QuTiP mesolve/sesolve result.

    Parameters
    ----------
    result : qutip.Result
        QuTiP solver result object.
    e_ops_index : int
        Index into result.expect list. Default 0.

    Returns
    -------
    np.ndarray
        1-D array of expectation values.

    Raises
    ------
    ValueError
        If result has no .expect attribute or index out of range.
    """
    if hasattr(result, "expect"):
        if isinstance(result.expect, list) and len(result.expect) > e_ops_index:
            return np.asarray(result.expect[e_ops_index])
        return np.asarray(result.expect)
    raise ValueError("Result has no .expect attribute")


def extract_population(result, level: int) -> np.ndarray:
    """Extract population of a Fock level from result.states.

    Parameters
    ----------
    result : qutip.Result
        QuTiP solver result (must have store_states=True).
    level : int
        Fock level to project onto.

    Returns
    -------
    np.ndarray
        1-D array of populations at each time step.

    Raises
    ------
    ValueError
        If result has no .states attribute.
    """
    from qutip import basis, expect

    if not hasattr(result, "states"):
        raise ValueError(
            "Result has no .states attribute (set store_states=True)"
        )
    populations = []
    for state in result.states:
        n_levels = state.dims[0][0]
        proj = basis(n_levels, level) * basis(n_levels, level).dag()
        populations.append(expect(proj, state))
    return np.asarray(populations)
