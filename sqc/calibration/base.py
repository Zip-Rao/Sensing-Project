"""sqc.calibration.base — Calibration ABC + CalibrationTable.

Calibration: abstract interface for qubit/system calibration workflows.
CalibrationTable: dataclass holding calibration results with interpolation.

Interpolation methods (evaluate / inverse) added in P3a per
_handbook.md §3.7.1.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np


@dataclass
class CalibrationTable:
    """Container for calibration results.

    Supports both legacy parameter-dict style (P1) and table-lookup style
    with evaluate/inverse methods (P3a+).

    Attributes
    ----------
    name : str
        Calibration name.
    parameters : dict[str, Any]
        Calibrated parameter values (legacy).
    uncertainties : dict[str, float]
        Uncertainty estimates (legacy).
    metadata : dict
        Additional metadata.
    qubit_name : str
        Qubit identifier (P3a+ table style). Defaults to name.
    kind : str
        Calibration kind, e.g. "phi_h", "f_phi" (P3a+ table style).
    inputs : np.ndarray or None
        Independent variable (e.g. flux) for interpolation.
    outputs : np.ndarray or None
        Dependent variable (e.g. phase, frequency) for interpolation.
    fit_params : dict
        Fit/interpolation parameters (P3a+ table style).
    """

    name: str
    parameters: dict[str, Any] = field(default_factory=dict)
    uncertainties: dict[str, float] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    # -- P3a+ table-style fields --
    qubit_name: str = ""
    kind: str = ""
    inputs: Optional[np.ndarray] = None
    outputs: Optional[np.ndarray] = None
    fit_params: dict = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Interpolation methods (P3a, per handbook §3.7.1)
    # ------------------------------------------------------------------

    def evaluate(self, x: np.ndarray) -> np.ndarray:
        """Interpolate outputs at query points x.

        Uses cubic spline interpolation with extrapolation.

        Parameters
        ----------
        x : np.ndarray
            Query points (same units as self.inputs).

        Returns
        -------
        np.ndarray
            Interpolated outputs at x.

        Raises
        ------
        ValueError
            If self.inputs / self.outputs are not set.
        """
        if self.inputs is None or self.outputs is None:
            raise ValueError(
                "CalibrationTable.inputs/outputs not set; "
                "call a calibration workflow first."
            )
        from scipy.interpolate import interp1d

        inp = np.asarray(self.inputs)
        out = np.asarray(self.outputs)

        # Fall back to linear/quadratic if too few points for cubic
        n_pts = len(inp)
        if n_pts >= 4:
            kind = "cubic"
        elif n_pts >= 3:
            kind = "quadratic"
        else:
            kind = "linear"

        f = interp1d(inp, out, kind=kind, fill_value="extrapolate")
        return np.asarray(f(x))

    def inverse(self, y: np.ndarray) -> np.ndarray:
        """Inverse interpolation: find inputs such that outputs ≈ y.

        Builds the inverse by sorting outputs and taking the monotonic
        chunk. Uses cubic spline with extrapolation.

        Parameters
        ----------
        y : np.ndarray
            Query points (same units as self.outputs).

        Returns
        -------
        np.ndarray
            Estimated inputs at y.

        Raises
        ------
        ValueError
            If self.inputs / self.outputs are not set, or if no
            monotonic region can be found.
        """
        if self.inputs is None or self.outputs is None:
            raise ValueError(
                "CalibrationTable.inputs/outputs not set; "
                "call a calibration workflow first."
            )
        from scipy.interpolate import interp1d

        inp = np.asarray(self.inputs)
        out = np.asarray(self.outputs)

        # Sort by outputs for inverse
        idx = np.argsort(out)
        x_sorted = inp[idx]
        y_sorted = out[idx]

        # Restrict to strictly monotonic segment
        mask = np.diff(y_sorted, prepend=-np.inf) > 1e-12
        if not np.any(mask):
            raise ValueError(
                "CalibrationTable outputs are not monotonic; "
                "cannot build inverse."
            )

        y_mono = y_sorted[mask]
        x_mono = x_sorted[mask]
        n_pts = len(y_mono)
        if n_pts >= 4:
            kind = "cubic"
        elif n_pts >= 3:
            kind = "quadratic"
        else:
            kind = "linear"

        f_inv = interp1d(
            y_mono, x_mono,
            kind=kind, fill_value="extrapolate",
        )
        return np.asarray(f_inv(y))


class Calibration(ABC):
    """Abstract base for qubit/system calibration.

    Subclasses implement calibrate() which returns a CalibrationTable.

    Full implementation in Phase 2-3.
    """

    @abstractmethod
    def calibrate(self) -> CalibrationTable:
        """Run the calibration workflow and return results."""
        raise NotImplementedError("implemented in Phase 2/3")
