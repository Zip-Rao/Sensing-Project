"""Calibration abstractions and data tables."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

import numpy as np


class Calibration(ABC):
    """Base class for calibration routines."""

    @abstractmethod
    def calibrate(self):
        """Run calibration and return calibration data."""


@dataclass(frozen=True)
class CalibrationTable:
    """Lookup table produced by a calibration routine."""

    qubit_name: str
    kind: str
    inputs: np.ndarray
    outputs: np.ndarray
    fit_params: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    def evaluate(self, x):
        """Evaluate by one-dimensional interpolation."""
        from scipy.interpolate import interp1d

        kind = "cubic" if len(self.inputs) >= 4 else "linear"
        interp = interp1d(
            self.inputs,
            self.outputs,
            kind=kind,
            fill_value="extrapolate",
            assume_sorted=False,
        )
        return interp(x)

    def inverse(self, y):
        """Evaluate an inverse lookup by sorted interpolation."""
        from scipy.interpolate import interp1d

        idx = np.argsort(self.outputs)
        x_sorted = np.asarray(self.inputs)[idx]
        y_sorted = np.asarray(self.outputs)[idx]
        mask = np.diff(y_sorted, prepend=-np.inf) > 1e-12
        if np.count_nonzero(mask) < 2:
            raise ValueError("CalibrationTable inverse requires monotonic outputs")
        kind = "cubic" if np.count_nonzero(mask) >= 4 else "linear"
        interp = interp1d(
            y_sorted[mask],
            x_sorted[mask],
            kind=kind,
            fill_value="extrapolate",
            assume_sorted=True,
        )
        return interp(y)
