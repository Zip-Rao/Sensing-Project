"""sqc.calibration.base — Calibration ABC + CalibrationTable.

Calibration: abstract interface for qubit/system calibration workflows.
CalibrationTable: dataclass holding calibration results.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class CalibrationTable:
    """Container for calibration results.

    Attributes
    ----------
    name : str
        Calibration name.
    parameters : dict[str, Any]
        Calibrated parameter values.
    uncertainties : dict[str, float]
        Uncertainty estimates.
    metadata : dict
        Additional metadata.
    """

    name: str
    parameters: dict[str, Any] = field(default_factory=dict)
    uncertainties: dict[str, float] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)


class Calibration(ABC):
    """Abstract base for qubit/system calibration.

    Subclasses implement calibrate() which returns a CalibrationTable.

    Full implementation in Phase 2-3.
    """

    @abstractmethod
    def calibrate(self) -> CalibrationTable:
        """Run the calibration workflow and return results."""
        raise NotImplementedError("implemented in Phase 2/3")
