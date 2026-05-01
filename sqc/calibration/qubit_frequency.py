"""Qubit-frequency calibration routines."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from sqc.calibration.base import Calibration, CalibrationTable


@dataclass
class QubitFrequencyCalibration(Calibration):
    """Calibrate f01 using the current transmon model."""

    qubit: object

    def calibrate(self) -> CalibrationTable:
        """Return a one-point frequency table for the current flux."""
        return CalibrationTable(
            qubit_name=self.qubit.spec().name,
            kind="f01",
            inputs=np.asarray([self.qubit.flux]),
            outputs=np.asarray([self.qubit.frequency]),
            fit_params={"method": "static_model"},
            metadata={},
        )


@dataclass
class TransientFrequencyCalibration(Calibration):
    """Transient-protocol frequency calibration interface."""

    qubit: object
    polynomial_order: int = 3
    test_signal_kind: Literal["ramp", "sine", "gaussian"] = "ramp"

    def calibrate(self) -> CalibrationTable:
        """Run transient frequency calibration once Track B case 8 exists."""
        raise NotImplementedError(
            "TransientFrequencyCalibration awaits Track B case 8 and "
            "calibrate_frequency_response implementation."
        )
