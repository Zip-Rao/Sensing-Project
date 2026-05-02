"""src_mirror.distortion — Re-export from sqc.hardware.distortion.

After P4 internalization, the canonical implementation lives in sqc/.
This file maintains backward-compatibility for any code importing
from src_mirror.distortion.

These were originally developed in src_mirror/ per R1 (never modify src/).
P4 internalized them into sqc/hardware/distortion.py.

Physical reference: Gao 2021 §III.D (control-line transfer functions),
§V.E (Cryoscope calibration of distortion tails).
"""
from __future__ import annotations

from sqc.hardware.distortion import (
    DistortionModel,
    SingleExponentialDistortion,
    MultiExponentialDistortion,
    FIRDistortion,
    IIRDistortion,
    CustomTransferDistortion,
)

__all__ = [
    "DistortionModel",
    "SingleExponentialDistortion",
    "MultiExponentialDistortion",
    "FIRDistortion",
    "IIRDistortion",
    "CustomTransferDistortion",
]
