"""sqc.hardware — Physical hardware models (control lines, transfer matrices, etc.)."""
from .distortion import (
    DistortionModel,
    SingleExponentialDistortion,
    MultiExponentialDistortion,
    FIRDistortion,
    IIRDistortion,
    CustomTransferDistortion,
    CascadeDistortion,
)
from .control_line import ControlLine
from .transfer_matrix import TransferMatrix

__all__ = [
    "DistortionModel",
    "SingleExponentialDistortion",
    "MultiExponentialDistortion",
    "FIRDistortion",
    "IIRDistortion",
    "CustomTransferDistortion",
    "CascadeDistortion",
    "ControlLine",
    "TransferMatrix",
]
