"""Physical control-line data structures."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

from sqc.control.waveform import Waveform
from sqc.hardware.distortion import DistortionModel

if TYPE_CHECKING:
    from sqc.calibration.predistortion import PredistortionDesigner


@dataclass
class ControlLine:
    """Physical xy/z/readout control line with optional transfer behavior."""

    name: str
    kind: Literal["xy", "z", "readout"]
    source: str
    target: str
    transfer_function: DistortionModel | None = None
    metadata: dict = field(default_factory=dict)

    def apply(self, waveform: Waveform) -> Waveform:
        """Forward model an AWG waveform as it arrives on chip."""
        if self.transfer_function is None:
            return waveform.copy()
        return self.transfer_function.apply(waveform)

    def predistort(
        self,
        target_waveform: Waveform,
        designer: "PredistortionDesigner",
    ) -> Waveform:
        """Compute the AWG waveform needed for the target on-chip waveform."""
        if self.transfer_function is None:
            return target_waveform.copy()
        return designer.predistort(target_waveform, self.transfer_function)
