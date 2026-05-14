"""sqc.hardware.control_line — ControlLine dataclass.

Models a physical control line (xy / z / readout) with
optional transfer function distortion.

Per _refactor_plan.md §5.3 and phase_4_handbook.md §3.1.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional

import numpy as np

from sqc.control.waveform import Waveform


@dataclass
class ControlLine:
    """Models a physical control line (xy / z / readout) with
    optional transfer function distortion.

    Attributes
    ----------
    name : str
        Identifier (e.g., "Z0", "XY0").
    kind : Literal["xy", "z", "readout"]
        Type of control line.
    source : str
        Source identifier (e.g., "AWG0:CH1").
    target : str
        Target identifier (e.g., "Q0", "C01").
    transfer_function : DistortionModel or None
        Optional distortion model for this line.
    metadata : dict
        Arbitrary key-value metadata.

    # Backward-compatible fields (from P1 stub)
    impedance : float
        Characteristic impedance (Ohm). Default 50.0.
    attenuation_db : float
        Total attenuation from room temp to chip (dB). Default 20.0.
    delay : float
        Propagation delay (ns). Default 0.0.
    filter_type : str or None
        Type of filter on the line (e.g., "lowpass", "bias_tee").
    cutoff_freq : float or None
        Filter cutoff frequency (GHz), if applicable.
    """

    name: str
    kind: Literal["xy", "z", "readout"]
    source: str
    target: str
    transfer_function: Optional["DistortionModel"] = None
    metadata: dict = field(default_factory=dict)

    # backward-compatible physical parameters
    impedance: float = 50.0
    attenuation_db: float = 20.0
    delay: float = 0.0
    filter_type: Optional[str] = None
    cutoff_freq: Optional[float] = None

    def apply(self, awg_waveform: Waveform) -> Waveform:
        """Forward-model: AWG waveform -> on-chip waveform.

        If no transfer_function is set, returns awg_waveform unchanged
        (possibly with a delay if delay > 0).

        Parameters
        ----------
        awg_waveform : Waveform
            Input waveform from AWG.

        Returns
        -------
        Waveform
            Waveform after propagation through the control line.
        """
        if self.transfer_function is None:
            wf = awg_waveform.copy()
        else:
            wf = self.transfer_function.apply_to_waveform(awg_waveform)

        # Apply propagation delay if nonzero (simple time shift)
        if self.delay > 0:
            dt = float(awg_waveform.t_list[1] - awg_waveform.t_list[0])
            shift_samples = int(np.round(self.delay / dt))
            if shift_samples > 0 and shift_samples < len(wf.samples):
                shifted = np.zeros_like(wf.samples)
                shifted[shift_samples:] = wf.samples[:-shift_samples]
                wf = Waveform(
                    t_list=wf.t_list.copy(),
                    samples=shifted,
                    metadata={**wf.metadata, "delay_ns": self.delay},
                )

        return wf

    def predistort(
        self,
        target_waveform: Waveform,
        designer: "PredistortionDesigner",
    ) -> Waveform:
        """Inverse: target on-chip waveform -> required AWG waveform.

        Parameters
        ----------
        target_waveform : Waveform
            Desired on-chip waveform.
        designer : PredistortionDesigner
            Designer that computes the predistortion filter.

        Returns
        -------
        Waveform
            Predistorted AWG waveform.
        """
        return designer.predistort(target_waveform, self.transfer_function)
