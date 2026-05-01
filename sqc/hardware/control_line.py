"""sqc.hardware.control_line — ControlLine dataclass.

Describes the physical properties of a flux/charge control line.
Full implementation in Phase 4.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class ControlLine:
    """Physical parameters of a control line.

    Attributes
    ----------
    name : str
        Identifier (e.g., "fast_z_1").
    impedance : float
        Characteristic impedance (Ohm).
    attenuation_db : float
        Total attenuation from room temp to chip (dB).
    delay : float
        Propagation delay (ns).
    filter_type : str or None
        Type of filter on the line (e.g., "lowpass", "bias_tee").
    cutoff_freq : float or None
        Filter cutoff frequency (GHz), if applicable.
    """

    name: str
    impedance: float = 50.0
    attenuation_db: float = 20.0
    delay: float = 0.0
    filter_type: Optional[str] = None
    cutoff_freq: Optional[float] = None

    def transfer_function(self, frequency: float) -> complex:
        """Compute the complex transfer function at a given frequency.

        Parameters
        ----------
        frequency : float
            Frequency (GHz).

        Returns
        -------
        complex
            Transfer function H(f).
        """
        raise NotImplementedError("implemented in Phase 4")
