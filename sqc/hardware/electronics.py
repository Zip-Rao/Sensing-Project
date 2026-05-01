"""sqc.hardware.electronics — AWG/ADC abstractions.

Full implementation in Phase 5.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class AWG:
    """Abstract AWG (arbitrary waveform generator).

    Full implementation in Phase 5.
    """

    name: str = ""
    sample_rate: float = 1.0  # GS/s
    resolution_bits: int = 14

    def generate(self, waveform, **kwargs):
        raise NotImplementedError("implemented in Phase 5")


@dataclass
class ADC:
    """Abstract ADC (analog-to-digital converter).

    Full implementation in Phase 5.
    """

    name: str = ""
    sample_rate: float = 1.0  # GS/s
    resolution_bits: int = 12

    def acquire(self, duration: float, **kwargs):
        raise NotImplementedError("implemented in Phase 5")
