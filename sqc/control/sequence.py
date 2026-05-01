"""Pulse sequence factories."""
from __future__ import annotations

from dataclasses import dataclass, field

from sqc.control.pulse import (
    CompositePulse,
    Pulse,
    create_cpmg_pulse,
    create_cryoscope_pulse,
    create_diff_echo_pulse,
    create_echo_pulse,
    create_pulse,
    create_ramsey_pulse,
)


@dataclass
class PulseSequence:
    """Container for an ordered list of pulse objects."""

    pulses: list[Pulse] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)


__all__ = [
    "Pulse",
    "CompositePulse",
    "PulseSequence",
    "create_pulse",
    "create_ramsey_pulse",
    "create_diff_echo_pulse",
    "create_echo_pulse",
    "create_cpmg_pulse",
    "create_cryoscope_pulse",
]
