"""src_mirror.pulse — Mirror of src/pulse.py using sqc.control.pulse + sqc.control.sequence.

New code should use sqc.control.* directly.
This file exists for backward-compatible imports from src_mirror.pulse.

See idea/refactor/_refactor_plan.md §8.
"""
from sqc.control.pulse import Pulse, CompositePulse
from sqc.control.sequence import (
    create_pulse,
    create_ramsey_pulse,
    create_diff_echo_pulse,
    create_echo_pulse,
    create_cpmg_pulse,
    create_cryoscope_pulse,
)

__all__ = [
    "Pulse",
    "CompositePulse",
    "create_pulse",
    "create_ramsey_pulse",
    "create_diff_echo_pulse",
    "create_echo_pulse",
    "create_cpmg_pulse",
    "create_cryoscope_pulse",
]
