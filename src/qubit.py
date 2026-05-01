"""Compatibility mirror for legacy `src.qubit` imports.

New code should import from `sqc.devices.*` and `sqc.control.gates`.
This module keeps historical notebooks and web_demo.py working.
"""
from __future__ import annotations

from sqc.control.gates import (
    ideal_CZ,
    ideal_iSWAP,
    simulate_CZ,
    simulate_iSWAP_simple as simulate_iSWAP,
)
from sqc.devices.chip import CoupledSystem as Coupled_System
from sqc.devices.resonator import Resonator as Cavity
from sqc.devices.transmon import QubitSpec, TransmonQubit


__all__ = [
    "TransmonQubit",
    "QubitSpec",
    "Cavity",
    "Coupled_System",
    "ideal_iSWAP",
    "simulate_iSWAP",
    "ideal_CZ",
    "simulate_CZ",
]
