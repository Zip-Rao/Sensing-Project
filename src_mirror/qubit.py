"""src_mirror.qubit — Mirror of src/qubit.py using sqc.devices.transmon.

New code should use sqc.devices.* directly.
This file exists for backward-compatible imports from src_mirror.qubit.

See idea/refactor/_refactor_plan.md §8.
"""
from sqc.devices.transmon import TransmonQubit, QubitSpec
from sqc.devices.resonator import Resonator as Cavity
from sqc.devices.chip import CoupledSystem as Coupled_System
from sqc.control.gates import (
    ideal_iSWAP,
    simulate_iSWAP as simulate_iSWAP_simple,
    ideal_CZ,
    simulate_CZ,
)

# Provide module-level simulate_iSWAP with same name
simulate_iSWAP = simulate_iSWAP_simple

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
