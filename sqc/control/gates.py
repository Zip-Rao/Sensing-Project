"""Single- and two-qubit gate helpers."""
from __future__ import annotations

from sqc.devices.transmon import ideal_CZ, ideal_iSWAP, simulate_CZ, simulate_iSWAP


simulate_iSWAP_simple = simulate_iSWAP

__all__ = [
    "ideal_iSWAP",
    "simulate_iSWAP",
    "simulate_iSWAP_simple",
    "ideal_CZ",
    "simulate_CZ",
]
