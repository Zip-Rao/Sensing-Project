"""src_mirror.signal — Mirror of src/signal.py using sqc.control.flux_signal.

New code should use sqc.control.* directly.
This file exists for backward-compatible imports from src_mirror.signal.

See idea/refactor/_refactor_plan.md §8.
"""
from sqc.control.flux_signal import FluxSignal as Signal, CompositeSignal

__all__ = ["Signal", "CompositeSignal"]
