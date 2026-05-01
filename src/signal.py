"""Compatibility mirror for legacy `src.signal` imports.

New code should use `sqc.control.flux_signal`.
"""
from __future__ import annotations

from sqc.control.flux_signal import CompositeSignal, FluxSignal as Signal


__all__ = ["Signal", "CompositeSignal"]
