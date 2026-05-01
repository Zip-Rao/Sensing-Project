"""Coupler abstractions."""
from __future__ import annotations

from abc import ABC, abstractmethod


class TunableCoupler(ABC):
    """Base class for tunable coupling elements."""

    @abstractmethod
    def coupling_strength(self, *args, **kwargs) -> float:
        """Return the effective coupling strength."""
