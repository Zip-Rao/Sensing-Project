"""Abstract base classes for physical devices."""
from __future__ import annotations

from abc import ABC, abstractmethod

from qutip import Qobj


class Device(ABC):
    """Base class for physical devices.

    Device subclasses describe physical parameters and operators. New code
    should keep experiment-specific state outside device parameter specs.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Device identifier."""

    @abstractmethod
    def hilbert_dim(self) -> int:
        """Return the Hilbert-space dimension."""

    @abstractmethod
    def hamiltonian_static(self) -> Qobj:
        """Return the time-independent Hamiltonian."""

    @abstractmethod
    def collapse_operators(self) -> list[Qobj]:
        """Return Lindblad collapse operators."""
