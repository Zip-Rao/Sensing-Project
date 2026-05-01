"""sqc.devices.base — Device ABC.

Every physical device (qubit, resonator, coupler, chip) inherits from Device.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

from qutip import Qobj


class Device(ABC):
    """Base class for all physical devices.

    Subclasses MUST be parameter-only objects: they describe physics but
    do not store experimental state. To represent a device under perturbation
    (e.g., qubit at non-zero flux), use a method that returns a NEW Device.

    Invariants
    ----------
    - hamiltonian_static() must return a time-independent H.
    - collapse_operators() must return valid Lindblad ops.
    - hilbert_dim() must return the dimension of the full Hilbert space.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable device identifier."""
        ...

    @abstractmethod
    def hilbert_dim(self) -> int:
        """Dimension of the full Hilbert space."""
        ...

    @abstractmethod
    def hamiltonian_static(self) -> Qobj:
        """Time-independent Hamiltonian without driving or external fields."""
        ...

    @abstractmethod
    def collapse_operators(self) -> list[Qobj]:
        """Lindblad collapse operators for dissipation and dephasing."""
        ...
