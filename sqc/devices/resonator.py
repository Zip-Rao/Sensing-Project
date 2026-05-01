"""sqc.devices.resonator — Multi-mode resonator model.

Verbatim port of src/qubit.py:Cavity → Resonator.
Renamed from Cavity to Resonator per _refactor_plan.md §7.3.
"""
from __future__ import annotations

import numpy as np
from qutip import Qobj, destroy, expect, qeye, tensor

from .base import Device


class Resonator(Device):
    """Multi-mode readout/coupling resonator.

    Parameters
    ----------
    frequencies : list[float]
        Mode frequencies (rad*GHz).
    n_levels : list[int]
        Fock truncation for each mode.
    """

    def __init__(
        self, frequencies: list[float], n_levels: list[int]
    ) -> None:
        self.frequencies: list[float] = list(frequencies)
        self.M: int = len(frequencies)
        self.n_levels: list[int] = list(n_levels)

        # Build mode operators
        self.a: list[Qobj] = []
        self.adag: list[Qobj] = []
        self.n: list[Qobj] = []
        self.I: Qobj = tensor([qeye(n) for n in self.n_levels])

        for M_idx in range(self.M):
            ops = []
            for m in range(self.M):
                if m == M_idx:
                    ops.append(destroy(self.n_levels[m]))
                else:
                    ops.append(qeye(self.n_levels[m]))
            a_M = tensor(ops)
            self.a.append(a_M)
            self.adag.append(a_M.dag())
            self.n.append(a_M.dag() * a_M)

        # Hamiltonian
        self.hamiltonian: Qobj = self.get_hamiltonian()

    # -- Device ABC ---------------------------------------------------------

    @property
    def name(self) -> str:
        return f"Resonator({self.M}modes)"

    def hilbert_dim(self) -> int:
        return int(np.prod(self.n_levels))

    def hamiltonian_static(self) -> Qobj:
        return self.hamiltonian

    def collapse_operators(self) -> list[Qobj]:
        # Not implemented yet for resonator
        return []

    # -- Hamiltonian --------------------------------------------------------

    def get_hamiltonian(self) -> Qobj:
        """Linear Hamiltonian: sum(f_i * n_i)."""
        H = 0 * self.I
        for i in range(self.M):
            H += self.frequencies[i] * self.n[i]
        return H

    def get_hamiltonian_rot(self, omega_d: float) -> Qobj:
        """Rotating-frame Hamiltonian.

        Parameters
        ----------
        omega_d : float
            Reference frequency (rad*GHz).

        Returns
        -------
        Qobj
        """
        H = 0 * self.I
        for i in range(self.M):
            H += (self.frequencies[i] - omega_d) * self.n[i]
        return H
