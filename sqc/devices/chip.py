"""Chip topology containers and coupled-system adapters."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from qutip import Qobj, destroy, qeye, tensor

from sqc.devices.transmon import Coupled_System as _LegacyCoupledSystem
from sqc.devices.transmon import TransmonQubit

if TYPE_CHECKING:
    from sqc.hardware.control_line import ControlLine
    from sqc.hardware.transfer_matrix import TransferMatrix


@dataclass
class ChipTopology:
    """Multi-device chip topology for qubits, resonators, and Z lines."""

    qubits: list[TransmonQubit]
    resonators: list[object] = field(default_factory=list)
    couplings: dict[tuple[str, str], float] = field(default_factory=dict)
    control_lines: dict[str, "ControlLine"] = field(default_factory=dict)
    transfer_matrix: "TransferMatrix | None" = None

    def __post_init__(self) -> None:
        names = [q.name for q in self.qubits]
        if len(set(names)) != len(names):
            raise ValueError("ChipTopology requires unique qubit names")

    def hilbert_dim(self) -> int:
        """Return the total Hilbert-space dimension."""
        dim = 1
        for qubit in self.qubits:
            dim *= int(qubit.n_levels)
        for resonator in self.resonators:
            for n_level in self._resonator_levels(resonator):
                dim *= int(n_level)
        return dim

    def lift_qubit_op(self, op: Qobj, qubit_name: str) -> Qobj:
        """Embed a single-qubit operator into the full chip Hilbert space."""
        factors: list[Qobj] = []
        found = False
        for qubit in self.qubits:
            if qubit.name == qubit_name:
                factors.append(op)
                found = True
            else:
                factors.append(qeye(qubit.n_levels))
        if not found:
            raise KeyError(f"unknown qubit: {qubit_name}")
        factors.extend(self._resonator_identity_factors())
        return tensor(*factors)

    def lift_resonator_mode_op(
        self,
        op: Qobj,
        resonator_name: str,
        mode_index: int = 0,
    ) -> Qobj:
        """Embed a single resonator-mode operator into the full chip space."""
        factors: list[Qobj] = [qeye(qubit.n_levels) for qubit in self.qubits]
        found = False
        for index, resonator in enumerate(self.resonators):
            name = self._resonator_name(index, resonator)
            levels = self._resonator_levels(resonator)
            for mode, n_level in enumerate(levels):
                if name == resonator_name and mode == mode_index:
                    factors.append(op)
                    found = True
                else:
                    factors.append(qeye(n_level))
        if not found:
            raise KeyError(f"unknown resonator mode: {resonator_name}[{mode_index}]")
        return tensor(*factors)

    def hamiltonian_static(self) -> Qobj:
        """Return the static Hamiltonian for the topology."""
        identity = tensor(*self._identity_factors())
        h_total = 0.0 * identity

        for qubit in self.qubits:
            h_total += self.lift_qubit_op(qubit.hamiltonian_static(), qubit.name)

        for r_index, resonator in enumerate(self.resonators):
            r_name = self._resonator_name(r_index, resonator)
            frequencies = np.asarray(getattr(resonator, "frequencies", []), dtype=float)
            for mode, n_level in enumerate(self._resonator_levels(resonator)):
                freq = float(frequencies[mode]) if mode < len(frequencies) else 0.0
                a_mode = destroy(n_level)
                h_total += freq * self.lift_resonator_mode_op(
                    a_mode.dag() * a_mode,
                    r_name,
                    mode,
                )

        for (qubit_name, resonator_name), g in self.couplings.items():
            qubit = self.qubit(qubit_name)
            a_qubit = self.lift_qubit_op(qubit.a, qubit_name)
            mode_index = 0
            r_name, parsed_mode = self._parse_resonator_key(resonator_name)
            if parsed_mode is not None:
                mode_index = parsed_mode
            resonator = self.resonator(r_name)
            n_level = self._resonator_levels(resonator)[mode_index]
            a_res = self.lift_resonator_mode_op(destroy(n_level), r_name, mode_index)
            h_total += float(g) * (a_qubit.dag() * a_res + a_qubit * a_res.dag())

        return h_total

    def collapse_operators(self) -> list[Qobj]:
        """Return lifted qubit collapse operators."""
        c_ops: list[Qobj] = []
        for qubit in self.qubits:
            for c_op in qubit.collapse_operators():
                c_ops.append(self.lift_qubit_op(c_op, qubit.name))
        return c_ops

    def qubit(self, name: str) -> TransmonQubit:
        """Return a qubit by name."""
        for qubit in self.qubits:
            if qubit.name == name:
                return qubit
        raise KeyError(f"unknown qubit: {name}")

    def resonator(self, name: str) -> object:
        """Return a resonator by name."""
        for index, resonator in enumerate(self.resonators):
            if self._resonator_name(index, resonator) == name:
                return resonator
        raise KeyError(f"unknown resonator: {name}")

    @classmethod
    def from_legacy_coupled_system(cls, coupled_system) -> "ChipTopology":
        """Build a topology wrapper around a legacy ``Coupled_System``."""
        qubits = [coupled_system.qubit1, coupled_system.qubit2]
        resonators = [coupled_system.cavity]
        couplings: dict[tuple[str, str], float] = {}
        if hasattr(coupled_system, "g"):
            for q_index, qubit in enumerate(qubits):
                for mode, g in enumerate(coupled_system.g[q_index]):
                    couplings[(qubit.name, f"R0:{mode}")] = float(g)
        return cls(qubits=qubits, resonators=resonators, couplings=couplings)

    def _identity_factors(self) -> list[Qobj]:
        factors = [qeye(qubit.n_levels) for qubit in self.qubits]
        factors.extend(self._resonator_identity_factors())
        return factors

    def _resonator_identity_factors(self) -> list[Qobj]:
        factors: list[Qobj] = []
        for resonator in self.resonators:
            for n_level in self._resonator_levels(resonator):
                factors.append(qeye(n_level))
        return factors

    @staticmethod
    def _resonator_levels(resonator) -> list[int]:
        levels = getattr(resonator, "n_levels", None)
        if levels is None:
            levels = [getattr(resonator, "n_levels", 0)]
        if isinstance(levels, (int, np.integer)):
            return [int(levels)]
        return [int(level) for level in levels]

    @staticmethod
    def _resonator_name(index: int, resonator) -> str:
        return str(getattr(resonator, "name", f"R{index}"))

    @staticmethod
    def _parse_resonator_key(name: str) -> tuple[str, int | None]:
        if ":" not in name:
            return name, None
        base, mode = name.split(":", 1)
        return base, int(mode)


CoupledSystem = _LegacyCoupledSystem
Coupled_System = _LegacyCoupledSystem

__all__ = ["ChipTopology", "CoupledSystem", "Coupled_System"]
