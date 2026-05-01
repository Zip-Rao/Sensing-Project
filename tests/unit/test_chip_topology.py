"""Unit tests for P5 chip topology containers."""
from __future__ import annotations

import numpy as np
import pytest
from qutip import destroy, qeye

from sqc.devices.chip import ChipTopology
from sqc.devices.transmon import Cavity, Coupled_System, TransmonQubit


pytestmark = pytest.mark.unit


def _make_qubit(name: str) -> TransmonQubit:
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        n_levels=2,
        name=name,
    )


def test_chip_topology_hilbert_dim_and_lift_qubit_op():
    qa = _make_qubit("QA")
    qb = _make_qubit("QB")
    chip = ChipTopology(qubits=[qa, qb])

    lifted = chip.lift_qubit_op(qa.n, "QA")

    assert chip.hilbert_dim() == 4
    assert lifted.shape == (4, 4)
    assert lifted == qa.n & qeye(qb.n_levels)


def test_chip_topology_static_hamiltonian_dimension_with_resonator():
    qa = _make_qubit("QA")
    qb = _make_qubit("QB")
    resonator = Cavity(frequencies=[2 * np.pi * 6.0], n_levels=[2])
    chip = ChipTopology(
        qubits=[qa, qb],
        resonators=[resonator],
        couplings={("QA", "R0"): 2 * np.pi * 0.01},
    )

    h_static = chip.hamiltonian_static()

    assert chip.hilbert_dim() == 8
    assert h_static.shape == (8, 8)
    assert len(chip.collapse_operators()) == 4


def test_chip_topology_rejects_duplicate_qubit_names():
    with pytest.raises(ValueError):
        ChipTopology(qubits=[_make_qubit("Q"), _make_qubit("Q")])


def test_chip_topology_from_legacy_coupled_system_preserves_dimensions():
    qa = _make_qubit("QA")
    qb = _make_qubit("QB")
    resonator = Cavity(frequencies=[2 * np.pi * 6.0], n_levels=[2])
    legacy = Coupled_System(qa, qb, resonator)

    chip = ChipTopology.from_legacy_coupled_system(legacy)

    assert chip.hilbert_dim() == 8
    assert ("QA", "R0:0") in chip.couplings
