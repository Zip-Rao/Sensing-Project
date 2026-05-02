"""tests.unit.test_chip_topology — Unit tests for ChipTopology."""
from __future__ import annotations

import numpy as np
import pytest
from qutip import Qobj, basis, destroy, qeye, tensor

from sqc.devices.chip import ChipTopology, CoupledSystem
from sqc.devices.transmon import TransmonQubit
from sqc.devices.resonator import Resonator


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def qubit_qa():
    return TransmonQubit(
        EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
        T1=10000, T2=8000, flux=0.0, n_levels=3, name="QA",
    )


@pytest.fixture
def qubit_qb():
    return TransmonQubit(
        EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
        T1=10000, T2=8000, flux=0.0, n_levels=3, name="QB",
    )


@pytest.fixture
def chip_2q(qubit_qa, qubit_qb) -> ChipTopology:
    """Two-qubit chip (no resonators)."""
    return ChipTopology(qubits=[qubit_qa, qubit_qb])


# ---------------------------------------------------------------------------
# hilbert_dim
# ---------------------------------------------------------------------------

class TestHilbertDim:
    """Tests for ChipTopology.hilbert_dim()."""

    def test_two_qubits_nlevels_3(self, chip_2q):
        assert chip_2q.hilbert_dim() == 9  # 3 * 3

    def test_single_qubit(self, qubit_qa):
        chip = ChipTopology(qubits=[qubit_qa])
        assert chip.hilbert_dim() == 3

    def test_three_qubits(self):
        qA = TransmonQubit(EC=0.4, EJ=10, T1=10000, T2=8000, n_levels=2, name="A")
        qB = TransmonQubit(EC=0.4, EJ=10, T1=10000, T2=8000, n_levels=2, name="B")
        qC = TransmonQubit(EC=0.4, EJ=10, T1=10000, T2=8000, n_levels=3, name="C")
        chip = ChipTopology(qubits=[qA, qB, qC])
        assert chip.hilbert_dim() == 2 * 2 * 3  # 12

    def test_with_resonator(self, qubit_qa):
        r = Resonator(frequencies=[2 * np.pi * 5.0], n_levels=[4])
        chip = ChipTopology(qubits=[qubit_qa], resonators=[r])
        assert chip.hilbert_dim() == 3 * 4  # 12

    def test_empty_chip(self):
        chip = ChipTopology()
        assert chip.hilbert_dim() == 1


# ---------------------------------------------------------------------------
# lift_qubit_op
# ---------------------------------------------------------------------------

class TestLiftQubitOp:
    """Tests for ChipTopology.lift_qubit_op()."""

    def test_lift_dims_correct(self, chip_2q):
        """Lifted op has correct dims: full chip space."""
        a_op = destroy(3)
        lifted = chip_2q.lift_qubit_op(a_op, "QA")
        assert lifted.dims == [[3, 3], [3, 3]]

    def test_lift_qa_gives_correct_operator(self, chip_2q):
        """Lifting destroy(3) onto QA yields a ⊗ I."""
        a = destroy(3)
        I = qeye(3)
        expected = tensor(a, I)
        lifted = chip_2q.lift_qubit_op(a, "QA")
        # Compare matrix representations
        np.testing.assert_array_almost_equal(lifted.full(), expected.full())

    def test_lift_qb_gives_I_tensor_a(self, chip_2q):
        """Lifting destroy(3) onto QB yields I ⊗ a."""
        a = destroy(3)
        I = qeye(3)
        expected = tensor(I, a)
        lifted = chip_2q.lift_qubit_op(a, "QB")
        np.testing.assert_array_almost_equal(lifted.full(), expected.full())

    def test_lift_unknown_qubit_raises(self, chip_2q):
        """Unknown qubit name raises ValueError."""
        with pytest.raises(ValueError, match="not found"):
            chip_2q.lift_qubit_op(destroy(3), "QC")

    def test_lift_preserves_hermiticity(self, chip_2q):
        """Hermitian op stays Hermitian after lifting."""
        n_op = destroy(3).dag() * destroy(3)  # number operator
        lifted = chip_2q.lift_qubit_op(n_op, "QA")
        np.testing.assert_array_almost_equal(
            lifted.full(), lifted.dag().full(),
        )

    def test_lift_three_qubits(self):
        """3-qubit chip: lift to correct position."""
        qA = TransmonQubit(EC=0.4, EJ=10, T1=10000, T2=8000, n_levels=2, name="A")
        qB = TransmonQubit(EC=0.4, EJ=10, T1=10000, T2=8000, n_levels=2, name="B")
        qC = TransmonQubit(EC=0.4, EJ=10, T1=10000, T2=8000, n_levels=2, name="C")
        chip = ChipTopology(qubits=[qA, qB, qC])
        a = destroy(2)
        I = qeye(2)
        # Lift onto B → I ⊗ a ⊗ I
        expected = tensor(I, a, I)
        lifted = chip.lift_qubit_op(a, "B")
        np.testing.assert_array_almost_equal(lifted.full(), expected.full())

    def test_lift_with_resonator(self, qubit_qa):
        """Lift onto qubit in presence of resonator."""
        r = Resonator(frequencies=[2 * np.pi * 5.0], n_levels=[2])
        chip = ChipTopology(qubits=[qubit_qa], resonators=[r])
        a = destroy(3)
        I_r = qeye(2)
        expected = tensor(a, I_r)
        lifted = chip.lift_qubit_op(a, "QA")
        np.testing.assert_array_almost_equal(lifted.full(), expected.full())


# ---------------------------------------------------------------------------
# hamiltonian_static
# ---------------------------------------------------------------------------

class TestHamiltonianStatic:
    """Tests for ChipTopology.hamiltonian_static()."""

    def test_is_qobj(self, chip_2q):
        H = chip_2q.hamiltonian_static()
        assert isinstance(H, Qobj)

    def test_correct_dims(self, chip_2q):
        H = chip_2q.hamiltonian_static()
        assert H.dims == [[3, 3], [3, 3]]

    def test_hermitian(self, chip_2q):
        H = chip_2q.hamiltonian_static()
        np.testing.assert_array_almost_equal(H.full(), H.dag().full())

    def test_equals_sum_of_individual(self, chip_2q, qubit_qa, qubit_qb):
        """H_chip = H_QA ⊗ I_QB + I_QA ⊗ H_QB."""
        H_chip = chip_2q.hamiltonian_static()
        H_QA = chip_2q.lift_qubit_op(
            qubit_qa.hamiltonian_static(), "QA",
        )
        H_QB = chip_2q.lift_qubit_op(
            qubit_qb.hamiltonian_static(), "QB",
        )
        expected = H_QA + H_QB
        np.testing.assert_array_almost_equal(H_chip.full(), expected.full())


# ---------------------------------------------------------------------------
# collapse_operators
# ---------------------------------------------------------------------------

class TestCollapseOperators:
    """Tests for ChipTopology.collapse_operators()."""

    def test_returns_list_of_qobj(self, chip_2q):
        c_ops = chip_2q.collapse_operators()
        assert isinstance(c_ops, list)
        # 2 qubits × 2 collapse ops each = 4
        assert len(c_ops) == 4
        for c in c_ops:
            assert isinstance(c, Qobj)

    def test_correct_dims(self, chip_2q):
        c_ops = chip_2q.collapse_operators()
        for c in c_ops:
            assert c.dims == [[3, 3], [3, 3]]


# ---------------------------------------------------------------------------
# from_legacy_coupled_system
# ---------------------------------------------------------------------------

class TestFromLegacyCoupledSystem:
    """Tests for ChipTopology.from_legacy_coupled_system()."""

    def test_extracts_qubits(self, qubit_qa, qubit_qb):
        """Build ChipTopology from CoupledSystem."""
        r = Resonator(frequencies=[2 * np.pi * 5.0], n_levels=[3])
        cs = CoupledSystem(qubit_qa, qubit_qb, r)
        chip = ChipTopology.from_legacy_coupled_system(cs)
        assert len(chip.qubits) == 2
        assert chip.qubits[0] is qubit_qa
        assert chip.qubits[1] is qubit_qb
        assert len(chip.resonators) == 1
        assert chip.resonators[0] is r


# ---------------------------------------------------------------------------
# get_qubit
# ---------------------------------------------------------------------------

class TestGetQubit:
    """Tests for ChipTopology.get_qubit()."""

    def test_get_existing(self, chip_2q):
        q = chip_2q.get_qubit("QA")
        assert q.name == "QA"

    def test_get_nonexistent_raises(self, chip_2q):
        with pytest.raises(ValueError, match="not found"):
            chip_2q.get_qubit("Nonexistent")
