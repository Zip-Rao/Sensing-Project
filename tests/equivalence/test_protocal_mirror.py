"""Equivalence tests: src_mirror.protocal.Protocal == src.protocal.Protocal.

For each protocol case (0, 1, 2, 4), verify that the old and new
Protocal implementations produce identical numeric results within
the physics regression tolerance (rtol=1e-6, atol=1e-9).
"""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import assert_array_close

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_qubit(n_levels=2):
    """Create a fresh TransmonQubit with default test parameters."""
    from src.qubit import TransmonQubit

    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=n_levels,
    )


# ---------------------------------------------------------------------------
# Case 0: Rabi
# ---------------------------------------------------------------------------

def test_protocal_mirror_case_0_rabi():
    """Case 0 (Rabi): old vs new expect values must match."""
    from src.protocal import Protocal as OldProtocal
    from src_mirror.protocal import Protocal as NewProtocal

    q_old = _make_qubit()
    old_p = OldProtocal(type=0)
    old_p.initialize(q_old, state=0)
    old_result = old_p.evolve(q_old)

    q_new = _make_qubit()
    new_p = NewProtocal(type=0)
    new_p.initialize(q_new, state=0)
    new_result = new_p.evolve(q_new)

    assert_array_close(
        np.array(old_result.expect),
        np.array(new_result.expect),
        name="case_0_expect",
    )


# ---------------------------------------------------------------------------
# Case 1: Ramsey
# ---------------------------------------------------------------------------

def test_protocal_mirror_case_1_ramsey():
    """Case 1 (Ramsey): old vs new p_e list must match."""
    from src.protocal import Protocal as OldProtocal
    from src_mirror.protocal import Protocal as NewProtocal

    q_old = _make_qubit()
    old_p = OldProtocal(type=1)
    old_p.initialize(q_old, state=0)
    Phi_old, tau_old, pe_old = old_p.evolve(q_old)

    q_new = _make_qubit()
    new_p = NewProtocal(type=1)
    new_p.initialize(q_new, state=0)
    Phi_new, tau_new, pe_new = new_p.evolve(q_new)

    assert_array_close(np.asarray(pe_old), np.asarray(pe_new), name="case_1_p_e")
    assert_array_close(np.asarray(tau_old), np.asarray(tau_new), name="case_1_tau")
    assert_array_close(
        np.asarray(Phi_old.signal),
        np.asarray(Phi_new.signal),
        name="case_1_Phi",
    )


# ---------------------------------------------------------------------------
# Case 2: DiffEcho
# ---------------------------------------------------------------------------

def test_protocal_mirror_case_2_diff_echo():
    """Case 2 (DiffEcho): old vs new p_e, k, t_int must match."""
    from src.protocal import Protocal as OldProtocal
    from src_mirror.protocal import Protocal as NewProtocal

    q_old = _make_qubit()
    old_p = OldProtocal(type=2)
    old_p.initialize(q_old, state=0)
    Phi_old, tau_old, pe_old, k_old, t_int_old = old_p.evolve(q_old)

    q_new = _make_qubit()
    new_p = NewProtocal(type=2)
    new_p.initialize(q_new, state=0)
    Phi_new, tau_new, pe_new, k_new, t_int_new = new_p.evolve(q_new)

    assert k_old == k_new
    assert t_int_old == pytest.approx(t_int_new)
    assert_array_close(np.asarray(pe_old), np.asarray(pe_new), name="case_2_p_e")
    assert_array_close(
        np.asarray(Phi_old.signal),
        np.asarray(Phi_new.signal),
        name="case_2_Phi",
    )


# ---------------------------------------------------------------------------
# Case 4: Transient sensing
# ---------------------------------------------------------------------------

def test_protocal_mirror_case_4_transient():
    """Case 4 (Transient): old vs new kernel, delta_p, p_e must match."""
    from src.protocal import Protocal as OldProtocal
    from src_mirror.protocal import Protocal as NewProtocal

    q_old = _make_qubit()
    old_p = OldProtocal(type=4)
    old_p.initialize(q_old, state=0)
    (t_s_old, k_old, scan_old, dp_old,
     pe_old, Phi_old, cp_old) = old_p.evolve(q_old)

    q_new = _make_qubit()
    new_p = NewProtocal(type=4)
    new_p.initialize(q_new, state=0)
    (t_s_new, k_new, scan_new, dp_new,
     pe_new, Phi_new, cp_new) = new_p.evolve(q_new)

    assert_array_close(np.asarray(k_old), np.asarray(k_new), name="case_4_kernel")
    assert_array_close(np.asarray(dp_old), np.asarray(dp_new), name="case_4_delta_p")
    assert_array_close(np.asarray(pe_old), np.asarray(pe_new), name="case_4_p_e")
    assert_array_close(np.asarray(scan_old), np.asarray(scan_new), name="case_4_scan")
