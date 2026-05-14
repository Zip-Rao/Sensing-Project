"""Equivalence test: src_mirror.protocal.Protocal case 5 vs src.protocal.Protocal.

Cryoscope case 5. Uses a reduced truncation set (5 points) to keep
test runtime manageable. Verifies that old and new produce identical
varphi values within physics regression tolerance.
"""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import assert_array_close


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


@pytest.mark.slow
@pytest.mark.xfail(
    reason="sqc IQ formula fixed (atan2(0.5-pI, pQ-0.5)); "
    "src/ has legacy buggy formula (atan2(pQ-0.5, pI-0.5)). "
    "src/ frozen by R1 — divergence is intentional."
)
def test_protocal_mirror_case_5_cryoscope_equivalence():
    """Case 5 (Cryoscope): old vs new varphi must match.

    Uses a reduced truncation set (5 points) for speed.
    Full equivalence with 120 truncation points would take ~10 min.

    Currently xfail: sqc uses corrected IQ phase formula, src/ uses
    legacy buggy version.  See 2026-05-13 IQ formula fix.
    """
    from src.protocal import Protocal as OldProtocal
    from src_mirror.protocal import Protocal as NewProtocal
    from src.signal import Signal

    # Create a short flux signal for faster testing
    # Match the same Signal type and amplitude as the legacy case 5 default
    Phi_old = Signal(
        type=2, t_list=np.linspace(0, 40, 80), amplitude=0.01,
    )

    # Use a small subset of truncation times (every ~15th point, reverse order)
    full_trunc = Phi_old.t_list[60:10:-1]
    step = max(1, len(full_trunc) // 5)
    trunc_subset = full_trunc[::step][:5]

    # --- Old Protocal path ---
    # We need to monkey-patch the trunc_list and IQ_readout parameters
    import src.protocal as old_mod

    q_old = _make_qubit()
    old_p = OldProtocal(type=5)
    old_p.initialize(q_old, state=0)

    # Run old protocol by directly simulating the case 5 logic with subset
    # The legacy case 5 hardcodes the trunc_list, so we run it manually
    psi_e_old = __import__("qutip", fromlist=["basis"]).basis(q_old.n_levels, 1)
    tau_old = 100.0
    t_rabi_old = np.linspace(0, 10, 20)
    p_e_list_old = [[], []]

    for trunc in trunc_subset:
        # Create a fresh Phi copy for each truncation (matching legacy in-place
        # semantics by copying then truncating)
        Phi_copy = Signal(
            type=2, t_list=np.linspace(0, 40, 80), amplitude=0.01,
        )
        Phi_copy.truncate(0, trunc)
        q_old.qubit_in_mag(Phi_copy, frame=1, omega_d=q_old.frequency)
        p_e_I, p_e_Q = old_mod.IQ_readout(q_old, type=3, tau=tau_old)
        p_e_list_old[0].append(p_e_I)
        p_e_list_old[1].append(p_e_Q)

    varphi_old = np.arctan2(
        np.array(p_e_list_old[1]) - 0.5,
        np.array(p_e_list_old[0]) - 0.5,
    )[::-1]

    # --- New Protocal path (via CryoscopeExperiment) ---
    q_new = _make_qubit()
    new_p = NewProtocal(type=5)
    new_p.initialize(q_new, state=0)

    # Create matching flux signal and trunc_list for new path
    from sqc.control.flux_signal import FluxSignal

    phi_new = FluxSignal(
        type=2, t_list=np.linspace(0, 40, 80), amplitude=0.01,
    )
    from sqc.experiments.cryoscope import CryoscopeExperiment

    exp = CryoscopeExperiment(
        qubit=q_new,
        flux_signal=phi_new,
        trunc_list=trunc_subset,
        t_rabi=t_rabi_old,
        tau=tau_old,
    )
    result = exp.run()
    varphi_new = result.data["varphi"]

    # Compare varphi values
    assert_array_close(
        varphi_old, varphi_new,
        name="case_5_varphi",
    )
