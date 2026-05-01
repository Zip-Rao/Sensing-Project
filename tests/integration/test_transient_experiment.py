"""Integration test: TransientSensingExperiment matches baseline."""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import assert_array_close, load_baseline


def _make_qubit():
    from src.qubit import TransmonQubit
    return TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, state=0, n_levels=2,
    )


@pytest.mark.slow
def test_transient_experiment_matches_baseline():
    """TransientSensingExperiment.run() matches transient_default baseline."""
    from sqc.experiments.transient import TransientSensingExperiment

    q = _make_qubit()
    exp = TransientSensingExperiment(qubit=q)
    result = exp.run()

    bl = load_baseline("transient_default")
    assert_array_close(result.data["kernel"], bl["kernel"], name="kernel")
    assert_array_close(result.data["delta_p"], bl["delta_p"], name="delta_p")
    assert_array_close(result.data["p_e"], bl["p_e"], name="p_e")


@pytest.mark.slow
def test_transient_experiment_matches_old_protocal():
    """TransientSensingExperiment matches src.Protocal(type=4).evolve."""
    from sqc.experiments.transient import TransientSensingExperiment
    from src.protocal import Protocal as OldProtocal

    q1 = _make_qubit()
    exp = TransientSensingExperiment(qubit=q1)
    result = exp.run()

    q2 = _make_qubit()
    old_p = OldProtocal(type=4)
    old_p.initialize(q2, state=0)
    (t_s, k_old, scan_old, dp_old,
     pe_old, Phi_old, cp_old) = old_p.evolve(q2)

    assert_array_close(result.data["kernel"], np.asarray(k_old), name="kernel")
    assert_array_close(result.data["delta_p"], np.asarray(dp_old), name="delta_p")
    assert_array_close(result.data["p_e"], np.asarray(pe_old), name="p_e")
