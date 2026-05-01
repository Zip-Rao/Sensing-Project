"""Integration test: RamseyExperiment matches baseline."""
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


def test_ramsey_experiment_matches_baseline():
    """RamseyExperiment.run() output matches ramsey_default baseline."""
    from sqc.experiments.ramsey import RamseyExperiment

    q = _make_qubit()
    exp = RamseyExperiment(qubit=q)
    result = exp.run()

    bl = load_baseline("ramsey_default")
    assert_array_close(result.data["p_e"], bl["p_e_list"], name="p_e")
    assert_array_close(result.axes["tau"], bl["tau_list"], name="tau")


def test_ramsey_experiment_matches_old_protocal():
    """RamseyExperiment p_e matches src.Protocal(type=1).evolve."""
    from sqc.experiments.ramsey import RamseyExperiment
    from src.protocal import Protocal as OldProtocal

    q1 = _make_qubit()
    exp = RamseyExperiment(qubit=q1)
    result = exp.run()

    q2 = _make_qubit()
    old_p = OldProtocal(type=1)
    old_p.initialize(q2, state=0)
    Phi_old, tau_old, pe_old = old_p.evolve(q2)

    assert_array_close(result.data["p_e"], np.asarray(pe_old), name="p_e_vs_old")
