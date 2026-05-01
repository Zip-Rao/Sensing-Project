"""Integration test: DiffEchoExperiment matches baseline."""
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


def test_diff_echo_experiment_matches_baseline():
    """DiffEchoExperiment.run() matches diff_echo_default baseline."""
    from sqc.experiments.echo import DiffEchoExperiment

    q = _make_qubit()
    exp = DiffEchoExperiment(qubit=q)
    result = exp.run()

    bl = load_baseline("diff_echo_default")
    assert exp.k == bl["k"]
    assert exp.t_int == pytest.approx(bl["t_int"])
    assert_array_close(result.data["p_e"], bl["p_e_list"], name="p_e")


def test_diff_echo_experiment_matches_old_protocal():
    """DiffEchoExperiment p_e matches src.Protocal(type=2).evolve."""
    from sqc.experiments.echo import DiffEchoExperiment
    from src.protocal import Protocal as OldProtocal

    q1 = _make_qubit()
    exp = DiffEchoExperiment(qubit=q1)
    result = exp.run()

    q2 = _make_qubit()
    old_p = OldProtocal(type=2)
    old_p.initialize(q2, state=0)
    Phi_old, tau_old, pe_old, k_old, t_int_old = old_p.evolve(q2)

    assert_array_close(result.data["p_e"], np.asarray(pe_old), name="p_e_vs_old")
