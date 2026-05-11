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


def test_diff_echo_experiment_self_consistent():
    """Two DiffEchoExperiment runs produce identical results."""
    from sqc.experiments.echo import DiffEchoExperiment

    q1 = _make_qubit()
    r1 = DiffEchoExperiment(qubit=q1).run()

    q2 = _make_qubit()
    r2 = DiffEchoExperiment(qubit=q2).run()

    assert_array_close(r1.data["p_e"], r2.data["p_e"], name="p_e_self")
