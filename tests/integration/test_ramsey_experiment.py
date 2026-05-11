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


def test_ramsey_experiment_self_consistent():
    """Two RamseyExperiment runs produce identical results."""
    from sqc.experiments.ramsey import RamseyExperiment

    q1 = _make_qubit()
    r1 = RamseyExperiment(qubit=q1).run()

    q2 = _make_qubit()
    r2 = RamseyExperiment(qubit=q2).run()

    assert_array_close(r1.data["p_e"], r2.data["p_e"], name="p_e_self")
    assert_array_close(r1.axes["tau"], r2.axes["tau"], name="tau_self")
