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
def test_transient_experiment_self_consistent():
    """Two TransientSensingExperiment runs produce identical results."""
    from sqc.experiments.transient import TransientSensingExperiment

    q1 = _make_qubit()
    r1 = TransientSensingExperiment(qubit=q1).run()

    q2 = _make_qubit()
    r2 = TransientSensingExperiment(qubit=q2).run()

    assert_array_close(r1.data["kernel"], r2.data["kernel"], name="kernel")
    assert_array_close(r1.data["delta_p"], r2.data["delta_p"], name="delta_p")
    assert_array_close(r1.data["p_e"], r2.data["p_e"], name="p_e")
