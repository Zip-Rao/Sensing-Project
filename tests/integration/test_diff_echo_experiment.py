"""Integration tests for DiffEchoExperiment."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.devices.transmon import TransmonQubit
from sqc.experiments.echo import DiffEchoExperiment
from tests.conftest import assert_array_close, load_baseline


pytestmark = pytest.mark.integration


def _make_qubit():
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000,
        T2=8000,
        n_levels=2,
    )


def test_diff_echo_experiment_matches_baseline():
    exp = DiffEchoExperiment(qubit=_make_qubit())
    result = exp.run()
    baseline = load_baseline("diff_echo_default")

    assert result.config["k"] == baseline["k"]
    assert result.config["t_int"] == pytest.approx(baseline["t_int"], rel=1e-12)
    assert_array_close(result.data["flux_samples"], baseline["Phi_signal"], name="Phi")
    assert_array_close(result.axes["tau"], baseline["tau_list"], name="tau")
    assert_array_close(result.data["p_e"], baseline["p_e_list"], name="p_e")
