"""Integration tests for TransientSensingExperiment."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.devices.transmon import TransmonQubit
from sqc.experiments.transient import TransientSensingExperiment
from tests.conftest import assert_array_close, load_baseline


pytestmark = [pytest.mark.integration, pytest.mark.slow]


def _make_qubit():
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000,
        T2=8000,
        n_levels=2,
    )


def test_transient_experiment_matches_baseline():
    result = TransientSensingExperiment(qubit=_make_qubit()).run()
    baseline = load_baseline("transient_default")

    assert_array_close(result.axes["t_samples"], baseline["t_samples"], name="t_samples")
    assert_array_close(result.data["kernel"], baseline["kernel"], name="kernel")
    assert_array_close(result.axes["scan"], baseline["scan_list"], name="scan")
    assert_array_close(result.data["delta_p"], baseline["delta_p"], name="delta_p")
    assert_array_close(result.data["p_e"], baseline["p_e"], name="p_e")
    assert_array_close(result.data["flux_samples"], baseline["Phi_signal"], name="Phi")
