"""Integration tests for CryoscopeExperiment."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.devices.transmon import TransmonQubit
from sqc.experiments.cryoscope import CryoscopeExperiment
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


def test_cryoscope_experiment_matches_baseline():
    result = CryoscopeExperiment(qubit=_make_qubit()).run()
    baseline = load_baseline("cryoscope_default")

    assert_array_close(result.axes["trunc"], baseline["trunc_list"], name="trunc")
    assert_array_close(result.data["varphi"], baseline["varphi"], name="varphi")
    assert_array_close(result.data["p_e_I"], baseline["p_e_I"], name="p_e_I")
    assert_array_close(result.data["p_e_Q"], baseline["p_e_Q"], name="p_e_Q")
