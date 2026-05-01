"""Regression baseline for P5 Z-crosstalk behavior."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import Waveform
from sqc.devices.chip import ChipTopology
from sqc.devices.transmon import TransmonQubit
from sqc.hardware.transfer_matrix import TransferMatrix
from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow
from tests.conftest import assert_array_close, load_baseline


pytestmark = pytest.mark.regression


def _make_qubit(name: str) -> TransmonQubit:
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        n_levels=2,
        name=name,
    )


def test_z_crosstalk_default_baseline():
    bl = load_baseline("z_crosstalk_default")
    pulse = Waveform(t_list=bl["t"], samples=bl["pulse_A"])
    transfer = TransferMatrix.from_dc_matrix(
        dc_matrix=np.array([[1.0, 0.0], [0.04, 1.0]]),
        source_names=["QA", "QB"],
        target_names=["QA", "QB"],
    )
    chip = ChipTopology(
        qubits=[_make_qubit("QA"), _make_qubit("QB")],
        transfer_matrix=transfer,
    )

    result = ZCrosstalkWorkflow(
        chip=chip,
        flux_pulse_on_A=pulse,
        true_transfer_matrix=transfer,
    ).run()

    assert_array_close(result["phi_B_true"].samples, bl["phi_B_true"], name="phi_B")
    assert_array_close(
        result["compensation_pulse"].samples,
        bl["compensation_pulse"],
        name="compensation_pulse",
    )
    assert_array_close(
        result["phi_B_after_compensation"].samples,
        bl["phi_B_after_compensation"],
        atol=1e-12,
        name="phi_B_after",
    )
    assert_array_close(
        result["H_BA_estimated"].real,
        bl["H_BA_estimated_real"],
        name="H_BA_estimated_real",
    )
    assert result["metrics"]["compensation_factor"] == pytest.approx(
        bl["metrics"]["compensation_factor"],
        rel=1e-9,
    )
    assert result["metrics"]["compensation_factor"] > 5.0
