"""Integration test for the P5 Z-crosstalk workflow."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import Waveform
from sqc.devices.chip import ChipTopology
from sqc.devices.transmon import TransmonQubit
from sqc.hardware.transfer_matrix import TransferMatrix
from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow


pytestmark = pytest.mark.integration


def _make_qubit(name: str) -> TransmonQubit:
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        n_levels=2,
        name=name,
    )


def test_z_crosstalk_workflow_extracts_and_compensates_crosstalk():
    t = np.linspace(0.0, 120.0, 1024)
    pulse = Waveform(
        t_list=t,
        samples=0.5
        * (np.tanh((t - 20.0) / 2.0) - np.tanh((t - 80.0) / 2.0)),
    )
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

    assert result["fit_error_dB"] < -20.0
    assert result["compensation_factor"] > 5.0
    assert result["phi_B_after_compensation"].samples.shape == pulse.samples.shape
