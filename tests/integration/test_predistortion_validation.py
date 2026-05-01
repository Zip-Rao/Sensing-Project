"""Integration test for the P4 predistortion validation workflow."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.calibration.predistortion import PredistortionDesigner
from sqc.control.waveform import Waveform
from sqc.devices.transmon import TransmonQubit
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.workflows.predistortion_validation import PredistortionValidationWorkflow


pytestmark = pytest.mark.integration


def _make_qubit():
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        n_levels=2,
    )


def test_predistortion_validation_improves_rmse_by_more_than_10x():
    t = np.linspace(0.0, 120.0, 1024)
    target = Waveform(
        t_list=t,
        samples=0.5
        * (np.tanh((t - 20.0) / 2.0) - np.tanh((t - 80.0) / 2.0)),
    )
    workflow = PredistortionValidationWorkflow(
        qubit=_make_qubit(),
        target_waveform=target,
        true_distortion=SingleExponentialDistortion(amplitude=0.15, tau=12.0),
        designer=PredistortionDesigner(method="iir_inverse"),
        fit_type="iir",
    )

    result = workflow.run()

    assert result["metrics"]["improvement_factor"] > 10.0
    assert result["metrics"]["rmse_corrected"] < 1e-10
    assert result["on_chip_corrected"].samples.shape == target.samples.shape
