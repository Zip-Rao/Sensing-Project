"""Regression baseline for P4 predistortion behavior."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.calibration.predistortion import PredistortionDesigner
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.workflows.predistortion_validation import PredistortionValidationWorkflow
from tests.conftest import assert_array_close, load_baseline


pytestmark = pytest.mark.regression


def test_predistortion_default_baseline(qubit_default):
    bl = load_baseline("predistortion_default")
    target = Waveform(t_list=bl["t"], samples=bl["target"])
    result = PredistortionValidationWorkflow(
        qubit=qubit_default,
        target_waveform=target,
        true_distortion=SingleExponentialDistortion(amplitude=0.15, tau=12.0),
        designer=PredistortionDesigner(method="iir_inverse"),
        fit_type="iir",
    ).run()

    assert_array_close(
        result["on_chip_uncorrected"].samples,
        bl["on_chip_uncorrected"],
        name="on_chip_uncorrected",
    )
    assert_array_close(
        result["awg_predistorted"].samples,
        bl["awg_predistorted"],
        name="awg_predistorted",
    )
    assert_array_close(
        result["on_chip_corrected"].samples,
        bl["on_chip_corrected"],
        atol=1e-12,
        name="on_chip_corrected",
    )
    assert result["metrics"]["improvement_factor"] == pytest.approx(
        bl["metrics"]["improvement_factor"], rel=1e-9
    )
    assert result["metrics"]["improvement_factor"] > 10.0
