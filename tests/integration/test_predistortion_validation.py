"""Integration test: PredistortionValidationWorkflow end-to-end.

Verifies that the full workflow achieves > 10x improvement
for a single-exponential distortion.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import Waveform
from sqc.hardware.distortion import (
    SingleExponentialDistortion,
    MultiExponentialDistortion,
)
from sqc.calibration.waveform import PredistortionDesigner
from sqc.workflows.predistortion_validation import PredistortionValidationWorkflow


class TestPredistortionValidationWorkflow:
    """Integration tests for the full predistortion validation workflow."""

    @pytest.fixture
    def target(self):
        t = np.linspace(0, 200, 2000)
        return Waveform(
            t_list=t,
            samples=np.where((t > 30) & (t < 100), 1.0, 0.0),
        )

    def test_improvement_greater_than_10x_single_exp(self, target):
        """Single exponential distortion: improvement > 10x."""
        dist = SingleExponentialDistortion(amplitude=0.05, tau=20.0)
        designer = PredistortionDesigner(method="auto")

        wf = PredistortionValidationWorkflow(
            target_waveform=target,
            true_distortion=dist,
            designer=designer,
        )
        result = wf.run()
        metrics = result["metrics"]

        assert metrics["improvement_factor"] > 10, (
            f"Improvement {metrics['improvement_factor']:.1f} <= 10"
        )
        assert metrics["rmse_corrected"] < metrics["rmse_uncorrected"]

    def test_improvement_greater_than_100x_strong_distortion(self, target):
        """Stronger distortion (amp=0.2): improvement > 100x."""
        dist = SingleExponentialDistortion(amplitude=0.2, tau=30.0)
        designer = PredistortionDesigner(method="auto")

        wf = PredistortionValidationWorkflow(
            target_waveform=target,
            true_distortion=dist,
            designer=designer,
        )
        result = wf.run()
        metrics = result["metrics"]

        assert metrics["improvement_factor"] > 100, (
            f"Improvement {metrics['improvement_factor']:.1f} <= 100"
        )

    def test_workflow_returns_expected_keys(self, target):
        """Workflow result dict should contain all expected keys."""
        dist = SingleExponentialDistortion(amplitude=0.05, tau=20.0)

        wf = PredistortionValidationWorkflow(
            target_waveform=target,
            true_distortion=dist,
        )
        result = wf.run()

        expected_keys = {
            "target", "on_chip_uncorrected", "on_chip_corrected",
            "awg_predistorted", "inverse_model", "measured_model", "metrics",
        }
        assert set(result.keys()) == expected_keys

        # Metrics sub-keys
        expected_metrics = {
            "rmse_uncorrected", "rmse_corrected", "improvement_factor",
            "settling_uncorrected_ns", "settling_corrected_ns",
            "calibration_fit_type", "inverse_model_type",
        }
        assert set(result["metrics"].keys()) >= expected_metrics

    def test_settling_time_reduces(self, target):
        """Corrected settling time should be less than uncorrected."""
        dist = SingleExponentialDistortion(amplitude=0.1, tau=25.0)
        designer = PredistortionDesigner(method="auto")

        wf = PredistortionValidationWorkflow(
            target_waveform=target,
            true_distortion=dist,
            designer=designer,
        )
        result = wf.run()
        metrics = result["metrics"]

        assert metrics["settling_corrected_ns"] <= metrics["settling_uncorrected_ns"]

    def test_multi_exp_workflow_runs(self, target):
        """MultiExp distortion: workflow should run without error."""
        dist = MultiExponentialDistortion(
            amplitudes=np.array([0.03, 0.02]),
            taus=np.array([15.0, 50.0]),
        )
        designer = PredistortionDesigner(method="auto")

        wf = PredistortionValidationWorkflow(
            target_waveform=target,
            true_distortion=dist,
            designer=designer,
        )
        result = wf.run()
        # Just verify it runs and returns valid metrics
        assert result["metrics"]["improvement_factor"] > 0
        assert result["metrics"]["rmse_uncorrected"] > 0
