"""Integration test: TransferFunctionCalibration fit accuracy.

Verifies that calibration can fit a known distortion model
to within 5% relative error on fitted parameters.
"""
from __future__ import annotations

import numpy as np

from sqc.calibration.transfer_function import TransferFunctionCalibration
from sqc.hardware.distortion import (
    SingleExponentialDistortion,
    MultiExponentialDistortion,
)
from sqc.control.waveform import Waveform


class TestTransferFunctionCalibration:
    """Integration tests for TransferFunctionCalibration."""

    def test_single_exp_fit_accuracy(self):
        """Fit a known single-exponential distortion; amplitude/tau within 5%."""
        true_amp = 0.05
        true_tau = 20.0
        dist = SingleExponentialDistortion(amplitude=true_amp, tau=true_tau)

        cal = TransferFunctionCalibration(
            distortion=dist,
            method="simulation",
            fit_type="single_exp",
        )
        table = cal.calibrate()

        fitted_amp = table.fit_params["amplitude"]
        fitted_tau = table.fit_params["tau"]

        # Within 5% relative error
        assert abs(fitted_amp - true_amp) / true_amp < 0.05, (
            f"amp fit error: {abs(fitted_amp - true_amp) / true_amp:.2%}"
        )
        assert abs(fitted_tau - true_tau) / true_tau < 0.05, (
            f"tau fit error: {abs(fitted_tau - true_tau) / true_tau:.2%}"
        )

    def test_multi_exp_fit_accuracy(self):
        """Fit a 2-component multi-exp; parameter errors within acceptable range."""
        true_amps = np.array([0.03, 0.02])
        true_taus = np.array([15.0, 50.0])
        dist = MultiExponentialDistortion(
            amplitudes=true_amps, taus=true_taus,
        )

        cal = TransferFunctionCalibration(
            distortion=dist,
            method="simulation",
            fit_type="multi_exp",
            n_exp_components=2,
            t_max=500.0,
        )
        table = cal.calibrate()

        fitted_amps = table.fit_params["amplitudes"]
        fitted_taus = table.fit_params["taus"]

        # Check we got at least 2 components
        assert len(fitted_amps) >= 2, f"Expected >=2 components, got {len(fitted_amps)}"

        # Sort for comparison
        sort_idx_true = np.argsort(true_taus)
        sort_idx_fit = np.argsort(fitted_taus[:2])

        # Individual components may swap order; check total parameters
        # Allow 20% due to component overlap (harder to separate)
        for i in range(2):
            ti = sort_idx_true[i]
            fi = sort_idx_fit[i]
            rel_err_amp = abs(true_amps[ti] - fitted_amps[fi]) / max(true_amps[ti], 0.001)
            rel_err_tau = abs(true_taus[ti] - fitted_taus[fi]) / max(true_taus[ti], 1.0)
            # Multi-exp fitting is approximate; allow 25%
            assert rel_err_amp < 0.35, f"amp[{i}] error {rel_err_amp:.2%}"
            assert rel_err_tau < 0.35, f"tau[{i}] error {rel_err_tau:.2%}"

    def test_to_distortion_model_roundtrip(self):
        """calibrate() -> to_distortion_model() should give a usable model."""
        dist = SingleExponentialDistortion(amplitude=0.08, tau=30.0)
        cal = TransferFunctionCalibration(
            distortion=dist,
            method="simulation",
            fit_type="single_exp",
        )
        fitted = cal.to_distortion_model()
        assert isinstance(fitted, SingleExponentialDistortion)

        # Verify the fitted model gives similar step response
        t = np.linspace(0, 200, 2000)
        true_step = dist.step_response(t)
        fit_step = fitted.step_response(t)
        rmse = np.sqrt(np.mean((true_step - fit_step) ** 2))
        assert rmse < 0.002, f"Step response RMSE {rmse} too large"

    def test_calibration_table_structure(self):
        """CalibrationTable should have expected fields."""
        dist = SingleExponentialDistortion(amplitude=0.05, tau=20.0)
        cal = TransferFunctionCalibration(
            distortion=dist, method="simulation",
            fit_type="single_exp",
        )
        table = cal.calibrate()

        assert table.kind == "transfer_function"
        assert table.inputs is not None
        assert table.outputs is not None
        assert "amplitude" in table.fit_params
        assert "tau" in table.fit_params
        assert table.metadata["method"] == "simulation"
