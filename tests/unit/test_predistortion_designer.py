"""Unit tests for PredistortionDesigner.

Tests the inverse filter design and predistortion application.
Minimum 6 test cases covering auto-detection, IIR inverse,
predistort convenience, FIR fallback, and pole stability.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import Waveform
from sqc.hardware.distortion import (
    SingleExponentialDistortion,
    MultiExponentialDistortion,
    IIRDistortion,
    FIRDistortion,
)
from sqc.calibration.predistortion import PredistortionDesigner


@pytest.fixture
def t_default():
    return np.linspace(0, 200, 2000)


@pytest.fixture
def step_waveform(t_default):
    return Waveform(
        t_list=t_default,
        samples=np.where((t_default > 30) & (t_default < 100), 1.0, 0.0),
    )


@pytest.fixture
def dist():
    return SingleExponentialDistortion(amplitude=0.05, tau=20.0)


class TestPredistortionDesigner:
    """Tests for PredistortionDesigner."""

    # ---- auto-detection ----

    def test_auto_method_single_exp(self, dist):
        """auto should detect SingleExponential and use IIR inverse."""
        designer = PredistortionDesigner(method="auto")
        dt = 0.1
        inverse = designer.design(dist, dt=dt)
        assert isinstance(inverse, IIRDistortion)

    def test_auto_method_multi_exp(self):
        """auto should detect MultiExponential and use freq inverse."""
        mdist = MultiExponentialDistortion(
            amplitudes=np.array([0.03, 0.02]),
            taus=np.array([10.0, 50.0]),
        )
        designer = PredistortionDesigner(method="auto")
        dt = 0.1
        inverse = designer.design(mdist, dt=dt)
        # MultiExp uses CustomTransferDistortion (frequency inverse)
        cls_name = type(inverse).__name__
        assert "CustomTransfer" in cls_name or "IIR" in cls_name

    # ---- IIR inverse (analytical) ----

    def test_iir_inverse_perfect_cancellation(self, step_waveform, dist):
        """IIR inverse should give nearly perfect cascade."""
        designer = PredistortionDesigner(method="iir_inverse")
        dt = float(step_waveform.t_list[1] - step_waveform.t_list[0])
        inverse = designer.design(dist, dt=dt)

        # Cascade: inverse -> forward
        predistorted = inverse.apply_to_waveform(step_waveform)
        corrected = dist.apply_to_waveform(predistorted)

        rmse = np.sqrt(
            np.mean((corrected.samples - step_waveform.samples) ** 2)
        )
        assert rmse < 1e-12, f"RMSE {rmse} should be near zero"

    def test_iir_inverse_stable(self, dist):
        """IIR inverse poles should be inside unit circle."""
        designer = PredistortionDesigner(method="iir_inverse")
        dt = 0.1
        inverse = designer.design(dist, dt=dt)
        assert isinstance(inverse, IIRDistortion)
        assert PredistortionDesigner.check_pole_stability(
            inverse.b_coeffs, inverse.a_coeffs,
        )

    # ---- predistort() convenience ----

    def test_predistort_method(self, step_waveform, dist):
        """predistort() should return a waveform that improves cascade."""
        designer = PredistortionDesigner(method="auto")
        predistorted = designer.predistort(step_waveform, transfer=dist)
        corrected = dist.apply_to_waveform(predistorted)

        # Without predistortion
        uncorrected = dist.apply_to_waveform(step_waveform)
        rmse_uncorrected = np.sqrt(
            np.mean((uncorrected.samples - step_waveform.samples) ** 2)
        )
        rmse_corrected = np.sqrt(
            np.mean((corrected.samples - step_waveform.samples) ** 2)
        )
        improvement = rmse_uncorrected / max(rmse_corrected, 1e-30)
        assert improvement > 1e9, f"Improvement {improvement} should be huge"

    def test_predistort_with_precomputed_inverse(self, step_waveform, dist):
        """predistort() with pre-computed inverse_model should skip design."""
        designer = PredistortionDesigner(method="auto")
        dt = float(step_waveform.t_list[1] - step_waveform.t_list[0])
        inverse = designer.design(dist, dt=dt)

        predistorted = designer.predistort(
            step_waveform, inverse_model=inverse,
        )
        assert predistorted.samples.shape == step_waveform.samples.shape

    # ---- FIR inverse ----

    def test_fir_inverse_produces_fir(self, dist):
        """fir_inverse method should produce FIRDistortion."""
        designer = PredistortionDesigner(method="fir_inverse", n_taps=32)
        dt = 0.1
        inverse = designer.design(dist, dt=dt)
        assert isinstance(inverse, FIRDistortion)
        assert len(inverse.taps) == 32

    # ---- pole stability ----

    def test_check_pole_stability_stable(self):
        """Stable filter should return True."""
        b = np.array([1.0, 0.5])
        a = np.array([1.0, -0.7])  # pole at 0.7
        assert bool(PredistortionDesigner.check_pole_stability(b, a)) is True

    def test_check_pole_stability_unstable(self):
        """Unstable filter should return False."""
        b = np.array([1.0])
        a = np.array([1.0, -1.5])  # pole at 1.5
        assert bool(PredistortionDesigner.check_pole_stability(b, a)) is False

    # ---- error handling ----

    def test_no_transfer_or_inverse_raises(self, step_waveform):
        """predistort() with no transfer and no inverse should raise."""
        designer = PredistortionDesigner(method="auto")
        with pytest.raises(ValueError, match="Either transfer or inverse"):
            designer.predistort(step_waveform)

    def test_design_no_dt_raises(self, dist):
        """design() with no dt should raise TypeError (required argument)."""
        designer = PredistortionDesigner(method="fir_inverse")
        with pytest.raises(TypeError):
            designer.design(dist)  # type: ignore[call-arg]
