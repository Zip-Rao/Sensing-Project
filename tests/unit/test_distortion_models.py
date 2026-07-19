"""Unit tests for DistortionModel subclasses.

Covers all 5 subclasses: SingleExponential, MultiExponential,
FIR, IIR, CustomTransfer. Minimum 3 test cases per subclass.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.signal import lfilter

from sqc.control.waveform import Waveform
from sqc.hardware.distortion import (
    DistortionModel,
    SingleExponentialDistortion,
    MultiExponentialDistortion,
    FIRDistortion,
    IIRDistortion,
    CustomTransferDistortion,
)


# ===========================================================================
# Helpers
# ===========================================================================

@pytest.fixture
def t_default():
    return np.linspace(0, 200, 2000)


@pytest.fixture
def step_waveform(t_default):
    """A step from 0 to 1 at t=10 ns."""
    return Waveform(
        t_list=t_default,
        samples=np.where(t_default > 10, 1.0, 0.0),
    )


# ===========================================================================
# SingleExponentialDistortion (4 cases)
# ===========================================================================

class TestSingleExponentialDistortion:
    """Tests for SingleExponentialDistortion."""

    def test_step_response_settles_to_one(self):
        """Step response must asymptote to exactly 1.0."""
        dist = SingleExponentialDistortion(amplitude=0.05, tau=20.0)
        t = np.linspace(0, 200, 1000)
        step = dist.step_response(t)
        assert step[0] == pytest.approx(1.0 - 0.05, rel=1e-6)
        assert step[-1] == pytest.approx(1.0, rel=1e-4)

    def test_impulse_response_dc_gain_one(self):
        """Impulse response integrated over time should give DC gain = 1."""
        dist = SingleExponentialDistortion(amplitude=0.1, tau=30.0)
        t = np.linspace(0, 500, 5000)
        h = dist.impulse_response(t)
        dt = t[1] - t[0]
        integ = np.sum(h) * dt
        assert integ == pytest.approx(1.0, rel=1e-3)

    def test_frequency_response_dc(self):
        """H(omega=0) must equal 1.0 (DC gain)."""
        dist = SingleExponentialDistortion(amplitude=0.2, tau=15.0)
        H_dc = dist.frequency_response(np.array([0.0]))
        assert abs(H_dc[0] - 1.0) < 1e-12
        assert abs(H_dc[0].imag) < 1e-12

    def test_apply_to_waveform_step(self, step_waveform):
        """Applying to a step waveform should produce expected distortion."""
        dist = SingleExponentialDistortion(amplitude=0.05, tau=20.0)
        out = dist.apply_to_waveform(step_waveform)
        # The step response should settle to 1.0
        assert out.samples[-1] == pytest.approx(1.0, rel=1e-3)
        # After the step edge (t > 10), the value should approach 1
        # At the step edge, value = 1 - amp = 0.95
        idx_step = np.searchsorted(step_waveform.t_list, 10)
        assert out.samples[idx_step + 1] == pytest.approx(0.95, abs=0.01)

    def test_iir_coeffs_positive_dt(self):
        """_iir_coeffs should produce valid filter coefficients."""
        dist = SingleExponentialDistortion(amplitude=0.05, tau=20.0)
        b, a = dist._iir_coeffs(dt=0.1)
        assert len(b) == 2
        assert len(a) == 2
        assert a[0] == pytest.approx(1.0)
        # Filter should be stable: pole inside unit circle
        roots = np.roots(a)
        assert np.all(np.abs(roots) < 1.0)


# ===========================================================================
# MultiExponentialDistortion (4 cases)
# ===========================================================================

class TestMultiExponentialDistortion:
    """Tests for MultiExponentialDistortion."""

    def test_step_response_settles_to_one(self):
        """Step response must asymptote to 1.0."""
        dist = MultiExponentialDistortion(
            amplitudes=np.array([0.03, 0.02]),
            taus=np.array([10.0, 50.0]),
        )
        t = np.linspace(0, 500, 2000)
        step = dist.step_response(t)
        assert step[0] == pytest.approx(1.0 - 0.03 - 0.02, rel=1e-6)
        assert step[-1] == pytest.approx(1.0, rel=1e-4)

    def test_frequency_response_dc(self):
        """H(omega=0) must equal 1.0."""
        dist = MultiExponentialDistortion(
            amplitudes=np.array([0.01, 0.02, 0.03]),
            taus=np.array([5.0, 20.0, 100.0]),
        )
        H_dc = dist.frequency_response(np.array([0.0]))
        assert abs(H_dc[0] - 1.0) < 1e-12

    def test_n_components(self):
        """n_components should match input array length."""
        dist = MultiExponentialDistortion(
            amplitudes=np.array([0.01, 0.02]),
            taus=np.array([10.0, 50.0]),
        )
        assert dist.n_components == 2

    def test_apply_to_waveform(self, step_waveform):
        """apply_to_waveform should produce valid output."""
        dist = MultiExponentialDistortion(
            amplitudes=np.array([0.03, 0.02]),
            taus=np.array([15.0, 100.0]),
        )
        out = dist.apply_to_waveform(step_waveform)
        assert out.samples.shape == step_waveform.samples.shape
        # Long tau=100ns tail: exp(-200/100) ≈ 0.135, residual ≈ 0.02*0.135 ≈ 0.003
        assert out.samples[-1] == pytest.approx(1.0, abs=0.005)

    def test_empty_components(self):
        """Zero components should pass through unchanged."""
        dist = MultiExponentialDistortion(
            amplitudes=np.array([]),
            taus=np.array([]),
        )
        wf = np.random.randn(100)
        out = dist.apply(wf, dt=0.1)
        np.testing.assert_array_almost_equal(out, wf)


# ===========================================================================
# FIRDistortion (3 cases)
# ===========================================================================

class TestFIRDistortion:
    """Tests for FIRDistortion."""

    def test_identity_taps(self, step_waveform):
        """Taps = [1.0] should pass through unchanged."""
        dist = FIRDistortion(taps=np.array([1.0]))
        out = dist.apply_to_waveform(step_waveform)
        np.testing.assert_array_almost_equal(out.samples, step_waveform.samples)

    def test_step_response_from_unit_taps(self):
        """Single unit tap: step response should be 1 everywhere."""
        dist = FIRDistortion(taps=np.array([1.0]))
        t = np.linspace(0, 100, 1000)
        step = dist.step_response(t)
        np.testing.assert_array_almost_equal(step, np.ones_like(t))

    def test_frequency_response_dc(self):
        """DC gain = sum of taps."""
        taps = np.array([0.2, 0.3, 0.5])
        dist = FIRDistortion(taps=taps)
        H_dc = dist.frequency_response(np.array([0.0]))
        assert abs(H_dc[0]) == pytest.approx(1.0, rel=1e-6)

    def test_two_tap_average(self):
        """Two-tap averaging filter [0.5, 0.5] should produce correct output."""
        dist = FIRDistortion(taps=np.array([0.5, 0.5]))
        wf = np.array([1.0, 0.0, 0.0, 0.0])
        out = dist.apply(wf, dt=1.0)
        expected = np.array([0.5, 0.5, 0.0, 0.0])
        np.testing.assert_array_almost_equal(out, expected)


# ===========================================================================
# IIRDistortion (3 cases)
# ===========================================================================

class TestIIRDistortion:
    """Tests for IIRDistortion."""

    def test_identity_coeffs(self, step_waveform):
        """b=[1], a=[1] should pass through unchanged."""
        dist = IIRDistortion(b_coeffs=np.array([1.0]), a_coeffs=np.array([1.0]))
        out = dist.apply_to_waveform(step_waveform)
        np.testing.assert_array_almost_equal(out.samples, step_waveform.samples)

    def test_step_response(self):
        """Step response of a leaky integrator."""
        # y[n] = x[n] + 0.5*y[n-1] (leaky integrator, DC gain = 2)
        # Normalize: b=[0.5], a=[1, -0.5], DC gain = 0.5/(1-0.5) = 1
        dist = IIRDistortion(
            b_coeffs=np.array([0.5]),
            a_coeffs=np.array([1.0, -0.5]),
        )
        t = np.linspace(0, 50, 500)
        step = dist.step_response(t)
        assert step[0] == pytest.approx(0.5)
        assert step[-1] == pytest.approx(1.0, rel=1e-4)

    def test_frequency_response_dc(self):
        """DC gain = sum(b) / sum(a) for normalized a[0]=1."""
        dist = IIRDistortion(
            b_coeffs=np.array([0.3, 0.2]),
            a_coeffs=np.array([1.0, -0.5]),
        )
        H_dc = dist.frequency_response(np.array([0.0]))
        expected_dc = np.sum([0.3, 0.2]) / np.sum([1.0, -0.5])
        assert abs(H_dc[0]) == pytest.approx(expected_dc, rel=1e-6)

    def test_normalize_a0(self):
        """If a[0] != 1, coefficients should be normalized."""
        dist = IIRDistortion(
            b_coeffs=np.array([2.0, 1.0]),
            a_coeffs=np.array([2.0, -1.0]),
        )
        assert dist.a_coeffs[0] == pytest.approx(1.0)
        assert dist.b_coeffs[0] == pytest.approx(1.0)
        assert dist.b_coeffs[1] == pytest.approx(0.5)
        assert dist.a_coeffs[1] == pytest.approx(-0.5)


# ===========================================================================
# CustomTransferDistortion (3 cases)
# ===========================================================================

class TestCustomTransferDistortion:
    """Tests for CustomTransferDistortion."""

    def test_identity_passthrough(self, step_waveform):
        """Identity H(omega)=1 should pass through unchanged."""
        omega_grid = np.linspace(-np.pi, np.pi, 1024)
        H_grid = np.ones(1024, dtype=complex)
        dist = CustomTransferDistortion(
            omega_grid=omega_grid, H_grid=H_grid,
        )
        out = dist.apply_to_waveform(step_waveform)
        np.testing.assert_array_almost_equal(
            out.samples, step_waveform.samples, decimal=10,
        )

    def test_step_response(self):
        """Step response of identity transfer function."""
        omega_grid = np.linspace(-np.pi, np.pi, 1024)
        H_grid = np.ones(1024, dtype=complex)
        dist = CustomTransferDistortion(
            omega_grid=omega_grid, H_grid=H_grid,
        )
        t = np.linspace(0, 100, 1000)
        step = dist.step_response(t)
        # Should be close to 1 for t >= 0
        mask = t >= 0
        assert np.all(np.abs(step[mask] - 1.0) < 0.01)

    def test_frequency_response_interpolation(self):
        """frequency_response should interpolate correctly."""
        omega_grid = np.linspace(-np.pi, np.pi, 1024)
        H_grid = omega_grid * 0.0 + (1.0 + 0j)  # identity
        dist = CustomTransferDistortion(
            omega_grid=omega_grid, H_grid=H_grid,
        )
        # Query at specific frequencies
        omega_query = np.array([-1.0, 0.0, 1.0])
        H = dist.frequency_response(omega_query)
        assert np.all(np.abs(H - 1.0) < 1e-12)

    def test_frequency_response_conjugate_symmetry(self):
        """For real-valued H, H(-omega) should equal conj(H(omega))."""
        n = 512
        omega_grid = np.linspace(-np.pi, np.pi, n)
        # Create a H(omega) with some structure
        H_grid = 1.0 + 0.5 * np.exp(-omega_grid**2 / 0.1) + 0j
        H_grid = H_grid.astype(complex)
        dist = CustomTransferDistortion(
            omega_grid=omega_grid, H_grid=H_grid,
        )
        # Check at omega and -omega
        omega_pos = np.array([0.5])
        omega_neg = np.array([-0.5])
        H_pos = dist.frequency_response(omega_pos)
        H_neg = dist.frequency_response(omega_neg)
        # For real H: both should be real, and equal
        assert abs(H_pos[0].imag) < 1e-12
        assert abs(H_neg[0].imag) < 1e-12
        assert H_pos[0].real == pytest.approx(H_neg[0].real, rel=1e-6)


# ===========================================================================
# DistortionModel ABC (2 cases)
# ===========================================================================

class TestDistortionModelABC:
    """Verify the ABC cannot be instantiated directly."""

    def test_cannot_instantiate_abc(self):
        """ABC should raise TypeError on instantiation."""
        with pytest.raises(TypeError):
            DistortionModel()  # type: ignore[abstract]

    def test_apply_to_waveform_convenience(self, step_waveform):
        """apply_to_waveform should work on any concrete subclass."""
        dist = SingleExponentialDistortion(amplitude=0.01, tau=50.0)
        out = dist.apply_to_waveform(step_waveform)
        assert out.samples.shape == step_waveform.samples.shape
        assert "distortion" in out.metadata
        assert out.metadata["distortion"] == "SingleExponentialDistortion"
