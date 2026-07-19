"""Unit tests for CascadeDistortion."""

import numpy as np
import pytest

from sqc.hardware.distortion import (
    CascadeDistortion,
    FIRDistortion,
    IIRDistortion,
    SingleExponentialDistortion,
)


class TestCascadeDistortion:
    """P9.A.1 — CascadeDistortion model correctness."""

    @staticmethod
    def _step(t):
        return np.where(t >= 0, 1.0, 0.0)

    def test_empty_stages_is_identity(self):
        """Cascade with empty stages list is identity."""
        cascade = CascadeDistortion(stages=[])
        t = np.arange(0, 100, 0.5)
        x = np.sin(2 * np.pi * 0.01 * t)
        out = cascade.apply(x.copy(), dt=0.5)
        np.testing.assert_allclose(out, x, atol=1e-12)

        # frequency response is unity
        omega = 2 * np.pi * np.fft.fftfreq(256, d=0.5)
        H = cascade.frequency_response(omega)
        np.testing.assert_allclose(H, 1.0 + 0j, atol=1e-12)

    def test_single_exp_cascade_vs_manual(self):
        """Two SingleExponential in cascade matches manual sequential apply."""
        s1 = SingleExponentialDistortion(amplitude=0.05, tau=50.0)
        s2 = SingleExponentialDistortion(amplitude=0.03, tau=200.0)
        cascade = CascadeDistortion(stages=[s1, s2])

        dt = 0.5
        t = np.arange(0, 200, dt)
        x = self._step(t)

        out_cascade = cascade.apply(x.copy(), dt)
        out_manual = s2.apply(s1.apply(x.copy(), dt), dt)
        np.testing.assert_allclose(out_cascade, out_manual, atol=1e-12)

    def test_dc_gain_product(self):
        """DC gain of cascade equals product of individual DC gains."""
        s1 = SingleExponentialDistortion(amplitude=0.05, tau=50.0)
        s2 = SingleExponentialDistortion(amplitude=0.03, tau=200.0)
        cascade = CascadeDistortion(stages=[s1, s2])

        omega_dc = np.array([0.0])
        H_cascade = cascade.frequency_response(omega_dc)
        H1 = s1.frequency_response(omega_dc)
        H2 = s2.frequency_response(omega_dc)
        np.testing.assert_allclose(H_cascade, H1 * H2, atol=1e-12)
        # Both single-exp have DC gain = 1.0
        np.testing.assert_allclose(np.abs(H_cascade), 1.0, atol=1e-12)

    def test_step_response_vs_apply(self):
        """step_response equals apply(unit_step)."""
        s1 = SingleExponentialDistortion(amplitude=0.1, tau=30.0)
        cascade = CascadeDistortion(stages=[s1])
        dt = 0.2
        t = np.arange(0, 100, dt)

        step_via_method = cascade.step_response(t)
        step_via_apply = cascade.apply(self._step(t), dt)
        np.testing.assert_allclose(step_via_method, step_via_apply, atol=1e-12)

    def test_three_stage_cascade(self):
        """Three-stage cascade: IIR + single-exp + FIR."""
        s1 = SingleExponentialDistortion(amplitude=0.05, tau=40.0)
        s2 = IIRDistortion(
            b_coeffs=np.array([0.2, 0.3, 0.2]),
            a_coeffs=np.array([1.0, -0.3, 0.1]),
        )
        s3 = FIRDistortion(taps=np.array([0.5, 0.3, 0.2]))
        cascade = CascadeDistortion(stages=[s1, s2, s3])

        dt = 0.5
        t = np.arange(0, 100, dt)
        x = self._step(t)

        out_cascade = cascade.apply(x.copy(), dt)
        out_manual = s3.apply(s2.apply(s1.apply(x.copy(), dt), dt), dt)
        np.testing.assert_allclose(out_cascade, out_manual, atol=1e-12)

    def test_frequency_response_product(self):
        """Cascade frequency response is pointwise product."""
        s1 = SingleExponentialDistortion(amplitude=0.05, tau=40.0)
        s2 = FIRDistortion(taps=np.array([0.7, 0.3]))
        cascade = CascadeDistortion(stages=[s1, s2])

        dt = 0.5
        omega = 2 * np.pi * np.fft.fftfreq(256, d=dt)

        H_cascade = cascade.frequency_response(omega)
        H_product = s1.frequency_response(omega) * s2.frequency_response(omega)
        np.testing.assert_allclose(H_cascade, H_product, atol=1e-12)

    def test_impulse_response_approximate(self):
        """Impulse response applied to step ≈ step response derivative."""
        s1 = SingleExponentialDistortion(amplitude=0.05, tau=40.0)
        cascade = CascadeDistortion(stages=[s1])
        dt = 0.2
        t = np.arange(0, 100, dt)

        h = cascade.impulse_response(t)
        # cumulative sum of h * dt ≈ step response
        step_from_imp = np.cumsum(h) * dt
        step_direct = cascade.step_response(t)
        # Allow some error from delta approximation
        np.testing.assert_allclose(step_from_imp, step_direct, atol=1e-3)
