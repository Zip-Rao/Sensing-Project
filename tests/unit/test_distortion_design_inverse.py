"""Unit tests for DistortionModel.design_inverse methods (P9.A.2)."""

import numpy as np
import pytest

from sqc.hardware.distortion import (
    CascadeDistortion,
    SingleExponentialDistortion,
)


class TestSingleExpDesignInverse:
    """P9.A.2 — SingleExponentialDistortion.design_inverse compensation."""

    @staticmethod
    def _step(t):
        return np.where(t >= 0, 1.0, 0.0)

    @pytest.mark.parametrize("formula", ["bilinear", "rol2020"])
    @pytest.mark.parametrize("amp,tau", [(0.05, 50.0), (0.10, 100.0), (0.02, 30.0)])
    def test_compensation_rmse(self, formula, amp, tau):
        """Inverse + forward cascade should recover unit step within 1e-3."""
        dt = 0.5
        dist = SingleExponentialDistortion(amplitude=amp, tau=tau)
        inv = dist.design_inverse(dt, formula=formula)
        cascade = CascadeDistortion(stages=[dist, inv])

        t = np.arange(0, 500, dt)
        step_in = self._step(t)
        step_out = cascade.apply(step_in, dt)

        rmse = float(np.sqrt(np.mean((step_out - step_in) ** 2)))
        assert rmse < 1e-3, f"{formula} RMSE={rmse:.2e} exceeds 1e-3"

    def test_bilinear_vs_rol2020_within_5pct(self):
        """bilinear and rol2020 produce similar inverse filters."""
        dt = 0.5
        dist = SingleExponentialDistortion(amplitude=0.05, tau=50.0)
        inv_b = dist.design_inverse(dt, formula="bilinear")
        inv_r = dist.design_inverse(dt, formula="rol2020")

        # Compare frequency responses
        omega = 2 * np.pi * np.fft.fftfreq(512, d=dt)
        H_b = inv_b.frequency_response(omega)
        H_r = inv_r.frequency_response(omega)

        diff = np.abs(H_b - H_r) / np.maximum(np.abs(H_b), 1e-12)
        assert np.max(diff) < 0.05, f"max relative diff {np.max(diff):.4f} >= 5%"

    def test_warning_for_large_amplitude(self):
        """|A| > 0.5 should trigger a warning."""
        import warnings

        dist = SingleExponentialDistortion(amplitude=0.6, tau=50.0)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            dist.design_inverse(0.5)
            assert len(w) > 0, "No warning raised for |A| > 0.5"

    @pytest.mark.parametrize("formula", ["bilinear", "rol2020"])
    def test_amplitude_near_one_returns_stable_inverse(self, formula):
        """A -> 1 collapses the forward DC gain; the inverse must degrade to
        a stable (identity) filter rather than an unbounded one."""
        dt = 0.5
        dist = SingleExponentialDistortion(amplitude=1.0, tau=50.0)
        inv = dist.design_inverse(dt, formula=formula)

        # Identity fallback: b == a == [1.0], and a finite frequency response.
        np.testing.assert_allclose(inv.b_coeffs, [1.0], atol=1e-12)
        np.testing.assert_allclose(inv.a_coeffs, [1.0], atol=1e-12)

        omega = 2 * np.pi * np.fft.fftfreq(256, d=dt)
        H = inv.frequency_response(omega)
        assert np.all(np.isfinite(H)), "inverse frequency response is not finite"

    def test_bilinear_matches_existing_wrapper(self):
        """Output matches PredistortionDesigner._single_exp_to_iir_inverse."""
        from sqc.calibration.waveform import PredistortionDesigner

        dt = 0.5
        dist = SingleExponentialDistortion(amplitude=0.05, tau=50.0)

        inv_new = dist.design_inverse(dt, formula="bilinear")
        inv_old = PredistortionDesigner._single_exp_to_iir_inverse(dist, dt)

        np.testing.assert_allclose(inv_new.b_coeffs, inv_old.b_coeffs, atol=1e-12)
        np.testing.assert_allclose(inv_new.a_coeffs, inv_old.a_coeffs, atol=1e-12)
