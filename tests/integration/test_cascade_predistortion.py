"""Integration test: cascade predistortion (P9.A.3).

Verifies that MultiExponentialDistortion.design_inverse() → CascadeDistortion
outperforms frequency_inverse for multi-exponential distortion.
"""

import numpy as np
import pytest

from sqc.hardware.distortion import (
    CascadeDistortion,
    CustomTransferDistortion,
    MultiExponentialDistortion,
)


def _step(t):
    return np.where(t >= 0, 1.0, 0.0)


def _rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


@pytest.mark.parametrize("formula", ["bilinear", "rol2020"])
def test_cascade_beats_frequency_inverse_2x(formula):
    """Cascade inverse RMSE ≤ frequency_inverse RMSE / 2."""
    dt = 0.5
    dist = MultiExponentialDistortion(
        amplitudes=np.array([0.06, 0.03, 0.01]),
        taus=np.array([200.0, 50.0, 10.0]),
    )

    # Cascade inverse
    inv_cascade = dist.design_inverse(
        dt, formula=formula, n_iir_stages=3, fir_taps=72,
        fir_threshold_ns=30.0,
    )
    assert isinstance(inv_cascade, CascadeDistortion)
    assert len(inv_cascade.stages) >= 2  # at least 2 IIR for τ=200, τ=50

    # Frequency inverse
    n_fft = 4096
    omega = 2.0 * np.pi * np.fft.fftfreq(n_fft, d=dt)
    H = dist.frequency_response(omega)
    H_inv = np.conj(H) / (np.abs(H) ** 2 + 1e-4)
    inv_freq = CustomTransferDistortion(
        omega_grid=omega.copy(), H_grid=H_inv.copy(),
    )

    t = np.arange(0, 400, dt)
    step_in = _step(t)

    cascade_corr = CascadeDistortion(stages=[dist, inv_cascade])
    out_cascade = cascade_corr.apply(step_in, dt)
    rmse_cascade = _rmse(out_cascade, step_in)

    freq_corr = CascadeDistortion(stages=[dist, inv_freq])
    out_freq = freq_corr.apply(step_in, dt)
    rmse_freq = _rmse(out_freq, step_in)

    assert rmse_cascade < rmse_freq / 2.0, (
        f"Cascade RMSE {rmse_cascade:.2e} not better than "
        f"freq_inv RMSE {rmse_freq:.2e} / 2 = {rmse_freq/2:.2e}"
    )


def test_cascade_no_iir_fallback_to_fir():
    """When no τ > threshold, falls back to FIR-only."""
    dt = 0.5
    dist = MultiExponentialDistortion(
        amplitudes=np.array([0.02]),
        taus=np.array([5.0]),  # fast — below 30ns threshold
    )
    inv = dist.design_inverse(
        dt, n_iir_stages=3, fir_taps=20, fir_threshold_ns=30.0,
    )
    assert isinstance(inv, CascadeDistortion)
    # Should have FIR stage even though no IIR stages qualified
    assert len(inv.stages) > 0


def test_cascade_empty_model():
    """Empty MultiExponentialDistortion returns identity cascade."""
    dist = MultiExponentialDistortion(
        amplitudes=np.array([]), taus=np.array([]),
    )
    inv = dist.design_inverse(0.5)
    assert isinstance(inv, CascadeDistortion)
    assert len(inv.stages) == 0


def test_predistortion_designer_auto_uses_cascade():
    """PredistortionDesigner(method='auto') delegates to model.design_inverse."""
    from sqc.calibration.waveform import PredistortionDesigner

    dt = 0.5
    dist = MultiExponentialDistortion(
        amplitudes=np.array([0.06, 0.03]),
        taus=np.array([200.0, 50.0]),
    )
    designer = PredistortionDesigner(method="auto")
    inv = designer.design(dist, dt=dt)
    assert isinstance(inv, CascadeDistortion)

    # Verify correction quality
    cascade = CascadeDistortion(stages=[dist, inv])
    t = np.arange(0, 300, dt)
    step_in = _step(t)
    out = cascade.apply(step_in, dt)
    rmse_val = _rmse(out, step_in)
    assert rmse_val < 0.01, f"RMSE {rmse_val:.2e} excessive"
