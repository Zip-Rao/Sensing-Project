"""Unit tests for P4 control-line distortion models."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import Waveform
from sqc.hardware.distortion import (
    CustomTransferDistortion,
    FIRDistortion,
    IIRDistortion,
    MultiExponentialDistortion,
    SingleExponentialDistortion,
)


pytestmark = pytest.mark.unit


def _time_grid():
    return np.linspace(0.0, 100.0, 1001)


def test_single_exponential_step_settles_to_dc_gain():
    t = _time_grid()
    waveform = Waveform(t_list=t, samples=np.where(t > 10.0, 1.0, 0.0))
    model = SingleExponentialDistortion(amplitude=0.05, tau=20.0)

    out = model.apply(waveform)

    assert out.samples.shape == waveform.samples.shape
    assert out.samples[-1] == pytest.approx(1.0, abs=1e-3)
    assert model.frequency_response(np.array([0.0]))[0] == pytest.approx(1.0)
    assert np.sum(model.impulse_response(t)) == pytest.approx(1.0, abs=1e-3)


def test_multi_exponential_has_unit_dc_gain():
    t = _time_grid()
    model = MultiExponentialDistortion(
        amplitudes=np.array([0.03, 0.02]),
        taus=np.array([8.0, 30.0]),
    )

    step = model.step_response(t)

    assert step[-1] == pytest.approx(1.0, abs=1e-3)
    assert model.frequency_response(np.array([0.0]))[0] == pytest.approx(1.0)
    assert np.sum(model.impulse_response(t)) == pytest.approx(1.0, abs=1e-3)


def test_fir_distortion_uses_discrete_taps():
    t = np.arange(0.0, 20.0, 0.5)
    model = FIRDistortion(taps=np.array([0.5, 0.3, 0.2]), dt=0.5)

    step = model.step_response(t)

    assert step[-1] == pytest.approx(1.0)
    assert model.frequency_response(np.array([0.0]))[0] == pytest.approx(1.0)


def test_iir_distortion_matches_single_exp_coefficients():
    t = _time_grid()
    waveform = Waveform(t_list=t, samples=np.sin(0.1 * t))
    single = SingleExponentialDistortion(amplitude=0.1, tau=12.0)
    b, a = single.iir_coefficients(t[1] - t[0])
    generic = IIRDistortion(b=b, a=a, dt=t[1] - t[0])

    assert np.allclose(generic.apply(waveform).samples, single.apply(waveform).samples)


def test_custom_transfer_scales_waveform_in_frequency_domain():
    t = _time_grid()
    waveform = Waveform(t_list=t, samples=np.sin(0.2 * t))
    omega = np.linspace(-100.0, 100.0, 401)
    model = CustomTransferDistortion(omega=omega, H=np.full_like(omega, 0.5))

    out = model.apply(waveform)

    assert np.max(np.abs(out.samples - 0.5 * waveform.samples)) < 1e-10


def test_distortion_rejects_non_uniform_time_grid():
    waveform = Waveform(
        t_list=np.array([0.0, 1.0, 2.2]),
        samples=np.array([0.0, 1.0, 0.0]),
    )
    model = SingleExponentialDistortion(amplitude=0.1, tau=5.0)

    with pytest.raises(ValueError):
        model.apply(waveform)
