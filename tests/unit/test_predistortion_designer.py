"""Unit tests for predistortion filter design."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.calibration.predistortion import PredistortionDesigner
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import (
    CustomTransferDistortion,
    FIRDistortion,
    IIRDistortion,
    SingleExponentialDistortion,
)


pytestmark = pytest.mark.unit


def _target_waveform() -> Waveform:
    t = np.linspace(0.0, 120.0, 1024)
    samples = 0.5 * (np.tanh((t - 20.0) / 2.0) - np.tanh((t - 80.0) / 2.0))
    return Waveform(t_list=t, samples=samples)


def test_iir_inverse_exactly_corrects_single_exponential_model():
    target = _target_waveform()
    transfer = SingleExponentialDistortion(amplitude=0.15, tau=12.0)
    designer = PredistortionDesigner(method="iir_inverse")

    awg = designer.predistort(target, transfer)
    corrected = transfer.apply(awg)

    assert isinstance(designer.design(transfer, dt=target.t_list[1] - target.t_list[0]), IIRDistortion)
    assert np.max(np.abs(corrected.samples - target.samples)) < 1e-10


def test_frequency_inverse_returns_custom_transfer_and_improves_rmse():
    target = _target_waveform()
    transfer = SingleExponentialDistortion(amplitude=0.15, tau=12.0)
    designer = PredistortionDesigner(method="frequency_inverse", regularization=1e-5)

    inverse = designer.design(transfer, dt=target.t_list[1] - target.t_list[0])
    corrected = transfer.apply(designer.predistort(target, transfer))
    uncorrected = transfer.apply(target)

    rmse_uncorrected = np.sqrt(np.mean((uncorrected.samples - target.samples) ** 2))
    rmse_corrected = np.sqrt(np.mean((corrected.samples - target.samples) ** 2))
    assert isinstance(inverse, CustomTransferDistortion)
    assert rmse_corrected < rmse_uncorrected / 10.0


def test_fir_inverse_returns_fir_and_improves_rmse_with_enough_taps():
    target = _target_waveform()
    transfer = SingleExponentialDistortion(amplitude=0.10, tau=4.0)
    designer = PredistortionDesigner(
        method="fir_inverse",
        n_taps=128,
        n_freq=4096,
        regularization=1e-6,
    )

    inverse = designer.design(transfer, dt=target.t_list[1] - target.t_list[0])
    corrected = transfer.apply(designer.predistort(target, transfer))
    uncorrected = transfer.apply(target)

    rmse_uncorrected = np.sqrt(np.mean((uncorrected.samples - target.samples) ** 2))
    rmse_corrected = np.sqrt(np.mean((corrected.samples - target.samples) ** 2))
    assert isinstance(inverse, FIRDistortion)
    assert rmse_corrected < rmse_uncorrected


def test_predistortion_without_transfer_returns_copy():
    target = _target_waveform()
    out = PredistortionDesigner().predistort(target, transfer=None)

    assert out is not target
    assert np.array_equal(out.samples, target.samples)
