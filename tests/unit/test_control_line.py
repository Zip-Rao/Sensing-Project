"""Unit tests for physical control-line behavior."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.calibration.predistortion import PredistortionDesigner
from sqc.control.waveform import Waveform
from sqc.hardware.control_line import ControlLine
from sqc.hardware.distortion import SingleExponentialDistortion


pytestmark = pytest.mark.unit


def test_control_line_without_transfer_returns_copy():
    t = np.linspace(0.0, 10.0, 101)
    waveform = Waveform(t_list=t, samples=np.sin(t), metadata={"name": "target"})
    line = ControlLine(name="Z0", kind="z", source="AWG0", target="Q0")

    out = line.apply(waveform)

    assert out is not waveform
    assert np.array_equal(out.samples, waveform.samples)
    assert out.metadata == waveform.metadata


def test_control_line_applies_transfer_function():
    t = np.linspace(0.0, 100.0, 1001)
    waveform = Waveform(t_list=t, samples=np.where(t > 10.0, 1.0, 0.0))
    transfer = SingleExponentialDistortion(amplitude=0.05, tau=20.0)
    line = ControlLine(
        name="Z0",
        kind="z",
        source="AWG0",
        target="Q0",
        transfer_function=transfer,
    )

    out = line.apply(waveform)

    assert out.samples.shape == waveform.samples.shape
    assert out.samples[-1] == pytest.approx(1.0, abs=1e-3)


def test_control_line_predistort_inverts_standard_iir_transfer():
    t = np.linspace(0.0, 80.0, 801)
    target = Waveform(t_list=t, samples=0.5 * (np.tanh((t - 20.0) / 2.0) + 1.0))
    transfer = SingleExponentialDistortion(amplitude=0.1, tau=10.0)
    line = ControlLine(
        name="Z0",
        kind="z",
        source="AWG0",
        target="Q0",
        transfer_function=transfer,
    )
    designer = PredistortionDesigner(method="iir_inverse")

    awg = line.predistort(target, designer)
    corrected = line.apply(awg)

    assert np.max(np.abs(corrected.samples - target.samples)) < 1e-10
