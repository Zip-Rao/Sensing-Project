"""Unit tests for ControlLine.

Tests apply() and predistort() behavior with and without distortion.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import Waveform
from sqc.hardware.control_line import ControlLine
from sqc.hardware.distortion import SingleExponentialDistortion


@pytest.fixture
def t_default():
    return np.linspace(0, 200, 2000)


@pytest.fixture
def step_waveform(t_default):
    return Waveform(
        t_list=t_default,
        samples=np.where(t_default > 10, 1.0, 0.0),
    )


@pytest.fixture
def dist():
    return SingleExponentialDistortion(amplitude=0.05, tau=20.0)


class TestControlLine:
    """Tests for ControlLine."""

    # ---- construction ----

    def test_construction_minimal(self):
        """Minimal constructor with required fields."""
        line = ControlLine(
            name="Z0", kind="z", source="AWG0", target="Q0",
        )
        assert line.name == "Z0"
        assert line.kind == "z"
        assert line.transfer_function is None
        assert line.impedance == 50.0  # default

    def test_construction_with_distortion(self, dist):
        """Construction with transfer function."""
        line = ControlLine(
            name="XY0", kind="xy", source="AWG1", target="Q0",
            transfer_function=dist,
        )
        assert line.transfer_function is not None
        assert isinstance(line.transfer_function, SingleExponentialDistortion)

    def test_backward_compat_fields(self):
        """Legacy fields (impedance, attenuation, etc.) still work."""
        line = ControlLine(
            name="Z1", kind="z", source="AWG0", target="Q1",
            impedance=75.0, attenuation_db=30.0, delay=5.0,
            filter_type="lowpass", cutoff_freq=0.5,
        )
        assert line.impedance == 75.0
        assert line.attenuation_db == 30.0
        assert line.delay == 5.0
        assert line.filter_type == "lowpass"
        assert line.cutoff_freq == 0.5

    # ---- apply() ----

    def test_apply_no_distortion(self, step_waveform):
        """apply() without distortion returns a copy of the input."""
        line = ControlLine(
            name="Z0", kind="z", source="AWG0", target="Q0",
        )
        out = line.apply(step_waveform)
        np.testing.assert_array_equal(out.samples, step_waveform.samples)
        # Should be a different object
        assert out is not step_waveform

    def test_apply_with_distortion(self, step_waveform, dist):
        """apply() with distortion modifies the waveform."""
        line = ControlLine(
            name="Z0", kind="z", source="AWG0", target="Q0",
            transfer_function=dist,
        )
        out = line.apply(step_waveform)
        # Should be different from input
        assert not np.allclose(out.samples, step_waveform.samples, rtol=1e-14)

    def test_apply_preserves_shape(self, step_waveform, dist):
        """apply() preserves waveform length."""
        line = ControlLine(
            name="Z0", kind="z", source="AWG0", target="Q0",
            transfer_function=dist,
        )
        out = line.apply(step_waveform)
        assert out.samples.shape == step_waveform.samples.shape
        assert out.t_list.shape == step_waveform.t_list.shape

    def test_apply_with_delay(self, step_waveform):
        """apply() with delay > 0 shifts the waveform."""
        line = ControlLine(
            name="Z0", kind="z", source="AWG0", target="Q0",
            delay=1.0,
        )
        out = line.apply(step_waveform)
        # First sample should be zero (delayed)
        assert out.samples[0] == 0.0
        # Shape preserved
        assert out.samples.shape == step_waveform.samples.shape

    # ---- predistort() ----

    def test_predistort_returns_waveform(self, step_waveform, dist):
        """predistort() should return a Waveform."""
        from sqc.calibration.waveform import PredistortionDesigner

        line = ControlLine(
            name="Z0", kind="z", source="AWG0", target="Q0",
            transfer_function=dist,
        )
        designer = PredistortionDesigner(method="auto")
        predistorted = line.predistort(step_waveform, designer)
        assert isinstance(predistorted, Waveform)
        assert predistorted.samples.shape == step_waveform.samples.shape

    def test_predistort_improves_cascade(self, step_waveform, dist):
        """Predistortion should improve the forward model cascade."""
        from sqc.calibration.waveform import PredistortionDesigner

        line = ControlLine(
            name="Z0", kind="z", source="AWG0", target="Q0",
            transfer_function=dist,
        )
        designer = PredistortionDesigner(method="auto")

        # Without predistortion
        uncorrected = line.apply(step_waveform)
        rmse_uncorrected = np.sqrt(
            np.mean((uncorrected.samples - step_waveform.samples) ** 2)
        )

        # With predistortion
        predistorted = line.predistort(step_waveform, designer)
        corrected = line.apply(predistorted)
        rmse_corrected = np.sqrt(
            np.mean((corrected.samples - step_waveform.samples) ** 2)
        )

        improvement = rmse_uncorrected / max(rmse_corrected, 1e-30)
        assert improvement > 10, f"Improvement {improvement} <= 10"

    def test_kind_literals(self):
        """All valid kind values should be accepted."""
        for kind in ["xy", "z", "readout"]:
            line = ControlLine(
                name=f"line_{kind}", kind=kind,  # type: ignore[arg-type]
                source="AWG0", target="Q0",
            )
            assert line.kind == kind
