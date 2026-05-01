"""Import and mirror checks for Phase 1."""
from __future__ import annotations

import numpy as np
import pytest


pytestmark = pytest.mark.unit


def test_all_phase1_abcs_import():
    from sqc.calibration.base import Calibration, CalibrationTable
    from sqc.control.pulse import PulseBase
    from sqc.devices.base import Device
    from sqc.experiments.base import Experiment
    from sqc.hardware.distortion import DistortionModel
    from sqc.hardware.readout import ReadoutModel
    from sqc.reconstruction.base import Reconstruction
    from sqc.simulation.runner import MesolveRunner, RunnerBase, SlidingMeasurementRunner
    from sqc.workflows.base import Workflow

    assert Device is not None
    assert DistortionModel is not None
    assert ReadoutModel is not None
    assert PulseBase is not None
    assert RunnerBase is not None
    assert MesolveRunner is not None
    assert SlidingMeasurementRunner is not None
    assert Experiment is not None
    assert Calibration is not None
    assert CalibrationTable("Q0", "kind", np.array([0]), np.array([0])) is not None
    assert Reconstruction is not None
    assert Workflow is not None


def test_src_qubit_mirrors_sqc():
    import src.qubit as legacy
    import sqc.devices.transmon as native

    assert legacy.TransmonQubit is native.TransmonQubit


def test_src_signal_mirrors_sqc():
    import src.signal as legacy
    import sqc.control.flux_signal as native

    assert legacy.Signal is native.FluxSignal
    assert legacy.CompositeSignal is native.CompositeSignal


def test_src_pulse_mirrors_sqc():
    import src.pulse as legacy
    import sqc.control.pulse as native

    assert legacy.Pulse is native.Pulse
    assert legacy.CompositePulse is native.CompositePulse
