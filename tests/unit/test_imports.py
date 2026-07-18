"""Smoke tests: verify all sqc/ and src_mirror/ modules import cleanly.

Tests: all ABCs, all data structures, mirror equivalence, legacy src/ imports.
"""
from __future__ import annotations

import pytest


class TestSqcImports:
    """All sqc sub-packages must be importable."""

    def test_import_sqc(self):
        import sqc
        assert sqc.__version__ == "1.0.0"

    def test_import_devices(self):
        from sqc.devices.base import Device
        from sqc.devices.transmon import QubitSpec, TransmonQubit
        from sqc.devices.resonator import Resonator
        from sqc.devices.chip import CoupledSystem
        assert Device is not None
        assert QubitSpec is not None

    def test_import_hardware(self):
        from sqc.hardware.distortion import DistortionModel
        from sqc.hardware.readout import ReadoutModel
        from sqc.hardware.control_line import ControlLine
        from sqc.hardware.transfer_matrix import TransferMatrix
        assert DistortionModel is not None

    def test_import_control(self):
        from sqc.control.waveform import Waveform, CompositeWaveform
        from sqc.control.flux_signal import FluxSignal, Signal, CompositeSignal
        from sqc.control.pulse import PulseBase, Pulse, CompositePulse
        from sqc.control.sequence import (
            PulseSequence,
            create_pulse, create_ramsey_pulse, create_diff_echo_pulse,
            create_echo_pulse, create_cpmg_pulse, create_cryoscope_pulse,
        )
        from sqc.control.gates import ideal_iSWAP, ideal_CZ
        assert Pulse is not None

    def test_import_simulation(self):
        from sqc.simulation.runner import (
            RunnerBase, MesolveRunner, SlidingMeasurementRunner,
        )
        from sqc.simulation.hamiltonian import HamiltonianBuilder
        from sqc.simulation.noise import generate_1f_noise
        from sqc.simulation.result import (
            ExperimentResult, MeasurementTrace,
            extract_expectation, extract_population,
        )
        assert HamiltonianBuilder is not None

    def test_import_experiments(self):
        from sqc.experiments.base import Experiment
        assert Experiment is not None

    def test_import_calibration(self):
        from sqc.calibration.base import Calibration, CalibrationTable
        assert Calibration is not None

    def test_import_reconstruction(self):
        from sqc.reconstruction.base import Reconstruction
        from sqc.reconstruction.kernel import KernelEstimator
        assert Reconstruction is not None

    def test_import_workflows(self):
        from sqc.workflows.base import Workflow
        assert Workflow is not None


class TestMirrorImports:
    """src_mirror/ aliases must point to the same sqc objects."""

    def test_mirror_qubit(self):
        import src_mirror.qubit as m
        import sqc.devices.transmon as t
        assert m.TransmonQubit is t.TransmonQubit
        assert m.QubitSpec is t.QubitSpec

    def test_mirror_signal(self):
        import src_mirror.signal as m
        import sqc.control.flux_signal as f
        assert m.Signal is f.FluxSignal
        assert m.Signal is f.Signal

    def test_mirror_pulse(self):
        import src_mirror.pulse as m
        import sqc.control.pulse as p
        assert m.Pulse is p.Pulse
        assert m.CompositePulse is p.CompositePulse

    def test_legacy_src_imports(self):
        """src/ imports still work (src/ is untouched)."""
        from src.qubit import TransmonQubit, Cavity, Coupled_System
        from src.signal import Signal, CompositeSignal
        from src.pulse import (
            Pulse, CompositePulse,
            create_pulse, create_ramsey_pulse,
        )
        assert TransmonQubit is not None
        assert Signal is not None
        assert Pulse is not None

    def test_transmon_qubit_instantiates(self):
        """Both old and new TransmonQubit construct correctly."""
        import numpy as np
        from src.qubit import TransmonQubit as OldQ
        from sqc.devices.transmon import TransmonQubit as NewQ

        qo = OldQ(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000, n_levels=2)
        qn = NewQ(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000, n_levels=2)
        assert qo.frequency == pytest.approx(qn.frequency, rel=1e-12)
        assert qo.anharmonicity == pytest.approx(qn.anharmonicity, rel=1e-12)
