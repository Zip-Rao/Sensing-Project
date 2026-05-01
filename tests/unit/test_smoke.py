"""Smoke tests for the legacy src package."""
from __future__ import annotations

import pytest


pytestmark = pytest.mark.unit


def test_import_qubit():
    from src.qubit import Cavity, Coupled_System, TransmonQubit
    from src.qubit import ideal_CZ, ideal_iSWAP

    assert TransmonQubit is not None
    assert Cavity is not None
    assert Coupled_System is not None
    assert ideal_iSWAP is not None
    assert ideal_CZ is not None


def test_import_signal():
    from src.signal import CompositeSignal, Signal

    assert Signal is not None
    assert CompositeSignal is not None


def test_import_pulse():
    from src.pulse import CompositePulse, Pulse
    from src.pulse import create_cpmg_pulse, create_cryoscope_pulse
    from src.pulse import create_diff_echo_pulse, create_echo_pulse
    from src.pulse import create_pulse, create_ramsey_pulse

    assert Pulse is not None
    assert CompositePulse is not None
    assert create_pulse is not None
    assert create_ramsey_pulse is not None
    assert create_diff_echo_pulse is not None
    assert create_echo_pulse is not None
    assert create_cpmg_pulse is not None
    assert create_cryoscope_pulse is not None


def test_import_protocal():
    from src.protocal import IQ_readout, Calibration, Protocal

    assert Protocal is not None
    assert Calibration is not None
    assert IQ_readout is not None


def test_import_analysis():
    from src.analysis import Analysis
    from src.analysis import R, basis_function_decomposition
    from src.analysis import compute_jacobian, compute_jacobian_finite_difference
    from src.analysis import forward_simulation, generate_basis_functions
    from src.analysis import levenberg_marquardt

    assert Analysis is not None
    assert generate_basis_functions is not None
    assert basis_function_decomposition is not None
    assert R is not None
    assert forward_simulation is not None
    assert compute_jacobian is not None
    assert compute_jacobian_finite_difference is not None
    assert levenberg_marquardt is not None


def test_qubit_construction(qubit_default):
    """The default qubit fixture constructs and exposes basic properties."""
    assert qubit_default.frequency > 0
    assert qubit_default.anharmonicity < 0
    assert qubit_default.n_levels == 2
