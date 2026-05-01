"""Smoke tests: verify all src/ modules import cleanly."""
from __future__ import annotations

import pytest


pytestmark = pytest.mark.unit


def test_import_qubit():
    from src.qubit import TransmonQubit, Cavity, Coupled_System
    from src.qubit import ideal_iSWAP, ideal_CZ
    assert TransmonQubit is not None


def test_import_signal():
    from src.signal import Signal, CompositeSignal
    assert Signal is not None


def test_import_pulse():
    from src.pulse import (
        Pulse, CompositePulse,
        create_pulse, create_ramsey_pulse, create_diff_echo_pulse,
        create_echo_pulse, create_cpmg_pulse, create_cryoscope_pulse,
    )
    assert Pulse is not None


def test_import_protocal():
    from src.protocal import Protocal, Calibration, IQ_readout
    assert Protocal is not None


def test_import_analysis():
    from src.analysis import (
        Analysis,
        generate_basis_functions, basis_function_decomposition, R,
        forward_simulation, compute_jacobian,
        compute_jacobian_finite_difference, levenberg_marquardt,
    )
    assert Analysis is not None


def test_qubit_construction(qubit_default):
    """qubit_default fixture must construct without error."""
    assert qubit_default.frequency > 0
    assert qubit_default.anharmonicity < 0
    assert qubit_default.n_levels == 2
