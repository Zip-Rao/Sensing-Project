"""Unit tests for IQReadoutModel and IdealProjectiveReadout."""
from __future__ import annotations

import numpy as np
import pytest
from qutip import basis


@pytest.fixture
def qubit():
    from src.qubit import TransmonQubit

    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=2,
    )


def test_ideal_projective_readout_ground():
    """IdealProjectiveReadout on |0> returns p_e=0."""
    from sqc.hardware.readout import IdealProjectiveReadout

    state = basis(2, 0)
    readout = IdealProjectiveReadout()
    result = readout.measure(state)
    assert "p_e" in result
    assert result["p_e"] == pytest.approx(0.0, abs=1e-12)


def test_ideal_projective_readout_excited():
    """IdealProjectiveReadout on |1> returns p_e=1."""
    from sqc.hardware.readout import IdealProjectiveReadout

    state = basis(2, 1)
    readout = IdealProjectiveReadout()
    result = readout.measure(state)
    assert "p_e" in result
    assert result["p_e"] == pytest.approx(1.0, abs=1e-12)


def test_ideal_projective_readout_superposition():
    """IdealProjectiveReadout on (|0>+|1>)/sqrt(2) returns p_e=0.5."""
    from sqc.hardware.readout import IdealProjectiveReadout

    state = (basis(2, 0) + basis(2, 1)).unit()
    readout = IdealProjectiveReadout()
    result = readout.measure(state)
    assert "p_e" in result
    assert result["p_e"] == pytest.approx(0.5, abs=1e-12)


def test_iq_readout_model_creation():
    """IQReadoutModel can be created with default parameters."""
    from sqc.hardware.readout import IQReadoutModel

    model = IQReadoutModel()
    assert model.tau == 20.0
    assert model.t_rabi is not None


def test_iq_readout_model_measure_cryoscope_type(qubit):
    """IQReadoutModel.measure works with a qubit having H_list populated."""
    from sqc.hardware.readout import IQReadoutModel
    from sqc.control.flux_signal import FluxSignal

    # Populate H_list via qubit_in_mag
    Phi = FluxSignal(type=0, t_list=np.linspace(0, 120, 240))
    qubit.qubit_in_mag(Phi, frame=1, omega_d=qubit.frequency)

    model = IQReadoutModel(tau=20, t_rabi=np.linspace(0, 10, 20))
    result = model.measure(qubit)

    assert "p_e_I" in result
    assert "p_e_Q" in result
    assert 0 <= result["p_e_I"] <= 1
    assert 0 <= result["p_e_Q"] <= 1


def test_iq_readout_legacy_function(qubit):
    """IQ_readout_legacy works as drop-in replacement."""
    from sqc.hardware.readout import IQ_readout_legacy
    from sqc.control.flux_signal import FluxSignal

    Phi = FluxSignal(type=0, t_list=np.linspace(0, 120, 240))
    qubit.qubit_in_mag(Phi, frame=1, omega_d=qubit.frequency)

    p_e_I, p_e_Q = IQ_readout_legacy(qubit, type=3, tau=20, t_rabi=np.linspace(0, 10, 20))
    assert isinstance(p_e_I, float)
    assert isinstance(p_e_Q, float)
    assert 0 <= p_e_I <= 1
    assert 0 <= p_e_Q <= 1
