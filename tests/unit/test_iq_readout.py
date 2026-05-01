"""Unit tests for IQ readout models."""
from __future__ import annotations

import numpy as np
import pytest
from qutip import basis

from sqc.control.flux_signal import FluxSignal
from sqc.devices.transmon import TransmonQubit
from sqc.hardware.readout import IQReadoutModel, IQ_readout_legacy, IdealProjectiveReadout


pytestmark = pytest.mark.unit


def _make_qubit():
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000,
        T2=8000,
        n_levels=2,
    )


def test_ideal_projective_readout_ground_and_excited():
    q = _make_qubit()
    readout = IdealProjectiveReadout()

    assert readout.measure(q.state, q)["p_e"] == pytest.approx(0.0)
    q.state = basis(q.n_levels, 1)
    assert readout.measure(q.state, q)["p_e"] == pytest.approx(1.0)


def test_iq_readout_model_matches_legacy_helper():
    q = _make_qubit()
    phi = FluxSignal(type=1, t_list=np.linspace(0, 40, 80), amplitude=0.0)
    q.qubit_in_mag(phi, frame=1, omega_d=q.frequency)

    model = IQReadoutModel(tau=10.0, t_rabi=np.linspace(0, 4, 8))
    result = model.measure(q)
    legacy = IQ_readout_legacy(q, type=3, tau=10.0, t_rabi=np.linspace(0, 4, 8))

    assert result["p_e_I"] == pytest.approx(legacy[0])
    assert result["p_e_Q"] == pytest.approx(legacy[1])
