"""Unit tests for simulation runners."""
from __future__ import annotations

import numpy as np
import pytest
from qutip import basis, qeye

from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.devices.transmon import TransmonQubit
from sqc.simulation.runner import MesolveRunner, SlidingMeasurementRunner


pytestmark = pytest.mark.unit


def _make_qubit():
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000,
        T2=8000,
        n_levels=2,
    )


def test_mesolve_runner_returns_experiment_result():
    t_list = np.linspace(0, 2, 5)
    psi0 = basis(2, 0)
    result = MesolveRunner().run(
        qeye(2),
        psi0,
        t_list,
        e_ops=[basis(2, 0) * basis(2, 0).dag()],
    )

    assert result.axes["t"].shape == (5,)
    assert result.data["expect"].shape == (1, 5)
    np.testing.assert_allclose(result.data["expect"][0], np.ones(5))


def test_sliding_measurement_runner_returns_scan_and_population():
    q = _make_qubit()
    flux = FluxSignal(type=1, t_list=np.linspace(0, 4, 5), amplitude=0.0)
    pulse = create_ramsey_pulse(np.linspace(0, 2, 4), tau=0.0, omega_d=q.frequency)
    scan = np.linspace(0, 2, 3)

    result = SlidingMeasurementRunner().run(q, flux, pulse, scan_list=scan)

    np.testing.assert_allclose(result.axes["scan"], scan)
    assert result.data["p_e"].shape == (3,)
    assert np.all(np.isfinite(result.data["p_e"]))
