"""Unit tests for HamiltonianBuilder."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.flux_signal import FluxSignal
from sqc.devices.transmon import TransmonQubit
from sqc.simulation.hamiltonian import HamiltonianBuilder


pytestmark = pytest.mark.unit


def test_builder_matches_legacy_qubit_in_mag_rotating_frame():
    q = TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000,
        T2=8000,
        n_levels=2,
    )
    phi = FluxSignal(
        type=2,
        t_list=np.linspace(0, 20, 50),
        amplitude=0.001,
        frequency=0.01,
        noise_level=0,
    )
    q.qubit_in_mag(phi, frame=1, omega_d=q.frequency)
    h_list, t_list = HamiltonianBuilder.build(
        q, flux_signal=phi, frame="rotating", omega_d=q.frequency
    )
    np.testing.assert_allclose(t_list, phi.t_list)
    np.testing.assert_allclose(h_list[1][1], q.H_list[1][1])


def test_builder_lab_frame_returns_frequency_coefficients():
    q = TransmonQubit(2 * np.pi * 0.2, 2 * np.pi * 15, 10000, 8000, n_levels=2)
    phi = FluxSignal(type=1, t_list=np.linspace(0, 5, 6), amplitude=0.0)
    h_list, _ = HamiltonianBuilder.build(q, flux_signal=phi, frame="lab")
    np.testing.assert_allclose(h_list[1][1], np.full(6, q.frequency))


def test_builder_accepts_qubit_spec():
    q = TransmonQubit(2 * np.pi * 0.2, 2 * np.pi * 15, 10000, 8000, n_levels=2)
    h_list, t_list = HamiltonianBuilder.build(q.spec(), frame="rotating")
    assert len(t_list) == 100
    assert len(h_list) == 2
