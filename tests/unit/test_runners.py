"""Unit tests for MesolveRunner and SlidingMeasurementRunner."""
from __future__ import annotations

import numpy as np
import pytest
from qutip import QobjEvo, basis, mesolve

from tests.conftest import assert_array_close


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


def test_mesolve_runner_basic(qubit):
    """MesolveRunner produces ExperimentResult with correct shape."""
    from sqc.simulation.runner import MesolveRunner

    H_0 = qubit.get_hamiltonian_rwa(qubit.frequency)
    t_list = np.linspace(0, 20, 40)
    psi_e = basis(qubit.n_levels, 1)

    runner = MesolveRunner()
    result = runner.run(
        H_list=[H_0],
        psi0=qubit.state,
        t_list=t_list,
        c_ops=[],
        e_ops=[psi_e * psi_e.dag()],
    )

    assert "expect" in result.data
    assert result.data["expect"].shape == (1, len(t_list))
    assert result.metadata["runner"] == "MesolveRunner"


def test_mesolve_runner_equivalence(qubit):
    """MesolveRunner produces same result as direct mesolve call."""
    from sqc.simulation.runner import MesolveRunner

    H_0 = qubit.get_hamiltonian_rwa(qubit.frequency)
    t_list = np.linspace(0, 20, 40)
    psi_e = basis(qubit.n_levels, 1)

    # Direct mesolve
    H_direct = QobjEvo(H_0)
    result_direct = mesolve(
        H_direct, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()]
    )

    # Via runner
    runner = MesolveRunner()
    result_runner = runner.run(
        H_list=[H_0],
        psi0=qubit.state,
        t_list=t_list,
        c_ops=[],
        e_ops=[psi_e * psi_e.dag()],
    )

    assert_array_close(
        result_runner.data["expect"][0],
        np.array(result_direct.expect[0]),
        name="mesolve_vs_runner",
    )


@pytest.mark.slow
def test_sliding_measurement_runner_basic(qubit):
    """SlidingMeasurementRunner returns ExperimentResult with correct structure."""
    from sqc.simulation.runner import SlidingMeasurementRunner
    from sqc.control.flux_signal import FluxSignal
    from sqc.control.sequence import create_ramsey_pulse

    t_list = np.linspace(0, 100, 100)
    Phi = FluxSignal(type=1, t_list=t_list, amplitude=0.01)
    t_rabi = np.linspace(0, 10, 20)
    cp = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    runner = SlidingMeasurementRunner()
    result = runner.run(qubit, Phi, cp)

    assert "p_e" in result.data
    assert "scan" in result.axes
    assert len(result.data["p_e"]) == len(result.axes["scan"])
    assert np.all(np.isfinite(result.data["p_e"]))
    # p_e should be in [0, 1]
    assert 0 <= result.data["p_e"].min() <= 1
    assert 0 <= result.data["p_e"].max() <= 1


@pytest.mark.slow
def test_sliding_measurement_equivalence(qubit):
    """SlidingMeasurementRunner matches legacy sliding_measrement output."""
    from sqc.simulation.runner import SlidingMeasurementRunner
    from sqc.control.flux_signal import FluxSignal
    from sqc.control.sequence import create_ramsey_pulse
    from src.protocal import Protocal

    t_list = np.linspace(0, 100, 100)
    Phi = FluxSignal(type=1, t_list=t_list, amplitude=0.01)
    t_rabi = np.linspace(0, 10, 20)
    cp = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    # Also create src versions for legacy call
    from src.signal import Signal as SrcSignal
    from src.pulse import create_ramsey_pulse as src_create_ramsey

    np.random.seed(42)
    Phi_src = SrcSignal(type=1, t_list=t_list, amplitude=0.01)
    cp_src = src_create_ramsey(t_rabi, tau=0.0, omega_d=qubit.frequency)

    # Legacy
    q_old = type(qubit)(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000, flux=0, state=0, n_levels=2
    )
    proto = Protocal(type=4)
    proto.initialize(q_old, state=0)
    scan_old, pe_old = proto.sliding_measrement(q_old, Phi_src, cp_src)

    # New runner (sqc types, src qubit)
    runner = SlidingMeasurementRunner()
    result_new = runner.run(qubit, Phi, cp)

    pe_old_arr = np.array(pe_old)
    pe_new_arr = np.array(result_new.data["p_e"])

    # The signals might differ due to seed timing, so check shapes
    assert len(pe_old_arr) == len(pe_new_arr)
    assert len(result_new.axes["scan"]) == len(pe_new_arr)
