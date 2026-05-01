"""Physics regression: refactored code must reproduce frozen baselines.

These tests are SLOW (each runs the full Protocal.evolve). Mark with
@pytest.mark.regression so they can be skipped during fast iterations:

    pytest -m "not regression"   # skip
    pytest -m regression         # only regression
"""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import (
    assert_array_close,
    load_baseline,
)


pytestmark = pytest.mark.regression


def test_qubit_static_properties(qubit_default):
    bl = load_baseline("qubit_static")
    assert qubit_default.frequency == pytest.approx(
        bl["frequency"], rel=1e-12
    )
    assert qubit_default.anharmonicity == pytest.approx(
        bl["anharmonicity"], rel=1e-12
    )
    assert qubit_default.frequency_sensitivity(0.0) == pytest.approx(
        bl["sensitivity_at_zero"], rel=1e-9
    )


def test_ramsey_default_baseline():
    """Protocal(type=1).evolve must match frozen Ramsey output."""
    from src.qubit import TransmonQubit
    from src.protocal import Protocal

    q = TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, state=0, n_levels=2,
    )
    proto = Protocal(type=1)
    proto.initialize(q, state=0)
    Phi, tau_list, p_e_list = proto.evolve(q)

    bl = load_baseline("ramsey_default")
    assert_array_close(Phi.signal, bl["Phi_signal"], name="Phi.signal")
    assert_array_close(np.asarray(tau_list), bl["tau_list"], name="tau_list")
    assert_array_close(np.asarray(p_e_list), bl["p_e_list"], name="p_e_list")


def test_diff_echo_default_baseline():
    from src.qubit import TransmonQubit
    from src.protocal import Protocal

    q = TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, state=0, n_levels=2,
    )
    proto = Protocal(type=2)
    proto.initialize(q, state=0)
    Phi, tau_list, p_e_list, k, t_int = proto.evolve(q)

    bl = load_baseline("diff_echo_default")
    assert k == bl["k"]
    assert t_int == pytest.approx(bl["t_int"], rel=1e-12)
    assert_array_close(np.asarray(p_e_list), bl["p_e_list"], name="p_e_list")


def test_transient_default_baseline():
    from src.qubit import TransmonQubit
    from src.protocal import Protocal

    q = TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, state=0, n_levels=2,
    )
    proto = Protocal(type=4)
    proto.initialize(q, state=0)
    t_samples, kernel, scan_list, delta_p, p_e, Phi, ctrl = proto.evolve(q)

    bl = load_baseline("transient_default")
    assert_array_close(kernel, bl["kernel"], name="kernel")
    assert_array_close(np.asarray(delta_p), bl["delta_p"], name="delta_p")
    assert_array_close(np.asarray(p_e), bl["p_e"], name="p_e")


def test_lm_default_baseline():
    """LM numerical inversion: b_opt matches frozen baseline."""
    from src.qubit import TransmonQubit
    from src.pulse import create_ramsey_pulse

    bl = load_baseline("lm_default")

    q = TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, state=0, n_levels=2,
    )
    cp = create_ramsey_pulse(
        t_rabi=np.linspace(0, 5, 6), tau=5.0, omega_d=q.frequency,
    )

    from sqc.reconstruction.numerical_inverse import LMReconstruction
    from sqc.simulation.result import ExperimentResult

    recon = LMReconstruction(
        qubit=q,
        control_pulse=cp,
        basis_type=bl["basis_type"],
        n_basis=bl["n_basis"],
        lambda_reg=100.0,
        max_iter=2,
        tol=1e-3,
    )
    meas = ExperimentResult(
        data={"p_meas": bl["p_meas"]},
        axes={"t_signal": bl["t_list"]},
    )
    B_opt_flux, history = recon.reconstruct(meas)
    b_new = B_opt_flux.params["b"]

    assert_array_close(b_new, bl["b_opt"], name="lm_b_opt")
    assert_array_close(
        np.asarray(history["res"][-1]),
        bl["history_res_final"],
        name="lm_res_final",
    )
