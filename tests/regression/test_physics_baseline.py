"""Physics regression tests for frozen legacy behavior."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from tests.conftest import assert_array_close, load_baseline


pytestmark = pytest.mark.regression
plt.show = lambda *args, **kwargs: None


def _make_qubit():
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


def test_qubit_static_properties(qubit_default):
    bl = load_baseline("qubit_static")
    assert qubit_default.frequency == pytest.approx(bl["frequency"], rel=1e-12)
    assert qubit_default.anharmonicity == pytest.approx(
        bl["anharmonicity"], rel=1e-12
    )
    assert qubit_default.frequency_sensitivity(0.0) == pytest.approx(
        bl["sensitivity_at_zero"], rel=1e-9
    )


def test_ramsey_default_baseline():
    from src.protocal import Protocal

    q = _make_qubit()
    proto = Protocal(type=1)
    proto.initialize(q, state=0)
    phi, tau_list, p_e_list = proto.evolve(q)
    plt.close("all")

    bl = load_baseline("ramsey_default")
    assert_array_close(phi.signal, bl["Phi_signal"], name="Phi.signal")
    assert_array_close(np.asarray(tau_list), bl["tau_list"], name="tau_list")
    assert_array_close(np.asarray(p_e_list), bl["p_e_list"], name="p_e_list")


def test_diff_echo_default_baseline():
    from src.protocal import Protocal

    q = _make_qubit()
    proto = Protocal(type=2)
    proto.initialize(q, state=0)
    _phi, _tau_list, p_e_list, k, t_int = proto.evolve(q)
    plt.close("all")

    bl = load_baseline("diff_echo_default")
    assert k == bl["k"]
    assert t_int == pytest.approx(bl["t_int"], rel=1e-12)
    assert_array_close(np.asarray(p_e_list), bl["p_e_list"], name="p_e_list")


def test_transient_default_baseline():
    from src.protocal import Protocal

    q = _make_qubit()
    proto = Protocal(type=4)
    proto.initialize(q, state=0)
    t_samples, kernel, scan_list, delta_p, p_e, _phi, _ctrl = proto.evolve(q)
    plt.close("all")

    bl = load_baseline("transient_default")
    assert_array_close(t_samples, bl["t_samples"], name="t_samples")
    assert_array_close(kernel, bl["kernel"], name="kernel")
    assert_array_close(np.asarray(scan_list), bl["scan_list"], name="scan_list")
    assert_array_close(np.asarray(delta_p), bl["delta_p"], name="delta_p")
    assert_array_close(np.asarray(p_e), bl["p_e"], name="p_e")


def test_cryoscope_default_baseline():
    from src.protocal import Protocal

    q = _make_qubit()
    proto = Protocal(type=5)
    proto.initialize(q, state=0)
    trunc_list, varphi, phi, p_e_list = proto.evolve(q)
    plt.close("all")

    bl = load_baseline("cryoscope_default")
    assert_array_close(np.asarray(trunc_list), bl["trunc_list"], name="trunc_list")
    assert_array_close(np.asarray(varphi), bl["varphi"], name="varphi")
    assert_array_close(np.asarray(phi.signal), bl["Phi_signal"], name="Phi.signal")
    assert_array_close(np.asarray(p_e_list[0]), bl["p_e_I"], name="p_e_I")
    assert_array_close(np.asarray(p_e_list[1]), bl["p_e_Q"], name="p_e_Q")
