"""Legacy API checks kept while src facade migration is deferred."""
from __future__ import annotations

import numpy as np
import pytest


pytestmark = pytest.mark.integration


def _make_qubit():
    from src.qubit import TransmonQubit

    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000,
        T2=8000,
        n_levels=2,
    )


def test_legacy_protocol_type_1_tuple_shape_still_available():
    from src.protocal import Protocal

    q = _make_qubit()
    proto = Protocal(type=1)
    proto.initialize(q, state=0)
    phi, tau_list, p_e_list = proto.evolve(q)

    assert hasattr(phi, "signal")
    assert len(tau_list) == len(p_e_list)
