"""Unit tests for sqc.devices.transmon.QubitSpec."""
from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from sqc.devices.transmon import QubitSpec, TransmonQubit


pytestmark = pytest.mark.unit


def test_qubit_spec_frequency_matches_legacy_default():
    spec = QubitSpec("Q0", 2 * np.pi * 0.2, 2 * np.pi * 15, 10000, 8000)
    legacy = TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000,
        T2=8000,
        n_levels=2,
    )
    assert spec.frequency() == pytest.approx(legacy.frequency)


def test_qubit_spec_is_frozen():
    spec = QubitSpec("Q0", 1.0, 2.0, 10000, 8000)
    with pytest.raises(FrozenInstanceError):
        spec.EC = 2.0


def test_optimal_work_point_unit_convention():
    assert QubitSpec.optimal_work_point() == pytest.approx(
        np.arctan(np.sqrt(2)) / np.pi
    )
    legacy = TransmonQubit(1.0, 20.0, 10000, 8000, n_levels=2)
    assert legacy.optimal_work_point() == pytest.approx(np.arctan(np.sqrt(2)))


def test_with_flux_returns_new_spec():
    spec = QubitSpec("Q0", 1.0, 20.0, 10000, 8000)
    shifted = spec.with_flux(0.1)
    assert shifted is not spec
    assert shifted.flux_bias == pytest.approx(0.1)
    assert spec.flux_bias == pytest.approx(0.0)
