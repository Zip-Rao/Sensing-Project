"""Unit tests for FluxSignal legacy compatibility."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.flux_signal import CompositeSignal, FluxSignal, Signal


pytestmark = pytest.mark.unit


def test_flux_signal_alias():
    assert Signal is FluxSignal


def test_type_2_signal_values_are_stable():
    t = np.linspace(0, 10, 11)
    sig = FluxSignal(type=2, t_list=t, amplitude=0.5, frequency=0.1, noise_level=0)
    expected = 0.5 * np.sin(2 * np.pi * 0.1 * t)
    np.testing.assert_allclose(sig.signal, expected)
    np.testing.assert_allclose(sig.samples, expected)


def test_type_8_custom_signal_and_update():
    t = np.linspace(0, 2, 3)
    sig = FluxSignal(type=8, t_list=t, signal=np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(sig.signal, np.array([1.0, 2.0, 3.0]))
    sig.update_signal(signal=np.array([3.0, 2.0, 1.0]))
    np.testing.assert_allclose(sig.signal, np.array([3.0, 2.0, 1.0]))


def test_composite_signal_concatenates_components():
    t = np.array([0.0, 1.0])
    a = FluxSignal(type=1, t_list=t, amplitude=1.0)
    b = FluxSignal(type=1, t_list=t, amplitude=2.0)
    comp = CompositeSignal([a, b])
    assert len(comp.t_list) == 4
    np.testing.assert_allclose(comp.signal, np.array([1.0, 1.0, 2.0, 2.0]))
