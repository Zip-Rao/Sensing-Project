"""Unit tests for FluxSignal — type-based construction.

Tests: types 0-8 equivalence with old src.signal.Signal, value_at,
truncate (in-place), update_signal, copy.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.flux_signal import FluxSignal, CompositeSignal


class TestFluxSignalTypes:
    """FluxSignal type-based construction matches old Signal."""

    @pytest.fixture
    def t_list(self) -> np.ndarray:
        return np.linspace(0, 200, 400)

    def test_type_0_zero(self, t_list):
        s = FluxSignal(type=0, t_list=t_list)
        expected = np.zeros(len(t_list))  # zero + offset(0)
        assert np.allclose(s.signal, expected, atol=1e-10)

    def test_type_1_constant(self, t_list):
        s = FluxSignal(type=1, t_list=t_list, amplitude=5.0, offset=1.0)
        expected = 6.0 * np.ones(len(t_list))  # 5 + 1
        assert np.allclose(s.signal, expected, atol=1e-10)

    def test_type_2_sinusoidal(self, t_list):
        s = FluxSignal(type=2, t_list=t_list, amplitude=2.0, frequency=0.01, phase=0.0, offset=0.0)
        omega = 2 * np.pi * 0.01
        expected = 2.0 * np.sin(omega * t_list)
        # Small noise from back_signal (noise_level=0); check close
        assert np.allclose(s.signal, expected, atol=1e-12)

    def test_type_3_gaussian(self, t_list):
        s = FluxSignal(type=3, t_list=t_list, amplitude=3.0, center=100.0, width=40.0, offset=0.0)
        sigma = 40.0 / 4.0
        expected = 3.0 * np.exp(-((t_list - 100.0) ** 2) / (2 * sigma**2))
        assert np.allclose(s.signal, expected, atol=1e-12)

    def test_type_8_custom(self, t_list):
        """Type 8: user-defined signal."""
        custom = np.arange(len(t_list), dtype=float)
        s = FluxSignal(type=8, t_list=t_list, signal=custom)
        assert np.allclose(s.signal, custom, atol=1e-12)

    def test_compare_with_old_signal(self, t_list):
        """FluxSignal(type=2) should match src.signal.Signal(type=2)."""
        from src.signal import Signal as OldSignal

        old = OldSignal(type=2, t_list=t_list, amplitude=1.5, frequency=0.005, phase=0.3)
        new = FluxSignal(type=2, t_list=t_list, amplitude=1.5, frequency=0.005, phase=0.3)
        assert np.allclose(old.signal, new.signal, atol=1e-10)


class TestFluxSignalMethods:
    """Legacy method behavior."""

    @pytest.fixture
    def sig(self) -> FluxSignal:
        t = np.linspace(0, 100, 101)
        return FluxSignal(type=3, t_list=t, amplitude=5.0, center=50.0, width=30.0)

    def test_value_at(self, sig):
        """value_at within range returns sample; outside returns 0."""
        assert sig.value_at(50.0) > 0.0
        assert sig.value_at(-1.0) == 0.0
        assert sig.value_at(200.0) == 0.0

    def test_truncate_in_place(self, sig):
        """truncate modifies samples in-place (legacy behavior)."""
        original = sig.signal.copy()
        sig.truncate(30, 70)
        # Signal should be zeroed outside [30, 70]
        assert sig.signal[0] == 0.0
        assert sig.signal[-1] == 0.0
        # Original values inside range should be preserved
        center_idx = 50
        assert sig.signal[center_idx] == original[center_idx]

    def test_update_signal(self, sig):
        """update_signal regenerates samples with new params."""
        old_max = np.max(sig.signal)
        sig.update_signal(amplitude=10.0)
        new_max = np.max(sig.signal)
        assert new_max == pytest.approx(2 * old_max, rel=0.01)

    def test_copy(self, sig):
        """copy produces independent object."""
        sig2 = sig.copy()
        assert sig2 is not sig
        assert np.array_equal(sig2.signal, sig.signal)
        sig2.signal[0] = 999
        assert sig.signal[0] != 999


class TestCompositeSignal:
    """CompositeSignal concatenation."""

    def test_concatenation(self):
        t = np.linspace(0, 50, 51)
        s1 = FluxSignal(type=1, t_list=t, amplitude=1.0)
        s2 = FluxSignal(type=1, t_list=t, amplitude=2.0)
        cs = CompositeSignal([s1, s2])
        # First half should be ~1, second half ~2
        mid = len(cs.t_list) // 2
        assert cs.signal[0] == pytest.approx(1.0, rel=0.01)
        assert cs.signal[-1] == pytest.approx(2.0, rel=0.01)
