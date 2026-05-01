"""Unit tests for Waveform and CompositeWaveform.

Tests: construction, shape validation, value_at, truncate (returns new),
CompositeWaveform concatenation.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import Waveform, CompositeWaveform


class TestWaveform:
    """Waveform dataclass correctness."""

    def test_construction(self):
        """Basic construction works."""
        t = np.linspace(0, 100, 101)
        s = np.sin(t)
        w = Waveform(t_list=t, samples=s)
        assert w.duration == pytest.approx(100.0, rel=1e-12)
        assert w.n_points == 101

    def test_shape_mismatch_raises(self):
        """t_list and samples must have same shape."""
        with pytest.raises(ValueError, match="shape mismatch"):
            Waveform(t_list=np.linspace(0, 10, 11), samples=np.zeros(10))

    def test_value_at(self):
        """value_at should return nearest-index value."""
        t = np.array([0, 10, 20, 30], dtype=float)
        s = np.array([0, 1, 2, 3], dtype=float)
        w = Waveform(t_list=t, samples=s)
        assert w.value_at(0) == 0.0
        assert w.value_at(10) == 1.0
        assert w.value_at(20) == 2.0

    def test_value_at_out_of_range(self):
        """value_at should return 0 for times outside range."""
        t = np.linspace(0, 100, 101)
        w = Waveform(t_list=t, samples=np.ones_like(t))
        assert w.value_at(-1) == 0.0
        assert w.value_at(200) == 0.0

    def test_truncate_returns_new_object(self):
        """truncate must return a new Waveform, not modify the original."""
        t = np.linspace(0, 100, 101)
        s = np.ones_like(t)
        w = Waveform(t_list=t.copy(), samples=s.copy())
        w2 = w.truncate(20, 80)
        assert w2 is not w
        # Original unchanged
        assert np.all(w.samples == 1.0)
        # New has zeros outside [20, 80]
        assert w2.samples[0] == 0.0
        assert w2.samples[-1] == 0.0
        assert np.all(w2.samples[20:81] == 1.0)

    def test_copy(self):
        """copy should produce an independent copy."""
        t = np.linspace(0, 50, 51)
        w = Waveform(t_list=t, samples=np.arange(51, dtype=float))
        w2 = w.copy()
        assert w2 is not w
        assert np.array_equal(w2.samples, w.samples)
        w2.samples[0] = 999
        assert w.samples[0] != 999  # original unaffected


class TestCompositeWaveform:
    """CompositeWaveform concatenation."""

    def test_from_components(self):
        """Concatenation merges time axes and samples."""
        w1 = Waveform(t_list=np.array([0, 1, 2], dtype=float), samples=np.array([1, 1, 1], dtype=float))
        w2 = Waveform(t_list=np.array([0, 1], dtype=float), samples=np.array([2, 2], dtype=float))
        cw = CompositeWaveform.from_components([w1, w2])
        assert len(cw.t_list) == 5
        # Check values: w1's samples first, then w2's
        assert cw.samples[0] == 1.0
        assert cw.samples[-1] == 2.0
        # Time axis should start near 0
        assert cw.t_list[0] == pytest.approx(0.0)

    def test_empty_components_handled(self):
        """Empty component list should produce an empty waveform."""
        cw = CompositeWaveform.from_components([])
        assert len(cw.t_list) == 0
        assert len(cw.samples) == 0
