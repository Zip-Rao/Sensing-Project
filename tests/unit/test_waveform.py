"""Unit tests for generic waveform structures."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import CompositeWaveform, Waveform


pytestmark = pytest.mark.unit


def test_waveform_rejects_shape_mismatch():
    with pytest.raises(ValueError):
        Waveform(t_list=np.array([0.0, 1.0]), samples=np.array([1.0]))


def test_waveform_value_at_and_out_of_range():
    w = Waveform(t_list=np.array([0.0, 1.0, 2.0]), samples=np.array([0.0, 2.0, 4.0]))
    assert w.value_at(1.1) == pytest.approx(2.0)
    assert w.value_at(-1.0) == pytest.approx(0.0)
    assert w.value_at(3.0) == pytest.approx(0.0)


def test_truncate_returns_new_object():
    w = Waveform(t_list=np.array([0.0, 1.0, 2.0]), samples=np.array([1.0, 2.0, 3.0]))
    truncated = w.truncate(0.5, 1.5)
    assert truncated is not w
    np.testing.assert_allclose(w.samples, np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(truncated.samples, np.array([0.0, 2.0, 0.0]))


def test_composite_waveform_from_components():
    a = Waveform(np.array([0.0, 1.0]), np.array([1.0, 2.0]))
    b = Waveform(np.array([0.0, 1.0]), np.array([3.0, 4.0]))
    c = CompositeWaveform.from_components([a, b])
    assert c.n_points == 4
    np.testing.assert_allclose(c.samples, np.array([1.0, 2.0, 3.0, 4.0]))
