"""Unit tests for KernelEstimator."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.sequence import create_ramsey_pulse
from sqc.devices.transmon import TransmonQubit
from sqc.reconstruction.kernel import KernelEstimator


pytestmark = pytest.mark.unit


def _make_qubit():
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000,
        T2=8000,
        n_levels=2,
    )


def test_kernel_estimator_sets_same_values_as_composite_pulse_get_kernel():
    q = _make_qubit()
    pulse = create_ramsey_pulse(np.linspace(0, 4, 6), tau=0.0, omega_d=q.frequency)
    estimator = KernelEstimator()

    t_direct, k_direct = estimator.estimate(pulse, q)
    t_method, k_method = pulse.get_kernel(q)

    np.testing.assert_allclose(t_method, t_direct)
    np.testing.assert_allclose(k_method, k_direct)
    np.testing.assert_allclose(pulse.t_samples, t_direct)
    np.testing.assert_allclose(pulse.kernel, k_direct)


def test_kernel_auto_calibrate_returns_finite_trace():
    q = _make_qubit()
    q.change_flux(q.optimal_work_point() / np.pi)
    pulse = create_ramsey_pulse(np.linspace(0, 4, 6), tau=0.0, omega_d=q.frequency)
    t_samples, kernel = KernelEstimator(auto_calibrate=True).estimate(pulse, q)

    assert t_samples.shape == kernel.shape
    assert np.all(np.isfinite(kernel))
