"""Unit tests for KernelEstimator."""
from __future__ import annotations

import numpy as np
import pytest

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


def test_kernel_estimator_default_params():
    """KernelEstimator with default parameters can be created."""
    from sqc.reconstruction.kernel import KernelEstimator

    ke = KernelEstimator()
    assert ke.stim_amplitude == 0.0215
    assert ke.stim_width == 3.0
    assert not ke.auto_calibrate


def test_kernel_estimator_estimate_ramsey_pulse(qubit):
    """KernelEstimator.estimate returns correct shape for a ramsey pulse."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse

    t_rabi = np.linspace(0, 10, 10)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    ke = KernelEstimator()
    t_samples, kernel = ke.estimate(pulse, qubit)

    assert isinstance(t_samples, np.ndarray)
    assert isinstance(kernel, np.ndarray)
    assert len(t_samples) == len(kernel)
    assert len(t_samples) > 0
    # Kernel values should be finite
    assert np.all(np.isfinite(kernel))


def test_kernel_estimator_estimate_nonzero_kernel(qubit):
    """KernelEstimator produces non-trivial (non-zero) kernel values."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse

    t_rabi = np.linspace(0, 10, 20)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    ke = KernelEstimator()
    t_samples, kernel = ke.estimate(pulse, qubit)

    # Kernel should not be all zeros
    assert np.max(np.abs(kernel)) > 0.0
    # There should be variation across time points
    assert np.max(kernel) != np.min(kernel) or len(kernel) == 1


def test_kernel_estimator_matches_legacy(qubit):
    """KernelEstimator output matches legacy CompositePulse.get_kernel."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    from sqc.control.pulse import CompositePulse

    t_rabi = np.linspace(0, 10, 20)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    # Legacy get_kernel
    pulse.get_kernel(qubit)
    legacy_ts = np.array(pulse.t_samples)
    legacy_k = np.array(pulse.kernel)

    # New KernelEstimator
    ke = KernelEstimator()
    new_ts, new_k = ke.estimate(pulse, qubit)

    assert_array_close(new_k, legacy_k, name="kernel_vs_legacy")
    assert_array_close(new_ts, legacy_ts, name="t_samples_vs_legacy")
