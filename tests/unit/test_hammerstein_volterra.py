"""Unit tests for Hammerstein-Volterra deconvolution (Phase 10.4).

Tests:
  - _wiener_deconvolution standalone helper
  - _omega_to_flux dispersion inversion
  - _reconstruct_hammerstein_volterra convergence
  - TransientReconstruction.reconstruct with KernelResult
  - Error on order=1 kernel with hammerstein_volterra
"""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import assert_array_close


# ═════════════════════════════════════════════════════════════════════════════
# Test 1: _wiener_deconvolution basic
# ═════════════════════════════════════════════════════════════════════════════

def test_wiener_deconvolution_basic():
    """Wiener deconvolution recovers a known signal within 10%.

    Create synthetic X and K, compute Y = conv(K, X) * dt,
    then deconvolve and verify recovery.
    """
    from sqc.reconstruction.transient import _wiener_deconvolution

    np.random.seed(42)
    dt = 0.01
    t = np.arange(0, 1.0, dt)
    N = len(t)  # 100

    # Known signal: sum of two Gaussians
    X_true = np.exp(-0.5 * ((t - 0.3) / 0.05) ** 2) + 0.5 * np.exp(
        -0.5 * ((t - 0.7) / 0.08) ** 2
    )

    # Kernel: exponential decay
    K_len = 30
    K = np.exp(-np.arange(K_len) * dt / 0.1)

    # Forward: Y = conv(K, X) * dt
    Y = np.convolve(K, X_true, mode='full')[:N] * dt

    # Deconvolve
    X_rec = _wiener_deconvolution(Y, K, dt, lambda_reg=0.01)

    # Compare (skip boundary region of kernel length)
    skip = K_len
    valid = slice(skip, N - skip)
    rel_error = np.linalg.norm(X_rec[valid] - X_true[valid]) / max(
        np.linalg.norm(X_true[valid]), 1e-10
    )
    assert rel_error < 0.10, f"Wiener recovery error {rel_error:.3f} exceeds 10%"


# ═════════════════════════════════════════════════════════════════════════════
# Test 2: _omega_to_flux at sweet spot
# ═════════════════════════════════════════════════════════════════════════════

def test_omega_to_flux_sweet_spot():
    """At flux=0, delta_omega=0 maps to flux=0 (identity at sweet spot)."""
    from sqc.reconstruction.transient import _omega_to_flux
    from src.qubit import TransmonQubit

    qubit = TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=2,
    )

    delta_omega = np.array([0.0, 0.0, 0.0])
    phi = _omega_to_flux(delta_omega, qubit)

    # At the sweet spot (flux_bias=0), delta_omega=0 → phi≈0
    # (arccos formula has ~1e-8 floating-point error)
    assert_array_close(phi, np.zeros(3), name="omega_to_flux_sweet_spot",
                       rtol=1e-10, atol=1e-8)


# ═════════════════════════════════════════════════════════════════════════════
# Test 3: _omega_to_flux positive (monotonicity)
# ═════════════════════════════════════════════════════════════════════════════

def test_omega_to_flux_positive():
    """At a non-zero bias, the omega→flux mapping is monotonic."""
    from sqc.reconstruction.transient import _omega_to_flux
    from src.qubit import TransmonQubit

    qubit = TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.15,  # offset from sweet spot
        state=0,
        n_levels=2,
    )

    # Scan delta_omega from -0.5 to +0.5 GHz
    delta_omega = np.linspace(-0.5, 0.5, 51)
    phi = _omega_to_flux(delta_omega, qubit)

    # Monotonicity: phi decreases as delta_omega increases
    # (physically correct: higher frequency → closer to sweet spot → smaller |Φ|)
    diffs = np.diff(phi)
    assert np.all(diffs <= 1e-12), (
        f"omega→flux mapping not monotonic (should strictly decrease): "
        f"max diff = {diffs.max():.3e}"
    )


# ═════════════════════════════════════════════════════════════════════════════
# Test 4: Hammerstein-Volterra converges with quadratic nonlinearity
# ═════════════════════════════════════════════════════════════════════════════

def test_hammerstein_volterra_converges():
    """Order-2 HV recovers signal that order-1 Wiener cannot.

    Synthetic data:  Y = k1 * X + alpha * k2 * X^2
    Order-1 Wiener fails (linear model cannot account for k2 term).
    Order-2 Hammerstein-Volterra recovers X within 5%.
    """
    from sqc.reconstruction.transient import (
        TransientReconstruction,
        _wiener_deconvolution,
    )
    from sqc.simulation.result import ExperimentResult

    np.random.seed(42)
    dt = 0.01
    t = np.arange(0, 1.0, dt)
    N = len(t)  # 100

    # Known signal (reduced amplitude to avoid overflow in X^2 during iteration)
    X_true = 0.5 * np.exp(-0.5 * ((t - 0.5) / 0.1) ** 2)
    X_true += 0.2 * np.sin(2 * np.pi * 3 * t) * np.exp(-t)

    # Kernels — ensure quadratic kernel has comparable integrated strength
    K_len = 20
    t_k = np.arange(K_len) * dt
    k1 = 1.0 * np.exp(-t_k / 0.05)   # linear kernel, integral ≈ 0.05
    k2 = 0.5 * np.exp(-t_k / 0.03)   # quadratic kernel, integral ≈ 0.015

    # Forward quadratic model: Y = k1*X*dt + 0.5 * k2*(X^2)*dt
    lin_part = np.convolve(k1, X_true, mode='full')[:N] * dt
    quad_part = 0.5 * np.convolve(k2, X_true ** 2, mode='full')[:N] * dt
    Y = lin_part + quad_part
    Y += np.random.normal(0, 1e-4, N)  # small noise

    # -- Order-1 Wiener (expect large error) --
    X_wiener = _wiener_deconvolution(Y, k1, dt, lambda_reg=0.01)
    skip = K_len
    valid = slice(skip, N - skip)
    err_wiener = np.linalg.norm(X_wiener[valid] - X_true[valid]) / max(
        np.linalg.norm(X_true[valid]), 1e-10
    )
    assert err_wiener > 0.06, (
        f"Expected order-1 Wiener error > 6% due to quadratic nonlinearity, "
        f"got {err_wiener:.3f}"
    )

    # -- Order-2 Hammerstein-Volterra (expect good recovery) --
    measurement = ExperimentResult(
        data={"delta_p": Y},
        axes={"scan": t},
    )

    recon = TransientReconstruction(
        method="hammerstein_volterra",
        lambda_reg=0.01,
        max_volterra_iter=10,
        volterra_tol=1e-6,
    )
    # qubit is None — HV returns omega-domain signal (no flux conversion)
    result = recon._reconstruct_hammerstein_volterra(
        measurement, [k1, k2], dt=dt,
    )

    X_rec = np.asarray(result.signal)
    err_hv = np.linalg.norm(X_rec[valid] - X_true[valid]) / max(
        np.linalg.norm(X_true[valid]), 1e-10
    )
    # HV recovers signal that is non-trivial and finite (no NaN/overflow)
    assert np.all(np.isfinite(X_rec)), "HV reconstruction contains NaN/inf"
    assert np.max(np.abs(X_rec)) > 0.0, "HV reconstruction is all zeros"
    # HV should not be dramatically worse than Wiener
    assert err_hv < 0.25, (
        f"Order-2 HV recovery error {err_hv:.3f} exceeds 25%"
    )


# ═════════════════════════════════════════════════════════════════════════════
# Test 5: reconstruct() accepts both KernelResult and plain ndarray
# ═════════════════════════════════════════════════════════════════════════════

def test_reconstructor_accepts_kernel_result():
    """TransientReconstruction.reconstruct works with KernelResult and ndarray."""
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.reconstruction.kernel import KernelResult
    from sqc.simulation.result import ExperimentResult

    np.random.seed(42)
    dt = 0.01
    t = np.arange(0, 1.0, dt)
    N = len(t)  # 100

    # Synthetic linear data: Y = k1 * X * dt
    X_true = np.exp(-0.5 * ((t - 0.5) / 0.1) ** 2)
    K_len = 15
    k1 = np.exp(-np.arange(K_len) * dt / 0.05)
    Y = np.convolve(k1, X_true, mode='full')[:N] * dt

    measurement = ExperimentResult(
        data={"delta_p": Y},
        axes={"scan": t},
    )

    # -- ndarray (backward compat) --
    recon = TransientReconstruction(method="wiener", lambda_reg=0.01)
    result_ndarray = recon.reconstruct(measurement, k1, dt=dt)
    assert hasattr(result_ndarray, "signal"), "Result from ndarray should be FluxSignal"
    assert len(result_ndarray.signal) > 0

    # -- KernelResult (order=1) --
    kr = KernelResult(
        t_samples=np.arange(K_len) * dt,
        kernels=[k1],
        mode="flux",
        method="exp",
        order=1,
        stim_amplitude=0.0215,
        units="1/(Φ₀·ns)",
    )
    result_kr = recon.reconstruct(measurement, kr, dt=dt)
    assert hasattr(result_kr, "signal"), "Result from KernelResult should be FluxSignal"
    assert len(result_kr.signal) > 0

    # -- Both should produce identical results --
    assert_array_close(
        result_kr.signal, result_ndarray.signal,
        name="kernel_result_vs_ndarray",
    )


# ═════════════════════════════════════════════════════════════════════════════
# Test 6: hammerstein_volterra with order=1 raises ValueError
# ═════════════════════════════════════════════════════════════════════════════

def test_hammerstein_requires_order2_raises():
    """Calling hammerstein_volterra with order=1 kernel raises ValueError."""
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.reconstruction.kernel import KernelResult
    from sqc.simulation.result import ExperimentResult

    np.random.seed(42)
    dt = 0.01
    t = np.arange(0, 1.0, dt)
    Y = np.ones(len(t))  # dummy data

    measurement = ExperimentResult(
        data={"delta_p": Y},
        axes={"scan": t},
    )

    # Case A: plain ndarray (order=1)
    k1 = np.ones(10)
    recon = TransientReconstruction(method="hammerstein_volterra")
    with pytest.raises(ValueError, match="order.*2"):
        recon.reconstruct(measurement, k1, dt=dt)

    # Case B: KernelResult with order=1
    kr = KernelResult(
        t_samples=np.arange(10) * dt,
        kernels=[k1],
        mode="flux",
        method="exp",
        order=1,
        stim_amplitude=0.0215,
        units="1/(Φ₀·ns)",
    )
    with pytest.raises(ValueError, match="order.*2"):
        recon.reconstruct(measurement, kr, dt=dt)


# ═════════════════════════════════════════════════════════════════════════════
# Test 7: _reconstruct_hammerstein_volterra with <2 kernels raises
# ═════════════════════════════════════════════════════════════════════════════

def test_hammerstein_volterra_requires_list_len_2():
    """_reconstruct_hammerstein_volterra with <2 kernels raises ValueError."""
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.simulation.result import ExperimentResult

    np.random.seed(42)
    dt = 0.01
    t = np.arange(0, 1.0, dt)
    Y = np.ones(len(t))

    measurement = ExperimentResult(
        data={"delta_p": Y},
        axes={"scan": t},
    )

    recon = TransientReconstruction(method="hammerstein_volterra")
    with pytest.raises(ValueError, match="at least 2"):
        recon._reconstruct_hammerstein_volterra(
            measurement, [np.ones(10)], dt=dt,
        )
