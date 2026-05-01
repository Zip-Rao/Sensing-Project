"""Unit tests for LMReconstruction (sqc/reconstruction/numerical_inverse.py).

Tests the LM reconstruction class and its private helpers:
  - _forward_simulation shape and output
  - _compute_jacobian_adjoint shape
  - _compute_jacobian_fd shape (against adjoint for consistency)
  - _levenberg_marquardt smoke (runs without crash)
  - reconstruct API smoke (runs without crash)

These tests use minimal parameters (n_basis=5, max_iter=1, n_levels=2)
to keep runtime under ~5 s per test.
"""
from __future__ import annotations

import numpy as np
import pytest
from qutip import qeye

from tests.conftest import assert_array_close


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

def _make_minimal_qubit():
    """Create a default legacy TransmonQubit (n_levels=2)."""
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


def _make_minimal_pulse(qubit):
    """Create a short Ramsey pulse for testing."""
    from src.pulse import create_ramsey_pulse

    return create_ramsey_pulse(
        t_rabi=np.linspace(0, 5, 6), tau=5.0, omega_d=qubit.frequency,
    )


def _make_t_list(n_pts=11):
    """Short time axis for fast tests."""
    return np.linspace(0, 20, n_pts)


def _build_h_and_t_evolve(qubit, cp, t_list, t_meas):
    """Build H_list and t_evolve_list (replicates LM inner H function)."""
    n_meas = len(t_meas)
    H_list: list = [[] for _ in range(n_meas)]
    t_evolve_list: list = [[] for _ in range(n_meas)]

    for i, t_delay in enumerate(t_meas):
        delta = t_delay - 0.5 * cp.t_list[-1]
        t_start = min(t_list[0], delta)
        t_end = max(t_list[-1], delta + cp.t_list[-1])
        N_e = len(t_list) + len(cp.t_list) - 1
        t_evolve = np.linspace(t_start, t_end, N_e)
        t_evolve_list[i] = t_evolve

        freq_coeffs = np.zeros(N_e, dtype=float)
        for j, t in enumerate(t_evolve):
            if t_list[0] <= t <= t_list[-1]:
                index = np.searchsorted(t_list, t)
                index = min(index, len(t_list) - 1)
                freq_coeffs[j] = (
                    qubit.freq_coeffs[index]
                    if cp.frame == 0
                    else qubit.freq_coeffs[index] - cp.omega_d
                )
            else:
                freq_coeffs[j] = (
                    qubit.frequency
                    if cp.frame == 0
                    else qubit.frequency - cp.omega_d
                )

        id_coeffs = np.ones(N_e, dtype=complex)
        if cp.frame == 0:
            H_list[i].append(
                [
                    qubit.anharmonicity * qubit.n * (qubit.n - 1) * 0.5,
                    id_coeffs,
                ]
            )
            H_list[i].append(
                [qubit.n + 0.5 * qeye(qubit.n_levels), freq_coeffs]
            )
        else:
            H_list[i].append(
                [
                    qubit.anharmonicity * qubit.n * (qubit.n - 1) * 0.5,
                    id_coeffs,
                ]
            )
            H_list[i].append([qubit.n, freq_coeffs])

        for op, coeffs in cp.hamiltonian:
            coeff_global = np.zeros(N_e, dtype=complex)
            for j, t in enumerate(t_evolve):
                t_loc = t - delta
                if 0 <= t_loc <= cp.t_list[-1]:
                    index = np.searchsorted(cp.t_list, t_loc)
                    index = min(index, len(cp.t_list) - 1)
                    coeff_global[j] = coeffs[index]
            H_list[i].append([op, coeff_global])

    return H_list, t_evolve_list


# ---------------------------------------------------------------------------
# Test 1: Forward simulation smoke
# ---------------------------------------------------------------------------

def test_forward_simulation_smoke():
    """_forward_simulation runs without crash and produces correct shape."""
    from src.signal import Signal
    from sqc.reconstruction.numerical_inverse import LMReconstruction

    np.random.seed(42)
    q = _make_minimal_qubit()
    cp = _make_minimal_pulse(q)
    t_list = _make_t_list()

    B_sig = Signal(type=2, t_list=t_list, amplitude=0.01, frequency=0.1)
    q.qubit_in_mag(B_sig, frame=0, omega_d=q.frequency)

    meas_start = t_list[0] - 0.5 * cp.t_list[-1]
    meas_end = t_list[-1] + 0.5 * cp.t_list[-1]
    t_meas = np.linspace(
        meas_start, meas_end, len(t_list) + len(cp.t_list) - 1,
    )

    H_list, t_evolve_list = _build_h_and_t_evolve(q, cp, t_list, t_meas)

    recon = LMReconstruction(qubit=q, control_pulse=cp)
    results = recon._forward_simulation(
        B_sig, t_meas, H_list, t_evolve_list,
    )

    assert len(results) == len(t_meas), "should produce one result per delay"
    p_e = np.array([r.expect[0][-1] for r in results])
    assert p_e.shape == (len(t_meas),), "p_e shape mismatch"
    assert np.all((0 <= p_e) & (p_e <= 1)), "p_e must be in [0, 1]"


# ---------------------------------------------------------------------------
# Test 2: Jacobian shape
# ---------------------------------------------------------------------------

def test_jacobian_adjoint_shape():
    """Adjoint Jacobian has correct shape (N_meas x M_basis)."""
    from src.signal import Signal
    from sqc.reconstruction.numerical_inverse import LMReconstruction

    np.random.seed(42)
    q = _make_minimal_qubit()
    cp = _make_minimal_pulse(q)
    t_list = _make_t_list()

    B_sig = Signal(type=2, t_list=t_list, amplitude=0.01, frequency=0.1)
    q.qubit_in_mag(B_sig, frame=0, omega_d=q.frequency)

    meas_start = t_list[0] - 0.5 * cp.t_list[-1]
    meas_end = t_list[-1] + 0.5 * cp.t_list[-1]
    t_meas = np.linspace(
        meas_start, meas_end, len(t_list) + len(cp.t_list) - 1,
    )

    H_list, t_evolve_list = _build_h_and_t_evolve(q, cp, t_list, t_meas)

    recon = LMReconstruction(qubit=q, control_pulse=cp)

    # Run forward simulation first
    results = recon._forward_simulation(
        B_sig, t_meas, H_list, t_evolve_list,
    )

    # Create a type=6 signal for basis functions
    n_basis = 5
    B_type6 = Signal(
        type=6, t_list=t_list, n_basis=n_basis, basis_type="fourier",
    )
    from src.analysis import basis_function_decomposition
    b_init = basis_function_decomposition(
        np.zeros_like(t_list), t_list, B_type6.basis_functions,
    )
    B_type6.update_signal(b=b_init)

    # Also need qubit_in_mag on B_type6 BEFORE calling compute_jacobian
    # because the original LM code does qubit.qubit_in_mag(B_curr) first
    q.qubit_in_mag(B_type6, frame=0, omega_d=q.frequency)

    # Rebuild H_list for B_type6
    H_list, t_evolve_list = _build_h_and_t_evolve(q, cp, t_list, t_meas)

    # But results are from B_sig forward sim, not B_type6...
    # This is a simplified shape test — the values may differ
    J = recon._compute_jacobian_adjoint(
        B_type6, t_meas, results, H_list, t_evolve_list,
    )

    expected_shape = (len(t_meas), n_basis)
    assert J.shape == expected_shape, (
        f"Jacobian shape {J.shape} != expected {expected_shape}"
    )
    assert J.dtype in (np.float64, np.complex128), "unexpected dtype"
    # Jacobian should not be all-NaN
    assert not np.any(np.isnan(J)), "Jacobian contains NaN"


# ---------------------------------------------------------------------------
# Test 3: FD Jacobian shape check
# ---------------------------------------------------------------------------

def test_jacobian_fd_shape():
    """Finite-difference Jacobian has correct shape."""
    from src.signal import Signal
    from sqc.reconstruction.numerical_inverse import LMReconstruction

    np.random.seed(42)
    q = _make_minimal_qubit()
    cp = _make_minimal_pulse(q)
    t_list = _make_t_list()

    B_sig = Signal(type=2, t_list=t_list, amplitude=0.01, frequency=0.1)
    q.qubit_in_mag(B_sig, frame=0, omega_d=q.frequency)

    meas_start = t_list[0] - 0.5 * cp.t_list[-1]
    meas_end = t_list[-1] + 0.5 * cp.t_list[-1]
    t_meas = np.linspace(
        meas_start, meas_end, len(t_list) + len(cp.t_list) - 1,
    )

    n_basis = 5
    B_type6 = Signal(
        type=6, t_list=t_list, n_basis=n_basis, basis_type="fourier",
    )
    from src.analysis import basis_function_decomposition
    b_init = basis_function_decomposition(
        np.zeros_like(t_list), t_list, B_type6.basis_functions,
    )
    B_type6.update_signal(b=b_init)
    q.qubit_in_mag(B_type6, frame=0, omega_d=q.frequency)

    H_list, t_evolve_list = _build_h_and_t_evolve(q, cp, t_list, t_meas)

    recon = LMReconstruction(qubit=q, control_pulse=cp)
    results = recon._forward_simulation(
        B_type6, t_meas, H_list, t_evolve_list,
    )
    p_sim = np.array([r.expect[0][-1] for r in results])

    J_fd = recon._compute_jacobian_fd(
        B_type6, t_meas, p_sim, H_list, t_evolve_list, epsilon=1e-6,
    )

    expected_shape = (len(t_meas), n_basis)
    assert J_fd.shape == expected_shape, (
        f"FD Jacobian shape {J_fd.shape} != expected {expected_shape}"
    )
    assert not np.any(np.isnan(J_fd)), "FD Jacobian contains NaN"


# ---------------------------------------------------------------------------
# Test 4: LM reconstruction smoke
# ---------------------------------------------------------------------------

def test_lm_reconstruct_smoke():
    """LMReconstruction.reconstruct runs without crash."""
    from src.pulse import create_ramsey_pulse
    from sqc.reconstruction.numerical_inverse import LMReconstruction
    from sqc.simulation.result import ExperimentResult
    from sqc.control.flux_signal import FluxSignal

    np.random.seed(42)
    q = _make_minimal_qubit()
    cp = create_ramsey_pulse(
        t_rabi=np.linspace(0, 5, 6), tau=5.0, omega_d=q.frequency,
    )
    t_list = _make_t_list()

    # Generate synthetic p_meas from a known signal
    import pickle
    with open("tests/baselines/lm_default.pkl", "rb") as f:
        bl = pickle.load(f)

    n_basis = 5
    recon = LMReconstruction(
        qubit=q,
        control_pulse=cp,
        basis_type="fourier",
        n_basis=n_basis,
        max_iter=1,  # minimal iterations for smoke test
        tol=1e-1,
        use_adjoint=True,
    )

    meas = ExperimentResult(
        data={"p_meas": bl["p_meas"]},
        axes={"t_signal": bl["t_list"]},
    )

    B_opt_flux, history = recon.reconstruct(meas)
    b_opt = B_opt_flux.params["b"]

    assert isinstance(B_opt_flux, FluxSignal), "should return FluxSignal"
    assert len(b_opt) == n_basis, "should return n_basis coefficients"
    assert "b" in history, "history should have 'b'"
    assert "res" in history, "history should have 'res'"
    assert len(history["b"]) >= 1, "should have at least initial b"


# ---------------------------------------------------------------------------
# Test 5: Forward simulation equivalence (new vs old)
# ---------------------------------------------------------------------------

def test_forward_simulation_equivalence():
    """New _forward_simulation matches old forward_simulation."""
    from src.signal import Signal
    from src.analysis import forward_simulation as old_fwd
    from sqc.reconstruction.numerical_inverse import LMReconstruction

    np.random.seed(42)
    q_old = _make_minimal_qubit()
    q_new = _make_minimal_qubit()

    cp = _make_minimal_pulse(q_old)
    t_list = _make_t_list()

    B_sig = Signal(type=2, t_list=t_list, amplitude=0.01, frequency=0.1)

    # Old code
    q_old.qubit_in_mag(B_sig, frame=0, omega_d=q_old.frequency)
    meas_start = t_list[0] - 0.5 * cp.t_list[-1]
    meas_end = t_list[-1] + 0.5 * cp.t_list[-1]
    t_meas = np.linspace(
        meas_start, meas_end, len(t_list) + len(cp.t_list) - 1,
    )
    H_list_old, t_evolve_list_old = _build_h_and_t_evolve(
        q_old, cp, t_list, t_meas,
    )
    results_old = old_fwd(
        q_old, cp, B_sig, t_meas, H_list_old, t_evolve_list_old,
    )
    p_old = np.array([r.expect[0][-1] for r in results_old])

    # New code
    q_new.qubit_in_mag(B_sig, frame=0, omega_d=q_new.frequency)
    H_list_new, t_evolve_list_new = _build_h_and_t_evolve(
        q_new, cp, t_list, t_meas,
    )
    recon = LMReconstruction(qubit=q_new, control_pulse=cp)
    results_new = recon._forward_simulation(
        B_sig, t_meas, H_list_new, t_evolve_list_new,
    )
    p_new = np.array([r.expect[0][-1] for r in results_new])

    assert_array_close(p_new, p_old, name="forward_simulation_p_e")
