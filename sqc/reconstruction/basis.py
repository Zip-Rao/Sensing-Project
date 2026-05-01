"""sqc.reconstruction.basis — basis function generation and decomposition.

Verbatim port from src/analysis.py:327-401.

Provides:
- generate_basis_functions: create a list of basis functions (bspline, fourier, legendre)
- basis_function_decomposition: decompose a signal onto basis functions via least squares
- regularization_matrix: smoothness-encouraging regularization matrix for basis coefficients
- R: backward-compatible alias for regularization_matrix
"""
from __future__ import annotations

from typing import Literal

import numpy as np

BasisType = Literal["bspline", "fourier", "legendre"]


def generate_basis_functions(
    basis_type: BasisType,
    n_basis: int,
    t_min: float,
    t_max: float,
) -> list:
    """Generate a list of basis functions on [t_min, t_max].

    Verbatim port from src/analysis.py:327-359.

    Parameters
    ----------
    basis_type : {"bspline", "fourier", "legendre"}
        Type of basis functions.
    n_basis : int
        Number of basis functions.
    t_min : float
        Start of time domain.
    t_max : float
        End of time domain.

    Returns
    -------
    list of callable
        Each callable takes a time array and returns the basis values.
    """
    if basis_type == "bspline":
        from scipy.interpolate import BSpline

        # degree=3 B-spline needs n_basis-2 interior knots
        knots = np.linspace(t_min, t_max, n_basis - 2)
        knots = np.r_[[t_min] * 3, knots, [t_max] * 3]
        basis_functions = []
        for i in range(n_basis):
            coeffs = np.zeros(n_basis)
            coeffs[i] = 1.0
            basis_functions.append(BSpline(knots, coeffs, 3))
        return basis_functions

    elif basis_type == "fourier":
        T = t_max - t_min
        basis_functions = []
        # constant term
        basis_functions.append(lambda t: np.ones_like(t) / np.sqrt(T))
        n_max = (n_basis - 1) // 2
        for n in range(1, n_max + 1):
            basis_functions.append(
                lambda t, n=n: np.sqrt(2 / T) * np.sin(
                    2 * np.pi * n * (t - t_min) / T
                )
            )
            basis_functions.append(
                lambda t, n=n: np.sqrt(2 / T) * np.cos(
                    2 * np.pi * n * (t - t_min) / T
                )
            )
        return basis_functions[:n_basis]

    elif basis_type == "legendre":
        from scipy.special import legendre

        basis_functions = [
            lambda t, n=n: legendre(n)(
                (2 * (t - t_min) / (t_max - t_min)) - 1
            )
            for n in range(n_basis)
        ]
        return basis_functions

    else:
        raise ValueError(f"Unsupported basis type: {basis_type}")


def basis_function_decomposition(
    sig,
    t_array: np.ndarray,
    basis_functions: list,
) -> np.ndarray:
    """Decompose a signal onto basis functions via least squares.

    Verbatim port from src/analysis.py:361-384.

    Parameters
    ----------
    sig : np.ndarray or object with .signal attribute
        The signal to decompose. If it has a .signal attribute (e.g.,
        src.signal.Signal or sqc.control.flux_signal.FluxSignal),
        sig.signal is used.
    t_array : np.ndarray
        Time axis for the signal.
    basis_functions : list of callable
        Basis functions over the time domain.

    Returns
    -------
    np.ndarray (n_basis,)
        Coefficients of the basis function decomposition (least-squares fit).
    """
    # Unwrap Signal-like objects
    if hasattr(sig, "signal"):
        signal = np.asarray(sig.signal, dtype=float)
    else:
        signal = np.asarray(sig, dtype=float)

    n_basis = len(basis_functions)
    n_points = len(t_array)

    # Build design matrix A (n_points x n_basis)
    A = np.zeros((n_points, n_basis))
    for i, basis_func in enumerate(basis_functions):
        A[:, i] = basis_func(t_array)

    # Solve least-squares problem A @ b = signal
    b, _, _, _ = np.linalg.lstsq(A, signal, rcond=None)

    return b


def regularization_matrix(n: int, basis_type: BasisType) -> np.ndarray:
    """Smoothness-encouraging regularization matrix for basis coefficients.

    Verbatim port from src/analysis.py:386-400 (function ``R``).

    For Fourier basis: penalizes higher frequencies by frequency order^2.
    For other bases: uses second-order finite-difference approximation
    to the second derivative.

    Parameters
    ----------
    n : int
        Number of basis functions.
    basis_type : {"bspline", "fourier", "legendre"}
        Type of basis.

    Returns
    -------
    np.ndarray (n, n)
        Regularization matrix.
    """
    if basis_type == "fourier":
        # frequency order: constant term=0, sin/cos pairs increase by 1 each
        freq_order = np.array([(i + 1) // 2 for i in range(n)])
        R_mat = np.diag(freq_order**2)
    else:
        D = np.zeros((n - 2, n))
        for i in range(n - 2):
            D[i, i] = 1
            D[i, i + 1] = -2
            D[i, i + 2] = 1
        R_mat = D.T @ D
    return R_mat


# Backward-compatible alias (matches src/analysis.py:386)
R = regularization_matrix
