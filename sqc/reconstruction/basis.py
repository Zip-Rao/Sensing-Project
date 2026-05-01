"""Basis-function utilities for waveform reconstruction."""
from __future__ import annotations

from typing import Literal

import numpy as np

from sqc.control.flux_signal import FluxSignal


BasisType = Literal["bspline", "fourier", "legendre"]


def generate_basis_functions(
    basis_type: BasisType, n_basis: int, t_min: float, t_max: float
) -> list:
    """Generate basis functions on ``[t_min, t_max]``."""
    if basis_type == "bspline":
        from scipy.interpolate import BSpline

        knots = np.linspace(t_min, t_max, n_basis - 2)
        knots = np.r_[[t_min] * 3, knots, [t_max] * 3]
        basis_functions = []
        for i in range(n_basis):
            coeffs = np.zeros(n_basis)
            coeffs[i] = 1.0
            basis_functions.append(BSpline(knots, coeffs, 3))
        return basis_functions

    if basis_type == "fourier":
        total_time = t_max - t_min
        basis_functions = [lambda t: np.ones_like(t) / np.sqrt(total_time)]
        n_max = (n_basis - 1) // 2
        for n in range(1, n_max + 1):
            basis_functions.append(
                lambda t, n=n: np.sqrt(2 / total_time)
                * np.sin(2 * np.pi * n * (t - t_min) / total_time)
            )
            basis_functions.append(
                lambda t, n=n: np.sqrt(2 / total_time)
                * np.cos(2 * np.pi * n * (t - t_min) / total_time)
            )
        return basis_functions[:n_basis]

    if basis_type == "legendre":
        from scipy.special import legendre

        return [
            lambda t, n=n: legendre(n)((2 * (t - t_min) / (t_max - t_min)) - 1)
            for n in range(n_basis)
        ]

    raise ValueError(f"Unsupported basis type: {basis_type}")


def basis_function_decomposition(
    signal_or_obj, t_array: np.ndarray, basis_functions: list
) -> np.ndarray:
    """Decompose a sampled signal onto basis functions via least squares."""
    signal = signal_or_obj.signal if isinstance(signal_or_obj, FluxSignal) else signal_or_obj
    t_array = np.asarray(t_array, dtype=float)
    signal = np.asarray(signal, dtype=float)
    design = np.zeros((len(t_array), len(basis_functions)))
    for i, basis_func in enumerate(basis_functions):
        design[:, i] = basis_func(t_array)
    coeffs, _, _, _ = np.linalg.lstsq(design, signal, rcond=None)
    return coeffs


def regularization_matrix(n: int, basis_type: BasisType) -> np.ndarray:
    """Return a smoothness regularization matrix."""
    if basis_type == "fourier":
        freq_order = np.array([(i + 1) // 2 for i in range(n)])
        return np.diag(freq_order**2)

    if n < 3:
        return np.eye(n)
    diff = np.zeros((n - 2, n))
    for i in range(n - 2):
        diff[i, i] = 1
        diff[i, i + 1] = -2
        diff[i, i + 2] = 1
    return diff.T @ diff


R = regularization_matrix
