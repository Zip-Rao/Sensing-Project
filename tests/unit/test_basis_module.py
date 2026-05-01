"""Unit tests for reconstruction basis helpers."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.reconstruction.basis import (
    R,
    basis_function_decomposition,
    generate_basis_functions,
    regularization_matrix,
)


pytestmark = pytest.mark.unit


def test_fourier_basis_count_and_shape():
    basis = generate_basis_functions("fourier", 5, 0.0, 10.0)
    t = np.linspace(0, 10, 11)

    assert len(basis) == 5
    assert basis[0](t).shape == t.shape


def test_basis_decomposition_recovers_constant_component():
    t = np.linspace(0, 10, 50)
    basis = generate_basis_functions("fourier", 3, 0.0, 10.0)
    signal = np.ones_like(t)
    coeffs = basis_function_decomposition(signal, t, basis)
    reconstructed = sum(c * f(t) for c, f in zip(coeffs, basis))

    np.testing.assert_allclose(reconstructed, signal, atol=1e-10)


def test_regularization_alias():
    np.testing.assert_allclose(R(5, "fourier"), regularization_matrix(5, "fourier"))
    assert regularization_matrix(5, "bspline").shape == (5, 5)
