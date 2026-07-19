"""Unit tests for sqc.reconstruction.basis.

Tests generate_basis_functions, basis_function_decomposition,
regularization_matrix, and backward-compat alias R.
"""
from __future__ import annotations

import numpy as np
import pytest


class TestGenerateBasisFunctions:
    """Tests for generate_basis_functions."""

    def test_fourier_basis_count(self):
        """Fourier basis produces exactly n_basis functions."""
        from sqc.reconstruction.basis import generate_basis_functions

        n = 7
        fns = generate_basis_functions("fourier", n, 0.0, 100.0)
        assert len(fns) == n

    def test_fourier_constant_term_normalized(self):
        """First Fourier basis function is 1/sqrt(T)."""
        from sqc.reconstruction.basis import generate_basis_functions

        T = 100.0
        fns = generate_basis_functions("fourier", 5, 0.0, T)
        t = np.linspace(0, T, 50)
        vals = fns[0](t)
        # constant term: 1/sqrt(T)
        np.testing.assert_allclose(vals, 1.0 / np.sqrt(T))

    def test_fourier_basis_orthogonality(self):
        """Fourier basis functions are approximately orthonormal."""
        from sqc.reconstruction.basis import generate_basis_functions

        T = 100.0
        n = 7
        fns = generate_basis_functions("fourier", n, 0.0, T)
        t = np.linspace(0, T, 2000)
        dt = t[1] - t[0]

        for i in range(n):
            for j in range(n):
                inner = np.sum(fns[i](t) * fns[j](t)) * dt
                if i == j:
                    np.testing.assert_allclose(inner, 1.0, atol=0.02)
                else:
                    np.testing.assert_allclose(inner, 0.0, atol=0.02)

    def test_bspline_basis_count(self):
        """B-spline basis produces exactly n_basis functions."""
        from sqc.reconstruction.basis import generate_basis_functions

        n = 5
        fns = generate_basis_functions("bspline", n, 0.0, 100.0)
        assert len(fns) == n

    def test_legendre_basis_count(self):
        """Legendre basis produces exactly n_basis functions."""
        from sqc.reconstruction.basis import generate_basis_functions

        n = 6
        fns = generate_basis_functions("legendre", n, 0.0, 100.0)
        assert len(fns) == n

    def test_unknown_basis_raises(self):
        """Unknown basis type raises ValueError."""
        from sqc.reconstruction.basis import generate_basis_functions

        with pytest.raises(ValueError, match="Unsupported basis type"):
            generate_basis_functions("chebyshev", 5, 0.0, 1.0)


class TestBasisFunctionDecomposition:
    """Tests for basis_function_decomposition."""

    def test_decompose_pure_basis_signal(self):
        """A signal that is exactly a basis function decomposes with one
        non-zero coefficient."""
        from sqc.reconstruction.basis import (
            generate_basis_functions,
            basis_function_decomposition,
        )

        T = 100.0
        n = 5
        fns = generate_basis_functions("fourier", n, 0.0, T)
        t = np.linspace(0, T, 200)

        # Signal = second basis function (first sine)
        sig = fns[1](t)
        coeffs = basis_function_decomposition(sig, t, fns)

        # Coefficient for index 1 should be ~1.0, others ~0
        assert abs(coeffs[1] - 1.0) < 0.01
        assert abs(coeffs[0]) < 0.01
        if n > 2:
            assert abs(coeffs[2]) < 0.01

    def test_decompose_signal_object(self):
        """Decomposition works with objects that have .signal attribute."""
        from sqc.reconstruction.basis import (
            generate_basis_functions,
            basis_function_decomposition,
        )

        class FakeSignal:
            def __init__(self, s):
                self.signal = s

        T = 100.0
        fns = generate_basis_functions("fourier", 5, 0.0, T)
        t = np.linspace(0, T, 200)
        sig_raw = fns[1](t)
        fake = FakeSignal(sig_raw)

        coeffs_raw = basis_function_decomposition(sig_raw, t, fns)
        coeffs_obj = basis_function_decomposition(fake, t, fns)
        np.testing.assert_allclose(coeffs_raw, coeffs_obj)

    def test_decompose_linear_combination(self):
        """Linear combination of basis functions is exactly recovered."""
        from sqc.reconstruction.basis import (
            generate_basis_functions,
            basis_function_decomposition,
        )

        T = 100.0
        n = 5
        fns = generate_basis_functions("fourier", n, 0.0, T)
        t = np.linspace(0, T, 200)

        # true coefficients
        b_true = np.array([0.5, 0.0, 0.3, -0.2, 0.1])
        sig = np.zeros(len(t))
        for i in range(n):
            sig += b_true[i] * fns[i](t)

        coeffs = basis_function_decomposition(sig, t, fns)
        np.testing.assert_allclose(coeffs, b_true, atol=1e-10)


class TestRegularizationMatrix:
    """Tests for regularization_matrix and alias R."""

    def test_fourier_regularization_diagonal(self):
        """Fourier regularisation matrix is diagonal."""
        from sqc.reconstruction.basis import regularization_matrix

        R = regularization_matrix(5, "fourier")
        assert R.shape == (5, 5)
        # Check it's diagonal
        off_diag = R - np.diag(np.diag(R))
        np.testing.assert_allclose(off_diag, 0.0, atol=1e-15)

    def test_fourier_frequency_ordering(self):
        """Higher frequency basis functions are penalized more strongly."""
        from sqc.reconstruction.basis import regularization_matrix

        R = regularization_matrix(9, "fourier")
        diag = np.diag(R)
        # constant term (index 0) should have smallest penalty
        assert diag[0] < diag[1]
        # penalties should be non-decreasing in groups of 2
        for i in range(2, len(diag)):
            assert diag[i] >= diag[i - 2]

    def test_non_fourier_regularization_shape(self):
        """Non-Fourier regularisation matrix has shape (n, n)."""
        from sqc.reconstruction.basis import regularization_matrix

        for bt in ("bspline", "legendre"):
            R = regularization_matrix(6, bt)
            assert R.shape == (6, 6)

    def test_R_alias(self):
        """R is identical to regularization_matrix."""
        from sqc.reconstruction.basis import regularization_matrix, R

        R1 = regularization_matrix(7, "fourier")
        R2 = R(7, "fourier")
        np.testing.assert_allclose(R1, R2)

    def test_regularization_semidefinite(self):
        """Regularisation matrix should be positive semidefinite."""
        from sqc.reconstruction.basis import regularization_matrix

        for bt in ("bspline", "fourier", "legendre"):
            R = regularization_matrix(8, bt)
            eigvals = np.linalg.eigvalsh(R)
            assert np.all(eigvals >= -1e-12), (
                f"{bt} regularization has negative eigenvalues"
            )
