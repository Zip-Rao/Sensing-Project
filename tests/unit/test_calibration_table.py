"""Unit tests for CalibrationTable.evaluate and inverse methods.

Tests added in P3a per phase_3_handbook.md §5.1.
"""
from __future__ import annotations

import numpy as np
import pytest


class TestCalibrationTableEvaluate:
    """Tests for CalibrationTable.evaluate()."""

    def test_evaluate_linear(self):
        """evaluate interpolates a simple linear relationship."""
        from sqc.calibration.base import CalibrationTable

        tbl = CalibrationTable(
            name="test_linear",
            kind="phi_h",
            qubit_name="Q1",
            inputs=np.array([0.0, 0.5, 1.0]),
            outputs=np.array([0.0, 1.0, 2.0]),
        )
        x = np.array([0.25, 0.75])
        y = tbl.evaluate(x)
        np.testing.assert_allclose(y, np.array([0.5, 1.5]), rtol=1e-6)

    def test_evaluate_cubic_interpolation(self):
        """evaluate handles cubic interpolation between known points."""
        from sqc.calibration.base import CalibrationTable

        tbl = CalibrationTable(
            name="test_cubic",
            kind="f_phi",
            qubit_name="Q1",
            inputs=np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
            outputs=np.sin(np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])),
        )
        x = np.array([0.1, 0.3, 0.5])
        y = tbl.evaluate(x)
        # Should be close to sin(x)
        np.testing.assert_allclose(y, np.sin(x), rtol=1e-2)

    def test_evaluate_raises_without_inputs(self):
        """evaluate raises ValueError if inputs/outputs not set."""
        from sqc.calibration.base import CalibrationTable

        tbl = CalibrationTable(name="empty")
        with pytest.raises(ValueError, match="inputs/outputs not set"):
            tbl.evaluate(np.array([0.0]))

    def test_evaluate_extrapolation(self):
        """evaluate supports extrapolation beyond input range."""
        from sqc.calibration.base import CalibrationTable

        tbl = CalibrationTable(
            name="test_extrap",
            kind="phi_h",
            qubit_name="Q1",
            inputs=np.array([0.0, 1.0, 2.0]),
            outputs=np.array([0.0, 1.0, 2.0]),
        )
        y = tbl.evaluate(np.array([-0.5, 2.5]))
        # Should extrapolate linearly-ish (cubic with fill_value='extrapolate')
        assert y.shape == (2,)
        assert not np.any(np.isnan(y))


class TestCalibrationTableInverse:
    """Tests for CalibrationTable.inverse()."""

    def test_inverse_monotonic(self):
        """inverse recovers inputs from outputs for monotonic data."""
        from sqc.calibration.base import CalibrationTable

        inputs = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
        outputs = np.array([0.0, 0.25, 1.0, 2.25, 4.0])  # y = x^2

        tbl = CalibrationTable(
            name="test_inv",
            kind="phi_h",
            qubit_name="Q1",
            inputs=inputs,
            outputs=outputs,
        )

        y_query = np.array([0.25, 1.0, 2.25])
        x_est = tbl.inverse(y_query)
        np.testing.assert_allclose(x_est, np.array([0.5, 1.0, 1.5]), rtol=0.1)

    def test_inverse_roundtrip(self):
        """evaluate(inverse(y)) ≈ y for the full range."""
        from sqc.calibration.base import CalibrationTable

        inputs = np.linspace(0.0, 1.0, 20)
        outputs = np.sin(inputs * np.pi)  # monotonic on [0, 1]

        tbl = CalibrationTable(
            name="test_rt",
            kind="phi_h",
            qubit_name="Q1",
            inputs=inputs,
            outputs=outputs,
        )

        y_test = np.linspace(0.0, np.sin(np.pi), 15)
        x_est = tbl.inverse(y_test)
        y_est = tbl.evaluate(x_est)
        np.testing.assert_allclose(y_est, y_test, rtol=0.05)

    def test_inverse_raises_without_inputs(self):
        """inverse raises ValueError if inputs/outputs not set."""
        from sqc.calibration.base import CalibrationTable

        tbl = CalibrationTable(name="empty")
        with pytest.raises(ValueError, match="inputs/outputs not set"):
            tbl.inverse(np.array([0.0]))

    def test_inverse_non_monotonic_fallback(self):
        """inverse handles flat regions by taking the first monotonic chunk."""
        from sqc.calibration.base import CalibrationTable

        # Non-strictly monotonic in the middle
        inputs = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        outputs = np.array([0.0, 1.0, 1.5, 1.5, 2.0, 3.0])  # flat at 2-3

        tbl = CalibrationTable(
            name="test_flat",
            kind="phi_h",
            qubit_name="Q1",
            inputs=inputs,
            outputs=outputs,
        )

        # Should not raise; takes monotonic subset
        y_query = np.array([0.5, 2.5])
        x_est = tbl.inverse(y_query)
        assert x_est.shape == (2,)
        assert not np.any(np.isnan(x_est))

    def test_legacy_backward_compat(self):
        """Legacy-style CalibrationTable (no inputs/outputs) still works."""
        from sqc.calibration.base import CalibrationTable

        tbl = CalibrationTable(
            name="legacy",
            parameters={"kappa": 1.5},
            uncertainties={"kappa": 0.1},
        )
        assert tbl.name == "legacy"
        assert tbl.parameters["kappa"] == 1.5
        # qubit_name defaults to empty string
        assert tbl.qubit_name == ""
