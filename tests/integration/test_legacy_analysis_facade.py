"""Integration tests for the legacy Analysis facade."""
from __future__ import annotations

import numpy as np
import pytest
from qutip import basis, mesolve, qeye

from sqc.simulation.result import ExperimentResult


pytestmark = pytest.mark.integration


def test_analysis_extracts_expectation_values():
    from src.analysis import Analysis

    t = np.linspace(0, 1, 5)
    psi0 = basis(2, 0)
    result = mesolve(qeye(2), psi0, t, [], e_ops=[basis(2, 0) * basis(2, 0).dag()])

    np.testing.assert_allclose(Analysis().get_expectation_values(result, 0), np.ones(5))


def test_analysis_wiener_facade_matches_native_result():
    from src.analysis import Analysis
    from sqc.reconstruction.wiener import WienerReconstruction

    delta_p = np.array([0.0, 1.0, 2.0, 1.0])
    kernel = np.array([1.0, 0.5])
    dt = 0.5
    t_legacy, b_legacy = Analysis().wiener_deconvolution(delta_p, kernel, dt, 0.1)
    measurement = ExperimentResult(
        data={"delta_p": delta_p},
        axes={"scan": np.arange(len(delta_p)) * dt},
    )
    native = WienerReconstruction(lambda_reg=0.1).reconstruct(measurement, kernel, dt=dt)

    np.testing.assert_allclose(t_legacy, native.t_list)
    np.testing.assert_allclose(b_legacy, native.signal)
