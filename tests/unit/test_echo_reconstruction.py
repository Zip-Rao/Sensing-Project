"""Unit tests for sqc.reconstruction.echo.EchoReconstruction.

Covers the arcsin domain guard: mesolve integrator noise can push p_e
marginally outside [0, 1], and an unclipped arcsin(2*p_e - 1) would then
emit silent NaNs into the reconstructed field.
"""

import numpy as np

from sqc.reconstruction.echo import EchoReconstruction
from sqc.simulation.result import ExperimentResult


class _StubQubit:
    """Minimal qubit exposing the flux_bias / sensitivity() interface."""

    flux_bias = 0.1

    def sensitivity(self) -> float:
        return 2.0


def _reconstruct(p_e):
    recon = EchoReconstruction(qubit=_StubQubit(), t_int=10.0, k=1)
    meas = ExperimentResult(data={"p_e": np.asarray(p_e, dtype=float)})
    return recon.reconstruct(meas)


def test_reconstruct_no_nan_when_p_e_slightly_out_of_range():
    """p_e marginally outside [0, 1] must not produce NaNs."""
    p_e = np.array([1.0 + 1e-9, -1e-9, 0.5, 1.0000002, -2e-7])
    B = _reconstruct(p_e)
    assert np.all(np.isfinite(B)), f"reconstruction produced non-finite values: {B}"


def test_reconstruct_clips_to_arcsin_domain_endpoints():
    """p_e = 1 (2*p_e-1 = 1) maps to arcsin(1) = pi/2, not NaN."""
    B = _reconstruct([1.0, 0.0])
    kappa, k, t_int = 2.0, 1, 10.0
    expected = -np.array([np.pi / 2, -np.pi / 2]) / (2 * k * kappa * t_int)
    np.testing.assert_allclose(B, expected, atol=1e-12)


def test_reconstruct_matches_formula_in_valid_range():
    """Within [0, 1] the clip is a no-op and the closed form holds."""
    p_e = np.array([0.5, 0.25, 0.75])
    B = _reconstruct(p_e)
    expected = -np.arcsin(2 * p_e - 1) / (2 * 1 * 2.0 * 10.0)
    np.testing.assert_allclose(B, expected, atol=1e-12)
