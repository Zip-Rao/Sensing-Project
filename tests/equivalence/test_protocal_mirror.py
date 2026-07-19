"""Equivalence tests: sqc experiments are deterministic and self-consistent.

With the global config refactor (sqc/config.py → arange time axes),
new sqc code intentionally uses different (more physical) time grids than
the legacy src/ code.  These tests verify that the new pipeline is
internally consistent and deterministic rather than comparing against old
linspace-based outputs.
"""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import assert_array_close


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_qubit():
    """Create a fresh TransmonQubit with default test parameters."""
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


# ---------------------------------------------------------------------------
# Case 0: Rabi — self-consistency (two runs → identical)
# ---------------------------------------------------------------------------

def test_rabi_experiment_self_consistent():
    """Two RabiExperiment runs with same params produce identical output."""
    from sqc.experiments.rabi import RabiExperiment

    q1 = _make_qubit()
    r1 = RabiExperiment(qubit=q1).run()

    q2 = _make_qubit()
    r2 = RabiExperiment(qubit=q2).run()

    assert_array_close(np.array(r1.expect), np.array(r2.expect), name="rabi_expect")


# ---------------------------------------------------------------------------
# Case 1: Ramsey — self-consistency
# ---------------------------------------------------------------------------

def test_ramsey_experiment_self_consistent():
    """Two RamseyExperiment runs with same params produce identical output."""
    from sqc.experiments.ramsey import RamseyExperiment

    q1 = _make_qubit()
    r1 = RamseyExperiment(qubit=q1).run()

    q2 = _make_qubit()
    r2 = RamseyExperiment(qubit=q2).run()

    assert_array_close(r1.data["p_e"], r2.data["p_e"], name="ramsey_p_e")
    assert_array_close(r1.axes["tau"], r2.axes["tau"], name="ramsey_tau")
    assert_array_close(r1.data["flux_samples"], r2.data["flux_samples"],
                       name="ramsey_flux")


# ---------------------------------------------------------------------------
# Case 2: DiffEcho — self-consistency
# ---------------------------------------------------------------------------

def test_diff_echo_experiment_self_consistent():
    """Two DiffEchoExperiment runs with same params produce identical output."""
    from sqc.experiments.echo import DiffEchoExperiment

    q1 = _make_qubit()
    r1 = DiffEchoExperiment(qubit=q1).run()

    q2 = _make_qubit()
    r2 = DiffEchoExperiment(qubit=q2).run()

    assert_array_close(r1.data["p_e"], r2.data["p_e"], name="diffecho_p_e")
    assert_array_close(r1.data["flux_samples"], r2.data["flux_samples"],
                       name="diffecho_flux")


# ---------------------------------------------------------------------------
# Case 4: Transient sensing — self-consistency
# ---------------------------------------------------------------------------

def test_transient_experiment_self_consistent():
    """Two TransientSensingExperiment runs with same params produce identical output."""
    from sqc.experiments.transient import TransientSensingExperiment

    q1 = _make_qubit()
    r1 = TransientSensingExperiment(qubit=q1).run()

    q2 = _make_qubit()
    r2 = TransientSensingExperiment(qubit=q2).run()

    assert_array_close(r1.data["kernel"], r2.data["kernel"], name="transient_kernel")
    assert_array_close(r1.data["delta_p"], r2.data["delta_p"],
                       name="transient_delta_p")
    assert_array_close(r1.data["p_e"], r2.data["p_e"], name="transient_p_e")
    assert_array_close(r1.axes["scan"], r2.axes["scan"], name="transient_scan")
