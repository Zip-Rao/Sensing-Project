"""Integration test: TransientSensingExperiment matches baseline."""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import assert_array_close, load_baseline


def _make_qubit():
    from src.qubit import TransmonQubit
    return TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, state=0, n_levels=2,
    )


@pytest.mark.slow
def test_transient_experiment_matches_baseline():
    """TransientSensingExperiment.run() matches transient_default baseline."""
    from sqc.experiments.transient import TransientSensingExperiment

    q = _make_qubit()
    exp = TransientSensingExperiment(qubit=q)
    result = exp.run()

    bl = load_baseline("transient_default")
    assert_array_close(result.data["kernel"], bl["kernel"], name="kernel")
    assert_array_close(result.data["delta_p"], bl["delta_p"], name="delta_p")
    assert_array_close(result.data["p_e"], bl["p_e"], name="p_e")


@pytest.mark.slow
def test_transient_experiment_self_consistent():
    """Two TransientSensingExperiment runs produce identical results."""
    from sqc.experiments.transient import TransientSensingExperiment

    q1 = _make_qubit()
    r1 = TransientSensingExperiment(qubit=q1).run()

    q2 = _make_qubit()
    r2 = TransientSensingExperiment(qubit=q2).run()

    assert_array_close(r1.data["kernel"], r2.data["kernel"], name="kernel")
    assert_array_close(r1.data["delta_p"], r2.data["delta_p"], name="delta_p")
    assert_array_close(r1.data["p_e"], r2.data["p_e"], name="p_e")


@pytest.mark.slow
def test_transient_sense_then_wiener_reconstruct_is_sane():
    """End-to-end sense -> Wiener reconstruct produces a finite, non-trivial,
    bounded signal.

    Guards the real failure modes of the reconstruction pipeline: a flatlined
    ~0 output (cf. handoff issue #24), NaNs, or a diverging inverse. It does
    NOT assert tight amplitude accuracy — recovery on the default (short) grid
    is intentionally modest (peak is a few % of the true amplitude).
    """
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.reconstruction.transient import TransientReconstruction

    q = _make_qubit()
    res = TransientSensingExperiment(qubit=q).run()
    rec = TransientReconstruction(method="wiener").reconstruct(
        res, kernel=res.data["kernel"]
    )
    sig = np.asarray(rec.signal, dtype=float)

    assert sig.size > 0
    assert np.all(np.isfinite(sig)), "reconstruction produced non-finite values"
    peak = float(np.max(np.abs(sig)))
    # true signal amplitude ~0.01 Phi_0: non-trivial but not diverging.
    assert 1e-5 < peak < 1e-1, f"reconstruction peak {peak:.2e} out of sane range"
