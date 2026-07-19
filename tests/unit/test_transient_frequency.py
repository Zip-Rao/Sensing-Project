"""Unit tests for transient single-point frequency measurement (order >= 3).

Covers the Phase-12 fixes:
  - adaptive delta_max in the Route A fit (unbiased G1),
  - Route B full off-diagonal kernel (g3_source='kernel_full'),
  - the order>=3 sign fix (δω = −Δ convention),
  - removal of the diagonal-kernel option.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.config import CONFIG

# Short global time axis (pulse is ~20 ns; post-pulse detuning does not
# change p_e) keeps these mesolve-heavy tests fast.
_TR = CONFIG.pulse.t_rabi.copy()
_TG = CONFIG.pulse.make_time(0, 40)


def _q(flux: float = 0.0):
    from src.qubit import TransmonQubit

    return TransmonQubit(
        EC=0.2 * 2 * np.pi, EJ=15 * 2 * np.pi,
        T1=10000, T2=5000, n_levels=2, flux=flux,
    )


def test_g3_source_unknown_raises():
    """The removed diagonal option (or any unknown source) raises early."""
    from sqc.calibration.frequency import _measure_frequency_transient

    q = _q()
    with pytest.raises(ValueError, match="unknown g3_source"):
        _measure_frequency_transient(
            q, q.frequency, _TR, _TG, flux=0.0, order=3, g3_source="diag_legacy",
        )


def test_adaptive_delta_max_recovers_unbiased_g1():
    """Adaptive scan recovers a G1 far larger than the wide-range (0.08 GHz)
    fit, which is biased low (~0.6x) by fringe saturation."""
    from sqc.calibration.frequency import (
        _calibrate_g3_taylor, _fit_g3_at_delta_max, _g3_cache,
    )

    q = _q()
    _g3_cache.clear()
    G1_adaptive, _ = _calibrate_g3_taylor(q, q.frequency, _TR, _TG)  # None = adaptive
    G1_wide, _ = _fit_g3_at_delta_max(q, q.frequency, _TR, _TG, 21, 0.08)

    # adaptive recovers the true (steeper) near-zero slope
    assert abs(G1_adaptive) > 1.3 * abs(G1_wide)


def test_kernel_full_g3_matches_fit():
    """Route B (∭k₃) and Route A (fit) measure the same constant-Δ cubic
    object: |G3| agree within ~20%, opposite sign (δω vs Δ convention)."""
    from sqc.calibration.frequency import (
        _calibrate_g3_taylor, _calibrate_g3_kernel_full, _g3_cache,
    )

    q = _q()
    _g3_cache.clear()
    _, G3_taylor = _calibrate_g3_taylor(q, q.frequency, _TR, _TG)   # Δ convention
    _, G3_ker = _calibrate_g3_kernel_full(q, q.frequency, _TR)      # δω = −Δ convention

    assert np.sign(G3_ker) == -np.sign(G3_taylor)
    rel = abs(abs(G3_ker) - abs(G3_taylor)) / abs(G3_taylor)
    assert rel < 0.2, f"G3 kernel_full vs fit mismatch: {rel:.2f}"


def test_g3_delta_max_override_used():
    """Passing g3_delta_max overrides the adaptive search (distinct result
    cached under a distinct key)."""
    from sqc.calibration.frequency import _calibrate_g3_taylor, _g3_cache

    q = _q()
    _g3_cache.clear()
    G1_wide, _ = _calibrate_g3_taylor(q, q.frequency, _TR, _TG, delta_max_ghz=0.08)
    G1_adaptive, _ = _calibrate_g3_taylor(q, q.frequency, _TR, _TG)
    # the override (biased) and adaptive results differ and coexist in cache
    assert abs(G1_wide) < abs(G1_adaptive)
    assert len(_g3_cache) >= 2


def test_order3_sign_correct_and_beats_linear():
    """In the valid sub-fold regime, order-3 (both sources) is sign-correct
    and the fit is at least as accurate as order-1."""
    from sqc.calibration.frequency import _measure_frequency_transient, _g3_cache

    q = _q()
    wd = q.frequency
    flux = 0.025                      # Δ ≈ −7.6 MHz, inside the fold
    w_true = _q(flux).frequency
    assert w_true < wd                # below sweet spot → Δ < 0

    _g3_cache.clear()
    f1 = _measure_frequency_transient(q, wd, _TR, _TG, flux=flux, order=1)
    _g3_cache.clear()
    f3 = _measure_frequency_transient(q, wd, _TR, _TG, flux=flux, order=3, g3_source="fit")
    _g3_cache.clear()
    fk = _measure_frequency_transient(q, wd, _TR, _TG, flux=flux, order=3, g3_source="kernel_full")

    # sign-correct: all land below ω_d (toward the true, lower frequency),
    # NOT above (the old sign bug returned ω_d − Δ, i.e. above ω_d).
    assert f1 < wd and f3 < wd and fk < wd

    e1 = abs(f1 - w_true)
    e3 = abs(f3 - w_true)
    ek = abs(fk - w_true)

    # order-3 fit improves on (or matches) linear
    assert e3 <= e1 + 1e-9 * abs(wd)
    # all within a few MHz of truth
    for e in (e1, e3, ek):
        assert e / (2 * np.pi) * 1e3 < 5.0
