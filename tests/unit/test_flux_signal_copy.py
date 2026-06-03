"""tests.unit.test_flux_signal_copy — regression test for FluxSignal.copy().

Ensures copy() preserves in-place mutations (truncate, etc.).
"""
import numpy as np
from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal


def test_copy_preserves_truncation():
    """copy() must preserve truncation state, not regenerate from _params."""
    t = CONFIG.pulse.make_time(0, 80)
    sig = np.ones(len(t))
    orig = FluxSignal(type=8, t_list=t, signal=sig)

    copied = orig.copy()
    copied.truncate(10.0, 50.0)

    # The copy should be truncated
    assert np.all(copied.samples[t < 10.0] == 0.0), "t<10 not zeroed"
    assert np.all(copied.samples[t > 50.0] == 0.0), "t>50 not zeroed"
    assert np.all(copied.samples[(t >= 10.0) & (t <= 50.0)] == 1.0), \
        "inner region corrupted"

    # The original must NOT be affected
    assert np.all(orig.samples == 1.0), "original was mutated by copy+truncate"


def test_copy_preserves_update_signal():
    """copy() must preserve update_signal() changes."""
    t = CONFIG.pulse.make_time(0, 80)
    sig = np.zeros(len(t))
    orig = FluxSignal(type=8, t_list=t, signal=sig)

    orig.update_signal(offset=3.0)
    copied = orig.copy()

    assert np.allclose(copied.samples, 3.0), "copy lost update_signal changes"
    assert np.allclose(orig.samples, 3.0), "original corrupted"
    # _params should still have original signal
    assert np.allclose(orig._params["signal"], 0.0), \
        "_params signal was mutated"


def test_chained_copy_truncate():
    """Multiple copy+truncate cycles must not interfere."""
    t = CONFIG.pulse.make_time(0, 80)
    sig = np.ones(len(t))
    orig = FluxSignal(type=8, t_list=t, signal=sig)

    a = orig.copy()
    a.truncate(0, 60.0)
    b = a.copy()
    b.truncate(0, 40.0)
    c = orig.copy()
    c.truncate(0, 20.0)

    # All should be independent
    assert np.all(orig.samples == 1.0), "original mutated"
    assert np.all(a.samples[t > 60.0] == 0.0), "a: t>60 not zeroed"
    assert np.all(a.samples[t <= 60.0] == 1.0), "a: t<=60 corrupted"
    assert np.all(b.samples[t > 40.0] == 0.0), "b: t>40 not zeroed"
    assert np.all(b.samples[t <= 40.0] == 1.0), "b: t<=40 corrupted"
    assert np.all(c.samples[t > 20.0] == 0.0), "c: t>20 not zeroed"
    assert np.all(c.samples[t <= 20.0] == 1.0), "c: t<=20 corrupted"
