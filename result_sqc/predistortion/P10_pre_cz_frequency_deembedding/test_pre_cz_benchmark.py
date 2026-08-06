from pathlib import Path

import numpy as np
import pytest

import core as B


REPO = Path(r"C:\Users\21034\Desktop\Workspace\Sensing project\Sensing-Project")


def test_noiseless_frequency_deembedding_cancels_protocol_response():
    kernel = B.load_transient_kernel(REPO)
    xs = B.probe_set(11, 0)
    response_rec = []
    line_rec = []
    h = B.lowpass_impulse(10.0)
    for x in xs:
        d0 = B.transient_forward_periodic(x, kernel)
        dl = B.transient_forward_periodic(B.apply_lti_periodic(x, h), kernel)
        response_rec.append(B.transient_reconstruct(d0, kernel, 5.0))
        line_rec.append(B.transient_reconstruct(dl, kernel, 5.0))
    response, reliable, *_ = B.estimate_transfer(xs, response_rec)
    measured, support, _ = B.estimate_transfer(xs, line_rec)
    rel = reliable & support & (np.abs(response) >= 0.1 * np.max(np.abs(response)))
    estimate = B.deembed_transfer(measured, response, rel)
    truth = B.truth_frequency_response(h)
    assert np.max(np.abs(estimate[rel] - truth[rel])) < 2e-5


def test_fir_rejects_insufficient_band_instead_of_fallback():
    h = np.ones(B.N_SAMPLES // 2 + 1, dtype=complex)
    band = np.zeros_like(h, dtype=bool)
    band[:B.MIN_QUALIFIED_BINS] = True
    with pytest.raises(ValueError, match="insufficient contiguous qualified band"):
        B.design_inverse_fir(h, band, 64, 1e-3)


def test_contiguous_mask_stops_at_first_failed_bin():
    mask = np.zeros(B.N_SAMPLES // 2 + 1, dtype=bool)
    reliable = np.ones_like(mask)
    mask[:12] = True
    mask[5] = False
    contiguous = B.contiguous_mask(mask, reliable)
    assert np.array_equal(np.flatnonzero(contiguous), np.arange(5))


def test_zero_mean_transfer_has_explicit_unit_dc_anchor():
    xs = B.probe_set(21, 0)
    response, support, _ = B.estimate_transfer(xs, xs)
    assert response[0] == 1.0 + 0.0j
    assert support[0]


def test_filter_selection_never_returns_constraint_violating_choice():
    h_impulse = B.lowpass_impulse(5.0)
    h = B.truth_frequency_response(h_impulse)
    band = np.zeros_like(h, dtype=bool)
    band[:16] = True
    validation = B.heldout_waveforms()["validation_dual_pulse"]
    choice = B.choose_filter(h, band, h_impulse, validation)
    if choice is not None:
        assert choice.awg_peak <= 0.25
        assert choice.awg_slew <= 0.08
        assert choice.noise_gain <= 3.0
