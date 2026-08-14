"""Shared frozen config + metrics for the D2 resolution/sensitivity block.

D2 answers two questions the existing reconstruction figures cannot:

  A) TIME RESOLUTION — how close can two field peaks be and still be
     resolved?  Not the 0.5 ns AWG sample spacing; the limit is set by the
     control kernel's effective width sigma_k.
  B) SENSITIVITY — how small a field amplitude survives binomial projection
     noise?  Reported as A_min under an explicit shot-noise model.

EVERYTHING IN THE "FROZEN" SECTION WAS FIXED BEFORE THE FORMAL SCANS from
pilot runs (see D2_RESOLUTION_SENSITIVITY_TASK_GUIDE.md sections 5.3 / 6.3).
Thresholds must not be re-tuned after seeing formal results; PARAM_HASH
below is written into every artifact so a mismatch is detectable.

Design notes
------------
* Double peaks use ``FluxSignal(type=8, signal=array)``: type=5 ties peak
  spacing and single-peak width to the same ``width`` parameter, so it cannot
  scan spacing at fixed width (guide section 5.1).
* Shot noise is sampled INDEPENDENTLY for the signal and zero-flux reference
  arms, then differenced -- not added as an unexplained Gaussian on delta_p.
* eta_phi is NOT computed: sqc.config has no readout or reset time, so the
  experimental duty cycle is unknown.  Saved as NaN (guide section 6.3).
"""
from __future__ import annotations

import hashlib
import io
import json
import contextlib
import os
import platform
import sys
from datetime import datetime

import numpy as np

# repo root + result_sqc/_common.py (two levels up)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import _common as C  # noqa: E402

SUBDIR = "reconstruction/D2_resolution_sensitivity"

# ===========================================================================
# FROZEN CONFIG — fixed from pilot runs BEFORE the formal scans.
# ===========================================================================

# --- subtask A: time resolution ---------------------------------------------
# Fixed single-peak shape; ONLY the spacing DT is scanned.
PEAK_AMPLITUDE = 0.01        # Phi_0, per-peak amplitude A
PEAK_SIGMA_NS = 5.0          # Gaussian sigma of EACH peak (never scanned)
PEAK_CENTER_NS = 100.0       # midpoint t_c between the two peaks

# Dense near the 11-16 ns boundary the pilot bracketed, sparse in the wings.
SEPARATIONS_NS = np.array(
    [6.0, 8.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
     16.0, 18.0, 20.0, 22.0, 25.0, 30.0, 35.0, 40.0])

# LM is ~300 s/point MEASURED (full density-matrix inversion; the 1-2 min in
# the project notes is optimistic), so it runs on a frozen subset rather than
# the full axis. These are DISCRETE VERIFICATION POINTS, plotted as markers
# only and never joined into a curve.
#
# Chosen from the pilot: LM at 14 ns returns a SINGLE peak (NRMSE 0.051) where
# Wiener resolves two (0.026) -- fourier n_basis=100 at lambda=100 over-smooths,
# so LM's crossing lies ABOVE Wiener's. A subset confined to 10-16 ns would sit
# entirely inside LM's unresolved regime and could only yield a lower bound, so
# the upper points exist to BRACKET LM's own crossing instead of extrapolating
# one (guide sections 4 and 5.3).
LM_SEPARATIONS_NS = np.array([10.0, 12.0, 13.0, 14.0, 16.0, 20.0, 25.0])

# resolvability criteria (guide section 5.3) — all three must hold
C_MIN = 0.20                 # min valley contrast Cv
NRMSE_MAX = 0.10             # loose sanity gate; NOT the binding criterion
PEAK_POSITION_TOL_NS = 3.0   # max |detected - true| peak position error

# peak detection
PEAK_REL_HEIGHT = 0.30       # min peak height as fraction of window max
PEAK_SEARCH_HALF_WINDOW_NS = 60.0   # search t_c +/- this; keeps peaks interior

# --- subtask B: sensitivity -------------------------------------------------
# Unit-PEAK asymmetric double-exponential template s(t); input is phi = A*s(t),
# so A is literally the peak flux excursion in Phi_0.
TEMPLATE_CENTER_NS = 60.0
TEMPLATE_FALL_NS = 20.0
TEMPLATE_RISE_NS = 10.0

# Log-spaced; pilot put the SNR=3 crossing between 1e-4 and 3e-4 at N=1e4.
AMPLITUDES = np.array([1e-5, 2e-5, 3e-5, 5e-5, 7e-5, 1e-4,
                       2e-4, 3e-4, 5e-4, 1e-3, 2e-3, 3e-3])
N_SHOTS_VALUES = np.array([1000, 10000, 100000])

N_REPEAT = 32                # seeds per (N_shot, A) for the Wiener sweep
N_REPEAT_ZERO = 64           # seeds for the zero-signal noise baseline
SEED_BASE = 20260801         # frozen master seed

SNR_THRESHOLD = 3.0          # primary detection definition
SNR_THRESHOLD_SECONDARY = 1.0    # noise-equivalent amplitude, saved as aux

# LM sensitivity subset: amplitudes straddling the Wiener crossing, at the
# middle shot count only. 3 amps x 4 seeds + 6 zero-signal seeds = 18 runs.
LM_AMPLITUDES = np.array([1e-4, 3e-4, 1e-3])
LM_N_SHOTS = 10000
LM_N_REPEAT = 4
LM_N_REPEAT_ZERO = 6

# Shot budget asymmetry, stated explicitly because it makes the two methods
# NOT directly comparable: Wiener consumes the DIFFERENTIAL observation (both
# arms), while sqc's LM reconstructor consumes the signal arm's absolute
# p_meas only. LM therefore sees half the circuit executions per point.
LM_USES_SIGNAL_ARM_ONLY = True

# eta_phi requires a full duty cycle. sqc.config defines t_rabi_duration but
# NO readout or reset time, so T_tot is unknown and eta_phi stays NaN.
ETA_PHI_COMPUTABLE = False
ETA_PHI_MISSING = ["readout_time_ns", "reset_time_ns"]


# ===========================================================================
# signal builders
# ===========================================================================

def double_peak(separation_ns, t=None, amplitude=PEAK_AMPLITUDE,
                sigma_ns=PEAK_SIGMA_NS, center_ns=PEAK_CENTER_NS):
    """Two equal Gaussians at t_c +/- separation/2, each of width sigma_ns.

    Peak WIDTH and peak SPACING are independent by construction -- the whole
    point of not using ``FluxSignal(type=5)``, which couples them.
    """
    t = C.T_LIST if t is None else np.asarray(t, dtype=float)
    c1 = center_ns - separation_ns / 2.0
    c2 = center_ns + separation_ns / 2.0
    return amplitude * (np.exp(-(t - c1) ** 2 / (2 * sigma_ns ** 2))
                        + np.exp(-(t - c2) ** 2 / (2 * sigma_ns ** 2)))


def true_peak_positions(separation_ns, center_ns=PEAK_CENTER_NS):
    """Analytic peak centres (t_c -/+ separation/2)."""
    return np.array([center_ns - separation_ns / 2.0,
                     center_ns + separation_ns / 2.0])


def template(t=None):
    """Unit-PEAK asymmetric double-exponential s(t) (max exactly 1.0)."""
    t = C.T_LIST if t is None else np.asarray(t, dtype=float)
    raw = np.exp(-(t - TEMPLATE_CENTER_NS) / TEMPLATE_FALL_NS
                 - np.exp(-(t - TEMPLATE_CENTER_NS) / TEMPLATE_RISE_NS))
    return raw / raw.max()


# ===========================================================================
# metrics — kernel width
# ===========================================================================

def kernel_effective_width(t_samples, kernel):
    """Energy-weighted centroid / width of the response kernel.

        w_i = |k_i|^2 / sum_j |k_j|^2
        t_c = sum_i t_i w_i
        sigma_k = sqrt(sum_i (t_i - t_c)^2 w_i)

    FWHM of |k(t)| is reported ONLY when |k| has an interior single maximum
    with half-crossings on BOTH sides. The sqc flux kernel is a monotone ramp
    over its 0-19 ns support (magnitude maximal at the right edge), so FWHM is
    genuinely undefined there and comes back NaN with a reason -- rather than
    quoting one arbitrarily-chosen lobe width as the effective width
    (guide section 5.2).

    Returns
    -------
    dict with keys t_centroid, sigma, fwhm, fwhm_status, n_sign_changes
    """
    t = np.asarray(t_samples, dtype=float)
    k = np.asarray(kernel, dtype=float)
    w = np.abs(k) ** 2
    total = w.sum()
    if total <= 0:
        return {"t_centroid": np.nan, "sigma": np.nan, "fwhm": np.nan,
                "fwhm_status": "zero-energy kernel", "n_sign_changes": 0}
    w = w / total
    t_c = float((t * w).sum())
    sigma = float(np.sqrt(((t - t_c) ** 2 * w).sum()))

    a = np.abs(k)
    i_max = int(np.argmax(a))
    half = a[i_max] / 2.0
    n_sign = int(np.sum(np.diff(np.sign(k[k != 0])) != 0))

    if i_max == 0 or i_max == len(a) - 1:
        return {"t_centroid": t_c, "sigma": sigma, "fwhm": np.nan,
                "fwhm_status": f"|k| maximal at support edge (index {i_max}) "
                               f"-> monotone, no interior peak; FWHM undefined",
                "n_sign_changes": n_sign}

    left = np.where(a[:i_max] <= half)[0]
    right = np.where(a[i_max:] <= half)[0]
    if len(left) == 0 or len(right) == 0:
        return {"t_centroid": t_c, "sigma": sigma, "fwhm": np.nan,
                "fwhm_status": "|k| does not cross half-max on both sides",
                "n_sign_changes": n_sign}

    t_l = np.interp(half, [a[left[-1]], a[left[-1] + 1]],
                    [t[left[-1]], t[left[-1] + 1]])
    j = i_max + right[0]
    t_r = np.interp(half, [a[j], a[j - 1]], [t[j], t[j - 1]])
    return {"t_centroid": t_c, "sigma": sigma, "fwhm": float(t_r - t_l),
            "fwhm_status": "ok (interior single maximum)",
            "n_sign_changes": n_sign}


# ===========================================================================
# metrics — double-peak resolvability
# ===========================================================================

def nrmse(recon, truth):
    """RMSE normalised by the peak-to-peak range of the truth."""
    n = min(len(recon), len(truth))
    r = np.asarray(recon, dtype=float)[:n]
    x = np.asarray(truth, dtype=float)[:n]
    rng = x.max() - x.min()
    if rng <= 0:
        return np.nan
    return float(np.sqrt(np.mean((r - x) ** 2)) / rng)


def valley_contrast(t_recon, recon, separation_ns,
                    center_ns=PEAK_CENTER_NS,
                    rel_height=PEAK_REL_HEIGHT,
                    half_window=PEAK_SEARCH_HALF_WINDOW_NS):
    """Detect the two peaks and compute valley contrast.

        Cv = 1 - x_valley / min(x_peak1, x_peak2)

    Primary Cv requires TWO detected peaks; with one peak it is NaN (an
    unresolved pair has no valley, and inventing a number there would
    manufacture a resolution claim). A secondary contrast measured at the
    KNOWN truth peak positions is also returned so panel (b) has a continuous
    curve through the unresolved regime.

    Returns
    -------
    dict: n_peaks, peak_positions (len-2, NaN-padded), peak_values,
          valley_value, cv, cv_at_truth, peak_position_error
    """
    from scipy.signal import find_peaks

    t = np.asarray(t_recon, dtype=float)
    x = np.asarray(recon, dtype=float)
    mask = np.abs(t - center_ns) <= half_window
    if not np.any(mask):
        raise ValueError("empty peak-search window")
    idx_off = int(np.argmax(mask))
    tw, xw = t[mask], x[mask]

    idx, _ = find_peaks(xw, height=rel_height * xw.max())
    out = {"n_peaks": int(len(idx)),
           "peak_positions": np.full(2, np.nan),
           "peak_values": np.full(2, np.nan),
           "valley_value": np.nan, "cv": np.nan,
           "peak_position_error": np.nan,
           "peak_indices_global": np.full(2, -1, dtype=int)}

    truth_pos = true_peak_positions(separation_ns, center_ns)

    if len(idx) >= 2:
        # outermost pair — the candidates for the two lobes
        i1, i2 = int(idx[0]), int(idx[-1])
        out["peak_positions"] = np.array([tw[i1], tw[i2]])
        out["peak_values"] = np.array([xw[i1], xw[i2]])
        out["peak_indices_global"] = np.array([i1 + idx_off, i2 + idx_off])
        valley = float(xw[i1:i2 + 1].min())
        out["valley_value"] = valley
        denom = min(xw[i1], xw[i2])
        out["cv"] = float(1.0 - valley / denom) if denom > 0 else np.nan
        out["peak_position_error"] = float(
            np.max(np.abs(out["peak_positions"] - truth_pos)))
    elif len(idx) == 1:
        i1 = int(idx[0])
        out["peak_positions"] = np.array([tw[i1], np.nan])
        out["peak_values"] = np.array([xw[i1], np.nan])
        out["peak_indices_global"] = np.array([i1 + idx_off, -1])

    # secondary: contrast between the TRUE peak locations (always defined)
    j1, j2 = (int(np.argmin(np.abs(t - p))) for p in truth_pos)
    if j2 > j1:
        v = float(x[j1:j2 + 1].min())
        d = min(x[j1], x[j2])
        out["cv_at_truth"] = float(1.0 - v / d) if d > 0 else np.nan
    else:
        out["cv_at_truth"] = 0.0
    return out


def contrast_for_crossing(cv, ran=None):
    """Turn raw Cv into a curve usable for the Delta_t_min crossing.

    A NaN Cv means "fewer than two peaks were detected", which is a genuine
    measurement -- the pair is DEFINITIVELY unresolved there -- not missing
    data. Mapping it to 0.0 is what lets the crossing be bracketed.

    Dropping those points instead (treating NaN as absent) mis-reports a method
    that fails at every small spacing: the first surviving point already sits
    above the criterion, so the crossing comes back as "always_above" at that
    point rather than bracketed between it and the last failing spacing.

    ``ran`` marks where the method actually ran; NaN is preserved there so
    never-run points stay excluded.
    """
    cv = np.asarray(cv, dtype=float)
    out = np.where(np.isfinite(cv), cv, 0.0)
    if ran is not None:
        out = np.where(np.asarray(ran, dtype=bool), out, np.nan)
    return out


def lm_ran_from_arrays(tr):
    """Which spacings LM actually produced a reconstruction for.

    Derived from the reconstruction rows themselves (a row with any finite
    sample means LM ran and returned a waveform) rather than trusting a stored
    flag. An early revision wrote ``lm_ran`` as ``isfinite(Cv)``, which wrongly
    excludes spacings where LM ran successfully but found FEWER THAN TWO PEAKS
    -- exactly the unresolved points the crossing needs. Deriving it here makes
    old and new caches read correctly.
    """
    if "lm_reconstructions" in tr:
        rec = np.asarray(tr["lm_reconstructions"], dtype=float)
        return np.any(np.isfinite(rec), axis=1)
    if "lm_ran" in tr:
        return np.asarray(tr["lm_ran"], dtype=bool)
    return None


def dt_min_from_arrays(tr):
    """Recompute all three Delta_t_min values from a loaded npz dict.

    Single source of truth shared by the generator and the plot/metrics, so
    every number in the figure is reproducible from the saved arrays alone
    (guide section 9.10) and a cache written by an older revision is corrected
    on read instead of silently carried forward.

    Returns {method: (value, status)} for 'wiener', 'lm', 'truth'.
    """
    sep = np.asarray(tr["separations_ns"], dtype=float)
    lm_ran = lm_ran_from_arrays(tr)
    return {
        "wiener": first_crossing(
            sep, contrast_for_crossing(tr["wiener_valley_contrast"]), C_MIN),
        "truth": first_crossing(
            sep, contrast_for_crossing(tr["truth_valley_contrast"]), C_MIN),
        "lm": _lm_crossing(sep, tr["lm_valley_contrast"], lm_ran),
    }


def _lm_crossing(sep, lm_cv, lm_ran):
    curve = contrast_for_crossing(lm_cv, ran=lm_ran)
    ok = np.isfinite(curve)
    if ok.sum() < 2:
        return (np.nan, "insufficient_lm_points")
    return first_crossing(sep[ok], curve[ok], C_MIN)


def is_resolved(m, nrmse_val):
    """All three frozen criteria: two peaks, Cv >= C_MIN, NRMSE <= max, pos tol."""
    return bool(
        m["n_peaks"] >= 2
        and np.isfinite(m["cv"]) and m["cv"] >= C_MIN
        and np.isfinite(nrmse_val) and nrmse_val <= NRMSE_MAX
        and np.isfinite(m["peak_position_error"])
        and m["peak_position_error"] <= PEAK_POSITION_TOL_NS)


def first_crossing(x_vals, y_vals, threshold, log=False):
    """Smallest x where y first reaches threshold, interpolated in the bracket.

    ``log=True`` interpolates in log10(x)-log10(y), used for the sensitivity
    axis. It is exact for ANY power law SNR ~ A^p, while linear-in-A is exact
    only for p == 1; since the transduction p_e(phi) is nonlinear in general,
    p need not be exactly 1. The two agree when p == 1, so nothing is lost by
    preferring log on the log-spaced amplitude grid. Left False for the
    resolution axis, where Cv(dt) is not a power law at all.

    Returns (value, status): "bracketed", "always_above" (report as an upper
    bound) or "never" (lower bound). NEVER extrapolates outside the scanned
    range (guide sections 5.3 / 6.3).
    """
    x = np.asarray(x_vals, dtype=float)
    y = np.asarray(y_vals, dtype=float)
    ok = np.isfinite(y)
    if not np.any(ok):
        return np.nan, "no_finite_data"
    x, y = x[ok], y[ok]
    if y[0] >= threshold:
        return float(x[0]), "always_above"
    above = np.where(y >= threshold)[0]
    if len(above) == 0:
        return np.nan, "never"
    i = int(above[0])
    x0, x1, y0, y1 = x[i - 1], x[i], y[i - 1], y[i]
    if y1 == y0:
        return float(x1), "bracketed"
    if log and x0 > 0 and y0 > 0 and x1 > 0 and y1 > 0:
        lx = (np.log10(x0) + (np.log10(threshold) - np.log10(y0))
              * (np.log10(x1) - np.log10(x0)) / (np.log10(y1) - np.log10(y0)))
        return float(10.0 ** lx), "bracketed"
    frac = (threshold - y0) / (y1 - y0)
    return float(x0 + frac * (x1 - x0)), "bracketed"


# ===========================================================================
# forward model
# ===========================================================================

def forward(flux_array, t=None, quiet=True):
    """One ideal-coherent transient sensing run on a custom flux array.

    Returns (exp, result). The runner prints timings; quiet=True swallows them.
    """
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.control.flux_signal import FluxSignal

    t = C.T_LIST if t is None else t
    q = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    fs = FluxSignal(type=8, t_list=t, signal=np.asarray(flux_array, float))
    exp = TransientSensingExperiment(qubit=q, flux_signal=fs)
    if quiet:
        with contextlib.redirect_stdout(io.StringIO()):
            res = exp.run()
    else:
        res = exp.run()
    return exp, res


def wiener(delta_p, kernel, scan_axis, lambda_reg=C.WIENER_LAMBDA):
    """Wiener deconvolution of a (possibly noisy) delta_p.

    Wraps the delta_p into a minimal ExperimentResult so the unmodified sqc
    reconstructor -- including its /dt normalisation, which is verified against
    src and must not be touched -- does the work.
    """
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.simulation.result import ExperimentResult

    res = ExperimentResult(
        data={"delta_p": np.asarray(delta_p, dtype=float)},
        axes={"scan": np.asarray(scan_axis, dtype=float)})
    rec = TransientReconstruction(
        method="wiener", lambda_reg=lambda_reg
    ).reconstruct(res, kernel=np.asarray(kernel, dtype=float))
    return np.asarray(rec.t_list, dtype=float), np.asarray(rec.signal, dtype=float)


def lm_reconstruct(exp, res, p_meas=None, lambda_reg=C.LM_LAMBDA,
                   n_basis=C.LM_N_BASIS, max_iter=C.LM_MAX_ITER):
    """LM full-density-matrix inversion. Returns (signal, history, status).

    On failure returns (None, None, "<exception text>") so the caller can save
    NaN plus the reason -- LM is never silently replaced by Wiener.
    ``p_meas`` overrides the ideal populations with a shot-sampled copy.
    ``exp.control_pulse`` only exists after ``exp.run()``.
    """
    from sqc.reconstruction.transient import TransientReconstruction

    adapted = C.adapt_for_lm(res)
    if p_meas is not None:
        adapted.data["p_meas"] = np.asarray(p_meas, dtype=float)
    try:
        rec = TransientReconstruction(
            method="lm", qubit=C.make_sqc_qubit(C.OPTIMAL_FLUX),
            control_pulse=exp.control_pulse, basis_type="fourier",
            n_basis=n_basis, lambda_reg=lambda_reg, max_iter=max_iter)
        with contextlib.redirect_stdout(io.StringIO()):
            out, history = rec.reconstruct(adapted, kernel=res.data["kernel"])
        return np.asarray(out.signal, dtype=float), history, "ok"
    except Exception as exc:  # noqa: BLE001 — status is recorded, not swallowed
        return None, None, f"{type(exc).__name__}: {exc}"


# ===========================================================================
# shot noise + amplitude estimation
# ===========================================================================

def make_seeds(n, offset=0):
    """Deterministic frozen seed list (saved to the npz for reproducibility)."""
    return SEED_BASE + offset + np.arange(n, dtype=np.int64)


def sample_shots(p_signal, p_reference, n_shot, seed):
    """Independent binomial sampling of the signal and reference arms.

        n_i ~ Binomial(N_shot, p_i),   p_hat_i = n_i / N_shot

    The two arms use separate child streams of one SeedSequence, so they are
    statistically independent (they are distinct circuit executions), and the
    stream assignment does not shift if sampling order changes.

    Returns (p_hat_signal, p_hat_reference, delta_p_hat).
    """
    ss = np.random.SeedSequence(int(seed))
    g_sig, g_ref = (np.random.default_rng(s) for s in ss.spawn(2))
    ps = np.clip(np.asarray(p_signal, dtype=float), 0.0, 1.0)
    pr = np.clip(np.asarray(p_reference, dtype=float), 0.0, 1.0)
    n_shot = int(n_shot)
    hs = g_sig.binomial(n_shot, ps) / n_shot
    hr = g_ref.binomial(n_shot, pr) / n_shot
    return hs, hr, hs - hr


def project_amplitude(t_recon, recon, t_template=None, tmpl=None):
    """Template-projection amplitude estimate over the common time window.

        A_hat = <phi_hat, s> / <s, s>

    Preferred over a single max sample, which is noise-dominated at low SNR.
    """
    r = np.asarray(recon, dtype=float)
    s = template(t_recon) if tmpl is None else np.asarray(tmpl, dtype=float)
    n = min(len(r), len(s))
    r, s = r[:n], s[:n]
    denom = float(np.dot(s, s))
    if denom <= 0:
        return np.nan
    return float(np.dot(r, s) / denom)


# ===========================================================================
# metadata / IO
# ===========================================================================

METADATA_NAME = "resolution_sensitivity_metadata.json"
METRICS_NAME = "resolution_sensitivity_metrics.txt"


def frozen_config(section):
    """The frozen parameter block for one section ('time_resolution'/'sensitivity')."""
    common = {
        "EC": C.EC, "EJ": C.EJ, "T1": C.T1, "T2": C.T2,
        "n_levels": C.N_LEVELS, "flux_work_point": C.OPTIMAL_FLUX,
        "t_list_start": float(C.T_LIST[0]), "t_list_stop": float(C.T_LIST[-1]),
        "t_list_n": int(len(C.T_LIST)),
        "wiener_lambda": C.WIENER_LAMBDA, "lm_lambda": C.LM_LAMBDA,
        "lm_n_basis": C.LM_N_BASIS, "lm_max_iter": C.LM_MAX_ITER,
        "lm_basis_type": "fourier",
    }
    if section == "time_resolution":
        common.update({
            "peak_amplitude": PEAK_AMPLITUDE, "peak_sigma_ns": PEAK_SIGMA_NS,
            "peak_center_ns": PEAK_CENTER_NS,
            "separations_ns": SEPARATIONS_NS.tolist(),
            "lm_separations_ns": LM_SEPARATIONS_NS.tolist(),
            "C_MIN": C_MIN, "NRMSE_MAX": NRMSE_MAX,
            "PEAK_POSITION_TOL_NS": PEAK_POSITION_TOL_NS,
            "peak_rel_height": PEAK_REL_HEIGHT,
            "peak_search_half_window_ns": PEAK_SEARCH_HALF_WINDOW_NS,
        })
    elif section == "sensitivity":
        common.update({
            "template_center_ns": TEMPLATE_CENTER_NS,
            "template_fall_ns": TEMPLATE_FALL_NS,
            "template_rise_ns": TEMPLATE_RISE_NS,
            "amplitudes": AMPLITUDES.tolist(),
            "n_shots_values": N_SHOTS_VALUES.tolist(),
            "n_repeat": N_REPEAT, "n_repeat_zero": N_REPEAT_ZERO,
            "seed_base": SEED_BASE,
            "SNR_THRESHOLD": SNR_THRESHOLD,
            "SNR_THRESHOLD_SECONDARY": SNR_THRESHOLD_SECONDARY,
            "lm_amplitudes": LM_AMPLITUDES.tolist(),
            "lm_n_shots": LM_N_SHOTS, "lm_n_repeat": LM_N_REPEAT,
            "lm_n_repeat_zero": LM_N_REPEAT_ZERO,
            "lm_uses_signal_arm_only": LM_USES_SIGNAL_ARM_ONLY,
            "eta_phi_computable": ETA_PHI_COMPUTABLE,
            "eta_phi_missing_inputs": ETA_PHI_MISSING,
        })
    else:
        raise ValueError(f"unknown section: {section}")
    return common


def param_hash(section, quick=False):
    """SHA-256 of the frozen config, so stale caches are detectable."""
    payload = {"section": section, "quick": bool(quick),
               "config": frozen_config(section)}
    blob = json.dumps(payload, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def environment():
    """Versions + platform, recorded in the metadata JSON."""
    import scipy
    import qutip
    return {"python": sys.version.split()[0], "numpy": np.__version__,
            "scipy": scipy.__version__, "qutip": qutip.__version__,
            "platform": platform.platform()}


def npz_name(stem, quick=False):
    """Formal '<stem>_sqc' vs quick '<stem>_sqc_quick' — quick never overwrites."""
    return f"{stem}_sqc_quick" if quick else f"{stem}_sqc"


def out_path(name):
    return os.path.join(C.out_dir(SUBDIR), name)


def write_metadata_section(section, payload):
    """Read-modify-write one section of the shared metadata JSON.

    Both generators write into the same file, so this merges instead of
    clobbering the other's section.
    """
    path = out_path(METADATA_NAME)
    doc = {}
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                doc = json.load(fh)
        except (json.JSONDecodeError, OSError):
            doc = {}
    doc[section] = payload
    doc["updated"] = datetime.now().isoformat(timespec="seconds")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2, sort_keys=True, default=str)
    return path


def read_metadata():
    path = out_path(METADATA_NAME)
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)







