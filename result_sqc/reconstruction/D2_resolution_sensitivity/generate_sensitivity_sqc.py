#!/usr/bin/env python3
"""D2 subtask B — sensitivity from a low-amplitude scan with projection noise.

Asks the smallest field amplitude that survives binomial measurement noise.
The chain, per (N_shot, A, seed):

  1. forward model once per A  -> IDEAL populations p_signal, p_reference
     (p_reference = p_e - delta_p, the zero-flux arm the experiment already ran)
  2. sample each arm INDEPENDENTLY:  n ~ Binomial(N_shot, p),  p_hat = n/N_shot
  3. differential observation  delta_p_hat = p_hat_signal - p_hat_reference
  4. Wiener deconvolution -> phi_hat
  5. template projection    A_hat = <phi_hat, s> / <s, s>

Sensitivity is then defined against a SEPARATELY RUN zero-signal baseline:

    SNR(A) = E[A_hat | A] / sigma_Ahat_0        A_min = min{A : SNR >= SNR_th}

sigma_Ahat_0 comes from its own A=0 repetitions -- where both arms have
identical ideal populations but are still sampled independently, so the
spread is pure shot noise. It is never scavenged from the tail of a
signal-bearing run.

Deliberately NOT done here:
  * no Gaussian noise pasted onto delta_p and relabelled "projection noise";
  * FluxSignal.noise_level is left at 0 -- it is an INPUT FLUX noise, a
    different physical thing from readout projection noise;
  * eta_phi is NOT computed. sqc.config has t_rabi_duration but no readout or
    reset time, so the duty cycle T_tot is unknown; eta_phi is saved as NaN and
    solver wall-clock is never substituted for acquisition time.

LM caveat (recorded in metadata, matters when comparing the two curves): the
Wiener path consumes the DIFFERENTIAL observation (both arms), while sqc's LM
reconstructor consumes the signal arm's absolute p_meas only -- so LM sees half
the circuit executions per point. LM runs on a frozen subset; failures are
saved as NaN plus the exception text.

Outputs:  sensitivity_sqc.npz / ..._quick.npz  +  metadata "sensitivity" section

Usage:
    python result_sqc/reconstruction/D2_resolution_sensitivity/generate_sensitivity_sqc.py --quick --no-lm
    python result_sqc/reconstruction/D2_resolution_sensitivity/generate_sensitivity_sqc.py
"""
from __future__ import annotations

import os
import sys
import time
from datetime import datetime

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _d2_common as D  # noqa: E402

QUICK_AMPLITUDES = np.array([3e-5, 1e-4, 3e-4, 1e-3])
QUICK_N_SHOTS = np.array([1000, 10000])
QUICK_N_REPEAT = 4
QUICK_N_REPEAT_ZERO = 8
QUICK_LM_AMPLITUDES = np.array([3e-4])
QUICK_LM_N_REPEAT = 1
QUICK_LM_N_REPEAT_ZERO = 2


def forward_ideal(amplitudes, tmpl):
    """One forward run per amplitude (plus A=0). Returns ideal populations.

    The A=0 entry is a real, separately executed forward run -- the zero-signal
    baseline is never inferred from a signal-bearing run.
    """
    out = {}
    exp0, res0 = D.forward(0.0 * tmpl)
    pe0 = np.asarray(res0.data["p_e"], dtype=float)
    out[0.0] = {"p_sig": pe0,
                "p_ref": pe0 - np.asarray(res0.data["delta_p"], dtype=float),
                "exp": exp0, "res": res0}
    for A in amplitudes:
        exp, res = D.forward(float(A) * tmpl)
        pe = np.asarray(res.data["p_e"], dtype=float)
        out[float(A)] = {"p_sig": pe,
                         "p_ref": pe - np.asarray(res.data["delta_p"], float),
                         "exp": exp, "res": res}
    return out


def _check_populations(ideal):
    """Guide section 9.3: every ideal population must lie in [0, 1]."""
    for A, d in ideal.items():
        for key in ("p_sig", "p_ref"):
            p = d[key]
            if not np.all((p >= 0.0) & (p <= 1.0)):
                raise ValueError(
                    f"ideal population out of [0,1] at A={A}, {key}: "
                    f"min={p.min():.6g} max={p.max():.6g}")


def wiener_sweep(ideal, amplitudes, n_shots_values, n_repeat, n_repeat_zero,
                 tmpl, verbose=True):
    """Shot-sample + Wiener + template projection over the full grid."""
    n_s, n_a = len(n_shots_values), len(amplitudes)
    est = np.full((n_s, n_a, n_repeat), np.nan)
    est_peak = np.full((n_s, n_a, n_repeat), np.nan)
    rmse_w = np.full((n_s, n_a, n_repeat), np.nan)
    est_zero = np.full((n_s, n_repeat_zero), np.nan)

    seeds = D.make_seeds(n_repeat)
    seeds_zero = D.make_seeds(n_repeat_zero, offset=500_000)
    t0 = time.time()

    for si, N in enumerate(n_shots_values):
        z = ideal[0.0]
        kern = z["res"].data["kernel"]
        scan = z["res"].axes["scan"]
        # --- zero-signal baseline: its own run, its own seeds ---
        for r, sd in enumerate(seeds_zero):
            # distinct stream per (N, seed) so shot counts never share draws
            _, _, dp = D.sample_shots(z["p_sig"], z["p_ref"], N,
                                      int(sd) + 7919 * int(si))
            tr, rc = D.wiener(dp, kern, scan)
            est_zero[si, r] = D.project_amplitude(tr, rc)
        sd0 = float(np.nanstd(est_zero[si], ddof=1))

        for ai, A in enumerate(amplitudes):
            d = ideal[float(A)]
            for r, sd in enumerate(seeds):
                _, _, dp = D.sample_shots(d["p_sig"], d["p_ref"], N,
                                          int(sd) + 7919 * int(si))
                tr, rc = D.wiener(dp, kern, scan)
                est[si, ai, r] = D.project_amplitude(tr, rc)
                est_peak[si, ai, r] = float(np.max(np.abs(rc)))
                truth = float(A) * tmpl
                rmse_w[si, ai, r] = D.C.rmse(rc, truth)
            if verbose:
                m = np.nanmean(est[si, ai])
                print(f"    N={N:>7d} A={A:.1e}  A_hat={m:.3e}  "
                      f"SNR={m / sd0:6.2f}")
        if verbose:
            print(f"  N={N:>7d}: sigma_Ahat_0 = {sd0:.3e} "
                  f"(zero-signal, {n_repeat_zero} seeds)")

    return {"est": est, "est_peak": est_peak, "est_zero": est_zero,
            "rmse": rmse_w, "seeds": seeds, "seeds_zero": seeds_zero,
            "runtime_s": time.time() - t0}


def analyse(est, est_zero, amplitudes, n_boot=2000, boot_seed=12345):
    """SNR(A) = E[A_hat|A] / sigma_Ahat_0, plus A_min and its bracket status.

    Uncertainty is reported two ways (both saved; SEM is the one drawn):
      * sem   — standard error of the mean over the frozen seeds
      * boot  — 95% percentile bootstrap CI of the mean, fixed bootstrap seed
    """
    n_s, n_a, _ = est.shape
    mean = np.nanmean(est, axis=2)
    sd_a = np.nanstd(est, axis=2, ddof=1)
    n_eff = np.sum(np.isfinite(est), axis=2)
    sem = sd_a / np.sqrt(np.maximum(n_eff, 1))

    sigma0 = np.nanstd(est_zero, axis=1, ddof=1)
    mean0 = np.nanmean(est_zero, axis=1)

    snr = mean / sigma0[:, None]
    snr_sem = sem / sigma0[:, None]

    rng = np.random.default_rng(boot_seed)
    snr_lo = np.full((n_s, n_a), np.nan)
    snr_hi = np.full((n_s, n_a), np.nan)
    for si in range(n_s):
        for ai in range(n_a):
            v = est[si, ai][np.isfinite(est[si, ai])]
            if len(v) < 2:
                continue
            bs = rng.choice(v, size=(n_boot, len(v)), replace=True).mean(axis=1)
            snr_lo[si, ai], snr_hi[si, ai] = np.percentile(
                bs / sigma0[si], [2.5, 97.5])

    a_min, a_min_status = [], []
    a_ne, a_ne_status = [], []
    for si in range(n_s):
        v, st = D.first_crossing(amplitudes, snr[si], D.SNR_THRESHOLD,
                                 log=True)
        a_min.append(v); a_min_status.append(st)
        v2, st2 = D.first_crossing(amplitudes, snr[si],
                                   D.SNR_THRESHOLD_SECONDARY, log=True)
        a_ne.append(v2); a_ne_status.append(st2)

    return {"mean": mean, "sd": sd_a, "sem": sem, "n_eff": n_eff,
            "sigma0": sigma0, "mean0": mean0, "snr": snr, "snr_sem": snr_sem,
            "snr_ci_lo": snr_lo, "snr_ci_hi": snr_hi,
            "a_min": np.array(a_min, dtype=float), "a_min_status": a_min_status,
            "a_ne": np.array(a_ne, dtype=float), "a_ne_status": a_ne_status}


def lm_subset(ideal, amplitudes, lm_amplitudes, lm_n_shots, lm_n_repeat,
              lm_n_repeat_zero, tmpl, n_shots_values, verbose=True):
    """LM on a frozen low-amplitude subset, at ONE shot count.

    LM consumes the signal arm's absolute p_meas (see module docstring), so it
    sees half the circuit executions of the Wiener differential path. Failures
    are stored as NaN + the exception text; Wiener is never substituted.
    """
    n_s, n_a = len(n_shots_values), len(amplitudes)
    est = np.full((n_s, n_a, max(lm_n_repeat, 1)), np.nan)
    est_zero = np.full((n_s, max(lm_n_repeat_zero, 1)), np.nan)
    rmse_lm = np.full((n_s, n_a, max(lm_n_repeat, 1)), np.nan)
    status: dict[str, str] = {}
    ran = np.zeros((n_s, n_a), dtype=bool)

    si_list = np.where(np.isclose(np.asarray(n_shots_values, float),
                                  float(lm_n_shots)))[0]
    if len(si_list) == 0:
        return {"est": est, "est_zero": est_zero, "rmse": rmse_lm,
                "status": {"all": f"lm_n_shots={lm_n_shots} not in sweep"},
                "ran": ran, "runtime_s": 0.0, "shot_index": -1}
    si = int(si_list[0])
    t0 = time.time()

    seeds = D.make_seeds(lm_n_repeat, offset=900_000)
    seeds_zero = D.make_seeds(lm_n_repeat_zero, offset=950_000)

    # zero-signal LM baseline — its own runs, its own seeds
    z = ideal[0.0]
    for r, sd in enumerate(seeds_zero):
        hs, _, _ = D.sample_shots(z["p_sig"], z["p_ref"], lm_n_shots, int(sd))
        sig, _h, st = D.lm_reconstruct(z["exp"], z["res"], p_meas=hs)
        status[f"zero_seed{int(sd)}"] = st
        if sig is not None:
            est_zero[si, r] = D.project_amplitude(D.C.T_LIST[:len(sig)], sig)
        if verbose:
            print(f"    LM zero seed={int(sd)}: {st}  "
                  f"A_hat={est_zero[si, r]:.3e}")

    for A in lm_amplitudes:
        ai = int(np.argmin(np.abs(np.asarray(amplitudes, float) - float(A))))
        d = ideal[float(A)]
        ran[si, ai] = True
        for r, sd in enumerate(seeds):
            hs, _, _ = D.sample_shots(d["p_sig"], d["p_ref"], lm_n_shots,
                                      int(sd))
            sig, _h, st = D.lm_reconstruct(d["exp"], d["res"], p_meas=hs)
            status[f"A{A:g}_seed{int(sd)}"] = st
            if sig is not None:
                est[si, ai, r] = D.project_amplitude(
                    D.C.T_LIST[:len(sig)], sig)
                rmse_lm[si, ai, r] = D.C.rmse(sig, float(A) * tmpl)
            if verbose:
                print(f"    LM A={A:.1e} seed={int(sd)}: {st}  "
                      f"A_hat={est[si, ai, r]:.3e}")

    return {"est": est, "est_zero": est_zero, "rmse": rmse_lm,
            "status": status, "ran": ran, "runtime_s": time.time() - t0,
            "shot_index": si}


def save(sweep, ana, lm, amplitudes, n_shots_values, tmpl, ideal,
         n_repeat, n_repeat_zero, quick=False):
    """Write the npz + the 'sensitivity' metadata section."""
    n_s = len(n_shots_values)
    lm_snr = np.full((n_s, len(amplitudes)), np.nan)
    lm_a_min = np.full(n_s, np.nan)
    lm_sigma0 = np.full(n_s, np.nan)
    lm_status_note = "not run"

    if lm is not None:
        lm_sigma0 = np.nanstd(lm["est_zero"], axis=1, ddof=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            lm_snr = np.nanmean(lm["est"], axis=2) / lm_sigma0[:, None]
        n_pts = int(np.sum(np.any(lm["ran"], axis=0)))
        lm_status_note = (
            f"discrete spot checks at {n_pts} amplitude(s), "
            f"N_shot={D.LM_N_SHOTS if not quick else int(n_shots_values[-1])}; "
            "too few seeds for a calibrated A_min, so lm_a_min stays NaN")

    arrays = {
        "t_flux": D.C.T_LIST,
        "template": tmpl,
        "amplitudes": np.asarray(amplitudes, dtype=float),
        "n_shots_values": np.asarray(n_shots_values, dtype=float),
        "seeds": sweep["seeds"],
        "seeds_zero": sweep["seeds_zero"],
        "amplitude_estimates_wiener": sweep["est"],
        "amplitude_estimates_wiener_zero": sweep["est_zero"],
        "amplitude_estimates_wiener_peak": sweep["est_peak"],
        "amplitude_estimates_lm": lm["est"] if lm else
            np.full((n_s, len(amplitudes), 1), np.nan),
        "amplitude_estimates_lm_zero": lm["est_zero"] if lm else
            np.full((n_s, 1), np.nan),
        "lm_ran": lm["ran"] if lm else np.zeros((n_s, len(amplitudes)), bool),
        "mean_wiener": ana["mean"], "sd_wiener": ana["sd"],
        "sem_wiener": ana["sem"], "n_eff_wiener": ana["n_eff"],
        "sigma_Ahat_zero_wiener": ana["sigma0"],
        "mean_Ahat_zero_wiener": ana["mean0"],
        "sigma_Ahat_zero_lm": lm_sigma0,
        "snr_wiener": ana["snr"], "snr_sem_wiener": ana["snr_sem"],
        "snr_ci_lo_wiener": ana["snr_ci_lo"],
        "snr_ci_hi_wiener": ana["snr_ci_hi"],
        "snr_lm": lm_snr,
        "a_min_wiener": ana["a_min"], "a_min_lm": lm_a_min,
        "a_noise_equivalent_wiener": ana["a_ne"],
        # eta_phi requires a duty cycle sqc does not define -> NaN, by design
        "eta_phi_wiener": np.full(n_s, np.nan),
        "eta_phi_lm": np.full(n_s, np.nan),
        "rmse_wiener": sweep["rmse"],
        "rmse_lm": lm["rmse"] if lm else
            np.full((n_s, len(amplitudes), 1), np.nan),
        "runtime_wiener_s": np.array(sweep["runtime_s"]),
        "runtime_lm_s": np.array(lm["runtime_s"] if lm else 0.0),
        "SNR_THRESHOLD": np.array(D.SNR_THRESHOLD),
        "p_ideal_zero_signal": ideal[0.0]["p_sig"],
        "p_ideal_zero_reference": ideal[0.0]["p_ref"],
        "scan_axis": np.asarray(ideal[0.0]["res"].axes["scan"], dtype=float),
        "kernel": np.asarray(ideal[0.0]["res"].data["kernel"], dtype=float),
    }
    path = D.C.save_npz(D.SUBDIR, D.npz_name("sensitivity", quick), **arrays)

    meta = {
        "is_quick_run": bool(quick),
        "param_hash": D.param_hash("sensitivity", quick),
        "generated": datetime.now().isoformat(timespec="seconds"),
        "config": D.frozen_config("sensitivity"),
        "amplitudes_used": np.asarray(amplitudes, float).tolist(),
        "n_shots_used": np.asarray(n_shots_values, float).tolist(),
        "n_repeat_used": int(n_repeat),
        "n_repeat_zero_used": int(n_repeat_zero),
        "seeds_used": sweep["seeds"].tolist(),
        "seeds_zero_used": sweep["seeds_zero"].tolist(),
        "noise_model": {
            "kind": "binomial projection (shot) noise",
            "formula": "n_i ~ Binomial(N_shot, p_i); p_hat_i = n_i / N_shot",
            "arms_sampled_independently": True,
            "differential": "delta_p_hat = p_hat_signal - p_hat_reference",
            "input_flux_noise_level": 0.0,
            "T1_ns": D.C.T1, "T2_ns": D.C.T2,
            "dissipation_note": ("T1/T2 act inside the ideal forward model; "
                                 "they are separate from the projection noise "
                                 "added at the measurement step"),
        },
        "uncertainty": {
            "primary": "standard error of the mean over frozen seeds",
            "also_saved": "95% percentile bootstrap CI of the mean",
            "bootstrap_resamples": 2000, "bootstrap_seed": 12345,
        },
        "a_min": {
            "threshold": D.SNR_THRESHOLD,
            "interpolation": ("log10-log10 within the bracketing interval: "
                              "exact for any power law SNR ~ A^p, whereas "
                              "linear-in-A is exact only for p==1. Never "
                              "extrapolated outside the scanned grid."),
            "estimator_gain_caveat": (
                "A_min is defined on the ESTIMATOR mean, and the Wiener "
                "estimator over-recovers amplitude (Ahat/A settles near 1.2 "
                "at the top of the scanned range; see mean_wiener). A_min is "
                "therefore the detection threshold of THIS estimator, not a "
                "gain-corrected physical amplitude -- dividing by the "
                "large-A gain would shift it by that factor."),
            "per_n_shot": [
                {"n_shot": float(N), "a_min": None if not np.isfinite(v) else v,
                 "status": st}
                for N, v, st in zip(n_shots_values, ana["a_min"],
                                    ana["a_min_status"])],
            "secondary_threshold": D.SNR_THRESHOLD_SECONDARY,
            "noise_equivalent_amplitude_per_n_shot": [
                {"n_shot": float(N), "a": None if not np.isfinite(v) else v,
                 "status": st}
                for N, v, st in zip(n_shots_values, ana["a_ne"],
                                    ana["a_ne_status"])],
        },
        "eta_phi": {
            "value": None, "computable": False,
            "missing_inputs": D.ETA_PHI_MISSING,
            "reason": ("sqc.config defines pulse.t_rabi_duration but no "
                       "readout or reset time, so the acquisition duty cycle "
                       "T_tot is unknown. Only A_min under the stated "
                       "shot-noise model is reported. Solver wall-clock is "
                       "NOT a substitute for acquisition time."),
            "known_inputs": {
                "t_rabi_duration_ns": 10.0,
                "n_delay_points": int(len(arrays["scan_axis"])),
                "circuit_branches_per_delay": 2,
                "branch_note": "signal arm + zero-flux reference arm",
                "n_shot": np.asarray(n_shots_values, float).tolist(),
            },
        },
        "lm": {"status_note": lm_status_note,
               "uses_signal_arm_only": D.LM_USES_SIGNAL_ARM_ONLY,
               "shot_budget_caveat": (
                   "Wiener uses the differential (both arms); sqc's LM "
                   "reconstructor consumes the signal arm's absolute p_meas "
                   "only, i.e. half the circuit executions per point. The two "
                   "SNR curves are therefore not a like-for-like comparison."),
               "status_per_point": lm["status"] if lm else {},
               "runtime_s": lm["runtime_s"] if lm else 0.0},
        "runtime_s": {"wiener": sweep["runtime_s"],
                      "lm": lm["runtime_s"] if lm else 0.0},
        "notes": [
            "Zero-signal baseline is a SEPARATE forward run with its own "
            "seeds; at A=0 both arms share the same ideal populations but are "
            "sampled independently, so the spread is pure shot noise.",
            "Amplitude is estimated by template projection, not by a single "
            "max sample; the peak estimate is saved as an auxiliary only.",
            "A is the PEAK flux excursion in Phi_0 (template has max 1.0).",
            "kappa_phi (dOmega/dPhi) is a transduction slope, NOT a "
            "sensitivity; it is deliberately not reported as one.",
        ],
        "environment": D.environment(),
    }
    D.write_metadata_section("sensitivity", meta)
    return path, meta


def main():
    quick = "--quick" in sys.argv
    with_lm = "--no-lm" not in sys.argv
    force = "--force" in sys.argv

    if quick:
        amps, shots = QUICK_AMPLITUDES, QUICK_N_SHOTS
        n_rep, n_rep0 = QUICK_N_REPEAT, QUICK_N_REPEAT_ZERO
        lm_amps, lm_rep, lm_rep0 = (QUICK_LM_AMPLITUDES, QUICK_LM_N_REPEAT,
                                    QUICK_LM_N_REPEAT_ZERO)
        lm_shots = int(QUICK_N_SHOTS[-1])
    else:
        amps, shots = D.AMPLITUDES, D.N_SHOTS_VALUES
        n_rep, n_rep0 = D.N_REPEAT, D.N_REPEAT_ZERO
        lm_amps, lm_rep, lm_rep0 = (D.LM_AMPLITUDES, D.LM_N_REPEAT,
                                    D.LM_N_REPEAT_ZERO)
        lm_shots = D.LM_N_SHOTS

    target = D.out_path(D.npz_name("sensitivity", quick) + ".npz")
    if os.path.exists(target) and not force:
        raise SystemExit(
            f"{os.path.basename(target)} already exists; pass --force to "
            f"overwrite (quick and formal caches never collide)")

    kind = "QUICK SMOKE (not a formal result)" if quick else "FORMAL"
    print(f"[D2-B sensitivity] {kind}: {len(amps)} amplitudes x "
          f"{len(shots)} shot counts x {n_rep} seeds "
          f"(+{n_rep0} zero-signal seeds)")
    print(f"  frozen: SNR_th={D.SNR_THRESHOLD}  seed_base={D.SEED_BASE}  "
          f"noise=binomial(N_shot,p), arms sampled independently")
    print(f"  eta_phi: NOT computed (missing {', '.join(D.ETA_PHI_MISSING)})")

    tmpl = D.template()
    t0 = time.time()
    print(f"  forward model: {len(amps) + 1} runs (incl. zero-signal)...")
    ideal = forward_ideal(amps, tmpl)
    _check_populations(ideal)
    print(f"    ideal populations in [0,1] OK ({time.time() - t0:.1f} s)")

    sweep = wiener_sweep(ideal, amps, shots, n_rep, n_rep0, tmpl)
    ana = analyse(sweep["est"], sweep["est_zero"], amps)

    lm = None
    if with_lm:
        print(f"  LM subset ({len(lm_amps)} amp x {lm_rep} seed + {lm_rep0} "
              f"zero) at N={lm_shots} — slow, ~5 min/call:")
        lm = lm_subset(ideal, amps, lm_amps, lm_shots, lm_rep, lm_rep0,
                       tmpl, shots)

    path, meta = save(sweep, ana, lm, amps, shots, tmpl, ideal,
                      n_rep, n_rep0, quick=quick)

    print("\n  A_min (SNR >= %.1f):" % D.SNR_THRESHOLD)
    for row in meta["a_min"]["per_n_shot"]:
        val = ("%.3e" % row["a_min"]) if row["a_min"] else "not bracketed"
        print(f"    N={row['n_shot']:>9.0f}  A_min={val}  ({row['status']})")
    print(f"  eta_phi = NaN by design ({', '.join(D.ETA_PHI_MISSING)} undefined)")
    print(f"  total {time.time() - t0:.1f} s -> {path}")


if __name__ == "__main__":
    main()



