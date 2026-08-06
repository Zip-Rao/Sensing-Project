#!/usr/bin/env python3
"""D2 subtask A — time resolution from a double-peak spacing scan.

Scans ONLY the spacing between two Gaussian field peaks, at fixed per-peak
width / amplitude / centre / work point / time grid, and asks at what spacing
the reconstruction still shows two peaks separated by a real valley.

Why this exists: the existing `different_signals/` double-peak result is a
single fixed spacing, which cannot support a resolution claim, and the 0.5 ns
AWG sample interval is NOT a time resolution. The limit here is set by the
control kernel's effective width sigma_k.

Two reference levels are reported, and keeping them apart is the point:
  * TRUTH ceiling — two sigma=5 ns Gaussians only develop a valley at all
    for spacing > 2*sigma = 10 ns. Below that the INPUT is single-peaked, so
    no sensor could resolve it. The truth's own Cv curve is saved so the
    figure can show this geometric floor.
  * WIENER / LM Δt_min — the spacing at which each reconstruction meets all
    three frozen criteria (two peaks, Cv >= C_MIN, NRMSE <= NRMSE_MAX,
    peak-position error <= tol).

Wiener runs on every spacing. LM runs on a frozen subset near the crossing
(LM_SEPARATIONS_NS) and is stored as DISCRETE verification points with NaN
elsewhere -- never back-filled from Wiener.

Outputs (in this directory):
    time_resolution_sqc.npz  /  ..._quick.npz
    resolution_sensitivity_metadata.json   ("time_resolution" section)

Usage:
    python result_sqc/reconstruction/D2_resolution_sensitivity/generate_time_resolution_sqc.py --quick --no-lm
    python result_sqc/reconstruction/D2_resolution_sensitivity/generate_time_resolution_sqc.py
    python result_sqc/reconstruction/D2_resolution_sensitivity/generate_time_resolution_sqc.py --force
"""
from __future__ import annotations

import os
import sys
import time
from datetime import datetime

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _d2_common as D  # noqa: E402

QUICK_SEPARATIONS_NS = np.array([10.0, 14.0, 20.0, 30.0])
QUICK_LM_SEPARATIONS_NS = np.array([14.0])


def run_scan(separations, lm_separations, with_lm=True, verbose=True):
    """Forward + Wiener at every spacing; LM on the frozen subset."""
    n_sep = len(separations)
    n_t = len(D.C.T_LIST)

    truth_all = np.full((n_sep, n_t), np.nan)
    wiener_all = np.full((n_sep, n_t), np.nan)
    lm_all = np.full((n_sep, n_t), np.nan)     # NaN where LM did not run
    t_recon = None
    kernel = kernel_t = None
    kw = None

    rows = []
    lm_status: dict[str, str] = {}
    t_wall = {"wiener_s": 0.0, "lm_s": 0.0}

    for i, sep in enumerate(separations):
        flux = D.double_peak(sep)
        t0 = time.time()
        exp, res = D.forward(flux)
        truth = np.asarray(res.data["flux_samples"], dtype=float)

        if kernel is None:
            kernel = np.asarray(res.data["kernel"], dtype=float)
            kernel_t = np.asarray(res.axes["t_samples"], dtype=float)
            kw = D.kernel_effective_width(kernel_t, kernel)

        tr, rc = D.wiener(res.data["delta_p"], res.data["kernel"],
                          res.axes["scan"])
        t_wall["wiener_s"] += time.time() - t0
        t_recon = tr if t_recon is None else t_recon

        truth_all[i, :len(truth)] = truth
        wiener_all[i, :len(rc)] = rc

        mw = D.valley_contrast(tr, rc, sep)
        ew = D.nrmse(rc, truth)
        # truth ceiling: what the noiseless INPUT itself can show
        mt = D.valley_contrast(np.asarray(res.axes["t_flux"], float), truth, sep)

        row = {"sep": float(sep), "w": mw, "w_nrmse": ew,
               "w_resolved": D.is_resolved(mw, ew), "truth": mt,
               "lm": None, "lm_nrmse": np.nan, "lm_resolved": False}

        # --- LM on the frozen subset (discrete verification points) ---
        if with_lm and np.any(np.isclose(lm_separations, sep)):
            t1 = time.time()
            sig, _hist, status = D.lm_reconstruct(exp, res)
            t_wall["lm_s"] += time.time() - t1
            lm_status[f"{sep:g}"] = status
            if sig is not None:
                lm_all[i, :len(sig)] = sig
                mlm = D.valley_contrast(
                    np.asarray(res.axes["t_flux"], float)[:len(sig)], sig, sep)
                elm = D.nrmse(sig, truth)
                row.update({"lm": mlm, "lm_nrmse": elm,
                            "lm_resolved": D.is_resolved(mlm, elm)})

        rows.append(row)
        if verbose:
            cv = row["w"]["cv"]
            extra = ""
            if row["lm"] is not None:
                extra = f"  | LM Cv={row['lm']['cv']:.3f} res={row['lm_resolved']}"
            elif f"{sep:g}" in lm_status:
                extra = f"  | LM FAILED: {lm_status[f'{sep:g}']}"
            print(f"  sep={sep:5.1f}  npk={row['w']['n_peaks']}  "
                  f"Cv={cv:.4f}" if np.isfinite(cv) else
                  f"  sep={sep:5.1f}  npk={row['w']['n_peaks']}  Cv=  nan ",
                  end="")
            print(f"  Cv_truth={row['truth']['cv_at_truth']:.4f}"
                  f"  NRMSE={ew:.4f}  resolved={row['w_resolved']}{extra}")

    return {"rows": rows, "truth_all": truth_all, "wiener_all": wiener_all,
            "lm_all": lm_all, "t_recon": t_recon, "kernel": kernel,
            "kernel_t": kernel_t, "kw": kw, "lm_status": lm_status,
            "runtime": t_wall}


def _stack(rows, getter, fill=np.nan, n=1):
    """Column-stack a per-row quantity into an array, NaN where absent."""
    out = np.full((len(rows), n) if n > 1 else (len(rows),), fill, dtype=float)
    for i, r in enumerate(rows):
        v = getter(r)
        if v is None:
            continue
        out[i] = v
    return out


def save(scan, separations, lm_separations, quick=False):
    """Write the npz + the 'time_resolution' metadata section."""
    rows = scan["rows"]
    seps = np.asarray(separations, dtype=float)

    w_cv = _stack(rows, lambda r: r["w"]["cv"])
    lm_cv = _stack(rows, lambda r: r["lm"]["cv"] if r["lm"] else None)
    w_nr = _stack(rows, lambda r: r["w_nrmse"])
    lm_nr = _stack(rows, lambda r: r["lm_nrmse"])

    # Delta_t_min per method, plus the truth's own geometric ceiling.
    # NaN Cv (<2 peaks) is mapped to 0 -- definitively unresolved, not missing
    # -- so the crossing brackets properly. See D.contrast_for_crossing.
    lm_ran = np.array([np.any(np.isclose(np.asarray(lm_separations, float), s))
                       for s in seps], dtype=bool)
    crossings = D.dt_min_from_arrays({
        "separations_ns": seps,
        "wiener_valley_contrast": w_cv,
        "truth_valley_contrast": _stack(
            rows, lambda r: r["truth"]["cv_at_truth"]),
        "lm_valley_contrast": lm_cv,
        "lm_ran": lm_ran,
    })
    dtmin_w, st_w = crossings["wiener"]
    dtmin_t, st_t = crossings["truth"]
    dtmin_lm, st_lm = crossings["lm"]
    ok = lm_ran

    kw = scan["kw"]
    arrays = {
        "t_flux": D.C.T_LIST,
        "t_recon": scan["t_recon"],
        "kernel_time": scan["kernel_t"],
        "kernel": scan["kernel"],
        "separations_ns": seps,
        "lm_separations_ns": np.asarray(lm_separations, dtype=float),
        "truth_waveforms": scan["truth_all"],
        "wiener_reconstructions": scan["wiener_all"],
        "lm_reconstructions": scan["lm_all"],
        "wiener_nrmse": w_nr,
        "lm_nrmse": lm_nr,
        "wiener_valley_contrast": w_cv,
        "lm_valley_contrast": lm_cv,
        "truth_valley_contrast": _stack(
            rows, lambda r: r["truth"]["cv_at_truth"]),
        "wiener_peak_positions": np.vstack(
            [r["w"]["peak_positions"] for r in rows]),
        "lm_peak_positions": np.vstack(
            [r["lm"]["peak_positions"] if r["lm"] else np.full(2, np.nan)
             for r in rows]),
        "truth_peak_positions": np.vstack(
            [D.true_peak_positions(s) for s in seps]),
        "wiener_n_peaks": _stack(rows, lambda r: r["w"]["n_peaks"]),
        "lm_n_peaks": _stack(rows, lambda r: r["lm"]["n_peaks"] if r["lm"] else None),
        "wiener_peak_position_error": _stack(
            rows, lambda r: r["w"]["peak_position_error"]),
        "lm_peak_position_error": _stack(
            rows, lambda r: r["lm"]["peak_position_error"] if r["lm"] else None),
        "wiener_resolved": np.array([r["w_resolved"] for r in rows], dtype=bool),
        "lm_resolved": np.array([r["lm_resolved"] for r in rows], dtype=bool),
        # "LM ran here", NOT "LM produced a finite Cv" — a successful LM run
        # that finds <2 peaks is still a run, and the crossing needs it.
        "lm_ran": ok,
        "kernel_sigma_ns": np.array(kw["sigma"]),
        "kernel_fwhm_ns": np.array(kw["fwhm"]),
        "kernel_centroid_ns": np.array(kw["t_centroid"]),
        "dt_min_wiener_ns": np.array(dtmin_w),
        "dt_min_lm_ns": np.array(dtmin_lm),
        "dt_min_truth_ns": np.array(dtmin_t),
        "C_MIN": np.array(D.C_MIN),
        "NRMSE_MAX": np.array(D.NRMSE_MAX),
    }
    path = D.C.save_npz(D.SUBDIR, D.npz_name("time_resolution", quick), **arrays)

    meta = {
        "is_quick_run": bool(quick),
        "param_hash": D.param_hash("time_resolution", quick),
        "generated": datetime.now().isoformat(timespec="seconds"),
        "config": D.frozen_config("time_resolution"),
        "separations_used_ns": seps.tolist(),
        "lm_separations_used_ns": np.asarray(lm_separations, float).tolist(),
        "kernel": {k: (None if isinstance(v, float) and not np.isfinite(v) else v)
                   for k, v in kw.items()},
        "dt_min": {
            "wiener_ns": None if not np.isfinite(dtmin_w) else dtmin_w,
            "wiener_status": st_w,
            "lm_ns": None if not np.isfinite(dtmin_lm) else dtmin_lm,
            "lm_status": st_lm,
            "truth_ceiling_ns": None if not np.isfinite(dtmin_t) else dtmin_t,
            "truth_ceiling_status": st_t,
            "criterion": (f"two peaks AND Cv>={D.C_MIN} AND "
                          f"NRMSE<={D.NRMSE_MAX} AND pos_err<="
                          f"{D.PEAK_POSITION_TOL_NS} ns"),
        },
        "lm_status_per_separation": scan["lm_status"],
        "runtime_s": scan["runtime"],
        "notes": [
            "Peak WIDTH is fixed at sigma=5 ns; ONLY the spacing is scanned "
            "(FluxSignal type=8 custom array; type=5 couples the two).",
            "Two sigma=5 ns Gaussians have no valley at all for spacing "
            "<= 2*sigma = 10 ns, so the truth itself is single-peaked there; "
            "dt_min.truth_ceiling_ns is that geometric floor under the same "
            "Cv criterion, and no sensor can beat it.",
            "The 0.5 ns AWG sample interval is NOT the time resolution.",
            "Truth t_flux is linspace(0,200,400) (dt=0.50125) while the "
            "reconstruction axis is arange dt=0.5, a <=0.5 ns skew at the "
            "window edge; peak-position errors include it and stay well "
            f"inside the {D.PEAK_POSITION_TOL_NS} ns tolerance.",
            "LM values are DISCRETE verification points on a frozen subset; "
            "NaN elsewhere means not run, never Wiener back-fill.",
        ],
        "environment": D.environment(),
    }
    D.write_metadata_section("time_resolution", meta)
    return path, meta


def main():
    quick = "--quick" in sys.argv
    with_lm = "--no-lm" not in sys.argv
    force = "--force" in sys.argv

    seps = QUICK_SEPARATIONS_NS if quick else D.SEPARATIONS_NS
    lm_seps = QUICK_LM_SEPARATIONS_NS if quick else D.LM_SEPARATIONS_NS

    target = D.out_path(D.npz_name("time_resolution", quick) + ".npz")
    if os.path.exists(target) and not force:
        raise SystemExit(
            f"{os.path.basename(target)} already exists; pass --force to overwrite "
            f"(quick and formal caches use separate filenames and never collide)")

    kind = "QUICK SMOKE (not a formal result)" if quick else "FORMAL"
    print(f"[D2-A time resolution] {kind}: {len(seps)} spacings, "
          f"LM on {len(lm_seps) if with_lm else 0}")
    print(f"  frozen: C_MIN={D.C_MIN}  NRMSE_MAX={D.NRMSE_MAX}  "
          f"pos_tol={D.PEAK_POSITION_TOL_NS} ns  sigma_peak={D.PEAK_SIGMA_NS} ns")

    t0 = time.time()
    scan = run_scan(seps, lm_seps, with_lm=with_lm)
    path, meta = save(scan, seps, lm_seps, quick=quick)

    kw = scan["kw"]
    print(f"\n  kernel: sigma_k={kw['sigma']:.3f} ns  FWHM="
          f"{kw['fwhm']:.3f} ns ({kw['fwhm_status']})")
    dm = meta["dt_min"]
    print(f"  dt_min Wiener = {dm['wiener_ns']} ns ({dm['wiener_status']})")
    print(f"  dt_min LM     = {dm['lm_ns']} ns ({dm['lm_status']})")
    print(f"  truth ceiling = {dm['truth_ceiling_ns']} ns "
          f"({dm['truth_ceiling_status']}) <- signal geometry, not sensor")
    print(f"  total {time.time() - t0:.1f} s -> {path}")


if __name__ == "__main__":
    main()


