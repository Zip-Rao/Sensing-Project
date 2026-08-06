#!/usr/bin/env python3
"""D2 report figure — time resolution and sensitivity (2x2).

    (a) normalised |k(t)|^2 with sigma_k (and FWHM when |k| is single-peaked)
    (b) valley contrast Cv and NRMSE vs peak spacing, with the frozen criteria
        and Delta_t_min  (stacked sub-panels sharing the x axis — no twin y
        axis, which would let one curve hide the criterion line)
    (c) SNR vs amplitude with SNR_THRESHOLD and A_min
    (d) A_min vs N_shot (eta_phi is NOT plotted: no duty cycle is defined)

READ-ONLY: loads the cached npz/json and never re-runs quantum dynamics, so
restyling is instant and cannot silently change the numbers. Verifies each
file's param_hash against the current frozen config and warns on mismatch.

sqc results only -- no frozen-src overlay (D2 has no src counterpart).
The report figure shows truth/threshold in black or grey and Wiener in red.
LM spot checks remain in the caches and metrics, but are omitted from the main
figure because their shot budget and sampling density are not like-for-like.
Use --show-lm to create a separate diagnostic figure with the LM overlays.

Usage:
    python result_sqc/reconstruction/D2_resolution_sensitivity/plot_resolution_sensitivity_sqc.py
    python result_sqc/reconstruction/D2_resolution_sensitivity/plot_resolution_sensitivity_sqc.py --quick
    python result_sqc/reconstruction/D2_resolution_sensitivity/plot_resolution_sensitivity_sqc.py --show-lm
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _d2_common as D  # noqa: E402

C_TRUTH, C_WIENER, C_LM, C_CRIT = "black", "#e41a1c", "#377eb8", "#555555"
FIG_NAME = "resolution_sensitivity_sqc"


def load(quick=False):
    """Load both npz caches + metadata; check param hashes."""
    out = {}
    for stem, section in (("time_resolution", "time_resolution"),
                          ("sensitivity", "sensitivity")):
        path = D.out_path(D.npz_name(stem, quick) + ".npz")
        if not os.path.exists(path):
            raise SystemExit(
                f"missing {os.path.basename(path)} — run "
                f"generate_{stem}_sqc.py"
                + (" --quick" if quick else ""))
        out[section] = dict(np.load(path))
    meta = D.read_metadata()
    for section in ("time_resolution", "sensitivity"):
        m = meta.get(section, {})
        want = D.param_hash(section, quick)
        got = m.get("param_hash")
        if got and got != want:
            print(f"  WARNING: {section} param_hash {got} != current config "
                  f"{want} — the cache was made with different frozen "
                  f"parameters; regenerate before trusting this figure.")
        if bool(m.get("is_quick_run")) != bool(quick):
            print(f"  WARNING: {section} cache is_quick_run="
                  f"{m.get('is_quick_run')} but --quick={quick}")
    return out, meta


def panel_a(ax, tr, meta):
    """(a) normalised |k(t)|^2, sigma_k band, FWHM when |k| is single-peaked."""
    t = tr["kernel_time"]
    k = tr["kernel"]
    k2 = np.abs(k) ** 2
    k2n = k2 / k2.max()
    sigma = float(tr["kernel_sigma_ns"])
    fwhm = float(tr["kernel_fwhm_ns"])
    tc = float(tr["kernel_centroid_ns"])

    ax.plot(t, k2n, color=C_TRUTH, lw=1.8)
    ax.fill_between(t, 0, k2n, color=C_TRUTH, alpha=0.08)
    ax.axvspan(tc - sigma, tc + sigma, color=C_WIENER, alpha=0.13,
               label=rf"$t_c\pm\sigma_k$, $\sigma_k$={sigma:.2f} ns")
    ax.axvline(tc, color=C_CRIT, ls=":", lw=1.0)

    kn = np.abs(k) / np.abs(k).max()
    ax.plot(t, kn, color=C_WIENER, lw=1.2, ls="--", alpha=0.8,
            label=r"$|k(t)|$ (normalised)")
    status = (meta.get("time_resolution", {}).get("kernel", {})
              .get("fwhm_status", ""))
    if np.isfinite(fwhm):
        ax.axhline(0.5, color=C_CRIT, ls="-.", lw=0.8, alpha=0.7)
        ax.annotate(f"FWHM({'|k|'})={fwhm:.2f} ns", xy=(0.03, 0.60),
                    xycoords="axes fraction", fontsize=8, color=C_CRIT)
    else:
        ax.annotate(f"FWHM undefined\n({status})", xy=(0.03, 0.58),
                    xycoords="axes fraction", fontsize=7, color=C_CRIT)

    ax.set_xlabel("Kernel time $t$ (ns)", fontsize=10)
    ax.set_ylabel(r"$|k|^2$, $|k|$ (normalised)", fontsize=10)
    ax.set_title("(a) Response-kernel effective width", fontsize=10.5)
    ax.legend(fontsize=7.5, loc="upper right")
    ax.grid(True, ls=":", alpha=0.4, lw=0.5)


def panel_b(ax_cv, ax_nr, tr, meta, show_lm=True):
    """(b) Cv (top) and NRMSE (bottom) vs spacing, sharing the x axis.

    Stacked rather than twin-y on purpose: a second y axis would let the
    NRMSE curve overlap and visually bury the Cv criterion line.

    show_lm=False drops the LM series entirely (cleaner, but loses the
    time-vs-amplitude trade-off and the guide's discrete-verification-point
    requirement -- see module docstring).
    """
    sep = tr["separations_ns"]
    cv_w = tr["wiener_valley_contrast"]
    cv_t = tr["truth_valley_contrast"]
    cv_lm = tr["lm_valley_contrast"]
    # Recomputed from the arrays, not read from metadata, so the figure is
    # self-consistent with the saved data even if the cache predates a fix.
    cr = D.dt_min_from_arrays(tr)
    dm = {"wiener_ns": cr["wiener"][0], "lm_ns": cr["lm"][0],
          "truth_ceiling_ns": cr["truth"][0]}
    dm = {k: (v if np.isfinite(v) else None) for k, v in dm.items()}

    # truth ceiling: what the noiseless input itself can show
    ax_cv.plot(sep, cv_t, color=C_TRUTH, ls="--", lw=1.4, marker="s", ms=3.5,
               label="Truth (input geometry)")
    ax_cv.plot(sep, cv_w, color=C_WIENER, ls="-", lw=1.8, marker="o", ms=4.5,
               label="Wiener")
    ok = np.isfinite(cv_lm)
    if show_lm and ok.any():
        ax_cv.plot(sep[ok], cv_lm[ok], color=C_LM, ls="none", marker="D",
                   ms=5.5, mfc="none", mew=1.6,
                   label="LM (discrete points)")
    # unresolved: fewer than two peaks detected -> Cv undefined, marked not zero
    # <2 peaks -> Cv undefined. Drawn BELOW zero so it cannot be misread as
    # "Cv measured and equal to 0".
    miss = ~np.isfinite(cv_w)
    if miss.any():
        ax_cv.plot(sep[miss], np.full(miss.sum(), -0.075), color=C_WIENER,
                   ls="none", marker="x", ms=6, mew=1.4,
                   label=r"Wiener: <2 peaks ($C_{\rm v}$ undefined)")
    miss_lm = (~np.isfinite(cv_lm)) & tr.get(
        "lm_ran", np.zeros(len(sep), bool)).astype(bool)
    if show_lm and miss_lm.any():
        ax_cv.plot(sep[miss_lm], np.full(miss_lm.sum(), -0.14), color=C_LM,
                   ls="none", marker="x", ms=6, mew=1.4,
                   label=r"LM: <2 peaks ($C_{\rm v}$ undefined)")

    ax_cv.axhline(D.C_MIN, color=C_CRIT, ls="-.", lw=1.2,
                  label=rf"$C_{{\min}}$={D.C_MIN}")
    # The three Delta_t_min lines fall within ~4 ns of each other, so the
    # labels are stacked at staggered heights instead of all sitting at one
    # level, where they would overprint each other and the C_min label.
    threshold_lines = [
        (("truth_ceiling_ns", C_TRUTH, "truth"), 0.97),
        (("wiener_ns", C_WIENER, "Wiener"), 0.78),
    ]
    if show_lm:
        threshold_lines.append((("lm_ns", C_LM, "LM"), 0.59))
    for (key, col, lab), yfrac in threshold_lines:
        v = dm.get(key)
        if v:
            ax_cv.axvline(v, color=col, ls=":", lw=1.3, alpha=0.75)
            ax_cv.annotate(rf"$\Delta t_{{\min}}^{{\rm {lab}}}$={v:.1f} ns",
                           xy=(v, yfrac), xycoords=("data", "axes fraction"),
                           xytext=(6, 0), textcoords="offset points",
                           fontsize=7, color=col, ha="left", va="top",
                           bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                     ec="none", alpha=0.75))
    ax_cv.set_ylabel(r"Valley contrast $C_{\rm v}$", fontsize=9.5)
    ax_cv.set_title("(b) Double-peak resolvability vs spacing", fontsize=10.5)
    ax_cv.legend(fontsize=6.2, loc="center right", ncol=1, framealpha=0.9)
    ax_cv.grid(True, ls=":", alpha=0.4, lw=0.5)
    ax_cv.tick_params(labelbottom=False)
    ax_cv.set_ylim(-0.24, 1.42)   # headroom for the stacked dt_min labels

    ax_nr.plot(sep, tr["wiener_nrmse"], color=C_WIENER, ls="-", lw=1.8,
               marker="o", ms=4.5, label="Wiener")
    nl = tr["lm_nrmse"]
    okn = np.isfinite(nl)
    if show_lm and okn.any():
        ax_nr.plot(sep[okn], nl[okn], color=C_LM, ls="none", marker="D",
                   ms=5.5, mfc="none", mew=1.6, label="LM")
    ax_nr.axhline(D.NRMSE_MAX, color=C_CRIT, ls="-.", lw=1.2,
                  label=rf"NRMSE$_{{\max}}$={D.NRMSE_MAX}")
    ax_nr.set_xlabel(r"Peak spacing $\Delta t$ (ns)", fontsize=10)
    ax_nr.set_ylabel("NRMSE", fontsize=9.5)
    ax_nr.legend(fontsize=6.8, loc="upper left")
    ax_nr.grid(True, ls=":", alpha=0.4, lw=0.5)
    ax_nr.set_ylim(0, max(D.NRMSE_MAX * 1.25,
                          np.nanmax(tr["wiener_nrmse"]) * 1.3))


def panel_c(ax, sn, meta, show_lm=False):
    """(c) SNR vs amplitude, error bars = SEM over frozen seeds."""
    amps = sn["amplitudes"]
    shots = sn["n_shots_values"]
    snr = sn["snr_wiener"]
    sem = sn["snr_sem_wiener"]
    rows = meta.get("sensitivity", {}).get("a_min", {}).get("per_n_shot", [])
    alphas = np.linspace(0.42, 1.0, len(shots))

    for si, N in enumerate(shots):
        ax.errorbar(amps, snr[si], yerr=sem[si], color=C_WIENER,
                    alpha=alphas[si], lw=1.6, marker="o", ms=4.2,
                    capsize=2.5, elinewidth=0.9,
                    label=rf"Wiener, $N_{{\rm shot}}$={int(N):,}")
    snr_lm = sn["snr_lm"]
    ok = np.isfinite(snr_lm)
    if show_lm and ok.any():
        si, ai = np.where(ok)
        ax.plot(amps[ai], snr_lm[si, ai], color=C_LM, ls="none", marker="D",
                ms=6, mfc="none", mew=1.6, label="LM (discrete points)")

    ax.axhline(D.SNR_THRESHOLD, color=C_CRIT, ls="-.", lw=1.2,
               label=rf"SNR$_{{\rm th}}$={D.SNR_THRESHOLD:g}")
    # A_min labels sit on the threshold line with a white backing box; placing
    # them along the bottom axis put them straight through the SNR curves.
    for si, row in enumerate(rows):
        v = row.get("a_min")
        if v:
            ax.axvline(v, color=C_WIENER, ls=":", lw=1.1,
                       alpha=alphas[si] if si < len(alphas) else 1.0)
            ax.annotate(rf"{v:.2e}", xy=(v, D.SNR_THRESHOLD),
                        xytext=(-2, -14 - 11 * si), textcoords="offset points",
                        fontsize=6.4, color=C_WIENER, ha="right",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                  ec="none", alpha=0.8))
    # (no "A_min" hint text here — the labelled dotted verticals already say
    # it, and the hint collided with the legend.)
    if show_lm and np.any(np.isfinite(snr_lm)):
        ax.annotate("LM: few seeds, signal arm only\n"
                    "(half the shots) — not like-for-like",
                    xy=(0.97, 0.06), xycoords="axes fraction", fontsize=6.3,
                    color=C_LM, ha="right", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.25", fc="#f5f5f5",
                              ec=C_LM, lw=0.6, alpha=0.9))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Signal amplitude $A$ ($\Phi_0$, peak)", fontsize=10)
    ax.set_ylabel(r"SNR $=\mathbb{E}[\hat A]/\sigma_{\hat A,0}$", fontsize=10)
    ax.set_title("(c) Sensitivity under binomial projection noise",
                 fontsize=10.5)
    ax.legend(fontsize=6.8, loc="upper left")
    ax.grid(True, which="both", ls=":", alpha=0.4, lw=0.5)


def panel_d(ax, sn, meta):
    """(d) A_min vs N_shot, with a 1/sqrt(N) guide. eta_phi is NOT plotted."""
    shots = np.asarray(sn["n_shots_values"], dtype=float)
    rows = meta.get("sensitivity", {}).get("a_min", {}).get("per_n_shot", [])
    a_min = np.array([r.get("a_min") if r.get("a_min") else np.nan
                      for r in rows], dtype=float)
    sigma0 = np.asarray(sn["sigma_Ahat_zero_wiener"], dtype=float)

    ax.plot(shots, a_min, color=C_WIENER, lw=1.8, marker="o", ms=7,
            label=rf"$A_{{\min}}$ (SNR$\geq${D.SNR_THRESHOLD:g}), Wiener")
    ax.plot(shots, sigma0, color=C_CRIT, lw=1.3, ls="--", marker="s", ms=5,
            mfc="none", label=r"noise floor $\sigma_{\hat A,0}$")

    ok = np.isfinite(a_min)
    if ok.sum() >= 2:
        i = int(np.where(ok)[0][0])
        guide = a_min[i] * np.sqrt(shots[i] / shots)
        ax.plot(shots, guide, color=C_TRUTH, ls=":", lw=1.2, alpha=0.75,
                label=r"$\propto 1/\sqrt{N_{\rm shot}}$")
        p = np.polyfit(np.log10(shots[ok]), np.log10(a_min[ok]), 1)
        ax.annotate(rf"fit slope {p[0]:.2f}  (ideal $-0.5$)",
                    xy=(0.97, 0.63), xycoords="axes fraction", fontsize=7.5,
                    color=C_WIENER, ha="right")

    ax.annotate(r"$\eta_\Phi$ not computable:" "\n"
                r"no readout/reset time defined" "\n"
                r"$\Rightarrow$ $A_{\min}$ only",
                xy=(0.97, 0.97), xycoords="axes fraction", fontsize=7,
                color=C_CRIT, ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.35", fc="#f5f5f5",
                          ec=C_CRIT, lw=0.7))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Shots per point $N_{\rm shot}$", fontsize=10)
    ax.set_ylabel(r"Amplitude ($\Phi_0$, peak)", fontsize=10)
    ax.set_title(r"(d) Detection floor vs shot budget", fontsize=10.5)
    ax.legend(fontsize=6.8, loc="lower left")
    ax.grid(True, which="both", ls=":", alpha=0.4, lw=0.5)


def write_metrics(tr, sn, meta, quick=False):
    """Metrics TXT recomputed FROM THE NPZ, so every figure number is traceable."""
    L = []
    add = L.append
    add("D2 — time resolution and sensitivity (sqc)")
    add("=" * 62)
    if quick:
        add("*** QUICK SMOKE RUN — NOT A FORMAL RESULT ***")
    add(f"generated from: {D.npz_name('time_resolution', quick)}.npz, "
        f"{D.npz_name('sensitivity', quick)}.npz")
    add(f"qubit: EC={D.C.EC:.6f} EJ={D.C.EJ:.6f} T1={D.C.T1:g} ns "
        f"T2={D.C.T2:g} ns  work point Phi={D.C.OPTIMAL_FLUX:.6f}")
    add("")
    add("A. TIME RESOLUTION")
    add("-" * 62)
    kmeta = meta.get("time_resolution", {}).get("kernel", {})
    add(f"  kernel sigma_k        = {float(tr['kernel_sigma_ns']):.4f} ns")
    fw = float(tr["kernel_fwhm_ns"])
    add(f"  kernel FWHM(|k|)      = "
        + (f"{fw:.4f} ns" if np.isfinite(fw) else "undefined")
        + f"   [{kmeta.get('fwhm_status', '')}]")
    add(f"  kernel centroid       = {float(tr['kernel_centroid_ns']):.4f} ns")
    add(f"  kernel sign changes   = {kmeta.get('n_sign_changes')}")
    add(f"  fixed peak sigma      = {D.PEAK_SIGMA_NS} ns (never scanned)")
    add(f"  criteria: two peaks AND Cv>={D.C_MIN} AND NRMSE<={D.NRMSE_MAX} "
        f"AND pos_err<={D.PEAK_POSITION_TOL_NS} ns")
    add("  (Delta_t_min below is RECOMPUTED from the npz arrays, so it is "
        "traceable\n   and independent of the metadata JSON)")
    cr = D.dt_min_from_arrays(tr)
    for key, lab in (("truth", "truth ceiling (input geometry)"),
                     ("wiener", "Wiener  Delta_t_min"),
                     ("lm", "LM      Delta_t_min")):
        v, st = cr[key]
        add(f"  {lab:<32s}= "
            + (f"{v:.3f} ns" if np.isfinite(v) else "not bracketed")
            + f"   ({st})")
    add("  NOTE: a NaN Cv means <2 peaks were detected, i.e. DEFINITIVELY "
        "unresolved;\n        it is treated as Cv=0 for the crossing, not as "
        "missing data.")
    add("")
    add(f"  {'dt(ns)':>7s} {'Cv_truth':>9s} {'Cv_W':>8s} {'NRMSE_W':>8s} "
        f"{'res_W':>6s} {'Cv_LM':>8s} {'NRMSE_LM':>9s} {'res_LM':>7s}")
    for i, s in enumerate(tr["separations_ns"]):
        f = lambda v: f"{v:8.4f}" if np.isfinite(v) else "     nan"  # noqa: E731
        add(f"  {s:7.1f} {tr['truth_valley_contrast'][i]:9.4f} "
            f"{f(tr['wiener_valley_contrast'][i])} "
            f"{tr['wiener_nrmse'][i]:8.4f} "
            f"{str(bool(tr['wiener_resolved'][i])):>6s} "
            f"{f(tr['lm_valley_contrast'][i])} "
            f"{f(tr['lm_nrmse'][i]):>9s} "
            f"{str(bool(tr['lm_resolved'][i])):>7s}")
    add("")
    add("B. SENSITIVITY (binomial projection noise)")
    add("-" * 62)
    smeta = meta.get("sensitivity", {})
    nm = smeta.get("noise_model", {})
    add(f"  noise model     : {nm.get('formula')}")
    add(f"  arms sampled independently: {nm.get('arms_sampled_independently')}")
    add(f"  input flux noise_level    : {nm.get('input_flux_noise_level')} "
        f"(separate from projection noise)")
    add(f"  amplitude estimator       : template projection "
        f"<phi_hat,s>/<s,s>")
    add(f"  uncertainty               : "
        f"{smeta.get('uncertainty', {}).get('primary')}")
    add(f"  SNR threshold             : {D.SNR_THRESHOLD:g} "
        f"(secondary {D.SNR_THRESHOLD_SECONDARY:g})")
    add("")
    for si, N in enumerate(sn["n_shots_values"]):
        s0 = sn["sigma_Ahat_zero_wiener"][si]
        m0 = sn["mean_Ahat_zero_wiener"][si]
        add(f"  N_shot={int(N):>8d}: sigma_Ahat_0={s0:.4e}  "
            f"mean_Ahat_0={m0:+.4e} (should be ~0)")
        for ai, A in enumerate(sn["amplitudes"]):
            add(f"      A={A:.3e}  Ahat={sn['mean_wiener'][si, ai]:.4e} "
                f"+/-{sn['sem_wiener'][si, ai]:.2e}  "
                f"SNR={sn['snr_wiener'][si, ai]:8.3f} "
                f"[{sn['snr_ci_lo_wiener'][si, ai]:.2f},"
                f"{sn['snr_ci_hi_wiener'][si, ai]:.2f}]")
    add("")
    for row in smeta.get("a_min", {}).get("per_n_shot", []):
        v = row.get("a_min")
        add(f"  A_min(N={int(row['n_shot']):>8d}) = "
            + (f"{v:.4e} Phi_0" if v else "not bracketed")
            + f"   ({row.get('status')})")
    add("")
    ep = smeta.get("eta_phi", {})
    add(f"  eta_phi = NaN  — {ep.get('reason', '')}")
    add(f"  missing inputs: {', '.join(ep.get('missing_inputs', []))}")
    add("")
    add(f"  LM: {smeta.get('lm', {}).get('status_note')}")
    add(f"      {smeta.get('lm', {}).get('shot_budget_caveat')}")
    add("")
    add("CAVEATS")
    add("-" * 62)
    for n in (meta.get("time_resolution", {}).get("notes", [])
              + smeta.get("notes", [])):
        add(f"  - {n}")

    path = D.out_path(D.METRICS_NAME if not quick
                      else D.METRICS_NAME.replace(".txt", "_quick.txt"))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return path


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    plt.style.use("seaborn-v0_8-whitegrid")

    quick = "--quick" in sys.argv
    show_lm = "--show-lm" in sys.argv
    data, meta = load(quick)
    tr, sn = data["time_resolution"], data["sensitivity"]

    fig = plt.figure(figsize=(11.5, 8.2))
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.26)
    ax_a = fig.add_subplot(gs[0, 0])
    # (b) is two stacked sub-panels sharing x — never a twin y axis
    sub = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 1],
                                  height_ratios=[2.1, 1.0], hspace=0.08)
    ax_b1 = fig.add_subplot(sub[0])
    ax_b2 = fig.add_subplot(sub[1], sharex=ax_b1)
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_a(ax_a, tr, meta)
    panel_b(ax_b1, ax_b2, tr, meta, show_lm=show_lm)
    panel_c(ax_c, sn, meta, show_lm=show_lm)
    panel_d(ax_d, sn, meta)

    title = ("Time resolution and sensitivity of transient flux sensing (sqc)")
    if quick:
        title += "  —  QUICK SMOKE, NOT A FORMAL RESULT"
    fig.suptitle(title, fontsize=12.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    stem = FIG_NAME + ("_quick" if quick else "")
    if show_lm:
        stem += "_with_lm"
    png = D.out_path(stem + ".png")
    fig.savefig(png, dpi=300)
    fig.savefig(D.out_path(stem + ".pdf"))   # vector
    plt.close(fig)

    mp = write_metrics(tr, sn, meta, quick=quick)
    print(f"figure  -> {png}")
    print(f"        -> {png.replace('.png', '.pdf')}")
    print(f"metrics -> {mp}")


if __name__ == "__main__":
    main()


