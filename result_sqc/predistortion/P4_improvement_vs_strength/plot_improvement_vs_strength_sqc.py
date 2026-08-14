#!/usr/bin/env python3
"""P4 ⭐ — Predistortion improvement vs distortion strength (sqc).

sqc source:
    sqc.workflows.PredistortionValidationWorkflow, swept over the injected
    distortion strength A. Each run() returns metrics {rmse_uncorrected,
    rmse_corrected, improvement_factor, settling_*_ns} plus the AWG waveform
    that had to be played.

Physics story (⭐ paper role): quantify WHERE predistortion helps, and where it
stops helping. THIS SCRIPT'S ORIGINAL PREMISE WAS WRONG and the figure is built
to show why:
    Originally assumed: "uncorrected rmse rises but corrected rmse stays flat
    → improvement_factor grows [with strength]".
    Measured: uncorrected rmse ∝ A^1.0000 (exact), corrected rmse ∝ A^2.06,
    hence improvement ∝ A^-1.06 — it DECREASES with strength. Predistortion
    helps most for WEAK distortion (×783 at A=0.01, ×12.7 at A=0.50).
Two scenarios are swept, for the reason established in P2/P3:
    "matched"  single_exp truth + single_exp fit → the inverse cancels
               algebraically, so improvement is pure floating-point noise
               (scatters 4e8..2e14 with NO trend). Sweeping only this, as
               originally planned, would have produced a meaningless figure.
    "mismatch" 2-component truth + 3-component fit → a real O(A²) residual
               and a clean 1/A improvement law.

Produces:
    result_sqc/predistortion/P4_improvement_vs_strength/improvement_vs_strength_sqc.png / .pdf
        (a) rmse vs A, log-log, with fitted power laws
        (b) improvement factor vs A: 1/A law vs floating-point scatter
        (c) settling time vs A: the 1e-3 spec breaks between A=0.05 and 0.10
        (d) AWG dynamic range demanded by the correction
    result_sqc/predistortion/P4_improvement_vs_strength/improvement_vs_strength_sqc.npz

Cost: medium — 2 scenarios × 12 strengths, analytical measurement path.
Cached to npz; pass --recompute to force a rebuild.

Usage:
    python result_sqc/predistortion/P4_improvement_vs_strength/plot_improvement_vs_strength_sqc.py
    python result_sqc/predistortion/P4_improvement_vs_strength/plot_improvement_vs_strength_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

from sqc.config import CONFIG  # noqa: E402

SUBDIR = "predistortion/P4_improvement_vs_strength"
_CACHE = "improvement_vs_strength_sqc"

# --- target waveform (R9: time axis derives from CONFIG.awg.dt) ----------
# Same flat-top target as P2 so the two figures are directly comparable.
DT = CONFIG.awg.dt
T_END = 500.0
FLAT_START, FLAT_STOP = 50.0, 250.0

# Injected total tail amplitude A. Extended past the skeleton's 0.10 ceiling to
# expose the breakdown: the settling spec fails above ~0.05 and the AWG range
# demand blows up. NOTE SingleExponentialDistortion.design_inverse warns that
# the Rol 2020 verified range is |A| <= 0.1, so A >= 0.15 is extrapolation —
# panels shade it.
STRENGTHS = np.array([0.005, 0.01, 0.02, 0.04, 0.06, 0.08,
                      0.10, 0.15, 0.20, 0.30, 0.40, 0.50])
ROL_VALID_MAX = 0.10
FIXED_TAU_NS = 100.0                 # matched scenario tail time constant
MISMATCH_TAUS = (80.0, 400.0)        # mismatch scenario poles
MISMATCH_SPLIT = (2.0 / 3.0, 1.0 / 3.0)   # how A is divided between the poles

SETTLE_TOL = 1e-3                    # the workflow's own settling tolerance

SCENARIOS = ("matched", "mismatch")
LABELS = {
    "matched": "matched: single_exp truth + single_exp fit",
    "mismatch": r"mismatch: 2-pole truth (80/400 ns) + 3-component fit",
}
STYLE = {
    "matched":  {"color": "#377eb8", "marker": "o", "ls": "-"},
    "mismatch": {"color": "#ff7f00", "marker": "s", "ls": "-"},
}


def make_target():
    """Desired on-chip flat-top waveform (identical to P2's)."""
    from sqc.control.waveform import Waveform

    t = np.arange(0.0, T_END, DT)
    return Waveform(t_list=t,
                    samples=np.where((t > FLAT_START) & (t < FLAT_STOP),
                                     1.0, 0.0))


def true_distortion(scenario, amp):
    """Injected ground truth at total tail amplitude `amp`.

    Both scenarios carry the SAME total amplitude so the uncorrected error is
    the same at each A (verified: rmse_uncorrected agrees to ~0.4%). Only the
    pole structure differs, which isolates model mismatch as the single cause
    of the corrected-error difference.
    """
    from sqc.hardware.distortion import (
        SingleExponentialDistortion, MultiExponentialDistortion,
    )
    if scenario == "matched":
        return SingleExponentialDistortion(amplitude=float(amp),
                                           tau=FIXED_TAU_NS)
    amps = np.array([amp * MISMATCH_SPLIT[0], amp * MISMATCH_SPLIT[1]])
    return MultiExponentialDistortion(amplitudes=amps,
                                      taus=np.array(MISMATCH_TAUS))


def compute(recompute=False):
    """Sweep both scenarios over STRENGTHS. Cached to npz."""
    from sqc.workflows import PredistortionValidationWorkflow

    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    target = make_target()
    out = {"strengths": STRENGTHS}
    for s in SCENARIOS:
        rows = {k: [] for k in ("unc", "cor", "impr", "set_unc", "set_cor",
                                "awg_lo", "awg_hi")}
        for amp in STRENGTHS:
            r = PredistortionValidationWorkflow(
                target_waveform=target,
                true_distortion=true_distortion(s, amp)).run()
            m = r["metrics"]
            awg = r["awg_predistorted"].samples
            rows["unc"].append(m["rmse_uncorrected"])
            rows["cor"].append(m["rmse_corrected"])
            rows["impr"].append(m["improvement_factor"])
            rows["set_unc"].append(m["settling_uncorrected_ns"])
            rows["set_cor"].append(m["settling_corrected_ns"])
            rows["awg_lo"].append(float(awg.min()))
            rows["awg_hi"].append(float(awg.max()))
        for k, v in rows.items():
            out[f"{s}_{k}"] = np.array(v, dtype=float)
        # power-law exponents from a log-log fit
        for k in ("unc", "cor"):
            sl = np.polyfit(np.log(STRENGTHS), np.log(out[f"{s}_{k}"]), 1)[0]
            out[f"{s}_{k}_slope"] = np.array([sl])

    C.save_npz(SUBDIR, _CACHE, **out)
    return out


def _shade_extrapolation(ax):
    """Grey out A > ROL_VALID_MAX (outside the Rol 2020 verified range)."""
    ax.axvspan(ROL_VALID_MAX, STRENGTHS[-1] * 1.12, color="gray", alpha=0.11,
               zorder=0)
    ax.axvline(ROL_VALID_MAX, color="gray", ls="--", lw=0.9, alpha=0.8)


def _panel_rmse(ax, d):
    """(a) rmse vs A, log-log, with the measured power laws annotated."""
    A = d["strengths"]
    _shade_extrapolation(ax)
    ax.loglog(A, d["mismatch_unc"], "k^-", ms=5, lw=1.4,
              label=r"uncorrected (both scenarios), $\propto A^{%.4f}$"
                    % d["mismatch_unc_slope"][0])
    for s in SCENARIOS:
        st = STYLE[s]
        # Only the mismatch residual obeys a power law. The matched "corrected"
        # values are roundoff (8e-16..2e-12); a log-log slope through them is an
        # artefact (it comes out A^-0.42) and must NOT be presented as a law.
        tag = (rf"$\propto A^{{{d[f'{s}_cor_slope'][0]:.2f}}}$"
               if s == "mismatch" else "floating-point noise, no trend")
        ax.loglog(A, d[f"{s}_cor"], st["marker"] + st["ls"], color=st["color"],
                  ms=5, lw=1.4, label=f"corrected, {s}: {tag}")
    ax.axhline(SETTLE_TOL, color="gray", lw=0.9, ls=":", alpha=0.9)
    # placed low-left: the uncorrected line passes just above 1e-3 at A[0]
    ax.text(A[0] * 1.05, SETTLE_TOL * 0.28, r"$10^{-3}$ tolerance", fontsize=8,
            color="gray")
    ax.set_xlabel("Injected tail amplitude $A$", fontsize=11)
    ax.set_ylabel("rmse vs target", fontsize=11)
    ax.set_title("(a) Error scaling: uncorrected is linear in $A$,\n"
                 "the mismatch residual is quadratic", fontsize=10.5)
    # center right: the gap between the mismatch residual and the roundoff
    # floor is empty there, and lower-right would sit on the matched curve
    ax.legend(fontsize=7.5, loc="center right", framealpha=0.92)
    ax.grid(True, which="both", ls=":", alpha=0.4)


def _panel_improvement(ax, d):
    """(b) the star panel: 1/A law vs floating-point scatter."""
    A = d["strengths"]
    _shade_extrapolation(ax)
    for s in SCENARIOS:
        st = STYLE[s]
        ax.loglog(A, d[f"{s}_impr"], st["marker"] + st["ls"], color=st["color"],
                  ms=6, lw=1.4, label=LABELS[s])
    # 1/A reference anchored on the mismatch curve
    ref = d["mismatch_impr"][0] * (A / A[0]) ** -1.0
    ax.loglog(A, ref, "k:", lw=1.2, label=r"$\propto 1/A$ reference")
    ax.set_xlabel("Injected tail amplitude $A$", fontsize=11)
    ax.set_ylabel("improvement factor  rmse$_{unc}$ / rmse$_{cor}$", fontsize=11)
    ax.set_title("(b) Improvement DECREASES with strength ($\\propto 1/A$).\n"
                 "Matched is floating-point noise, not a trend", fontsize=10.5)
    # center left: the decades between the roundoff-driven matched curve and the
    # real mismatch law are empty, so the legend cannot cover either
    ax.legend(fontsize=7.5, loc="center left", framealpha=0.92)
    ax.grid(True, which="both", ls=":", alpha=0.4)


def _panel_settling(ax, d):
    """(c) settling time vs A — where the 1e-3 spec stops being met."""
    A = d["strengths"]
    _shade_extrapolation(ax)
    ax.semilogx(A, d["mismatch_set_unc"], "k^-", ms=5, lw=1.4,
                label="uncorrected (never settles in window)")
    # matched stays at 0 ns across the WHOLE sweep and mismatch coincides with
    # it for A <= 0.06, so draw matched first with larger markers: its ring then
    # remains visible around the mismatch squares instead of being hidden.
    for s, ms in (("matched", 9), ("mismatch", 6)):
        st = STYLE[s]
        ax.semilogx(A, d[f"{s}_set_cor"], st["marker"] + st["ls"],
                    color=st["color"], ms=ms, lw=1.4, mfc="none" if ms == 9
                    else st["color"], label=f"corrected, {s}")
    # locate the first strength where the mismatch correction stops settling
    broke = np.nonzero(d["mismatch_set_cor"] > 0.0)[0]
    if len(broke):
        i = broke[0]
        ax.annotate(f"spec breaks\nbetween A={A[i-1]:.2f} and {A[i]:.2f}",
                    xy=(A[i], d["mismatch_set_cor"][i]),
                    xytext=(A[i] * 0.22, d["mismatch_set_cor"][i] + 130),
                    fontsize=8, color="#ff7f00",
                    arrowprops=dict(arrowstyle="->", lw=0.9, color="#ff7f00"))
    ax.set_xlabel("Injected tail amplitude $A$", fontsize=11)
    ax.set_ylabel(rf"Settling time to ${SETTLE_TOL:g}$ (ns)", fontsize=11)
    ax.set_title("(c) Settling spec: correction holds to $A\\approx0.06$,\n"
                 "then degrades toward the uncorrected line", fontsize=10.5)
    ax.legend(fontsize=7.5, loc="center left", framealpha=0.92)
    ax.grid(True, which="both", ls=":", alpha=0.4)


def _panel_headroom(ax, d):
    """(d) AWG dynamic range the correction demands."""
    A = d["strengths"]
    _shade_extrapolation(ax)
    # One axis, one quantity: this panel plots ABSOLUTE AWG amplitude only.
    # (An earlier version also drew the peak-to-peak span here, which conflates
    # a difference with an absolute level on the same axis.) The span is read
    # off as the band height and annotated at the right edge.
    for s in SCENARIOS:
        st = STYLE[s]
        ax.semilogx(A, d[f"{s}_awg_hi"], st["marker"] + "-", color=st["color"],
                    ms=5, lw=1.4, label=f"{s}: max")
        ax.semilogx(A, d[f"{s}_awg_lo"], st["marker"] + "--", color=st["color"],
                    ms=5, lw=1.4, alpha=0.75, label=f"{s}: min")
        ax.fill_between(A, d[f"{s}_awg_lo"], d[f"{s}_awg_hi"],
                        color=st["color"], alpha=0.10)
    for y, lb in ((1.0, "ideal target max = 1"), (0.0, "ideal target min = 0")):
        ax.axhline(y, color="k", lw=0.9, ls="--", alpha=0.55)
        ax.text(A[0] * 1.05, y + 0.05, lb, fontsize=7.5, color="dimgray")
    span_end = d["mismatch_awg_hi"][-1] - d["mismatch_awg_lo"][-1]
    ax.annotate(f"span {span_end:.2f}$\\times$ ideal\nat $A$={A[-1]:.2f}",
                xy=(A[-1], d["mismatch_awg_hi"][-1]),
                xytext=(A[-1] * 0.16, d["mismatch_awg_hi"][-1] - 0.42),
                fontsize=8, color="#ff7f00",
                arrowprops=dict(arrowstyle="->", lw=0.9, color="#ff7f00"))
    ax.set_xlabel("Injected tail amplitude $A$", fontsize=11)
    ax.set_ylabel("AWG amplitude the DAC must reach (a.u.)", fontsize=11)
    ax.set_title("(d) Cost in AWG dynamic range: band height is the swing\n"
                 "needed to place a unit-amplitude pulse on chip",
                 fontsize=10.5)
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.92, ncol=2)
    ax.grid(True, which="both", ls=":", alpha=0.4)


def plot(d, out):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.8))
    (ax_a, ax_b), (ax_c, ax_e) = axes
    _panel_rmse(ax_a, d)
    _panel_improvement(ax_b, d)
    _panel_settling(ax_c, d)
    _panel_headroom(ax_e, d)

    fig.suptitle("Predistortion improvement vs distortion strength (sqc)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(os.path.join(out, f"{_CACHE}.png"), dpi=300)
    fig.savefig(os.path.join(out, f"{_CACHE}.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    d = compute(recompute="--recompute" in sys.argv)

    for s in SCENARIOS:
        if s == "mismatch":
            note = (f"rmse_cor ~ A^{d[f'{s}_cor_slope'][0]:.4f}   "
                    f"=> improvement ~ A^"
                    f"{d[f'{s}_unc_slope'][0] - d[f'{s}_cor_slope'][0]:.4f}")
        else:
            note = ("rmse_cor = roundoff (no power law; the fitted slope "
                    f"A^{d[f'{s}_cor_slope'][0]:.2f} is an artefact)")
        print(f"  [{s}] rmse_unc ~ A^{d[f'{s}_unc_slope'][0]:.4f}   {note}")
        for i, a in enumerate(d["strengths"]):
            print(f"     A={a:.3f}  unc={d[f'{s}_unc'][i]:.3e} "
                  f"cor={d[f'{s}_cor'][i]:.3e} x{d[f'{s}_impr'][i]:<10.4g} "
                  f"set_cor={d[f'{s}_set_cor'][i]:6.1f} ns  "
                  f"AWG {d[f'{s}_awg_lo'][i]:+.4f}..{d[f'{s}_awg_hi'][i]:+.4f}")

    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
