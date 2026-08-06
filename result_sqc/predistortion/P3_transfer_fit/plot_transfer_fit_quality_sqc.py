#!/usr/bin/env python3
"""P3 — Transfer-function fit quality and model-order selection (sqc).

sqc source:
    sqc.calibration.waveform.WaveformCalibration(method="transfer_function")
    → measures the step response, fits a DistortionModel via curve_fit
    (waveform.py:203/229), and to_distortion_model() returns the fitted model.
    Ground truth is known, so both the step-response fit and the recovered
    (amplitude, tau) parameters can be scored.

Physics story (paper role): calibration accuracy, and specifically MODEL-ORDER
SELECTION. The informative axis is the number of fitted exponential components
K (n_exp_components), NOT the calibration window:
    The calibration window is a DEAD axis on the analytical path — an analytic
    step response is noiseless, and a sum of exponentials is uniquely
    determined by any finite arc of it, so the fit is exact even at
    t_max/tau = 0.2 (verified: rel. err < 1e-9 for t_max in 20..2000 ns).
    Fit quality is therefore governed by whether K matches the true pole count.

Produces:
    result_sqc/predistortion/P3_transfer_fit/transfer_fit_quality_sqc.png / .pdf
        (a) measured vs fitted step, 2-pole truth, K=1 (underfit) vs K=2
        (b) fit residual |fitted - measured|, log-y, K=1/2/3
        (c) fit rms vs K for three truths, true pole count marked
        (d) parameter recovery: relative error of each recovered (A_k, tau_k)
    result_sqc/predistortion/P3_transfer_fit/transfer_fit_quality_sqc.npz

Cost: LOW — analytic step response, seconds. Cached to npz; --recompute rebuilds.

Usage:
    python result_sqc/predistortion/P3_transfer_fit/plot_transfer_fit_quality_sqc.py
    python result_sqc/predistortion/P3_transfer_fit/plot_transfer_fit_quality_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "predistortion/P3_transfer_fit"
_CACHE = "transfer_fit_quality_sqc"

# Calibration measurement grid. WaveformCalibration._measure_step_response uses
# np.linspace(0, t_max, n_points) internally — this is NOT a control time axis
# (no waveform is played), so R9's arange rule does not apply here; we just
# mirror the calibration defaults.
T_MAX = 500.0
N_POINTS = 2000
K_SWEEP = (1, 2, 3, 4, 5)

# --- ground truths ------------------------------------------------------
# "single" reuses P1/P2's canonical single_exp so all three figures share physics.
TRUE_AMP, TRUE_TAU_NS = 0.05, 100.0
MULTI_AMPS, MULTI_TAUS = (0.04, 0.02), (80.0, 400.0)
CASCADE = ((0.05, 50.0), (0.03, 200.0))    # (amplitude, tau) per stage

# Reference parameters for the recovery panel. For the cascade the honest
# reference is NOT the per-stage amplitudes: cascading two single-exp stages
# gives H = H1*H2, whose step response is a 2-pole sum whose weights are the
# PARTIAL-FRACTION RESIDUES
#     r1 = A1(1-A2) + A1*A2*t1/(t1-t2),  r2 = (1-A1)A2 - A1*A2*t2/(t1-t2)
# = (0.0480, 0.0305) here, summing to the total step deficit 1-(1-A1)(1-A2).
# Poles are unchanged (50, 200 ns).
TRUTHS = ("single", "multi2", "cascade")
TRUE_ORDER = {"single": 1, "multi2": 2, "cascade": 2}
LABELS = {
    "single": r"single_exp  $A$=0.05, $\tau$=100 ns  (1 pole)",
    "multi2": r"multi_exp  (0.04, 80 ns)+(0.02, 400 ns)  (2 poles)",
    "cascade": r"cascade  (0.05, 50 ns)$\rightarrow$(0.03, 200 ns)  (2 poles)",
}
STYLE = {
    "single":  {"color": "#377eb8", "marker": "o"},
    "multi2":  {"color": "#e41a1c", "marker": "s"},
    "cascade": {"color": "#4daf4a", "marker": "^"},
}
K_STYLE = {1: {"color": "#e41a1c", "ls": "--"},
           2: {"color": "#377eb8", "ls": "-"},
           3: {"color": "#984ea3", "ls": ":"}}
FOCUS = "multi2"           # truth used in panels (a) and (b)


def true_distortion(name):
    """Ground-truth distortion model for one truth label."""
    from sqc.hardware.distortion import (
        SingleExponentialDistortion, MultiExponentialDistortion,
        CascadeDistortion,
    )
    if name == "single":
        return SingleExponentialDistortion(amplitude=TRUE_AMP, tau=TRUE_TAU_NS)
    if name == "multi2":
        return MultiExponentialDistortion(amplitudes=np.array(MULTI_AMPS),
                                          taus=np.array(MULTI_TAUS))
    return CascadeDistortion(stages=[
        SingleExponentialDistortion(amplitude=a, tau=t) for a, t in CASCADE])


def reference_params(name):
    """(amplitudes, taus) the fit should recover, largest-amplitude first.

    See the CASCADE comment above: the cascade reference is the analytic
    partial-fraction residues, NOT the per-stage amplitudes.
    """
    if name == "single":
        return np.array([TRUE_AMP]), np.array([TRUE_TAU_NS])
    if name == "multi2":
        return np.array(MULTI_AMPS), np.array(MULTI_TAUS)
    (a1, t1), (a2, t2) = CASCADE
    r1 = a1 * (1 - a2) + a1 * a2 * t1 / (t1 - t2)
    r2 = (1 - a1) * a2 - a1 * a2 * t2 / (t1 - t2)
    return np.array([r1, r2]), np.array([t1, t2])


def fit_once(name, K):
    """Fit truth `name` with K exponential components.

    Returns (t, measured_step, fitted_step, amplitudes, taus) with the
    recovered components sorted largest-amplitude-first so they can be
    compared against reference_params() element-wise.
    """
    from sqc.calibration.waveform import WaveformCalibration

    cal = WaveformCalibration(distortion=true_distortion(name),
                              method="transfer_function",
                              fit_type="multi_exp", n_exp_components=K,
                              t_max=T_MAX, n_points=N_POINTS)
    table = cal.calibrate()
    fitted = np.asarray(cal.to_distortion_model().step_response(table.inputs))
    fp = table.fit_params
    amps = np.atleast_1d(np.asarray(fp["amplitudes"], dtype=float))
    taus = np.atleast_1d(np.asarray(fp["taus"], dtype=float))
    order = np.argsort(amps)[::-1]
    return table.inputs, table.outputs, fitted, amps[order], taus[order]


def compute(recompute=False):
    """Fit every (truth, K) pair; cache steps, rms and recovered params."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    out = {}
    for name in TRUTHS:
        ref_a, ref_t = reference_params(name)
        out[f"{name}_ref_amps"] = ref_a
        out[f"{name}_ref_taus"] = ref_t
        rms_vs_k = []
        for K in K_SWEEP:
            t, meas, fitted, amps, taus = fit_once(name, K)
            rms_vs_k.append(float(np.sqrt(np.mean((fitted - meas) ** 2))))
            if "t" not in out:
                out["t"] = t
            if name == FOCUS and K in K_STYLE:
                out[f"focus_fit_K{K}"] = fitted
            if K == TRUE_ORDER[name]:
                # at the correct order the components line up with the
                # reference one-to-one -> relative recovery error
                out[f"{name}_fit_amps"] = amps
                out[f"{name}_fit_taus"] = taus
                o_f, o_r = np.argsort(taus), np.argsort(ref_t)
                out[f"{name}_relerr_amps"] = np.abs(
                    amps[o_f] / ref_a[o_r] - 1.0)
                out[f"{name}_relerr_taus"] = np.abs(
                    taus[o_f] / ref_t[o_r] - 1.0)
            if K == 1:
                out[f"{name}_k1_amp"] = np.array([amps[0]])
                out[f"{name}_k1_tau"] = np.array([taus[0]])
        out[f"{name}_rms"] = np.array(rms_vs_k)
        if name == FOCUS:
            out["focus_measured"] = meas

    out["k_sweep"] = np.array(K_SWEEP, dtype=float)
    C.save_npz(SUBDIR, _CACHE, **out)
    return out


_RES_FLOOR = 1e-16   # log-y clip: exact fits reach machine precision


def _panel_step(ax, d):
    """(a) measured vs fitted step for the 2-pole focus truth."""
    t, meas = d["t"], d["focus_measured"]
    ax.plot(t, meas, color="k", lw=2.2, alpha=0.45,
            label="measured step (ground truth)")
    for K in sorted(K_STYLE):
        ax.plot(t, d[f"focus_fit_K{K}"], lw=1.5, **K_STYLE[K],
                label=f"fitted, K={K}")
    ax.set_xlabel("Time (ns)", fontsize=11)
    ax.set_ylabel("Step response $s(t)$", fontsize=11)
    ax.set_title("(a) Measured vs fitted step, 2-pole truth\n"
                 "K=1 underfits; K=2 lands on the data", fontsize=10.5)
    ax.set_xlim(0, T_MAX)
    ax.legend(fontsize=8, loc="lower right", framealpha=0.92)
    ax.grid(True, ls=":", alpha=0.4)
    # No zoom inset here: the step only spans 0.94..1.0, so this panel's
    # natural y-range already IS the zoom, and an inset would just cover the
    # curves it is meant to magnify.


def _panel_residual(ax, d):
    """(b) |fitted - measured| on log-y for K=1/2/3."""
    t, meas = d["t"], d["focus_measured"]
    for K in sorted(K_STYLE):
        res = np.maximum(np.abs(d[f"focus_fit_K{K}"] - meas), _RES_FLOOR)
        ax.semilogy(t, res, lw=1.4, **K_STYLE[K], label=f"K={K}")
    ax.set_xlabel("Time (ns)", fontsize=11)
    ax.set_ylabel(r"Fit residual $|s_{fit} - s_{meas}|$", fontsize=11)
    ax.set_title("(b) Fit residual (log-y): the underfit floor is\n"
                 "8 orders above the correct order", fontsize=10.5)
    ax.set_xlim(0, T_MAX)
    ax.legend(fontsize=8, loc="center right", framealpha=0.92)
    ax.grid(True, which="both", ls=":", alpha=0.4)


def _panel_order(ax, d):
    """(c) fit rms vs K, with each truth's real pole count marked."""
    k = d["k_sweep"]
    for name in TRUTHS:
        st = STYLE[name]
        rms = np.maximum(d[f"{name}_rms"], _RES_FLOOR)
        ax.semilogy(k, rms, st["marker"] + "-", color=st["color"],
                    ms=6, lw=1.4, label=LABELS[name])
        # ring the true order
        i = list(K_SWEEP).index(TRUE_ORDER[name])
        ax.semilogy([k[i]], [rms[i]], "o", ms=13, mfc="none",
                    mec=st["color"], mew=1.8)
        # a truly exact fit gives rms == 0, which a log axis cannot show;
        # say so rather than letting the clipped marker imply "very small"
        if d[f"{name}_rms"][i] == 0.0:
            ax.annotate("exact (rms = 0)", xy=(k[i], rms[i]),
                        xytext=(k[i] + 0.30, rms[i] * 300), fontsize=8,
                        color=st["color"],
                        arrowprops=dict(arrowstyle="->", lw=0.9,
                                        color=st["color"]))
    ax.set_xlabel("Fitted components $K$ (n_exp_components)", fontsize=11)
    ax.set_ylabel("Fit rms over the step window", fontsize=11)
    ax.set_title("(c) Model-order selection (circled = true pole count)\n"
                 "rms collapses at $K=K_{true}$, then scatters", fontsize=10.5)
    ax.set_xticks(list(K_SWEEP))
    # upper right: no curve reaches above 1e-3 for K>=3, and the lower left is
    # taken by the exact-fit marker + its annotation
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.92)
    ax.grid(True, which="both", ls=":", alpha=0.4)


def _panel_recovery(ax, d):
    """(d) relative recovery error of each (A_k, tau_k) at the correct K."""
    names, vals, cols, hatches = [], [], [], []
    for name in TRUTHS:
        st = STYLE[name]
        ra, rt = d[f"{name}_relerr_amps"], d[f"{name}_relerr_taus"]
        taus = np.sort(d[f"{name}_ref_taus"])
        for j, tau in enumerate(taus):
            names.append(rf"{name}" "\n" rf"$\tau$={tau:.0f} ns")
            vals.append(max(ra[j], _RES_FLOOR)); cols.append(st["color"])
            hatches.append("")
            names.append("")
            vals.append(max(rt[j], _RES_FLOOR)); cols.append(st["color"])
            hatches.append("///")

    x = np.arange(len(vals))
    bars = ax.bar(x, vals, color=cols, edgecolor="white", linewidth=0.6)
    for b, h in zip(bars, hatches):
        b.set_hatch(h)
    ax.set_yscale("log")
    # Same honesty fix as panel (c): an exactly-recovered parameter has
    # relative error 0, unplottable on log-y. Label it instead of showing a
    # clipped stub that reads as "small but nonzero".
    for xi, v in zip(x, vals):
        if v <= _RES_FLOOR:
            ax.annotate("exact", xy=(xi, _RES_FLOOR), rotation=90,
                        ha="center", va="bottom", fontsize=7.5, color="dimgray")
    ax.set_xticks(x[::2] + 0.5)
    ax.set_xticklabels([n for n in names if n], fontsize=7.5)
    ax.set_ylabel("Relative recovery error $|fit/true - 1|$", fontsize=11)
    ax.set_title("(d) Parameter recovery at $K=K_{true}$\n"
                 "solid = amplitude, hatched = $\\tau$", fontsize=10.5)
    ax.grid(True, which="both", axis="y", ls=":", alpha=0.4)


def plot(d, out):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.8))
    (ax_a, ax_b), (ax_c, ax_e) = axes
    _panel_step(ax_a, d)
    _panel_residual(ax_b, d)
    _panel_order(ax_c, d)
    _panel_recovery(ax_e, d)

    fig.suptitle("Transfer-function fit quality and model-order selection "
                 "(sqc)", fontsize=13)
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

    for name in TRUTHS:
        rms = d[f"{name}_rms"]
        best = int(d["k_sweep"][int(np.argmin(rms))])
        print(f"  {name:8s} true order={TRUE_ORDER[name]}  "
              f"rms(K=1..{K_SWEEP[-1]})=" +
              " ".join(f"{v:.2e}" for v in rms) + f"  argmin K={best}")
        print(f"           K=1 compromise: A={d[f'{name}_k1_amp'][0]:.5f} "
              f"tau={d[f'{name}_k1_tau'][0]:.2f} ns")
        print(f"           at K={TRUE_ORDER[name]}: relerr A=" +
              " ".join(f"{v:.2e}" for v in d[f"{name}_relerr_amps"]) +
              "  tau=" +
              " ".join(f"{v:.2e}" for v in d[f"{name}_relerr_taus"]))

    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
