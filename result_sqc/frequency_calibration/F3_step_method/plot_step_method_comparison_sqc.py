#!/usr/bin/env python3
"""F3 — Closed-loop step-method convergence comparison (sqc).

sqc source:
    sqc.calibration.frequency.SinglePointFrequencyCalibration, run three times
    with step_method in {"secant", "bisection", "gradient"} toward the same
    f_target. Each returns a CalibrationTable whose fit_params["history"] is
    a per-iteration list of {iter, V, f, residual, ...}.

Physics story (paper role): root-finding algorithm comparison. Secant is
superlinear (1–3 iters); bisection gives a clean exponential bracket-width
trend (ideal for the log-scale panel); gradient needs no pre-bracketing.

Produces:
    result_sqc/frequency_calibration/F3_step_method/step_method_comparison_sqc.png / .pdf
        left: |residual| (MHz, log-y) vs iteration, three methods;
        right: bisection bracket-width (log-y) vs iteration (exponential O(log)).
    result_sqc/frequency_calibration/F3_step_method/step_method_comparison_sqc.npz

Cost: medium~high — each method: several iters × (Ramsey double-sweep ≈ 8.6 s
per measurement). Cached to npz; pass --recompute to force a fresh run.

Usage:
    python result_sqc/frequency_calibration/F3_step_method/plot_step_method_comparison_sqc.py
    python result_sqc/frequency_calibration/F3_step_method/plot_step_method_comparison_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "frequency_calibration/F3_step_method"
_CACHE = "step_method_comparison_sqc"
TWO_PI = 2.0 * np.pi
METHODS = ["secant", "bisection", "gradient"]
# f_target 20 MHz below the sweet spot. IMPORTANT: V here is a flux OFFSET from
# the qubit construction point (OPTIMAL_FLUX), NOT absolute flux —
# SinglePointFrequencyCalibration forwards V to FrequencyMeasurement.measure(flux=V)
# which adds it on top of the base work point (the [[sqc-frequency-measurement-api]]
# "flux is an OFFSET" convention). f(Φ) is monotone here, so the −20 MHz target
# (root at offset ≈ −0.016) sits strictly inside f(offset=−0.05)≈−55 MHz ..
# f(offset=+0.01)≈+7 MHz → the bracket straddles the root (secant/bisection converge).
F_TARGET_OFFSET_MHZ = -20.0
V_BRACKET = (-0.05, 0.01)  # V_a, V_b (OFFSETS) bracketing the target
V_SEED = 0.0               # gradient start (offset): sweet spot

STYLE = {
    "secant":    {"color": "#377eb8", "marker": "o"},
    "bisection": {"color": "#e41a1c", "marker": "s"},
    "gradient":  {"color": "#4daf4a", "marker": "^"},
}


def run_method(step_method):
    """Run one closed-loop calibration. Returns fit_params['history']."""
    from sqc.calibration.frequency import SinglePointFrequencyCalibration

    q = C.make_sqc_qubit()
    f_target = q.frequency + TWO_PI * (F_TARGET_OFFSET_MHZ * 1e-3)
    kw = dict(qubit=q, f_target=f_target, step_method=step_method,
              V_a=V_BRACKET[0], V_b=V_BRACKET[1], measure_method="ramsey")
    if step_method == "gradient":
        kw["V_seed"] = V_SEED
    return SinglePointFrequencyCalibration(**kw).calibrate().fit_params["history"]


def compute(recompute=False):
    """Run all three methods, extract per-iter |residual| (MHz). Cached to npz.

    History iteration index differs per method (gradient logs iter 0, the
    bracketing methods start at iter 1); we re-index each on its own 0-based
    axis so the log-y panel compares convergence *rate*, not absolute count.
    """
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    out = {}
    for m in METHODS:
        hist = run_method(m)
        res_mhz = np.array([abs(h["residual"]) / TWO_PI * 1e3 for h in hist])
        out[f"{m}_res"] = res_mhz
        out[f"{m}_iter"] = np.arange(len(res_mhz), dtype=float)
        if m == "bisection":
            out["bisection_bw"] = np.array(
                [h.get("bracket_width", np.nan) for h in hist])
        print(f"{m}: {len(hist)} iters, final |res|={res_mhz[-1]:.4f} MHz")

    out["f_target_offset_mhz"] = np.array([F_TARGET_OFFSET_MHZ])
    C.save_npz(SUBDIR, _CACHE, **out)
    return out


def plot(d, out):
    import matplotlib.pyplot as plt

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    # -- left: |residual| log-y vs iteration --
    for m in METHODS:
        st = STYLE[m]
        ax0.semilogy(d[f"{m}_iter"], d[f"{m}_res"] + 1e-4, st["marker"] + "-",
                     color=st["color"], ms=6, lw=1.4, label=m)
    eps_mhz = 1e-4 / TWO_PI * 1e3  # epsilon_f=1e-4 rad·GHz → MHz
    ax0.axhline(eps_mhz, color="gray", ls="--", lw=0.9,
                label=f"tol ε={eps_mhz:.3f} MHz")
    ax0.set_xlabel("Iteration", fontsize=11)
    ax0.set_ylabel("|residual| $|f_q - f_{target}|$ (MHz)", fontsize=11)
    ax0.set_title(f"Closed-loop convergence (target {F_TARGET_OFFSET_MHZ:+.0f} MHz)",
                  fontsize=11)
    ax0.legend(fontsize=9); ax0.grid(True, which="both", ls=":", alpha=0.4)

    # -- right: bisection bracket-width log-y (exponential trend) --
    bw = d.get("bisection_bw")
    if bw is not None:
        it = np.arange(1, len(bw) + 1)
        ax1.semilogy(it, bw, "s-", color="#e41a1c", ms=6, lw=1.4,
                     label="measured bracket width")
        # ideal halving reference from the initial width
        ideal = bw[0] * 0.5 ** (it - 1)
        ax1.semilogy(it, ideal, "k:", lw=1.2, label=r"ideal $2^{-n}$ halving")
        ax1.set_xlabel("Iteration", fontsize=11)
        ax1.set_ylabel(r"Bracket width $V_{hi}-V_{lo}$ ($\Phi_0$)", fontsize=11)
        ax1.set_title("Bisection bracket contraction", fontsize=11)
        ax1.legend(fontsize=9); ax1.grid(True, which="both", ls=":", alpha=0.4)

    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{_CACHE}.png"), dpi=300)
    fig.savefig(os.path.join(out, f"{_CACHE}.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    recompute = "--recompute" in sys.argv
    d = compute(recompute=recompute)

    for m in METHODS:
        res = d[f"{m}_res"]
        print(f"  {m:10s}: {len(res)} iters, final |res|={res[-1]:.4f} MHz")

    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
