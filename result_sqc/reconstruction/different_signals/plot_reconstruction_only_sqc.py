#!/usr/bin/env python3
"""sqc-only reconstruction panel — original vs sqc Wiener/LM, NO src overlay.

Companion to plot_signal_comparison_panel_sqc.py (which overlays the frozen
src baseline). This figure is the clean "does sqc recover the field" story
for a report/paper: each cell shows only the original signal and sqc's own
Wiener + LM reconstructions.

Two data paths:
  - default: load the arrays cached by generate_signal_data_sqc.py (fast).
  - --recompute: re-run the full sqc API for each signal (slow; LM is
    ~344 s/signal full-rho inversion). The API calls are spelled out in
    recompute_one() so this script is self-contained.

Produces:
    result_sqc/reconstruction/different_signals/reconstruction_only_sqc.png (+ .pdf)

Usage:
    python result_sqc/reconstruction/different_signals/plot_reconstruction_only_sqc.py
    python result_sqc/reconstruction/different_signals/plot_reconstruction_only_sqc.py --recompute
    python result_sqc/reconstruction/different_signals/plot_reconstruction_only_sqc.py --recompute --no-lm
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "reconstruction/different_signals"
SIGNAL_ORDER = ["sine", "step", "double_peak", "complex"]
DESCRIPTIONS = {"sine": "Sinusoidal", "step": "Step-like",
                "double_peak": "Double-peak", "complex": "Complex"}

# src default signal params per type (must match generate_signal_data_sqc.py)
SIGNAL_PARAMS = {
    2: {"amplitude": 0.01, "frequency": 0.01},
    4: {"amplitude": 0.01, "center": 100, "rise": 10, "fall": 10},
    5: {"amplitude": 0.01, "center": 100, "width": 40},
    7: {"amplitude": 0.01},
}

C_ORIG, C_WIENER, C_LM = "black", "#e41a1c", "#377eb8"


def load_cached(name):
    """Load the arrays cached by generate_signal_data_sqc.py, or None."""
    path = os.path.join(C.out_dir(SUBDIR), f"signal_{name}_sqc.npz")
    if not os.path.exists(path):
        return None
    d = dict(np.load(path, allow_pickle=True))
    return {
        "t": d["original_time"], "orig": d["original_signal"],
        "wiener": d["wiener_recon"], "wiener_ratio": float(d["wiener_ratio"]),
        "lm": d.get("lm_recon"),
        "lm_ratio": float(d["lm_ratio"]) if "lm_ratio" in d else None,
    }


def recompute_one(signal_type, with_lm=True):
    """Re-run the full sqc API for one signal (self-contained). Returns dict.

    This is the same pipeline generate_signal_data_sqc.py uses, inlined here
    so the figure can be reproduced from scratch without the cache.
    """
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.control.flux_signal import FluxSignal

    q = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    fs = FluxSignal(type=signal_type, t_list=C.T_LIST, **SIGNAL_PARAMS[signal_type])
    exp = TransientSensingExperiment(qubit=q, flux_signal=fs)
    res = exp.run()
    orig = res.data["flux_samples"]

    # Wiener (lambda=5)
    wiener = TransientReconstruction(method="wiener", lambda_reg=C.WIENER_LAMBDA)
    rec_w = wiener.reconstruct(res, kernel=res.data["kernel"])

    out = {"t": res.axes["t_flux"], "orig": orig,
           "wiener": rec_w.signal, "wiener_ratio": C.peak_ratio(rec_w.signal, orig),
           "lm": None, "lm_ratio": None}

    # LM (fourier) — needs adapt_for_lm; returns (FluxSignal, history)
    if with_lm:
        lm = TransientReconstruction(
            method="lm", qubit=q, control_pulse=exp.control_pulse,
            basis_type="fourier", n_basis=C.LM_N_BASIS,
            lambda_reg=C.LM_LAMBDA, max_iter=C.LM_MAX_ITER)
        rec_lm, _ = lm.reconstruct(C.adapt_for_lm(res), kernel=res.data["kernel"])
        out["lm"] = rec_lm.signal
        out["lm_ratio"] = C.peak_ratio(rec_lm.signal, orig)
    return out


def _cell(ax, d, col):
    """One column cell: original + sqc Wiener + sqc LM (no src)."""
    ax.plot(d["t"], d["orig"], color=C_ORIG, ls="--", lw=1.5,
            label="Original" if col == 0 else "")
    ax.plot(d["t"], d["wiener"], color=C_WIENER, ls="-", lw=1.8,
            label="Wiener (sqc)" if col == 0 else "")
    if d.get("lm") is not None:
        ax.plot(d["t"], d["lm"], color=C_LM, ls="-", lw=1.8,
                label="LM (sqc)" if col == 0 else "")
    bits = [f"W r={d['wiener_ratio']:.2f}"]
    if d.get("lm_ratio") is not None:
        bits.append(f"LM r={d['lm_ratio']:.2f}")
    ax.text(0.97, 0.05, "\n".join(bits), transform=ax.transAxes,
            fontsize=8, ha="right", va="bottom")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    recompute = "--recompute" in sys.argv
    with_lm = "--no-lm" not in sys.argv

    loaded = {}
    for name in SIGNAL_ORDER:
        st = {v: k for k, v in C.SIGNAL_NAMES.items()}[name]
        d = recompute_one(st, with_lm=with_lm) if recompute else load_cached(name)
        if d is None:
            raise SystemExit(f"no cache for {name}; run generate_signal_data_sqc.py all "
                             f"or pass --recompute")
        loaded[name] = d

    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), sharex=True)
    for col, name in enumerate(SIGNAL_ORDER):
        ax = axes[col]
        _cell(ax, loaded[name], col)
        ax.set_title(DESCRIPTIONS[name], fontsize=11, pad=8)
        ax.set_xlabel("Time (ns)", fontsize=10)
        if col == 0:
            ax.set_ylabel("Field (arb.)", fontsize=10)
        ax.grid(True, ls=":", alpha=0.4, lw=0.5)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fig.legend(handles, labels, loc="lower center", ncol=len(labels),
               fontsize=9, frameon=True, bbox_to_anchor=(0.5, 0.0))

    out = os.path.join(C.out_dir(SUBDIR), "reconstruction_only_sqc.png")
    fig.savefig(out, dpi=300)
    fig.savefig(out.replace(".png", ".pdf"))
    plt.close(fig)
    print(f"panel -> {out}")


if __name__ == "__main__":
    main()
