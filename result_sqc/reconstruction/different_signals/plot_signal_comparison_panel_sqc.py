#!/usr/bin/env python3
"""Build the 4-signal comparison panel from saved sqc .npz files.

Baseline: result/different_signals/plot_signal_comparison_panel.py
Produces: result_sqc/reconstruction/different_signals/comparison_sqc.png (+ .pdf)

Reads signal_{sine,step,double_peak,complex}_sqc.npz produced by
generate_signal_data_sqc.py and lays out a grid of
original vs Wiener vs LM (with src baseline overlaid).
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

# ColorBrewer Set1 (colour-blind safe), matching src panel conventions.
C_ORIG, C_WIENER, C_LM = "black", "#e41a1c", "#377eb8"


def load(name):
    """Load one sqc npz as a dict, or None if missing."""
    path = os.path.join(C.out_dir(SUBDIR), f"signal_{name}_sqc.npz")
    if not os.path.exists(path):
        return None
    return dict(np.load(path, allow_pickle=True))


def _metrics(data):
    """Per-signal metric row: sqc vs src rmse/peak-ratio for wiener + LM."""
    m = {"wiener_rmse": float(data["wiener_rmse"]),
         "wiener_ratio": float(data["wiener_ratio"])}
    if "lm_rmse" in data:
        m["lm_rmse"] = float(data["lm_rmse"])
        m["lm_ratio"] = float(data["lm_ratio"])
    return m


def _cell(ax, data, row, col):
    """Draw one panel cell. row 0 = Wiener, row 1 = LM."""
    t = data["original_time"]
    ax.plot(t, data["original_signal"], color=C_ORIG, ls="--", lw=1.5,
            label="Original" if col == 0 else "")
    if row == 0:
        ax.plot(data["wiener_time"], data["wiener_recon"], color=C_WIENER,
                ls="-", lw=1.8, label="Wiener (sqc)" if col == 0 else "")
        if "src_wiener" in data:
            ax.plot(data["src_wiener_time"], data["src_wiener"], color=C_WIENER,
                    ls=":", lw=1.0, alpha=0.6, label="Wiener (src)" if col == 0 else "")
        ratio = float(data["wiener_ratio"])
    else:
        if "lm_recon" in data:
            ax.plot(data["lm_time"], data["lm_recon"], color=C_LM, ls="-",
                    lw=1.8, label="LM (sqc)" if col == 0 else "")
        if "src_lm" in data:
            ax.plot(data["src_lm_time"], data["src_lm"], color=C_LM, ls=":",
                    lw=1.0, alpha=0.6, label="LM (src)" if col == 0 else "")
        ratio = float(data.get("lm_ratio", np.nan))
    tag = chr(97 + col + row * 4)
    ax.text(0.02, 0.95, f"({tag})", transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")
    if np.isfinite(ratio):
        ax.text(0.97, 0.05, f"r={ratio:.2f}", transform=ax.transAxes,
                fontsize=8, ha="right", va="bottom")
    if row == 0:
        ax.set_title(DESCRIPTIONS.get(data_name(data), ""), fontsize=10, pad=8)
    else:
        ax.set_xlabel("Time (ns)", fontsize=10)
    if col == 0:
        ax.set_ylabel(("Wiener" if row == 0 else "LM") + "\nField (arb.)", fontsize=10)
    ax.grid(True, ls=":", alpha=0.4, lw=0.5)


def data_name(data):
    """Recover the signal name from the stored signal_type."""
    return C.SIGNAL_NAMES[int(data["signal_type"])]


def build_panel(loaded):
    """2x4 panel: rows = Wiener/LM, cols = signal types. Returns fig."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(12, 6.4), sharex="col")
    for col, name in enumerate(SIGNAL_ORDER):
        data = loaded[name]
        for row in range(2):
            _cell(axes[row, col], data, row, col)
    handles, labels = [], []
    for ax in axes.flat:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l and l not in labels:
                handles.append(h); labels.append(l)
    # reserve a bottom band for the shared legend so it never clips
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    fig.legend(handles, labels, loc="lower center", ncol=len(labels),
               fontsize=9, frameon=True, bbox_to_anchor=(0.5, 0.0))
    return fig


def write_metrics(loaded, path):
    """Write a compact sqc metrics table (rmse + peak-ratio per signal)."""
    lines = ["sqc reconstruction metrics (different_signals)", "=" * 60,
             f"{'signal':<14}{'method':<10}{'rmse':<14}{'peak-ratio':<12}", "-" * 60]
    for name in SIGNAL_ORDER:
        m = _metrics(loaded[name])
        lines.append(f"{name:<14}{'Wiener':<10}{m['wiener_rmse']:<14.4e}{m['wiener_ratio']:<12.3f}")
        if "lm_rmse" in m:
            lines.append(f"{'':<14}{'LM':<10}{m['lm_rmse']:<14.4e}{m['lm_ratio']:<12.3f}")
        lines.append("-" * 60)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    loaded = {n: load(n) for n in SIGNAL_ORDER}
    missing = [n for n, d in loaded.items() if d is None]
    if missing:
        raise SystemExit(f"missing npz for {missing}; run generate_signal_data_sqc.py all")

    out = C.out_dir(SUBDIR)
    fig = build_panel(loaded)
    png = os.path.join(out, "comparison_sqc.png")
    fig.savefig(png, dpi=300)
    fig.savefig(png.replace(".png", ".pdf"))
    plt.close(fig)

    metrics_path = os.path.join(out, "comparison_metrics_sqc.txt")
    write_metrics(loaded, metrics_path)
    print(f"panel -> {png}\nmetrics -> {metrics_path}")


if __name__ == "__main__":
    main()
