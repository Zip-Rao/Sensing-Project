#!/usr/bin/env python3
"""Generate Supplemental Fig. S2 from frozen response and estimator evidence."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fig2_simulation_common import ESTIMATOR_DATA, RESPONSE_DATA, TWO_PI, load_npz


HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[3]
STYLE = ROOT.parent / "scholaraio" / ".claude" / "skills" / "draw" / "nature_pub.mplstyle"
STEM = "figS2-response-simulation"


def main() -> None:
    response = load_npz(RESPONSE_DATA)
    estimator = load_npz(ESTIMATOR_DATA)
    if STYLE.exists():
        plt.style.use(STYLE)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 7.0})

    delta = response["dense_delta_mhz"]
    x = TWO_PI * delta * 1e-3
    fit_fold = float(estimator["delta_fold_mhz"])
    kernel_fold = float(estimator["delta_fold_kernel_mhz"])
    interval = float(estimator["delta_val_mhz"])
    fit_model = float(estimator["G1_fit_ns"]) * x + float(estimator["G3_fit_ns3"]) * x**3 / 6.0
    kernel_model = (
        float(estimator["G1_kernel_ns"]) * x
        + float(estimator["G3_kernel_ns3"]) * x**3 / 6.0
    )
    fit_mask = np.abs(delta) <= fit_fold
    kernel_mask = np.abs(delta) <= kernel_fold

    fig, ax = plt.subplots(figsize=(3.50, 2.85))
    ax.axvspan(-interval, interval, color="#D9EAD3", alpha=0.75, lw=0)
    for fold, color in ((fit_fold, "#B65A55"), (kernel_fold, "#16705A")):
        for sign in (-1.0, 1.0):
            ax.axvline(sign * fold, color=color, ls="--", lw=0.7, alpha=0.8)
    ax.plot(delta, response["dense_pdiff_corrected"], color="#222222", lw=1.2,
            label="Dissipative master equation")
    ax.plot(delta[fit_mask], fit_model[fit_mask], color="#B65A55", lw=1.05,
            label="Cubic fit")
    ax.plot(delta[kernel_mask], kernel_model[kernel_mask], color="#16705A", lw=1.05,
            label="Full-kernel cubic")
    ax.scatter(response["calibration_delta_mhz"], response["calibration_pdiff_corrected"],
               s=15, color="#3572B0", edgecolor="white", linewidth=0.4,
               label="Response calibration", zorder=4)
    ax.scatter(response["validation_delta_mhz"], response["validation_pdiff_corrected"],
               s=18, marker="D", facecolor="white", edgecolor="#B65A55", linewidth=0.75,
               label="Simulated held-out", zorder=5)
    fig.text(0.18, 0.055,
             rf"fit: $G_1={float(estimator['G1_fit_ns']):.3f}$ ns, "
             rf"$G_3={float(estimator['G3_fit_ns3']):.1f}$ ns$^3$, fold {fit_fold:.2f} MHz"
             "\n"
             rf"kernel: $G_1={float(estimator['G1_kernel_ns']):.3f}$ ns, "
             rf"$G_3={float(estimator['G3_kernel_ns3']):.1f}$ ns$^3$, fold {kernel_fold:.2f} MHz",
             ha="left", va="bottom", fontsize=4.8, color="#444444")
    ax.set(xlabel=r"True detuning $\Delta/(2\pi)$ (MHz)",
           ylabel=r"Corrected differential response $\widetilde p_{\rm d}$",
           xlim=(-28.0, 28.0))
    fig.text(0.02, 0.96, "S2", ha="left", va="top",
             fontsize=8.5, fontweight="bold")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
              frameon=False, fontsize=5.35, columnspacing=0.8)
    fig.subplots_adjust(left=0.18, right=0.97, bottom=0.27, top=0.78)
    for suffix, kwargs in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
        fig.savefig(HERE / f"{STEM}.{suffix}", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
