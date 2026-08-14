#!/usr/bin/env python3
"""Generate final Fig. 2(c): held-out validation accuracy summary."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from fig2_simulation_common import ESTIMATOR_DATA, SUMMARY_DATA, load_npz, save_npz


HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[3]
SOURCE = ESTIMATOR_DATA
STEM = "fig2c-validation-summary"
STYLE = ROOT.parent / "scholaraio" / ".claude" / "skills" / "draw" / "nature_pub.mplstyle"


def compute(source: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    delta_val = float(source["delta_val_mhz"])
    mask = np.abs(source["validation_delta_mhz"]) <= delta_val
    errors = np.asarray(source["validation_error_mhz"][:, mask], dtype=float)
    if not np.all(np.isfinite(errors)):
        raise ValueError("all estimators must be finite on the frozen validation set")
    bias = np.mean(errors, axis=1)
    mae = np.mean(np.abs(errors), axis=1)
    rmse = np.sqrt(np.mean(errors**2, axis=1))
    max_abs = np.max(np.abs(errors), axis=1)
    return {
        "labels": source["labels"],
        "bias_mhz": bias,
        "mae_mhz": mae,
        "rmse_mhz": rmse,
        "max_abs_error_mhz": max_abs,
        "metrics_mhz": np.column_stack((mae, rmse, max_abs)),
        "metric_labels": np.array(("MAE", "RMSE", "Maximum")),
        "validation_delta_mhz": source["validation_delta_mhz"][mask],
        "validation_count": np.array(int(np.sum(mask))),
        "delta_val_mhz": np.array(delta_val),
        "delta_fold_fit_mhz": source["delta_fold_mhz"],
        "delta_fold_kernel_mhz": source["delta_fold_kernel_mhz"],
        "sigma": source["sigma"],
        "uncertainty_source": source["uncertainty_source"],
        "validation_role": source["validation_role"],
        "truth_role": source["truth_role"],
        "source_npz": np.array(SOURCE.name),
    }


def plot(data: dict[str, np.ndarray], stem: str = STEM, panel_label: str = "c") -> None:
    if STYLE.exists():
        plt.style.use(STYLE)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 7.0})
    colors = ("#3572B0", "#C43C39", "#2E8B57")
    markers = ("o", "s", "^")
    labels = [str(x) for x in data["labels"]]

    fig, ax_metric = plt.subplots(figsize=(3.50, 2.70))
    x = np.arange(4)
    offsets = (-0.16, 0.0, 0.16)
    for idx, (label, color, marker, offset) in enumerate(
        zip(labels, colors, markers, offsets)
    ):
        values = np.concatenate(([data["bias_mhz"][idx]], data["metrics_mhz"][idx]))
        ax_metric.plot(x + offset, values, linestyle="none", marker=marker,
                       markersize=4.8, markerfacecolor="white",
                       markeredgecolor=color, markeredgewidth=0.9,
                       label=label, zorder=3)
    ax_metric.set_yscale("symlog", linthresh=1e-5, linscale=0.75, base=10)
    ax_metric.axhline(0.0, color="#777777", lw=0.65, zorder=1)
    ax_metric.set_xticks(x, np.concatenate((["Bias"], data["metric_labels"])))
    ax_metric.set_ylabel("Frequency error (MHz)")
    ax_metric.set_ylim(-1.2e-5, 3.2)
    ax_metric.grid(axis="y", which="both", color="#DDDDDD", lw=0.45, alpha=0.8)
    ax_metric.text(
        0.98, 0.66,
        rf"held-out $n={int(data['validation_count'])}$, "
        rf"$|\Delta|/(2\pi)\leq {float(data['delta_val_mhz']):.0f}$ MHz",
        transform=ax_metric.transAxes, ha="right", va="top",
        fontsize=5.45, color="#555555",
    )
    handles, legend_labels = ax_metric.get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="upper center", bbox_to_anchor=(0.55, 0.97),
               ncol=3, frameon=False, handletextpad=0.3, columnspacing=0.8,
               fontsize=5.6)
    fig.text(0.02, 0.955, panel_label, ha="left", va="top",
             fontsize=8.5, fontweight="bold")
    ax_metric.text(0.08, 0.18, "bias: numerical residual",
                   transform=ax_metric.transAxes, ha="left", va="bottom",
                   fontsize=5.3, color="#555555")
    fig.subplots_adjust(left=0.20, right=0.97, bottom=0.20, top=0.79)
    fig.savefig(HERE / f"{stem}.pdf")
    fig.savefig(HERE / f"{stem}.svg")
    fig.savefig(HERE / f"{stem}.png", dpi=600)
    plt.close(fig)


def main() -> None:
    source = load_npz(SOURCE)
    data = compute(source)
    save_npz(SUMMARY_DATA, data)
    plot(data)
    plot(data, "figS7-validation-simulation", "S7")
    for idx, label in enumerate(data["labels"]):
        print(
            f"{label}: bias={data['bias_mhz'][idx]:.9f}, "
            f"MAE={data['mae_mhz'][idx]:.6f}, "
            f"RMSE={data['rmse_mhz'][idx]:.6f}, "
            f"max={data['max_abs_error_mhz'][idx]:.6f} MHz"
        )


if __name__ == "__main__":
    main()
