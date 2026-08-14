#!/usr/bin/env python3
"""Build report figures RF1 and RF2 from frozen frequency-calibration data.

RF1 is the single-panel f(Phi) comparison used in the completion report.
RF2 combines the drive-tracking state machine with deterministic static-loop
performance. Until per-iteration hardware or wall-clock timestamps are
available, RF2 uses counted solver calls as the auditable resource axis.

The script only reads the existing SRC-F1 and SRC-F4 NPZ files. It does not
rerun the physics simulations or modify their cached results.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np


HERE = Path(__file__).resolve().parent
FC_ROOT = HERE.parent
RF1_SOURCE = FC_ROOT / "F1_flux_curve" / "f_phi_curve_sqc.npz"
RF2_SOURCE = FC_ROOT / "F4_hybrid_convergence" / "freq_calibration_sqc.npz"

COLORS = {
    "ramsey": "#356AA0",
    "hybrid": "#D98324",
    "tracking": "#178F6A",
    "analytic": "#202124",
    "trusted": "#8BCF9B",
    "locked": "#7561A8",
    "muted": "#667085",
}

THRESHOLDS_MHZ = np.array([10.0, 5.0, 1.0, 0.1, 0.016])
ACCURATE_HALFWIDTH = 0.015


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7.0,
        "axes.labelsize": 7.0,
        "axes.titlesize": 7.5,
        "legend.fontsize": 6.3,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.25,
        "lines.markersize": 4.0,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
    })


def _load(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"Missing frozen source data: {path}")
    with np.load(path, allow_pickle=True) as data:
        return {key: data[key] for key in data.files}


def _save(fig: plt.Figure, stem: str) -> None:
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / f"{stem}.png", dpi=600)
    fig.savefig(HERE / f"{stem}.pdf")
    fig.savefig(HERE / f"{stem}.svg")
    plt.close(fig)


def plot_rf1(data: dict[str, np.ndarray]) -> None:
    flux = data["flux"]
    optimum = float(flux[np.argmin(np.abs(data["offsets"]))])
    win_lo = optimum - ACCURATE_HALFWIDTH
    win_hi = optimum + ACCURATE_HALFWIDTH

    fig, ax = plt.subplots(figsize=(3.50, 2.45))
    ax.axvspan(
        win_lo,
        win_hi,
        color=COLORS["trusted"],
        alpha=0.20,
        lw=0,
        label="Trusted range",
        zorder=0,
    )
    ax.axvline(optimum, color="#8A8F98", lw=0.8, ls="--", zorder=1)
    ax.plot(
        flux,
        data["analytic"],
        color=COLORS["analytic"],
        lw=1.45,
        label="Analytic",
        zorder=2,
    )
    ax.plot(
        flux,
        data["ramsey"],
        linestyle="none",
        marker="o",
        color=COLORS["ramsey"],
        markeredgecolor="white",
        markeredgewidth=0.45,
        label="Ramsey",
        zorder=3,
    )
    ax.plot(
        flux,
        data["transient"],
        linestyle="none",
        marker="s",
        markerfacecolor="white",
        markeredgecolor=COLORS["hybrid"],
        markeredgewidth=0.9,
        label="Short sequence",
        zorder=4,
    )
    ax.text(
        optimum + 0.0007,
        0.985,
        r"$\Phi_{\rm ref}$",
        transform=ax.get_xaxis_transform(),
        color=COLORS["muted"],
        fontsize=5.8,
        ha="left",
        va="top",
    )
    ax.set_xlabel(r"Flux bias $\Phi/\Phi_0$")
    ax.set_ylabel(r"Qubit frequency $f_{01}$ (GHz)")
    ax.set_xlim(float(flux.min()) - 0.002, float(flux.max()) + 0.002)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.grid(axis="y", color="#D9DDE3", lw=0.55, alpha=0.75)
    ax.legend(
        loc="lower right",
        ncol=2,
        fontsize=5.6,
        frameon=True,
        framealpha=0.96,
        edgecolor="#D9DDE3",
        borderpad=0.35,
        handlelength=1.45,
        columnspacing=0.9,
    )
    fig.tight_layout()
    _save(fig, "RF1_flux_curve_report")


def _node(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    detail: str,
    color: str,
) -> None:
    x, y = xy
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        facecolor="white",
        edgecolor=color,
        linewidth=1.15,
    )
    ax.add_patch(box)
    ax.text(
        x + width / 2,
        y + height * 0.66,
        title,
        ha="center",
        va="center",
        fontsize=6.3,
        fontweight="bold",
        color=color,
    )
    ax.text(
        x + width / 2,
        y + height * 0.28,
        detail,
        ha="center",
        va="center",
        fontsize=5.5,
        color="#3F4650",
    )


def _arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    label: str | None = None,
    connectionstyle: str = "arc3",
    color: str = "#667085",
) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=7,
        linewidth=0.8,
        color=color,
        connectionstyle=connectionstyle,
        shrinkA=2,
        shrinkB=2,
    )
    ax.add_patch(arrow)
    if label:
        x = (start[0] + end[0]) / 2
        y = (start[1] + end[1]) / 2
        ax.text(
            x,
            y,
            label,
            ha="center",
            va="center",
            fontsize=5.0,
            color=COLORS["muted"],
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.7},
        )


def _plot_state_machine(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    width, height = 0.56, 0.16
    x = 0.04
    _node(ax, (x, 0.79), width, height, "Wide-range acquisition", "Ramsey estimate", COLORS["ramsey"])
    _node(ax, (x, 0.54), width, height, "Drive tracking", "Short sequence", COLORS["tracking"])
    _node(ax, (x, 0.29), width, height, "Target lock", r"$\omega_d=\omega_{\rm tar}$", COLORS["locked"])
    _node(ax, (0.67, 0.54), 0.30, height, "Re-acquire", "Ramsey", COLORS["hybrid"])

    _arrow(ax, (0.32, 0.79), (0.32, 0.70), "local range")
    _arrow(ax, (0.32, 0.54), (0.32, 0.45), "target tolerance")
    _arrow(ax, (0.60, 0.64), (0.67, 0.64), None)
    ax.text(
        0.635,
        0.705,
        "out of range",
        ha="center",
        va="center",
        fontsize=5.0,
        color=COLORS["muted"],
    )
    _arrow(
        ax,
        (0.82, 0.54),
        (0.60, 0.56),
        None,
        connectionstyle="arc3,rad=-0.34",
    )
    _arrow(
        ax,
        (0.60, 0.34),
        (0.78, 0.54),
        None,
        connectionstyle="arc3,rad=-0.20",
        color=COLORS["hybrid"],
    )

    ax.text(
        0.04,
        0.05,
        r"$V$ update: target convergence" "\n" r"$\omega_d$ update: estimator validity",
        fontsize=5.7,
        color="#3F4650",
        ha="left",
        va="bottom",
    )
    ax.text(0.0, 1.01, "a", fontsize=8.0, fontweight="bold", ha="left", va="bottom")


def _cost_to_reach(residual: np.ndarray, cost: np.ndarray, threshold: float) -> float:
    hit = np.asarray(residual) <= threshold
    if not hit.any():
        return np.nan
    return float(np.asarray(cost)[np.flatnonzero(hit)[0]])


def plot_rf2(data: dict[str, np.ndarray]) -> None:
    strategies = {
        "Ramsey gradient": (
            data["base_verified_res_mhz"],
            data["base_time_true"],
            COLORS["ramsey"],
            "o",
        ),
        "Fixed-drive hybrid": (
            data["verified_res_mhz"],
            data["time_true"],
            COLORS["hybrid"],
            "s",
        ),
        "Drive tracking": (
            data["trk_verified_res_mhz"],
            data["trk_time_true"],
            COLORS["tracking"],
            "^",
        ),
    }

    fig = plt.figure(figsize=(7.20, 2.75))
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.14, 1.0, 1.0],
        left=0.025,
        right=0.985,
        bottom=0.22,
        top=0.94,
        wspace=0.36,
    )
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    _plot_state_machine(ax0)

    for label, (residual, elapsed, color, marker) in strategies.items():
        ax1.loglog(
            elapsed,
            residual,
            marker=marker,
            color=color,
            markerfacecolor="white",
            markeredgewidth=0.8,
            label=label,
        )
    ax1.axhline(0.016, color="#4B5563", ls="--", lw=0.8)
    ax1.text(
        0.98,
        0.016,
        " stop threshold",
        transform=ax1.get_yaxis_transform(),
        fontsize=5.4,
        color=COLORS["muted"],
        ha="right",
        va="bottom",
    )
    ax1.set_xlabel("Cumulative simulation time (s)")
    ax1.set_ylabel(r"Verified residual $|f_q-f_{\rm tar}|$ (MHz)")
    ax1.grid(which="both", color="#D9DDE3", lw=0.5, alpha=0.7)
    ax1.legend(frameon=False, loc="upper right", handlelength=1.5)
    ax1.text(-0.22, 1.01, "b", transform=ax1.transAxes, fontsize=8.0, fontweight="bold", va="bottom")

    for label, (residual, elapsed, color, marker) in strategies.items():
        first_time = np.array([
            _cost_to_reach(residual, elapsed, threshold)
            for threshold in THRESHOLDS_MHZ
        ])
        ax2.loglog(
            THRESHOLDS_MHZ,
            first_time,
            marker=marker,
            color=color,
            markerfacecolor="white",
            markeredgewidth=0.8,
            label=label,
        )
    ax2.invert_xaxis()
    ax2.set_xlabel("Residual threshold (MHz)")
    ax2.set_ylabel("Time to first reach threshold (s)")
    ax2.set_xticks(THRESHOLDS_MHZ)
    ax2.set_xticklabels(["10", "5", "1", "0.1", "0.016"])
    ax2.grid(which="both", color="#D9DDE3", lw=0.5, alpha=0.7)
    ax2.text(-0.22, 1.01, "c", transform=ax2.transAxes, fontsize=8.0, fontweight="bold", va="bottom")

    _save(fig, "RF2_closed_loop_report")


def main() -> None:
    _style()
    rf1 = _load(RF1_SOURCE)
    rf2 = _load(RF2_SOURCE)
    plot_rf1(rf1)
    plot_rf2(rf2)
    print(f"RF1 source: {RF1_SOURCE}")
    print(f"RF2 source: {RF2_SOURCE}")
    print(f"Output directory: {HERE}")


if __name__ == "__main__":
    main()
