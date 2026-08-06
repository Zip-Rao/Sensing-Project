#!/usr/bin/env python3
"""F2 ⭐ — Transient frequency-measurement accuracy vs detuning (sqc).

sqc source:
    sqc.calibration.frequency.FrequencyMeasurement(method="transient")
      - order=1  → linear estimate  Δω = p_diff / G_freq
      - order=3  → cubic Newton correction  p_diff = G₁·δω + (G₃/6)·δω³
    with the cubic coefficient G₃ from TWO equivalent sources (cross-check):
      - g3_source="fit"          (odd-poly fit of p_diff(Δ), adaptive range)
      - g3_source="kernel_full"  (full off-diagonal ∭k₃ dt³, no Δ scan)

Physics story (⭐ paper centrepiece): the linear estimate biases low at
large |Δ| (the fringe folds); the cubic Newton correction widens the
accurate range. Also demonstrates the two G₃ routes agree — the resolution
of the G1_fit≈0.6× bias documented in [[transient-g3-sim-fullkernel]].

True detuning is injected by a DC flux offset from OPTIMAL_FLUX while the
drive stays at f(OPTIMAL_FLUX) (omega_d not overridden), so |Δ| grows with
|offset|. The G₃ calibration is cached by (pulse, omega_d) — fixed drive ⇒
it runs ONCE per order-3 variant, then is reused across all offset points.

Produces:
    result_sqc/frequency_calibration/F2_transient_accuracy/transient_accuracy_sqc.png / .pdf
        top: signed (measured − true) f01 vs true Δ, three curves;
        bottom: |error| (log-y) vs true Δ, with the 1 MHz line + cubic fold.
    result_sqc/frequency_calibration/F2_transient_accuracy/transient_accuracy_sqc.npz

Cost: HIGH — G₃ calibration ~42 mesolve (fit) once; ~2 mesolve per point per
variant. Cached to npz; pass --recompute to force a fresh scan.

Usage:
    python result_sqc/frequency_calibration/F2_transient_accuracy/plot_transient_accuracy_sqc.py
    python result_sqc/frequency_calibration/F2_transient_accuracy/plot_transient_accuracy_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "frequency_calibration/F2_transient_accuracy"
_CACHE = "transient_accuracy_sqc"
TWO_PI = 2.0 * np.pi

# Flux offsets from OPTIMAL_FLUX (Φ₀). Chosen so the induced true detuning Δ
# spans roughly ±25 MHz — past the cubic fold (|δω_fold|≈17.6 MHz) so the plot
# captures both the accurate window and the beyond-fold linear fallback.
# Non-time array ⇒ linspace OK (R9).
FLUX_OFFSETS = np.linspace(-0.028, 0.028, 13)

# Cubic-fold half-width (MHz), from _solve_cubic_detuning's turnover
# |δω_fold| = sqrt(-2·G_lin/G3) for the default 10 ns π/2 pulse. Drawn as a
# guide; beyond it the order-3 Newton solve gracefully falls back to linear.
FOLD_MHZ = 17.6

VARIANTS = [
    {"label": "order=1 (linear)", "order": 1, "g3_source": "fit",
     "color": "#377eb8", "marker": "o", "mfc": "#377eb8"},
    {"label": "order=3 (fit)", "order": 3, "g3_source": "fit",
     "color": "#e41a1c", "marker": "s", "mfc": "none"},
    {"label": "order=3 (kernel_full)", "order": 3, "g3_source": "kernel_full",
     "color": "#4daf4a", "marker": "^", "mfc": "none"},
]


def analytic_f(offset):
    """Analytic f01 (GHz) at OPTIMAL_FLUX + offset, from the qubit model."""
    q = C.make_sqc_qubit()
    q.change_flux(C.OPTIMAL_FLUX + float(offset))
    return q.frequency / TWO_PI


def measure_curve(order, g3_source):
    """Measured f01 (GHz) across FLUX_OFFSETS for one variant.

    Drive omega_d is left at qubit.frequency (sweet spot) so the offset shows
    up as a genuine detuning; the G₃ calibration is therefore cached and reused
    across all offsets.
    """
    from sqc.calibration.frequency import FrequencyMeasurement

    fm = FrequencyMeasurement(qubit=C.make_sqc_qubit(), method="transient",
                              order=order, g3_source=g3_source)
    return np.array([fm.measure(flux=float(o)) / TWO_PI for o in FLUX_OFFSETS])


def compute(recompute=False):
    """Measured f01 for all variants + analytic truth. Cached to npz."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    analytic = np.array([analytic_f(o) for o in FLUX_OFFSETS])  # GHz
    f0 = float(analytic[np.argmin(np.abs(FLUX_OFFSETS))])       # drive f_d (GHz)
    detuning_mhz = (analytic - f0) * 1e3                        # true Δ (MHz)

    measured = np.array([measure_curve(v["order"], v["g3_source"])
                         for v in VARIANTS])                    # (n_var, n_pt) GHz
    error_mhz = (measured - analytic[None, :]) * 1e3           # (n_var, n_pt) MHz

    data = {
        "offsets": FLUX_OFFSETS,
        "analytic": analytic,
        "detuning_mhz": detuning_mhz,
        "f_drive": np.array([f0]),
        "measured": measured,
        "error_mhz": error_mhz,
        "labels": np.array([v["label"] for v in VARIANTS]),
    }
    C.save_npz(SUBDIR, _CACHE, **data)
    return data


def plot(d, out):
    import matplotlib.pyplot as plt

    delta = d["detuning_mhz"]
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(8, 6.8), sharex=True,
        gridspec_kw={"height_ratios": [3, 2]})

    # cubic accurate window shading + fold lines (both panels)
    for ax in (ax0, ax1):
        ax.axvspan(-FOLD_MHZ, FOLD_MHZ, color="#4daf4a", alpha=0.10, zorder=0)
        ax.axvline(-FOLD_MHZ, color="#4daf4a", ls="--", lw=0.8, alpha=0.7, zorder=1)
        ax.axvline(FOLD_MHZ, color="#4daf4a", ls="--", lw=0.8, alpha=0.7, zorder=1)

    # -- top: signed error --
    ax0.axhline(0, color="k", lw=0.8, ls=":")
    for i, v in enumerate(VARIANTS):
        ax0.plot(delta, d["error_mhz"][i], v["marker"] + "-", color=v["color"],
                 ms=6, lw=1.3, mfc=v["mfc"], mew=1.4, label=v["label"], zorder=3)
    ax0.set_ylabel("Measured − true $f_{01}$ (MHz)", fontsize=11)
    ax0.set_title("F2 — Transient frequency-measurement accuracy vs detuning",
                  fontsize=12)
    ax0.annotate(r"cubic accurate window $|\Delta|\lesssim$"
                 f"{FOLD_MHZ:.0f} MHz",
                 xy=(0, 0), xytext=(0.5, 0.06), textcoords="axes fraction",
                 ha="center", fontsize=8.5, color="#2e7d32")
    ax0.legend(fontsize=9, loc="upper center"); ax0.grid(True, ls=":", alpha=0.4)

    # -- bottom: |error| log-y --
    for i, v in enumerate(VARIANTS):
        ax1.semilogy(delta, np.abs(d["error_mhz"][i]) + 1e-4, v["marker"] + "-",
                     color=v["color"], ms=6, lw=1.3, mfc=v["mfc"], mew=1.4,
                     label=v["label"], zorder=3)
    ax1.axhline(1.0, color="gray", lw=0.9, ls="--", label="1 MHz")
    ax1.set_xlabel(r"True detuning $\Delta = f(\Phi) - f(\Phi_{opt})$ (MHz)",
                   fontsize=11)
    ax1.set_ylabel("|error| (MHz)", fontsize=10)
    ax1.legend(fontsize=8, ncol=2, loc="upper center"); ax1.grid(True, ls=":", alpha=0.4)

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

    delta = d["detuning_mhz"]
    inside = np.abs(delta) <= FOLD_MHZ
    print(f"Δ range: {delta.min():.1f}–{delta.max():.1f} MHz, "
          f"{len(delta)} points ({int(inside.sum())} inside fold)")
    for i, v in enumerate(VARIANTS):
        err = d["error_mhz"][i]
        print(f"  {v['label']:22s}: MAE all={np.mean(np.abs(err)):6.2f} MHz  "
              f"MAE inside-fold={np.mean(np.abs(err[inside])):6.2f} MHz")

    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
