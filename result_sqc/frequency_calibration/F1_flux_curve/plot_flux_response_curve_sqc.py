#!/usr/bin/env python3
"""F1 — Flux-response calibration curve f(Φ): Ramsey vs transient (sqc).

sqc source:
    sqc.calibration.frequency.FrequencyMeasurement, one instance per method
    ("ramsey", "transient"), measured across a DC flux offset scan. Each
    .measure(flux=offset) returns f01 (angular) at optimal_work_point + offset.

Physics story (report/paper role): the transmon flux-dispersion curve f(Φ) is
the device calibration primitive. We recover it two ways and overlay both on
the analytic curve f(Φ) = √(8·E_J(Φ)·E_C) − E_C (sampled from the qubit model):
  - Ramsey  — τ-sweep FFT, robust, costlier (the standard method).
  - transient — τ=0 orthogonal readout + kernel sensitivity, cheap, relies on
    the weak-signal linear (+cubic) approximation.
The residual subpanel shows both track the analytic dispersion to <~0.1 MHz
over the scanned range, i.e. the cheap transient probe reproduces the Ramsey
calibration.

Produces:
    result_sqc/frequency_calibration/F1_flux_curve/f_phi_curve_sqc.png / .pdf
        top: analytic f(Φ) + Ramsey pts + transient pts;
        bottom: residual (measured − analytic) in MHz for both methods.
    result_sqc/frequency_calibration/F1_flux_curve/f_phi_curve_sqc.npz

Cost: medium — 2 methods × len(OFFSETS) single-point solves. Cached to npz;
pass --recompute to force a fresh scan.

Usage:
    python result_sqc/frequency_calibration/F1_flux_curve/plot_flux_response_curve_sqc.py
    python result_sqc/frequency_calibration/F1_flux_curve/plot_flux_response_curve_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "frequency_calibration/F1_flux_curve"
_CACHE = "f_phi_curve_sqc"
# flux offsets from optimal work point (Φ₀); non-time array, linspace OK.
OFFSETS = np.linspace(-0.03, 0.03, 21)
TWO_PI = 2.0 * np.pi
# transient accurate window (half-width in flux offset). Inside |offset|<=this,
# the detuning stays below the cubic fold turnover (δω*≈0.111 rad/ns ≈ 17.6 MHz)
# so the transient weak-signal inversion reproduces Ramsey to <~1 MHz.
ACCURATE_HALFWIDTH = 0.015


def analytic_f(offset):
    """Analytic f01 (GHz) at optimal_work_point + offset, from the qubit model."""
    q = C.make_sqc_qubit()
    q.change_flux(C.OPTIMAL_FLUX + float(offset))
    return q.frequency / TWO_PI


def compute(recompute=False):
    """Measure f(Φ) with both methods across OFFSETS. Cached to npz."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path)
        return {k: d[k] for k in d.files}

    from sqc.calibration.frequency import FrequencyMeasurement

    # base qubit sits at optimal; .measure(flux=offset) applies the offset.
    fm_ramsey = FrequencyMeasurement(qubit=C.make_sqc_qubit(), method="ramsey",
                                     f_artificial=0.1)
    fm_transient = FrequencyMeasurement(qubit=C.make_sqc_qubit(), method="transient",
                                        order=3, f_artificial=0.1)

    analytic = np.array([analytic_f(o) for o in OFFSETS])
    ramsey = np.array([fm_ramsey.measure(flux=float(o)) / TWO_PI for o in OFFSETS])
    transient = np.array([fm_transient.measure(flux=float(o)) / TWO_PI for o in OFFSETS])

    data = {"offsets": OFFSETS, "flux": C.OPTIMAL_FLUX + OFFSETS,
            "analytic": analytic, "ramsey": ramsey, "transient": transient}
    C.save_npz(SUBDIR, _CACHE, **data)
    return data


def plot(d, out):
    import matplotlib.pyplot as plt

    flux = d["flux"]
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(8, 6.4), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]})

    # transient accurate window (flux coords) + drive point
    win_lo = C.OPTIMAL_FLUX - ACCURATE_HALFWIDTH
    win_hi = C.OPTIMAL_FLUX + ACCURATE_HALFWIDTH
    f_drive = float(d["analytic"][np.argmin(np.abs(d["offsets"]))])

    # -- top: f(Φ) --
    ax0.axvspan(win_lo, win_hi, color="#4daf4a", alpha=0.12, zorder=0,
                label=r"Transient accurate window ($|\Delta|\lesssim$17 MHz)")
    ax0.axvline(C.OPTIMAL_FLUX, color="gray", ls="--", lw=0.9, zorder=1)
    ax0.axhline(f_drive, color="gray", ls=":", lw=0.9, zorder=1)
    ax0.plot(flux, d["analytic"], "k-", lw=1.6, label="Analytic f(Φ)", zorder=2)
    ax0.plot(flux, d["ramsey"], "o", color="#377eb8", ms=6, label="Ramsey", zorder=3)
    ax0.plot(flux, d["transient"], "s", color="#e41a1c", ms=5, mfc="none",
             mew=1.4, label="Transient", zorder=4)
    ax0.annotate(f"drive $f_d$={f_drive:.4f} GHz\n"
                 r"($\Delta=0$ at $\Phi_{opt}$)",
                 xy=(C.OPTIMAL_FLUX, f_drive), xytext=(0.04, 0.72),
                 textcoords="axes fraction", fontsize=8.5, color="dimgray",
                 arrowprops=dict(arrowstyle="->", color="gray", lw=1))
    ax0.set_ylabel("Qubit frequency $f_{01}$ (GHz)", fontsize=11)
    ax0.set_title("F1 — Flux-response calibration f(Φ): Ramsey vs transient",
                  fontsize=12)
    ax0.legend(fontsize=9, loc="lower right"); ax0.grid(True, ls=":", alpha=0.4)

    # -- bottom: residual (MHz) --
    res_r = (d["ramsey"] - d["analytic"]) * 1e3
    res_t = (d["transient"] - d["analytic"]) * 1e3
    ax1.axvspan(win_lo, win_hi, color="#4daf4a", alpha=0.12, zorder=0)
    ax1.axvline(C.OPTIMAL_FLUX, color="gray", ls="--", lw=0.9, zorder=1)
    ax1.axhline(0, color="k", lw=0.8, ls=":")
    ax1.axhline(1, color="#4daf4a", lw=0.7, ls="--", alpha=0.7)
    ax1.axhline(-1, color="#4daf4a", lw=0.7, ls="--", alpha=0.7)
    ax1.plot(flux, res_r, "o-", color="#377eb8", ms=4, lw=1, label="Ramsey")
    ax1.plot(flux, res_t, "s-", color="#e41a1c", ms=4, lw=1, mfc="none",
             label="Transient")
    ax1.set_xlabel(r"Flux bias $\Phi$ ($\Phi_0$)", fontsize=11)
    ax1.set_ylabel("Residual (MHz)", fontsize=10)
    ax1.legend(fontsize=8, ncol=2); ax1.grid(True, ls=":", alpha=0.4)

    fig.tight_layout()
    fig.savefig(os.path.join(out, "f_phi_curve_sqc.png"), dpi=300)
    fig.savefig(os.path.join(out, "f_phi_curve_sqc.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    recompute = "--recompute" in sys.argv
    d = compute(recompute=recompute)

    mae_r = np.mean(np.abs(d["ramsey"] - d["analytic"])) * 1e3
    mae_t = np.mean(np.abs(d["transient"] - d["analytic"])) * 1e3
    print(f"f(Φ): {len(d['offsets'])} points, "
          f"f range {d['analytic'].min():.4f}–{d['analytic'].max():.4f} GHz")
    print(f"MAE vs analytic:  ramsey={mae_r:.3f} MHz  transient={mae_t:.3f} MHz")
    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
