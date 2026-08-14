#!/usr/bin/env python3
"""NEW figure (no src counterpart): work-point effect on reconstruction.

Motivation
----------
The frozen src scripts all sit at optimal_work_point()=0.9553 (linear regime).
This figure makes the work-point dependence explicit — the core constraint of
transient sensing — by reconstructing the SAME signal at three regimes:

    - sweet spot   (flux=0.0):   dω/dΦ = 0, kernel ill-conditioned -> fails
    - linear       (flux=0.9553): φ_max <~ 1 rad -> faithful (ratio ~1)
    - nonlinear    (flux=0.25):  φ_max >> 1, phase-wrap -> Wiener under-recovers

Criterion: φ_max = κ · A · t_free must stay <~ 1 rad for linear Wiener.

Produces:
    result_sqc/reconstruction/workpoint_comparison/workpoint_comparison_sqc.png
    result_sqc/reconstruction/workpoint_comparison/workpoint_metrics_sqc.npz

Usage:
    python result_sqc/reconstruction/workpoint_comparison/generate_workpoint_comparison_sqc.py
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "reconstruction/workpoint_comparison"

# (flux, label, expected regime) — the three representative snapshot panels.
WORK_POINTS = [
    (0.0,    "sweet spot",  "kappa=0, ill-conditioned"),
    (0.9553, "linear",      "phi_max<~1, faithful"),
    (0.25,   "nonlinear",   "phase-wrap, under-recovery"),
]

# Dense scan for the summary panel (d). f01(Φ) is even about the sweet spot,
# so scanning symmetric ±Φ makes κ=dω/dΦ change sign across Φ=0 — the summary
# curve shows BOTH the slope sign (via signed κ) and the |κ| dependence of
# recovery. NON-UNIFORM in flux: a fine cluster near the sweet spot resolves
# the κ→0 rise-to-a-peak (uniform-in-flux sampling leaves κ∈(0,5.6) empty and
# makes the sweet-spot RMSE maximum look like a needle artifact).
_FINE = np.array([0.008, 0.016, 0.025, 0.035])          # near-sweet-spot cluster
_COARSE = np.array([0.045, 0.09, 0.135, 0.18, 0.225,
                    0.27, 0.315, 0.36, 0.405, 0.45])
_HALF = np.concatenate([_FINE, _COARSE])
SCAN_FLUX = np.concatenate([-_HALF[::-1], [0.0], _HALF])  # symmetric, 29 pts
_CACHE = "workpoint_scan_sqc"


def signed_kappa(flux, h=1e-4):
    """Signed dω/dΦ (rad/ns per Φ₀) — sign encodes the dispersion slope."""
    qp = C.make_sqc_qubit(flux + h)
    qm = C.make_sqc_qubit(flux - h)
    return (qp.frequency - qm.frequency) / (2 * h)


def kappa(flux, h=1e-4):
    """|dω/dΦ| at a work point (sensitivity magnitude)."""
    return abs(signed_kappa(flux, h))


def _sense_wiener(flux):
    """Sense + Wiener-reconstruct the shared Gaussian at one flux. Returns
    (recon, truth, delta_p)."""
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.control.flux_signal import FluxSignal

    q = C.make_sqc_qubit(flux)
    fs = FluxSignal(type=3, t_list=C.T_LIST, amplitude=0.01,
                    rise=10, fall=10, center=100, noise_level=0.0)
    res = TransientSensingExperiment(qubit=q, flux_signal=fs).run()
    rec = TransientReconstruction(
        method="wiener", lambda_reg=C.WIENER_LAMBDA
    ).reconstruct(res, kernel=res.data["kernel"])
    return rec.signal, res.data["flux_samples"], res.data["delta_p"]


def run_one(flux):
    """Full snapshot at one work point (waveforms + metrics)."""
    recon, gt, delta_p = _sense_wiener(flux)
    return {
        "flux": flux, "kappa": kappa(flux), "signed_kappa": signed_kappa(flux),
        "delta_p": delta_p, "recon": recon, "truth": gt,
        "peak_ratio": C.peak_ratio(recon, gt), "rmse": C.rmse(recon, gt),
    }


def run_scan(recompute=False):
    """Dense flux scan for panel (d). Returns (flux, signed_kappa, peak_ratio,
    rmse) arrays. Cached — Wiener per point ~10 s, 21 pts ~3-4 min.

    NOTE: peak_ratio is NOT monotonic-in-quality at large |κ| — phase-wrap
    oscillations inflate the reconstructed peak (r>1 = overshoot, not recovery),
    so RMSE is the honest quality axis for the summary panel (clean valley,
    minimum at the optimal |κ|)."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path)
        if "rmse" in d.files:
            return d["flux"], d["signed_kappa"], d["peak_ratio"], d["rmse"]

    sk, pr, er = [], [], []
    for phi in SCAN_FLUX:
        recon, gt, _ = _sense_wiener(float(phi))
        sk.append(signed_kappa(float(phi)))
        pr.append(C.peak_ratio(recon, gt))
        er.append(C.rmse(recon, gt))
    sk, pr, er = np.array(sk), np.array(pr), np.array(er)
    C.save_npz(SUBDIR, _CACHE, flux=SCAN_FLUX, signed_kappa=sk,
               peak_ratio=pr, rmse=er)
    return SCAN_FLUX, sk, pr, er


def plot(rows, scan, out):
    """3 snapshot columns + a signed-κ recovery summary (M-curve).

    scan = (flux, signed_kappa, peak_ratio) dense arrays for panel (d).
    """
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(13, 4.2))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.15])
    for i, ((fb, label, note), r) in enumerate(zip(WORK_POINTS, rows)):
        ax = fig.add_subplot(gs[0, i])
        t = C.T_LIST
        ax.plot(t, r["truth"], "k--", lw=1.5, label="Truth")
        ax.plot(t, r["recon"], "#e41a1c", lw=1.8, label="Wiener (sqc)")
        ax.set_title(f"{label}\n" + rf"$\Phi$={fb:.3f}, $\kappa$={r['kappa']:.1f}",
                     fontsize=10)
        ax.set_xlabel("Time (ns)", fontsize=9)
        if i == 0:
            ax.set_ylabel("Field (arb.)", fontsize=10)
        ax.text(0.97, 0.05, f"r={r['peak_ratio']:.2f}", transform=ax.transAxes,
                fontsize=8, ha="right", va="bottom")
        ax.text(0.03, 0.97, f"({chr(97+i)})", transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")
        ax.grid(True, ls=":", alpha=0.4)
        if i == 0:
            ax.legend(fontsize=8, loc="upper right")

    # summary panel (d): RMSE vs SIGNED κ (dense scan) + named markers.
    # RMSE (not peak-ratio) is the honest quality axis — peak-ratio overshoots
    # (>1) at large |κ| due to phase-wrap oscillations. RMSE gives a clean
    # valley: ill-conditioned at κ=0, minimum at optimal |κ|, phase-wrap at large.
    _flux, sk, pr, er = scan
    order = np.argsort(sk)
    ax = fig.add_subplot(gs[0, 3])
    ax.semilogy(sk[order], er[order], "-", color="#377eb8", lw=1.6, zorder=2,
                label="dense scan")
    ax.axvline(0.0, color="gray", ls="--", lw=0.8)
    for r, (_, label, _) in zip(rows, WORK_POINTS):
        ax.semilogy(r["signed_kappa"], r["rmse"], "o", color="#e41a1c",
                    ms=8, zorder=3)
        ax.annotate(label, (r["signed_kappa"], r["rmse"]), fontsize=7,
                    textcoords="offset points", xytext=(5, 4))
    ax.set(xlabel=r"Signed sensitivity $\kappa=d\omega/d\Phi$",
           ylabel="RMSE", title="Error vs sensitivity")
    ax.text(0.03, 0.97, "(d)", transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")
    ax.legend(fontsize=7, loc="upper center"); ax.grid(True, which="both",
                                                       ls=":", alpha=0.4)

    fig.tight_layout()
    fig.savefig(os.path.join(out, "workpoint_comparison_sqc.png"), dpi=300)
    fig.savefig(os.path.join(out, "workpoint_comparison_sqc.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    recompute = "--recompute" in sys.argv
    rows = [run_one(fb) for fb, _, _ in WORK_POINTS]
    print("flux    kappa    peak_ratio   rmse")
    for (fb, label, _), r in zip(WORK_POINTS, rows):
        print(f"{fb:.4f}  {r['kappa']:6.2f}   {r['peak_ratio']:8.3f}   "
              f"{r['rmse']:.2e}   [{label}]")

    scan = run_scan(recompute=recompute)
    print(f"dense scan: {len(scan[0])} pts, signed κ "
          f"[{scan[1].min():.1f}, {scan[1].max():.1f}], "
          f"RMSE min={scan[3].min():.2e} at κ={scan[1][np.argmin(scan[3])]:.1f}")

    out = C.out_dir(SUBDIR)
    plot(rows, scan, out)
    C.save_npz(SUBDIR, "workpoint_metrics_sqc",
               fluxes=np.array([r["flux"] for r in rows]),
               kappas=np.array([r["kappa"] for r in rows]),
               signed_kappas=np.array([r["signed_kappa"] for r in rows]),
               peak_ratios=np.array([r["peak_ratio"] for r in rows]),
               rmses=np.array([r["rmse"] for r in rows]),
               truths=np.array([r["truth"] for r in rows]),
               recons=np.array([r["recon"] for r in rows]),
               scan_flux=scan[0], scan_signed_kappa=scan[1], scan_peak_ratio=scan[2])
    print(f"-> {out}")


if __name__ == "__main__":
    main()
