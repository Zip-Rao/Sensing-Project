#!/usr/bin/env python3
"""Reproduce result/basis_comparison/ with the sqc API.

Baseline (frozen src):
    result/basis_comparison/generate_basis_comparison_data.py
    -> Protocal(4).evolve, then numerical_inverse over basis types
       {bspline, fourier, legendre} at n_basis=20, lambdas=100.

sqc reproduction (this file):
    One TransientSensingExperiment.run(), then for each basis type:
      TransientReconstruction(method="lm", basis_type=<b>, n_basis=20,
          lambda_reg=100, control_pulse=exp.control_pulse, qubit=q)
    sqc supports the same three bases (sqc/reconstruction/basis.py).

Produces:
    result_sqc/reconstruction/basis_comparison/basis_comparison_sqc.png
    result_sqc/reconstruction/basis_comparison/basis_comparison_sqc.npz  (per-basis recon+rmse)

Usage:
    python result_sqc/reconstruction/basis_comparison/generate_basis_comparison_sqc.py
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "reconstruction/basis_comparison"
BASIS_TYPES = ["bspline", "fourier", "legendre"]
N_BASIS = 20
LAMBDA = 100.0
MAX_ITER = 30  # match src generate_basis_comparison_data.py
_CACHE = "basis_comparison_sqc"

# NOTE: result/basis_comparison/ ships only scripts (no committed npz), so this
# section is sqc-only — no frozen src overlay available.


def compute(recompute=False):
    """One transient run, then LM per basis type. Cached (each LM ~few min)."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.control.flux_signal import FluxSignal

    q = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    fs = FluxSignal(type=3, t_list=C.T_LIST, amplitude=0.01,
                    rise=10, fall=10, center=100, noise_level=0.0)
    exp = TransientSensingExperiment(qubit=q, flux_signal=fs)
    res = exp.run()
    truth = res.data["flux_samples"]
    t = res.axes["t_flux"]

    recons, ratios, rmses = [], [], []
    for b in BASIS_TYPES:
        # LM requires adapt_for_lm (p_e->p_meas); see _common.adapt_for_lm.
        lm = TransientReconstruction(
            method="lm", qubit=q, control_pulse=exp.control_pulse,
            basis_type=b, n_basis=N_BASIS, lambda_reg=LAMBDA, max_iter=MAX_ITER)
        sig, _ = lm.reconstruct(C.adapt_for_lm(res), kernel=res.data["kernel"])
        recons.append(sig.signal)
        ratios.append(C.peak_ratio(sig.signal, truth))
        rmses.append(C.rmse(sig.signal, truth))

    data = {"t": t, "truth": truth, "basis_types": np.array(BASIS_TYPES),
            "recons": np.array(recons), "ratios": np.array(ratios),
            "rmses": np.array(rmses)}
    C.save_npz(SUBDIR, _CACHE, **data)
    return data


def plot(d, out):
    import matplotlib.pyplot as plt

    t, truth = d["t"], d["truth"]
    bases = [str(b) for b in d["basis_types"]]
    colors = ["#e41a1c", "#377eb8", "#4daf4a"]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t, truth, "k--", lw=2, label="Truth", zorder=5)
    for b, rec, c, r, e in zip(bases, d["recons"], colors, d["ratios"], d["rmses"]):
        ax.plot(t, rec, color=c, lw=1.6,
                label=f"{b} (r={float(r):.2f}, rmse={float(e):.1e})")
    ax.set(xlabel="Time (ns)", ylabel="Field (arb.)",
           title=f"LM basis comparison (sqc, n_basis={N_BASIS}, λ={LAMBDA:g})")
    ax.legend(fontsize=9); ax.grid(True, ls=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "basis_comparison_sqc.png"), dpi=300)
    fig.savefig(os.path.join(out, "basis_comparison_sqc.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    recompute = "--recompute" in sys.argv
    d = compute(recompute=recompute)
    for b, r, e in zip(d["basis_types"], d["ratios"], d["rmses"]):
        print(f"{str(b):10s} peak-ratio={float(r):.3f}  rmse={float(e):.2e}")
    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
