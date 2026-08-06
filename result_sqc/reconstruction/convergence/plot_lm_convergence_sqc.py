#!/usr/bin/env python3
"""Reproduce result/convergence/ with the sqc API.

Baseline (frozen src):
    result/convergence/plot_lm_convergence.py
    -> numerical_inverse returns history{res, b, mu}; plots convergence,
       residual heatmap, coefficient-norm evolution.

sqc reproduction (this file):
    TransientReconstruction(method="lm", ...).reconstruct(...) returns
    (FluxSignal, history) where history has:
      - "res": per-iteration residual VECTORS (p_meas - p_sim, length N_meas)
      - "b":   per-iteration coefficient vectors (length n_basis)
      - "mu":  per-iteration LM damping
    Same three figures as the src baseline.

The LM run is ~344 s (full-rho inversion), so the history is cached to
lm_history_sqc.npz on first run; plotting reads the cache. Pass --recompute
to force a fresh LM run.

Produces:
    result_sqc/reconstruction/convergence/lm_convergence_sqc.png
    result_sqc/reconstruction/convergence/lm_residual_heatmap_sqc.png
    result_sqc/reconstruction/convergence/lm_coefficient_norm_sqc.png
    result_sqc/reconstruction/convergence/lm_history_sqc.npz  (cache)

Usage:
    python result_sqc/reconstruction/convergence/plot_lm_convergence_sqc.py
    python result_sqc/reconstruction/convergence/plot_lm_convergence_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "reconstruction/convergence"
_CACHE = "lm_history_sqc"


def run_lm():
    """Run LM once via the sqc API. Returns (res_norms, res_matrix, b_matrix, mu).

    res_matrix: (n_iter, N_meas) residual vectors.
    b_matrix:   (n_iter+1, n_basis) coefficient vectors.
    """
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.control.flux_signal import FluxSignal

    q = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    fs = FluxSignal(type=3, t_list=C.T_LIST, amplitude=0.01,
                    rise=10, fall=10, center=100, noise_level=0.0)
    exp = TransientSensingExperiment(qubit=q, flux_signal=fs)
    res = exp.run()

    # LM requires the adapted result (p_e->p_meas, t_flux->t_signal); otherwise
    # reconstruct() raises KeyError: 'p_meas'. See _common.adapt_for_lm.
    lm = TransientReconstruction(
        method="lm", qubit=q, control_pulse=exp.control_pulse,
        basis_type="fourier", n_basis=C.LM_N_BASIS,
        lambda_reg=C.LM_LAMBDA, max_iter=C.LM_MAX_ITER)
    _, history = lm.reconstruct(C.adapt_for_lm(res), kernel=res.data["kernel"])

    res_matrix = np.asarray(history["res"], dtype=float)       # (n_iter, N_meas)
    b_matrix = np.asarray(history["b"], dtype=float)           # (n_iter+1, n_basis)
    mu = np.asarray(history["mu"], dtype=float)
    res_norms = np.linalg.norm(res_matrix, axis=1)
    return res_norms, res_matrix, b_matrix, mu


def get_history(recompute=False):
    """Load cached history, or run LM and cache it."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path)
        return d["res_norms"], d["res_matrix"], d["b_matrix"], d["mu"]
    res_norms, res_matrix, b_matrix, mu = run_lm()
    C.save_npz(SUBDIR, _CACHE, res_norms=res_norms, res_matrix=res_matrix,
               b_matrix=b_matrix, mu=mu)
    return res_norms, res_matrix, b_matrix, mu


def plot_convergence(res_norms, out):
    """Residual norm ||r|| vs iteration (log-y), src-style."""
    import matplotlib.pyplot as plt
    it = np.arange(1, len(res_norms) + 1)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(it, res_norms, "bo-", lw=2, ms=6, mfc="white", mew=1.5)
    ax.set(xlabel="Iteration Number", ylabel=r"Residual Norm $\|r\|$",
           title="LM Algorithm Convergence (sqc)")
    ax.set_yscale("log")
    ax.set_xlim(0, len(it) + 1)
    ax.grid(True, alpha=0.3, ls="--")
    ax.annotate(f"Final: {res_norms[-1]:.2e}", xy=(it[-1], res_norms[-1]),
                xytext=(it[-1] - 2, res_norms[-1] * 5),
                arrowprops=dict(arrowstyle="->", color="red", lw=1.5),
                fontsize=11, color="red")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "lm_convergence_sqc.png"), dpi=300)
    plt.close(fig)


def plot_residual_heatmap(res_matrix, out):
    """Residual vector magnitude across measurement index vs iteration."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5))
    im = ax.imshow(np.abs(res_matrix), aspect="auto", origin="lower",
                   cmap="viridis",
                   extent=[0, res_matrix.shape[1], 1, res_matrix.shape[0]])
    ax.set(xlabel="Measurement index", ylabel="Iteration",
           title="LM Residual |r| Heatmap (sqc)")
    fig.colorbar(im, ax=ax, label=r"$|r_i|$")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "lm_residual_heatmap_sqc.png"), dpi=300)
    plt.close(fig)


def plot_coefficient_norm(b_matrix, out):
    """||b|| evolution over iterations."""
    import matplotlib.pyplot as plt
    b_norms = np.linalg.norm(b_matrix, axis=1)
    it = np.arange(len(b_norms))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(it, b_norms, "rs-", lw=2, ms=6, mfc="white", mew=1.5)
    ax.set(xlabel="Iteration Number", ylabel=r"Coefficient Norm $\|b\|$",
           title="LM Coefficient-Norm Evolution (sqc)")
    ax.grid(True, alpha=0.3, ls="--")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "lm_coefficient_norm_sqc.png"), dpi=300)
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    recompute = "--recompute" in sys.argv
    res_norms, res_matrix, b_matrix, mu = get_history(recompute=recompute)
    out = C.out_dir(SUBDIR)

    plot_convergence(res_norms, out)
    plot_residual_heatmap(res_matrix, out)
    plot_coefficient_norm(b_matrix, out)
    print(f"{len(res_norms)} iters, final ||r||={res_norms[-1]:.3e} -> {out}")


if __name__ == "__main__":
    main()
