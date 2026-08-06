#!/usr/bin/env python3
"""Reproduce result/lambda_scan/ with the sqc API.

Baseline (frozen src):
    result/lambda_scan/plot_lambda_scan.py
    -> one Protocal(4).evolve, then wiener_deconvolution / numerical_inverse
       over a list of regularization lambdas.

sqc reproduction (this file):
    Run TransientSensingExperiment once, then reconstruct at each lambda:
      - Wiener:  TransientReconstruction(method="wiener", lambda_reg=lam)
      - LM:      TransientReconstruction(method="lm", ..., lambda_reg=lam)
    (Equivalently SensingWorkflow.sweep("reconstruction.lambda_reg", LAMBDAS).)

Produces:
    result_sqc/reconstruction/lambda_scan/wiener_lambda_scan_sqc.png
    result_sqc/reconstruction/lambda_scan/lm_lambda_scan_sqc.png

Usage:
    python result_sqc/reconstruction/lambda_scan/plot_lambda_scan_sqc.py
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "reconstruction/lambda_scan"
# Match the frozen src baseline lambda grids (result/lambda_scan/lambda_scan_data.npz).
WIENER_LAMBDAS = [2.0, 3.0, 4.0, 5.0, 10.0, 50.0, 100.0]
LM_LAMBDAS = [0.1, 1.0, 10.0, 100.0]
_SRC = os.path.join(C._REPO_ROOT, "result", "lambda_scan", "lambda_scan_data.npz")
_CACHE = "lambda_scan_sqc"


def run_once():
    """Single transient measurement reused across all lambdas."""
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.control.flux_signal import FluxSignal

    q = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    fs = FluxSignal(type=3, t_list=C.T_LIST, amplitude=0.01,
                    rise=10, fall=10, center=100, noise_level=0.0)
    exp = TransientSensingExperiment(qubit=q, flux_signal=fs)
    return q, exp, exp.run()


def compute(recompute=False):
    """Wiener over WIENER_LAMBDAS + LM over LM_LAMBDAS. Cached to npz (LM slow)."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    from sqc.reconstruction.transient import TransientReconstruction

    q, exp, res = run_once()
    truth = res.data["flux_samples"]
    t = res.axes["t_flux"]

    wiener_recons = []
    for lam in WIENER_LAMBDAS:
        rec = TransientReconstruction(method="wiener", lambda_reg=lam
                                      ).reconstruct(res, kernel=res.data["kernel"])
        wiener_recons.append(rec.signal)

    lm_recons = []
    for lam in LM_LAMBDAS:  # each ~6 min (full-rho)
        lm = TransientReconstruction(
            method="lm", qubit=q, control_pulse=exp.control_pulse,
            basis_type="fourier", n_basis=C.LM_N_BASIS,
            lambda_reg=lam, max_iter=C.LM_MAX_ITER)
        sig, _ = lm.reconstruct(C.adapt_for_lm(res), kernel=res.data["kernel"])
        lm_recons.append(sig.signal)

    data = {"t": t, "truth": truth,
            "wiener_lambdas": np.array(WIENER_LAMBDAS),
            "wiener_recons": np.array(wiener_recons),
            "lm_lambdas": np.array(LM_LAMBDAS),
            "lm_recons": np.array(lm_recons)}
    C.save_npz(SUBDIR, _CACHE, **data)
    return data


def _load_src():
    if not os.path.exists(_SRC):
        return None
    return dict(np.load(_SRC, allow_pickle=True))


def _scan_panel(t, truth, lambdas, recons, src_lam, src_recon, title, fname, out,
                match_lambda):
    """One lambda-scan panel. If src_lam is None -> sqc-only report candidate;
    otherwise the frozen src recon at the lambda closest to `match_lambda` is
    overlaid for validation."""
    import matplotlib.pyplot as plt
    from matplotlib import cm

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t, truth, "k--", lw=2, label="Truth", zorder=5)
    colors = cm.viridis(np.linspace(0, 0.9, len(lambdas)))
    for lam, rec, c in zip(lambdas, recons, colors):
        ax.plot(t, rec, color=c, lw=1.5, label=f"λ={float(lam):g}")
    if src_lam is not None:
        j = int(np.argmin(np.abs(np.asarray(src_lam) - match_lambda)))
        ax.plot(t, src_recon[j], "r:", lw=1.2, alpha=0.7,
                label=f"src λ={float(src_lam[j]):g}")
    ax.set(xlabel="Time (ns)", ylabel="Field (arb.)", title=title)
    ax.legend(fontsize=8, ncol=2); ax.grid(True, ls=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{fname}.png"), dpi=300)
    fig.savefig(os.path.join(out, f"{fname}.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    recompute = "--recompute" in sys.argv
    d = compute(recompute=recompute)
    src = _load_src()
    out = C.out_dir(SUBDIR)

    sw_lam = None if src is None else src["wiener_lambdas"]
    sw_rec = None if src is None else src["wiener_reconstructions"]
    slm_lam = None if src is None else src["lm_lambdas"]
    slm_rec = None if src is None else src["lm_reconstructions"]

    # --- Wiener: sqc-only report candidate + vs-src validation ---
    _scan_panel(d["t"], d["truth"], d["wiener_lambdas"], d["wiener_recons"],
                None, None, "Wiener λ scan (sqc)",
                "wiener_lambda_scan_sqc", out, C.WIENER_LAMBDA)
    _scan_panel(d["t"], d["truth"], d["wiener_lambdas"], d["wiener_recons"],
                sw_lam, sw_rec, "Wiener λ scan (sqc vs src)",
                "wiener_lambda_scan_vs_src_sqc", out, C.WIENER_LAMBDA)

    # --- LM: sqc-only report candidate + vs-src validation (match λ=100) ---
    _scan_panel(d["t"], d["truth"], d["lm_lambdas"], d["lm_recons"],
                None, None, "LM λ scan (sqc)",
                "lm_lambda_scan_sqc", out, C.LM_LAMBDA)
    _scan_panel(d["t"], d["truth"], d["lm_lambdas"], d["lm_recons"],
                slm_lam, slm_rec, "LM λ scan (sqc vs src)",
                "lm_lambda_scan_vs_src_sqc", out, C.LM_LAMBDA)
    print(f"wiener {len(d['wiener_lambdas'])} λ, LM {len(d['lm_lambdas'])} λ -> {out}")


if __name__ == "__main__":
    main()
