#!/usr/bin/env python3
"""Reproduce result/different_signals/ with the sqc API.

Baseline (frozen src):
    result/different_signals/generate_single_signal_data.py
    -> Protocal(4).evolve + Analysis.wiener_deconvolution(lambdas=5)
       + Analysis.numerical_inverse(fourier, n_basis=100, lambdas=100, 20 it)
    across signal types 2/sine, 4/step, 5/double_peak, 7/complex.

sqc reproduction (this file):
    TransientSensingExperiment(qubit@0.9553, FluxSignal(type=...)).run()
    -> TransientReconstruction(method="wiener", lambda_reg=5)
    -> TransientReconstruction(method="lm", control_pulse=..., qubit=...,
                               basis_type="fourier", n_basis=100,
                               lambda_reg=100, max_iter=20)

Produces per signal:
    result_sqc/reconstruction/different_signals/signal_<name>_sqc.npz   (all arrays + src overlay)
    result_sqc/reconstruction/different_signals/check_<name>_sqc.png    (quick-look)
Panel across all four is built by plot_signal_comparison_panel_sqc.py.

Usage:
    python result_sqc/reconstruction/different_signals/generate_signal_data_sqc.py sine
    python result_sqc/reconstruction/different_signals/generate_signal_data_sqc.py all
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "reconstruction/different_signals"

# src default signal params per type (from generate_single_signal_data.py)
SIGNAL_PARAMS = {
    2: {"amplitude": 0.01, "frequency": 0.01},
    4: {"amplitude": 0.01, "center": 100, "rise": 10, "fall": 10},
    5: {"amplitude": 0.01, "center": 100, "width": 40},
    7: {"amplitude": 0.01},
}

# Frozen src baseline (R1 ground truth). We LOAD these for overlay rather than
# re-running src live (~14 min/signal) — see memory result-sqc-src-overlay-frozen.
_SRC_DIR = os.path.join(C._REPO_ROOT, "result", "different_signals")


def load_src_baseline(name):
    """Load frozen src arrays for overlay. Returns dict or None if absent."""
    path = os.path.join(_SRC_DIR, f"signal_{name}.npz")
    if not os.path.exists(path):
        return None
    d = np.load(path, allow_pickle=True)
    return {
        "src_original_time": d["original_time"], "src_original": d["original_signal"],
        "src_wiener_time": d["wiener_time"], "src_wiener": d["wiener_recon"],
        "src_lm_time": d["lm_time"], "src_lm": d["lm_recon"],
    }


def run_one(signal_type, with_lm=True):
    """Run sqc sensing + wiener (+ lm) for one signal type. Returns data dict."""
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.control.flux_signal import FluxSignal

    name = C.SIGNAL_NAMES[signal_type]
    params = SIGNAL_PARAMS[signal_type]
    q = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    fs = FluxSignal(type=signal_type, t_list=C.T_LIST, **params)

    exp = TransientSensingExperiment(qubit=q, flux_signal=fs)
    res = exp.run()
    orig = res.data["flux_samples"]

    # --- Wiener (lambda=5) ---
    wiener = TransientReconstruction(method="wiener", lambda_reg=C.WIENER_LAMBDA)
    rec_w = wiener.reconstruct(res, kernel=res.data["kernel"])

    data = {
        "signal_type": signal_type,
        "original_time": res.axes["t_flux"], "original_signal": orig,
        "scan_time": res.axes["scan"], "delta_p": res.data["delta_p"],
        "kernel": res.data["kernel"],
        "wiener_time": res.axes["t_flux"], "wiener_recon": rec_w.signal,
        "wiener_ratio": C.peak_ratio(rec_w.signal, orig),
        "wiener_rmse": C.rmse(rec_w.signal, orig),
    }

    # --- LM (fourier, n_basis=100, lambda=100, 20 it) — full-rho, ~6 min ---
    # LM needs the adapted result (p_e->p_meas, t_flux->t_signal); returns
    # (FluxSignal, history). See _common.adapt_for_lm / lm-transient-p-meas-gap.
    if with_lm:
        lm = TransientReconstruction(
            method="lm", qubit=q, control_pulse=exp.control_pulse,
            basis_type="fourier", n_basis=C.LM_N_BASIS,
            lambda_reg=C.LM_LAMBDA, max_iter=C.LM_MAX_ITER)
        rec_lm, history = lm.reconstruct(C.adapt_for_lm(res), kernel=res.data["kernel"])
        data.update({
            "lm_time": res.axes["t_flux"], "lm_recon": rec_lm.signal,
            "lm_ratio": C.peak_ratio(rec_lm.signal, orig),
            "lm_rmse": C.rmse(rec_lm.signal, orig),
            "lm_res_history": np.asarray(history["res"], dtype=float),
        })

    # --- frozen src baseline overlay ---
    src = load_src_baseline(name)
    if src is not None:
        data.update(src)

    C.save_npz(SUBDIR, f"signal_{name}_sqc", **data)
    return data


def plot_check(data, name):
    """Quick-look 2-panel figure (original + reconstruction overlay).

    Mirrors src check_<name>.png: left = original signal, right = sqc
    wiener/LM reconstructions with the frozen src baseline overlaid.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = data["original_time"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].plot(t, data["original_signal"], "b-", label="Original (sqc)")
    if "src_original" in data:
        axes[0].plot(data["src_original_time"], data["src_original"],
                     "k:", alpha=0.6, label="Original (src)")
    axes[0].set(xlabel="Time (ns)", ylabel="Magnetic field",
                title=f"{name} Signal")
    axes[0].grid(True, alpha=0.3); axes[0].legend()

    ax = axes[1]
    ax.plot(t, data["original_signal"], "k-", alpha=0.7, label="Original")
    ax.plot(data["wiener_time"], data["wiener_recon"], "r--", lw=1.5,
            label=f"Wiener (sqc, r={data['wiener_ratio']:.2f})")
    if "src_wiener" in data:
        ax.plot(data["src_wiener_time"], data["src_wiener"], "r:", alpha=0.5,
                lw=1.0, label="Wiener (src)")
    if "lm_recon" in data:
        ax.plot(data["lm_time"], data["lm_recon"], "b-.", lw=1.5,
                label=f"LM (sqc, r={data['lm_ratio']:.2f})")
    if "src_lm" in data:
        ax.plot(data["src_lm_time"], data["src_lm"], "b:", alpha=0.5,
                lw=1.0, label="LM (src)")
    ax.set(xlabel="Time (ns)", ylabel="Magnetic field",
           title=f"{name} Reconstruction")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)

    plt.tight_layout()
    path = os.path.join(C.out_dir(SUBDIR), f"check_{name}_sqc.png")
    plt.savefig(path, dpi=150); plt.close(fig)
    return path


def _resolve_type(arg):
    """Map a CLI arg (name like 'sine' or int like '2') to a signal_type."""
    name_to_type = {v: k for k, v in C.SIGNAL_NAMES.items()}
    if arg in name_to_type:
        return name_to_type[arg]
    return int(arg)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    with_lm = "--no-lm" not in sys.argv
    arg = args[0] if args else "all"
    types = list(C.SIGNAL_NAMES) if arg == "all" else [_resolve_type(arg)]
    for st in types:
        name = C.SIGNAL_NAMES[st]
        d = run_one(st, with_lm=with_lm)
        lm_str = f", LM r={d['lm_ratio']:.3f}" if "lm_ratio" in d else ""
        print(f"[{name}] wiener r={d['wiener_ratio']:.3f} "
              f"rmse={d['wiener_rmse']:.2e}{lm_str}")
        print(f"    check -> {plot_check(d, name)}")


if __name__ == "__main__":
    main()
