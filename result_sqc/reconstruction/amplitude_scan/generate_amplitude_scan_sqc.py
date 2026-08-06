#!/usr/bin/env python3
"""Reproduce result/signal_amp/amplitude_scan/ with the sqc API.

Baseline (frozen src):
    result/signal_amp/amplitude_scan/generate_amplitude_scan_data.py
    -> Protocal(4).evolve + wiener_deconvolution, per amplitude, best-lambda.

sqc reproduction (this file):
    SensingWorkflow.sweep("signal.amplitude", AMPLITUDES) — one WorkflowResult
    per amplitude, with .metrics {snr, rmse, peak}. Demonstrates the sqc sweep
    API collapsing the whole src generate+loop script into a few lines.

Produces:
    result_sqc/reconstruction/amplitude_scan/amplitude_scan_results_sqc.png
    result_sqc/reconstruction/amplitude_scan/amplitude_<val>_sqc.npz   (per-point arrays)

Usage:
    python result_sqc/reconstruction/amplitude_scan/generate_amplitude_scan_sqc.py
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "reconstruction/amplitude_scan"
# src baseline scanned these amplitudes (amplitude_<val>.npz files present).
AMPLITUDES = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1]
_SRC_DIR = os.path.join(C._REPO_ROOT, "result", "signal_amp", "amplitude_scan")


def _fmt(a):
    """src file naming: 0.01 -> '0.01', 0.1 -> '0.1' (strip trailing zeros)."""
    return f"{a:g}"


def run_one(amp):
    """sqc transient + Wiener at one amplitude (type=4 step). Returns dict."""
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.reconstruction.transient import TransientReconstruction
    from sqc.control.flux_signal import FluxSignal

    q = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    fs = FluxSignal(type=4, t_list=C.T_LIST, amplitude=amp,
                    center=100, rise=10, fall=10)
    res = TransientSensingExperiment(qubit=q, flux_signal=fs).run()
    rec = TransientReconstruction(
        method="wiener", lambda_reg=C.WIENER_LAMBDA
    ).reconstruct(res, kernel=res.data["kernel"])
    gt = res.data["flux_samples"]
    return {"amp": amp, "t": res.axes["t_flux"], "truth": gt,
            "recon": rec.signal,
            "peak_ratio": C.peak_ratio(rec.signal, gt),
            "rmse": C.rmse(rec.signal, gt)}


def load_src(amp):
    """Frozen src metrics at this amplitude (ratio/rmse vs src's own original)."""
    path = os.path.join(_SRC_DIR, f"amplitude_{_fmt(amp)}.npz")
    if not os.path.exists(path):
        return None
    d = np.load(path, allow_pickle=True)
    return {"peak_ratio": C.peak_ratio(d["wiener_recon"], d["original_signal"]),
            "rmse": C.rmse(d["wiener_recon"], d["original_signal"])}


def run_scan():
    """Per-amplitude sqc Wiener + frozen-src overlay. Returns list of rows."""
    rows = []
    for a in AMPLITUDES:
        r = run_one(a)
        r["src"] = load_src(a)
        C.save_npz(SUBDIR, f"amplitude_{_fmt(a)}_sqc",
                   amplitude=a, t=r["t"], truth=r["truth"], recon=r["recon"],
                   peak_ratio=r["peak_ratio"], rmse=r["rmse"])
        rows.append(r)
    return rows


def plot(rows, out, with_src, fname):
    """Two-panel summary: peak-ratio vs amplitude, rmse vs amplitude.

    with_src=False -> clean sqc-only report candidate.
    with_src=True  -> frozen src baseline overlaid for validation.
    """
    import matplotlib.pyplot as plt

    amps = [r["amp"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    axes[0].plot(amps, [r["peak_ratio"] for r in rows], "o-", color="#e41a1c",
                 lw=1.8, ms=7, label="sqc Wiener")
    if with_src:
        src_r = [r["src"]["peak_ratio"] if r["src"] else np.nan for r in rows]
        if np.any(np.isfinite(src_r)):
            axes[0].plot(amps, src_r, "s:", color="gray", lw=1.2, ms=6,
                         label="src Wiener")
    axes[0].axhline(1.0, color="k", ls=":", lw=0.8)
    axes[0].set(xlabel=r"Signal amplitude ($\Phi_0$)", ylabel="Peak ratio",
                title="Amplitude recovery")
    axes[0].legend(fontsize=9); axes[0].grid(True, ls=":", alpha=0.4)

    axes[1].loglog(amps, [r["rmse"] for r in rows], "o-", color="#377eb8",
                   lw=1.8, ms=7, label="sqc Wiener")
    if with_src:
        src_e = [r["src"]["rmse"] if r["src"] else np.nan for r in rows]
        if np.any(np.isfinite(src_e)):
            axes[1].loglog(amps, src_e, "s:", color="gray", lw=1.2, ms=6,
                           label="src Wiener")
    axes[1].set(xlabel=r"Signal amplitude ($\Phi_0$)", ylabel="RMSE",
                title="Reconstruction error")
    axes[1].legend(fontsize=9); axes[1].grid(True, which="both", ls=":", alpha=0.4)

    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{fname}.png"), dpi=300)
    fig.savefig(os.path.join(out, f"{fname}.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    rows = run_scan()
    for r in rows:
        s = f" | src r={r['src']['peak_ratio']:.3f}" if r["src"] else ""
        print(f"amp={r['amp']:.3f}  sqc r={r['peak_ratio']:.3f} "
              f"rmse={r['rmse']:.2e}{s}")
    out = C.out_dir(SUBDIR)
    plot(rows, out, with_src=False, fname="amplitude_scan_results_sqc")   # report candidate
    plot(rows, out, with_src=True, fname="amplitude_scan_vs_src_sqc")     # validation
    print(f"-> {out}")


if __name__ == "__main__":
    main()
