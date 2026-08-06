#!/usr/bin/env python3
"""Generate D1 with a superconducting-qubit Gaussian control envelope.

All three measurements and their response kernels are generated with adjacent
Gaussian pi/2 rotations.  The double-peak and complex cases also rerun LM with
the same control pulse; no square-pulse kernel or LM result is reused.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
RESULT_SQC = os.path.abspath(os.path.join(HERE, "..", ".."))
REPO_ROOT = os.path.abspath(os.path.join(RESULT_SQC, ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
if RESULT_SQC not in sys.path:
    sys.path.insert(0, RESULT_SQC)

import _common as C

STEM = "d1_waveform_regularization_sqc"
CONTROL_ENVELOPE = "gaussian"
CONTROL_ANGLE = np.pi / 2
GAUSSIAN_SIGMA_NS = 2.0
CACHE_VERSION = 2
SIGNAL_PARAMS = {
    "double_peak": (5, {"amplitude": 0.01, "center": 100, "width": 40}),
    "complex": (7, {"amplitude": 0.01}),
    "step": (4, {"amplitude": 0.01, "center": 100, "rise": 10, "fall": 10}),
}
LAMBDA_GRID = np.array(
    [0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0,
     15.0, 20.0, 30.0, 50.0],
    dtype=float,
)
DISPLAY_LAMBDAS = np.array([0.1, 5.0, 20.0, 50.0], dtype=float)
TOP_WIENER_LAMBDA = 10.0
BASELINE_LAMBDA = 5.0

COL_TRUTH = "#202020"
COL_WIENER = "#d92523"
COL_LM = "#2778b5"
COL_METRIC = "#2a8c62"
COL_NOTE = "#666666"


def load_npz(path: str) -> dict[str, np.ndarray]:
    with np.load(path) as data:
        return {key: data[key] for key in data.files}


def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def nrmse(signal: np.ndarray, truth: np.ndarray) -> float:
    scale = float(np.max(truth) - np.min(truth))
    return float(np.sqrt(np.mean((signal - truth) ** 2)) / scale)


def peak_ratio(signal: np.ndarray, truth: np.ndarray) -> float:
    return float(np.max(np.abs(signal)) / np.max(np.abs(truth)))


def cache_path(name: str) -> str:
    return os.path.join(HERE, f"d1_{name}_gaussian_cache_sqc.npz")


def generate_measurement_cache(name: str, with_lm: bool) -> dict[str, np.ndarray]:
    """Run dynamics, estimate the Gaussian-pulse kernel, and optionally LM."""
    from sqc.control.flux_signal import FluxSignal
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.reconstruction.transient import TransientReconstruction

    signal_type, params = SIGNAL_PARAMS[name]
    qubit = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    flux = FluxSignal(type=signal_type, t_list=C.T_LIST, **params)
    experiment = TransientSensingExperiment(
        qubit=qubit,
        flux_signal=flux,
        rotation_angle=CONTROL_ANGLE,
        envelope=CONTROL_ENVELOPE,
        envelope_sigma=GAUSSIAN_SIGMA_NS,
    )
    result = experiment.run()
    truth = np.asarray(result.data["flux_samples"], dtype=float)
    data = {
        "cache_version": np.array(CACHE_VERSION),
        "signal_name": np.array(name),
        "original_time": np.asarray(result.axes["t_flux"], dtype=float),
        "original_signal": truth,
        "scan_time": np.asarray(result.axes["scan"], dtype=float),
        "delta_p": np.asarray(result.data["delta_p"], dtype=float),
        "kernel_time": np.asarray(result.axes["t_samples"], dtype=float),
        "kernel": np.asarray(result.data["kernel"], dtype=float),
        "pulse_time": np.asarray(experiment.t_rabi, dtype=float),
        "pulse_envelope": np.asarray(
            experiment.control_pulse.pulses[0].Omega.signal, dtype=float
        ),
        "control_angle": np.array(result.metadata["rotation_angle"]),
        "envelope_sigma_ns": np.array(GAUSSIAN_SIGMA_NS),
    }

    if with_lm:
        lm = TransientReconstruction(
            method="lm",
            qubit=qubit,
            control_pulse=experiment.control_pulse,
            basis_type="fourier",
            n_basis=C.LM_N_BASIS,
            lambda_reg=C.LM_LAMBDA,
            max_iter=C.LM_MAX_ITER,
        )
        reconstructed, history = lm.reconstruct(
            C.adapt_for_lm(result), kernel=result.data["kernel"]
        )
        lm_signal = np.asarray(reconstructed.signal, dtype=float)
        data.update({
            "lm_time": np.asarray(result.axes["t_flux"], dtype=float),
            "lm_recon": lm_signal,
            "lm_ratio": np.array(peak_ratio(lm_signal, truth)),
            "lm_rmse": np.array(
                float(np.sqrt(np.mean((lm_signal - truth) ** 2)))
            ),
            "lm_res_history": np.asarray(history["res"], dtype=float),
        })

    np.savez(cache_path(name), **data)
    return data


def load_or_generate(name: str, with_lm: bool, force: bool) -> dict[str, np.ndarray]:
    path = cache_path(name)
    if not force and os.path.exists(path):
        cached = load_npz(path)
        version = int(np.asarray(cached.get("cache_version", -1)).item())
        sigma = float(np.asarray(cached.get("envelope_sigma_ns", np.nan)).item())
        if version == CACHE_VERSION and np.isclose(sigma, GAUSSIAN_SIGMA_NS):
            if not with_lm or "lm_recon" in cached:
                return cached
    return generate_measurement_cache(name, with_lm=with_lm)


def recompute_wiener(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    from sqc.reconstruction.transient import TransientReconstruction

    dt = float(data["scan_time"][1] - data["scan_time"][0])
    reconstructed = TransientReconstruction(
        method="wiener", lambda_reg=TOP_WIENER_LAMBDA
    ).reconstruct(data["delta_p"], kernel=data["kernel"], dt=dt)
    updated = dict(data)
    signal = np.asarray(reconstructed.signal, dtype=float)
    truth = np.asarray(data["original_signal"], dtype=float)
    updated["wiener_time"] = np.asarray(data["original_time"], dtype=float)
    updated["wiener_recon"] = signal
    updated["wiener_ratio"] = np.array(peak_ratio(signal, truth))
    updated["wiener_rmse"] = np.array(
        float(np.sqrt(np.mean((signal - truth) ** 2)))
    )
    return updated


def compute_step_scan(step: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    from sqc.reconstruction.transient import TransientReconstruction

    dt = float(step["scan_time"][1] - step["scan_time"][0])
    recons = []
    for lam in LAMBDA_GRID:
        reconstructor = TransientReconstruction(
            method="wiener", lambda_reg=float(lam)
        )
        reconstructed = reconstructor.reconstruct(
            step["delta_p"], kernel=step["kernel"], dt=dt
        )
        recons.append(np.asarray(reconstructed.signal, dtype=float))
    recons = np.asarray(recons)
    truth = np.asarray(step["original_signal"], dtype=float)
    return {
        "step_time": np.asarray(step["original_time"], dtype=float),
        "step_truth": truth,
        "step_lambdas": LAMBDA_GRID.copy(),
        "step_reconstructions": recons,
        "step_nrmse": np.array([nrmse(row, truth) for row in recons]),
        "step_peak_ratio": np.array([peak_ratio(row, truth) for row in recons]),
        "step_correlation": np.array(
            [np.corrcoef(row, truth)[0, 1] for row in recons]
        ),
        "step_dt_ns": np.array(dt),
    }


def plot_reconstruction(ax, data, title: str, panel: str) -> None:
    ax.plot(data["original_time"], data["original_signal"], color=COL_TRUTH,
            ls="--", lw=1.6, label="Truth")
    ax.plot(data["wiener_time"], data["wiener_recon"], color=COL_WIENER,
            lw=1.5,
            label=rf"Wiener ($\lambda_W={TOP_WIENER_LAMBDA:g}$)")
    ax.plot(data["lm_time"], data["lm_recon"], color=COL_LM,
            lw=1.5, label=r"LM ($\lambda_{LM}=100$)")
    ax.text(0.97, 0.06,
            f"W r={float(data['wiener_ratio']):.2f}\n"
            f"LM r={float(data['lm_ratio']):.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.7,
            color="#333333",
            bbox=dict(fc="white", ec="none", alpha=0.78, pad=1.5))
    ax.set_title(f"({panel}) {title}", fontsize=9.5)
    ax.set(xlim=(0, 200), xlabel="Time (ns)", ylabel=r"Flux ($\Phi_0$)")
    ax.tick_params(labelsize=7.5)
    ax.xaxis.label.set_size(8.5)
    ax.yaxis.label.set_size(8.5)
    ax.grid(True, ls=":", alpha=0.35, lw=0.5)


def plot_figure(double_peak, complex_signal, scan) -> tuple[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    plt.style.use("seaborn-v0_8-whitegrid")
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.05],
                  hspace=0.38, wspace=0.28)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    plot_reconstruction(ax_a, double_peak, "Double-peak reconstruction", "a")
    plot_reconstruction(ax_b, complex_signal,
                        "Complex-waveform reconstruction", "b")
    handles, labels = ax_a.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.95),
               ncol=3, frameon=False, fontsize=8, handlelength=2.6,
               columnspacing=1.8, handletextpad=0.7)

    ax_c = fig.add_subplot(gs[1, 0])
    t = scan["step_time"]
    truth = scan["step_truth"]
    ax_c.plot(t, truth, color=COL_TRUTH, ls="--", lw=1.7, label="Truth")
    styles = [
        ("#6a3d9a", ":"),
        ("#2778b5", "-"),
        ("#238b57", "--"),
        ("#e76f2e", "-."),
    ]
    for lam, (color, linestyle) in zip(DISPLAY_LAMBDAS, styles):
        idx = int(np.argmin(np.abs(scan["step_lambdas"] - lam)))
        ax_c.plot(t, scan["step_reconstructions"][idx], color=color,
                  ls=linestyle, lw=1.45, label=rf"$\lambda_W={lam:g}$")
    ax_c.set_title(r"(c) Step-like signal: Wiener $\lambda_W$ scan",
                   fontsize=9.5)
    ax_c.set(xlim=(45, 165), xlabel="Time (ns)",
             ylabel=r"Flux ($\Phi_0$)")
    ax_c.tick_params(labelsize=7.5)
    ax_c.xaxis.label.set_size(8.5)
    ax_c.yaxis.label.set_size(8.5)
    legend_c = ax_c.legend(fontsize=6.5, ncol=2, loc="lower left",
                           frameon=True, fancybox=False, framealpha=0.92,
                           borderpad=0.35, columnspacing=1.0,
                           handlelength=2.2)
    legend_c.get_frame().set_facecolor("white")
    legend_c.get_frame().set_edgecolor("#bbbbbb")
    legend_c.get_frame().set_linewidth(0.5)
    ax_c.grid(True, ls=":", alpha=0.35, lw=0.5)

    sub = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 1], hspace=0.12)
    ax_d1 = fig.add_subplot(sub[0])
    ax_d2 = fig.add_subplot(sub[1], sharex=ax_d1)
    lambdas = scan["step_lambdas"]
    errors = scan["step_nrmse"]
    ratios = scan["step_peak_ratio"]
    best_idx = int(np.argmin(errors))
    best_lambda = float(lambdas[best_idx])

    ax_d1.semilogx(lambdas, errors, "o-", color=COL_WIENER, lw=1.4, ms=3.5)
    ax_d1.axvline(BASELINE_LAMBDA, color=COL_NOTE, ls=":", lw=1.0)
    ax_d1.axvline(best_lambda, color=COL_METRIC, ls="--", lw=1.0)
    ax_d1.plot(best_lambda, errors[best_idx], marker="*", ms=7,
               color=COL_METRIC, zorder=4)
    ax_d1.annotate(rf"min at $\lambda_W={best_lambda:g}$",
                   xy=(best_lambda, errors[best_idx]), xytext=(-7, 7),
                   textcoords="offset points", fontsize=6.5, color=COL_METRIC,
                   ha="right")
    ax_d1.set_title("(d) Regularization-selection metrics", fontsize=9.5)
    ax_d1.set_ylabel("NRMSE", fontsize=8)
    ax_d1.tick_params(labelsize=7, labelbottom=False)
    ax_d1.grid(True, which="both", ls=":", alpha=0.35, lw=0.5)

    ax_d2.semilogx(lambdas, ratios, "s-", color=COL_METRIC, lw=1.4, ms=3.5)
    ax_d2.axhline(1.0, color=COL_TRUTH, ls="--", lw=0.9)
    ax_d2.axvline(BASELINE_LAMBDA, color=COL_NOTE, ls=":", lw=1.0,
                  label=rf"baseline $\lambda_W={BASELINE_LAMBDA:g}$")
    ax_d2.axvline(best_lambda, color=COL_METRIC, ls="--", lw=1.0,
                  label=rf"min NRMSE $\lambda_W={best_lambda:g}$")
    ax_d2.set_xlabel(r"Regularization $\lambda_W$", fontsize=8.5)
    ax_d2.set_ylabel("Peak ratio", fontsize=8)
    ax_d2.tick_params(labelsize=7)
    legend_d = ax_d2.legend(fontsize=5.9, frameon=True, fancybox=False,
                            framealpha=0.90, loc="lower left", ncol=1,
                            borderpad=0.3, handlelength=2.2)
    legend_d.get_frame().set_facecolor("white")
    legend_d.get_frame().set_edgecolor("#bbbbbb")
    legend_d.get_frame().set_linewidth(0.5)
    ax_d2.grid(True, which="both", ls=":", alpha=0.35, lw=0.5)

    fig.suptitle("Waveform reconstruction and regularization dependence",
                 fontsize=11, y=0.99)
    fig.text(0.98, 0.975,
             rf"Gaussian $\pi/2$ control, $\sigma={GAUSSIAN_SIGMA_NS:g}$ ns",
             ha="right", va="top",
             fontsize=6.5, color=COL_NOTE)
    fig.subplots_adjust(top=0.885, bottom=0.09, left=0.09, right=0.98)
    png = os.path.join(HERE, STEM + ".png")
    pdf = os.path.join(HERE, STEM + ".pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def save_outputs(double_peak, complex_signal, scan) -> tuple[str, str]:
    npz_path = os.path.join(HERE, STEM + ".npz")
    np.savez(
        npz_path,
        double_peak_time=double_peak["original_time"],
        double_peak_truth=double_peak["original_signal"],
        double_peak_wiener=double_peak["wiener_recon"],
        double_peak_lm=double_peak["lm_recon"],
        complex_time=complex_signal["original_time"],
        complex_truth=complex_signal["original_signal"],
        complex_wiener=complex_signal["wiener_recon"],
        complex_lm=complex_signal["lm_recon"],
        display_lambdas=DISPLAY_LAMBDAS,
        top_wiener_lambda=np.array(TOP_WIENER_LAMBDA),
        baseline_lambda=np.array(BASELINE_LAMBDA),
        **scan,
    )
    best_idx = int(np.argmin(scan["step_nrmse"]))
    metadata = {
        "figure": "D1",
        "status": "formal report candidate",
        "control": {
            "envelope": CONTROL_ENVELOPE,
            "rotation_angle_rad": CONTROL_ANGLE,
            "gaussian_sigma_ns": GAUSSIAN_SIGMA_NS,
            "pulse_window_ns": float(
                double_peak["pulse_time"][-1] - double_peak["pulse_time"][0]
            ),
            "kernel_recomputed_for_each_signal": True,
        },
        "source_caches": {
            name: os.path.relpath(cache_path(name), REPO_ROOT)
            for name in ("double_peak", "complex", "step")
        },
        "lambda_grid": LAMBDA_GRID.tolist(),
        "display_lambdas": DISPLAY_LAMBDAS.tolist(),
        "top_wiener_lambda": TOP_WIENER_LAMBDA,
        "baseline_lambda": BASELINE_LAMBDA,
        "step_min_nrmse_lambda": float(scan["step_lambdas"][best_idx]),
        "step_min_nrmse": float(scan["step_nrmse"][best_idx]),
        "notes": [
            "All quantum dynamics and kernels use Gaussian pi/2 control pulses.",
            "Top-row Wiener reconstructions are recomputed with lambda_W=10.",
            "The step scan reuses one Gaussian-control noiseless measurement.",
            "The step-specific NRMSE minimum is not claimed as a global optimum.",
            "Top-row LM results are rerun with the same Gaussian control pulse.",
        ],
    }
    meta_path = os.path.join(HERE, STEM + "_metadata.json")
    with open(meta_path, "w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2, ensure_ascii=True)

    metrics_path = os.path.join(HERE, STEM + "_metrics.txt")
    lines = [
        "D1 waveform reconstruction and regularization metrics",
        "=" * 72,
        f"control envelope = {CONTROL_ENVELOPE}",
        f"control rotation angle = {CONTROL_ANGLE:.10f} rad (pi/2)",
        f"Gaussian sigma = {GAUSSIAN_SIGMA_NS:g} ns",
        f"kernel points = {len(double_peak['kernel'])}",
        f"top-row Wiener lambda_W = {TOP_WIENER_LAMBDA:g}",
        f"double-peak NRMSE = "
        f"{nrmse(double_peak['wiener_recon'], double_peak['original_signal']):.5e}",
        f"double-peak peak ratio = {float(double_peak['wiener_ratio']):.5f}",
        f"double-peak LM NRMSE = "
        f"{nrmse(double_peak['lm_recon'], double_peak['original_signal']):.5e}",
        f"double-peak LM peak ratio = {float(double_peak['lm_ratio']):.5f}",
        f"complex-waveform NRMSE = "
        f"{nrmse(complex_signal['wiener_recon'], complex_signal['original_signal']):.5e}",
        f"complex-waveform peak ratio = {float(complex_signal['wiener_ratio']):.5f}",
        f"complex-waveform LM NRMSE = "
        f"{nrmse(complex_signal['lm_recon'], complex_signal['original_signal']):.5e}",
        f"complex-waveform LM peak ratio = {float(complex_signal['lm_ratio']):.5f}",
        "",
        "Step-like waveform regularization scan",
        "lambda_W      NRMSE       peak_ratio   correlation",
    ]
    for lam, err, ratio, corr in zip(
        scan["step_lambdas"], scan["step_nrmse"],
        scan["step_peak_ratio"], scan["step_correlation"]
    ):
        lines.append(f"{lam:8.3g}  {err:11.5e}  {ratio:11.5f}  {corr:11.7f}")
    lines.extend([
        "",
        f"baseline lambda_W = {BASELINE_LAMBDA:g}",
        f"step-specific minimum NRMSE lambda_W = "
        f"{float(scan['step_lambdas'][best_idx]):g}",
        "The minimum is specific to this cached noiseless step-like signal.",
    ])
    with open(metrics_path, "w", encoding="utf-8") as stream:
        stream.write("\n".join(lines) + "\n")
    return npz_path, metrics_path


def main() -> None:
    os.makedirs(HERE, exist_ok=True)
    force = "--force" in sys.argv
    double_peak = recompute_wiener(
        load_or_generate("double_peak", with_lm=True, force=force)
    )
    complex_signal = recompute_wiener(
        load_or_generate("complex", with_lm=True, force=force)
    )
    step = load_or_generate("step", with_lm=False, force=force)
    scan = compute_step_scan(step)
    npz_path, metrics_path = save_outputs(double_peak, complex_signal, scan)
    png, pdf = plot_figure(double_peak, complex_signal, scan)
    print(f"figure  -> {png}")
    print(f"        -> {pdf}")
    print(f"data    -> {npz_path}")
    print(f"metrics -> {metrics_path}")
    print(f"png sha256 -> {sha256(png)}")


if __name__ == "__main__":
    main()
