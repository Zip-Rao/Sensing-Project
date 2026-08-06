#!/usr/bin/env python3
"""Generate report figure D3 with Gaussian pi/2 sensing control.

All amplitude, work-point, Wiener, and LM results are generated with the same
truncated Gaussian control used by D1.  D3 owns and validates its caches, so
legacy square-pulse data cannot be mixed into the report figure.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import _common as C


SUBDIR = "reconstruction/D3_nonlinearity_boundary"
AMPLITUDES = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10])
LM_AMPLITUDES = np.array([0.01, 0.04, 0.08, 0.10])
PULSE_WINDOW_NS = 10.0
GAUSSIAN_SIGMA_NS = 2.0
ROTATION_ANGLE = np.pi / 2.0
AWG_DT_NS = 0.5
CACHE_VERSION = 2

_FINE = np.array([0.008, 0.016, 0.025, 0.035])
_COARSE = np.array([0.045, 0.09, 0.135, 0.18, 0.225,
                    0.27, 0.315, 0.36, 0.405, 0.45])
_HALF = np.concatenate([_FINE, _COARSE])
WORKPOINT_FLUX = np.concatenate([-_HALF[::-1], [0.0], _HALF])


def _fmt(value: float) -> str:
    return f"{value:g}"


def omega01(flux):
    """Transmon dispersion in the angular-frequency units used by sqc."""
    flux = np.asarray(flux, dtype=float)
    ej_eff = C.EJ * np.abs(np.cos(np.pi * flux))
    return np.sqrt(8.0 * ej_eff * C.EC) - C.EC


def signed_kappa(flux: float, step: float = 1e-6) -> float:
    return float((omega01(flux + step) - omega01(flux - step)) / (2.0 * step))


def _experiment(qubit, signal):
    """Transient experiment under the D1 Gaussian-control convention."""
    from sqc.experiments.transient import TransientSensingExperiment

    return TransientSensingExperiment(
        qubit=qubit,
        flux_signal=signal,
        t_rabi=np.arange(0.0, PULSE_WINDOW_NS, AWG_DT_NS),
        rotation_angle=ROTATION_ANGLE,
        envelope="gaussian",
        envelope_sigma=GAUSSIAN_SIGMA_NS,
    )


def _cache_matches(data) -> bool:
    """Reject legacy square-pulse or differently parameterised caches."""
    try:
        return (
            int(np.asarray(data["cache_version"]).item()) == CACHE_VERSION
            and str(np.asarray(data["control_envelope"]).item()) == "gaussian"
            and np.isclose(float(np.asarray(data["rotation_angle_rad"])),
                           ROTATION_ANGLE)
            and np.isclose(float(np.asarray(data["envelope_sigma_ns"])),
                           GAUSSIAN_SIGMA_NS)
            and np.isclose(float(np.asarray(data["pulse_window_ns"])),
                           PULSE_WINDOW_NS)
        )
    except (KeyError, ValueError, TypeError):
        return False


def run_wiener(amplitude: float, recompute: bool = False) -> dict:
    """Full nonlinear forward run plus Wiener reconstruction at one amplitude."""
    from sqc.control.flux_signal import FluxSignal
    from sqc.reconstruction.transient import TransientReconstruction

    out = C.out_dir(SUBDIR)
    path = os.path.join(out, f"gaussian_amplitude_{_fmt(amplitude)}_sqc.npz")
    if os.path.exists(path) and not recompute:
        data = np.load(path)
        if _cache_matches(data):
            return {key: data[key] for key in data.files}

    q = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    signal = FluxSignal(type=4, t_list=C.T_LIST, amplitude=float(amplitude),
                        center=100, rise=10, fall=10, noise_level=0.0)
    exp = _experiment(q, signal)
    result = exp.run()
    truth = np.asarray(result.data["flux_samples"], dtype=float)
    reconstructed = TransientReconstruction(
        method="wiener", lambda_reg=C.WIENER_LAMBDA
    ).reconstruct(result, kernel=result.data["kernel"])
    recon = np.asarray(reconstructed.signal, dtype=float)
    payload = {
        "cache_version": np.array(CACHE_VERSION),
        "control_envelope": np.array("gaussian"),
        "rotation_angle_rad": np.array(ROTATION_ANGLE),
        "envelope_sigma_ns": np.array(GAUSSIAN_SIGMA_NS),
        "pulse_window_ns": np.array(PULSE_WINDOW_NS),
        "amplitude": np.array(amplitude),
        "t": np.asarray(result.axes["t_flux"], dtype=float),
        "truth": truth,
        "recon": recon,
        "kernel": np.asarray(result.data["kernel"], dtype=float),
        "peak_ratio": np.array(C.peak_ratio(recon, truth)),
        "rmse": np.array(C.rmse(recon, truth)),
    }
    np.savez(path, **payload)
    return payload


def load_wiener_rows(recompute: bool = False) -> list[dict]:
    rows = []
    for amplitude in AMPLITUDES:
        data = run_wiener(float(amplitude), recompute=recompute)
        truth = np.asarray(data["truth"], dtype=float)
        recon = np.asarray(data["recon"], dtype=float)
        rows.append({
            "amplitude": float(amplitude),
            "t": np.asarray(data["t"], dtype=float),
            "truth": truth,
            "recon": recon,
            "peak_excursion": float(np.max(np.abs(truth))),
            "peak_ratio": C.peak_ratio(recon, truth),
            "rmse": C.rmse(recon, truth),
        })
    return rows


def linear_model_error(signal: np.ndarray, bias: float) -> float:
    """Relative L2 error of the local linear frequency-shift model."""
    exact = omega01(bias + signal) - omega01(bias)
    linear = signed_kappa(bias) * signal
    denominator = np.linalg.norm(exact)
    return float(np.linalg.norm(linear - exact) / denominator) if denominator > 1e-15 else np.nan


def load_workpoint_scan(recompute: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Dense work-point scan under the same Gaussian sensing control."""
    from sqc.control.flux_signal import FluxSignal
    from sqc.reconstruction.transient import TransientReconstruction

    path = os.path.join(C.out_dir(SUBDIR), "gaussian_workpoint_scan_sqc.npz")
    if os.path.exists(path) and not recompute:
        data = np.load(path)
        if _cache_matches(data):
            return data["flux"], data["signed_kappa"], data["rmse"]

    slopes, errors = [], []
    for i, flux in enumerate(WORKPOINT_FLUX):
        print(f"work point {i + 1:02d}/{len(WORKPOINT_FLUX)}: Phi={flux:.4f}",
              flush=True)
        q = C.make_sqc_qubit(float(flux))
        signal = FluxSignal(type=3, t_list=C.T_LIST, amplitude=0.01,
                            rise=10, fall=10, center=100, noise_level=0.0)
        result = _experiment(q, signal).run()
        reconstructed = TransientReconstruction(
            method="wiener", lambda_reg=C.WIENER_LAMBDA
        ).reconstruct(result, kernel=result.data["kernel"])
        truth = np.asarray(result.data["flux_samples"], dtype=float)
        slopes.append(signed_kappa(float(flux)))
        errors.append(C.rmse(reconstructed.signal, truth))

    payload = {
        "cache_version": np.array(CACHE_VERSION),
        "control_envelope": np.array("gaussian"),
        "rotation_angle_rad": np.array(ROTATION_ANGLE),
        "envelope_sigma_ns": np.array(GAUSSIAN_SIGMA_NS),
        "pulse_window_ns": np.array(PULSE_WINDOW_NS),
        "flux": WORKPOINT_FLUX,
        "signed_kappa": np.asarray(slopes, dtype=float),
        "rmse": np.asarray(errors, dtype=float),
    }
    np.savez(path, **payload)
    return payload["flux"], payload["signed_kappa"], payload["rmse"]


def run_lm(amplitude: float, recompute: bool = False) -> dict:
    """Run one LM point and cache success or failure without fallback."""
    out = C.out_dir(SUBDIR)
    path = os.path.join(out, f"lm_amplitude_{_fmt(amplitude)}_sqc.npz")
    if os.path.exists(path) and not recompute:
        data = np.load(path)
        if _cache_matches(data):
            return {key: data[key] for key in data.files}

    from sqc.control.flux_signal import FluxSignal
    from sqc.reconstruction.transient import TransientReconstruction

    started = time.perf_counter()
    q = C.make_sqc_qubit(C.OPTIMAL_FLUX)
    signal = FluxSignal(type=4, t_list=C.T_LIST, amplitude=float(amplitude),
                        center=100, rise=10, fall=10, noise_level=0.0)
    exp = _experiment(q, signal)

    try:
        result = exp.run()
        truth = np.asarray(result.data["flux_samples"], dtype=float)
        wiener = TransientReconstruction(
            method="wiener", lambda_reg=C.WIENER_LAMBDA
        ).reconstruct(result, kernel=result.data["kernel"])
        solver = TransientReconstruction(
            method="lm", qubit=q, control_pulse=exp.control_pulse,
            basis_type="fourier", n_basis=C.LM_N_BASIS,
            lambda_reg=C.LM_LAMBDA, max_iter=C.LM_MAX_ITER,
        )
        reconstructed, history = solver.reconstruct(
            C.adapt_for_lm(result), kernel=result.data["kernel"],
            initial_guess=np.asarray(wiener.signal, dtype=float),
        )
        lm_signal = np.asarray(reconstructed.signal, dtype=float)
        payload = {
            "cache_version": np.array(CACHE_VERSION),
            "control_envelope": np.array("gaussian"),
            "rotation_angle_rad": np.array(ROTATION_ANGLE),
            "envelope_sigma_ns": np.array(GAUSSIAN_SIGMA_NS),
            "pulse_window_ns": np.array(PULSE_WINDOW_NS),
            "amplitude": np.array(amplitude),
            "t": np.asarray(result.axes["t_flux"], dtype=float),
            "truth": truth,
            "lm_recon": lm_signal,
            "lm_rmse": np.array(C.rmse(lm_signal, truth)),
            "lm_peak_ratio": np.array(C.peak_ratio(lm_signal, truth)),
            "residual_norm": np.asarray([
                np.linalg.norm(item) for item in history.get("res", [])
            ], dtype=float),
            "runtime_s": np.array(time.perf_counter() - started),
            "status": np.array("success"),
            "error": np.array(""),
        }
    except Exception as exc:
        payload = {
            "cache_version": np.array(CACHE_VERSION),
            "control_envelope": np.array("gaussian"),
            "rotation_angle_rad": np.array(ROTATION_ANGLE),
            "envelope_sigma_ns": np.array(GAUSSIAN_SIGMA_NS),
            "pulse_window_ns": np.array(PULSE_WINDOW_NS),
            "amplitude": np.array(amplitude),
            "t": C.T_LIST.copy(),
            "truth": np.full_like(C.T_LIST, np.nan, dtype=float),
            "lm_recon": np.full_like(C.T_LIST, np.nan, dtype=float),
            "lm_rmse": np.array(np.nan),
            "lm_peak_ratio": np.array(np.nan),
            "residual_norm": np.array([], dtype=float),
            "runtime_s": np.array(time.perf_counter() - started),
            "status": np.array("failed"),
            "error": np.array(f"{type(exc).__name__}: {exc}"),
        }
    np.savez(path, **payload)
    return payload


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=C._REPO_ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def save_outputs(wiener_rows, lm_rows, workpoint_scan, model_errors):
    out = C.out_dir(SUBDIR)
    flux_scan, kappa_scan, workpoint_rmse = workpoint_scan
    lm_status = np.array([str(np.asarray(row["status"]).item()) for row in lm_rows])
    lm_error = np.array([str(np.asarray(row["error"]).item()) for row in lm_rows])
    lm_rmse = np.array([float(np.asarray(row["lm_rmse"])) for row in lm_rows])
    lm_peak_ratio = np.array([float(np.asarray(row["lm_peak_ratio"])) for row in lm_rows])
    lm_runtime = np.array([float(np.asarray(row["runtime_s"])) for row in lm_rows])

    np.savez(
        os.path.join(out, "nonlinearity_boundary_sqc.npz"),
        cache_version=np.array(CACHE_VERSION),
        control_envelope=np.array("gaussian"),
        rotation_angle_rad=np.array(ROTATION_ANGLE),
        envelope_sigma_ns=np.array(GAUSSIAN_SIGMA_NS),
        pulse_window_ns=np.array(PULSE_WINDOW_NS),
        amplitudes=AMPLITUDES,
        peak_excursions=np.array([row["peak_excursion"] for row in wiener_rows]),
        truth_waveforms=np.array([row["truth"] for row in wiener_rows]),
        wiener_reconstructions=np.array([row["recon"] for row in wiener_rows]),
        wiener_rmse=np.array([row["rmse"] for row in wiener_rows]),
        wiener_peak_ratio=np.array([row["peak_ratio"] for row in wiener_rows]),
        linear_model_relative_error=np.asarray(model_errors),
        workpoint_flux=flux_scan,
        workpoint_signed_kappa=kappa_scan,
        workpoint_rmse=workpoint_rmse,
        lm_amplitudes=LM_AMPLITUDES,
        lm_rmse=lm_rmse,
        lm_peak_ratio=lm_peak_ratio,
        lm_runtime_s=lm_runtime,
        lm_status=lm_status,
        lm_error=lm_error,
    )

    metadata = {
        "figure": "D3",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(),
        "work_point_phi0": C.OPTIMAL_FLUX,
        "control": {
            "envelope": "gaussian",
            "rotation_angle_rad": ROTATION_ANGLE,
            "pulse_window_ns": PULSE_WINDOW_NS,
            "envelope_sigma_ns": GAUSSIAN_SIGMA_NS,
            "cache_version": CACHE_VERSION,
        },
        "wiener_lambda": C.WIENER_LAMBDA,
        "lm_lambda": C.LM_LAMBDA,
        "lm_basis": "fourier",
        "lm_n_basis": C.LM_N_BASIS,
        "lm_max_iter": C.LM_MAX_ITER,
        "amplitudes": AMPLITUDES.tolist(),
        "lm_amplitudes": LM_AMPLITUDES.tolist(),
        "amplitude_axis_note": "Figure uses realized peak flux_samples, not the type=4 parameter",
        "model_error": "Relative L2 error of exact vs local-linear frequency shift",
        "source_caches": [
            "result_sqc/reconstruction/D3_nonlinearity_boundary/gaussian_amplitude_*_sqc.npz",
            "result_sqc/reconstruction/D3_nonlinearity_boundary/gaussian_workpoint_scan_sqc.npz",
            "result_sqc/reconstruction/D3_nonlinearity_boundary/lm_amplitude_*_sqc.npz",
        ],
    }
    with open(os.path.join(out, "nonlinearity_boundary_metadata.json"),
              "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)

    lines = ["D3 Gaussian-control nonlinearity boundary metrics", "=" * 78,
             ("control: Gaussian pi/2, pulse window=10 ns, "
              "envelope sigma=2 ns"),
             "",
             "amp_param  peak_flux  model_err  wiener_rmse  wiener_ratio"]
    for row, model_error in zip(wiener_rows, model_errors):
        lines.append(
            f"{row['amplitude']:9.3f}  {row['peak_excursion']:9.5f}  "
            f"{model_error:9.4f}  {row['rmse']:11.4e}  {row['peak_ratio']:12.4f}"
        )
    lines.extend(["", "LM boundary points",
                  "amp_param  rmse         ratio    runtime_s  status"])
    for amplitude, rmse, ratio, runtime, status, error in zip(
        LM_AMPLITUDES, lm_rmse, lm_peak_ratio, lm_runtime, lm_status, lm_error
    ):
        lines.append(
            f"{amplitude:9.3f}  {rmse:11.4e}  {ratio:8.4f}  "
            f"{runtime:9.1f}  {status} {error}"
        )
    with open(os.path.join(out, "nonlinearity_boundary_metrics.txt"),
              "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def plot_figure(wiener_rows, lm_rows, workpoint_scan, model_errors):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 8,
        "axes.titlesize": 9, "axes.labelsize": 8,
        "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    red, blue, green, gray = "#D62728", "#1F77B4", "#2E8B57", "#666666"
    bias = C.OPTIMAL_FLUX
    peak_excursions = np.array([row["peak_excursion"] for row in wiener_rows])
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), constrained_layout=True)

    ax = axes[0, 0]
    phi = np.linspace(bias - 0.07, bias + 0.07, 700)
    exact = omega01(phi) / (2.0 * np.pi)
    tangent = (omega01(bias) + signed_kappa(bias) * (phi - bias)) / (2.0 * np.pi)
    ax.plot(phi, exact, color="black", lw=1.7, label="Exact dispersion")
    ax.plot(phi, tangent, color=blue, ls="--", lw=1.4, label="Local linear model")
    ax.axvline(bias, color=gray, ls=":", lw=1.0)
    ax.axvspan(bias, bias + peak_excursions[-1], color="#F2C14E", alpha=0.24,
               label="largest tested excursion")
    ax.scatter([bias], [omega01(bias) / (2.0 * np.pi)], color=red, s=24, zorder=4)
    ax.set(xlabel=r"Flux bias $\Phi/\Phi_0$", ylabel=r"$\omega_{01}/2\pi$ (GHz)",
           title="Transmon dispersion near the work point")
    ax.legend(frameon=False, loc="best")

    ax = axes[0, 1]
    ax.plot(peak_excursions, 100.0 * np.asarray(model_errors), "o-",
            color=green, lw=1.5, ms=4)
    ax.set(xlabel=r"Peak flux excursion ($\Phi_0$)",
           ylabel="Relative model error (%)",
           title="Breakdown of the local linear model")
    ax.grid(True, ls=":", alpha=0.35)

    ax = axes[1, 0]
    _flux, kappa_scan, workpoint_rmse = workpoint_scan
    order = np.argsort(kappa_scan)
    ax.semilogy(kappa_scan[order], workpoint_rmse[order], color=blue, lw=1.5)
    kappa_bias = signed_kappa(bias)
    nearest = int(np.argmin(np.abs(kappa_scan - kappa_bias)))
    ax.scatter([kappa_scan[nearest]], [workpoint_rmse[nearest]], color=red,
               s=27, zorder=4, label="selected work point")
    ax.axvline(0.0, color=gray, ls=":", lw=1.0)
    ax.set(xlabel=r"Signed slope $\kappa_\phi=d\omega_{01}/d\Phi$",
           ylabel="Wiener RMSE", title="Reconstruction error versus work point")
    ax.legend(frameon=False)
    ax.grid(True, which="both", ls=":", alpha=0.35)

    ax = axes[1, 1]
    wiener_rmse = np.array([row["rmse"] for row in wiener_rows])
    ax.loglog(peak_excursions, wiener_rmse, "o-", color=red, lw=1.5,
              ms=4, label="Wiener")
    lm_rmse = np.array([float(np.asarray(row["lm_rmse"])) for row in lm_rows])
    lm_peak = [next(row for row in wiener_rows
                    if np.isclose(row["amplitude"], amplitude))["peak_excursion"]
               for amplitude in LM_AMPLITUDES]
    finite = np.isfinite(lm_rmse)
    if np.any(finite):
        ax.plot(np.asarray(lm_peak)[finite], lm_rmse[finite], "s--", color=blue,
                lw=1.3, ms=4.5, label="LM boundary points")
    ax.set(xlabel=r"Peak flux excursion ($\Phi_0$)", ylabel="Reconstruction RMSE",
           title="High-amplitude reconstruction boundary")
    ax.legend(frameon=False)
    ax.grid(True, which="both", ls=":", alpha=0.35)

    for label, panel in zip("abcd", axes.flat):
        panel.text(0.02, 0.98, f"({label})", transform=panel.transAxes,
                   ha="left", va="top", fontweight="bold")

    out = C.out_dir(SUBDIR)
    fig.savefig(os.path.join(out, "nonlinearity_boundary_sqc.png"), dpi=300)
    fig.savefig(os.path.join(out, "nonlinearity_boundary_sqc.pdf"))
    plt.close(fig)


def main():
    no_lm = "--no-lm" in sys.argv
    recompute_lm = "--recompute-lm" in sys.argv
    force = "--force" in sys.argv
    recompute_wiener = force or "--recompute-wiener" in sys.argv
    recompute_workpoint = force or "--recompute-workpoint" in sys.argv
    recompute_lm = force or recompute_lm
    wiener_rows = load_wiener_rows(recompute=recompute_wiener)
    workpoint_scan = load_workpoint_scan(recompute=recompute_workpoint)
    model_errors = np.array([
        linear_model_error(row["truth"], C.OPTIMAL_FLUX) for row in wiener_rows
    ])
    if no_lm:
        lm_rows = [{
            "amplitude": np.array(amplitude), "lm_rmse": np.array(np.nan),
            "lm_peak_ratio": np.array(np.nan), "runtime_s": np.array(0.0),
            "status": np.array("not_run"), "error": np.array("--no-lm"),
        } for amplitude in LM_AMPLITUDES]
    else:
        lm_rows = []
        for amplitude in LM_AMPLITUDES:
            print(f"LM amplitude={amplitude:g} ...", flush=True)
            row = run_lm(float(amplitude), recompute=recompute_lm)
            print(f"  status={np.asarray(row['status']).item()} "
                  f"rmse={float(np.asarray(row['lm_rmse'])):.4e} "
                  f"runtime={float(np.asarray(row['runtime_s'])):.1f}s", flush=True)
            lm_rows.append(row)
    save_outputs(wiener_rows, lm_rows, workpoint_scan, model_errors)
    plot_figure(wiener_rows, lm_rows, workpoint_scan, model_errors)
    print(f"-> {C.out_dir(SUBDIR)}")


if __name__ == "__main__":
    main()
