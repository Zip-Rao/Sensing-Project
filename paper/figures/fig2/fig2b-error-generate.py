#!/usr/bin/env python3
"""Generate final Fig. 2(b): signed local-frequency estimation error."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sqc.control.sequence import create_ramsey_pulse
from sqc.devices.transmon import TransmonQubit
from sqc.reconstruction.kernel import KernelEstimator
from fig2_simulation_common import (
    CANDIDATE_INTERVAL_MHZ,
    ESTIMATOR_DATA,
    REFERENCE_FLUX_PHI0,
    RESPONSE_DATA,
    TWO_PI,
    local_flux_from_detuning_mhz,
    load_npz,
    save_npz,
)


HERE = Path(__file__).resolve().parent
SOURCE = RESPONSE_DATA
STEM = "fig2b-frequency-error"
STYLE = ROOT.parent / "scholaraio" / ".claude" / "skills" / "draw" / "nature_pub.mplstyle"


def kernel_coefficients(source: dict[str, np.ndarray]) -> tuple[float, float]:
    """Return full-kernel G1/G3 in the physical Delta convention."""
    t_rabi = source["t_rabi_ns"]
    envelope = source["envelope_samples"]
    qubit = TransmonQubit(
        EC=TWO_PI * float(source["EC_over_h_GHz"]),
        EJ=TWO_PI * float(source["EJ_over_h_GHz"]),
        T1=float(source["T1_ns"]), T2=float(source["T2_ns"]),
        flux=0.0, state=0, n_levels=2,
    )
    omega_d = float(source["omega_d_rad_GHz"])
    common = dict(
        t_rabi=t_rabi, tau=0.0, omega_d=omega_d,
        phase1=np.pi / 2.0, qubit=qubit,
        rotation_angle=float(source["rotation_angle_rad"]), envelope=envelope,
    )
    pulse_plus = create_ramsey_pulse(phase2=0.0, **common)
    pulse_minus = create_ramsey_pulse(phase2=np.pi, **common)
    estimator = KernelEstimator(
        mode="omega", method="sim", order=3, extract_off_diagonal=True,
    )
    result_plus = estimator.estimate_full(pulse_plus, qubit)
    result_minus = estimator.estimate_full(pulse_minus, qubit)
    tk = np.asarray(result_plus.t_samples, dtype=float)
    k1 = (np.asarray(result_plus.kernels[0]) - np.asarray(result_minus.kernels[0])) / 2.0
    k3 = (np.asarray(result_plus.kernels[2]) - np.asarray(result_minus.kernels[2])) / 2.0
    g1_raw = float(np.trapezoid(k1, tk))
    g3_raw = float(
        np.trapezoid(
            np.trapezoid(np.trapezoid(k3, tk, axis=0), tk, axis=0),
            tk, axis=0,
        )
    )
    # The Virtual-Z kernel uses delta_omega=-Delta. Odd coefficients change sign.
    return -g1_raw, -g3_raw


def invert_cubic(
    response: np.ndarray, g1: float, g3: float, branch_limit_mhz: float,
) -> np.ndarray:
    """Invert the root continuously connected to zero; reject past the fold."""
    p = np.asarray(response, dtype=float)
    estimate = p / g1
    for _ in range(50):
        residual = g1 * estimate + (g3 / 6.0) * estimate**3 - p
        derivative = g1 + (g3 / 2.0) * estimate**2
        good = np.abs(derivative) > 1e-12
        update = np.zeros_like(estimate)
        update[good] = residual[good] / derivative[good]
        estimate -= update
        if np.max(np.abs(update[good])) < 1e-12:
            break
    estimate_mhz = estimate / TWO_PI * 1e3
    estimate_mhz[np.abs(estimate_mhz) > branch_limit_mhz] = np.nan
    return estimate_mhz


def estimator_errors(
    delta_mhz: np.ndarray,
    response: np.ndarray,
    g1_fit: float,
    g3_fit: float,
    g1_kernel: float,
    g3_kernel: float,
    fold_fit_mhz: float,
    fold_kernel_mhz: float,
) -> np.ndarray:
    linear = response / g1_fit / TWO_PI * 1e3
    cubic_fit = invert_cubic(response, g1_fit, g3_fit, fold_fit_mhz)
    cubic_kernel = invert_cubic(response, g1_kernel, g3_kernel, fold_kernel_mhz)
    cubic_fit[np.abs(delta_mhz) > fold_fit_mhz] = np.nan
    cubic_kernel[np.abs(delta_mhz) > fold_kernel_mhz] = np.nan
    estimates = np.vstack((linear, cubic_fit, cubic_kernel))
    return estimates - np.asarray(delta_mhz)[None, :]


def compute(source: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    g1_fit = float(source["G1_ns"])
    g3_fit = float(source["G3_ns3"])
    g1_kernel, g3_kernel = kernel_coefficients(source)
    fold_fit_mhz = float(source["delta_fold_mhz"])
    fold_kernel_mhz = (
        np.sqrt(-2.0 * g1_kernel / g3_kernel) / TWO_PI * 1e3
        if g1_kernel * g3_kernel < 0.0 else np.nan
    )
    dense_error = estimator_errors(
        source["dense_delta_mhz"], source["dense_pdiff_corrected"],
        g1_fit, g3_fit, g1_kernel, g3_kernel, fold_fit_mhz, fold_kernel_mhz,
    )
    validation_error = estimator_errors(
        source["validation_delta_mhz"], source["validation_pdiff_corrected"],
        g1_fit, g3_fit, g1_kernel, g3_kernel, fold_fit_mhz, fold_kernel_mhz,
    )
    calibration_error = estimator_errors(
        source["calibration_delta_mhz"], source["calibration_pdiff_corrected"],
        g1_fit, g3_fit, g1_kernel, g3_kernel, fold_fit_mhz, fold_kernel_mhz,
    )
    ec = TWO_PI * float(source["EC_over_h_GHz"])
    ej = TWO_PI * float(source["EJ_over_h_GHz"])
    omega_ref = float(source["omega_d_rad_GHz"])
    dense_flux = local_flux_from_detuning_mhz(
        source["dense_delta_mhz"], ec, ej, omega_ref,
    )
    validation_flux = local_flux_from_detuning_mhz(
        source["validation_delta_mhz"], ec, ej, omega_ref,
    )
    calibration_flux = local_flux_from_detuning_mhz(
        source["calibration_delta_mhz"], ec, ej, omega_ref,
    )
    max_frequency_error = float(source.get("validation_max_frequency_error_mhz", 0.5))
    delta_val = CANDIDATE_INTERVAL_MHZ
    mask = np.abs(source["validation_delta_mhz"]) <= delta_val + 1e-12
    cubic_errors = validation_error[1:, mask]
    if not (
        np.all(np.isfinite(cubic_errors))
        and np.max(np.abs(cubic_errors)) <= max_frequency_error
        and delta_val < min(fold_fit_mhz, fold_kernel_mhz)
    ):
        raise ValueError("frozen candidate interval fails its estimator qualification")
    return {
        "dense_delta_mhz": source["dense_delta_mhz"],
        "dense_error_mhz": dense_error,
        "validation_delta_mhz": source["validation_delta_mhz"],
        "validation_error_mhz": validation_error,
        "calibration_delta_mhz": source["calibration_delta_mhz"],
        "calibration_error_mhz": calibration_error,
        "dense_flux_phi0": dense_flux,
        "validation_flux_phi0": validation_flux,
        "calibration_flux_phi0": calibration_flux,
        "analytic_truth_offset_mhz": source["dense_delta_mhz"],
        "dense_estimate_offset_mhz": dense_error + source["dense_delta_mhz"][None, :],
        "validation_estimate_offset_mhz": (
            validation_error + source["validation_delta_mhz"][None, :]
        ),
        "reference_flux_phi0": np.array(REFERENCE_FLUX_PHI0),
        "reference_frequency_ghz": np.array(omega_ref / TWO_PI),
        "flux_truth_role": np.array("withheld_analytic_dispersion_scoring_only"),
        "labels": np.array(("First-order", "Cubic fit", "Cubic kernel")),
        "G1_fit_ns": np.array(g1_fit),
        "G3_fit_ns3": np.array(g3_fit),
        "G1_kernel_ns": np.array(g1_kernel),
        "G3_kernel_ns3": np.array(g3_kernel),
        "delta_fold_mhz": source["delta_fold_mhz"],
        "delta_fold_kernel_mhz": np.array(fold_kernel_mhz),
        "delta_val_mhz": np.array(delta_val),
        "validation_max_frequency_error_mhz": np.array(max_frequency_error),
        "candidate_interval_role": np.array("simulation_defined_candidate_interval"),
        "sigma": source["sigma"],
        "uncertainty_source": source["uncertainty_source"],
        "truth_role": np.array("withheld_analytic_simulation_ground_truth"),
        "validation_role": source["validation_role"],
        "source_npz": np.array(SOURCE.name),
    }


def plot(data: dict[str, np.ndarray], stem: str = STEM, panel_label: str = "b") -> None:
    if STYLE.exists():
        plt.style.use(STYLE)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 7.0})
    colors = ("#8A9AA8", "#C98D8A", "#16705A")
    markers = ("o", "s", "^")
    fig = plt.figure(figsize=(3.50, 3.25))
    grid = fig.add_gridspec(2, 1, height_ratios=(1.15, 1.0), hspace=0.08)
    ax_top = fig.add_subplot(grid[0])
    ax = fig.add_subplot(grid[1], sharex=ax_top)
    dense_x = data["dense_flux_phi0"]
    val_x = data["validation_flux_phi0"]
    valid = float(data["delta_val_mhz"])
    display_limit = min(18.0, float(data["delta_fold_kernel_mhz"]))
    dense_mask = np.isfinite(dense_x) & (np.abs(data["dense_delta_mhz"]) <= display_limit)
    val_mask = np.isfinite(val_x) & (np.abs(data["validation_delta_mhz"]) <= display_limit)
    candidate_flux = local_flux_from_detuning_mhz(
        np.array([-valid, valid]),
        TWO_PI * 0.2, TWO_PI * 10.0,
        TWO_PI * float(data["reference_frequency_ghz"]),
    )
    display_flux = local_flux_from_detuning_mhz(
        np.array([-display_limit, display_limit]),
        TWO_PI * 0.2, TWO_PI * 10.0,
        TWO_PI * float(data["reference_frequency_ghz"]),
    )
    for target in (ax_top, ax):
        target.axvspan(candidate_flux[0], candidate_flux[1], color="#D9EAD3",
                       alpha=0.75, lw=0, zorder=0)
    ax_top.plot(dense_x[dense_mask], data["analytic_truth_offset_mhz"][dense_mask], color="#222222",
                lw=1.2, label="Withheld analytic truth", zorder=2)
    kernel = 2
    ax_top.plot(dense_x[dense_mask], data["dense_estimate_offset_mhz"][kernel, dense_mask],
                color=colors[kernel], lw=1.15, label="Cubic-kernel estimate", zorder=3)
    ax_top.scatter(val_x[val_mask], data["validation_estimate_offset_mhz"][kernel, val_mask], s=20,
                   marker=markers[kernel], facecolor="white", edgecolor=colors[kernel],
                   linewidth=0.8, zorder=4)
    ax_top.set_ylabel(r"$\widehat f_{01}-f_{\rm ref}$ (MHz)")
    ax_top.tick_params(labelbottom=False)
    ax_top.legend(loc="upper left", frameon=False, fontsize=5.6)
    ax_top.text(0.98, 0.04, "analytic truth used for scoring only",
                transform=ax_top.transAxes, ha="right", va="bottom",
                fontsize=4.9, color="#555555")
    ax.axhline(0.0, color="#777777", lw=0.65, zorder=1)
    for idx, (label, color, marker) in enumerate(zip(data["labels"], colors, markers)):
        alpha = 1.0 if idx == kernel else 0.65
        width = 1.2 if idx == kernel else 0.8
        ax.plot(dense_x[dense_mask], data["dense_error_mhz"][idx, dense_mask], color=color, lw=width,
                alpha=alpha, label=str(label), zorder=2 + idx)
        ax.scatter(val_x[val_mask], data["validation_error_mhz"][idx, val_mask], s=18, marker=marker,
                   facecolor="white", edgecolor=color, linewidth=0.8,
                   alpha=alpha, zorder=6 + idx)
    ax.plot(
        data["calibration_flux_phi0"][np.isfinite(data["calibration_flux_phi0"])],
        np.full(np.count_nonzero(np.isfinite(data["calibration_flux_phi0"])), 0.018),
        linestyle="none", marker="|", markersize=4.5, color="#777777",
        transform=ax.get_xaxis_transform(), label="Calibration grid", zorder=8,
    )
    ax.set_xlabel(r"Local flux bias $\Phi/\Phi_0$")
    ax.set_ylabel("Signed residual (MHz)")
    ax.set_xlim(float(display_flux[0]), float(display_flux[1]))
    finite_local = data["dense_error_mhz"][:, dense_mask]
    bound = max(0.5, 1.15 * float(np.nanmax(np.abs(finite_local))))
    ax.set_ylim(-bound, bound)
    ax_top.text(-0.16, 1.03, panel_label, transform=ax_top.transAxes, ha="left", va="bottom",
            fontsize=8.5, fontweight="bold")
    ax.text(0.98, 0.94, r"deterministic simulation, $\sigma=0$",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5,
            color="#555555")
    ax.legend(loc="upper left", ncol=2, frameon=False, handlelength=1.5,
              columnspacing=0.8, fontsize=4.9)
    fig.subplots_adjust(left=0.21, right=0.97, bottom=0.15, top=0.96)
    fig.savefig(HERE / f"{stem}.pdf")
    fig.savefig(HERE / f"{stem}.svg")
    fig.savefig(HERE / f"{stem}.png", dpi=600)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    source = load_npz(SOURCE)
    data_path = ESTIMATOR_DATA
    if args.plot_only:
        data = load_npz(data_path)
        required = {"dense_flux_phi0", "dense_estimate_offset_mhz"}
        if not required.issubset(data):
            data = compute(source)
        else:
            data["candidate_interval_role"] = np.array("simulation_defined_candidate_interval")
            data["source_npz"] = np.array(SOURCE.name)
        save_npz(data_path, data)
    else:
        data = compute(source)
        save_npz(data_path, data)
    plot(data)
    plot(data, "figS3-frequency-recovery-simulation", "S3")
    valid = np.abs(data["validation_delta_mhz"]) <= float(data["delta_val_mhz"])
    print(f"G1 kernel={float(data['G1_kernel_ns']):.8g} ns")
    print(f"G3 kernel={float(data['G3_kernel_ns3']):.8g} ns^3")
    print(f"kernel fold={float(data['delta_fold_kernel_mhz']):.4f} MHz")
    for idx, label in enumerate(data["labels"]):
        errors = data["validation_error_mhz"][idx, valid]
        finite = errors[np.isfinite(errors)]
        print(f"{label}: bias={np.mean(finite):.6f}, MAE={np.mean(np.abs(finite)):.6f}, "
              f"RMSE={np.sqrt(np.mean(finite**2)):.6f}, max={np.max(np.abs(finite)):.6f} MHz")


if __name__ == "__main__":
    main()
