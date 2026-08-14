#!/usr/bin/env python3
"""Generate a D2 preview for Gaussian-pulse parameter dependence.

The preview keeps the existing D2 definitions of double-peak resolvability,
projection noise, and Wiener reconstruction, but recomputes the first-order
response kernel for every Gaussian pulse width and rotation angle.  It uses
the linear convolution model associated with each recomputed kernel; selected
points can be promoted to full nonlinear forward runs after the layout is
approved.

Outputs are deliberately named ``preview`` and do not replace the formal D2
figure or the report image.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, HERE)

import _d2_common as D  # noqa: E402


STEM = "d2_pulse_parameter_preview_sqc"
PULSE_WIDTHS_NS = np.array([4.0, 6.0, 8.0, 10.0])
ROTATION_ANGLES = np.array([np.pi / 6.0, np.pi / 4.0, np.pi / 2.0])
ANGLE_LABELS = [r"$\pi/6$", r"$\pi/4$", r"$\pi/2$"]
SEPARATIONS_NS = np.array(
    [8.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
     17.0, 18.0, 20.0, 22.0, 25.0, 30.0, 35.0, 40.0]
)
AMPLITUDES = np.array(
    [1e-5, 2e-5, 3e-5, 5e-5, 7e-5, 1e-4,
     2e-4, 3e-4, 5e-4, 1e-3, 2e-3, 3e-3]
)
N_SHOTS = np.array([1_000, 10_000, 100_000])
N_REPEAT = 96
N_REPEAT_ZERO = 192
SENSITIVITY_PULSE_WIDTH_NS = 10.0
AWG_DT_NS = 0.5


def angle_key(angle: float) -> str:
    return f"a{angle / np.pi:.8f}pi"


def kernel_key(width: float, angle: float) -> str:
    return f"T{width:g}_{angle_key(angle)}"


def compute_kernel(width_ns: float, angle: float) -> dict[str, np.ndarray | float]:
    """Build a Gaussian Ramsey pair and estimate its flux-response kernel."""
    from sqc.control.sequence import create_ramsey_pulse
    from sqc.reconstruction.kernel import KernelEstimator

    q = D.C.make_sqc_qubit(D.C.OPTIMAL_FLUX)
    t_rabi = np.arange(0.0, float(width_ns), AWG_DT_NS)
    sigma = float(width_ns) / 5.0
    pulse = create_ramsey_pulse(
        t_rabi,
        tau=0.0,
        omega_d=q.frequency,
        phase1=np.pi / 2.0,
        phase2=0.0,
        qubit=q,
        rotation_angle=float(angle),
        envelope="gaussian",
        envelope_sigma=sigma,
    )
    result = KernelEstimator(
        mode="flux", method="exp", order=1
    ).estimate_full(pulse, q)
    kernel = np.asarray(result.k1, dtype=float)
    t_kernel = np.asarray(result.t_samples, dtype=float)
    width = D.kernel_effective_width(t_kernel, kernel)
    return {
        "kernel": kernel,
        "kernel_time": t_kernel,
        "kernel_sigma_ns": float(width["sigma"]),
        "pulse_time": t_rabi,
        "pulse_sigma_ns": sigma,
        "control_angle_rad": float(pulse.pulses[0].get_angle_simple()),
    }


def build_kernels() -> dict[str, dict[str, np.ndarray | float]]:
    kernels = {}
    for width in PULSE_WIDTHS_NS:
        for angle in ROTATION_ANGLES:
            key = kernel_key(width, angle)
            print(
                f"kernel: T_p={width:g} ns, alpha={angle / np.pi:.3f} pi"
            )
            kernels[key] = compute_kernel(width, angle)
    return kernels


def linear_measurement(signal: np.ndarray, kernel: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Riemann-sum first-order response and its scan axis."""
    dt = float(D.C.T_LIST[1] - D.C.T_LIST[0])
    delta_p = np.convolve(np.asarray(signal, float), kernel, mode="full") * dt
    scan = np.arange(len(delta_p), dtype=float) * dt
    return scan, delta_p


def resolution_scan(kernels):
    shape = (len(ROTATION_ANGLES), len(PULSE_WIDTHS_NS), len(SEPARATIONS_NS))
    cv = np.full(shape, np.nan)
    nr = np.full(shape, np.nan)
    pos = np.full(shape, np.nan)
    resolved = np.zeros(shape, dtype=bool)
    dt_min = np.full((len(ROTATION_ANGLES), len(PULSE_WIDTHS_NS)), np.nan)
    kernel_sigma = np.full_like(dt_min, np.nan)

    for ai, angle in enumerate(ROTATION_ANGLES):
        for wi, width in enumerate(PULSE_WIDTHS_NS):
            item = kernels[kernel_key(width, angle)]
            kernel = np.asarray(item["kernel"], float)
            kernel_sigma[ai, wi] = float(item["kernel_sigma_ns"])
            for si, sep in enumerate(SEPARATIONS_NS):
                truth = D.double_peak(float(sep))
                scan, delta_p = linear_measurement(truth, kernel)
                tr, rec = D.wiener(delta_p, kernel, scan)
                metric = D.valley_contrast(tr, rec, float(sep))
                cv[ai, wi, si] = metric["cv"]
                nr[ai, wi, si] = D.nrmse(rec, truth)
                pos[ai, wi, si] = metric["peak_position_error"]
                resolved[ai, wi, si] = D.is_resolved(metric, nr[ai, wi, si])

            ok = np.where(resolved[ai, wi])[0]
            if len(ok):
                dt_min[ai, wi] = float(SEPARATIONS_NS[ok[0]])
            print(
                f"resolution: T_p={width:g} ns, alpha={angle / np.pi:.3f} pi, "
                f"sigma_k={kernel_sigma[ai, wi]:.3f} ns, "
                f"dt_min={dt_min[ai, wi]:g} ns"
            )

    return {
        "valley_contrast": cv,
        "nrmse": nr,
        "peak_position_error_ns": pos,
        "resolved": resolved,
        "dt_min_ns": dt_min,
        "kernel_sigma_ns": kernel_sigma,
    }


def zero_population(angle: float) -> float:
    """Ideal baseline for two orthogonal equal-angle Ramsey rotations."""
    return 0.5 * float(np.sin(angle) ** 2)


def sensitivity_scan(kernels):
    n_a = len(ROTATION_ANGLES)
    n_s = len(N_SHOTS)
    n_amp = len(AMPLITUDES)
    snr = np.full((n_a, n_s, n_amp), np.nan)
    snr_sem = np.full_like(snr, np.nan)
    sigma0 = np.full((n_a, n_s), np.nan)
    a_min = np.full((n_a, n_s), np.nan)
    a_min_status = np.full((n_a, n_s), "never", dtype="U24")
    mean_est = np.full_like(snr, np.nan)
    template = D.template()

    for ai, angle in enumerate(ROTATION_ANGLES):
        kernel = np.asarray(
            kernels[kernel_key(SENSITIVITY_PULSE_WIDTH_NS, angle)]["kernel"],
            float,
        )
        p0_value = zero_population(float(angle))
        _, delta_unit = linear_measurement(template, kernel)
        scan = np.arange(len(delta_unit), dtype=float) * (
            D.C.T_LIST[1] - D.C.T_LIST[0]
        )
        p_ref = np.full(len(delta_unit), p0_value, dtype=float)

        for ni, shots in enumerate(N_SHOTS):
            zero_est = np.empty(N_REPEAT_ZERO, dtype=float)
            for ri in range(N_REPEAT_ZERO):
                seed = D.SEED_BASE + 1_000_000 * ai + 10_000 * ni + ri
                _, _, dp = D.sample_shots(p_ref, p_ref, int(shots), seed)
                tr, rec = D.wiener(dp, kernel, scan)
                zero_est[ri] = D.project_amplitude(tr, rec)
            sigma0[ai, ni] = float(np.std(zero_est, ddof=1))

            for mi, amplitude in enumerate(AMPLITUDES):
                p_sig = np.clip(p_ref + amplitude * delta_unit, 0.0, 1.0)
                estimates = np.empty(N_REPEAT, dtype=float)
                for ri in range(N_REPEAT):
                    seed = (
                        D.SEED_BASE + 1_000_000 * ai + 10_000 * ni
                        + 100 * mi + ri
                    )
                    _, _, dp = D.sample_shots(p_sig, p_ref, int(shots), seed)
                    tr, rec = D.wiener(dp, kernel, scan)
                    estimates[ri] = D.project_amplitude(tr, rec)
                mean_est[ai, ni, mi] = float(np.mean(estimates))
                sem = float(np.std(estimates, ddof=1) / np.sqrt(N_REPEAT))
                snr[ai, ni, mi] = mean_est[ai, ni, mi] / sigma0[ai, ni]
                snr_sem[ai, ni, mi] = sem / sigma0[ai, ni]

            value, status = D.first_crossing(
                AMPLITUDES, snr[ai, ni], D.SNR_THRESHOLD, log=True
            )
            a_min[ai, ni] = value
            a_min_status[ai, ni] = status
            print(
                f"sensitivity: alpha={angle / np.pi:.3f} pi, N={shots}, "
                f"sigma0={sigma0[ai, ni]:.3e}, A_min={value:.3e} ({status})"
            )

    return {
        "template": template,
        "snr": snr,
        "snr_sem": snr_sem,
        "sigma0": sigma0,
        "a_min": a_min,
        "a_min_status": a_min_status,
        "mean_estimate": mean_est,
    }


def save_data(kernels, resolution, sensitivity):
    arrays = {
        "pulse_widths_ns": PULSE_WIDTHS_NS,
        "rotation_angles_rad": ROTATION_ANGLES,
        "separations_ns": SEPARATIONS_NS,
        "amplitudes": AMPLITUDES,
        "n_shots": N_SHOTS,
        "resolution_valley_contrast": resolution["valley_contrast"],
        "resolution_nrmse": resolution["nrmse"],
        "resolution_peak_position_error_ns": resolution[
            "peak_position_error_ns"
        ],
        "resolution_resolved": resolution["resolved"],
        "dt_min_ns": resolution["dt_min_ns"],
        "kernel_sigma_ns": resolution["kernel_sigma_ns"],
        "sensitivity_template": sensitivity["template"],
        "snr": sensitivity["snr"],
        "snr_sem": sensitivity["snr_sem"],
        "sigma_Ahat_zero": sensitivity["sigma0"],
        "a_min": sensitivity["a_min"],
        "a_min_status": sensitivity["a_min_status"],
        "mean_amplitude_estimate": sensitivity["mean_estimate"],
    }
    for width in PULSE_WIDTHS_NS:
        for angle in ROTATION_ANGLES:
            item = kernels[kernel_key(width, angle)]
            key = kernel_key(width, angle)
            arrays[f"kernel_{key}"] = item["kernel"]
            arrays[f"kernel_time_{key}"] = item["kernel_time"]
            arrays[f"pulse_time_{key}"] = item["pulse_time"]
            arrays[f"control_angle_{key}"] = np.array(item["control_angle_rad"])
    path = os.path.join(HERE, STEM + ".npz")
    np.savez(path, **arrays)

    metadata = {
        "figure": "D2 pulse-parameter preview",
        "status": "layout and parameter-scan preview; not formal report data",
        "generated": datetime.now().isoformat(timespec="seconds"),
        "control": {
            "envelope": "gaussian",
            "sigma_rule": "sigma = T_p / 5",
            "pulse_widths_ns": PULSE_WIDTHS_NS.tolist(),
            "rotation_angles_rad": ROTATION_ANGLES.tolist(),
            "sensitivity_pulse_width_ns": SENSITIVITY_PULSE_WIDTH_NS,
        },
        "model": {
            "kernel": "recomputed by KernelEstimator for every (T_p, alpha)",
            "measurement": "first-order Riemann-sum convolution with each recomputed kernel",
            "reconstruction": f"Wiener, lambda_W={D.C.WIENER_LAMBDA:g}",
            "projection_noise": "independent binomial sampling of signal and reference arms",
            "zero_population": "ideal two orthogonal equal-angle rotations: p0=sin(alpha)^2/2",
        },
        "resolution_criteria": {
            "C_v_min": D.C_MIN,
            "NRMSE_max": D.NRMSE_MAX,
            "peak_position_tolerance_ns": D.PEAK_POSITION_TOL_NS,
            "dt_min_grid_note": "first passing scanned separation; no interpolation in preview",
        },
        "sensitivity": {
            "SNR_threshold": D.SNR_THRESHOLD,
            "n_repeat": N_REPEAT,
            "n_repeat_zero": N_REPEAT_ZERO,
            "seed_base": D.SEED_BASE,
        },
        "notes": [
            "This preview does not overwrite the existing formal D2 artifacts.",
            "A formal run should validate selected parameter points with the full nonlinear forward model.",
            "The 0.5 ns AWG grid is a numerical control grid, not the reported time resolution.",
        ],
    }
    meta_path = os.path.join(HERE, STEM + "_metadata.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    lines = [
        "D2 Gaussian-pulse parameter preview",
        "=" * 72,
        "Control: Gaussian envelope, sigma=T_p/5, Wiener lambda_W=5",
        "Resolution: first scanned spacing satisfying Cv, NRMSE and peak-position criteria",
        "",
    ]
    for ai, label in enumerate(ANGLE_LABELS):
        for wi, width in enumerate(PULSE_WIDTHS_NS):
            lines.append(
                f"alpha={label:>7s}, T_p={width:>4g} ns: "
                f"sigma_k={resolution['kernel_sigma_ns'][ai, wi]:.4f} ns, "
                f"Delta_t_min={resolution['dt_min_ns'][ai, wi]:.4f} ns"
            )
    lines.extend(["", "Sensitivity at T_p=10 ns (A_min for SNR>=3)"])
    for ai, label in enumerate(ANGLE_LABELS):
        for ni, shots in enumerate(N_SHOTS):
            lines.append(
                f"alpha={label:>7s}, N_shot={shots:>6d}: "
                f"A_min={sensitivity['a_min'][ai, ni]:.6e}, "
                f"sigma0={sensitivity['sigma0'][ai, ni]:.6e}"
            )
    metrics_path = os.path.join(HERE, STEM + "_metrics.txt")
    with open(metrics_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return path, meta_path, metrics_path


def plot(kernels, resolution, sensitivity):
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.linewidth": 0.8,
        "mathtext.fontset": "dejavusans",
    })
    colors = ["#6a3d9a", "#1f78b4", "#238b45"]
    markers = ["o", "s", "^"]
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.2))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    # (a) Representative kernels: pulse-width dependence plus one angle check.
    representative = [
        (4.0, np.pi / 2.0, "#d73027", "-", r"$T_p=4$ ns, $\alpha=\pi/2$"),
        (6.0, np.pi / 2.0, "#fc8d59", "-", r"$T_p=6$ ns, $\alpha=\pi/2$"),
        (10.0, np.pi / 2.0, "#252525", "-", r"$T_p=10$ ns, $\alpha=\pi/2$"),
        (10.0, np.pi / 4.0, "#1f78b4", "--", r"$T_p=10$ ns, $\alpha=\pi/4$"),
    ]
    centered_limits = []
    for width, angle, color, style, label in representative:
        item = kernels[kernel_key(width, angle)]
        k = np.abs(np.asarray(item["kernel"], float))
        t_kernel = np.asarray(item["kernel_time"], float)
        centroid = D.kernel_effective_width(t_kernel, k)["t_centroid"]
        t_centered = t_kernel - centroid
        centered_limits.extend([float(t_centered.min()), float(t_centered.max())])
        ax_a.plot(
            t_centered, k / k.max(), color=color, ls=style,
            lw=2.0 if (width == 10.0 and np.isclose(angle, np.pi / 2.0)) else 1.7,
            zorder=3 if (width == 10.0 and np.isclose(angle, np.pi / 2.0)) else 2,
            label=label + rf", $\sigma_k={item['kernel_sigma_ns']:.2f}$ ns",
        )
    extent = max(abs(min(centered_limits)), abs(max(centered_limits)))
    ax_a.axvline(0.0, color="#777777", ls=":", lw=0.9, zorder=1)
    ax_a.set_title("(a) Gaussian-control response kernels", fontsize=10)
    ax_a.set_xlabel(r"Time relative to kernel centroid $t-t_c$ (ns)")
    ax_a.set_ylabel(r"Normalised $|k(t)|$")
    ax_a.set_xlim(-extent, extent)
    ax_a.set_ylim(-0.03, 1.06)
    ax_a.legend(fontsize=6.4, loc="lower right", framealpha=0.82)
    ax_a.grid(True, ls=":", alpha=0.35, lw=0.55)

    # (b) Operational resolution against pulse width for each rotation angle.
    finite_dt = resolution["dt_min_ns"][np.isfinite(resolution["dt_min_ns"])]
    fail_y = float(np.max(finite_dt) + 0.55)
    for ai, label in enumerate(ANGLE_LABELS):
        values = resolution["dt_min_ns"][ai]
        ax_b.plot(
            PULSE_WIDTHS_NS, values,
            color=colors[ai], marker=markers[ai], lw=1.7, ms=5.2,
            label=rf"$\alpha={label.strip('$')}$",
        )
        failed = ~np.isfinite(values)
        if np.any(failed):
            ax_b.plot(
                PULSE_WIDTHS_NS[failed], np.full(np.sum(failed), fail_y),
                color=colors[ai], marker="x", ls="none", ms=6.5, mew=1.5,
                label="No pass through 40 ns" if ai == 0 else None,
            )
    ax_b.axhline(10.0, color="#555555", ls=":", lw=1.0,
                 label=r"Input-geometry floor $2\sigma_{\rm peak}$")
    ax_b.set_title("(b) Minimum resolvable peak spacing", fontsize=10)
    ax_b.set_xlabel(r"Pulse window $T_p$ (ns)")
    ax_b.set_ylabel(r"$\Delta t_{\min}$ (ns)")
    ax_b.set_xticks(PULSE_WIDTHS_NS)
    ax_b.set_ylim(9.7, fail_y + 0.65)
    ax_b.legend(fontsize=6.2, loc="lower center", ncol=2, framealpha=0.92)
    ax_b.grid(True, ls=":", alpha=0.35, lw=0.55)

    # (c) SNR at the middle shot budget, comparing rotation angles.
    shot_index = int(np.where(N_SHOTS == 10_000)[0][0])
    for ai, label in enumerate(ANGLE_LABELS):
        y = sensitivity["snr"][ai, shot_index]
        e = sensitivity["snr_sem"][ai, shot_index]
        good = y > 0
        ax_c.errorbar(
            AMPLITUDES[good], y[good], yerr=e[good], color=colors[ai],
            marker=markers[ai], lw=1.5, ms=4.2, capsize=2,
            label=rf"$\alpha={label.strip('$')}$",
        )
    ax_c.axhline(D.SNR_THRESHOLD, color="#555555", ls="--", lw=1.1,
                 label=rf"SNR$_{{\rm th}}={D.SNR_THRESHOLD:g}$")
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_title(r"(c) SNR at $N_{\rm shot}=10^4$", fontsize=10)
    ax_c.set_xlabel(r"Signal amplitude $A$ ($\Phi_0$)")
    ax_c.set_ylabel(r"SNR $=\mathbb{E}[\hat A]/\sigma_{\hat A,0}$")
    ax_c.legend(fontsize=6.7, loc="upper left", framealpha=0.92)
    ax_c.grid(True, which="both", ls=":", alpha=0.35, lw=0.55)

    # (d) Detection threshold against shot count for each rotation angle.
    for ai, label in enumerate(ANGLE_LABELS):
        ax_d.plot(
            N_SHOTS, sensitivity["a_min"][ai], color=colors[ai],
            marker=markers[ai], lw=1.7, ms=5.2,
            label=rf"$\alpha={label.strip('$')}$",
        )
    base = sensitivity["a_min"][-1, 0]
    guide = base * np.sqrt(N_SHOTS[0] / N_SHOTS)
    ax_d.plot(N_SHOTS, guide, color="#333333", ls=":", lw=1.2,
              label=r"$\propto N_{\rm shot}^{-1/2}$")
    ax_d.set_xscale("log")
    ax_d.set_yscale("log")
    ax_d.set_title(r"(d) Detection floor at $T_p=10$ ns", fontsize=10)
    ax_d.set_xlabel(r"Shots per point $N_{\rm shot}$")
    ax_d.set_ylabel(r"$A_{\min}$ ($\Phi_0$), SNR$\geq3$")
    ax_d.legend(fontsize=6.7, loc="lower left", framealpha=0.92)
    ax_d.grid(True, which="both", ls=":", alpha=0.35, lw=0.55)

    fig.suptitle("Pulse-parameter dependence of time resolution and sensitivity",
                 fontsize=11.5, y=0.985)
    fig.text(
        0.985, 0.96,
        rf"Gaussian envelope: $\sigma=T_p/5$; Wiener $\lambda_W={D.C.WIENER_LAMBDA:g}$",
        ha="right", va="top", fontsize=6.8, color="#555555",
    )
    fig.subplots_adjust(top=0.90, bottom=0.10, left=0.09, right=0.98,
                        hspace=0.36, wspace=0.28)
    png = os.path.join(HERE, STEM + ".png")
    pdf = os.path.join(HERE, STEM + ".pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def main():
    cache = os.path.join(HERE, STEM + ".npz")
    force = "--force" in sys.argv
    if os.path.exists(cache) and not force:
        z = dict(np.load(cache))
        kernels = {}
        for width in PULSE_WIDTHS_NS:
            for angle in ROTATION_ANGLES:
                key = kernel_key(width, angle)
                kernels[key] = {
                    "kernel": z[f"kernel_{key}"],
                    "kernel_time": z[f"kernel_time_{key}"],
                    "pulse_time": z[f"pulse_time_{key}"],
                    "control_angle_rad": float(z[f"control_angle_{key}"]),
                    "kernel_sigma_ns": float(
                        z["kernel_sigma_ns"][
                            np.where(np.isclose(ROTATION_ANGLES, angle))[0][0],
                            np.where(np.isclose(PULSE_WIDTHS_NS, width))[0][0],
                        ]
                    ),
                }
        resolution = {
            "valley_contrast": z["resolution_valley_contrast"],
            "nrmse": z["resolution_nrmse"],
            "peak_position_error_ns": z["resolution_peak_position_error_ns"],
            "resolved": z["resolution_resolved"],
            "dt_min_ns": z["dt_min_ns"],
            "kernel_sigma_ns": z["kernel_sigma_ns"],
        }
        sensitivity = {
            "template": z["sensitivity_template"],
            "snr": z["snr"],
            "snr_sem": z["snr_sem"],
            "sigma0": z["sigma_Ahat_zero"],
            "a_min": z["a_min"],
            "a_min_status": z["a_min_status"],
            "mean_estimate": z["mean_amplitude_estimate"],
        }
        data_path = cache
        meta_path = os.path.join(HERE, STEM + "_metadata.json")
        metrics_path = os.path.join(HERE, STEM + "_metrics.txt")
        print(f"loaded preview cache: {cache}")
    else:
        kernels = build_kernels()
        resolution = resolution_scan(kernels)
        sensitivity = sensitivity_scan(kernels)
        data_path, meta_path, metrics_path = save_data(
            kernels, resolution, sensitivity
        )

    png, pdf = plot(kernels, resolution, sensitivity)
    digest = hashlib.sha256(open(png, "rb").read()).hexdigest().upper()
    print(f"figure  -> {png}\n        -> {pdf}")
    print(f"data    -> {data_path}")
    print(f"metadata-> {meta_path}\nmetrics -> {metrics_path}")
    print(f"png sha256 -> {digest}")


if __name__ == "__main__":
    main()
