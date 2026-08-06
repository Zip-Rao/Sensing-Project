#!/usr/bin/env python3
"""Run P9 measurement and held-out predistortion benchmark."""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
import benchmark_core as B


def _json_ready(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def _bootstrap_ci(values: np.ndarray, seed: int = 17) -> list[float]:
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    draws = np.median(rng.choice(values, (4000, len(values)), replace=True), axis=1)
    return [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))]


def run() -> tuple[dict[str, np.ndarray], dict]:
    kernel = B.load_transient_kernel(REPO)
    cases = B.make_cases()
    protocols = ("transient", "cryoscope")
    n_case, n_seed = len(cases), len(B.SEEDS)
    global_rmse = np.full((2, n_case, n_seed), np.nan)
    edge_rmse = np.full_like(global_rmse, np.nan)
    bandwidth = np.full_like(global_rmse, np.nan)
    amp_error = np.full_like(global_rmse, np.nan)
    phase_error = np.full_like(global_rmse, np.nan)
    gains = np.full((2, n_seed), np.nan)
    main_index = next(i for i, c in enumerate(cases) if c.family == "D4")

    representative = {}
    compensation_rmse = np.full((2, n_seed, 2), np.nan)
    compensation_edge_rmse = np.full_like(compensation_rmse, np.nan)
    awg_peak = np.full((2, n_seed), np.nan)
    awg_slew = np.full_like(awg_peak, np.nan)
    noise_gain = np.full_like(awg_peak, np.nan)
    constraints_pass = np.zeros((2, n_seed), dtype=bool)
    uncorrected_rmse = np.full(2, np.nan)
    heldout = B.heldout_waveforms()
    resource_grid = np.array([4096, 16384, 65536, 262144, 1048576], dtype=int)
    resource_amp_error = np.full((2, len(resource_grid), n_seed), np.nan)
    resource_phase_error = np.full_like(resource_amp_error, np.nan)
    resource_qualified_bins = np.zeros_like(resource_amp_error, dtype=int)

    for pi, protocol in enumerate(protocols):
        for si, seed in enumerate(B.SEEDS):
            gains[pi, si] = B.calibrate_protocol_gain(protocol, kernel, int(seed + 70000))

    for ci, case in enumerate(cases):
        freq_ghz, h_true_f = B.true_frequency_response(case.impulse)
        for si, seed in enumerate(B.SEEDS):
            xs = [B.multisine(int(seed + 1000 * j), amp)
                  for j, amp in enumerate(B.TRAIN_AMPLITUDES)]
            y_truth = [B.apply_lti(x, case.impulse) for x in xs]
            for pi, protocol in enumerate(protocols):
                recovered = []
                raw_obs = []
                for j, y in enumerate(y_truth):
                    obs_seed = int(seed + 10000 * ci + 100 * pi + j)
                    if protocol == "transient":
                        obs, _ = B.transient_observe(y, kernel, obs_seed)
                        rec = B.transient_reconstruct(obs, kernel, gains[pi, si])
                    else:
                        p_i, p_q, _ = B.cryoscope_observe(y, obs_seed)
                        obs = np.stack([p_i, p_q])
                        rec = B.cryoscope_reconstruct(p_i, p_q, gains[pi, si])
                    recovered.append(rec)
                    raw_obs.append(obs)
                h_est, support, power = B.estimate_h(xs, recovered)
                qualified = B.qualified_mask(h_est, h_true_f, support)
                global_rmse[pi, ci, si] = B.rmse(recovered[1], y_truth[1])
                edge_rmse[pi, ci, si] = B.edge_rmse(recovered[1], y_truth[1])
                bandwidth[pi, ci, si] = B.effective_bandwidth_mhz(qualified, support)
                common = support & (freq_ghz <= 0.35)
                amp_error[pi, ci, si] = np.median(
                    np.abs(np.abs(h_est[common]) - np.abs(h_true_f[common]))
                    / np.maximum(np.abs(h_true_f[common]), 1e-9)
                )
                phase_error[pi, ci, si] = np.median(
                    np.abs(np.angle(h_est[common] * np.conj(h_true_f[common])))
                )

                if ci == main_index:
                    taps = B.design_inverse_fir(h_est, qualified, support)
                    for wi, target in enumerate(heldout.values()):
                        _, chip = B.compensate(target, case.impulse, taps)
                        compensation_rmse[pi, si, wi] = B.rmse(chip, target)
                        compensation_edge_rmse[pi, si, wi] = B.edge_rmse(chip, target)
                    awg, _ = B.compensate(next(iter(heldout.values())), case.impulse, taps)
                    hw = B.hardware_metrics(awg, taps)
                    awg_peak[pi, si] = hw["awg_peak"]
                    awg_slew[pi, si] = hw["awg_slew_phi0_per_ns"]
                    noise_gain[pi, si] = hw["filter_noise_gain"]
                    constraints_pass[pi, si] = hw["constraints_pass"]
                    if si == 0:
                        representative[protocol] = {
                            "x_train": xs[1], "y_truth": y_truth[1],
                            "y_reconstructed": recovered[1], "raw": raw_obs[1],
                            "h_est": h_est, "support": support,
                            "qualified": qualified, "taps": taps,
                        }

    main_case = cases[main_index]
    for wi, target in enumerate(heldout.values()):
        uncorrected_rmse[wi] = B.rmse(B.apply_lti(target, main_case.impulse), target)

    # Pre-registered accuracy-resource diagnostic on the frozen main case.
    freq_main, h_main = B.true_frequency_response(main_case.impulse)
    for ri, shots in enumerate(resource_grid):
        for si, seed in enumerate(B.SEEDS):
            xs = [B.multisine(int(seed + 3000 * j), amp)
                  for j, amp in enumerate(B.TRAIN_AMPLITUDES)]
            ys = [B.apply_lti(x, main_case.impulse) for x in xs]
            for pi, protocol in enumerate(protocols):
                gain = B.calibrate_protocol_gain(
                    protocol, kernel, int(seed + 90000 + ri), shots=int(shots)
                )
                recovered = []
                for j, y in enumerate(ys):
                    obs_seed = int(seed + 500000 * ri + 100 * pi + j)
                    if protocol == "transient":
                        obs, _ = B.transient_observe(y, kernel, obs_seed, int(shots))
                        rec = B.transient_reconstruct(obs, kernel, gain)
                    else:
                        p_i, p_q, _ = B.cryoscope_observe(y, obs_seed, int(shots))
                        rec = B.cryoscope_reconstruct(p_i, p_q, gain)
                    recovered.append(rec)
                h_est, support, _ = B.estimate_h(xs, recovered)
                qualified = B.qualified_mask(h_est, h_main, support)
                common = support & (freq_main <= 0.35)
                resource_amp_error[pi, ri, si] = np.median(
                    np.abs(np.abs(h_est[common]) - np.abs(h_main[common]))
                    / np.maximum(np.abs(h_main[common]), 1e-9)
                )
                resource_phase_error[pi, ri, si] = np.median(
                    np.abs(np.angle(h_est[common] * np.conj(h_main[common])))
                )
                resource_qualified_bins[pi, ri, si] = int(np.sum(qualified & support))

    arrays = {
        "t_ns": np.arange(B.N_SAMPLES) * B.DT_NS,
        "freq_mhz": np.fft.rfftfreq(B.N_SAMPLES, d=B.DT_NS) * 1000.0,
        "case_names": np.array([c.name for c in cases]),
        "case_families": np.array([c.family for c in cases]),
        "case_parameters": np.array([c.parameter for c in cases]),
        "seed_values": B.SEEDS,
        "protocol_names": np.array(protocols),
        "transient_kernel": kernel,
        "gain_values": gains,
        "measurement_rmse": global_rmse,
        "edge_rmse": edge_rmse,
        "qualified_bandwidth_mhz": bandwidth,
        "h_amplitude_error": amp_error,
        "h_phase_error_rad": phase_error,
        "compensation_rmse": compensation_rmse,
        "compensation_edge_rmse": compensation_edge_rmse,
        "heldout_names": np.array(list(heldout)),
        "uncorrected_rmse": uncorrected_rmse,
        "awg_peak": awg_peak,
        "awg_slew_phi0_per_ns": awg_slew,
        "filter_noise_gain": noise_gain,
        "constraints_pass": constraints_pass,
        "resource_grid_shots_per_arm": resource_grid,
        "resource_h_amplitude_error": resource_amp_error,
        "resource_h_phase_error_rad": resource_phase_error,
        "resource_qualified_bins": resource_qualified_bins,
        "main_impulse": main_case.impulse,
        "main_h_truth": B.true_frequency_response(main_case.impulse)[1],
    }
    for protocol, data in representative.items():
        for key, value in data.items():
            arrays[f"representative_{protocol}_{key}"] = value
    for name, target in heldout.items():
        arrays[f"heldout_{name}"] = target
        arrays[f"uncorrected_{name}"] = B.apply_lti(target, main_case.impulse)
        for protocol in protocols:
            taps = representative[protocol]["taps"]
            awg, chip = B.compensate(target, main_case.impulse, taps)
            arrays[f"representative_awg_{protocol}_{name}"] = awg
            arrays[f"representative_chip_{protocol}_{name}"] = chip

    short_mask = np.array([c.family in {"D1", "D2", "D3"} for c in cases])
    edge_improvement = 1.0 - edge_rmse[0, short_mask] / edge_rmse[1, short_mask]
    family_improvements = {}
    family_passes = 0
    for family in ("D1", "D2", "D3"):
        mask = np.array([c.family == family for c in cases])
        values = (1.0 - edge_rmse[0, mask] / edge_rmse[1, mask]).ravel()
        median = float(np.median(values))
        ci = _bootstrap_ci(values)
        passed = median >= 0.30 and ci[0] > 0.0
        family_passes += int(passed)
        family_improvements[family] = {"median": median, "bootstrap_95_ci": ci,
                                       "passed": passed}
    comp_ratio = 1.0 - compensation_rmse[0, :, 1] / compensation_rmse[1, :, 1]
    bw_t = float(np.median(bandwidth[0, main_index]))
    bw_c = float(np.median(bandwidth[1, main_index]))
    metrics = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "protocol-level noisy LTI benchmark; gate validation is separate",
        "resource": {
            "shots_per_arm": B.N_SHOTS, "arms_per_delay": 2,
            "samples_per_probe": B.N_SAMPLES, "training_probes": 3,
            "executions_per_training_set": 2 * B.N_SHOTS * B.N_SAMPLES * 3,
            "random_seeds": B.N_SEEDS,
        },
        "configuration": {
            "dt_ns": B.DT_NS, "fir_taps": B.FIR_TAPS,
            "fir_delay_samples": B.FIR_DELAY, "fir_ridge": B.FIR_RIDGE,
            "transient_wiener_lambda": B.TRANSIENT_LAMBDA,
            "cryoscope_sg_window": B.CRYO_SG_WINDOW,
            "main_distortion": "10 ns low-pass + 6 ns 4% echo + 200 ns 1% tail",
        },
        "G0": {
            "transient_identity_rmse_median": float(np.median(global_rmse[0, 0])),
            "cryoscope_identity_rmse_median": float(np.median(global_rmse[1, 0])),
            "passed": bool(np.median(global_rmse[:, 0]) < 0.003),
        },
        "G1": {"family_results": family_improvements,
               "families_passing": family_passes, "passed": family_passes >= 2},
        "G2": {
            "main_bandwidth_transient_mhz": bw_t,
            "main_bandwidth_cryoscope_mhz": bw_c,
            "bandwidth_ratio": bw_t / max(bw_c, 1e-12),
            "passed": bool(bw_t >= 1.5 * max(bw_c, 1e-12)),
        },
        "G3": {
            "flat_top_rmse_uncorrected": float(uncorrected_rmse[1]),
            "flat_top_rmse_transient_median": float(np.median(compensation_rmse[0, :, 1])),
            "flat_top_rmse_cryoscope_median": float(np.median(compensation_rmse[1, :, 1])),
            "transient_relative_advantage_median": float(np.median(comp_ratio)),
            "bootstrap_95_ci": _bootstrap_ci(comp_ratio),
            "hardware_constraints_all_pass": bool(np.all(constraints_pass)),
            "passed": bool(np.median(comp_ratio) >= 0.30
                           and _bootstrap_ci(comp_ratio)[0] > 0.0
                           and np.all(constraints_pass)),
        },
        "claim_gate": "G4 pending full-Hamiltonian CZ validation",
        "resource_diagnostic": {
            "shot_grid": resource_grid.tolist(),
            "qualified_bins_at_max_shots_transient_median": float(
                np.median(resource_qualified_bins[0, -1])
            ),
            "qualified_bins_at_max_shots_cryoscope_median": float(
                np.median(resource_qualified_bins[1, -1])
            ),
            "interpretation": (
                "Increasing shots reduces variance but does not establish a "
                "continuous qualified band under the frozen reconstruction rules."
            ),
        },
    }
    return arrays, metrics


def plot(arrays: dict[str, np.ndarray], metrics: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 7.2, "pdf.fonttype": 42, "axes.linewidth": 0.7})
    t = arrays["t_ns"]
    f = arrays["freq_mhz"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.7), constrained_layout=True)
    ax = axes.ravel()
    colors = {"transient": "#237a73", "cryoscope": "#2f6ea3"}
    ax[0].plot(t, 1e3 * arrays["representative_transient_x_train"], color="#999999",
               lw=0.8, label="AWG probe")
    ax[0].plot(t, 1e3 * arrays["representative_transient_y_truth"], color="black",
               lw=1.2, label="on-chip truth")
    for protocol in ("transient", "cryoscope"):
        ax[0].plot(t, 1e3 * arrays[f"representative_{protocol}_y_reconstructed"],
                   color=colors[protocol], ls="--" if protocol == "cryoscope" else "-.",
                   lw=1.0, label=protocol.capitalize())
    ax[0].set(xlabel="Time (ns)", ylabel=r"Flux ($10^{-3}\Phi_0$)",
              title="(a) Main mixed distortion")
    ax[0].legend(frameon=False, ncol=2)

    families = ("D1", "D2", "D3")
    positions = np.arange(3)
    for pi, protocol in enumerate(("transient", "cryoscope")):
        vals = []
        errs = []
        for family in families:
            mask = arrays["case_families"] == family
            samples = arrays["edge_rmse"][pi, mask].ravel()
            vals.append(np.median(samples))
            errs.append([np.median(samples) - np.quantile(samples, 0.16),
                         np.quantile(samples, 0.84) - np.median(samples)])
        ax[1].errorbar(positions + (pi - 0.5) * 0.08, vals,
                       yerr=np.array(errs).T, marker="o", capsize=2,
                       color=colors[protocol], label=protocol.capitalize())
    ax[1].set_xticks(positions, ["fast edge", "echo", "ringing"])
    ax[1].set(ylabel=r"Edge RMSE ($\Phi_0$)", title="(b) Short-time scan, 20 seeds")
    ax[1].set_yscale("log")
    ax[1].legend(frameon=False)

    truth = arrays["main_h_truth"]
    ax[2].plot(f, np.abs(truth), color="black", lw=1.2, label="truth")
    for protocol in ("transient", "cryoscope"):
        est = arrays[f"representative_{protocol}_h_est"]
        qual = arrays[f"representative_{protocol}_qualified"].astype(bool)
        ax[2].plot(f[qual], np.abs(est[qual]), color=colors[protocol], lw=1.0,
                   label=protocol.capitalize())
        ax[2].plot(f[~qual], np.abs(est[~qual]), color="#c9c9c9", lw=0.55)
    ax[2].set(xlim=(0, 400), xlabel="Frequency (MHz)", ylabel=r"$|H(f)|$",
              title="(c) Qualified transfer-function bins")
    ax[2].legend(frameon=False)

    target = arrays["heldout_short_flat_top"]
    ax[3].plot(t, target, color="black", lw=1.2, label="target")
    ax[3].plot(t, arrays["uncorrected_short_flat_top"], color="#c54b3c", lw=0.9,
               label="uncorrected")
    for protocol in ("transient", "cryoscope"):
        ax[3].plot(t, arrays[f"representative_chip_{protocol}_short_flat_top"],
                   color=colors[protocol], lw=1.0, label=f"{protocol}-FIR")
    ax[3].set(xlabel="Time (ns)", ylabel=r"On-chip flux ($\Phi_0$)",
              title="(d) Held-out flat-top compensation")
    ax[3].legend(frameon=False, ncol=2)
    for axis in ax:
        axis.grid(alpha=0.18, lw=0.45)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"fig_E1_short_time_identification.{suffix}",
                    dpi=400 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)


def write_environment(script_path: Path) -> None:
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO,
                                         text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"],
                                             cwd=REPO, text=True).strip())
    except Exception:
        commit, dirty = "unknown", True
    env = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "git_commit": commit, "git_dirty": dirty,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    (HERE / "environment.json").write_text(json.dumps(env, indent=2), encoding="utf-8")


def main() -> None:
    arrays, metrics = run()
    np.savez_compressed(HERE / "short_time_measurement_data.npz", **arrays)
    # This benchmark stores measurement and predistortion arrays together; the
    # second name keeps the plan's public artifact contract explicit.
    np.savez_compressed(HERE / "predistortion_validation_data.npz", **arrays)
    (HERE / "metrics.json").write_text(
        json.dumps(metrics, indent=2, default=_json_ready), encoding="utf-8"
    )
    plot(arrays, metrics)
    write_environment(Path(__file__))
    print(json.dumps({k: metrics[k] for k in ("G0", "G1", "G2", "G3")}, indent=2))


if __name__ == "__main__":
    main()
