#!/usr/bin/env python3
"""Generate the corrected frequency-deembedded benchmark through G3 only."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("SENSING_PROJECT_ROOT", HERE.parents[2])).resolve()
sys.path.insert(0, str(HERE))
import core as B


PROTOCOLS = ("transient", "cryoscope")
RESOURCE_GRID = np.array([4096, 16384, 65536, 131072, 262144, 524288, 1048576])
MAIN_LINE_DIAGNOSTIC_GRID = np.array([262144, 524288, 1048576])


def _bootstrap_ci(values: np.ndarray, seed: int = 901) -> list[float]:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if not len(values):
        return [float("nan"), float("nan")]
    rng = np.random.default_rng(seed)
    draws = np.median(rng.choice(values, (4000, len(values)), replace=True), axis=1)
    return [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))]


def _single_transfer(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    xf = np.fft.rfft(x)
    yf = np.fft.rfft(y)
    out = np.full_like(yf, np.nan + 1j * np.nan)
    support = np.abs(xf) >= 0.05 * np.max(np.abs(xf))
    out[support] = yf[support] / xf[support]
    out[0] = 1.0 + 0.0j
    return out


def _empirical_line_band(xs, reconstructed, response, reliable, d0_mask):
    estimates = []
    for x, y in zip(xs, reconstructed):
        measured = _single_transfer(x, y)
        estimates.append(B.deembed_transfer(measured, response, reliable))
    stack = np.stack(estimates)
    amp_cv = np.full(stack.shape[1], np.inf)
    phase_std = np.full(stack.shape[1], np.inf)
    for idx in np.flatnonzero(reliable):
        vals = stack[:, idx]
        vals = vals[np.isfinite(vals)]
        if len(vals):
            amp_cv[idx] = np.std(np.abs(vals)) / max(np.mean(np.abs(vals)), 1e-12)
            phase_std[idx] = np.std(np.unwrap(np.angle(vals)))
    stable = reliable & d0_mask & (amp_cv <= 0.05) & (phase_std <= 0.10)
    stable[0] = True
    return B.contiguous_mask(stable, reliable), amp_cv, phase_std


def _d0_validation(protocol, parameter, kernel, seed, shots, response, reliable):
    xs = B.probe_set(seed, split=1)
    reconstructed = []
    for j, x in enumerate(xs):
        _, rec = B.observe_and_reconstruct(
            protocol, x, kernel, int(seed + 5000 + j), shots, parameter
        )
        reconstructed.append(rec)
    measured, support, _ = B.estimate_transfer(xs, reconstructed)
    rel = reliable & support
    h_id = B.deembed_transfer(measured, response, rel)
    truth = np.ones_like(h_id)
    qualified = B.qualified_mask(h_id, truth, rel)
    return h_id, B.contiguous_mask(qualified, rel), xs, reconstructed


def _main_line_resource_diagnostic(protocol, parameter, kernel, case, case_index):
    """Test shot scaling on D4 without feeding the result into model selection."""
    bandwidth = np.full((len(MAIN_LINE_DIAGNOSTIC_GRID), len(B.SEEDS)), np.nan)
    qualified = np.zeros_like(bandwidth, dtype=bool)
    for ri, shots_np in enumerate(MAIN_LINE_DIAGNOSTIC_GRID):
        shots = int(shots_np)
        for si, seed_np in enumerate(B.SEEDS):
            seed = int(seed_np)
            response, reliable, *_ = B.calibrate_protocol_response(
                protocol, parameter, kernel, seed, shots, split=0
            )
            _, d0_mask, _, _ = _d0_validation(
                protocol, parameter, kernel, seed, shots, response, reliable
            )
            xs = B.probe_set(seed + 30000 * case_index, split=2)
            reconstructed = []
            for j, x in enumerate(xs):
                y = B.apply_lti_periodic(x, case.impulse)
                _, rec = B.observe_and_reconstruct(
                    protocol, y, kernel,
                    seed + 100000 * case_index + 1000 * PROTOCOLS.index(protocol) + j,
                    shots, parameter
                )
                reconstructed.append(rec)
            _, support, _ = B.estimate_transfer(xs, reconstructed)
            rel = reliable & support
            design_band, _, _ = _empirical_line_band(
                xs, reconstructed, response, rel, d0_mask
            )
            bandwidth[ri, si] = B.bandwidth_mhz(design_band)
            qualified[ri, si] = np.sum(design_band[1:]) >= B.MIN_QUALIFIED_BINS
    return bandwidth, qualified


def run():
    kernel = B.load_transient_kernel(REPO)
    cases = B.make_cases()
    main_index = next(i for i, case in enumerate(cases) if case.family == "D4")
    choices = {
        p: B.choose_protocol_parameter(p, kernel, B.MAIN_SHOTS) for p in PROTOCOLS
    }

    n_p, n_c, n_s = len(PROTOCOLS), len(cases), len(B.SEEDS)
    measurement_rmse = np.full((n_p, n_c, n_s), np.nan)
    edge_rmse = np.full_like(measurement_rmse, np.nan)
    bandwidth = np.full_like(measurement_rmse, np.nan)
    h_amp_error = np.full_like(measurement_rmse, np.nan)
    h_phase_error = np.full_like(measurement_rmse, np.nan)
    design_bandwidth = np.full_like(measurement_rmse, np.nan)
    d0_bandwidth = np.full((n_p, n_s), np.nan)
    d0_error = np.full((n_p, n_s), np.nan)

    validation_rmse = np.full((n_p, n_s), np.nan)
    test_rmse = np.full((n_p, n_s), np.nan)
    filter_taps_count = np.full((n_p, n_s), np.nan)
    filter_ridge = np.full((n_p, n_s), np.nan)
    filter_noise_gain = np.full((n_p, n_s), np.nan)
    awg_peak = np.full((n_p, n_s), np.nan)
    awg_slew = np.full((n_p, n_s), np.nan)
    filter_valid = np.zeros((n_p, n_s), dtype=bool)
    # 0=valid, 1=insufficient empirical band, 2=no constraint-satisfying FIR
    filter_failure_code = np.full((n_p, n_s), 1, dtype=int)

    heldout = B.heldout_waveforms()
    validation = heldout["validation_dual_pulse"]
    test = heldout["test_short_flat_top"]
    main_impulse = cases[main_index].impulse
    uncorrected_validation = B.apply_lti(validation, main_impulse)
    uncorrected_test = B.apply_lti(test, main_impulse)
    uncorrected_validation_rmse = B.rmse(uncorrected_validation, validation)
    uncorrected_test_rmse = B.rmse(uncorrected_test, test)
    representative = {}
    filter_representative = {}

    for pi, protocol in enumerate(PROTOCOLS):
        parameter = choices[protocol].parameter
        for si, seed_np in enumerate(B.SEEDS):
            seed = int(seed_np)
            response, reliable, cal_xs, cal_rec, cal_raw, _ = (
                B.calibrate_protocol_response(
                    protocol, parameter, kernel, seed, B.MAIN_SHOTS, split=0
                )
            )
            h_id, d0_mask, _, d0_rec = _d0_validation(
                protocol, parameter, kernel, seed, B.MAIN_SHOTS,
                response, reliable
            )
            d0_bandwidth[pi, si] = B.bandwidth_mhz(d0_mask)
            d0_error[pi, si] = B.transfer_error(
                h_id, np.ones_like(h_id), reliable
            )

            for ci, case in enumerate(cases):
                xs = B.probe_set(seed + 30000 * ci, split=2)
                y_truth = [B.apply_lti_periodic(x, case.impulse) for x in xs]
                reconstructed, raw = [], []
                for j, y in enumerate(y_truth):
                    obs, rec = B.observe_and_reconstruct(
                        protocol, y, kernel,
                        seed + 100000 * ci + 1000 * pi + j,
                        B.MAIN_SHOTS, parameter
                    )
                    raw.append(obs)
                    reconstructed.append(rec)
                measured, support, _ = B.estimate_transfer(xs, reconstructed)
                rel = reliable & support
                h_est = B.deembed_transfer(measured, response, rel)
                h_true = B.truth_frequency_response(case.impulse)
                eval_qualified = B.qualified_mask(h_est, h_true, rel)
                eval_band = B.contiguous_mask(eval_qualified, rel)
                design_band, amp_cv, phase_std = _empirical_line_band(
                    xs, reconstructed, response, rel, d0_mask
                )

                freq_mhz = np.fft.rfftfreq(B.N_SAMPLES, d=B.DT_NS) * 1000.0
                common_wave_band = freq_mhz <= B.COMMON_WAVEFORM_BAND_MHZ
                common_wave_band[0] = True
                common_wave_band &= rel
                calibrated_wave = B.deembed_waveform(
                    reconstructed[1], response, common_wave_band
                )
                truth_bandlimited = B.deembed_waveform(
                    y_truth[1], np.ones_like(response), common_wave_band
                )
                measurement_rmse[pi, ci, si] = B.rmse(
                    calibrated_wave, truth_bandlimited
                )
                edge_rmse[pi, ci, si] = B.edge_rmse(
                    calibrated_wave, truth_bandlimited
                )
                bandwidth[pi, ci, si] = B.bandwidth_mhz(eval_band)
                design_bandwidth[pi, ci, si] = B.bandwidth_mhz(design_band)
                valid = rel & np.isfinite(h_est)
                valid[0] = False
                h_amp_error[pi, ci, si] = np.median(
                    np.abs(np.abs(h_est[valid]) - np.abs(h_true[valid]))
                    / np.maximum(np.abs(h_true[valid]), 1e-12)
                )
                h_phase_error[pi, ci, si] = np.median(
                    np.abs(np.angle(h_est[valid] * np.conj(h_true[valid])))
                )

                if ci == main_index:
                    selected = B.choose_filter(
                        h_est, design_band, case.impulse, validation
                    )
                    if selected is not None:
                        filter_valid[pi, si] = True
                        filter_failure_code[pi, si] = 0
                        validation_rmse[pi, si] = selected.validation_rmse
                        filter_taps_count[pi, si] = selected.n_taps
                        filter_ridge[pi, si] = selected.ridge
                        filter_noise_gain[pi, si] = selected.noise_gain
                        awg_peak[pi, si] = selected.awg_peak
                        awg_slew[pi, si] = selected.awg_slew
                        _, chip_test = B.compensate(
                            test, case.impulse, selected.taps
                        )
                        test_rmse[pi, si] = B.rmse(chip_test, test)
                    elif np.sum(design_band[1:]) >= B.MIN_QUALIFIED_BINS:
                        filter_failure_code[pi, si] = 2
                    if si == 0:
                        representative[protocol] = {
                            "response": response, "reliable": reliable,
                            "d0_mask": d0_mask, "x_train": xs[1],
                            "chip_truth_train": y_truth[1],
                            "raw_observation": raw[1],
                            "raw_reconstruction": reconstructed[1],
                            "calibrated_waveform": calibrated_wave,
                            "h_est": h_est, "h_true": h_true,
                            "eval_band": eval_band, "design_band": design_band,
                            "amp_cv": amp_cv, "phase_std": phase_std,
                        }
                    if selected is not None and protocol not in filter_representative:
                        filter_representative[protocol] = {
                            "seed": np.array([seed]), "filter_taps": selected.taps,
                        }

    resource_tuning_bandwidth = np.full((len(PROTOCOLS), len(RESOURCE_GRID)), np.nan)
    resource_tuning_parameter = np.full_like(resource_tuning_bandwidth, np.nan)
    for ri, shots in enumerate(RESOURCE_GRID):
        for pi, protocol in enumerate(PROTOCOLS):
            choice = B.choose_protocol_parameter(protocol, kernel, int(shots))
            resource_tuning_bandwidth[pi, ri] = choice.validation_bandwidth_mhz
            resource_tuning_parameter[pi, ri] = choice.parameter

    main_line_resource_bandwidth = np.full(
        (len(PROTOCOLS), len(MAIN_LINE_DIAGNOSTIC_GRID), n_s), np.nan
    )
    main_line_resource_qualified = np.zeros_like(
        main_line_resource_bandwidth, dtype=bool
    )
    for pi, protocol in enumerate(PROTOCOLS):
        diag_bandwidth, diag_qualified = _main_line_resource_diagnostic(
            protocol, choices[protocol].parameter, kernel,
            cases[main_index], main_index
        )
        main_line_resource_bandwidth[pi] = diag_bandwidth
        main_line_resource_qualified[pi] = diag_qualified

    arrays = {
        "t_ns": np.arange(B.N_SAMPLES) * B.DT_NS,
        "freq_mhz": np.fft.rfftfreq(B.N_SAMPLES, d=B.DT_NS) * 1000.0,
        "protocol_names": np.array(PROTOCOLS),
        "case_names": np.array([c.name for c in cases]),
        "case_families": np.array([c.family for c in cases]),
        "case_parameters": np.array([c.parameter for c in cases]),
        "seed_values": B.SEEDS,
        "transient_kernel": kernel,
        "protocol_parameters": np.array([choices[p].parameter for p in PROTOCOLS]),
        "d0_bandwidth_mhz": d0_bandwidth,
        "d0_transfer_error": d0_error,
        "measurement_rmse": measurement_rmse,
        "edge_rmse": edge_rmse,
        "qualified_bandwidth_mhz": bandwidth,
        "design_bandwidth_mhz": design_bandwidth,
        "h_amplitude_error": h_amp_error,
        "h_phase_error_rad": h_phase_error,
        "validation_rmse": validation_rmse,
        "test_rmse": test_rmse,
        "filter_valid": filter_valid,
        "filter_failure_code": filter_failure_code,
        "filter_taps_count": filter_taps_count,
        "filter_ridge": filter_ridge,
        "filter_noise_gain": filter_noise_gain,
        "awg_peak": awg_peak,
        "awg_slew_phi0_per_ns": awg_slew,
        "resource_grid_shots_per_arm": RESOURCE_GRID,
        "resource_tuning_bandwidth_mhz": resource_tuning_bandwidth,
        "resource_tuning_parameter": resource_tuning_parameter,
        "main_line_diagnostic_grid_shots_per_arm": MAIN_LINE_DIAGNOSTIC_GRID,
        "main_line_diagnostic_design_bandwidth_mhz": main_line_resource_bandwidth,
        "main_line_diagnostic_qualified": main_line_resource_qualified,
        "validation_target": validation,
        "test_target": test,
        "validation_uncorrected": uncorrected_validation,
        "test_uncorrected": uncorrected_test,
        "main_impulse": main_impulse,
        "main_h_truth": B.truth_frequency_response(main_impulse),
    }
    for protocol, data in representative.items():
        for key, value in data.items():
            arrays[f"representative_{protocol}_{key}"] = value
        taps = filter_representative.get(protocol, {}).get("filter_taps", np.array([]))
        if len(taps):
            arrays[f"filter_representative_{protocol}_seed"] = (
                filter_representative[protocol]["seed"]
            )
            arrays[f"filter_representative_{protocol}_taps"] = taps
            awg_v, chip_v = B.compensate(validation, main_impulse, taps)
            awg_t, chip_t = B.compensate(test, main_impulse, taps)
            arrays[f"representative_{protocol}_validation_awg"] = awg_v
            arrays[f"representative_{protocol}_validation_chip"] = chip_v
            arrays[f"representative_{protocol}_test_awg"] = awg_t
            arrays[f"representative_{protocol}_test_chip"] = chip_t

    family_results = {}
    family_passes = 0
    for family in ("D1", "D2", "D3"):
        mask = np.array([c.family == family for c in cases])
        improvement = (1.0 - edge_rmse[0, mask] / edge_rmse[1, mask]).ravel()
        median = float(np.nanmedian(improvement))
        ci = _bootstrap_ci(improvement)
        passed = bool(median >= 0.30 and ci[0] > 0.0)
        family_passes += int(passed)
        family_results[family] = {
            "median_relative_improvement": median,
            "bootstrap_95_ci": ci, "passed": passed,
        }

    valid_both = filter_valid[0] & filter_valid[1]
    test_advantage = np.full(n_s, np.nan)
    test_advantage[valid_both] = (
        1.0 - test_rmse[0, valid_both] / test_rmse[1, valid_both]
    )
    paired_values = test_advantage[np.isfinite(test_advantage)]
    g3_ci = _bootstrap_ci(paired_values) if len(paired_values) else [None, None]
    g3_median = float(np.median(paired_values)) if len(paired_values) else None
    bw_t = float(np.median(bandwidth[0, main_index]))
    bw_c = float(np.median(bandwidth[1, main_index]))
    g1_pass = family_passes >= 2
    g2_pass = bool(bw_t >= 1.5 * max(bw_c, 1e-12))
    g3_pass = bool(
        np.all(filter_valid)
        and g3_median is not None and g3_median >= 0.30
        and g3_ci[0] is not None and g3_ci[0] > 0.0
    )
    metrics = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "frequency-deembedded protocol-level LTI benchmark through G3",
        "resource": {
            "shots_per_arm": B.MAIN_SHOTS, "arms_per_delay": 2,
            "samples_per_probe": B.N_SAMPLES, "probes_per_split": B.N_PROBES,
            "calibration_splits": 2, "random_test_seeds": B.N_SEEDS,
            "main_resource_selection": (
                "lowest D0-only grid point with >=8 validated continuous "
                "non-zero bins for both protocols"
            ),
            "main_line_resource_diagnostic": {
                "role": (
                    "post-selection D4 diagnostic only; excluded from resource, "
                    "protocol-parameter, and FIR selection"
                ),
                "shots_per_arm": MAIN_LINE_DIAGNOSTIC_GRID.tolist(),
                "frozen_protocol_parameters": {
                    p: float(choices[p].parameter) for p in PROTOCOLS
                },
                "median_design_bandwidth_mhz": {
                    p: np.median(main_line_resource_bandwidth[pi], axis=1).tolist()
                    for pi, p in enumerate(PROTOCOLS)
                },
                "fraction_with_minimum_continuous_band": {
                    p: np.mean(main_line_resource_qualified[pi], axis=1).tolist()
                    for pi, p in enumerate(PROTOCOLS)
                },
            },
        },
        "protocol_choices": {
            p: {
                "parameter": choices[p].parameter,
                "d0_tuning_bandwidth_mhz": choices[p].validation_bandwidth_mhz,
                "d0_tuning_error": choices[p].validation_error,
            } for p in PROTOCOLS
        },
        "G0": {
            "transient_d0_bandwidth_median_mhz": float(np.median(d0_bandwidth[0])),
            "cryoscope_d0_bandwidth_median_mhz": float(np.median(d0_bandwidth[1])),
            "transient_d0_error_median": float(np.median(d0_error[0])),
            "cryoscope_d0_error_median": float(np.median(d0_error[1])),
            "passed": bool(np.median(d0_bandwidth, axis=1).min() >= 31.25),
        },
        "G1": {"family_results": family_results,
               "families_passing": family_passes, "passed": g1_pass},
        "G2": {
            "main_bandwidth_transient_mhz": bw_t,
            "main_bandwidth_cryoscope_mhz": bw_c,
            "bandwidth_ratio": bw_t / max(bw_c, 1e-12),
            "passed": g2_pass,
        },
        "G3": {
            "uncorrected_validation_rmse": uncorrected_validation_rmse,
            "uncorrected_test_rmse": uncorrected_test_rmse,
            "valid_filters_transient": int(np.sum(filter_valid[0])),
            "valid_filters_cryoscope": int(np.sum(filter_valid[1])),
            "paired_valid_seeds": int(np.sum(valid_both)),
            "failure_counts": {
                p: {
                    "insufficient_band": int(np.sum(filter_failure_code[pi] == 1)),
                    "hardware_constraints": int(np.sum(filter_failure_code[pi] == 2)),
                } for pi, p in enumerate(PROTOCOLS)
            },
            "test_rmse_transient_median": float(np.nanmedian(test_rmse[0])),
            "test_rmse_cryoscope_median": float(np.nanmedian(test_rmse[1])),
            "transient_relative_advantage_median": g3_median,
            "bootstrap_95_ci": g3_ci,
            "passed": g3_pass,
        },
        "cz_entry": {
            "eligible": bool(g1_pass and g2_pass and g3_pass),
            "reason": (
                "G1-G3 did not all pass; CZ simulation and plotting were not entered."
            ),
        },
    }
    return arrays, metrics


def _style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7.0, "axes.labelsize": 7.0,
        "axes.titlesize": 7.5, "legend.fontsize": 6.1,
        "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.65,
    })


def plot_identification(d):
    import matplotlib.pyplot as plt
    colors = {"transient": "#23877b", "cryoscope": "#276b9a"}
    t, f = d["t_ns"], d["freq_mhz"]
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 4.6), constrained_layout=True)
    ax = ax.ravel()
    ax[0].plot(t, 1e3*d["representative_transient_chip_truth_train"],
               color="black", label="on-chip truth")
    for p in PROTOCOLS:
        ax[0].plot(t, 1e3*d[f"representative_{p}_calibrated_waveform"],
                   color=colors[p], label=p.capitalize(), alpha=0.9)
    ax[0].set(title="(a) Frequency-calibrated reconstruction", xlabel="Time (ns)",
              ylabel=r"Flux ($10^{-3}\Phi_0$)")
    ax[0].legend(frameon=False)

    for p in PROTOCOLS:
        r = d[f"representative_{p}_response"]
        reliable = d[f"representative_{p}_reliable"].astype(bool)
        ax[1].plot(f, np.abs(r), color=colors[p], label=p.capitalize())
        ax[1].plot(f[~reliable], np.abs(r[~reliable]), color="#cccccc", lw=0.6)
    ax[1].set(title="(b) Calibrated protocol response $R(f)$", xlabel="Frequency (MHz)",
              ylabel=r"$|R(f)|$", xlim=(0, 400))
    ax[1].legend(frameon=False)

    truth = d["main_h_truth"]
    ax[2].plot(f, np.abs(truth), color="black", label="line truth")
    for p in PROTOCOLS:
        h = d[f"representative_{p}_h_est"]
        band = d[f"representative_{p}_eval_band"].astype(bool)
        ax[2].plot(f[band], np.abs(h[band]), color=colors[p], label=p.capitalize())
        ax[2].plot(f[~band], np.abs(h[~band]), color="#cccccc", lw=0.55)
    ax[2].set(title="(c) Deembedded line transfer $H(f)$", xlabel="Frequency (MHz)",
              ylabel=r"$|H(f)|$", xlim=(0, 400))
    ax[2].legend(frameon=False)

    for pi, p in enumerate(PROTOCOLS):
        med = np.median(d["qualified_bandwidth_mhz"][pi], axis=1)
        ax[3].plot(np.arange(len(med)), med, marker="o", ms=3,
                   color=colors[p], label=p.capitalize())
    ax[3].set_xticks(np.arange(len(d["case_names"])), d["case_families"], rotation=45)
    ax[3].set(title="(d) Qualified bandwidth across distortion cases",
              xlabel="Distortion case", ylabel="Bandwidth (MHz)")
    ax[3].legend(frameon=False)
    for a in ax: a.grid(alpha=0.18, lw=0.45)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"fig_E1_frequency_deembedded_identification.{ext}",
                    dpi=400 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def plot_predistortion(d):
    import matplotlib.pyplot as plt
    colors = {"transient": "#23877b", "cryoscope": "#276b9a"}
    t = d["t_ns"]
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 4.6), constrained_layout=True)
    ax = ax.ravel()
    ax[0].plot(t, d["test_target"], color="black", label="target")
    ax[0].plot(t, d["test_uncorrected"], color="#c23b33", label="uncorrected")
    for p in PROTOCOLS:
        key = f"representative_{p}_test_chip"
        if key in d:
            ax[0].plot(t, d[key], color=colors[p], label=f"{p}-FIR")
    ax[0].set(title="(a) Frozen-filter test waveform", xlabel="Time (ns)",
              ylabel=r"On-chip flux ($\Phi_0$)")
    ax[0].legend(frameon=False, ncol=2)

    for p in PROTOCOLS:
        key = f"representative_{p}_test_chip"
        if key in d:
            ax[1].plot(t, d[key]-d["test_target"], color=colors[p], label=p.capitalize())
    ax[1].plot(t, d["test_uncorrected"]-d["test_target"], color="#c23b33",
               label="Uncorrected")
    ax[1].set(title="(b) Test residual", xlabel="Time (ns)",
              ylabel=r"Residual ($\Phi_0$)")
    ax[1].legend(frameon=False)

    x = np.arange(len(B.SEEDS))
    for pi, p in enumerate(PROTOCOLS):
        ax[2].scatter(x, d["test_rmse"][pi], s=11, color=colors[p], label=p.capitalize())
    ax[2].axhline(B.rmse(d["test_uncorrected"], d["test_target"]),
                  color="#c23b33", ls="--", label="Uncorrected")
    ax[2].set(title="(c) Independent test RMSE, 20 seeds", xlabel="Seed index",
              ylabel=r"RMSE ($\Phi_0$)", yscale="log")
    ax[2].legend(frameon=False)

    for pi, p in enumerate(PROTOCOLS):
        values = d["main_line_diagnostic_design_bandwidth_mhz"][pi]
        med = np.median(values, axis=1)
        lo, hi = np.quantile(values, [0.25, 0.75], axis=1)
        ax[3].plot(d["main_line_diagnostic_grid_shots_per_arm"], med, marker="o",
                   color=colors[p], label=p.capitalize())
        ax[3].fill_between(d["main_line_diagnostic_grid_shots_per_arm"], lo, hi,
                           color=colors[p], alpha=0.15, linewidth=0)
    ax[3].axvline(B.MAIN_SHOTS, color="black", ls=":", label="Frozen resource")
    ax[3].set(title="(d) D4 design-band resource diagnostic", xlabel="Shots per arm",
              ylabel="Design bandwidth (MHz)", xscale="log")
    ax[3].legend(frameon=False)
    for a in ax: a.grid(alpha=0.18, lw=0.45)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"fig_E2_pre_cz_predistortion_validation.{ext}",
                    dpi=400 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def write_environment():
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO,
                                         text=True, stderr=subprocess.DEVNULL).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"],
                                             cwd=REPO, text=True,
                                             stderr=subprocess.DEVNULL).strip())
    except Exception:
        commit, dirty = "unknown", True
    script = Path(__file__)
    env = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "numpy": np.__version__,
        "platform": platform.platform(), "git_commit": commit,
        "git_dirty": dirty,
        "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
    }
    (HERE / "environment.json").write_text(json.dumps(env, indent=2), encoding="utf-8")


def main():
    arrays, metrics = run()
    np.savez_compressed(HERE / "pre_cz_measurement_data.npz", **arrays)
    np.savez_compressed(HERE / "pre_cz_predistortion_data.npz", **arrays)
    (HERE / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    (HERE / "cz_entry_status.json").write_text(
        json.dumps(metrics["cz_entry"], indent=2), encoding="utf-8"
    )
    _style()
    plot_identification(arrays)
    plot_predistortion(arrays)
    write_environment()
    print(json.dumps({k: metrics[k] for k in ("protocol_choices", "G0", "G1", "G2", "G3", "cz_entry")}, indent=2))


if __name__ == "__main__":
    main()
