#!/usr/bin/env python3
"""Generate final Fig. 2(d): local accuracy versus instrumented solver cost."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from qutip import QobjEvo, basis, mesolve

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sqc.calibration.frequency import _fft_peak
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.devices.transmon import TransmonQubit
from fig2_simulation_common import (
    COST_CSV,
    COST_DATA,
    ESTIMATOR_DATA,
    RESPONSE_DATA,
    SUMMARY_DATA,
    load_npz,
    save_npz,
)


TWO_PI = 2.0 * np.pi
HERE = Path(__file__).resolve().parent
SOURCE_A = RESPONSE_DATA
SOURCE_B = ESTIMATOR_DATA
SOURCE_C = SUMMARY_DATA
STEM = "fig2d-accuracy-cost"
STYLE = ROOT.parent / "scholaraio" / ".claude" / "skills" / "draw" / "nature_pub.mplstyle"
TAU_MAX_NS = 200.0
ARTIFICIAL_GHZ = 0.05
RAMSEY_SETTINGS = (
    ("Single sweep, 2 ns", 2.0, False),
    ("Single sweep, 1 ns", 1.0, False),
    ("Single sweep, 0.5 ns", 0.5, False),
    ("Double sweep, 0.5 ns", 0.5, True),
)


def make_qubit(a: dict[str, np.ndarray]) -> TransmonQubit:
    return TransmonQubit(
        EC=TWO_PI * float(a["EC_over_h_GHz"]),
        EJ=TWO_PI * float(a["EJ_over_h_GHz"]),
        T1=float(a["T1_ns"]), T2=float(a["T2_ns"]),
        flux=0.0, state=0, n_levels=2,
    )


def ramsey_sweep(
    delta_mhz: float,
    tau_list: np.ndarray,
    artificial_ghz: float,
    a: dict[str, np.ndarray],
) -> tuple[np.ndarray, int]:
    """Run one dissipative Ramsey sweep and return populations plus call count."""
    qubit = make_qubit(a)
    omega_d = float(a["omega_d_rad_GHz"])
    t_rabi = np.asarray(a["t_rabi_ns"], dtype=float)
    envelope = np.asarray(a["envelope_samples"], dtype=float)
    dt = float(a["dt_ns"])
    t_global = np.arange(0.0, 2.0 * float(a["pulse_window_ns"]) + TAU_MAX_NS + dt, dt)
    flux = FluxSignal(type=0, t_list=t_global)
    qubit.qubit_in_mag(flux, frame=1, omega_d=omega_d)
    target_omega_q = omega_d + TWO_PI * float(delta_mhz) * 1e-3
    detuning_shift = target_omega_q - qubit.frequency
    h_base = QobjEvo(qubit.H_list, tlist=t_global, order=1) + detuning_shift * qubit.n
    projector_e = basis(2, 1) * basis(2, 1).dag()
    populations = np.zeros(len(tau_list), dtype=float)
    calls = 0
    for idx, tau in enumerate(tau_list):
        phase2 = TWO_PI * artificial_ghz * float(tau)
        ctrl = create_ramsey_pulse(
            t_rabi, tau=float(tau), omega_d=omega_d,
            phase1=0.0, phase2=phase2, qubit=qubit,
            rotation_angle=float(a["rotation_angle_rad"]), envelope=envelope,
        )
        h_total = h_base + QobjEvo(ctrl.hamiltonian_on(t_global), tlist=t_global, order=1)
        result = mesolve(
            h_total, qubit.state, t_global, qubit.get_collapse_operators(),
            e_ops=[projector_e],
            options={"atol": 1e-10, "rtol": 1e-8, "max_step": dt / 2.0},
        )
        populations[idx] = float(result.expect[0][-1])
        calls += 1
    return populations, calls


def ramsey_estimate(
    delta_mhz: float,
    spacing_ns: float,
    double_sweep: bool,
    a: dict[str, np.ndarray],
) -> tuple[float, int]:
    tau = np.arange(0.0, TAU_MAX_NS, spacing_ns)
    plus, calls = ramsey_sweep(delta_mhz, tau, +ARTIFICIAL_GHZ, a)
    f_plus = _fft_peak(plus, spacing_ns)
    if f_plus is None:
        raise RuntimeError("Ramsey + artificial-detuning sweep has insufficient contrast")
    if not double_sweep:
        estimate_ghz = ARTIFICIAL_GHZ - f_plus
        return estimate_ghz * 1e3, calls
    minus, minus_calls = ramsey_sweep(delta_mhz, tau, -ARTIFICIAL_GHZ, a)
    f_minus = _fft_peak(minus, spacing_ns)
    if f_minus is None:
        raise RuntimeError("Ramsey - artificial-detuning sweep has insufficient contrast")
    estimate_ghz = (f_minus**2 - f_plus**2) / (4.0 * ARTIFICIAL_GHZ)
    return estimate_ghz * 1e3, calls + minus_calls


def compute() -> dict[str, np.ndarray]:
    a = load_npz(SOURCE_A)
    b = load_npz(SOURCE_B)
    c = load_npz(SOURCE_C)
    delta = np.asarray(c["validation_delta_mhz"], dtype=float)
    estimates = np.zeros((len(RAMSEY_SETTINGS), len(delta)))
    calls = np.zeros_like(estimates, dtype=int)
    for col, truth in enumerate(delta):
        print(f"running Delta={truth:+.0f} MHz", flush=True)
        plus_peaks: dict[float, float] = {}
        for spacing in (2.0, 1.0, 0.5):
            tau = np.arange(0.0, TAU_MAX_NS, spacing)
            populations, _ = ramsey_sweep(truth, tau, +ARTIFICIAL_GHZ, a)
            peak = _fft_peak(populations, spacing)
            if peak is None:
                raise RuntimeError("Ramsey + artificial-detuning sweep has insufficient contrast")
            plus_peaks[spacing] = peak
        tau_fine = np.arange(0.0, TAU_MAX_NS, 0.5)
        minus_populations, _ = ramsey_sweep(truth, tau_fine, -ARTIFICIAL_GHZ, a)
        minus_peak = _fft_peak(minus_populations, 0.5)
        if minus_peak is None:
            raise RuntimeError("Ramsey - artificial-detuning sweep has insufficient contrast")
        for row, (label, spacing, double_sweep) in enumerate(RAMSEY_SETTINGS):
            if double_sweep:
                estimate_ghz = (minus_peak**2 - plus_peaks[spacing]**2) / (4.0 * ARTIFICIAL_GHZ)
                calls[row, col] = 2 * len(np.arange(0.0, TAU_MAX_NS, spacing))
            else:
                estimate_ghz = ARTIFICIAL_GHZ - plus_peaks[spacing]
                calls[row, col] = len(np.arange(0.0, TAU_MAX_NS, spacing))
            estimates[row, col] = estimate_ghz * 1e3
            print(f"  {label}: estimate={estimates[row, col]:+.6f} calls={calls[row, col]}", flush=True)
    errors = estimates - delta[None, :]
    ramsey_mae = np.mean(np.abs(errors), axis=1)
    ramsey_calls = np.max(calls, axis=1).astype(float)

    kernel_mae = float(c["mae_mhz"][2])
    n_validation = int(c["validation_count"])
    steady_calls = 2.0
    response_calibration_calls = 2.0 * len(a["calibration_delta_mhz"])
    response_validation_calls = 2.0 * len(a["validation_delta_mhz"])
    kernel_calibration_calls = 2.0
    qualification_calls = (
        response_calibration_calls + response_validation_calls + kernel_calibration_calls
    )
    amortized_calls = steady_calls + qualification_calls / n_validation
    cold_start_calls = steady_calls + qualification_calls
    return {
        "validation_delta_mhz": delta,
        "ramsey_labels": np.array([item[0] for item in RAMSEY_SETTINGS]),
        "ramsey_spacing_ns": np.array([item[1] for item in RAMSEY_SETTINGS]),
        "ramsey_double_sweep": np.array([item[2] for item in RAMSEY_SETTINGS]),
        "ramsey_estimate_mhz": estimates,
        "ramsey_error_mhz": errors,
        "ramsey_mae_mhz": ramsey_mae,
        "ramsey_solver_calls_per_measurement": ramsey_calls,
        "ramsey_solver_calls_by_point": calls,
        "short_pulse_mae_mhz": np.array(kernel_mae),
        "short_pulse_steady_calls": np.array(steady_calls),
        "short_pulse_response_calibration_calls": np.array(response_calibration_calls),
        "short_pulse_response_validation_calls": np.array(response_validation_calls),
        "short_pulse_kernel_calibration_calls": np.array(kernel_calibration_calls),
        "short_pulse_qualification_calls": np.array(qualification_calls),
        "short_pulse_amortized_calls": np.array(amortized_calls),
        "short_pulse_cold_start_calls": np.array(cold_start_calls),
        "short_pulse_amortization_count": np.array(n_validation),
        "delta_val_mhz": c["delta_val_mhz"],
        "sigma": c["sigma"],
        "cost_unit": np.array("counted_solver_calls_mesolve_plus_sesolve"),
        "truth_role": c["truth_role"],
        "validation_role": c["validation_role"],
        "source_a": np.array(SOURCE_A.name),
        "source_b": np.array(SOURCE_B.name),
        "source_c": np.array(SOURCE_C.name),
    }


def refresh_short_pulse_costs(data: dict[str, np.ndarray]) -> None:
    """Recompute cost roles from frozen response and validation metadata."""
    a = load_npz(SOURCE_A)
    c = load_npz(SOURCE_C)
    steady = 2.0
    response_calibration = 2.0 * len(a["calibration_delta_mhz"])
    response_validation = 2.0 * len(a["validation_delta_mhz"])
    kernel_calibration = 2.0
    qualification = response_calibration + response_validation + kernel_calibration
    n = int(c["validation_count"])
    data.update({
        "short_pulse_steady_calls": np.array(steady),
        "short_pulse_response_calibration_calls": np.array(response_calibration),
        "short_pulse_response_validation_calls": np.array(response_validation),
        "short_pulse_kernel_calibration_calls": np.array(kernel_calibration),
        "short_pulse_qualification_calls": np.array(qualification),
        "short_pulse_amortized_calls": np.array(steady + qualification / n),
        "short_pulse_cold_start_calls": np.array(steady + qualification),
        "short_pulse_amortization_count": np.array(n),
        "cost_unit": np.array("counted_solver_calls_mesolve_plus_sesolve"),
        "source_a": np.array(SOURCE_A.name),
        "source_b": np.array(SOURCE_B.name),
        "source_c": np.array(SOURCE_C.name),
    })


def write_cost_csv(data: dict[str, np.ndarray]) -> None:
    rows = ["method,variant,mae_mhz,solver_calls,cost_scope"]
    for idx, label in enumerate(data["ramsey_labels"]):
        rows.append(
            f'Ramsey,"{label}",{float(data["ramsey_mae_mhz"][idx]):.9g},'
            f'{float(data["ramsey_solver_calls_per_measurement"][idx]):g},per_estimate'
        )
    short_mae = float(data["short_pulse_mae_mhz"])
    for variant, key in (
        ("steady", "short_pulse_steady_calls"),
        ("amortized_8", "short_pulse_amortized_calls"),
        ("cold_start", "short_pulse_cold_start_calls"),
    ):
        rows.append(
            f"cubic short pulse,{variant},{short_mae:.9g},"
            f"{float(data[key]):g},counted_generation_chain"
        )
    COST_CSV.write_text("\n".join(rows) + "\n", encoding="utf-8")


def plot(data: dict[str, np.ndarray], stem: str = STEM, panel_label: str = "d") -> None:
    if STYLE.exists():
        plt.style.use(STYLE)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 7.0})
    fig, ax = plt.subplots(figsize=(3.50, 2.70))
    ramsey_x = data["ramsey_solver_calls_per_measurement"]
    ramsey_y = data["ramsey_mae_mhz"]
    ax.plot(ramsey_x, ramsey_y, color="#777777", lw=0.85, zorder=1)
    for idx, label in enumerate(data["ramsey_labels"]):
        filled = bool(data["ramsey_double_sweep"][idx])
        ax.scatter(
            ramsey_x[idx], ramsey_y[idx], s=27, marker="o",
            facecolor="#3572B0" if filled else "white", edgecolor="#3572B0",
            linewidth=0.9, zorder=3,
        )
        if not filled:
            spacing = float(data["ramsey_spacing_ns"][idx])
            offsets_by_spacing = {2.0: (4, 7), 1.0: (4, 4), 0.5: (4, 1)}
            ax.annotate(
                rf"{spacing:g} ns",
                (ramsey_x[idx], ramsey_y[idx]), xytext=offsets_by_spacing[spacing],
                textcoords="offset points", fontsize=5.1, color="#444444",
            )
        else:
            ax.annotate(
                rf"double, $\Delta\tau={float(data['ramsey_spacing_ns'][idx]):g}$ ns",
                (ramsey_x[idx], ramsey_y[idx]), xytext=(-72, 9),
                textcoords="offset points", fontsize=5.1, color="#285D91",
            )

    short_y = float(data["short_pulse_mae_mhz"])
    short_x = np.array([
        float(data["short_pulse_steady_calls"]),
        float(data["short_pulse_amortized_calls"]),
        float(data["short_pulse_cold_start_calls"]),
    ])
    ax.plot(short_x, np.full(3, short_y), color="#2E8B57", lw=0.9, zorder=2)
    ax.scatter(short_x, np.full(3, short_y), s=30, marker="^", facecolor="white",
               edgecolor="#2E8B57", linewidth=0.95, zorder=4)
    ax.annotate(
        "cubic short pulse\n2 steady / 8.5 amortized\n54 cold-start calls",
        (short_x[1], short_y), xytext=(8, -48), textcoords="offset points",
        fontsize=5.2, color="#2E6F48", ha="left",
    )
    ax.text(0.62, 0.70, r"single-sweep Ramsey ($\Delta\tau$)",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=5.1, color="#444444")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Counted solver calls")
    ax.set_ylabel("Held-out local-window MAE (MHz)")
    ax.grid(which="both", color="#DDDDDD", lw=0.45, alpha=0.8)
    ax.text(-0.16, 1.03, panel_label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.5, fontweight="bold")
    ax.text(0.98, 0.96, r"solver calls only; not acquisition time",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.4,
            color="#555555")
    fig.subplots_adjust(left=0.20, right=0.97, bottom=0.20, top=0.92)
    fig.savefig(HERE / f"{stem}.pdf")
    fig.savefig(HERE / f"{stem}.svg")
    fig.savefig(HERE / f"{stem}.png", dpi=600)
    plt.close(fig)


def smoke() -> None:
    a = load_npz(SOURCE_A)
    for truth in (-14.0, 14.0):
        estimate, calls = ramsey_estimate(truth, 2.0, False, a)
        print(f"truth={truth:+.3f} estimate={estimate:+.6f} error={estimate-truth:+.6f} calls={calls}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        smoke()
        return
    data_path = COST_DATA
    if args.plot_only:
        data = load_npz(data_path)
        refresh_short_pulse_costs(data)
        save_npz(data_path, data)
    else:
        data = compute()
        save_npz(data_path, data)
    write_cost_csv(data)
    plot(data)
    plot(data, "figS6-solver-cost-simulation", "S6")
    for idx, label in enumerate(data["ramsey_labels"]):
        print(f"{label}: MAE={data['ramsey_mae_mhz'][idx]:.6f} MHz, "
              f"calls={data['ramsey_solver_calls_per_measurement'][idx]:.0f}")
    print(f"Cubic short pulse: MAE={float(data['short_pulse_mae_mhz']):.6f} MHz, "
          f"calls steady/amortized/cold={float(data['short_pulse_steady_calls']):g}/"
          f"{float(data['short_pulse_amortized_calls']):g}/"
          f"{float(data['short_pulse_cold_start_calls']):g}")


if __name__ == "__main__":
    main()
