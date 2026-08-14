#!/usr/bin/env python3
"""Generate final Fig. 2(a): short-pulse response and local models.

The default protocol follows the manuscript parameters and the sin-squared
envelope shown in Fig. 1(b).  ``make_envelope`` is the single protocol hook:
it also accepts square, Gaussian, or a custom array without changing the
calibration/validation and plotting logic.
"""
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

from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.config import CONFIG
from sqc.devices.transmon import TransmonQubit
from fig2_simulation_common import (
    CANDIDATE_INTERVAL_MHZ,
    MAX_FREQUENCY_ERROR_MHZ,
    MAX_MODEL_ERROR,
    RESPONSE_DATA,
    load_npz,
    save_npz,
)


TWO_PI = 2.0 * np.pi
DT_NS = float(CONFIG.awg.dt)
PULSE_WINDOW_NS = 10.0
EC = TWO_PI * 0.2
EJ = TWO_PI * 10.0
T1_NS = 100_000.0
T2_NS = 50_000.0
N_LEVELS = 2
OMEGA_D = TWO_PI * 3.7802780210349702
ENVELOPE_DEFAULT = "sin2"
CAL_DELTA_MHZ = np.arange(-24.0, 24.01, 4.0)
VAL_DELTA_MHZ = np.arange(-22.0, 22.01, 4.0)
DENSE_DELTA_MHZ = np.arange(-28.0, 28.01, 0.5)

HERE = Path(__file__).resolve().parent
STEM = "fig2a-short-pulse-response"
STYLE = ROOT.parent / "scholaraio" / ".claude" / "skills" / "draw" / "nature_pub.mplstyle"


def make_envelope(kind: str | np.ndarray, t_rabi: np.ndarray) -> np.ndarray:
    """Return a dimensionless pulse profile on ``t_rabi``."""
    if not isinstance(kind, str):
        profile = np.asarray(kind, dtype=float)
        if profile.shape != t_rabi.shape:
            raise ValueError("custom envelope must match t_rabi")
        return profile
    if kind == "sin2":
        u = (t_rabi - t_rabi[0]) / PULSE_WINDOW_NS
        return np.sin(np.pi * u) ** 2
    if kind == "square":
        return np.ones_like(t_rabi)
    if kind == "gaussian":
        center = 0.5 * PULSE_WINDOW_NS
        return np.exp(-0.5 * ((t_rabi - center) / (PULSE_WINDOW_NS / 4.0)) ** 2)
    raise ValueError(f"unsupported envelope: {kind}")


def make_qubit() -> TransmonQubit:
    return TransmonQubit(
        EC=EC, EJ=EJ, T1=T1_NS, T2=T2_NS, flux=0.0,
        state=0, n_levels=N_LEVELS,
    )


def branch_population(delta_mhz: float, phase2: float, envelope: np.ndarray) -> float:
    """Run one dissipative master-equation branch at fixed detuning."""
    qubit = make_qubit()
    t_rabi = np.arange(0.0, PULSE_WINDOW_NS + 0.5 * DT_NS, DT_NS)
    t_global = np.arange(0.0, 2.0 * PULSE_WINDOW_NS + DT_NS, DT_NS)
    omega_q = OMEGA_D + TWO_PI * float(delta_mhz) * 1e-3
    flux = FluxSignal(type=0, t_list=t_global)
    qubit.qubit_in_mag(flux, frame=1, omega_d=OMEGA_D)
    # Set the constant rotating-frame detuning explicitly while retaining the
    # project qubit, pulse, and collapse-operator implementations.
    detuning_shift = omega_q - qubit.frequency
    h_base = QobjEvo(qubit.H_list, tlist=t_global, order=1)
    # Delta/2*sigma_z equals Delta*n up to an irrelevant scalar offset in
    # the {|0>, |1>} basis used by the project Hamiltonian.
    h_base += detuning_shift * qubit.n
    ctrl = create_ramsey_pulse(
        t_rabi, tau=0.0, omega_d=OMEGA_D,
        phase1=np.pi / 2.0, phase2=phase2, qubit=qubit,
        rotation_angle=np.pi / 2.0, envelope=envelope,
    )
    h_total = h_base + QobjEvo(ctrl.hamiltonian_on(t_global), tlist=t_global, order=1)
    projector_e = basis(N_LEVELS, 1) * basis(N_LEVELS, 1).dag()
    result = mesolve(
        h_total, qubit.state, t_global, qubit.get_collapse_operators(),
        e_ops=[projector_e],
        options={"atol": 1e-10, "rtol": 1e-8, "max_step": DT_NS / 2.0},
    )
    return float(result.expect[0][-1])


def response(delta_mhz: np.ndarray, envelope: np.ndarray) -> np.ndarray:
    out = []
    for delta in np.asarray(delta_mhz, dtype=float):
        p_plus = branch_population(delta, 0.0, envelope)
        p_minus = branch_population(delta, np.pi, envelope)
        out.append(0.5 * (p_plus - p_minus))
    return np.asarray(out)


def compute(envelope_name: str) -> dict[str, np.ndarray]:
    t_rabi = np.arange(0.0, PULSE_WINDOW_NS + 0.5 * DT_NS, DT_NS)
    envelope = make_envelope(envelope_name, t_rabi)
    cal_raw = response(CAL_DELTA_MHZ, envelope)
    val_raw = response(VAL_DELTA_MHZ, envelope)
    dense_raw = response(DENSE_DELTA_MHZ, envelope)
    b0 = float(cal_raw[np.argmin(np.abs(CAL_DELTA_MHZ))])
    cal = cal_raw - b0
    val = val_raw - b0
    dense = dense_raw - b0
    x_cal = TWO_PI * CAL_DELTA_MHZ * 1e-3
    design = np.column_stack((x_cal, x_cal**3))
    c1, c3 = np.linalg.lstsq(design, cal, rcond=None)[0]
    g1, g3 = float(c1), float(6.0 * c3)
    fold_rad = np.sqrt(-2.0 * g1 / g3) if g1 * g3 < 0.0 else np.nan
    fold_mhz = float(fold_rad / TWO_PI * 1e3) if np.isfinite(fold_rad) else np.nan
    x_val = TWO_PI * VAL_DELTA_MHZ * 1e-3
    val_model = g1 * x_val + (g3 / 6.0) * x_val**3
    delta_val = CANDIDATE_INTERVAL_MHZ
    candidate = np.abs(VAL_DELTA_MHZ) <= delta_val + 1e-12
    if np.max(np.abs(val[candidate] - val_model[candidate])) > MAX_MODEL_ERROR:
        raise ValueError("frozen candidate interval fails its response-model qualification")
    if not delta_val < fold_mhz:
        raise ValueError("frozen candidate interval crosses the fitted cubic fold")
    x_dense = TWO_PI * DENSE_DELTA_MHZ * 1e-3
    return {
        "calibration_delta_mhz": CAL_DELTA_MHZ,
        "calibration_pdiff_corrected": cal,
        "validation_delta_mhz": VAL_DELTA_MHZ,
        "validation_pdiff_corrected": val,
        "dense_delta_mhz": DENSE_DELTA_MHZ,
        "dense_pdiff_corrected": dense,
        "linear_model": g1 * x_dense,
        "cubic_model": g1 * x_dense + (g3 / 6.0) * x_dense**3,
        "t_rabi_ns": t_rabi,
        "envelope_samples": envelope,
        "envelope_name": np.array(envelope_name),
        "pulse_window_ns": np.array(PULSE_WINDOW_NS),
        "dt_ns": np.array(DT_NS),
        "EC_over_h_GHz": np.array(EC / TWO_PI),
        "EJ_over_h_GHz": np.array(EJ / TWO_PI),
        "T1_ns": np.array(T1_NS),
        "T2_ns": np.array(T2_NS),
        "rotation_angle_rad": np.array(np.pi / 2.0),
        "omega_d_rad_GHz": np.array(OMEGA_D),
        "b0": np.array(b0),
        "G1_ns": np.array(g1),
        "G3_ns3": np.array(g3),
        "delta_fold_mhz": np.array(fold_mhz),
        "delta_val_mhz": np.array(delta_val),
        "validation_max_model_error_threshold": np.array(MAX_MODEL_ERROR),
        "validation_max_frequency_error_mhz": np.array(MAX_FREQUENCY_ERROR_MHZ),
        "candidate_interval_role": np.array("simulation_defined_candidate_interval"),
        "sigma": np.array(0.0),
        "uncertainty_source": np.array("deterministic_zero"),
        "calibration_role": np.array("response_calibration"),
        "validation_role": np.array("simulated_held_out_response_validation"),
    }


def plot(data: dict[str, np.ndarray], stem: str = STEM, panel_label: str = "a") -> None:
    if STYLE.exists():
        plt.style.use(STYLE)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 7.0})
    fig, ax = plt.subplots(figsize=(3.50, 2.70))
    delta = data["dense_delta_mhz"]
    fold = float(data["delta_fold_mhz"])
    valid = float(data["delta_val_mhz"])
    branch = np.ones_like(delta, dtype=bool)
    if np.isfinite(fold):
        branch = np.abs(delta) <= fold
    ax.axvspan(-valid, valid, color="#D9EAD3", alpha=0.75, lw=0, zorder=0)
    if np.isfinite(fold):
        for sign in (-1.0, 1.0):
            ax.axvline(sign * fold, color="#777777", ls="--", lw=0.8, zorder=1)
    ax.plot(delta, data["dense_pdiff_corrected"], color="#222222", lw=1.25,
            label="Master equation", zorder=2)
    ax.plot(delta[branch], data["linear_model"][branch], color="#3572B0", lw=1.0,
            label="Linear response", zorder=3)
    ax.plot(delta[branch], data["cubic_model"][branch], color="#C43C39", lw=1.15,
            label="Cubic response", zorder=4)
    ax.scatter(data["calibration_delta_mhz"], data["calibration_pdiff_corrected"],
               s=17, marker="o", color="#3572B0", edgecolor="white", linewidth=0.45,
               label="Response calibration", zorder=5)
    ax.scatter(data["validation_delta_mhz"], data["validation_pdiff_corrected"],
               s=21, marker="D", facecolor="white", edgecolor="#C43C39", linewidth=0.8,
               label="Held-out validation", zorder=6)
    ax.axhline(0.0, color="#999999", lw=0.55, zorder=0)
    ax.set_xlabel(r"True detuning $\Delta/(2\pi)$ (MHz)")
    ax.set_ylabel(r"Corrected differential population $\widetilde{p}_{\rm d}$")
    ax.set_xlim(-28.0, 28.0)
    ax.text(-0.16, 1.03, panel_label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.5, fontweight="bold")
    ax.text(0.50, 0.04, rf"$|\Delta|/(2\pi)\leq {valid:.0f}$ MHz",
            transform=ax.transAxes, ha="center", va="bottom", color="#3E6B3A", fontsize=6.2)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
              frameon=False, handlelength=1.6, columnspacing=0.9, fontsize=5.8)
    fig.subplots_adjust(left=0.18, right=0.97, bottom=0.19, top=0.76)
    fig.savefig(HERE / f"{stem}.pdf")
    fig.savefig(HERE / f"{stem}.svg")
    fig.savefig(HERE / f"{stem}.png", dpi=600)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--envelope", choices=("sin2", "square", "gaussian"),
                        default=ENVELOPE_DEFAULT)
    parser.add_argument("--plot-only", action="store_true",
                        help="reuse the existing NPZ and only rebuild exports")
    args = parser.parse_args()
    data_path = RESPONSE_DATA
    if args.plot_only:
        data = load_npz(data_path)
        data["candidate_interval_role"] = np.array("simulation_defined_candidate_interval")
        save_npz(data_path, data)
    else:
        data = compute(args.envelope)
        save_npz(data_path, data)
    plot(data)
    plot(data, "figS2-response-simulation", "S2")
    print(f"envelope={args.envelope}")
    print(f"G1={float(data['G1_ns']):.8g} ns")
    print(f"G3={float(data['G3_ns3']):.8g} ns^3")
    print(f"Delta_fold/(2pi)={float(data['delta_fold_mhz']):.4f} MHz")
    print(f"Delta_val/(2pi)={float(data['delta_val_mhz']):.4f} MHz")
    print(f"output={HERE / STEM}")


if __name__ == "__main__":
    main()
