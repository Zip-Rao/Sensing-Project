#!/usr/bin/env python3
"""Build numerical placeholder panels for Fig. 1(a) and Fig. 1(c).

The panels intentionally use the currently frozen numerical protocol:
10 ns square pi/2 pulses. They are layout/evidence placeholders until the
experimental data sets and the final physical pulse envelope are frozen.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from qutip import Qobj, QobjEvo, basis, mesolve
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE
DATA = REPO / "paper" / "data" / "fig1-simulation"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "result_sqc"))

from _common import OPTIMAL_FLUX, make_sqc_qubit  # noqa: E402
from sqc.calibration.frequency import CONFIG  # noqa: E402
from sqc.control.flux_signal import FluxSignal  # noqa: E402
from sqc.control.sequence import create_ramsey_pulse  # noqa: E402


MM = 1.0 / 25.4
TWO_PI = 2.0 * np.pi
CHARCOAL = "#24272B"
BLUE = "#2C6FA3"
LIGHT_BLUE = "#DCEAF3"
RED = "#B84A3A"
ORANGE = "#D98424"
GREEN = "#DCE9DF"
MID_GREY = "#747A80"


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Latin Modern Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 7.2,
            "axes.labelsize": 7.2,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.2,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.55,
            "ytick.major.width": 0.55,
            "xtick.major.size": 2.6,
            "ytick.major.size": 2.6,
            "lines.linewidth": 1.0,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )


def save_all(fig: mpl.figure.Figure, stem: str) -> None:
    for suffix, kwargs in (
        ("svg", {}),
        ("pdf", {}),
        ("png", {"dpi": 600}),
    ):
        fig.savefig(OUT / f"{stem}.{suffix}", facecolor="white", **kwargs)
    plt.close(fig)


def analytic_frequency_at_offset(offset: np.ndarray) -> np.ndarray:
    values = []
    for value in np.asarray(offset, dtype=float):
        q = make_sqc_qubit()
        q.change_flux(OPTIMAL_FLUX + float(value))
        values.append(q.frequency / TWO_PI)
    return np.asarray(values)


def build_panel_a() -> None:
    f1 = np.load(
        REPO
        / "result_sqc/frequency_calibration/F1_flux_curve/f_phi_curve_sqc.npz",
        allow_pickle=True,
    )
    loop = np.load(
        REPO
        / "result_sqc/frequency_calibration/F4_hybrid_convergence/freq_calibration_sqc.npz",
        allow_pickle=True,
    )

    flux_dense = np.linspace(float(f1["flux"].min()), float(f1["flux"].max()), 401)
    f_dense = analytic_frequency_at_offset(flux_dense - OPTIMAL_FLUX)
    track_flux = OPTIMAL_FLUX + np.asarray(loop["trk_V"], dtype=float)
    track_freq = analytic_frequency_at_offset(track_flux - OPTIMAL_FLUX)
    f_ref = float(analytic_frequency_at_offset(np.array([0.0]))[0])
    f_target = f_ref + float(loop["f_target_offset_mhz"][0]) * 1e-3

    order = np.argsort(f_dense)
    phi_target = float(np.interp(f_target, f_dense[order], flux_dense[order]))

    fig, ax = plt.subplots(figsize=(88 * MM, 44 * MM))
    fig.subplots_adjust(left=0.155, right=0.975, bottom=0.22, top=0.92)

    ax.plot(flux_dense, f_dense, color=CHARCOAL, lw=1.25, zorder=1)
    ax.text(
        0.22,
        0.87,
        "withheld dispersion",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=MID_GREY,
        fontsize=6.2,
    )

    ax.axhline(f_target, color=RED, lw=0.75, ls=(0, (3, 2)), zorder=0)
    ax.plot(phi_target, f_target, marker="*", ms=7.2, color=RED, zorder=5)
    ax.annotate(
        "target",
        xy=(phi_target, f_target),
        xytext=(5, -13),
        textcoords="offset points",
        color=RED,
        ha="left",
        va="top",
        fontsize=6.4,
    )

    ax.plot(track_flux, track_freq, color=BLUE, lw=1.05, zorder=3)
    ax.scatter(
        track_flux,
        track_freq,
        s=18,
        facecolor="white",
        edgecolor=BLUE,
        linewidth=0.9,
        zorder=4,
    )
    ax.scatter(track_flux[0], track_freq[0], s=25, color=ORANGE, zorder=5)
    ax.annotate(
        "seed",
        xy=(track_flux[0], track_freq[0]),
        xytext=(5, 7),
        textcoords="offset points",
        color=ORANGE,
        fontsize=6.4,
    )
    ax.annotate(
        "",
        xy=(track_flux[2], track_freq[2]),
        xytext=(track_flux[1], track_freq[1]),
        arrowprops={"arrowstyle": "-|>", "color": BLUE, "lw": 0.8, "mutation_scale": 7},
        zorder=5,
    )
    ax.text(
        track_flux[1] - 0.0006,
        track_freq[1] + 0.0015,
        "Track updates",
        color=BLUE,
        ha="right",
        va="bottom",
        fontsize=6.2,
    )

    ax.set_xlabel(r"Flux bias $\Phi/\Phi_0$")
    ax.set_ylabel(r"Qubit frequency $f_{01}$ (GHz)")
    ax.set_xlim(0.932, 0.9615)
    ax.set_ylim(f_target - 0.0045, f_ref + 0.007)
    ax.set_xticks([0.935, 0.945, 0.955])
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.13, 1.02, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=8.5)

    save_all(fig, "fig1a-numerical-placeholder-v1")


def build_panel_a_global_inset() -> None:
    loop = np.load(
        REPO
        / "result_sqc/frequency_calibration/F4_hybrid_convergence/freq_calibration_sqc.npz",
        allow_pickle=True,
    )
    flux_global = np.linspace(0.82, 1.18, 501)
    f_global = analytic_frequency_at_offset(flux_global - OPTIMAL_FLUX)
    track_flux = OPTIMAL_FLUX + np.asarray(loop["trk_V"], dtype=float)
    track_freq = analytic_frequency_at_offset(track_flux - OPTIMAL_FLUX)
    f_ref = float(analytic_frequency_at_offset(np.array([0.0]))[0])
    f_target = f_ref + float(loop["f_target_offset_mhz"][0]) * 1e-3
    order = np.argsort(f_global)
    phi_target = float(np.interp(f_target, f_global[order], flux_global[order]))

    fig, ax = plt.subplots(figsize=(88 * MM, 56 * MM))
    fig.subplots_adjust(left=0.16, right=0.975, bottom=0.18, top=0.94)

    ax.plot(flux_global, f_global, color=CHARCOAL, lw=1.25, zorder=2)
    ax.text(
        0.07,
        0.94,
        "withheld dispersion",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=MID_GREY,
        fontsize=6.2,
    )
    ax.set_xlabel(r"Flux bias $\Phi/\Phi_0$")
    ax.set_ylabel(r"Qubit frequency $f_{01}$ (GHz)")
    ax.set_xlim(0.82, 1.18)
    ax.set_ylim(3.455, 3.815)
    ax.set_xticks([0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15])
    ax.set_yticks([3.5, 3.6, 3.7, 3.8])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.15, 1.015, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=8.5)

    axins = ax.inset_axes([0.39, 0.13, 0.44, 0.38])
    axins.set_facecolor("white")
    local_flux = np.linspace(0.9335, 0.9595, 151)
    local_freq = analytic_frequency_at_offset(local_flux - OPTIMAL_FLUX)
    axins.plot(local_flux, local_freq, color=CHARCOAL, lw=0.85, zorder=1)
    axins.axhline(f_target, color=RED, lw=0.6, ls=(0, (3, 2)), zorder=0)
    axins.plot(track_flux, track_freq, color=BLUE, lw=0.85, zorder=3)
    axins.scatter(
        track_flux,
        track_freq,
        s=10,
        facecolor="white",
        edgecolor=BLUE,
        linewidth=0.7,
        zorder=4,
    )
    axins.scatter(track_flux[0], track_freq[0], s=15, color=ORANGE, zorder=5)
    axins.plot(phi_target, f_target, marker="*", ms=5.2, color=RED, zorder=5)
    axins.annotate(
        "seed",
        xy=(track_flux[0], track_freq[0]),
        xytext=(-3, -4),
        textcoords="offset points",
        ha="right",
        va="top",
        fontsize=5.1,
        color=ORANGE,
    )
    axins.text(0.965, 0.055, r"$f_{\mathrm{tar}}$", transform=axins.transAxes,
               ha="right", va="bottom", fontsize=5.1, color=RED)
    axins.annotate(
        "",
        xy=(track_flux[2], track_freq[2]),
        xytext=(track_flux[1], track_freq[1]),
        arrowprops={"arrowstyle": "-|>", "color": BLUE, "lw": 0.65, "mutation_scale": 5.5},
        zorder=5,
    )
    axins.set_xlim(0.9335, 0.9595)
    axins.set_ylim(f_target - 0.0025, f_ref + 0.003)
    axins.set_xticks([0.94, 0.95])
    axins.set_yticks([3.76, 3.78])
    axins.tick_params(labelsize=5.0, length=1.8, width=0.45, pad=1.3)
    for spine in axins.spines.values():
        spine.set_linewidth(0.55)
        spine.set_color(MID_GREY)
    axins.text(0.04, 0.93, "local calibration", transform=axins.transAxes,
               ha="left", va="top", fontsize=5.2, color=CHARCOAL)

    patch, connector1, connector2 = mark_inset(
        ax,
        axins,
        loc1=2,
        loc2=1,
        fc="#F5D6A0",
        ec=ORANGE,
        alpha=0.28,
        lw=0.7,
    )
    patch.set_zorder(1)
    connector1.set_color(MID_GREY)
    connector2.set_color(MID_GREY)
    connector1.set_linewidth(0.55)
    connector2.set_linewidth(0.55)

    save_all(fig, "fig1a-global-inset-placeholder-v2")


def simulate_pdiff(detuning_mhz: np.ndarray) -> np.ndarray:
    qubit = make_sqc_qubit()
    omega_d = float(qubit.frequency)
    n_levels = qubit.n_levels
    psi_e = basis(n_levels, 1)
    t_rabi = CONFIG.pulse.t_rabi.copy()
    t_global = CONFIG.pulse.t_global.copy()
    t_sig = CONFIG.pulse.make_time(0, 300)
    sigma_z = Qobj(np.diag([1.0, -1.0] + [0.0] * (n_levels - 2)))

    qubit.qubit_in_mag(FluxSignal(type=0, t_list=t_sig), frame=1, omega_d=omega_d)
    ctrl_x = create_ramsey_pulse(
        t_rabi,
        0.0,
        omega_d=omega_d,
        phase1=np.pi / 2,
        phase2=0.0,
        qubit=qubit,
        rotation_angle=np.pi / 2,
        envelope="square",
    )
    ctrl_mx = create_ramsey_pulse(
        t_rabi,
        0.0,
        omega_d=omega_d,
        phase1=np.pi / 2,
        phase2=np.pi,
        qubit=qubit,
        rotation_angle=np.pi / 2,
        envelope="square",
    )
    h_base = QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1)
    h_x = QobjEvo(ctrl_x.hamiltonian_on(t_global), tlist=t_global, order=1)
    h_mx = QobjEvo(ctrl_mx.hamiltonian_on(t_global), tlist=t_global, order=1)
    e_op = psi_e * psi_e.dag()

    values = []
    for delta_mhz in np.asarray(detuning_mhz, dtype=float):
        delta_rad_ghz = TWO_PI * delta_mhz * 1e-3
        coeff = -0.5 * delta_rad_ghz * np.ones_like(t_global)
        h_detuning = QobjEvo([[sigma_z, coeff]], tlist=t_global, order=1)
        res_x = mesolve(
            h_base + h_detuning + h_x,
            qubit.state,
            t_global,
            [],
            e_ops=[e_op],
            options={"max_step": float(CONFIG.awg.dt)},
        )
        res_mx = mesolve(
            h_base + h_detuning + h_mx,
            qubit.state,
            t_global,
            [],
            e_ops=[e_op],
            options={"max_step": float(CONFIG.awg.dt)},
        )
        values.append((float(res_x.expect[0][-1]) - float(res_mx.expect[0][-1])) / 2.0)
    return np.asarray(values)


def build_panel_c() -> None:
    dense_delta = np.linspace(-28.0, 28.0, 113)
    validation_delta = np.array([-20.5, -15.5, -10.5, -5.5, 0.0, 5.5, 10.5, 15.5, 20.5])
    dense_p = simulate_pdiff(dense_delta)
    validation_p = simulate_pdiff(validation_delta)

    calibration_mask = np.abs(dense_delta) <= 10.0
    coeff = np.polyfit(dense_delta[calibration_mask], dense_p[calibration_mask], 3)
    slope = float(np.polyval(np.polyder(coeff), 0.0))
    tangent_delta = np.array([-7.5, 7.5])
    tangent_p = float(np.polyval(coeff, 0.0)) + slope * tangent_delta

    np.savez(
        DATA / "fig1c-v1-data.npz",
        dense_delta_mhz=dense_delta,
        dense_pdiff=dense_p,
        validation_delta_mhz=validation_delta,
        validation_pdiff=validation_p,
        pulse_window_ns=np.array([10.0]),
        rotation_angle_rad=np.array([np.pi / 2]),
        envelope=np.array(["square"]),
    )

    fig, ax = plt.subplots(figsize=(88 * MM, 44 * MM))
    fig.subplots_adjust(left=0.155, right=0.975, bottom=0.22, top=0.92)

    ax.axvspan(-12.0, 12.0, color=GREEN, alpha=0.78, lw=0, zorder=0)
    ax.text(
        0.5,
        0.955,
        "candidate validated interval",
        transform=ax.transAxes,
        ha="center",
        va="top",
        color="#526C59",
        fontsize=6.1,
    )
    ax.axhline(0.0, color="#A8ADB2", lw=0.55, zorder=0)
    ax.axvline(0.0, color="#A8ADB2", lw=0.55, ls=(0, (2, 2)), zorder=0)
    ax.plot(dense_delta, dense_p, color=BLUE, lw=1.35, zorder=2)
    ax.scatter(
        validation_delta,
        validation_p,
        s=15,
        facecolor="white",
        edgecolor=ORANGE,
        linewidth=0.95,
        zorder=3,
    )
    ax.plot(tangent_delta, tangent_p, color=CHARCOAL, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_xlabel(r"Detuning $\Delta/(2\pi)$ (MHz)")
    ax.set_ylabel(r"Differential response $p_{\mathrm{d}}$")
    ax.set_xlim(-28, 28)
    pad = 0.07 * float(np.ptp(dense_p))
    ax.set_ylim(float(dense_p.min() - pad), float(dense_p.max() + pad))
    ax.set_xticks([-20, -10, 0, 10, 20])
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.13, 1.02, "(c)", transform=ax.transAxes, fontweight="bold", fontsize=8.5)

    save_all(fig, "fig1c-numerical-placeholder-v1")


def main() -> None:
    configure_matplotlib()
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    build_panel_a()
    build_panel_a_global_inset()
    build_panel_c()


if __name__ == "__main__":
    main()
