#!/usr/bin/env python3
"""Generate compact report figures for transient characterization and CZ validation.

The script uses only simulated data produced by the project models.  E1 runs a
finite-duration probe through the control-line model and reconstructs the
on-chip response with the transient protocol.  E2 reuses the frozen P7
Cryoscope calibration and evaluates the same two-qubit model before and after
predistortion on an amplitude-duration grid.
"""
from __future__ import annotations

import contextlib
import io
import json
import os
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "result_sqc"))

E1_CACHE = HERE / "fig_E1_transient_characterization_data.npz"
E2_CACHE = HERE / "fig_E2_predistortion_cz_validation_data.npz"
P7_CACHE = REPO / "result_sqc/predistortion/P7_cz_gate_impact/cz_gate_impact_sqc.npz"

INK = "#202124"
GRAY = "#8b9097"
RED = "#c23b33"
BLUE = "#276b9a"
TEAL = "#23877b"
LIGHT = "#eceff1"


def _quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def _save_npz(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **data)


def _make_probe(t):
    """Finite-band probe with both slow and fast time-domain structure."""
    g1 = np.exp(-0.5 * ((t - 55.0) / 9.0) ** 2)
    g2 = np.exp(-0.5 * ((t - 103.0) / 15.0) ** 2)
    g3 = np.exp(-0.5 * ((t - 148.0) / 7.0) ** 2)
    return 0.0045 * (g1 - 0.72 * g2 + 0.50 * g3)


def compute_e1(recompute=False):
    if E1_CACHE.exists() and not recompute:
        z = np.load(E1_CACHE)
        return {k: z[k] for k in z.files}

    from sqc.config import CONFIG
    from sqc.control.flux_signal import FluxSignal
    from sqc.control.waveform import Waveform
    from sqc.devices.transmon import TransmonQubit
    from sqc.experiments.transient import TransientSensingExperiment
    from sqc.hardware.control_line import ControlLine
    from sqc.hardware.distortion import MultiExponentialDistortion
    from sqc.reconstruction.transient import TransientReconstruction

    import _common as C

    t = CONFIG.pulse.make_time(0, 200)
    x = _make_probe(t)
    line = ControlLine(
        name="Z0", kind="z", source="AWG0", target="Q0",
        transfer_function=MultiExponentialDistortion(
            amplitudes=[0.04, 0.02], taus=[80.0, 400.0]
        ),
    )
    y = line.apply(Waveform(t_list=t, samples=x)).samples
    signal = FluxSignal(type=8, t_list=t, signal=y)
    qubit = TransmonQubit(
        EC=C.EC, EJ=C.EJ, T1=C.T1, T2=C.T2,
        flux=C.OPTIMAL_FLUX, n_levels=3,
    )
    result = _quiet(TransientSensingExperiment(qubit=qubit, flux_signal=signal).run)
    recon = TransientReconstruction(method="wiener", lambda_reg=C.WIENER_LAMBDA)
    recovered_raw = recon.reconstruct(result, kernel=result.data["kernel"]).signal
    # The first-order transient kernel fixes waveform shape but leaves one
    # scalar measurement gain.  Estimate that gain from this known test probe,
    # exactly as an experimental amplitude calibration would do.
    transient_gain = float(
        np.dot(recovered_raw, y) / np.dot(recovered_raw, recovered_raw)
    )
    recovered = transient_gain * recovered_raw

    # Cross-spectral estimate with light regularization. Only supported bins
    # are retained for plotting; unsupported bins are never interpreted.
    dt = float(t[1] - t[0])
    window = np.hanning(len(t))
    x0 = (x - np.mean(x)) * window
    y0 = (y - np.mean(y)) * window
    yr0 = (recovered - np.mean(recovered)) * window
    xf = np.fft.rfft(x0)
    yf = np.fft.rfft(y0)
    yrf = np.fft.rfft(yr0)
    freq = np.fft.rfftfreq(len(t), d=dt)
    reg = 1e-4 * np.max(np.abs(xf) ** 2)
    h_true = yf * np.conj(xf) / (np.abs(xf) ** 2 + reg)
    h_rec = yrf * np.conj(xf) / (np.abs(xf) ** 2 + reg)
    support = (np.abs(xf) >= 0.08 * np.max(np.abs(xf))) & (freq > 0)

    data = {
        "t": t, "probe": x, "chip_true": y, "chip_reconstructed": recovered,
        "scan": result.axes["scan"], "delta_p": result.data["delta_p"],
        "kernel": result.data["kernel"], "freq": freq,
        "h_true": h_true, "h_reconstructed": h_rec, "support": support,
        "transient_gain": np.array([transient_gain]),
        "waveform_rmse": np.array([np.sqrt(np.mean((recovered - y) ** 2))]),
    }
    _save_npz(E1_CACHE, data)
    return data


def _load_p7_helpers():
    p7_dir = REPO / "result_sqc/predistortion/P7_cz_gate_impact"
    sys.path.insert(0, str(p7_dir))
    import generate_cz_gate_impact_sqc as p7
    return p7


def _cz_waveform(duration, amplitude, rise_frac, dt):
    t = np.arange(0.0, duration + 0.5 * dt, dt)
    edge = rise_frac * duration
    phi = np.full_like(t, amplitude)
    left = t < edge
    right = t > duration - edge
    phi[left] = amplitude * np.sin(0.5 * np.pi * t[left] / edge) ** 2
    phi[right] = amplitude * np.sin(
        0.5 * np.pi * (duration - t[right]) / edge
    ) ** 2
    return t, phi


def compute_e2(recompute=False):
    if E2_CACHE.exists() and not recompute:
        z = np.load(E2_CACHE)
        return {k: z[k] for k in z.files}
    if not P7_CACHE.exists():
        raise FileNotFoundError(f"Run P7 first: {P7_CACHE}")

    from sqc.calibration.waveform import PredistortionDesigner
    from sqc.hardware.distortion import SingleExponentialDistortion

    p7 = _load_p7_helpers()
    base = np.load(P7_CACHE)
    dt = float(base["t_cz"][1] - base["t_cz"][0])
    fitted = SingleExponentialDistortion(
        amplitude=float(base["fitted_amplitude"][0]),
        tau=float(base["fitted_tau"][0]),
    )
    inverse = PredistortionDesigner(method="auto", regularization=1e-4).design(
        fitted, dt=dt
    )
    line = p7._make_control_line()

    amplitudes = np.linspace(0.148, 0.188, 17)
    durations = np.linspace(126.0, 168.0, 17)
    pop_unc = np.zeros((len(amplitudes), len(durations)))
    pop_cor = np.zeros_like(pop_unc)
    infid_unc = np.zeros_like(pop_unc)
    infid_cor = np.zeros_like(pop_unc)

    for ia, amp in enumerate(amplitudes):
        for it, duration in enumerate(durations):
            t, target = _cz_waveform(duration, amp, p7.CZ_RISE_FRAC, dt)
            unc = p7._apply_control_line(t, target, line)
            awg = inverse.apply_to_waveform(p7._flux_to_waveform(t, target)).samples
            cor = p7._apply_control_line(t, awg, line)
            ru = p7._simulate_cz(t, unc)[0]
            rc = p7._simulate_cz(t, cor)[0]
            pop_unc[ia, it] = ru.populations["pop_20"]
            pop_cor[ia, it] = rc.populations["pop_20"]
            infid_unc[ia, it] = ru.infidelity
            infid_cor[ia, it] = rc.infidelity

    data = {
        "t_cz": base["t_cz"], "phi_target": base["phi_target"],
        "phi_chip_uncorrected": base["phi_chip_uncorrected"],
        "phi_chip_corrected": base["phi_chip_corrected"],
        "awg_predistorted": base["awg_predistorted"],
        "amplitudes": amplitudes, "durations": durations,
        "pop20_uncorrected": pop_unc, "pop20_corrected": pop_cor,
        "infidelity_uncorrected": infid_unc,
        "infidelity_corrected": infid_cor,
        "nominal_amplitude": base["cz_flux_amp"],
        "nominal_duration": base["cz_duration"],
        "nominal_infid_uncorrected": base["infidelity_uncorrected"],
        "nominal_infid_corrected": base["infidelity_corrected"],
    }
    _save_npz(E2_CACHE, data)
    return data


def _style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 7.0,
        "axes.labelsize": 7.0,
        "axes.titlesize": 7.5,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "legend.fontsize": 6.2,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "lines.linewidth": 1.15,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })


def _panel(ax, label):
    ax.text(-0.16, 1.06, label, transform=ax.transAxes, ha="left",
            va="bottom", fontsize=8, fontweight="bold")


def _save_figure(fig, stem):
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"{stem}.{suffix}", dpi=400 if suffix == "png" else None,
                    bbox_inches="tight", facecolor="white")


def plot_e1(d):
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, Rectangle

    fig = plt.figure(figsize=(7.2, 4.15), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[0.72, 1.36, 1.14],
                          height_ratios=[0.88, 1.12])
    ax_flow = fig.add_subplot(gs[:, 0])
    ax_raw = fig.add_subplot(gs[0, 1])
    ax_wave = fig.add_subplot(gs[1, 1])
    ax_mag = fig.add_subplot(gs[0, 2])
    ax_phase = fig.add_subplot(gs[1, 2], sharex=ax_mag)

    ax_flow.set_xlim(0, 1); ax_flow.set_ylim(0, 1); ax_flow.axis("off")
    labels = ["Test probe", "Control line", "Qubit sensing", "Reconstruction"]
    ys = [0.83, 0.61, 0.39, 0.17]
    for label, y in zip(labels, ys):
        rect = Rectangle((0.10, y - 0.065), 0.80, 0.13,
                         facecolor="white", edgecolor=INK, linewidth=0.8)
        ax_flow.add_patch(rect)
        ax_flow.text(0.50, y, label, ha="center", va="center", fontsize=6.8)
    for y0, y1 in zip(ys[:-1], ys[1:]):
        ax_flow.add_patch(FancyArrowPatch((0.50, y0 - 0.07), (0.50, y1 + 0.07),
                                         arrowstyle="-|>", mutation_scale=7,
                                         linewidth=0.75, color=INK))
    _panel(ax_flow, "a")

    ax_raw.plot(d["scan"], d["delta_p"], color=TEAL, lw=1.0)
    ax_raw.axhline(0, color=GRAY, lw=0.5)
    ax_raw.set(xlabel="Delay (ns)", ylabel=r"$\Delta P_e$", title="Transient observable")
    ax_raw.grid(alpha=0.18, lw=0.45)
    inset = ax_raw.inset_axes([0.66, 0.58, 0.30, 0.34])
    inset.plot(np.arange(len(d["kernel"])) * 0.5, d["kernel"], color=INK, lw=0.8)
    inset.set_title("Kernel", fontsize=5.8, pad=1)
    inset.tick_params(labelsize=5.0, length=2)
    _panel(ax_raw, "b")

    ax_wave.plot(d["t"], 1e3 * d["probe"], color=GRAY, ls=":", label="AWG input")
    ax_wave.plot(d["t"], 1e3 * d["chip_true"], color=INK, label="On-chip truth")
    ax_wave.plot(d["t"], 1e3 * d["chip_reconstructed"], color=BLUE, ls="--",
                 label="Transient reconstruction")
    ax_wave.set(xlabel="Time (ns)", ylabel=r"Flux ($10^{-3}\Phi_0$)",
                title="On-chip waveform reconstruction")
    ax_wave.legend(frameon=False, ncol=1, loc="upper right")
    ax_wave.grid(alpha=0.18, lw=0.45)
    _panel(ax_wave, "c")

    mask = d["support"].astype(bool)
    f_mhz = 1e3 * d["freq"][mask]
    qualified = f_mhz <= 15.0
    ax_mag.plot(f_mhz, np.abs(d["h_true"][mask]), color=INK, label="Model truth")
    ax_mag.plot(f_mhz[qualified], np.abs(d["h_reconstructed"][mask])[qualified],
                color=BLUE, ls="--", label="Transient estimate")
    ax_mag.plot(f_mhz[~qualified], np.abs(d["h_reconstructed"][mask])[~qualified],
                color=GRAY, ls=":", label="Outside qualified band")
    ax_mag.axvline(15.0, color=GRAY, lw=0.65)
    ax_mag.set(ylabel=r"$|H(f)|$", title="Transfer-function identification")
    ax_mag.legend(frameon=False, loc="best")
    ax_mag.grid(alpha=0.18, lw=0.45)
    _panel(ax_mag, "d")

    phase_true = np.unwrap(np.angle(d["h_true"][mask]))
    phase_rec = np.unwrap(np.angle(d["h_reconstructed"][mask]))
    ax_phase.plot(f_mhz, phase_true, color=INK)
    ax_phase.plot(f_mhz[qualified], phase_rec[qualified], color=BLUE, ls="--")
    ax_phase.plot(f_mhz[~qualified], phase_rec[~qualified], color=GRAY, ls=":")
    ax_phase.axvline(15.0, color=GRAY, lw=0.65)
    ax_phase.set(xlabel="Frequency (MHz)", ylabel="Phase (rad)")
    ax_phase.grid(alpha=0.18, lw=0.45)

    _save_figure(fig, "fig_E1_transient_characterization")
    plt.close(fig)


def plot_e2(d):
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(7.2, 4.35), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.82, 1.18])
    ax_wave = fig.add_subplot(gs[0, :])
    ax_unc = fig.add_subplot(gs[1, 0])
    ax_cor = fig.add_subplot(gs[1, 1], sharex=ax_unc, sharey=ax_unc)

    t = d["t_cz"]
    target = d["phi_target"]
    unc = d["phi_chip_uncorrected"]
    cor = d["phi_chip_corrected"]
    tol = 0.0025
    ax_wave.fill_between(t, target - tol, target + tol, color=LIGHT,
                         label=r"Target $\pm 2.5\times10^{-3}\Phi_0$")
    ax_wave.plot(t, target, color=INK, label="Target")
    ax_wave.plot(t, unc, color=RED, label="Uncorrected")
    ax_wave.plot(t, cor, color=BLUE, ls="--", label="Cryoscope-corrected")
    ax_wave.set(xlabel="Time (ns)", ylabel=r"On-chip flux ($\Phi_0$)",
                title="Reference predistortion recovery of the CZ flux pulse")
    ax_wave.legend(frameon=False, ncol=4, loc="lower center")
    ax_wave.grid(alpha=0.18, lw=0.45)
    _panel(ax_wave, "a")

    extent = [d["durations"][0], d["durations"][-1],
              d["amplitudes"][0], d["amplitudes"][-1]]
    vmax = max(float(np.max(d["pop20_uncorrected"])),
               float(np.max(d["pop20_corrected"])))
    kwargs = dict(cmap="magma", vmin=0.0, vmax=vmax, shading="nearest")
    im = ax_unc.pcolormesh(d["durations"], d["amplitudes"],
                           d["pop20_uncorrected"], **kwargs)
    ax_cor.pcolormesh(d["durations"], d["amplitudes"],
                      d["pop20_corrected"], **kwargs)
    for ax, title, label in ((ax_unc, "Uncorrected", "b"),
                             (ax_cor, "Cryoscope-corrected", "c")):
        ax.plot(float(d["nominal_duration"][0]), float(d["nominal_amplitude"][0]),
                marker="+", ms=7, mew=1.0, color="white")
        ax.set(xlabel="Pulse duration (ns)", title=title)
        _panel(ax, label)
    ax_unc.set_ylabel(r"Flux amplitude ($\Phi_0$)")
    plt.setp(ax_cor.get_yticklabels(), visible=False)
    cbar = fig.colorbar(im, ax=[ax_unc, ax_cor], pad=0.025, fraction=0.035)
    cbar.set_label(r"Final $|20\rangle$ population")

    _save_figure(fig, "fig_E2_predistortion_cz_validation")
    plt.close(fig)


def write_notes(e1, e2):
    metrics = {
        "e1_waveform_rmse_phi0": float(e1["waveform_rmse"][0]),
        "e1_transient_gain": float(e1["transient_gain"][0]),
        "e2_nominal_infidelity_uncorrected": float(e2["nominal_infid_uncorrected"][0]),
        "e2_nominal_infidelity_corrected": float(e2["nominal_infid_corrected"][0]),
        "e2_grid_shape": list(e2["pop20_uncorrected"].shape),
    }
    (HERE / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    readme = """# Report figure draft

## E1 — Transient characterization

Finite-duration probe, control-line response, raw transient observable, waveform
reconstruction, and the transfer function estimated only in frequency bins
supported by the input spectrum. The black transfer curve is the simulator's
control-line truth; the dashed blue curve is inferred from transient reconstruction.

## E2 — Predistortion and CZ validation

The top panel compares target, distorted, and corrected on-chip CZ waveforms.
The lower panels use the same Hamiltonian, parameter grid, normalization, and
color scale. The white cross marks the frozen nominal CZ operating point.

This first draft uses Cryoscope-derived predistortion from P7. It validates the
application-side figure structure, but must not be captioned as a transient-derived
compensation result until a transient-derived filter passes independent validation.
"""
    (HERE / "README.md").write_text(readme, encoding="utf-8")


def main():
    import matplotlib
    matplotlib.use("Agg")
    _style()
    recompute = "--recompute" in sys.argv
    e1 = compute_e1(recompute=recompute)
    plot_e1(e1)
    e2 = compute_e2(recompute=recompute)
    plot_e2(e2)
    write_notes(e1, e2)
    print(f"Generated report figures in {HERE}")


if __name__ == "__main__":
    main()
