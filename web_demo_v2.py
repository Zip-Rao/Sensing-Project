#!/usr/bin/env python3
"""
Sensing-Project v2 — Interactive cQED Full-Stack Visualization & Protocol Demo

Uses the refactored sqc/ framework to provide a visual, interactive
interface for exploring superconducting qubit sensing protocols.

Run:
    python web_demo_v2.py

Requires: gradio, qutip, numpy, scipy, matplotlib
"""
from __future__ import annotations

import sys
from pathlib import Path

# Ensure project root on path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for Gradio
import matplotlib.pyplot as plt
import gradio as gr

# ── sqc/ imports ──────────────────────────────────────────────────
from sqc.devices.transmon import TransmonQubit, QubitSpec
from sqc.control.flux_signal import FluxSignal
from sqc.control.waveform import Waveform
from sqc.experiments.ramsey import RamseyExperiment
from sqc.experiments.rabi import RabiExperiment
from sqc.experiments.echo import DiffEchoExperiment
from sqc.experiments.transient import TransientSensingExperiment
from sqc.experiments.cryoscope import CryoscopeExperiment
from sqc.reconstruction.kernel import KernelEstimator
from sqc.reconstruction.wiener import WienerReconstruction
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.hardware.control_line import ControlLine
from sqc.simulation.result import ExperimentResult

# ═══════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════

def _make_qubit(EC, EJ, T1, T2, flux, n_levels):
    return TransmonQubit(
        EC=2 * np.pi * EC,
        EJ=2 * np.pi * EJ,
        T1=T1, T2=T2, flux=flux,
        n_levels=int(n_levels),
    )


def _fig_to_array(fig):
    """Convert matplotlib Figure to numpy array for Gradio."""
    fig.canvas.draw()
    data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    w, h = fig.canvas.get_width_height()
    return data.reshape(h, w, 3)


def _empty_plot(msg="No data"):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes, fontsize=14, color="gray")
    ax.set_axis_off()
    return fig


# ═══════════════════════════════════════════════════════════════════
# Tab 0: Architecture Diagram
# ═══════════════════════════════════════════════════════════════════

def render_architecture():
    """Render the cQED full-stack architecture diagram."""
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 11)
    ax.set_axis_off()
    ax.set_title("Sensing-Project — cQED Full-Stack Architecture (Gao 2021)", fontsize=16, fontweight="bold", pad=20)

    layers = [
        ("workflows/\nTop-level Pipelines", 9.5, "#e74c3c",
         "PredistortionValidation . ZCrosstalk"),
        ("calibration/\nCalibration Layer", 8.0, "#e67e22",
         "QubitFrequency . FluxResponse . TransferFunction . PredistortionDesigner"),
        ("reconstruction/\nWaveform Recovery", 6.5, "#f1c40f",
         "Wiener . Hammerstein . LM . Cryoscope . RamseyIQ"),
        ("experiments/\nExperiment Protocols", 5.0, "#2ecc71",
         "Rabi . Ramsey . Echo . Transient . Cryoscope"),
        ("simulation/\nQuTiP Interface", 3.5, "#3498db",
         "HamiltonianBuilder . MesolveRunner . SlidingMeasurementRunner"),
        ("control/\nControl Pulses", 2.0, "#9b59b6",
         "Waveform . FluxSignal . Pulse . Sequence . Gates"),
        ("hardware/\nControl Electronics", 0.5, "#1abc9c",
         "ControlLine . DistortionModel . TransferMatrix . Readout"),
    ]

    for name, y, color, details in layers:
        rect = plt.Rectangle((1, y), 12, 1.2, facecolor=color, edgecolor="white",
                              linewidth=2, alpha=0.85, zorder=2)
        ax.add_patch(rect)
        ax.text(1.3, y + 0.85, name, fontsize=12, fontweight="bold", color="white",
                va="top", zorder=3)
        ax.text(1.3, y + 0.35, details, fontsize=9, color="white",
                va="center", alpha=0.9, zorder=3)

    # devices at bottom
    rect = plt.Rectangle((1, -1.0), 12, 1.2, facecolor="#34495e", edgecolor="white",
                          linewidth=2, alpha=0.85, zorder=2)
    ax.add_patch(rect)
    ax.text(1.3, -0.15, "devices/\nPhysical Device Layer", fontsize=12, fontweight="bold",
            color="white", va="top", zorder=3)
    ax.text(1.3, -0.65, "QubitSpec . TransmonQubit . Resonator . Coupler . ChipTopology",
            fontsize=9, color="white", va="center", alpha=0.9, zorder=3)

    # dependency arrows
    for i in range(len(layers)):
        y_top = layers[i][1] + 1.2
        y_bot = layers[i][1]
        if i < len(layers) - 1:
            y_next = layers[i + 1][1] + 1.2
            ax.annotate("", xy=(6.5, y_next + 0.05), xytext=(6.5, y_top - 0.05),
                        arrowprops=dict(arrowstyle="->", color="white", lw=2))
    ax.annotate("", xy=(6.5, -1.0 + 1.2 + 0.05), xytext=(6.5, layers[-1][1] - 0.05),
                arrowprops=dict(arrowstyle="->", color="white", lw=2))

    # side labels
    ax.text(13.5, 10.1, "^ High-level", fontsize=9, color="gray", ha="center")
    ax.text(13.5, -1.5, "v Low-level", fontsize=9, color="gray", ha="center")

    return fig


# ═══════════════════════════════════════════════════════════════════
# Tab 1: Qubit Configuration
# ═══════════════════════════════════════════════════════════════════

def qubit_info(EC, EJ, T1, T2, flux, n_levels):
    try:
        q = _make_qubit(EC, EJ, T1, T2, flux, n_levels)
        f01_ghz = q.frequency / (2 * np.pi)
        alpha_mhz = q.anharmonicity * 1000 / (2 * np.pi)
        kappa = q.frequency_sensitivity(flux)
        EJ_EC = q.EJ_0 / q.EC

        text = f"""### Qubit Parameters
| Parameter | Value |
|---|---|
| f₀₁ | {f01_ghz:.4f} GHz |
| α (anharmonicity) | {alpha_mhz:.1f} MHz |
| EJ/EC | {EJ_EC:.1f} |
| κ (sensitivity) | {kappa:.4f} rad·GHz/Φ₀ |
| T₁ | {q.T1:.0f} ns |
| T₂ | {q.T2:.0f} ns |
| flux bias | {q.flux:.4f} Φ₀ |
| n_levels | {q.n_levels} |

**Sanity checks** (Gao 2021 §III.B):
- f₀₁ ∈ [4, 8] GHz: **{"✓" if 4 <= f01_ghz <= 8 else "✗"}**
- EJ/EC ∈ [40, 80]: **{"✓" if 40 <= EJ_EC <= 80 else "✗"}**
- α/h ∈ [200, 300] MHz: **{"✓" if 200 <= abs(alpha_mhz) <= 300 else "✗"}**
"""
        return text
    except Exception as e:
        return f"**Error**: {e}"


# ═══════════════════════════════════════════════════════════════════
# Tab 2: Waveform Designer
# ═══════════════════════════════════════════════════════════════════

def plot_waveform(wf_type, amplitude, frequency, center, width, rise, fall,
                  duration, n_points, noise_level):
    try:
        t = np.linspace(0, duration, int(n_points))
        type_map = {
            "Zero": 0, "Constant": 1, "Sinusoidal": 2, "Gaussian": 3,
            "Asymmetric Pulse": 4, "Double Peak": 5, "Wavepacket": 7,
        }
        phi = FluxSignal(
            type=type_map[wf_type], t_list=t,
            amplitude=amplitude, frequency=frequency,
            center=center, width=width, rise=rise, fall=fall,
            noise_level=noise_level,
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(phi.t_list, phi.samples, "b-", linewidth=1.5)
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Φ (Φ₀)")
        ax.set_title(f"Flux Signal — {wf_type}")
        ax.grid(True, alpha=0.3)
        stats = f"max={phi.samples.max():.4f}, min={phi.samples.min():.4f}, mean={phi.samples.mean():.6f}"
        ax.text(0.02, 0.95, stats, transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
        return fig
    except Exception as e:
        return _empty_plot(f"Error: {e}")


# ═══════════════════════════════════════════════════════════════════
# Tab 3: Protocol Runner
# ═══════════════════════════════════════════════════════════════════

def run_protocol(protocol, EC, EJ, T1, T2, flux, n_levels,
                 wf_type, wf_amp, wf_freq, wf_center, wf_width,
                 wf_dur, wf_npts, tau_max, n_tau):
    try:
        q = _make_qubit(EC, EJ, T1, T2, flux, n_levels)
        t = np.linspace(0, wf_dur, int(wf_npts))
        tau_list = np.linspace(0, tau_max, int(n_tau))

        type_map = {"Zero": 0, "Constant": 1, "Sinusoidal": 2, "Gaussian": 3,
                    "Asymmetric Pulse": 4, "Wavepacket": 7}
        phi = FluxSignal(type=type_map[wf_type], t_list=t, amplitude=wf_amp,
                         frequency=wf_freq, center=wf_center, width=wf_width)

        if protocol == "Rabi":
            exp = RabiExperiment(qubit=q, t_rabi=np.linspace(0, tau_max, int(n_tau)))
            result = exp.run()
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(result.axes.get("t", result.axes.get("tau", np.arange(len(result.data["p_e"])))),
                    result.data["p_e"], "b.-", markersize=3)
            ax.set_xlabel("Pulse duration (ns)"); ax.set_ylabel("p_e")
            ax.set_title("Rabi Oscillations"); ax.grid(True, alpha=0.3)

        elif protocol == "Ramsey":
            exp = RamseyExperiment(qubit=q, flux_signal=phi, tau_list=tau_list)
            result = exp.run()
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(result.axes["tau"], result.data["p_e"], "b-", linewidth=1)
            ax.set_xlabel("τ (ns)"); ax.set_ylabel("p_e")
            ax.set_title("Ramsey Interference"); ax.grid(True, alpha=0.3)

        elif protocol == "Differential Echo":
            exp = DiffEchoExperiment(qubit=q, flux_signal=phi,
                                     tau_list=tau_list, k=5)
            result = exp.run()
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(result.axes["tau"], result.data["p_e"], "r-", linewidth=1)
            ax.set_xlabel("τ (ns)"); ax.set_ylabel("p_e")
            ax.set_title("Differential Echo"); ax.grid(True, alpha=0.3)

        elif protocol == "Transient Sensing":
            exp = TransientSensingExperiment(qubit=q, flux_signal=phi)
            result = exp.run()
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4))
            ax1.plot(result.data.get("t_samples", result.axes.get("scan", [])),
                     result.data["kernel"] if "kernel" in result.data else
                     result.data.get("p_e", []), "g-", linewidth=1)
            ax1.set_title("Control Kernel"); ax1.set_xlabel("Time (ns)")
            ax1.grid(True, alpha=0.3)
            ax2.plot(result.axes.get("scan", []), result.data.get("delta_p", result.data.get("p_e", [])), "m-", linewidth=1)
            ax2.set_title("Δp vs Scan Position"); ax2.set_xlabel("Scan position")
            ax2.grid(True, alpha=0.3)

        elif protocol == "Cryoscope":
            exp = CryoscopeExperiment(qubit=q, flux_signal=phi)
            result = exp.run()
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(result.axes["trunc"], result.data["varphi"], "b.-", markersize=4)
            ax.set_xlabel("Truncation time (ns)"); ax.set_ylabel("φ (rad)")
            ax.set_title("Cryoscope: φ vs trunc"); ax.grid(True, alpha=0.3)

        else:
            return _empty_plot(f"Unknown protocol: {protocol}")

        plt.tight_layout()
        return fig

    except Exception as e:
        import traceback
        return _empty_plot(f"Error:\n{traceback.format_exc()}")


# ═══════════════════════════════════════════════════════════════════
# Tab 4: Reconstruction
# ═══════════════════════════════════════════════════════════════════

def run_reconstruction(EC, EJ, T1, T2, flux, n_levels, recon_method, lambda_reg):
    try:
        q = _make_qubit(EC, EJ, T1, T2, flux, n_levels)

        # Create a test signal and run transient sensing to get data
        t = np.linspace(0, 100, 200)
        phi_true = FluxSignal(type=3, t_list=t, amplitude=0.01, center=50, width=5)

        exp = TransientSensingExperiment(qubit=q, flux_signal=phi_true)
        result = exp.run()

        if recon_method == "Wiener Deconvolution":
            recon = WienerReconstruction(lambda_reg=lambda_reg)
            dt = result.axes["scan"][1] - result.axes["scan"][0]
            phi_rec = recon.reconstruct(
                measurement=ExperimentResult(
                    data={"delta_p": result.data["delta_p"]},
                    axes={"scan": result.axes["scan"]},
                ),
                kernel=result.data["kernel"],
                dt=dt,
            )
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4))
            ax1.plot(phi_true.t_list, phi_true.samples, "b-", label="True Φ(t)", linewidth=1.5)
            ax1.plot(phi_rec.t_list[:len(phi_rec.samples)], phi_rec.samples, "r--", label="Reconstructed", linewidth=1.5)
            ax1.set_xlabel("Time (ns)"); ax1.set_ylabel("Φ (Φ₀)")
            ax1.set_title("Wiener Reconstruction"); ax1.legend(); ax1.grid(True, alpha=0.3)
            ax2.plot(result.axes["scan"], result.data["delta_p"], "g-", linewidth=1)
            ax2.set_xlabel("Scan position"); ax2.set_ylabel("Δp")
            ax2.set_title("Measurement (Δp)"); ax2.grid(True, alpha=0.3)

        elif recon_method == "Kernel Estimation Only":
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(result.data["t_samples"], result.data["kernel"], "g-", linewidth=1)
            ax.set_xlabel("Time (ns)"); ax.set_ylabel("Kernel")
            ax.set_title("Control Kernel"); ax.grid(True, alpha=0.3)

        else:
            return _empty_plot(f"Method '{recon_method}' not available")

        plt.tight_layout()
        return fig

    except Exception as e:
        import traceback
        return _empty_plot(f"Error:\n{traceback.format_exc()}")


# ═══════════════════════════════════════════════════════════════════
# Tab 5: Distortion / Hardware Demo
# ═══════════════════════════════════════════════════════════════════

def run_distortion_demo(dist_amp, dist_tau, pulse_amp, pulse_width, pulse_center):
    try:
        t = np.linspace(0, 200, 2000)
        target = Waveform(
            t_list=t,
            samples=pulse_amp * np.exp(-0.5 * ((t - pulse_center) / pulse_width) ** 2),
        )
        dist = SingleExponentialDistortion(amplitude=dist_amp, tau=dist_tau)

        # Forward distortion
        distorted = dist.apply_to_waveform(target)

        # Predistortion via analytical inverse
        from sqc.calibration.predistortion import PredistortionDesigner
        designer = PredistortionDesigner(method="fir_inverse", n_taps=64, regularization=1e-4)
        inverse_model = designer.design(dist)
        predistorted = inverse_model.apply_to_waveform(target)
        corrected = dist.apply_to_waveform(predistorted)

        fig, axes = plt.subplots(2, 2, figsize=(14, 8))

        # Step response
        ax = axes[0, 0]
        s = dist.step_response(t - 10)
        ax.plot(t, s, "b-", linewidth=1.5)
        ax.set_title("Step Response s(t)"); ax.set_xlabel("Time (ns)")
        ax.grid(True, alpha=0.3)

        # Frequency response
        ax = axes[0, 1]
        omega = np.linspace(0, 0.5, 500)
        H = dist.frequency_response(omega)
        ax.semilogy(omega, np.abs(H), "r-", linewidth=1.5)
        ax.set_title("|H(ω)|"); ax.set_xlabel("ω (rad/ns)")
        ax.grid(True, alpha=0.3)

        # AWG → Chip (without predistortion)
        ax = axes[1, 0]
        ax.plot(t, target.samples, "b-", label="Target (desired)", linewidth=1.5, alpha=0.7)
        ax.plot(t, distorted.samples, "r--", label="On-chip (distorted)", linewidth=1.5)
        ax.set_title("Without Predistortion"); ax.set_xlabel("Time (ns)")
        ax.legend(); ax.grid(True, alpha=0.3)

        # With predistortion
        ax = axes[1, 1]
        ax.plot(t, target.samples, "b-", label="Target", linewidth=1.5, alpha=0.7)
        ax.plot(t, corrected.samples, "g--", label="On-chip (corrected)", linewidth=1.5)
        rmse = np.sqrt(np.mean((corrected.samples - target.samples) ** 2))
        ax.set_title(f"With Predistortion (RMSE={rmse:.2e})"); ax.set_xlabel("Time (ns)")
        ax.legend(); ax.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig

    except Exception as e:
        import traceback
        return _empty_plot(f"Error:\n{traceback.format_exc()}")


# ═══════════════════════════════════════════════════════════════════
# Gradio UI
# ═══════════════════════════════════════════════════════════════════

THEME = gr.themes.Soft(primary_hue="blue", secondary_hue="cyan")

with gr.Blocks(theme=THEME, title="Sensing-Project v2 — cQED Full-Stack Demo") as demo:
    gr.Markdown("""
    # 🔬 Sensing-Project v2 — Superconducting Quantum Control Platform

    **Full-stack cQED simulation framework** for flux-tunable Transmon qubit sensing,
    waveform reconstruction, calibration, and predistortion.
    Built on the Gao 2021 six-layer architecture.
    """)

    with gr.Tabs():
        # ── Tab 0: Architecture ──
        with gr.Tab("🏗️ Architecture"):
            gr.Markdown("## cQED Full-Stack Architecture")
            arch_plot = gr.Plot(render_architecture(), label="architecture_diagram")
            gr.Markdown("""
            ### Layer Descriptions
            | Layer | Module | Responsibility |
            |---|---|---|
            | **Workflows** | `sqc/workflows/` | End-to-end research pipelines |
            | **Calibration** | `sqc/calibration/` | Qubit & transfer function calibration |
            | **Reconstruction** | `sqc/reconstruction/` | Waveform recovery from measurements |
            | **Experiments** | `sqc/experiments/` | Sensing protocol objects |
            | **Simulation** | `sqc/simulation/` | QuTiP mesolve interface |
            | **Control** | `sqc/control/` | Waveforms, pulses, sequences |
            | **Hardware** | `sqc/hardware/` | Control lines, distortion, crosstalk |
            | **Devices** | `sqc/devices/` | Qubit, resonator, chip models |

            *Reference: Gao, Rol, Touzard, Wang, PRX Quantum 2, 040202 (2021)*
            """)

        # ── Tab 1: Qubit Config ──
        with gr.Tab("⚛️ Qubit Configuration"):
            gr.Markdown("## Transmon Qubit Parameters")
            with gr.Row():
                with gr.Column(scale=1):
                    ec = gr.Slider(0.1, 0.5, 0.2, step=0.01, label="EC (GHz·2π)")
                    ej = gr.Slider(5.0, 25.0, 15.0, step=0.5, label="EJ (GHz·2π)")
                    t1 = gr.Slider(1000, 50000, 10000, step=1000, label="T₁ (ns)")
                    t2 = gr.Slider(1000, 50000, 8000, step=1000, label="T₂ (ns)")
                    flux_bias = gr.Slider(-0.5, 0.5, 0.0, step=0.01, label="Flux bias (Φ₀)")
                    nlev = gr.Radio([2, 3, 4, 5], value=3, label="n_levels (Fock truncation)")
                with gr.Column(scale=1):
                    qubit_output = gr.Markdown("*(Configure parameters to see qubit info)*")

            gr.Button("Refresh Qubit Info").click(
                qubit_info, [ec, ej, t1, t2, flux_bias, nlev], qubit_output
            )

        # ── Tab 2: Waveform Designer ──
        with gr.Tab("📈 Waveform Designer"):
            gr.Markdown("## Flux Signal Designer")
            with gr.Row():
                with gr.Column(scale=1):
                    wf_type = gr.Dropdown(
                        ["Zero", "Constant", "Sinusoidal", "Gaussian",
                         "Asymmetric Pulse", "Double Peak", "Wavepacket"],
                        value="Gaussian", label="Signal Type"
                    )
                    wf_amp = gr.Slider(0.0, 0.1, 0.01, step=0.001, label="Amplitude (Φ₀)")
                    wf_freq = gr.Slider(0.0, 0.1, 0.01, step=0.001, label="Frequency (for sinusoidal)")
                    wf_center = gr.Slider(0, 200, 50, step=1, label="Center (ns)")
                    wf_width = gr.Slider(1, 50, 10, step=1, label="Width (ns)")
                    wf_rise = gr.Slider(1, 50, 10, step=1, label="Rise time (ns)")
                    wf_fall = gr.Slider(1, 50, 10, step=1, label="Fall time (ns)")
                    wf_dur = gr.Slider(50, 500, 250, step=10, label="Duration (ns)")
                    wf_npts = gr.Slider(100, 2000, 500, step=50, label="N points")
                    wf_noise = gr.Slider(0.0, 0.01, 0.0, step=0.0001, label="Noise level")
                with gr.Column(scale=2):
                    wf_plot = gr.Plot(label="waveform_plot")

            gr.Button("Generate Waveform").click(
                plot_waveform,
                [wf_type, wf_amp, wf_freq, wf_center, wf_width, wf_rise, wf_fall, wf_dur, wf_npts, wf_noise],
                wf_plot,
            )

        # ── Tab 3: Protocol Runner ──
        with gr.Tab("🧪 Protocol Runner"):
            gr.Markdown("## Execute Sensing Protocols")
            with gr.Row():
                with gr.Column(scale=1):
                    protocol = gr.Dropdown(
                        ["Rabi", "Ramsey", "Differential Echo", "Transient Sensing", "Cryoscope"],
                        value="Ramsey", label="Protocol"
                    )
                    tau_max = gr.Slider(20, 500, 250, step=10, label="τ max / scan range (ns)")
                    n_tau = gr.Slider(20, 1000, 200, step=20, label="N points (τ/scan)")
                    gr.Markdown("*Qubit and waveform params from Tabs 1 & 2 are shared*")
                    run_btn = gr.Button("▶ Run Protocol", variant="primary")
                with gr.Column(scale=2):
                    proto_plot = gr.Plot(label="protocol_result")

            run_btn.click(
                run_protocol,
                [protocol, ec, ej, t1, t2, flux_bias, nlev,
                 wf_type, wf_amp, wf_freq, wf_center, wf_width, wf_dur, wf_npts, tau_max, n_tau],
                proto_plot,
            )

        # ── Tab 4: Reconstruction ──
        with gr.Tab("🔄 Reconstruction"):
            gr.Markdown("## Waveform Reconstruction")
            with gr.Row():
                with gr.Column(scale=1):
                    recon_method = gr.Dropdown(
                        ["Wiener Deconvolution", "Kernel Estimation Only"],
                        value="Wiener Deconvolution", label="Method"
                    )
                    lambda_reg = gr.Slider(0.001, 100.0, 1.0, step=0.1, label="Regularization λ")
                    recon_btn = gr.Button("▶ Run Reconstruction", variant="primary")
                with gr.Column(scale=2):
                    recon_plot = gr.Plot(label="reconstruction_result")

            recon_btn.click(
                run_reconstruction,
                [ec, ej, t1, t2, flux_bias, nlev, recon_method, lambda_reg],
                recon_plot,
            )

        # ── Tab 5: Distortion / Hardware ──
        with gr.Tab("🔧 Distortion & Predistortion"):
            gr.Markdown("## Control Line Distortion & Predistortion Demo")
            with gr.Row():
                with gr.Column(scale=1):
                    dist_amp = gr.Slider(0.0, 0.2, 0.05, step=0.005, label="Distortion amplitude")
                    dist_tau = gr.Slider(5, 200, 50, step=5, label="Distortion τ (ns)")
                    pulse_amp2 = gr.Slider(0.0, 2.0, 1.0, step=0.1, label="Pulse amplitude")
                    pulse_w = gr.Slider(1, 30, 8, step=1, label="Pulse width (ns)")
                    pulse_c = gr.Slider(0, 200, 80, step=5, label="Pulse center (ns)")
                    dist_btn = gr.Button("▶ Run Distortion Demo", variant="primary")
                with gr.Column(scale=2):
                    dist_plot = gr.Plot(label="distortion_result")

            dist_btn.click(
                run_distortion_demo,
                [dist_amp, dist_tau, pulse_amp2, pulse_w, pulse_c],
                dist_plot,
            )

        # ── Tab 6: About ──
        with gr.Tab("📖 About"):
            gr.Markdown("""
            ## Sensing-Project v2

            **Version**: sqc v0.1.0
            **Framework**: QuTiP + NumPy + SciPy + Matplotlib
            **Reference**: Gao et al., *PRX Quantum* 2, 040202 (2021)

            ### Key Features
            - Six-layer cQED full-stack architecture
            - Flux-tunable Transmon qubit model with SQUID modulation
            - Five sensing protocols: Rabi, Ramsey, Differential Echo, Transient, Cryoscope
            - Waveform reconstruction: Wiener, Hammerstein-Wiener, Levenberg-Marquardt
            - Control line distortion modeling & predistortion
            - Multi-qubit Z-crosstalk transfer matrix
            - Physics regression testing with frozen baselines

            ### Documentation
            - Technical docs: `docs/architecture.md`
            - Refactoring plan: `idea/refactor/_refactor_plan.md`
            - Test suite: `pytest tests/ -v`

            ### Environment
            Conda env: `qutip-env` | Python 3.x | qutip ≥ 5.0
            """)

# ── Entry point ──
if __name__ == "__main__":
    print("=" * 60)
    print("Sensing-Project v2 — cQED Full-Stack Demo")
    print(f"sqc version: 0.1.0")
    print("=" * 60)
    demo.launch(
        server_name="127.0.0.1",
        server_port=7860,
        share=False,
        show_error=True,
    )
