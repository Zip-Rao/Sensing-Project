#!/usr/bin/env python3
"""
Sensing-Project — Interactive cQED Full-Stack Visualization & Protocol Demo

A polished Gradio interface for the refactored sqc/ framework, providing
visual exploration of the six-layer cQED architecture and interactive
sensing-protocol experimentation.

Run:
    python web_demo_v2.py

Requires: gradio>=4.0, qutip>=5.0, numpy, scipy, matplotlib
"""
from __future__ import annotations

import sys
import io
import contextlib
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


@contextlib.contextmanager
def suppress_stdout():
    """Silence noisy debug prints from legacy src/protocal.py during runs."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        yield

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
import gradio as gr

# ─── sqc/ imports ────────────────────────────────────────────────────────
from sqc.devices.transmon import TransmonQubit, QubitSpec
from sqc.control.flux_signal import FluxSignal
from sqc.control.waveform import Waveform
from sqc.experiments.ramsey import RamseyExperiment
from sqc.experiments.rabi import RabiExperiment
from sqc.experiments.echo import DiffEchoExperiment
from sqc.experiments.transient import TransientSensingExperiment
from sqc.experiments.cryoscope import CryoscopeExperiment
from sqc.reconstruction.kernel import KernelEstimator
from sqc.reconstruction.transient import TransientReconstruction
from sqc.hardware.distortion import (
    SingleExponentialDistortion, MultiExponentialDistortion, FIRDistortion,
)
from sqc.hardware.control_line import ControlLine
from sqc.hardware.transfer_matrix import TransferMatrix
from sqc.simulation.result import ExperimentResult
from sqc.simulation.hamiltonian import HamiltonianBuilder

# ═══════════════════════════════════════════════════════════════════════════
# Global styling
# ═══════════════════════════════════════════════════════════════════════════

# Color palette (semantic, keyed to cQED stack layers)
COLORS = {
    "workflows":     "#E63946",  # red
    "calibration":   "#F77F00",  # orange
    "reconstruction":"#FFB627",  # amber
    "experiments":   "#06A77D",  # green
    "simulation":    "#1D70A2",  # blue
    "control":       "#7B2CBF",  # purple
    "hardware":      "#00B4A6",  # teal
    "devices":       "#3D405B",  # navy
    # Plot accents
    "primary":       "#1D70A2",
    "secondary":     "#E63946",
    "accent":        "#06A77D",
    "warning":       "#F77F00",
    "neutral":       "#6C757D",
    "bg":            "#F8F9FA",
    "fg":            "#1A202C",
}

# Configure matplotlib for cleaner plots
plt.rcParams.update({
    "figure.facecolor":   "white",
    "axes.facecolor":     "#FAFBFC",
    "axes.edgecolor":     "#CBD5E0",
    "axes.labelcolor":    "#2D3748",
    "axes.titlecolor":    "#1A202C",
    "axes.titlesize":     12,
    "axes.titleweight":   "bold",
    "axes.titlepad":      10,
    "axes.labelsize":     10,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.color":         "#E2E8F0",
    "grid.linestyle":     "--",
    "grid.alpha":         0.7,
    "xtick.color":        "#4A5568",
    "ytick.color":        "#4A5568",
    "xtick.labelsize":    9,
    "ytick.labelsize":    9,
    "legend.frameon":     True,
    "legend.facecolor":   "white",
    "legend.edgecolor":   "#CBD5E0",
    "legend.fontsize":    9,
    "font.family":        ["DejaVu Sans", "Arial", "sans-serif"],
    "lines.linewidth":    1.8,
})

# ═══════════════════════════════════════════════════════════════════════════
# Custom CSS for Gradio
# ═══════════════════════════════════════════════════════════════════════════

CUSTOM_CSS = """
/* --- Global container --- */
.gradio-container {
    max-width: 1400px !important;
    margin: 0 auto;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
}

/* --- Header --- */
#header {
    background: linear-gradient(135deg, #1D70A2 0%, #06A77D 100%);
    color: white;
    padding: 28px 32px;
    border-radius: 12px;
    margin-bottom: 16px;
    box-shadow: 0 4px 12px rgba(0,0,0,0.08);
}
#header h1 {
    margin: 0 0 8px 0;
    font-size: 28px;
    font-weight: 700;
}
#header p {
    margin: 0;
    opacity: 0.95;
    font-size: 14px;
    line-height: 1.6;
}

/* --- Layer badges --- */
.layer-badge {
    display: inline-block;
    padding: 3px 10px;
    border-radius: 12px;
    font-size: 11px;
    font-weight: 600;
    color: white;
    margin: 2px 4px 2px 0;
    letter-spacing: 0.3px;
}
.badge-workflows     { background: #E63946; }
.badge-calibration   { background: #F77F00; }
.badge-reconstruction{ background: #FFB627; color: #1A202C; }
.badge-experiments   { background: #06A77D; }
.badge-simulation    { background: #1D70A2; }
.badge-control       { background: #7B2CBF; }
.badge-hardware      { background: #00B4A6; }
.badge-devices       { background: #3D405B; }

/* --- Cards --- */
.info-card {
    background: white;
    border: 1px solid #E2E8F0;
    border-radius: 10px;
    padding: 16px 20px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04);
    margin: 8px 0;
}
.info-card h3 {
    margin-top: 0;
    color: #1D70A2;
    border-bottom: 2px solid #E2E8F0;
    padding-bottom: 8px;
}

/* --- Tab navigation --- */
button.svelte-1uw5tnk {
    font-weight: 600 !important;
}

/* --- Sliders --- */
input[type="range"] {
    accent-color: #1D70A2;
}

/* --- Primary buttons --- */
.primary-btn button {
    background: linear-gradient(135deg, #1D70A2 0%, #06A77D 100%) !important;
    color: white !important;
    font-weight: 600 !important;
    border: none !important;
    box-shadow: 0 2px 6px rgba(29,112,162,0.25) !important;
}
.primary-btn button:hover {
    box-shadow: 0 4px 12px rgba(29,112,162,0.4) !important;
    transform: translateY(-1px);
}

/* --- Status pills --- */
.status-ok    { color: #06A77D; font-weight: 600; }
.status-warn  { color: #F77F00; font-weight: 600; }
.status-fail  { color: #E63946; font-weight: 600; }

/* --- Footer --- */
#footer {
    text-align: center;
    padding: 16px;
    color: #6C757D;
    font-size: 12px;
    border-top: 1px solid #E2E8F0;
    margin-top: 24px;
}
"""

# ═══════════════════════════════════════════════════════════════════════════
# Helper utilities
# ═══════════════════════════════════════════════════════════════════════════

def _make_qubit(EC, EJ, T1, T2, flux, n_levels):
    """Construct a TransmonQubit from UI parameters."""
    return TransmonQubit(
        EC=2 * np.pi * EC,
        EJ=2 * np.pi * EJ,
        T1=T1, T2=T2, flux=flux,
        n_levels=int(n_levels),
    )


def _empty_plot(msg="No data"):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes,
            fontsize=13, color=COLORS["neutral"])
    ax.set_axis_off()
    return fig


def _layer_badges(*layers):
    """Render HTML badges for layers used by current operation."""
    html = '<div style="margin: 8px 0;">📦 <b>Stack layers used:</b> '
    for layer in layers:
        html += f'<span class="layer-badge badge-{layer}">{layer}</span>'
    html += '</div>'
    return html


def _status_pill(label, ok):
    cls = "status-ok" if ok else "status-fail"
    icon = "✓" if ok else "✗"
    return f'<span class="{cls}">{icon} {label}</span>'


# ═══════════════════════════════════════════════════════════════════════════
# TAB 0 — Architecture (layered diagram + data flow)
# ═══════════════════════════════════════════════════════════════════════════

def render_architecture_layered():
    """Detailed layered view of the cQED stack."""
    fig, ax = plt.subplots(figsize=(13, 9))
    ax.set_xlim(0, 13)
    ax.set_ylim(-0.5, 9.5)
    ax.set_axis_off()
    ax.set_title("cQED Full-Stack Architecture (Gao 2021 Fig. 1a)",
                 fontsize=15, fontweight="bold", pad=12, color=COLORS["fg"])

    layers = [
        ("workflows",     "Top-level Pipelines",   8.4,
         "PredistortionValidation  •  ZCrosstalk",                        "🔬"),
        ("calibration",   "Calibration Workflows", 7.2,
         "FluxResponse  •  SinglePointFrequency  •  Waveform  •  PredistortionDesigner", "🎯"),
        ("reconstruction","Waveform Recovery",     6.0,
         "Wiener  •  Hammerstein  •  LM  •  Cryoscope  •  RamseyIQ  •  KernelEstimator", "🔄"),
        ("experiments",   "Sensing Protocols",     4.8,
         "Rabi  •  Ramsey  •  DiffEcho  •  TransientSensing  •  Cryoscope", "🧪"),
        ("simulation",    "QuTiP Interface",       3.6,
         "HamiltonianBuilder  •  MesolveRunner  •  SlidingMeasurementRunner  •  noise", "⚙️"),
        ("control",       "Control Pulses",        2.4,
         "Waveform  •  FluxSignal  •  Pulse  •  Sequence  •  Gates",      "📈"),
        ("hardware",      "Control Electronics",   1.2,
         "ControlLine  •  DistortionModel  •  TransferMatrix  •  Readout", "🔧"),
        ("devices",       "Physical Devices",      0.0,
         "QubitSpec  •  TransmonQubit  •  Resonator  •  Coupler  •  ChipTopology", "⚛️"),
    ]

    box_w, box_h = 10.5, 1.05
    x_box = 1.2

    for layer_key, name, y, members, icon in layers:
        color = COLORS[layer_key]
        # rounded box
        box = FancyBboxPatch(
            (x_box, y), box_w, box_h,
            boxstyle="round,pad=0.02,rounding_size=0.12",
            facecolor=color, edgecolor="white", linewidth=2.5, alpha=0.95,
            zorder=2,
        )
        ax.add_patch(box)

        # icon
        ax.text(x_box + 0.4, y + box_h/2, icon, fontsize=22,
                ha="center", va="center", zorder=3)

        # layer name
        ax.text(x_box + 0.95, y + box_h - 0.3,
                f"sqc/{layer_key}/", fontsize=11, fontweight="bold",
                color="white", va="center", zorder=3)
        ax.text(x_box + 0.95, y + box_h - 0.65,
                name, fontsize=9, color="white", alpha=0.95,
                va="center", zorder=3, style="italic")

        # members
        ax.text(x_box + 0.95, y + 0.22, members,
                fontsize=8, color="white", alpha=0.9, va="center", zorder=3)

    # dependency arrow on the left
    arrow = FancyArrowPatch(
        (0.5, 8.9), (0.5, 0.1),
        arrowstyle="->,head_width=8,head_length=10",
        color=COLORS["neutral"], lw=2.5, zorder=1,
    )
    ax.add_patch(arrow)
    ax.text(0.18, 4.5, "depends on →", rotation=90, fontsize=10,
            color=COLORS["neutral"], ha="center", va="center", style="italic")

    # Reference label
    ax.text(12.6, -0.4, "Ref: Gao et al., PRX Quantum 2, 040202 (2021)",
            fontsize=8, color=COLORS["neutral"], ha="right", style="italic")

    fig.tight_layout()
    return fig


def render_data_flow():
    """Visual data flow diagram showing how a measurement becomes a recovered waveform."""
    fig, ax = plt.subplots(figsize=(13, 7))
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 7)
    ax.set_axis_off()
    ax.set_title("Example Data Flow: Sensing → Reconstruction Pipeline",
                 fontsize=14, fontweight="bold", pad=10, color=COLORS["fg"])

    nodes = [
        # (x, y, w, h, title, subtitle, color_key)
        (0.4,  4.7, 2.2, 1.1, "FluxSignal Φ(t)",  "control/flux_signal",        "control"),
        (3.0,  4.7, 2.2, 1.1, "Pulse Sequence",    "control/sequence",           "control"),
        (5.6,  4.7, 2.4, 1.1, "QubitSpec",         "devices/transmon",           "devices"),
        (8.4,  4.7, 2.4, 1.1, "ControlLine",       "hardware/control_line",      "hardware"),

        (3.5,  2.8, 2.4, 1.1, "HamiltonianBuilder","simulation/hamiltonian",     "simulation"),
        (6.5,  2.8, 2.4, 1.1, "MesolveRunner",     "simulation/runner",          "simulation"),

        (1.0,  0.6, 2.4, 1.1, "ExperimentResult",  "{p_e, axes, metadata}",      "experiments"),
        (4.2,  0.6, 2.4, 1.1, "KernelEstimator",   "reconstruction/kernel",      "reconstruction"),
        (7.4,  0.6, 2.4, 1.1, "WienerReconstr.",   "reconstruction/wiener",      "reconstruction"),
        (10.6, 0.6, 2.0, 1.1, "Recovered Φ̂(t)",   "FluxSignal output",           "reconstruction"),
    ]

    for x, y, w, h, title, sub, key in nodes:
        color = COLORS[key]
        box = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.02,rounding_size=0.1",
            facecolor=color, edgecolor="white", linewidth=2, alpha=0.92,
            zorder=2,
        )
        ax.add_patch(box)
        ax.text(x + w/2, y + h - 0.32, title, fontsize=10, fontweight="bold",
                color="white", ha="center", va="center", zorder=3)
        ax.text(x + w/2, y + 0.28, sub, fontsize=8, color="white",
                alpha=0.85, ha="center", va="center", style="italic", zorder=3)

    # arrows
    def arrow(x1, y1, x2, y2, **kw):
        a = FancyArrowPatch((x1, y1), (x2, y2),
                            arrowstyle="->,head_width=6,head_length=8",
                            color=COLORS["neutral"], lw=1.5, zorder=1, **kw)
        ax.add_patch(a)

    # row 1 → row 2
    arrow(1.5, 4.7, 4.7, 3.9)   # FluxSignal → HamiltonianBuilder
    arrow(4.1, 4.7, 4.7, 3.9)   # Pulse → HamiltonianBuilder
    arrow(6.8, 4.7, 5.0, 3.9)   # QubitSpec → HamiltonianBuilder
    arrow(9.6, 4.7, 7.7, 3.9)   # ControlLine → MesolveRunner

    # row 2 connection
    arrow(5.9, 3.35, 6.5, 3.35)  # HamiltonianBuilder → MesolveRunner

    # row 2 → row 3
    arrow(4.7, 2.8, 2.2, 1.7)  # HamiltonianBuilder → ExperimentResult
    arrow(7.7, 2.8, 5.4, 1.7)  # MesolveRunner → KernelEstimator (kernel)
    arrow(7.7, 2.8, 8.6, 1.7)  # MesolveRunner → ExperimentResult (Δp)

    # row 3 chain
    arrow(3.4, 1.15, 4.2, 1.15)
    arrow(6.6, 1.15, 7.4, 1.15)
    arrow(9.4, 1.15, 10.6, 1.15)

    # data type labels
    ax.text(2.8, 3.95, "H_list, t_list", fontsize=8, color=COLORS["neutral"],
            ha="center", style="italic")
    ax.text(8.5, 2.5, "QuTiP\nResult", fontsize=8, color=COLORS["neutral"],
            ha="center", style="italic")
    ax.text(3.8, 1.4, "Δp(scan)", fontsize=8, color=COLORS["neutral"], style="italic")
    ax.text(7.0, 1.4, "kernel(t)", fontsize=8, color=COLORS["neutral"], style="italic")
    ax.text(10.2, 1.4, "deconvolve", fontsize=8, color=COLORS["neutral"], style="italic")

    # legend at bottom
    layers_to_show = ["control", "devices", "hardware", "simulation", "experiments", "reconstruction"]
    handles = [Line2D([], [], marker="s", linestyle="", markersize=12,
                       markerfacecolor=COLORS[k], markeredgecolor="white", label=k)
               for k in layers_to_show]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              ncol=6, frameon=False, fontsize=9)

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# TAB 1 — Qubit Configuration (live updates)
# ═══════════════════════════════════════════════════════════════════════════

def qubit_info_html(EC, EJ, T1, T2, flux, n_levels):
    """Render qubit info as styled HTML (auto-refresh)."""
    try:
        q = _make_qubit(EC, EJ, T1, T2, flux, n_levels)
        f01_ghz = q.frequency / (2 * np.pi)
        alpha_mhz = abs(q.anharmonicity) * 1000 / (2 * np.pi)
        kappa = q.frequency_sensitivity(flux)
        EJ_EC = q.EJ_0 / q.EC

        ok_f = 4 <= f01_ghz <= 8
        ok_ec = 40 <= EJ_EC <= 80
        ok_a = 200 <= alpha_mhz <= 300

        all_ok = ok_f and ok_ec and ok_a
        status_color = COLORS["accent"] if all_ok else COLORS["warning"]
        status_text = "✓ All parameters in recommended range" if all_ok else "⚠ Some parameters out of range"

        html = f"""
<div class="info-card">
  <h3>⚛️ Derived Quantities</h3>
  <table style="width:100%; border-collapse: collapse; font-size:13px;">
    <tr style="background:#F1F5F9;">
      <th style="text-align:left; padding:8px; border-bottom:2px solid #CBD5E0;">Parameter</th>
      <th style="text-align:right; padding:8px; border-bottom:2px solid #CBD5E0;">Value</th>
      <th style="text-align:left; padding:8px; border-bottom:2px solid #CBD5E0;">Units</th>
    </tr>
    <tr><td style="padding:6px;">f₀₁ (transition frequency)</td>
        <td style="text-align:right; padding:6px; font-family:monospace;"><b>{f01_ghz:.4f}</b></td>
        <td style="padding:6px; color:#6C757D;">GHz</td></tr>
    <tr style="background:#F8FAFC;"><td style="padding:6px;">α (anharmonicity)</td>
        <td style="text-align:right; padding:6px; font-family:monospace;"><b>−{alpha_mhz:.1f}</b></td>
        <td style="padding:6px; color:#6C757D;">MHz</td></tr>
    <tr><td style="padding:6px;">EJ/EC ratio</td>
        <td style="text-align:right; padding:6px; font-family:monospace;"><b>{EJ_EC:.1f}</b></td>
        <td style="padding:6px; color:#6C757D;">—</td></tr>
    <tr style="background:#F8FAFC;"><td style="padding:6px;">κ = dω/dΦ (sensitivity)</td>
        <td style="text-align:right; padding:6px; font-family:monospace;"><b>{kappa:.4f}</b></td>
        <td style="padding:6px; color:#6C757D;">rad·GHz/Φ₀</td></tr>
    <tr><td style="padding:6px;">EJ(Φ) at current flux</td>
        <td style="text-align:right; padding:6px; font-family:monospace;"><b>{q.EJ/(2*np.pi):.3f}</b></td>
        <td style="padding:6px; color:#6C757D;">GHz·2π</td></tr>
    <tr style="background:#F8FAFC;"><td style="padding:6px;">Hilbert dim</td>
        <td style="text-align:right; padding:6px; font-family:monospace;"><b>{q.n_levels}</b></td>
        <td style="padding:6px; color:#6C757D;">levels</td></tr>
  </table>
</div>

<div class="info-card" style="border-left: 4px solid {status_color};">
  <h3>📋 Sanity Checks <span style="color:{status_color}; font-size:14px;">— {status_text}</span></h3>
  <ul style="line-height:1.8; margin:8px 0;">
    <li>f₀₁ ∈ [4, 8] GHz: {_status_pill(f"{f01_ghz:.3f} GHz", ok_f)}</li>
    <li>EJ/EC ∈ [40, 80]: {_status_pill(f"{EJ_EC:.1f}", ok_ec)}</li>
    <li>|α|/h ∈ [200, 300] MHz: {_status_pill(f"{alpha_mhz:.0f} MHz", ok_a)}</li>
  </ul>
  <p style="font-size:11px; color:#6C757D; margin:8px 0 0 0;">
    Recommended ranges from Gao 2021 §III.B for Transmon qubit design.
  </p>
</div>
"""
        return html
    except Exception as e:
        return f'<div class="info-card" style="border-left:4px solid #E63946;"><h3>Error</h3><pre>{e}</pre></div>'


def qubit_spectrum_plot(EC, EJ, T1, T2, flux, n_levels):
    """Plot f01(Φ) and κ(Φ) showing the current operating point."""
    try:
        flux_range = np.linspace(-0.5, 0.5, 401)
        spec = QubitSpec(name="Q", EC=2*np.pi*EC, EJ=2*np.pi*EJ,
                         T1=T1, T2=T2, n_levels=int(n_levels))
        f_curve = np.array([spec.frequency(f)/(2*np.pi) for f in flux_range])
        kappa_curve = np.array([spec.sensitivity(f) for f in flux_range])

        f_now = spec.frequency(flux)/(2*np.pi)
        kappa_now = spec.sensitivity(flux)

        fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
        ax = axes[0]
        ax.plot(flux_range, f_curve, color=COLORS["primary"], linewidth=2)
        ax.axvline(flux, color=COLORS["secondary"], linestyle="--", linewidth=1.2, alpha=0.8)
        ax.axhline(f_now, color=COLORS["secondary"], linestyle=":", linewidth=1, alpha=0.5)
        ax.scatter([flux], [f_now], color=COLORS["secondary"], s=80, zorder=5,
                   edgecolor="white", linewidth=1.5,
                   label=f"current: f={f_now:.3f} GHz")
        ax.set_xlabel("Flux Φ (Φ₀)")
        ax.set_ylabel("f₀₁ (GHz)")
        ax.set_title("Qubit Frequency vs Flux  —  f₀₁(Φ) = √(8 EJ(Φ) EC) − EC")
        ax.legend(loc="best")

        ax = axes[1]
        ax.plot(flux_range, kappa_curve, color=COLORS["accent"], linewidth=2)
        ax.axvline(flux, color=COLORS["secondary"], linestyle="--", linewidth=1.2, alpha=0.8)
        ax.scatter([flux], [kappa_now], color=COLORS["secondary"], s=80, zorder=5,
                   edgecolor="white", linewidth=1.5,
                   label=f"current: κ={kappa_now:.3f}")
        ax.set_xlabel("Flux Φ (Φ₀)")
        ax.set_ylabel("κ = dω₀₁/dΦ  (rad·GHz/Φ₀)")
        ax.set_title("Flux Sensitivity Spectrum")
        ax.legend(loc="best")

        fig.tight_layout()
        return fig
    except Exception as e:
        return _empty_plot(f"Error: {e}")


# ═══════════════════════════════════════════════════════════════════════════
# TAB 2 — Waveform Designer (live updates)
# ═══════════════════════════════════════════════════════════════════════════

WF_TYPE_MAP = {
    "Zero": 0, "Constant": 1, "Sinusoidal": 2, "Gaussian": 3,
    "Asymmetric Pulse": 4, "Double Peak": 5, "Wavepacket": 7,
}


def plot_waveform(wf_type, amplitude, frequency, center, width, rise, fall,
                  duration, n_points, noise_level):
    try:
        t = np.linspace(0, duration, int(n_points))
        phi = FluxSignal(
            type=WF_TYPE_MAP[wf_type], t_list=t,
            amplitude=amplitude, frequency=frequency,
            center=center, width=width, rise=rise, fall=fall,
            noise_level=noise_level,
        )

        fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))

        # Time domain
        ax = axes[0]
        ax.fill_between(phi.t_list, phi.samples, alpha=0.18, color=COLORS["primary"])
        ax.plot(phi.t_list, phi.samples, color=COLORS["primary"], linewidth=1.6)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Φ (Φ₀)")
        ax.set_title(f"Time Domain  —  {wf_type}")
        stats = (f"max = {phi.samples.max():.4f} Φ₀\n"
                 f"min = {phi.samples.min():.4f} Φ₀\n"
                 f"mean = {phi.samples.mean():.6f} Φ₀\n"
                 f"RMS = {np.sqrt(np.mean(phi.samples**2)):.4f} Φ₀")
        ax.text(0.97, 0.97, stats, transform=ax.transAxes, fontsize=8.5,
                va="top", ha="right", family="monospace",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor=COLORS["neutral"], alpha=0.92))

        # Frequency domain (FFT magnitude)
        ax = axes[1]
        if len(t) > 4 and (t[1] - t[0]) > 0:
            dt = t[1] - t[0]
            spectrum = np.abs(np.fft.rfft(phi.samples)) * dt
            freqs = np.fft.rfftfreq(len(t), d=dt) * 1000  # MHz
            ax.semilogy(freqs[1:], spectrum[1:] + 1e-12,
                        color=COLORS["accent"], linewidth=1.4)
            ax.set_xlabel("Frequency (MHz)")
            ax.set_ylabel("|FFT|  (Φ₀·ns)")
            ax.set_title("Frequency Spectrum")
            ax.set_xlim(0, min(freqs[-1], 100))
        fig.tight_layout()
        return fig
    except Exception as e:
        return _empty_plot(f"Error: {e}")


# ═══════════════════════════════════════════════════════════════════════════
# TAB 3 — Protocol Runner (richer panels)
# ═══════════════════════════════════════════════════════════════════════════

PROTOCOL_INFO = {
    "Rabi": {
        "desc": "Resonant driving sweep — scans pulse duration, observes p_e oscillations.",
        "layers": ["devices", "control", "simulation", "experiments"],
        "physics": "Ω·t rotation on the Bloch sphere; period = 2π/Ω.",
    },
    "Ramsey": {
        "desc": "π/2 — τ — π/2 sequence; measures phase accumulated under flux signal.",
        "layers": ["devices", "control", "simulation", "experiments"],
        "physics": "p_e(τ) = ½(1 − cos[∫κ·Φ(t')dt']); decays via T₂*.",
    },
    "Differential Echo": {
        "desc": "π/2 — [Hahn echo]^k — π/2 (k=5 default). Cancels low-frequency noise.",
        "layers": ["devices", "control", "simulation", "experiments"],
        "physics": "B = −φ / (2k·κ·t_int) — direct closed-form retrieval.",
    },
    "Transient Sensing": {
        "desc": "Sliding π/2-pulse measurement; extracts kernel + Δp signature.",
        "layers": ["devices", "control", "simulation", "experiments", "reconstruction"],
        "physics": "Δp(t) = ∫kernel(t−t')·Φ(t')dt'  → Wiener-decoded.",
    },
    "Cryoscope": {
        "desc": "Truncated flux-pulse + IQ readout; recovers φ(t_d) for cryoscope inversion.",
        "layers": ["devices", "control", "simulation", "experiments"],
        "physics": "φ(t_d) = ∫₀^t_d κ·Φ(t')dt'; differentiate → Φ(t).",
    },
}


def protocol_info_html(protocol):
    info = PROTOCOL_INFO.get(protocol, {})
    layers_html = "".join(f'<span class="layer-badge badge-{l}">{l}</span>'
                           for l in info.get("layers", []))
    return f"""
<div class="info-card">
  <h3>🧪 {protocol}</h3>
  <p style="margin:6px 0;">{info.get('desc','')}</p>
  <p style="margin:6px 0; color:#4A5568;"><b>Physics:</b> <code>{info.get('physics','')}</code></p>
  <p style="margin:6px 0;">{layers_html}</p>
</div>
"""


def run_protocol(protocol, EC, EJ, T1, T2, flux, n_levels,
                 wf_type, wf_amp, wf_freq, wf_center, wf_width,
                 wf_dur, wf_npts, tau_max, n_tau):
    try:
        q = _make_qubit(EC, EJ, T1, T2, flux, n_levels)
        t = np.linspace(0, wf_dur, int(wf_npts))
        tau_list = np.linspace(0, tau_max, int(n_tau))
        phi = FluxSignal(type=WF_TYPE_MAP.get(wf_type, 3), t_list=t,
                         amplitude=wf_amp, frequency=wf_freq,
                         center=wf_center, width=wf_width)

        # Two-panel layout: input + result
        fig, axes = plt.subplots(1, 2, figsize=(13, 4))

        # Panel 1 — Input flux signal
        ax = axes[0]
        ax.fill_between(phi.t_list, phi.samples, alpha=0.18, color=COLORS["control"])
        ax.plot(phi.t_list, phi.samples, color=COLORS["control"], linewidth=1.6,
                label="Input Φ(t)")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Φ (Φ₀)")
        ax.set_title("Input Flux Signal")
        ax.legend(loc="best")

        # Panel 2 — Protocol result
        ax = axes[1]

        if protocol == "Rabi":
            t_rabi = np.linspace(0, tau_max, int(n_tau))
            exp = RabiExperiment(qubit=q, t_rabi=t_rabi)
            with suppress_stdout():
                result = exp.run()
            # RabiExperiment returns raw qutip.Result
            p_e = np.asarray(result.expect[0])
            ax.plot(t_rabi, p_e, color=COLORS["primary"], linewidth=1.5,
                    marker="o", markersize=3, markevery=max(1, len(t_rabi)//50))
            ax.set_xlabel("Pulse duration (ns)")
            ax.set_ylabel("p_e")
            ax.set_title("Rabi Oscillations")
            ax.set_ylim(-0.05, 1.05)

        elif protocol == "Ramsey":
            exp = RamseyExperiment(qubit=q, flux_signal=phi, tau_list=tau_list)
            with suppress_stdout():
                result = exp.run()
            ax.plot(result.axes["tau"], result.data["p_e"],
                    color=COLORS["primary"], linewidth=1.4)
            ax.set_xlabel("Free precession τ (ns)")
            ax.set_ylabel("p_e")
            ax.set_title("Ramsey Interference")
            ax.set_ylim(-0.05, 1.05)

        elif protocol == "Differential Echo":
            exp = DiffEchoExperiment(qubit=q, flux_signal=phi,
                                     tau_list=tau_list, k=5)
            with suppress_stdout():
                result = exp.run()
            ax.plot(result.axes["tau"], result.data["p_e"],
                    color=COLORS["secondary"], linewidth=1.4)
            ax.set_xlabel("τ (ns)")
            ax.set_ylabel("p_e")
            ax.set_title("Differential Echo (k=5 cycles)")
            ax.set_ylim(-0.05, 1.05)

        elif protocol == "Transient Sensing":
            exp = TransientSensingExperiment(qubit=q, flux_signal=phi)
            with suppress_stdout():
                result = exp.run()
            scan = result.axes.get("scan", np.arange(len(result.data.get("p_e", []))))
            ax.plot(scan, result.data.get("delta_p", result.data.get("p_e", [])),
                    color=COLORS["accent"], linewidth=1.4)
            ax.set_xlabel("Scan position (ns)")
            ax.set_ylabel("Δp")
            ax.set_title("Transient Sensing — Δp(t)")

        elif protocol == "Cryoscope":
            exp = CryoscopeExperiment(qubit=q, flux_signal=phi)
            with suppress_stdout():
                result = exp.run()
            ax.plot(result.axes["trunc"], result.data["varphi"],
                    color=COLORS["primary"], linewidth=0,
                    marker="o", markersize=3.5)
            ax.set_xlabel("Truncation time (ns)")
            ax.set_ylabel("φ (rad)")
            ax.set_title("Cryoscope Phase  —  φ vs t_d")

        else:
            return _empty_plot(f"Unknown protocol: {protocol}")

        fig.tight_layout()
        return fig

    except Exception as e:
        return _empty_plot(f"Error:\n{traceback.format_exc()}")


# ═══════════════════════════════════════════════════════════════════════════
# TAB 4 — Reconstruction
# ═══════════════════════════════════════════════════════════════════════════

def run_reconstruction(EC, EJ, T1, T2, flux, n_levels, recon_method, lambda_reg):
    try:
        q = _make_qubit(EC, EJ, T1, T2, flux, n_levels)
        t = np.linspace(0, 100, 200)
        phi_true = FluxSignal(type=3, t_list=t, amplitude=0.01,
                              center=50, width=5)
        exp = TransientSensingExperiment(qubit=q, flux_signal=phi_true)
        with suppress_stdout():
            result = exp.run()

        if recon_method == "Wiener Deconvolution":
            recon = TransientReconstruction(method="wiener", lambda_reg=lambda_reg)
            dt = result.axes["scan"][1] - result.axes["scan"][0]
            phi_rec = recon.reconstruct(
                measurement=ExperimentResult(
                    data={"delta_p": result.data["delta_p"]},
                    axes={"scan": result.axes["scan"]},
                ),
                kernel=result.data["kernel"],
                dt=dt,
            )
            fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))

            # 1. Measured Δp
            ax = axes[0]
            ax.fill_between(result.axes["scan"], result.data["delta_p"],
                            color=COLORS["accent"], alpha=0.18)
            ax.plot(result.axes["scan"], result.data["delta_p"],
                    color=COLORS["accent"], linewidth=1.4)
            ax.set_xlabel("Scan position (ns)")
            ax.set_ylabel("Δp")
            ax.set_title("1. Measurement (Δp)")

            # 2. Kernel
            ax = axes[1]
            ax.plot(result.axes["t_samples"], result.data["kernel"],
                    color=COLORS["warning"], linewidth=1.4)
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("Kernel")
            ax.set_title("2. Control Kernel")

            # 3. True vs Reconstructed
            ax = axes[2]
            ax.plot(phi_true.t_list, phi_true.samples,
                    color=COLORS["primary"], linewidth=2.2, label="True Φ(t)")
            ax.plot(phi_rec.t_list[:len(phi_rec.samples)], phi_rec.samples,
                    color=COLORS["secondary"], linewidth=1.6, linestyle="--",
                    label="Reconstructed Φ̂(t)")
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("Φ (Φ₀)")
            ax.set_title("3. Wiener Reconstruction")
            ax.legend(loc="best")

        elif recon_method == "Kernel Estimation Only":
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(result.axes["t_samples"], result.data["kernel"],
                    color=COLORS["warning"], linewidth=1.6)
            ax.set_xlabel("Time (ns)"); ax.set_ylabel("Kernel value")
            ax.set_title("Estimated Control Kernel")

        else:
            return _empty_plot(f"Method '{recon_method}' not available")

        fig.tight_layout()
        return fig

    except Exception as e:
        return _empty_plot(f"Error:\n{traceback.format_exc()}")


# ═══════════════════════════════════════════════════════════════════════════
# TAB 5 — Distortion & Predistortion
# ═══════════════════════════════════════════════════════════════════════════

def run_distortion_demo(dist_amp, dist_tau, pulse_amp, pulse_width, pulse_center):
    try:
        t = np.linspace(0, 200, 2000)
        target = Waveform(
            t_list=t,
            samples=pulse_amp * np.exp(-0.5 * ((t - pulse_center) / pulse_width) ** 2),
        )
        dist = SingleExponentialDistortion(amplitude=dist_amp, tau=dist_tau)
        distorted = dist.apply_to_waveform(target)

        from sqc.calibration.waveform import PredistortionDesigner
        designer = PredistortionDesigner(method="fir_inverse", n_taps=64,
                                          regularization=1e-4)
        dt_step = float(t[1] - t[0])
        inverse_model = designer.design(dist, dt=dt_step)
        predistorted = inverse_model.apply_to_waveform(target)
        corrected = dist.apply_to_waveform(predistorted)

        fig, axes = plt.subplots(2, 2, figsize=(13, 7))

        # Step response
        ax = axes[0, 0]
        s = dist.step_response(t - 10)
        ax.fill_between(t, s, alpha=0.16, color=COLORS["primary"])
        ax.plot(t, s, color=COLORS["primary"], linewidth=1.8)
        ax.axhline(1.0, color=COLORS["neutral"], linestyle=":", linewidth=1, alpha=0.7)
        ax.set_title("Step Response  s(t)")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Output")
        ax.set_ylim(min(s.min()*0.95, 0.85), 1.02)

        # Frequency response
        ax = axes[0, 1]
        omega = np.linspace(1e-4, 0.5, 500)
        H = dist.frequency_response(omega)
        ax.semilogy(omega, np.abs(H), color=COLORS["secondary"], linewidth=1.8)
        ax.set_title("Magnitude Response  |H(ω)|")
        ax.set_xlabel("ω  (rad/ns)"); ax.set_ylabel("|H(ω)|")

        # Without predistortion
        ax = axes[1, 0]
        ax.plot(t, target.samples, color=COLORS["primary"], linewidth=2.2,
                label="Target", alpha=0.85)
        ax.plot(t, distorted.samples, color=COLORS["secondary"], linewidth=1.4,
                linestyle="--", label="On-chip (distorted)")
        rmse_un = np.sqrt(np.mean((distorted.samples - target.samples)**2))
        ax.set_title(f"Without Predistortion  —  RMSE = {rmse_un:.3e}")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Amplitude")
        ax.legend(loc="best")

        # With predistortion
        ax = axes[1, 1]
        ax.plot(t, target.samples, color=COLORS["primary"], linewidth=2.2,
                label="Target", alpha=0.85)
        ax.plot(t, corrected.samples, color=COLORS["accent"], linewidth=1.4,
                linestyle="--", label="On-chip (corrected)")
        rmse_co = np.sqrt(np.mean((corrected.samples - target.samples)**2))
        improvement = rmse_un / max(rmse_co, 1e-30)
        ax.set_title(f"With Predistortion  —  RMSE = {rmse_co:.3e}  ({improvement:.1f}× better)")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Amplitude")
        ax.legend(loc="best")

        fig.tight_layout()
        return fig

    except Exception as e:
        return _empty_plot(f"Error:\n{traceback.format_exc()}")


# ═══════════════════════════════════════════════════════════════════════════
# TAB 6 — Z-Crosstalk (multi-qubit)
# ═══════════════════════════════════════════════════════════════════════════

def run_zcrosstalk_demo(h_aa, h_ab, h_ba, h_bb, pulse_amp, pulse_width):
    try:
        t = np.linspace(0, 200, 2000)
        v_a = Waveform(
            t_list=t,
            samples=pulse_amp * np.exp(-0.5 * ((t - 80) / pulse_width) ** 2),
        )
        v_b = Waveform(t_list=t, samples=np.zeros_like(t))

        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[h_aa, h_ab], [h_ba, h_bb]]),
            source_names=["QA", "QB"],
            target_names=["QA", "QB"],
        )
        fluxes = tm.apply({"QA": v_a, "QB": v_b})
        phi_a = fluxes["QA"].samples
        phi_b = fluxes["QB"].samples

        crosstalk_db = 20*np.log10(max(abs(phi_b).max(), 1e-12)
                                    / max(abs(phi_a).max(), 1e-12))

        fig, axes = plt.subplots(2, 2, figsize=(13, 7))

        # Source pulse on QA
        ax = axes[0, 0]
        ax.fill_between(t, v_a.samples, alpha=0.18, color=COLORS["primary"])
        ax.plot(t, v_a.samples, color=COLORS["primary"], linewidth=1.6,
                label="V_A applied")
        ax.set_title("Source: AWG pulse on QA Z-line")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("V_A (a.u.)")
        ax.legend(loc="best")

        # Source on QB (zero)
        ax = axes[0, 1]
        ax.plot(t, v_b.samples, color=COLORS["neutral"], linewidth=1.6)
        ax.set_title("Source: QB Z-line (no input)")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("V_B (a.u.)")

        # On-chip flux on QA
        ax = axes[1, 0]
        ax.fill_between(t, phi_a, alpha=0.18, color=COLORS["accent"])
        ax.plot(t, phi_a, color=COLORS["accent"], linewidth=1.6,
                label=f"Φ_A = H_AA·V_A (peak={abs(phi_a).max():.4f})")
        ax.set_title("On-chip Φ_A  (intended channel)")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Φ_A (Φ₀)")
        ax.legend(loc="best")

        # Crosstalk on QB
        ax = axes[1, 1]
        ax.fill_between(t, phi_b, alpha=0.18, color=COLORS["secondary"])
        ax.plot(t, phi_b, color=COLORS["secondary"], linewidth=1.6,
                label=f"Φ_B = H_BA·V_A (peak={abs(phi_b).max():.4f})")
        ax.set_title(f"Parasitic Φ_B  —  Crosstalk = {crosstalk_db:.1f} dB")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Φ_B (Φ₀)")
        ax.legend(loc="best")

        fig.tight_layout()
        return fig

    except Exception as e:
        return _empty_plot(f"Error:\n{traceback.format_exc()}")


# ═══════════════════════════════════════════════════════════════════════════
# TAB 7 — Stack Inspector (live code references)
# ═══════════════════════════════════════════════════════════════════════════

STACK_REFERENCES = {
    "Run a Ramsey experiment": [
        ("devices",       "from sqc.devices.transmon import TransmonQubit"),
        ("control",       "from sqc.control.flux_signal import FluxSignal"),
        ("control",       "from sqc.control.sequence import create_ramsey_pulse"),
        ("simulation",    "from sqc.simulation.hamiltonian import HamiltonianBuilder"),
        ("simulation",    "from sqc.simulation.runner import MesolveRunner"),
        ("experiments",   "from sqc.experiments.ramsey import RamseyExperiment"),
    ],
    "Wiener-deconvolve a transient signal": [
        ("devices",       "from sqc.devices.transmon import TransmonQubit"),
        ("control",       "from sqc.control.flux_signal import FluxSignal"),
        ("experiments",   "from sqc.experiments.transient import TransientSensingExperiment"),
        ("reconstruction","from sqc.reconstruction.kernel import KernelEstimator"),
        ("reconstruction","from sqc.reconstruction.transient import TransientReconstruction"),
    ],
    "Predistort a flux pulse": [
        ("control",       "from sqc.control.waveform import Waveform"),
        ("hardware",      "from sqc.hardware.distortion import SingleExponentialDistortion"),
        ("hardware",      "from sqc.hardware.control_line import ControlLine"),
        ("calibration",   "from sqc.calibration.waveform import PredistortionDesigner"),
        ("workflows",     "from sqc.workflows.predistortion_validation import PredistortionValidationWorkflow"),
    ],
    "Quantify Z-crosstalk between two qubits": [
        ("devices",       "from sqc.devices.chip import ChipTopology"),
        ("hardware",      "from sqc.hardware.transfer_matrix import TransferMatrix"),
        ("control",       "from sqc.control.waveform import Waveform"),
        ("experiments",   "from sqc.experiments.transient import TransientSensingExperiment"),
        ("workflows",     "from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow"),
    ],
    "LM full-density-matrix waveform inversion": [
        ("devices",       "from sqc.devices.transmon import TransmonQubit"),
        ("control",       "from sqc.control.pulse import CompositePulse"),
        ("reconstruction","from sqc.reconstruction.basis import generate_basis_functions"),
        ("reconstruction","from sqc.reconstruction.transient import TransientReconstruction"),
    ],
}


def show_stack_for_task(task):
    if task not in STACK_REFERENCES:
        return "<p>No reference for this task.</p>"
    refs = STACK_REFERENCES[task]
    layers_used = sorted(set(layer for layer, _ in refs),
                         key=lambda x: ["devices","hardware","control","simulation",
                                        "experiments","reconstruction","calibration","workflows"].index(x))

    html = f'<div class="info-card"><h3>📦 Task: {task}</h3>'
    html += '<p><b>Layers involved:</b> '
    html += "".join(f'<span class="layer-badge badge-{l}">{l}</span>' for l in layers_used)
    html += '</p><h3>📜 Required imports</h3>'
    html += '<pre style="background:#1A202C; color:#F8F9FA; padding:14px; border-radius:8px; overflow-x:auto; font-size:13px;">'
    for layer, imp in refs:
        html += f'<span style="color:#A78BFA;"># {layer}/</span>\n{imp}\n'
    html += '</pre></div>'
    return html


# ═══════════════════════════════════════════════════════════════════════════
# UI Construction
# ═══════════════════════════════════════════════════════════════════════════

# Feature flag for experimental / post-v1 tabs (e.g. Z-Crosstalk). These are
# hidden from the v1 public demo (decision D3): the implementation stays in
# sqc/ and is importable via its deep path, only the interactive tab is gated.
# Flip to True to restore the experimental tabs.
SHOW_EXPERIMENTAL = False


def build_app():
    theme = gr.themes.Soft(
        primary_hue="blue",
        secondary_hue="cyan",
        neutral_hue="slate",
        font=[gr.themes.GoogleFont("Inter"), "system-ui", "sans-serif"],
    )

    with gr.Blocks(theme=theme, css=CUSTOM_CSS,
                   title="Sensing-Project — cQED Full-Stack Demo",
                   analytics_enabled=False) as app:

        # --- Header ---
        gr.HTML("""
<div id="header">
  <h1>🔬 Sensing-Project — Superconducting Quantum Control Platform</h1>
  <p>Interactive exploration of the <b>six-layer cQED stack</b> (Gao 2021):
  Transmon qubit sensing • Waveform reconstruction • Calibration • Predistortion • Multi-qubit crosstalk.<br>
  Powered by the refactored <code>sqc/</code> framework — 230 tests passing • src/ untouched.</p>
</div>
""")

        with gr.Tabs():

            # ─────────────────── Tab: Architecture ───────────────────
            with gr.Tab("🏗️ Architecture"):
                with gr.Tabs():
                    with gr.Tab("Layered View"):
                        gr.Markdown("### The six-layer cQED stack")
                        gr.HTML('''
<div class="info-card">
  <p>Each layer is a Python sub-package under <code>sqc/</code>. <b>Higher layers depend on lower
  layers, never the reverse.</b> This mirrors the real cQED control-electronics stack from
  Gao et al., <i>PRX Quantum</i> 2, 040202 (2021).</p>
</div>
''')
                        gr.Plot(render_architecture_layered(), label="cqed_stack",
                                show_label=False)

                    with gr.Tab("Data Flow"):
                        gr.Markdown("### How a measurement becomes a recovered waveform")
                        gr.HTML('''
<div class="info-card">
  <p>A typical sensing experiment threads through <b>~7 modules</b>: control signals are
  combined with device parameters in <code>HamiltonianBuilder</code>, evolved by <code>MesolveRunner</code>,
  postprocessed by an <code>Experiment</code>, and finally inverted by a <code>Reconstruction</code> class.</p>
</div>
''')
                        gr.Plot(render_data_flow(), label="data_flow", show_label=False)

                    with gr.Tab("Module Index"):
                        gr.HTML('''
<div class="info-card">
  <h3>📚 Module Reference</h3>
  <table style="width:100%; border-collapse: collapse; font-size:13px;">
    <thead style="background:#F1F5F9;">
      <tr>
        <th style="text-align:left; padding:8px;">Layer</th>
        <th style="text-align:left; padding:8px;">Path</th>
        <th style="text-align:left; padding:8px;">Key Classes / Functions</th>
      </tr>
    </thead>
    <tbody>
      <tr><td><span class="layer-badge badge-workflows">workflows</span></td>
          <td><code>sqc/workflows/</code></td>
          <td>PredistortionValidationWorkflow • ZCrosstalkWorkflow</td></tr>
      <tr style="background:#F8FAFC;"><td><span class="layer-badge badge-calibration">calibration</span></td>
          <td><code>sqc/calibration/</code></td>
          <td>FluxResponseCalibration • SinglePointFrequencyCalibration • WaveformCalibration • PredistortionDesigner • CalibrationScheduler</td></tr>
      <tr><td><span class="layer-badge badge-reconstruction">reconstruction</span></td>
          <td><code>sqc/reconstruction/</code></td>
          <td>RamseyReconstruction • EchoReconstruction • TransientReconstruction • CryoscopeReconstruction • DelayRamseyReconstruction • PiPulseCompReconstruction</td></tr>
      <tr style="background:#F8FAFC;"><td><span class="layer-badge badge-experiments">experiments</span></td>
          <td><code>sqc/experiments/</code></td>
          <td>RabiExperiment • RamseyExperiment • DiffEchoExperiment • TransientSensingExperiment • CryoscopeExperiment</td></tr>
      <tr><td><span class="layer-badge badge-simulation">simulation</span></td>
          <td><code>sqc/simulation/</code></td>
          <td>HamiltonianBuilder • MesolveRunner • SlidingMeasurementRunner • ExperimentResult • noise</td></tr>
      <tr style="background:#F8FAFC;"><td><span class="layer-badge badge-control">control</span></td>
          <td><code>sqc/control/</code></td>
          <td>Waveform • FluxSignal • Pulse • CompositePulse • create_*_pulse • gates</td></tr>
      <tr><td><span class="layer-badge badge-hardware">hardware</span></td>
          <td><code>sqc/hardware/</code></td>
          <td>ControlLine • SingleExponentialDistortion • TransferMatrix • IQReadoutModel</td></tr>
      <tr style="background:#F8FAFC;"><td><span class="layer-badge badge-devices">devices</span></td>
          <td><code>sqc/devices/</code></td>
          <td>QubitSpec • TransmonQubit • Resonator • Coupler • ChipTopology</td></tr>
    </tbody>
  </table>
</div>
''')

            # ─────────────────── Tab: Qubit ───────────────────
            with gr.Tab("⚛️ Qubit Configuration"):
                gr.HTML(_layer_badges("devices"))
                with gr.Row():
                    with gr.Column(scale=1):
                        gr.Markdown("### 🔧 Transmon parameters")
                        ec = gr.Slider(0.10, 0.50, 0.20, step=0.01,
                                       label="EC  (GHz / 2π)",
                                       info="Charging energy — sets anharmonicity")
                        ej = gr.Slider(5.0, 25.0, 15.0, step=0.5,
                                       label="EJ  (GHz / 2π)",
                                       info="Josephson energy at zero flux")
                        flux_bias = gr.Slider(-0.5, 0.5, 0.0, step=0.005,
                                              label="Flux bias  (Φ₀)",
                                              info="Static flux through SQUID loop")
                        nlev = gr.Radio([2, 3, 4, 5], value=3,
                                        label="n_levels  (Fock truncation)")
                        with gr.Accordion("⏱ Coherence times", open=False):
                            t1 = gr.Slider(1000, 50000, 10000, step=1000, label="T₁  (ns)")
                            t2 = gr.Slider(1000, 50000, 8000, step=1000, label="T₂  (ns)")
                    with gr.Column(scale=2):
                        qubit_html = gr.HTML(qubit_info_html(0.2, 15, 10000, 8000, 0.0, 3))

                gr.Markdown("### 📊 Frequency-flux dispersion")
                spectrum_plot = gr.Plot(qubit_spectrum_plot(0.2, 15, 10000, 8000, 0.0, 3),
                                         show_label=False)

                # Live updates on slider change
                inputs = [ec, ej, t1, t2, flux_bias, nlev]
                for slider in inputs:
                    slider.change(qubit_info_html, inputs, qubit_html, show_progress="hidden")
                    slider.change(qubit_spectrum_plot, inputs, spectrum_plot, show_progress="hidden")

            # ─────────────────── Tab: Waveform ───────────────────
            with gr.Tab("📈 Waveform Designer"):
                gr.HTML(_layer_badges("control"))
                gr.Markdown("Design a flux signal **Φ(t)** with one of 7 parametric shapes; "
                            "see both the time waveform and its FFT spectrum live.")
                with gr.Row():
                    with gr.Column(scale=1):
                        wf_type = gr.Dropdown(
                            list(WF_TYPE_MAP.keys()),
                            value="Gaussian", label="Signal type",
                        )
                        wf_amp = gr.Slider(0.0, 0.1, 0.01, step=0.001,
                                           label="Amplitude  (Φ₀)")
                        wf_freq = gr.Slider(0.0, 0.1, 0.01, step=0.001,
                                            label="Frequency  (1/ns) — for sinusoid")
                        wf_center = gr.Slider(0, 250, 100, step=1, label="Center  (ns)")
                        wf_width = gr.Slider(1, 50, 10, step=1, label="Width  (ns)")
                        with gr.Accordion("⚙ Advanced", open=False):
                            wf_rise = gr.Slider(1, 50, 10, step=1, label="Rise time  (ns)")
                            wf_fall = gr.Slider(1, 50, 10, step=1, label="Fall time  (ns)")
                            wf_dur = gr.Slider(50, 500, 250, step=10, label="Duration  (ns)")
                            wf_npts = gr.Slider(100, 2000, 500, step=50, label="N points")
                            wf_noise = gr.Slider(0.0, 0.01, 0.0, step=0.0001,
                                                  label="Noise level")
                    with gr.Column(scale=2):
                        wf_plot = gr.Plot(
                            plot_waveform("Gaussian", 0.01, 0.01, 100, 10, 10, 10, 250, 500, 0.0),
                            show_label=False,
                        )

                wf_inputs = [wf_type, wf_amp, wf_freq, wf_center, wf_width,
                              wf_rise, wf_fall, wf_dur, wf_npts, wf_noise]
                for ctrl in wf_inputs:
                    ctrl.change(plot_waveform, wf_inputs, wf_plot, show_progress="hidden")

            # ─────────────────── Tab: Protocols ───────────────────
            with gr.Tab("🧪 Protocol Lab"):
                gr.HTML(_layer_badges("devices", "control", "simulation", "experiments"))
                gr.Markdown("Combine the qubit + waveform from Tabs ⚛️/📈 and run a sensing protocol.")

                with gr.Row():
                    with gr.Column(scale=1):
                        # "Differential Echo" is hidden from the v1 demo
                        # (SHOW_EXPERIMENTAL=False): its closed-form retrieval
                        # B = -phi/(2k·kappa·t_int) is not yet reliable
                        # (t_int mis-derived from gap timing, off-scale B).
                        # See RELEASE_TODO.md. Implementation stays importable
                        # via sqc.experiments.echo.DiffEchoExperiment.
                        _proto_choices = [
                            p for p in PROTOCOL_INFO.keys()
                            if SHOW_EXPERIMENTAL or p != "Differential Echo"
                        ]
                        protocol = gr.Dropdown(
                            _proto_choices,
                            value="Ramsey", label="Protocol",
                        )
                        proto_info = gr.HTML(protocol_info_html("Ramsey"))
                        tau_max = gr.Slider(20, 500, 250, step=10,
                                            label="τ max / scan range  (ns)")
                        n_tau = gr.Slider(20, 1000, 200, step=20,
                                          label="N points  (τ/scan)")
                        run_btn = gr.Button("▶ Run Protocol", variant="primary",
                                            elem_classes="primary-btn")
                    with gr.Column(scale=2):
                        proto_plot = gr.Plot(show_label=False)

                protocol.change(protocol_info_html, protocol, proto_info, show_progress="hidden")

                run_btn.click(
                    run_protocol,
                    [protocol, ec, ej, t1, t2, flux_bias, nlev,
                     wf_type, wf_amp, wf_freq, wf_center, wf_width, wf_dur, wf_npts,
                     tau_max, n_tau],
                    proto_plot,
                )

            # ─────────────────── Tab: Reconstruction ───────────────────
            with gr.Tab("🔄 Reconstruction Bench"):
                gr.HTML(_layer_badges("experiments", "reconstruction"))
                gr.Markdown("Recover an unknown flux waveform from measurement data using "
                            "the **transient sensing → kernel → Wiener** pipeline.")
                with gr.Row():
                    with gr.Column(scale=1):
                        recon_method = gr.Dropdown(
                            ["Wiener Deconvolution", "Kernel Estimation Only"],
                            value="Wiener Deconvolution", label="Reconstruction method",
                        )
                        lambda_reg = gr.Slider(0.001, 100.0, 1.0, step=0.1,
                                               label="Regularization  λ")
                        gr.Markdown(
                            "*True signal: Gaussian @ t=50 ns, σ=5 ns, A=0.01 Φ₀.*\n\n"
                            "*Increase λ for noisy measurements; decrease for sharp signals.*"
                        )
                        recon_btn = gr.Button("▶ Run Reconstruction", variant="primary",
                                              elem_classes="primary-btn")
                    with gr.Column(scale=2):
                        recon_plot = gr.Plot(show_label=False)

                recon_btn.click(
                    run_reconstruction,
                    [ec, ej, t1, t2, flux_bias, nlev, recon_method, lambda_reg],
                    recon_plot,
                )

            # ─────────────────── Tab: Distortion ───────────────────
            with gr.Tab("🔧 Distortion & Predistortion"):
                gr.HTML(_layer_badges("hardware", "calibration"))
                gr.Markdown("Inject a single-exponential flux-line distortion and validate "
                            "**FIR predistortion** that cancels it.")
                with gr.Row():
                    with gr.Column(scale=1):
                        gr.Markdown("### Distortion model")
                        dist_amp = gr.Slider(0.0, 0.2, 0.05, step=0.005,
                                              label="Tail amplitude")
                        dist_tau = gr.Slider(5, 200, 50, step=5,
                                              label="Tail τ  (ns)")
                        gr.Markdown("### Target pulse")
                        pulse_amp2 = gr.Slider(0.0, 2.0, 1.0, step=0.1,
                                                label="Pulse amplitude")
                        pulse_w = gr.Slider(1, 30, 8, step=1, label="Pulse width  (ns)")
                        pulse_c = gr.Slider(0, 200, 80, step=5, label="Pulse center  (ns)")
                        dist_btn = gr.Button("▶ Run Distortion + Predistortion",
                                             variant="primary", elem_classes="primary-btn")
                    with gr.Column(scale=2):
                        dist_plot = gr.Plot(show_label=False)

                dist_btn.click(
                    run_distortion_demo,
                    [dist_amp, dist_tau, pulse_amp2, pulse_w, pulse_c],
                    dist_plot,
                )

            # ─────────── Tab: Z-Crosstalk (experimental — gated, D3) ───────────
            # Hidden from v1 (SHOW_EXPERIMENTAL=False). Flip the flag to restore.
            if SHOW_EXPERIMENTAL:
                with gr.Tab("🌐 Z-Crosstalk (2-qubit)"):
                    gr.HTML(_layer_badges("hardware", "devices"))
                    gr.Markdown("Apply a flux pulse on **QA**'s Z-line; observe the parasitic flux "
                                "on **QB** through the **TransferMatrix** model.")
                    with gr.Row():
                        with gr.Column(scale=1):
                            gr.Markdown("### Transfer Matrix  H[target, source]")
                            with gr.Row():
                                h_aa = gr.Number(value=1.00, label="H_AA (self)", precision=3)
                                h_ab = gr.Number(value=0.00, label="H_AB (B → A)", precision=3)
                            with gr.Row():
                                h_ba = gr.Number(value=0.05, label="H_BA  (A → B)  ⚡", precision=3)
                                h_bb = gr.Number(value=1.00, label="H_BB (self)", precision=3)
                            gr.Markdown("### Pulse on QA")
                            zc_amp = gr.Slider(0.0, 2.0, 1.0, step=0.1, label="Pulse amplitude")
                            zc_width = gr.Slider(1, 30, 8, step=1, label="Pulse width  (ns)")
                            zc_btn = gr.Button("▶ Run Z-Crosstalk", variant="primary",
                                                elem_classes="primary-btn")
                        with gr.Column(scale=2):
                            zc_plot = gr.Plot(show_label=False)

                    zc_btn.click(
                        run_zcrosstalk_demo,
                        [h_aa, h_ab, h_ba, h_bb, zc_amp, zc_width],
                        zc_plot,
                    )

            # ─────────────────── Tab: Stack Inspector ───────────────────
            with gr.Tab("🔍 Stack Inspector"):
                gr.Markdown("### See exactly which `sqc/` layers your task uses")
                with gr.Row():
                    with gr.Column(scale=1):
                        task = gr.Dropdown(
                            list(STACK_REFERENCES.keys()),
                            value=list(STACK_REFERENCES.keys())[0],
                            label="Choose a task",
                        )
                    with gr.Column(scale=2):
                        stack_html = gr.HTML(show_stack_for_task(list(STACK_REFERENCES.keys())[0]))

                task.change(show_stack_for_task, task, stack_html, show_progress="hidden")

            # ─────────────────── Tab: About ───────────────────
            with gr.Tab("📖 About"):
                gr.HTML(f"""
<div class="info-card">
  <h3>Sensing-Project v2 — sqc v0.1.0</h3>
  <p><b>What it is:</b> a research-grade simulation framework for flux-tunable Transmon
  qubit sensing. The codebase has been refactored into a clean six-layer architecture
  (per Gao 2021) with full backward compatibility (<code>src_mirror/</code>) and a
  comprehensive regression-test suite.</p>

  <h3>Stats</h3>
  <ul>
    <li><b>~50 modules</b> across 8 sub-packages in <code>sqc/</code></li>
    <li><b>230 tests passing</b> (185 unit + 14 equivalence + 24 integration + 7 regression)</li>
    <li><b>0 modifications to <code>src/</code></b> — old code untouched (R1 compliance)</li>
    <li><b>7 frozen physics baselines</b> (rtol = 1e-6, atol = 1e-9)</li>
  </ul>

  <h3>Documentation</h3>
  <ul>
    <li>Technical docs: <code>docs/architecture.md</code></li>
    <li>Refactoring plan: <code>idea/refactor/_refactor_plan.md</code></li>
    <li>Phase handbooks: <code>idea/refactor/phase_*_handbook.md</code></li>
    <li>Session log: <code>idea/refactor/_session_log_2026-05-01.md</code></li>
  </ul>

  <h3>References</h3>
  <ul>
    <li>Gao, Rol, Touzard, Wang, "Practical Guide for Building Superconducting Quantum Devices",
        <i>PRX Quantum</i> <b>2</b>, 040202 (2021)</li>
    <li>Koch et al., <i>Phys. Rev. A</i> <b>76</b>, 042319 (2007) — Transmon qubit</li>
    <li>Motzoi et al., <i>Phys. Rev. Lett.</i> <b>103</b>, 110501 (2009) — DRAG pulses</li>
  </ul>
</div>
""")

        # Footer
        gr.HTML('<div id="footer">Sensing-Project · sqc v0.1.0 · Built with Gradio · '
                'Conda env: <code>qutip-env</code></div>')

    return app


# ═══════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 64)
    print("Sensing-Project — cQED Full-Stack Demo  (sqc v0.1.0)")
    print("=" * 64)
    app = build_app()
    app.launch(
        server_name="127.0.0.1",
        server_port=7860,
        share=False,
        show_error=True,
        inbrowser=True,
    )
