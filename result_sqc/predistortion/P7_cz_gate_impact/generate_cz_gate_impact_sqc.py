#!/usr/bin/env python3
"""P7 — Impact of predistortion on CZ gate performance (sqc).

Takes the cryoscope-derived predistortion filter (from P6) and applies it
to a CZ flux pulse.  Three groups are compared:
  - Ideal baseline: target flux reaches qubit directly (no control line)
  - Uncorrected: target AWG → control_line distortion → chip flux → gate
  - Corrected: predistorted AWG → control_line → chip flux → gate

All groups use the SAME gate parameters, Hamiltonian, basis, and solver
settings.  The predistortion filter is designed from the cryoscope step
response and frozen before any gate metric is computed.

Produces Figure E2 and associated data files.

Cost: HIGH — cryoscope measurement (~3s) + 3× CZ gate mesolve (~10s each)
       + robustness scan.  Cached to npz; pass --recompute to force rebuild.

Usage:
    python result_sqc/predistortion/P7_cz_gate_impact/generate_cz_gate_impact_sqc.py
    python result_sqc/predistortion/P7_cz_gate_impact/generate_cz_gate_impact_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os
import io
import contextlib
import json
import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "predistortion/P7_cz_gate_impact"
_CACHE = "cz_gate_impact_sqc"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Qubit parameters (same as P6 / _common.py)
EC = C.EC          # rad·GHz
EJ = C.EJ          # rad·GHz
T1 = C.T1          # ns
T2 = C.T2          # ns
N_LEVELS = 3

# Coupling
G_COUPLING = 0.05   # GHz

# CZ pulse calibrated once on the ideal control line, then frozen for every
# comparison group.  Q2 is statically detuned from the flux-pulsed Q1.
Q2_FLUX_BIAS = 0.2166779812  # Φ₀
CZ_FLUX_AMP = 0.1691783618   # Φ₀
CZ_DURATION = 147.3225170    # ns
CZ_RISE_FRAC = 0.2761768316

# Control-line ground truth (same as P6: dual-exponential, model mismatch)
TRUE_AMPS = [0.04, 0.02]
TRUE_TAUS = [80.0, 400.0]

# Cryoscope measurement (for predistortion)
T_MAX_STEP = 200.0      # ns
STEP_AMP = 0.05         # Φ₀
FIT_TYPE = "single_exp"

# Robustness scan: vary the dominant distortion amplitude
SCAN_AMPS = np.linspace(0.032, 0.048, 9)
SCAN_TAU_FIXED = 80.0   # ns — vary amp of the fast component

# Styles
STYLE_IDEAL = {"color": "k", "ls": "-"}
STYLE_UNCORRECTED = {"color": "#e41a1c", "ls": "-"}
STYLE_CORRECTED = {"color": "#377eb8", "ls": "--"}
LW_MAIN = 1.6
LW_SEC = 1.2


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _quiet(fn, *a, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


def _make_qubit_pair(flux1=0.0, flux2=Q2_FLUX_BIAS):
    from sqc.devices.transmon import TransmonQubit
    q1 = TransmonQubit(EC=EC, EJ=EJ, T1=T1, T2=T2, flux=flux1,
                        n_levels=N_LEVELS)
    q2 = TransmonQubit(EC=EC, EJ=EJ, T1=T1, T2=T2, flux=flux2,
                        n_levels=N_LEVELS)
    return q1, q2


def _make_control_line(amps=None, taus=None):
    from sqc.hardware.control_line import ControlLine
    from sqc.hardware.distortion import MultiExponentialDistortion
    amps = TRUE_AMPS if amps is None else amps
    taus = TRUE_TAUS if taus is None else taus
    return ControlLine(
        name="Z0", kind="z", source="AWG0", target="Q0",
        transfer_function=MultiExponentialDistortion(
            amplitudes=list(amps), taus=list(taus),
        ),
    )


def _design_cz_waveform(flux_amp=CZ_FLUX_AMP):
    """Design a CZ flux pulse with sin² edges.

    Returns (t_axis, phi_waveform) where phi_waveform is in Φ₀.
    """
    from sqc.config import CONFIG

    T = CZ_DURATION
    T1 = CZ_RISE_FRAC * T
    T2 = T1

    t_axis = CONFIG.pulse.make_time(0, T)
    phi = np.zeros(len(t_axis))
    for i, t in enumerate(t_axis):
        if t < 0:
            v = 0.0
        elif t < T1:
            v = flux_amp * np.sin(0.5 * np.pi * t / T1) ** 2
        elif t < T - T2:
            v = flux_amp
        else:
            v = flux_amp * np.sin(0.5 * np.pi * (T - t) / T2) ** 2
        phi[i] = v
    return t_axis, phi


def _measure_cryoscope_step(control_line):
    """Measure step response via cryoscope → fit → design predistortion.

    Returns (fitted_model, inverse_model).
    """
    from sqc.calibration.waveform import WaveformCalibration, PredistortionDesigner
    from sqc.hardware.distortion import SingleExponentialDistortion

    q = _make_qubit_pair(flux1=0.0)[0]  # cryoscope needs sweet spot

    cal = WaveformCalibration(
        qubit=q,
        control_line=control_line,
        measurement_protocol="cryoscope",
        method="transfer_function",
        fit_type=FIT_TYPE,
        t_max=T_MAX_STEP,
        step_amplitude=STEP_AMP,
    )
    table = _quiet(cal.calibrate)
    fp = table.fit_params
    fitted = SingleExponentialDistortion(
        amplitude=float(fp.get("amplitude", 0.01)),
        tau=float(fp.get("tau", 100.0)),
    )
    from sqc.config import CONFIG as _CFG
    dt = float(_CFG.awg.dt)
    designer = PredistortionDesigner(method="auto", regularization=1e-4)
    inverse = designer.design(fitted, dt=dt)
    return fitted, inverse


def _flux_to_waveform(t_axis, phi):
    """Convert flux array to Waveform object for ControlLine.apply()."""
    from sqc.control.waveform import Waveform
    return Waveform(t_list=t_axis.copy(), samples=np.asarray(phi, dtype=float))


def _apply_control_line(t_axis, phi_awg, control_line):
    """Pass AWG waveform through control line → on-chip flux."""
    wf_in = _flux_to_waveform(t_axis, phi_awg)
    wf_out = control_line.apply(wf_in)
    return np.asarray(wf_out.samples, dtype=float)


# ---------------------------------------------------------------------------
# Gate simulation wrapper
# ---------------------------------------------------------------------------

def _simulate_cz(t_axis, phi_chip):
    """Simulate CZ gate with the given on-chip flux waveform.

    Returns (CZResult, q1, q2).
    """
    from sqc.control.gates import simulate_cz_from_flux
    q1, q2 = _make_qubit_pair(flux1=0.0, flux2=Q2_FLUX_BIAS)
    result = simulate_cz_from_flux(
        t_axis, phi_chip, q1, q2, G_COUPLING,
        solver_options={"nsteps": 100000, "max_step": 0.25},
    )
    return result, q1, q2


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------

def compute(recompute=False):
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    out: dict = {}

    # -- 1. Design CZ target waveform --
    t_cz, phi_target = _design_cz_waveform()
    out.update({
        "t_cz": t_cz,
        "phi_target": phi_target,
        "cz_flux_amp": np.array([CZ_FLUX_AMP]),
        "g_coupling": np.array([G_COUPLING]),
        "q2_flux_bias": np.array([Q2_FLUX_BIAS]),
        "cz_duration": np.array([CZ_DURATION]),
        "cz_rise_frac": np.array([CZ_RISE_FRAC]),
    })

    # -- 2. Control line --
    cl = _make_control_line()

    # -- 3. Cryoscope step response → predistortion --
    fitted_model, inverse_model = _measure_cryoscope_step(cl)
    dt = float(t_cz[1] - t_cz[0])
    out.update({
        "fitted_amplitude": np.array([fitted_model.amplitude]),
        "fitted_tau": np.array([fitted_model.tau]),
    })

    # -- 4. Three-group comparison --
    # (a) Ideal: target flux → qubit directly
    result_ideal, q1, q2 = _simulate_cz(t_cz, phi_target)
    out.update(_result_to_dict("ideal", result_ideal))

    # (b) Uncorrected: target AWG → control_line → chip flux
    phi_uncorrected_chip = _apply_control_line(t_cz, phi_target, cl)
    result_unc, _, _ = _simulate_cz(t_cz, phi_uncorrected_chip)
    out.update({
        "phi_chip_uncorrected": phi_uncorrected_chip,
        **_result_to_dict("uncorrected", result_unc),
    })

    # (c) Corrected: predistorted AWG → control_line → chip flux
    phi_predistorted_awg = np.asarray(
        inverse_model.apply_to_waveform(_flux_to_waveform(t_cz, phi_target)).samples,
        dtype=float,
    )
    phi_corrected_chip = _apply_control_line(t_cz, phi_predistorted_awg, cl)
    result_cor, _, _ = _simulate_cz(t_cz, phi_corrected_chip)
    out.update({
        "awg_predistorted": phi_predistorted_awg,
        "phi_chip_corrected": phi_corrected_chip,
        **_result_to_dict("corrected", result_cor),
    })

    # -- 5. Robustness scan: vary distortion amplitude --
    scan = {key: [] for key in [
        "scan_amps", "phase_err_ideal", "phase_err_unc", "phase_err_cor",
        "leak_ideal", "leak_unc", "leak_cor",
        "infid_ideal", "infid_unc", "infid_cor",
    ]}
    for amp in SCAN_AMPS:
        cl_scan = _make_control_line(amps=[amp, TRUE_AMPS[1]], taus=[SCAN_TAU_FIXED, TRUE_TAUS[1]])
        # Robustness means the calibration is frozen while the line drifts.
        inv_scan = inverse_model

        # Ideal (same target, independent of distortion)
        r_ideal, _, _ = _simulate_cz(t_cz, phi_target)

        # Uncorrected
        phi_unc_chip = _apply_control_line(t_cz, phi_target, cl_scan)
        r_unc, _, _ = _simulate_cz(t_cz, phi_unc_chip)

        # Corrected
        phi_pred_awg = np.asarray(
            inv_scan.apply_to_waveform(_flux_to_waveform(t_cz, phi_target)).samples,
            dtype=float,
        )
        phi_cor_chip = _apply_control_line(t_cz, phi_pred_awg, cl_scan)
        r_cor, _, _ = _simulate_cz(t_cz, phi_cor_chip)

        scan["scan_amps"].append(amp)
        scan["phase_err_ideal"].append(r_ideal.phase_error)
        scan["phase_err_unc"].append(r_unc.phase_error)
        scan["phase_err_cor"].append(r_cor.phase_error)
        scan["leak_ideal"].append(r_ideal.leakage)
        scan["leak_unc"].append(r_unc.leakage)
        scan["leak_cor"].append(r_cor.leakage)
        scan["infid_ideal"].append(r_ideal.infidelity)
        scan["infid_unc"].append(r_unc.infidelity)
        scan["infid_cor"].append(r_cor.infidelity)

    for k, v in scan.items():
        out[k] = np.array(v, dtype=float)

    C.save_npz(SUBDIR, _CACHE, **out)
    return out


def _result_to_dict(prefix, result):
    return {
        f"cond_phase_{prefix}": np.array([result.conditional_phase]),
        f"phase_error_{prefix}": np.array([result.phase_error]),
        f"leakage_{prefix}": np.array([result.leakage]),
        f"infidelity_{prefix}": np.array([result.infidelity]),
        f"pop_20_{prefix}": np.array([result.populations["pop_20"]]),
        f"pop_02_{prefix}": np.array([result.populations["pop_02"]]),
    }


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def _compute_metrics(d):
    m = {}
    for prefix in ["ideal", "uncorrected", "corrected"]:
        m[f"cond_phase_{prefix}"] = float(d[f"cond_phase_{prefix}"][0])
        m[f"phase_error_{prefix}"] = float(d[f"phase_error_{prefix}"][0])
        m[f"leakage_{prefix}"] = float(d[f"leakage_{prefix}"][0])
        m[f"infidelity_{prefix}"] = float(d[f"infidelity_{prefix}"][0])
    # Additional error from distortion
    m["phase_error_added_unc"] = m["phase_error_uncorrected"] - m["phase_error_ideal"]
    m["phase_error_added_cor"] = m["phase_error_corrected"] - m["phase_error_ideal"]
    m["infidelity_added_unc"] = m["infidelity_uncorrected"] - m["infidelity_ideal"]
    m["infidelity_added_cor"] = m["infidelity_corrected"] - m["infidelity_ideal"]
    # Waveform RMSEs
    target = d["phi_target"]
    m["rmse_chip_unc"] = float(np.sqrt(np.mean((d["phi_chip_uncorrected"] - target)**2)))
    m["rmse_chip_cor"] = float(np.sqrt(np.mean((d["phi_chip_corrected"] - target)**2)))
    return m


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(d, out_dir):
    import matplotlib.pyplot as plt

    t_cz = d["t_cz"]
    target = d["phi_target"]
    phi_unc = d["phi_chip_uncorrected"]
    phi_cor = d["phi_chip_corrected"]
    awg_pred = d["awg_predistorted"]

    metrics = _compute_metrics(d)

    fig, ax = plt.subplots(2, 2, figsize=(12.8, 8.8))
    ((a, b), (c, d_ax)) = ax

    # -- (a) On-chip flux trajectories --
    a.plot(t_cz, target, **STYLE_IDEAL, lw=LW_MAIN, label="target (ideal)")
    a.plot(t_cz, phi_unc, **STYLE_UNCORRECTED, lw=LW_SEC,
           label="uncorrected")
    a.plot(t_cz, phi_cor, **STYLE_CORRECTED, lw=LW_MAIN,
           label="corrected")
    a.set_xlabel("t (ns)", fontsize=11)
    a.set_ylabel("On-chip flux $\\Phi_{\\rm chip}$ ($\\Phi_0$)", fontsize=11)
    a.set_title("(a) On-chip flux trajectories", fontsize=11)
    a.legend(fontsize=8.5, loc="upper right")
    a.grid(True, ls=":", alpha=0.4)

    # Inset: AWG predistorted waveform
    ai = a.inset_axes([0.55, 0.12, 0.42, 0.28])
    ai.plot(t_cz, target, **STYLE_IDEAL, lw=0.8, label="target AWG")
    ai.plot(t_cz, awg_pred, **STYLE_CORRECTED, lw=0.8, label="predistorted")
    ai.set_title("AWG waveforms", fontsize=7)
    ai.tick_params(labelsize=6)
    ai.legend(fontsize=6, loc="upper right")
    ai.grid(True, ls=":", alpha=0.3)

    # -- (b) Gate metrics: three compact axes, one physical unit each --
    groups = ["ideal", "uncorrected", "corrected"]
    labels = ["Ideal", "Uncorr.", "Corrected"]
    colors = ["#6b6b6b", "#d73027", "#2878b5"]
    metric_specs = [
        ("phase_error", 1e3, "Phase error", "mrad"),
        ("leakage", 100.0, "Leakage", "%"),
        ("infidelity", 100.0, "Infidelity", "%"),
    ]
    b.set_axis_off()
    b.set_title("(b) Gate metrics", fontsize=11, pad=10)
    for j, (key, scale, title, unit) in enumerate(metric_specs):
        metric_ax = b.inset_axes([0.02 + 0.335 * j, 0.08, 0.30, 0.82])
        values = [metrics[f"{key}_{g}"] * scale for g in groups]
        bars = metric_ax.bar(np.arange(3), values, width=0.64,
                             color=colors, edgecolor="black", lw=0.45)
        metric_ax.set_title(f"{title} ({unit})", fontsize=8.5)
        metric_ax.set_xticks(np.arange(3), labels, rotation=25,
                             ha="right", fontsize=7)
        metric_ax.tick_params(axis="y", labelsize=7)
        metric_ax.grid(True, axis="y", ls=":", alpha=0.35)
        ymax = max(values) if max(values) > 0 else 1.0
        metric_ax.set_ylim(0, 1.24 * ymax)
        for bar, value in zip(bars, values):
            metric_ax.text(bar.get_x() + bar.get_width() / 2,
                           value + 0.035 * ymax, f"{value:.2f}",
                           ha="center", va="bottom", fontsize=6.5)

    # -- (c) Flux residuals --
    res_unc = np.abs(phi_unc - target)
    res_cor = np.abs(phi_cor - target)
    c.semilogy(t_cz, res_unc, **STYLE_UNCORRECTED, lw=LW_SEC,
               label=f"uncorrected (RMSE={metrics['rmse_chip_unc']:.2e})")
    c.semilogy(t_cz, res_cor, **STYLE_CORRECTED, lw=LW_MAIN,
               label=f"corrected (RMSE={metrics['rmse_chip_cor']:.2e})")
    c.set_xlabel("t (ns)", fontsize=11)
    c.set_ylabel("Absolute flux residual ($\\Phi_0$)", fontsize=11)
    c.set_title("(c) On-chip flux residuals", fontsize=11)
    c.legend(fontsize=8.5, loc="lower left")
    c.grid(True, which="both", ls=":", alpha=0.4)

    # -- (d) Robustness to drift around the calibrated line --
    d_ax.plot(d["scan_amps"], 100 * d["infid_ideal"], marker="s", **STYLE_IDEAL,
              lw=LW_MAIN, label="ideal baseline")
    d_ax.plot(d["scan_amps"], 100 * d["infid_unc"], marker="o", **STYLE_UNCORRECTED,
              lw=LW_SEC, label="uncorrected")
    d_ax.plot(d["scan_amps"], 100 * d["infid_cor"], marker="^", **STYLE_CORRECTED,
              lw=LW_MAIN, label="corrected")
    d_ax.set_xlabel("Distortion amplitude $A$", fontsize=11)
    d_ax.set_ylabel("Gate infidelity (%)", fontsize=11)
    d_ax.set_title("(d) Robustness to control-line drift\n"
                   f"(calibrated at A={TRUE_AMPS[0]:.3f})", fontsize=11)
    d_ax.legend(fontsize=8.5)
    d_ax.grid(True, ls=":", alpha=0.4)

    # Shade the "improvement zone"
    d_ax.fill_between(d["scan_amps"], 100 * d["infid_cor"],
                      100 * d["infid_unc"],
                      alpha=0.15, color="#377eb8",
                      label="improvement")

    fig.suptitle("Predistortion recovery of a calibrated CZ gate", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.965))

    for fmt in ["png", "pdf"]:
        fig.savefig(os.path.join(out_dir, f"{_CACHE}.{fmt}"),
                    dpi=300 if fmt == "png" else None)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Environment recording
# ---------------------------------------------------------------------------

def _record_environment(d, metrics):
    import sys as _sys
    import numpy as _np
    import scipy as _sp
    import qutip as _qt
    import matplotlib as _mpl

    env = {
        "timestamp": datetime.datetime.now().isoformat(),
        "python": _sys.version,
        "numpy": _np.__version__,
        "scipy": _sp.__version__,
        "qutip": _qt.__version__,
        "matplotlib": _mpl.__version__,
        "config": {
            "cz_flux_amp": CZ_FLUX_AMP,
            "cz_duration": CZ_DURATION,
            "cz_rise_frac": CZ_RISE_FRAC,
            "q2_flux_bias": Q2_FLUX_BIAS,
            "g_coupling": G_COUPLING,
            "true_amps": TRUE_AMPS,
            "true_taus": TRUE_TAUS,
            "t_max_step": T_MAX_STEP,
            "step_amp": STEP_AMP,
            "fit_type": FIT_TYPE,
            "scan_amps": SCAN_AMPS.tolist(),
            "scan_tau_fixed": SCAN_TAU_FIXED,
        },
        "metrics_summary": {k: float(v) if isinstance(v, (np.floating, float))
                            else v for k, v in metrics.items()},
    }
    json_path = os.path.join(C.out_dir(SUBDIR),
                             "cz_gate_impact_metrics.json")
    with open(json_path, "w") as f:
        json.dump(env, f, indent=2, default=str)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    d = compute(recompute="--recompute" in sys.argv)
    metrics = _compute_metrics(d)

    # Print report
    print("=== P7 CZ Gate Impact Metrics ===")
    print(f"  CZ flux amp: {CZ_FLUX_AMP:.4f} Phi_0, coupling g={G_COUPLING:.3f} GHz")
    print(f"  Cryoscope fit: A={d['fitted_amplitude'][0]:.5f}, "
          f"τ={d['fitted_tau'][0]:.2f} ns")
    for label, prefix in [("Ideal baseline  ", "ideal"),
                           ("Uncorrected     ", "uncorrected"),
                           ("Corrected       ", "corrected")]:
        print(f"  {label}: cond_phase={metrics[f'cond_phase_{prefix}']:.4f} rad, "
              f"phase_err={metrics[f'phase_error_{prefix}']:.4f}, "
              f"leak={metrics[f'leakage_{prefix}']:.4f}, "
              f"infid={metrics[f'infidelity_{prefix}']:.4f}")
    print(f"  Additional phase error: unc +{metrics['phase_error_added_unc']:.4f}, "
          f"cor +{metrics['phase_error_added_cor']:.4f} rad")
    print(f"  Additional infidelity:  unc +{metrics['infidelity_added_unc']:.4f}, "
          f"cor +{metrics['infidelity_added_cor']:.4f}")
    print(f"  Flux RMSE: uncorrected={metrics['rmse_chip_unc']:.2e}, "
          f"corrected={metrics['rmse_chip_cor']:.2e}")

    _record_environment(d, metrics)
    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
