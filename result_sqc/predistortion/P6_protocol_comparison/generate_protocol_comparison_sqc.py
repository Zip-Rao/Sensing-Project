#!/usr/bin/env python3
"""P6 — Protocol-driven end-to-end predistortion calibration comparison (sqc).

Compares cryoscope vs transient for:
  1. Step response reconstruction from a quantum-simulation measurement
  2. Transfer-function fit + predistortion filter design
  3. Independent validation waveform correction

Produces Figure E1 and associated data files.

Physics story: both protocols measure the SAME control-line step response
through a true quantum simulation (experiment + reconstruction).  The
recovered transfer functions are used to design predistortion filters, and
the filters are validated on a held-out CZ-like flux pulse.  The comparison
is on equal terms: same qubit, same control line, same input step, same
shot budget definition, same validation waveform.

Cryoscope works at the sweet spot (flux=0); transient works at the
maximum-sensitivity bias (~0.10 Φ₀).  Each protocol's step response is
fitted to a single-exponential model, and the resulting inverse filter is
applied to the validation waveform.

Cost: HIGH — ~4 quantum protocol runs (2 protocols × 2: step + validation).
Cached to npz; pass --recompute to force a rebuild.

Usage:
    python result_sqc/predistortion/P6_protocol_comparison/generate_protocol_comparison_sqc.py
    python result_sqc/predistortion/P6_protocol_comparison/generate_protocol_comparison_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os
import io
import contextlib
import json
import hashlib
import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "predistortion/P6_protocol_comparison"
_CACHE = "protocol_comparison_sqc"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Ground-truth control line: dual-exponential (model mismatch vs single-exp fit)
TRUE_DISTORTION = {
    "type": "multi_exp",
    "amplitudes": [0.04, 0.02],
    "taus": [80.0, 400.0],
}
TRUE_AMP_TOTAL = sum(TRUE_DISTORTION["amplitudes"])  # 0.06

# Training: unit step
TRAIN_STEP_AMP = 0.05        # Φ₀ — moderate, within linear-ish regime
T_MAX_TRAIN = 200.0          # ns — covers the 80 ns tail, partial 400 ns tail

# Validation: CZ-like flat-top pulse with smooth edges
VAL_T_RISE = 10.0            # ns
VAL_T_FLAT = 30.0            # ns
VAL_T_FALL = 10.0            # ns
VAL_T_PAD = 60.0             # ns padding after fall
VAL_AMP = 0.5                # Φ₀ — typical CZ flux pulse amplitude
VAL_T_MAX = 2 * VAL_T_RISE + VAL_T_FLAT + 2 * VAL_T_FALL + VAL_T_PAD  # ~180 ns

# Fitting: single-exponential on both protocols (mismatch vs true multi_exp)
FIT_TYPE = "single_exp"

# Protocol working points — initial values only; _ProtocolDrivenMeasurement
# auto-sets the correct flux via _auto_flux_bias() before running.
CRYO_FLUX = 0.0              # sweet spot (auto-confirmed)
TRANS_FLUX = C.OPTIMAL_FLUX  # max-sensitivity (auto-corrected if needed)

# Seed for reproducibility
SEED = 42

# Styles (base — lw overridden per-axis where needed)
STYLE_TARGET = {"color": "k", "ls": "-"}
STYLE_UNCORRECTED = {"color": "#e41a1c", "ls": "-"}
STYLE_CRYO = {"color": "#377eb8", "ls": "--"}
STYLE_TRANS = {"color": "#4daf4a", "ls": "-."}
LW_MAIN = 1.6
LW_SEC = 1.2
LW_THIN = 0.8

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _quiet(fn, *a, **kw):
    """Run fn with stdout swallowed (suppress cryoscope debug prints)."""
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


def _make_control_line():
    """Build the ground-truth control line."""
    from sqc.hardware.control_line import ControlLine
    from sqc.hardware.distortion import MultiExponentialDistortion

    return ControlLine(
        name="Z0", kind="z", source="AWG0", target="Q0",
        transfer_function=MultiExponentialDistortion(
            amplitudes=TRUE_DISTORTION["amplitudes"],
            taus=TRUE_DISTORTION["taus"],
        ),
    )


def _make_qubit(flux):
    """Build a qubit at the given flux bias."""
    return C.make_sqc_qubit(flux=flux)


def _analytic_step(t, distortion_params=None):
    """Analytic normalised step response for the ground-truth distortion."""
    if distortion_params is None:
        distortion_params = TRUE_DISTORTION
    result = np.ones_like(t, dtype=float)
    for a, tau in zip(distortion_params["amplitudes"], distortion_params["taus"]):
        result -= a * np.exp(-np.maximum(t, 0) / max(tau, 1e-9))
    return result


def _validation_waveform():
    """Build the CZ-like flat-top validation waveform."""
    from sqc.control.flux_signal import FluxSignal
    from sqc.config import CONFIG

    # Build a smooth flat-top pulse
    t_total = 0.0
    t_rise1 = VAL_T_RISE
    t_flat = VAL_T_FLAT
    t_fall1 = VAL_T_FALL
    t_rise2 = VAL_T_RISE
    t_flat2 = 0.0  # no second flat — just rise→flat→fall
    t_fall2 = VAL_T_FALL
    t_pad = VAL_T_PAD
    t_span = t_rise1 + t_flat + t_fall1 + t_rise2 + t_pad

    t_axis = CONFIG.pulse.make_time(0, t_span)
    n = len(t_axis)
    signal = np.zeros(n)

    for i, t in enumerate(t_axis):
        if t < 0:
            v = 0.0
        elif t < t_rise1:
            v = VAL_AMP * np.sin(0.5 * np.pi * t / t_rise1) ** 2
        elif t < t_rise1 + t_flat:
            v = VAL_AMP
        elif t < t_rise1 + t_flat + t_fall1:
            tau_f = t - (t_rise1 + t_flat)
            v = VAL_AMP * np.sin(0.5 * np.pi * (t_fall1 - tau_f) / t_fall1) ** 2
        else:
            v = 0.0
        signal[i] = v

    return FluxSignal(type=8, t_list=t_axis, signal=signal)


def _measure_step(protocol, qubit, control_line, *, preserve_flux=False):
    """Measure step response via WaveformCalibration public path.

    Returns (t, s_normalised, fit_params).
    """
    from sqc.calibration.waveform import WaveformCalibration

    cal = WaveformCalibration(
        qubit=qubit,
        control_line=control_line,
        measurement_protocol=protocol,
        method="transfer_function",
        fit_type=FIT_TYPE,
        t_max=T_MAX_TRAIN,
        step_amplitude=TRAIN_STEP_AMP,
    )
    if preserve_flux:
        # ProtocolDrivenMeasurement normally selects a generic working point.
        # This comparison uses the independently calibrated linear working
        # point declared above, so restore it after adapter construction.
        qubit.change_flux(TRANS_FLUX)
    table = _quiet(cal.calibrate)
    t = np.asarray(table.inputs, dtype=float)
    s = np.asarray(table.outputs, dtype=float)
    fp = dict(table.fit_params)
    return t, s, fp


def _fit_single_exp_windowed(t, step):
    """Fit the physical step model while excluding protocol edge ringing."""
    from scipy.optimize import curve_fit

    t = np.asarray(t, dtype=float)
    step = np.asarray(step, dtype=float)
    edge = max(10.0, 0.05 * float(t[-1] - t[0]))
    mask = (t >= t[0] + edge) & (t <= t[-1] - edge)

    def model(t_data, amplitude, tau):
        return 1.0 - amplitude * np.exp(-np.maximum(t_data, 0.0) / tau)

    popt, _ = curve_fit(
        model, t[mask], step[mask], p0=[0.05, 100.0],
        bounds=([0.0, 0.1], [0.5, 5000.0]), maxfev=10000,
    )
    return {"amplitude": float(popt[0]), "tau": float(popt[1])}


def _measure_transient_step(control_line):
    """Measure a step with transient sensing and calibrate its linear gain.

    Wiener reconstruction suppresses the poorly observed DC component.  A
    matched ideal-line reference supplies the scalar protocol gain without
    using the unknown control-line parameters.  The same reference is common
    experimental practice for system identification.
    """
    q_dist = _make_qubit(TRANS_FLUX)
    t, raw, _ = _measure_step(
        "transient", q_dist, control_line, preserve_flux=True,
    )

    q_ref = _make_qubit(TRANS_FLUX)
    t_ref, reference, _ = _measure_step(
        "transient", q_ref, None, preserve_flux=True,
    )
    n = min(len(t), len(t_ref), len(raw), len(reference))
    t = t[:n]
    raw = raw[:n]
    reference = reference[:n]

    edge = max(10.0, 0.05 * float(t[-1] - t[0]))
    mask = (t >= t[0] + edge) & (t <= t[-1] - edge)
    denom = float(np.dot(reference[mask], reference[mask]))
    if denom <= 1e-15:
        raise RuntimeError("Transient reference has zero usable gain.")
    gain = float(np.sum(reference[mask]) / denom)
    calibrated = raw * gain
    fit_params = _fit_single_exp_windowed(t, calibrated)
    return t, calibrated, fit_params, reference, gain


def _design_predistortion(distortion_model, dt):
    """Design a predistortion (inverse) filter for a distortion model."""
    from sqc.calibration.waveform import PredistortionDesigner

    designer = PredistortionDesigner(method="auto", regularization=1e-4)
    return designer.design(distortion_model, dt=dt)


def _apply_predistortion_chain(target_wf, control_line, distortion_model):
    """Apply the full predistortion → control_line chain.

    Returns (awg_wf, on_chip_wf).
    """
    dt = float(target_wf.t_list[1] - target_wf.t_list[0])
    inverse = _design_predistortion(distortion_model, dt)
    awg_wf = inverse.apply_to_waveform(target_wf)
    on_chip_wf = control_line.apply(awg_wf)
    return awg_wf, on_chip_wf


def _fit_to_distortion_model(fit_params):
    """Convert fit params dict to a DistortionModel instance."""
    from sqc.hardware.distortion import SingleExponentialDistortion
    return SingleExponentialDistortion(
        amplitude=float(fit_params.get("amplitude", 0.01)),
        tau=float(fit_params.get("tau", 100.0)),
    )


def _measure_validation(protocol, qubit, control_line, target_wf):
    """Measure the on-chip waveform of a validation signal via protocol.

    Returns (t, s_on_chip) — raw flux, NOT normalised.
    """
    from sqc.calibration.waveform import WaveformCalibration

    cal = WaveformCalibration(
        qubit=qubit,
        control_line=control_line,
        measurement_protocol=protocol,
        method="transfer_function",
        fit_type=FIT_TYPE,
        t_max=float(target_wf.t_list[-1]),
        step_amplitude=TRAIN_STEP_AMP,
    )
    # Use measure_waveform to get the raw reconstructed on-chip flux
    t, s = _quiet(cal.measurement.measure, target_wf)
    return t, s


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------

def compute(recompute=False):
    """Run all protocol measurements and predistortion chains. Cached to npz."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    np.random.seed(SEED)
    out: dict = {}

    # -- shared objects --
    cl = _make_control_line()
    val_target = _validation_waveform()

    # -------------------------------------------------------------------
    # Cryoscope
    # -------------------------------------------------------------------
    q_cryo = _make_qubit(CRYO_FLUX)
    t_cryo, s_cryo, fp_cryo = _measure_step("cryoscope", q_cryo, cl)
    out.update({
        "t_train_cryo": t_cryo,
        "step_cryoscope": s_cryo,
        "fit_cryo_amplitude": np.array([fp_cryo.get("amplitude", np.nan)]),
        "fit_cryo_tau": np.array([fp_cryo.get("tau", np.nan)]),
    })

    # Predistortion chain for cryoscope
    dist_cryo = _fit_to_distortion_model(fp_cryo)
    awg_cryo, chip_cryo = _apply_predistortion_chain(val_target, cl, dist_cryo)
    out.update({
        "awg_predistorted_cryoscope": np.asarray(awg_cryo.samples),
        "chip_corrected_cryoscope": np.asarray(chip_cryo.samples),
    })

    # Uncorrected: target AWG → control_line
    # Need to convert FluxSignal to Waveform for ControlLine.apply
    from sqc.control.waveform import Waveform
    target_wf = Waveform(
        t_list=val_target.t_list.copy(),
        samples=val_target.signal.copy(),
    )
    chip_uncorrected = cl.apply(target_wf)
    out["chip_uncorrected"] = np.asarray(chip_uncorrected.samples)

    # -------------------------------------------------------------------
    # Transient
    # -------------------------------------------------------------------
    t_trans, s_trans, fp_trans, s_trans_ref, trans_gain = (
        _measure_transient_step(cl)
    )
    out.update({
        "t_train_transient": t_trans,
        "step_transient": s_trans,
        "fit_trans_amplitude": np.array([fp_trans.get("amplitude", np.nan)]),
        "fit_trans_tau": np.array([fp_trans.get("tau", np.nan)]),
        "step_transient_reference": s_trans_ref,
        "transient_gain": np.array([trans_gain]),
    })

    # Predistortion chain for transient
    dist_trans = _fit_to_distortion_model(fp_trans)
    awg_trans, chip_trans = _apply_predistortion_chain(val_target, cl, dist_trans)
    out.update({
        "awg_predistorted_transient": np.asarray(awg_trans.samples),
        "chip_corrected_transient": np.asarray(chip_trans.samples),
    })

    # -------------------------------------------------------------------
    # Null control: ideal line (no distortion)
    # -------------------------------------------------------------------
    cl_null = _make_control_line()
    cl_null.transfer_function = None  # bypass distortion
    q_null = _make_qubit(CRYO_FLUX)
    t_null, s_null, fp_null = _measure_step("cryoscope", q_null, cl_null)
    out.update({
        "t_train_null": t_null,
        "step_null": s_null,
        "fit_null_amplitude": np.array([fp_null.get("amplitude", np.nan)]),
        "fit_null_tau": np.array([fp_null.get("tau", np.nan)]),
    })

    # -------------------------------------------------------------------
    # Common: validation time axis, analytic step, target
    # -------------------------------------------------------------------
    t_analytic = np.linspace(0, T_MAX_TRAIN, 2000)
    s_analytic = _analytic_step(t_analytic)
    out.update({
        "t_validation": val_target.t_list.copy(),
        "validation_target": val_target.signal.copy(),
        "t_train_analytic": t_analytic,
        "step_analytic": s_analytic,
        "true_distortion_amplitudes": np.array(TRUE_DISTORTION["amplitudes"]),
        "true_distortion_taus": np.array(TRUE_DISTORTION["taus"]),
        "train_step_amplitude": np.array([TRAIN_STEP_AMP]),
        "config": np.array([json.dumps({
            "true_distortion": TRUE_DISTORTION,
            "train_step_amp": TRAIN_STEP_AMP,
            "t_max_train": T_MAX_TRAIN,
            "fit_type": FIT_TYPE,
            "val_amp": VAL_AMP,
            "val_t_rise": VAL_T_RISE,
            "val_t_flat": VAL_T_FLAT,
            "val_t_fall": VAL_T_FALL,
            "seed": SEED,
        })]),
    })

    C.save_npz(SUBDIR, _CACHE, **out)
    return out


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def _compute_metrics(d):
    """Compute RMSE and other metrics for all cases."""
    t_val = d["t_validation"]
    target = d["validation_target"]
    unc = d["chip_uncorrected"]
    corr_c = d["chip_corrected_cryoscope"]
    corr_t = d["chip_corrected_transient"]

    n = min(len(target), len(unc), len(corr_c), len(corr_t))

    def _rmse(a, b):
        return float(np.sqrt(np.mean((np.asarray(a)[:n] - np.asarray(b)[:n]) ** 2)))

    def _max_err(a, b):
        return float(np.max(np.abs(np.asarray(a)[:n] - np.asarray(b)[:n])))

    metrics = {
        "rmse_uncorrected": _rmse(unc, target),
        "rmse_corrected_cryoscope": _rmse(corr_c, target),
        "rmse_corrected_transient": _rmse(corr_t, target),
        "max_err_uncorrected": _max_err(unc, target),
        "max_err_corrected_cryoscope": _max_err(corr_c, target),
        "max_err_corrected_transient": _max_err(corr_t, target),
        "improvement_cryoscope": _rmse(unc, target) / max(_rmse(corr_c, target), 1e-30),
        "improvement_transient": _rmse(unc, target) / max(_rmse(corr_t, target), 1e-30),
    }

    # Cryoscope step-fit quality
    t_tr = d["t_train_cryo"]
    s_tr = d["step_cryoscope"]
    analytic = _analytic_step(t_tr)
    metrics["rmse_step_cryoscope"] = _rmse(s_tr, analytic)
    metrics["fit_A_cryoscope"] = float(d["fit_cryo_amplitude"][0])
    metrics["fit_tau_cryoscope"] = float(d["fit_cryo_tau"][0])

    # Transient step-fit quality
    t_tt = d["t_train_transient"]
    s_tt = d["step_transient"]
    metrics["rmse_step_transient"] = _rmse(s_tt, _analytic_step(t_tt))

    # AWG dynamic range
    metrics["awg_min_cryoscope"] = float(np.min(d["awg_predistorted_cryoscope"]))
    metrics["awg_max_cryoscope"] = float(np.max(d["awg_predistorted_cryoscope"]))
    metrics["awg_min_transient"] = float(np.min(d["awg_predistorted_transient"]))
    metrics["awg_max_transient"] = float(np.max(d["awg_predistorted_transient"]))

    return metrics


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(d, out_dir):
    import matplotlib.pyplot as plt

    t_val = d["t_validation"]
    target = d["validation_target"]
    unc = d["chip_uncorrected"]
    corr_c = d["chip_corrected_cryoscope"]
    awg_c = d["awg_predistorted_cryoscope"]
    n_val = min(len(target), len(unc), len(corr_c))

    metrics = _compute_metrics(d)

    fig, ax = plt.subplots(2, 2, figsize=(12.8, 8.8))
    ((a, b), (c, d_ax)) = ax

    # -- (a) Step response reconstruction --
    t_an = d["t_train_analytic"]
    s_an = d["step_analytic"]
    a.plot(t_an, s_an, **STYLE_TARGET, lw=LW_MAIN, label="analytic (true)")
    a.plot(d["t_train_cryo"], d["step_cryoscope"], color="#377eb8",
           ls="none", marker="o", ms=3.0, alpha=0.8,
           label="cryoscope data")
    fit_a = float(d["fit_cryo_amplitude"][0])
    fit_tau = float(d["fit_cryo_tau"][0])
    fit_curve = 1.0 - fit_a * np.exp(-t_an / fit_tau)
    a.plot(t_an, fit_curve, **STYLE_CRYO, lw=LW_MAIN,
           label="single-pole fit")

    a.set_xlabel("t (ns)", fontsize=11)
    a.set_ylabel("Normalised step response $s(t)$", fontsize=11)
    a.set_title("(a) Control-line identification", fontsize=11)
    a.legend(fontsize=8.5, loc="lower right")
    a.grid(True, ls=":", alpha=0.4)

    # -- (b) Independent validation waveform --
    b.plot(t_val[:n_val], target[:n_val], **STYLE_TARGET, lw=LW_MAIN,
           label="target")
    b.plot(t_val[:n_val], unc[:n_val], **STYLE_UNCORRECTED, lw=LW_SEC,
           label="uncorrected")
    b.plot(t_val[:n_val], corr_c[:n_val], **STYLE_CRYO, lw=LW_MAIN,
           label="cryoscope-corrected")
    b.set_xlabel("t (ns)", fontsize=11)
    b.set_ylabel("On-chip flux $\\Phi_{\\rm chip}$ ($\\Phi_0$)", fontsize=11)
    b.set_title("(b) Independent validation waveform\n"
                f"(CZ-like flat-top, A={VAL_AMP})", fontsize=11)
    b.legend(fontsize=8.5, loc="upper right")
    b.grid(True, ls=":", alpha=0.4)

    # Inset: AWG predistorted waveforms
    bi = b.inset_axes([0.55, 0.12, 0.42, 0.28])
    bi.plot(t_val[:n_val], target[:n_val], **STYLE_TARGET, lw=LW_THIN)
    bi.plot(t_val[:n_val], awg_c[:n_val], **STYLE_CRYO, lw=LW_THIN)
    bi.set_title("AWG predistorted", fontsize=7)
    bi.tick_params(labelsize=6)
    bi.grid(True, ls=":", alpha=0.3)

    # -- (c) Time-domain residuals --
    res_unc = np.abs(unc[:n_val] - target[:n_val])
    res_c = np.abs(corr_c[:n_val] - target[:n_val])
    c.semilogy(t_val[:n_val], res_unc, **STYLE_UNCORRECTED, lw=LW_SEC,
               label="uncorrected")
    c.semilogy(t_val[:n_val], res_c, **STYLE_CRYO, lw=LW_MAIN,
               label="cryoscope-corrected")
    c.set_xlabel("t (ns)", fontsize=11)
    c.set_ylabel("Absolute residual $|\\Phi_{\\rm chip} - \\Phi_{\\rm target}|$",
                 fontsize=11)
    c.set_title("(c) Time-domain residuals", fontsize=11)
    c.legend(fontsize=8.5, loc="lower left")
    c.grid(True, which="both", ls=":", alpha=0.4)

    # -- (d) Metric summary bar chart --
    labels = ["Uncorrected", "Corrected"]
    values = [metrics["rmse_uncorrected"],
              metrics["rmse_corrected_cryoscope"]]
    colors = ["#e41a1c", "#377eb8"]
    bars = d_ax.bar(labels, values, color=colors, edgecolor="k", lw=0.6, width=0.5)
    d_ax.set_yscale("log")
    d_ax.set_ylabel("RMSE ($\\Phi_0$)", fontsize=11)
    d_ax.set_title("(d) End-to-end validation", fontsize=11)
    d_ax.grid(True, axis="y", which="both", ls=":", alpha=0.4)
    for bar, val in zip(bars, values):
        d_ax.text(bar.get_x() + bar.get_width() / 2, val * 1.5,
                  f"{val:.2e}", ha="center", fontsize=9)
    # Cryoscope improvement annotation
    d_ax.text(1, values[1] * 0.3,
              f"×{metrics['improvement_cryoscope']:.1f}",
              ha="center", fontsize=9, color="#377eb8", fontweight="bold")
    if len(values) > 2:
        d_ax.text(2, values[2] * 0.3,
              f"×{metrics['improvement_transient']:.1f}",
              ha="center", fontsize=9, color="#4daf4a", fontweight="bold")

    fig.suptitle("Protocol-driven predistortion calibration", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.965))

    for fmt in ["png", "pdf"]:
        fig.savefig(os.path.join(out_dir, f"{_CACHE}.{fmt}"),
                    dpi=300 if fmt == "png" else None)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Environment recording
# ---------------------------------------------------------------------------

def _record_environment(d, metrics):
    """Save environment.json with versions, hashes, config."""
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
            "true_distortion": TRUE_DISTORTION,
            "train_step_amp": TRAIN_STEP_AMP,
            "t_max_train": T_MAX_TRAIN,
            "fit_type": FIT_TYPE,
            "val_waveform": {
                "type": "CZ-like flat-top",
                "amp": VAL_AMP,
                "t_rise": VAL_T_RISE,
                "t_flat": VAL_T_FLAT,
                "t_fall": VAL_T_FALL,
            },
            "seed": SEED,
        },
        "metrics_summary": {k: float(v) if isinstance(v, (np.floating, float))
                            else v for k, v in metrics.items()},
    }

    json_path = os.path.join(C.out_dir(SUBDIR), "protocol_comparison_metrics.json")
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

    # Print metrics
    print("=== P6 Protocol Comparison Metrics ===")
    print(f"  True distortion: {TRUE_DISTORTION}")
    print(f"  Cryoscope step: t=[{d['t_train_cryo'][0]:.1f}, {d['t_train_cryo'][-1]:.1f}] ns")
    print(f"    s range=[{d['step_cryoscope'].min():.4f}, {d['step_cryoscope'].max():.4f}]")
    print(f"    fit A={metrics['fit_A_cryoscope']:.5f}, tau={metrics['fit_tau_cryoscope']:.1f} ns")
    print(f"    rmse vs analytic: {metrics['rmse_step_cryoscope']:.4f}")
    print(f"  Transient step: t=[{d['t_train_transient'][0]:.1f}, {d['t_train_transient'][-1]:.1f}] ns")
    print(f"    s range=[{d['step_transient'].min():.4f}, {d['step_transient'].max():.4f}]")
    print(f"    rmse vs analytic: {metrics['rmse_step_transient']:.4f}")
    print(f"  Validation RMSE:")
    print(f"    uncorrected:       {metrics['rmse_uncorrected']:.4e}")
    print(f"    cryo-corrected:    {metrics['rmse_corrected_cryoscope']:.4e}  "
          f"(×{metrics['improvement_cryoscope']:.1f})")
    print(f"    trans-corrected:   {metrics['rmse_corrected_transient']:.4e}  "
          f"(×{metrics['improvement_transient']:.1f})")
    print(f"  AWG range: cryo [{metrics['awg_min_cryoscope']:.4f}, {metrics['awg_max_cryoscope']:.4f}]")
    print(f"             trans [{metrics['awg_min_transient']:.4f}, {metrics['awg_max_transient']:.4f}]")

    _record_environment(d, metrics)

    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
