#!/usr/bin/env python3
"""P5 ⭐ — Protocol-driven vs analytic transfer-function measurement (sqc).

sqc source:
    sqc.calibration.waveform.WaveformCalibration(measurement_protocol=...)
    measures the step response through a ControlLine via a TRUE
    quantum-simulation protocol (experiment + reconstruction), versus the
    analytical distortion.step_response(t).

Physics story (⭐ paper role): connect predistortion back to real qubit
protocols — and report honestly how far that connection currently reaches.
The cryoscope-measured step response tracks the analytic transfer function to
rmse 3.7e-3 and recovers (A, tau) to (-8.6%, +11.4%) at tau=100 ns, so the
P2–P4 loop IS driveable from a simulated measurement rather than a formula.
But the accompanying control experiment (panel b) shows that residual is
dominated by a protocol artifact, not by distortion-measurement error, so
3.7e-3 must NOT be quoted as the protocol's accuracy. An earlier draft of this
script called it a "measurement accuracy floor"; that was wrong and the control
experiment is what disproved it.

CRYOSCOPE ONLY, and that is a measured result, not a shortcut. Of the four
protocols registered in _build_protocol_registry, only cryoscope returns a
usable normalised step response through this public path (all verified here at
A=0.05, tau=20 ns, t_max=200 ns):
    cryoscope     3.2 s   window 20.5..90 ns    rmse vs analytic 4.2e-3   OK
    delay_ramsey  9.6 s   window 0..199.5 ns    rmse 0.73..0.94 — the
                  reconstruction decays to ~0 instead of staying near 1, and
                  t_max has NO effect (identical output for t_max=30/60/100)
    transient     —       crashes: TypeError: object of type 'NoneType' has
                  no len()
    pi_pulse      47.9 s  window 0..99.5 ns     rmse 0.82, range 0..0.2
See FIGURES.md P5 for the full write-up. Those three are sqc defects, out of
scope for this figure; do NOT "fix" them from here.

TWO INDEPENDENT ERROR SOURCES, which is the physics this figure isolates:
    1. a PROTOCOL ARTIFACT, not a noise floor. rmse vs analytic is
       ~3.6e-3..4.2e-3 and barely moves with tau — because it is not really
       measuring the tail at all. Control experiment (panel b): run the same
       cryoscope on an IDEAL line with NO distortion, where the true step
       response is exactly s(t)=1. It still comes back oscillating, rms
       5.4e-3 at step_amplitude=0.05 — LARGER than the 3.7e-3 "error" against
       the distorted line. The residual is therefore dominated by the protocol,
       not by distortion-measurement error, and calling it a noise floor would
       be wrong: it is coherent, with a well-defined period.
       Decomposing the artifact on the ideal line (all measured here):
         step_amp   DC bias     oscillation std   FFT period   phase/pi over tau
           0.002    +3.45e-3    4.20e-4           (unresolved)   0.008
           0.005    +3.77e-3    5.43e-4           (unresolved)   0.049
           0.010    +4.32e-3    2.88e-4           (unresolved)   0.197
           0.020    +3.46e-3    1.14e-3           (unresolved)   0.790
           0.050    +5.51e-4    5.39e-3            35 ns         4.940
           0.100    (rms 9.6e-1 — reconstruction collapses)   5 ns   19.823
       Two regimes, split where the phase accumulated over the cryoscope delay
       tau=100 ns crosses pi (between step_amp 0.02 and 0.05):
         - DC-dominated below: a roughly amplitude-INDEPENDENT bias ~+3.5e-3.
         - oscillation-dominated above, and the period matches the phase-winding
           prediction 2*pi/dw(a), where dw(a) = w(a) - w(0) is QUADRATIC in flux
           because dw/dPhi = 0 at the sweet spot (verified: dw scales x4.0 per
           doubling of a). Predicted 40.5 ns vs observed 35 ns at a=0.05;
           predicted 10.1 ns vs observed 5.0 ns at a=0.10 — agreement to within
           the coarse trunc grid.
       IMPORTANT — what the mechanism is NOT: the model-guided unwrap
       (unwrap_phase_with_model, cryoscope.py:175) does NOT fail at a=0.05.
       Checked on the intermediate varphi: it recovers the phase slope to +0.1%
       across a 3.44*pi span. Only a=0.10 (13.96*pi span) breaks it
       (nonlinearity rms 1.63 rad vs 1.0e-2 rad at a=0.05). So at the working
       point of this figure the artifact comes from the ~1e-2 rad NONLINEARITY
       left in varphi being amplified by the d/dt gradient, not from a wrap
       failure. (Also not from dividing by a vanishing sensitivity:
       invert_frequency_to_flux solves the exact transmon dispersion by
       bisection, it never divides by dw/dPhi.)
       Consequence: at the step_amplitude=0.05 used here the artifact peak
       (+-7.6e-3) is ~37% of the true in-window step deficit
       (0.05*(exp(-20.5/100) - exp(-90/100)) = 0.0204), which is what biases the
       fitted A low by 8.6%. step_amplitude=0.10 is unusable outright.
    2. parameter recoverability — set by how much tail amplitude is left at the
       window start (exp(-20.5/tau)): A error goes -49.8% (tau=20) -> -5.4%
       (tau=200) while rmse stays flat. Orthogonal to (1).

Produces:
    result_sqc/predistortion/P5_protocol_driven/protocol_driven_measurement_sqc.png / .pdf
        (a) quantum-measured vs analytic step (+ fitted model), tau=100 ns
        (b) residual vs the IDEAL-line control -> it is a protocol artifact,
            with an inset of artifact rms vs step_amplitude (operating range)
        (c) window validity: the reconstruction is only valid on 20.5..90 ns
        (d) tau sweep: recovery error vs surviving tail fraction
    result_sqc/predistortion/P5_protocol_driven/protocol_driven_measurement_sqc.npz

Cost: HIGH-ish — 8 quantum protocol runs (mesolve), ~40 s total. Cached to npz;
pass --recompute to force a rebuild.

Usage:
    python result_sqc/predistortion/P5_protocol_driven/plot_protocol_driven_measurement_sqc.py
    python result_sqc/predistortion/P5_protocol_driven/plot_protocol_driven_measurement_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os
import io
import contextlib

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "predistortion/P5_protocol_driven"
_CACHE = "protocol_driven_measurement_sqc"

PROTOCOL = "cryoscope"
# Cryoscope requires the SWEET SPOT (flux=0) — _warn_flux_bias warns otherwise.
# This is the one place in result_sqc/ that does NOT sit at C.OPTIMAL_FLUX.
CRYO_FLUX = 0.0

TRUE_AMP = 0.05
TRUE_TAU_NS = 100.0        # main panel; matches P1–P4's canonical tail
TAU_SWEEP = (20.0, 50.0, 100.0, 200.0)
T_MAX = 200.0

# The default trunc_list is flux_signal.t_list[180:40:-1] (cryoscope.py:80),
# a HARDCODED index slice -> at dt=0.5 ns the measured window is always
# 20.5..90 ns regardless of t_max. Verified that this is not an arbitrary
# default but the reconstruction's true validity range: pushing the upper edge
# out breaks it (rmse 4.2e-3 -> 1.6e-1 at 120 ns -> 6.7e-1 at 150 ns, where the
# response inversion blows up to 6.2). Panel (c) shows this.
WINDOW_SLICES = ((180, 40), (240, 40), (300, 40))
WINDOW_T_MAX = 250.0       # signal must outlast the widest slice tested

# step_amplitude scan for the artifact control experiment. 0.10 is included
# deliberately: that is where the reconstruction collapses (artifact rms 0.96,
# i.e. ~100% error) and it sets the protocol's usable operating range.
STEP_AMP_SCAN = (0.002, 0.005, 0.01, 0.02, 0.05, 0.10)

STYLE_MEAS = {"color": "#377eb8", "marker": "o"}
STYLE_ANA = {"color": "k", "ls": "-"}
STYLE_FIT = {"color": "#e41a1c", "ls": "--"}


def _quiet(fn, *a, **kw):
    """Run fn with stdout swallowed.

    sqc.experiments.cryoscope has two stray debug prints — a fake
    'Warning: trunc_list not provided...' (cryoscope.py:79, whose message is
    also stale: it says [140:20:-1] while the code does [180:40:-1], and the
    string is unterminated) and print("trunc_list:", ...) (cryoscope.py:116).
    They dump the whole trunc_list array into stdout and would bury this
    script's own output. Suppressed here rather than patched: editing
    sqc/experiments/ is Track B work, out of scope for a figure script.
    """
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


def _line(tau):
    """ControlLine carrying the ground-truth single-exp distortion."""
    from sqc.hardware.control_line import ControlLine
    from sqc.hardware.distortion import SingleExponentialDistortion

    return ControlLine(name="Z0", kind="z", source="AWG0", target="Q0",
                       transfer_function=SingleExponentialDistortion(
                           amplitude=TRUE_AMP, tau=tau))


def measure_public(tau):
    """Cryoscope step response via the public WaveformCalibration path.

    Returns (t, measured_step, fit_amplitude, fit_tau). The step is already
    normalised by step_amplitude inside _ProtocolDrivenMeasurement.measure().
    """
    from sqc.calibration.waveform import WaveformCalibration

    cal = WaveformCalibration(
        qubit=C.make_sqc_qubit(flux=CRYO_FLUX), control_line=_line(tau),
        measurement_protocol=PROTOCOL, method="transfer_function",
        fit_type="single_exp", t_max=T_MAX, step_amplitude=TRUE_AMP)
    table = _quiet(cal.calibrate)
    fp = table.fit_params
    return (np.asarray(table.inputs, dtype=float),
            np.asarray(table.outputs, dtype=float),
            float(fp["amplitude"]), float(fp["tau"]))


def measure_window(tau, hi, lo):
    """Same measurement with an EXPLICIT trunc_list window.

    Needed for panel (c): _ProtocolDrivenMeasurement._build_exp_kwargs does not
    forward trunc_list, so the window can only be varied by driving
    CryoscopeExperiment / CryoscopeReconstruction directly. Same machinery,
    just with the window under our control.
    """
    from sqc.experiments.cryoscope import CryoscopeExperiment
    from sqc.reconstruction.cryoscope import CryoscopeReconstruction
    from sqc.control.flux_signal import FluxSignal
    from sqc.config import CONFIG

    tau_c = CONFIG.reconstruction.cryoscope_tau
    t_ax = CONFIG.pulse.make_time(0, WINDOW_T_MAX)
    q = C.make_sqc_qubit(flux=CRYO_FLUX)
    sig = FluxSignal(type=1, t_list=t_ax, amplitude=TRUE_AMP)

    def _run():
        res = CryoscopeExperiment(qubit=q, control_line=_line(tau),
                                  flux_signal=sig, tau=tau_c,
                                  trunc_list=t_ax[hi:lo:-1]).run()
        return CryoscopeReconstruction(inversion="response", qubit=q,
                                       tau=tau_c).reconstruct(res)

    fr = _quiet(_run)
    return (np.asarray(fr.t_list, dtype=float),
            np.asarray(fr.samples, dtype=float) / TRUE_AMP)


def measure_public_ideal(step_amp):
    """Cryoscope measurement with NO control_line -> true s(t) == 1 exactly.

    The control experiment for panel (b): anything that comes back other than a
    flat 1.0 is protocol + reconstruction artifact, independent of any
    distortion. Verified identical to passing a control_line with amplitude=0.
    """
    from sqc.calibration.waveform import WaveformCalibration

    cal = WaveformCalibration(
        qubit=C.make_sqc_qubit(flux=CRYO_FLUX),
        measurement_protocol=PROTOCOL, method="transfer_function",
        fit_type="single_exp", t_max=T_MAX, step_amplitude=step_amp)
    table = _quiet(cal.calibrate)
    return (np.asarray(table.inputs, dtype=float),
            np.asarray(table.outputs, dtype=float))


def _dom_period(t, r):
    """Dominant oscillation period (ns) of a residual, via rFFT peak."""
    r = np.asarray(r, dtype=float)
    r = r - r.mean()
    freq = np.fft.rfftfreq(len(r), float(t[1] - t[0]))
    power = np.abs(np.fft.rfft(r))
    k = int(np.argmax(power[1:])) + 1
    return float(1.0 / freq[k]) if freq[k] > 0 else np.nan


def analytic_step(t, tau):
    """Ground-truth analytic step response on the measured time axis."""
    from sqc.hardware.distortion import SingleExponentialDistortion

    return np.asarray(SingleExponentialDistortion(
        amplitude=TRUE_AMP, tau=tau).step_response(t), dtype=float)


def compute(recompute=False):
    """Run all cryoscope measurements. Cached to npz (~40 s cold)."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    out = {}

    # -- main: tau=100 through the public path --
    t, meas, fa, ft = measure_public(TRUE_TAU_NS)
    ana = analytic_step(t, TRUE_TAU_NS)
    out.update({"t": t, "meas": meas, "ana": ana,
                "fit_amp": np.array([fa]), "fit_tau": np.array([ft])})
    out["fit_curve"] = 1.0 - fa * np.exp(-np.maximum(t, 0.0) / ft)

    # -- tau sweep: recovery error vs surviving tail fraction --
    sw = {k: [] for k in ("rmse", "amp", "tau", "frac")}
    for tau in TAU_SWEEP:
        ts, ms, a_i, t_i = measure_public(tau)
        sw["rmse"].append(float(np.sqrt(np.mean(
            (ms - analytic_step(ts, tau)) ** 2))))
        sw["amp"].append(a_i)
        sw["tau"].append(t_i)
        # tail amplitude still alive at the window start
        sw["frac"].append(float(np.exp(-ts[0] / tau)))
    out["sweep_tau_true"] = np.array(TAU_SWEEP, dtype=float)
    for k, v in sw.items():
        out[f"sweep_{k}"] = np.array(v, dtype=float)

    # -- protocol artifact: measure an IDEAL (undistorted) line --
    # With no control_line the true step response is exactly s(t)=1, so whatever
    # comes back is pure protocol/reconstruction artifact. This is the control
    # experiment that shows panel (b)'s residual is NOT a distortion-measurement
    # error and NOT noise.
    t_id, s_id = measure_public_ideal(TRUE_AMP)
    out["ideal_t"] = t_id
    out["ideal_dev"] = s_id - 1.0

    # artifact vs step_amplitude: DC bias is ~amplitude-independent while the
    # oscillation grows and its period collapses (70 -> 5 ns) -> phase wrapping.
    sc = {k: [] for k in ("rms", "dc", "osc", "period")}
    for sa in STEP_AMP_SCAN:
        _t, _s = measure_public_ideal(sa)
        dev = _s - 1.0
        sc["rms"].append(float(np.sqrt(np.mean(dev ** 2))))
        sc["dc"].append(float(dev.mean()))
        sc["osc"].append(float(dev.std()))
        sc["period"].append(_dom_period(_t, dev))
    out["scan_step_amp"] = np.array(STEP_AMP_SCAN, dtype=float)
    for k, v in sc.items():
        out[f"scan_{k}"] = np.array(v, dtype=float)

    # Falsifiable check on the phase-winding mechanism: the artifact period
    # should equal 2*pi/dw(a). Computed analytically here so the script's own
    # output can be compared against scan_period. Also record the phase
    # accumulated over the cryoscope delay, whose crossing of pi separates the
    # DC-dominated from the oscillation-dominated regime.
    from sqc.config import CONFIG as _CFG
    tau_c = float(_CFG.reconstruction.cryoscope_tau)
    w0 = float(C.make_sqc_qubit(flux=CRYO_FLUX).frequency)
    dw = np.array([float(C.make_sqc_qubit(flux=a).frequency) - w0
                   for a in STEP_AMP_SCAN])
    out["scan_dw"] = dw
    out["scan_period_pred"] = 2.0 * np.pi / np.abs(dw)
    out["scan_phase_over_pi"] = np.abs(dw) * tau_c / np.pi
    # the FFT cannot resolve a period longer than the window span
    out["window_span_ns"] = np.array([float(t[-1] - t[0])])

    # -- window validity --
    wi = {k: [] for k in ("hi_ns", "rmse", "smin", "smax")}
    for hi, lo in WINDOW_SLICES:
        tw, sw_ = measure_window(TRUE_TAU_NS, hi, lo)
        wi["hi_ns"].append(float(tw[-1]))
        wi["rmse"].append(float(np.sqrt(np.mean(
            (sw_ - analytic_step(tw, TRUE_TAU_NS)) ** 2))))
        wi["smin"].append(float(sw_.min()))
        wi["smax"].append(float(sw_.max()))
    for k, v in wi.items():
        out[f"win_{k}"] = np.array(v, dtype=float)

    C.save_npz(SUBDIR, _CACHE, **out)
    return out


def plot(d, out):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(12.8, 8.8))
    ((a, b), (c, e)) = ax

    # -- (a) quantum-measured vs analytic step, tau=100 --
    t = d["t"]
    a.plot(t, d["ana"], **STYLE_ANA, lw=1.6, label="analytic (true)")
    a.plot(t, d["meas"], **STYLE_MEAS, ms=5, lw=1.2, ls="-",
           label="cryoscope-measured")
    a.plot(t, d["fit_curve"], **STYLE_FIT, lw=1.4,
           label=f"fit: A={d['fit_amp'][0]:.4f}, "
                 f"$\\tau$={d['fit_tau'][0]:.1f} ns")
    a.set_xlabel("t (ns)", fontsize=11)
    a.set_ylabel("Normalised step response $s(t)$", fontsize=11)
    a.set_title(f"(a) Cryoscope vs analytic, true A={TRUE_AMP}, "
                f"$\\tau$={TRUE_TAU_NS:.0f} ns", fontsize=11)
    a.legend(fontsize=9, loc="lower right")
    a.grid(True, ls=":", alpha=0.4)

    # -- (b) the residual IS the protocol artifact (control experiment) --
    res = d["meas"] - d["ana"]
    rms = float(np.sqrt(np.mean(res ** 2)))
    ideal = d["ideal_dev"]
    ideal_rms = float(np.sqrt(np.mean(ideal ** 2)))
    b.axhline(0.0, color="k", lw=0.8, alpha=0.5)
    b.plot(t, res, color="#377eb8", lw=1.4,
           label=f"residual, distorted line (rms {rms:.1e})")
    b.plot(d["ideal_t"], ideal, color="#984ea3", ls="--", lw=1.4,
           label=f"IDEAL line, no distortion (rms {ideal_rms:.1e})")
    b.set_xlabel("t (ns)", fontsize=11)
    b.set_ylabel("measured $-$ true", fontsize=11)
    b.set_title("(b) The residual is a PROTOCOL ARTIFACT, not noise:\n"
                "an undistorted line shows the same structure",
                fontsize=11)
    b.legend(fontsize=8, loc="lower left")
    b.grid(True, ls=":", alpha=0.4)

    # inset: artifact vs step_amplitude -> the operating range
    # NOTE: the step_amplitude scan lives in panel (c), not as an inset here.
    # Both traces in this panel sweep the full vertical range, so any inset
    # placed inside these axes covers the curves it is meant to annotate.

    # -- (c) window validity --
    hi, wr = d["win_hi_ns"], d["win_rmse"]
    cols = ["#4daf4a" if v < 1e-2 else "#e41a1c" for v in wr]
    c.bar([f"{TRUE_TAU_NS*0+20.5:.1f}–{h:.0f}" for h in hi], wr,
          color=cols, edgecolor="k", lw=0.6, width=0.55)
    c.set_yscale("log")
    c.set_xlabel("Reconstruction window (ns)", fontsize=11)
    c.set_ylabel("rmse vs analytic", fontsize=11)
    c.set_title("(c) Operating range: window validity (bars) and\n"
                "step-amplitude artifact (inset)", fontsize=11)
    for i, (h, v) in enumerate(zip(hi, wr)):
        c.text(i, v * 1.4, f"{v:.1e}\nmax s={d['win_smax'][i]:.2f}",
               ha="center", fontsize=8)
    # headroom raised to 3e2 so the inset clears the tallest bar's label
    c.set_ylim(1e-3, 3e2)
    c.grid(True, axis="y", which="both", ls=":", alpha=0.4)

    # The other half of the operating envelope: artifact vs step_amplitude.
    # Lives here because panel (b)'s curves leave no free space.
    # raised to 0.70: lower placement put the inset's x-label on top of the
    # tallest bar's "max s=" annotation
    ci = c.inset_axes([0.55, 0.70, 0.42, 0.26])
    ci.loglog(d["scan_step_amp"], d["scan_rms"], "o-", color="#984ea3",
              ms=4, lw=1.2)
    ci.axvline(TRUE_AMP, color="#377eb8", ls="--", lw=0.9)
    ci.text(TRUE_AMP * 0.92, 3e-1, "used here", fontsize=6, rotation=90,
            ha="right", va="center", color="#377eb8")
    ci.set_ylabel("artifact rms", fontsize=6.5, labelpad=1)
    ci.set_xlabel("step amp", fontsize=6.5, labelpad=0)
    ci.tick_params(labelsize=6)
    ci.grid(True, which="both", ls=":", alpha=0.4)

    # -- (d) tau sweep: two orthogonal error sources --
    frac = d["sweep_frac"]
    amp_err = 100.0 * (d["sweep_amp"] / TRUE_AMP - 1.0)
    tau_err = 100.0 * (d["sweep_tau"] / d["sweep_tau_true"] - 1.0)
    e.plot(frac, amp_err, "o-", color="#e41a1c", ms=7, lw=1.5,
           label="fitted $A$ error")
    e.plot(frac, tau_err, "s-", color="#984ea3", ms=7, lw=1.5,
           label=r"fitted $\tau$ error")
    e.axhline(0.0, color="k", lw=0.8, alpha=0.5)
    e2 = e.twinx()
    e2.plot(frac, d["sweep_rmse"], "^--", color="#377eb8", ms=7, lw=1.5,
            label="rmse vs analytic (right axis)")
    e2.set_ylim(0.0, 0.006)
    e2.set_ylabel("rmse vs analytic", fontsize=11, color="#377eb8")
    e2.tick_params(axis="y", labelcolor="#377eb8")
    # tau labels along the TOP of the axes: at y=0 they sat on the A-error curve
    for f, tt in zip(frac, d["sweep_tau_true"]):
        e.annotate(rf"$\tau$={tt:.0f}", (f, 1.0),
                   xycoords=("data", "axes fraction"), xytext=(0, -13),
                   textcoords="offset points", ha="center", fontsize=8,
                   color="dimgray")
    # widen x a little so the leftmost tau label is not clipped by the spine
    e.set_xlim(frac.min() - 0.045, frac.max() + 0.035)
    # headroom for the tau labels along the top edge (the tau=20 point sits at
    # +36.8%, which would otherwise collide with its own label)
    _lo = min(amp_err.min(), tau_err.min()) - 6.0
    _hi = max(amp_err.max(), tau_err.max()) + 12.0
    e.set_ylim(_lo, _hi)
    e.set_xlabel(r"Tail fraction surviving at window start, "
                 r"$e^{-t_0/\tau}$", fontsize=11)
    e.set_ylabel("Parameter recovery error (%)", fontsize=11)
    e.set_title("(d) Two orthogonal error sources: rmse is flat (artifact-\n"
                "dominated), recoverability tracks the surviving tail",
                fontsize=11)
    h1, l1 = e.get_legend_handles_labels()
    h2, l2 = e2.get_legend_handles_labels()
    e.legend(h1 + h2, l1 + l2, fontsize=8.5, loc="lower right")
    e.grid(True, ls=":", alpha=0.4)

    fig.suptitle("P5 — Protocol-driven (cryoscope) vs analytic "
                 "transfer-function measurement", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(os.path.join(out, f"{_CACHE}.png"), dpi=300)
    fig.savefig(os.path.join(out, f"{_CACHE}.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    d = compute(recompute="--recompute" in sys.argv)

    res = d["meas"] - d["ana"]
    ideal = d["ideal_dev"]
    print(f"  window      {d['t'][0]:.1f}..{d['t'][-1]:.1f} ns "
          f"({len(d['t'])} pts)")
    print(f"  ARTIFACT (ideal line, true s=1): rms "
          f"{np.sqrt(np.mean(ideal ** 2)):.3e}  "
          f"DC {ideal.mean():+.3e}  osc_std {ideal.std():.3e}  "
          f"period {_dom_period(d['ideal_t'], ideal):.1f} ns")
    span = float(d["window_span_ns"][0])
    print(f"  artifact vs step_amplitude  (window span {span:.1f} ns; a period "
          f"longer than that is NOT resolvable — shown as '-'):")
    print("     step_amp     rms         DC         osc      phase/pi  "
          "period_pred  period_obs")
    for i, sa in enumerate(d["scan_step_amp"]):
        pred = d["scan_period_pred"][i]
        resolvable = pred <= span
        obs = f"{d['scan_period'][i]:9.1f}" if resolvable else "        -"
        print(f"     {sa:8.3f}  {d['scan_rms'][i]:.3e}  "
              f"{d['scan_dc'][i]:+.3e}  {d['scan_osc'][i]:.3e}  "
              f"{d['scan_phase_over_pi'][i]:7.3f}  {pred:10.1f}  {obs}")
    print(f"  main tau={TRUE_TAU_NS:.0f}  rmse vs analytic "
          f"{np.sqrt(np.mean(res ** 2)):.3e}   "
          f"fit A={d['fit_amp'][0]:.5f} ({100*(d['fit_amp'][0]/TRUE_AMP-1):+.1f}%)"
          f"  tau={d['fit_tau'][0]:.2f} "
          f"({100*(d['fit_tau'][0]/TRUE_TAU_NS-1):+.1f}%)")
    print("  tau sweep:")
    for i, tt in enumerate(d["sweep_tau_true"]):
        print(f"     tau={tt:6.0f}  tail_frac={d['sweep_frac'][i]:.3f}  "
              f"rmse={d['sweep_rmse'][i]:.3e}  "
              f"A={d['sweep_amp'][i]:.5f} "
              f"({100*(d['sweep_amp'][i]/TRUE_AMP-1):+6.1f}%)  "
              f"tau_fit={d['sweep_tau'][i]:7.2f} "
              f"({100*(d['sweep_tau'][i]/tt-1):+6.1f}%)")
    print("  window validity:")
    for i, h in enumerate(d["win_hi_ns"]):
        print(f"     20.5..{h:5.0f} ns  rmse={d['win_rmse'][i]:.3e}  "
              f"s range {d['win_smin'][i]:.4f}..{d['win_smax'][i]:.4f}")

    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
