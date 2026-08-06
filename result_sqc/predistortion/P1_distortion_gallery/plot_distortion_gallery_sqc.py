#!/usr/bin/env python3
"""P1 — Distortion-model gallery: step + frequency response (sqc, analytic).

sqc source:
    sqc.hardware.distortion.{SingleExponentialDistortion,
    MultiExponentialDistortion, CascadeDistortion, FIRDistortion,
    IIRDistortion} — each exposes step_response(t) and
    frequency_response(omega). CustomTransferDistortion is omitted: it needs a
    user-supplied H(omega) grid, so it has no representative parameter set.

Physics story (paper role): control-line distortion model overview — the model
zoo behind predistortion. Six representative models spanning the exponential
family (direct-path vs pure-lowpass, single vs multi vs cascaded stages) and
the filter family (AWG reconstruction FIR, cryo-line bandwidth IIR).

UNITS PITFALL (verified, do not "fix"): the exponential family defines
frequency_response(omega) in PHYSICAL rad/ns, but FIRDistortion /
IIRDistortion define it on the NORMALISED grid dt=1 (polyval on exp(-1j*omega),
see distortion.py:597 / :745, both documented "nominal dt=1"). To put all six on
one Bode axis the filter family must be evaluated at omega*DT. Cross-check:
the 200 MHz first-order IIR then lands at |H|=0.695 (-3.15 dB) at f=0.2 GHz,
i.e. its design corner — correct.

Produces:
    result_sqc/predistortion/P1_distortion_gallery/distortion_gallery_sqc.png / .pdf
        (a) step response, full 500 ns window
        (b) step response, first 30 ns (direct-path jump vs smooth ramp)
        (c) Bode magnitude (dB, log-f)
        (d) Bode phase (deg, log-f)
    result_sqc/predistortion/P1_distortion_gallery/distortion_gallery_sqc.npz

Cost: LOW — fully analytic, seconds. Cached to npz anyway (block convention);
pass --recompute to force a rebuild.

Usage:
    python result_sqc/predistortion/P1_distortion_gallery/plot_distortion_gallery_sqc.py
    python result_sqc/predistortion/P1_distortion_gallery/plot_distortion_gallery_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "predistortion/P1_distortion_gallery"
_CACHE = "distortion_gallery_sqc"

# --- grids (R9: time axes derive from CONFIG.awg.dt via np.arange) -------
from sqc.config import CONFIG  # noqa: E402

DT = CONFIG.awg.dt                 # 0.5 ns
T_END = 500.0                      # ns — plotted window (panel a)
TAIL_FLOOR = 1e-5                  # panel (b) log-y floor (echo FIR hits 0)
# Settling time is a METRIC, not a plotted curve, and the slowest tail
# (multi_exp, tau=400 ns, A=0.02) needs 400*ln(0.02/1e-3) ~ 1200 ns to reach
# the 1e-3 band. Evaluate it on its own long grid so T_END can stay short
# enough that panel (a) still resolves the leading-edge structure.
T_SETTLE_END = 3000.0              # ns
# frequency grid: 1e-4..1 GHz. NOT a time axis, so linspace/logspace is fine
# (R9 constrains time axes only). Upper end = Nyquist 1/(2*DT) = 1 GHz.
# 2000 pts (not 600): the echo FIR ripples with period 1/FIR_ECHO_DELAY_NS =
# 0.25 GHz, and a coarser log grid under-samples that near 1 GHz -> visible
# aliasing jaggies in panel (c).
F_GHZ = np.logspace(-4, 0, 2000)

# --- model zoo ----------------------------------------------------------
# Exponential family: physical rad/ns frequency response.
EXP_FAMILY = ("single_exp", "single_exp_smooth", "multi_exp", "cascade")
# Filter family: frequency_response on the NORMALISED dt=1 grid (see module
# docstring) -> must be evaluated at omega*DT.
FILTER_FAMILY = ("fir_awg", "iir_cryo")

# FIR: single reflection echo on the control line (impedance mismatch), the
# textbook FIR-shaped distortion. A centred Gaussian was tried first and
# rejected: truncating it leaves -35 dB sidelobes that shred the Bode phase,
# and its group delay swamps the exponential family on the shared phase axis.
FIR_ECHO_DELAY_NS = 4.0            # round-trip delay of the reflection
FIR_ECHO_AMP = 0.05                # reflection coefficient
IIR_FC_GHZ = 0.2                   # cryo-line first-order corner (200 MHz)

LABELS = {
    "single_exp": r"single_exp  $A$=0.05, $\tau$=100 ns",
    "single_exp_smooth": r"single_exp smooth  $\tau$=100 ns",
    "multi_exp": r"multi_exp  (0.04, 80 ns)+(0.02, 400 ns)",
    "cascade": r"cascade  (0.05, 50 ns)$\rightarrow$(0.03, 200 ns)",
    "fir_awg": rf"FIR  echo {FIR_ECHO_AMP:.2f} @ {FIR_ECHO_DELAY_NS:.0f} ns",
    "iir_cryo": rf"IIR  1st-order lowpass $f_c$={IIR_FC_GHZ*1e3:.0f} MHz",
}
STYLE = {
    "single_exp":        {"color": "#377eb8", "ls": "-"},
    "single_exp_smooth": {"color": "#80b1d3", "ls": "--"},
    "multi_exp":         {"color": "#e41a1c", "ls": "-"},
    "cascade":           {"color": "#4daf4a", "ls": "-"},
    "fir_awg":           {"color": "#984ea3", "ls": "-."},
    "iir_cryo":          {"color": "#ff7f00", "ls": ":"},
}
MODELS = EXP_FAMILY + FILTER_FAMILY


def build_models():
    """The six representative distortion models, keyed by short name.

    Parameter choices are the canonical ones already exercised by the phase-4
    unit tests (tests/unit/test_distortion_models.py,
    tests/unit/test_cascade_distortion.py) and by P2's validation harness
    (A=0.05, tau=100 ns), so the gallery and the correction figures show the
    same physics.
    """
    from scipy.signal import bilinear
    from sqc.hardware.distortion import (
        SingleExponentialDistortion, MultiExponentialDistortion,
        CascadeDistortion, FIRDistortion, IIRDistortion,
    )

    # FIR: causal reflection echo — main tap plus one delayed tap.
    taps = np.zeros(int(round(FIR_ECHO_DELAY_NS / DT)) + 1)
    taps[0] = 1.0
    taps[-1] = FIR_ECHO_AMP
    taps /= taps.sum()             # unit DC gain, like every other model here

    # IIR: bilinear-discretised first-order lowpass at IIR_FC_GHZ.
    tau_c = 1.0 / (2.0 * np.pi * IIR_FC_GHZ)   # rad/ns -> ns
    b, a = bilinear([0.0, 1.0], [tau_c, 1.0], fs=1.0 / DT)

    return {
        "single_exp": SingleExponentialDistortion(amplitude=0.05, tau=100.0),
        "single_exp_smooth": SingleExponentialDistortion(tau=100.0, smooth=True),
        "multi_exp": MultiExponentialDistortion(
            amplitudes=np.array([0.04, 0.02]), taus=np.array([80.0, 400.0])),
        "cascade": CascadeDistortion(stages=[
            SingleExponentialDistortion(amplitude=0.05, tau=50.0),
            SingleExponentialDistortion(amplitude=0.03, tau=200.0)]),
        "fir_awg": FIRDistortion(taps=taps),
        "iir_cryo": IIRDistortion(b_coeffs=b, a_coeffs=a),
    }


def settling_time(t, step, tol=1e-3):
    """First time after which |s(t) - s(inf)| stays below tol (ns).

    s(inf) is taken as the analytic DC gain 1.0 (all six models are unit-DC by
    construction). Returns nan if the tail never settles inside the window.
    """
    bad = np.abs(np.asarray(step) - 1.0) > tol
    idx = np.nonzero(bad)[0]
    return float(t[idx[-1]]) if len(idx) and idx[-1] + 1 < len(t) else (
        np.nan if len(idx) else 0.0)


def compute(recompute=False):
    """Evaluate step + frequency response for all six models. Cached to npz."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    models = build_models()
    t = np.arange(0.0, T_END, DT)
    t_long = np.arange(0.0, T_SETTLE_END, DT)   # settling metric only
    omega = 2.0 * np.pi * F_GHZ            # physical rad/ns

    out = {"t": t, "f_ghz": F_GHZ, "dt": np.array([DT])}
    for name in MODELS:
        m = models[name]
        step = np.asarray(m.step_response(t), dtype=float)
        step_long = np.asarray(m.step_response(t_long), dtype=float)
        # UNITS: filter family lives on the normalised dt=1 grid -> omega*DT.
        w = omega * DT if name in FILTER_FAMILY else omega
        H = np.asarray(m.frequency_response(w), dtype=complex)
        out[f"{name}_step"] = step
        out[f"{name}_H_mag"] = np.abs(H)
        out[f"{name}_H_phase_deg"] = np.degrees(np.unwrap(np.angle(H)))
        out[f"{name}_settling"] = np.array([settling_time(t_long, step_long)])

    C.save_npz(SUBDIR, _CACHE, **out)
    return out


_MAG_FLOOR = 1e-4   # clip: the bilinear IIR has |H| -> 0 exactly at Nyquist


def plot(d, out):
    import matplotlib.pyplot as plt

    t, f = d["t"], d["f_ghz"]
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.2))
    (ax_a, ax_b), (ax_c, ax_e) = axes

    for name in MODELS:
        st, lb = STYLE[name], LABELS[name]
        step = d[f"{name}_step"]
        ax_a.plot(t, step, lw=1.5, label=lb, **st)
        # tail deviation on log-y: an exponential tail is a STRAIGHT LINE of
        # slope -1/(tau*ln10), so the panel reads off tau directly.
        ax_b.semilogy(t, np.maximum(np.abs(step - 1.0), TAIL_FLOOR),
                      lw=1.5, **st)
        mag = np.maximum(d[f"{name}_H_mag"], _MAG_FLOOR)
        ax_c.semilogx(f, 20.0 * np.log10(mag), lw=1.5, **st)
        ax_e.semilogx(f, d[f"{name}_H_phase_deg"], lw=1.5, **st)

    # -- (a) full-window step --
    ax_a.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax_a.set_xlabel("Time (ns)", fontsize=11)
    ax_a.set_ylabel("Step response $s(t)$", fontsize=11)
    ax_a.set_title("(a) Step response, full window", fontsize=11)
    ax_a.set_xlim(0, T_END)
    ax_a.legend(fontsize=8, loc="lower right", framealpha=0.9)

    # -- (b) tail deviation, log-y (separates what (a) cannot resolve) --
    ax_b.axhline(1e-3, color="gray", lw=0.9, ls=":", alpha=0.9)
    ax_b.text(T_END * 0.03, 1.3e-3, "$10^{-3}$ settling band",
              fontsize=8, color="gray")
    ax_b.set_xlabel("Time (ns)", fontsize=11)
    ax_b.set_ylabel(r"Tail deviation $|1 - s(t)|$", fontsize=11)
    ax_b.set_title(r"(b) Tail decay (log-y: slope $=-1/\tau$)", fontsize=11)
    ax_b.set_xlim(0, T_END)
    ax_b.set_ylim(TAIL_FLOOR * 0.7, 2.0)

    # -- (c) Bode magnitude, plus an inset for the direct-path shelf --
    # The four direct-path models never fall below ~-0.9 dB, so on the axis
    # range needed by smooth/IIR (-45 dB) their (1-A) shelf is invisible.
    # Inset x is LINEAR (not log): the echo FIR combs with period
    # 1/FIR_ECHO_DELAY_NS, which a log axis crushes into the last half-decade.
    ax_ci = ax_c.inset_axes([0.11, 0.13, 0.46, 0.40])
    for name in ("single_exp", "multi_exp", "cascade", "fir_awg"):
        ax_ci.plot(f, 20.0 * np.log10(d[f"{name}_H_mag"]),
                   lw=1.3, **STYLE[name])
    ax_ci.set_xlim(0, f[-1])
    ax_ci.set_ylim(-1.0, 0.12)
    ax_ci.set_title(r"direct-path shelf $20\log_{10}(1-A)$, linear $f$",
                    fontsize=7.5, pad=2)
    ax_ci.set_xlabel("GHz", fontsize=7, labelpad=1)
    ax_ci.tick_params(labelsize=7)
    ax_ci.grid(True, which="both", ls=":", alpha=0.4)

    ax_c.axhline(0.0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax_c.axhline(-3.0, color="gray", lw=0.8, ls=":", alpha=0.8)
    ax_c.text(f[0] * 1.4, -3.0, "-3 dB", fontsize=8, color="gray",
              va="bottom")
    ax_c.set_xlabel("Frequency (GHz)", fontsize=11)
    ax_c.set_ylabel(r"$|H(\omega)|$ (dB)", fontsize=11)
    ax_c.set_title("(c) Bode magnitude", fontsize=11)
    ax_c.set_ylim(-45, 5)

    # -- (d) Bode phase --
    ax_e.axhline(0.0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax_e.set_xlabel("Frequency (GHz)", fontsize=11)
    ax_e.set_ylabel(r"$\arg H(\omega)$ (deg)", fontsize=11)
    ax_e.set_title("(d) Bode phase", fontsize=11)

    for ax in (ax_a, ax_b, ax_c, ax_e):
        ax.grid(True, which="both", ls=":", alpha=0.4)

    fig.suptitle("Control-line distortion model gallery (sqc, analytic)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(out, f"{_CACHE}.png"), dpi=300)
    fig.savefig(os.path.join(out, f"{_CACHE}.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    d = compute(recompute="--recompute" in sys.argv)

    f = d["f_ghz"]
    for name in MODELS:
        step, mag = d[f"{name}_step"], d[f"{name}_H_mag"]
        # -3 dB corner: first crossing of |H| = 1/sqrt(2)
        below = np.nonzero(mag < 2.0 ** -0.5)[0]
        fc = f"{f[below[0]] * 1e3:8.1f} MHz" if len(below) else "      none"
        print(f"  {name:18s} s(0)={step[0]:.4f}  s(end)={step[-1]:.6f}  "
              f"settling={d[f'{name}_settling'][0]:6.1f} ns  f_-3dB={fc}")

    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
