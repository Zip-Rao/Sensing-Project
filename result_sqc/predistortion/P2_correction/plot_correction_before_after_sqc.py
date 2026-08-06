#!/usr/bin/env python3
"""P2 — Predistortion before/after, matched vs mismatched model (sqc).

This block has NO frozen `result/` counterpart — closed-loop predistortion is a
new sqc capability. Figures demonstrate sqc's own correction quality, not a
per-point overlay against src.

sqc source:
    sqc.workflows.PredistortionValidationWorkflow — inject known distortion →
    measure step response → fit H(omega) → design inverse → apply → verify.
    run() returns {target, on_chip_uncorrected, on_chip_corrected,
    awg_predistorted, inverse_model, measured_model, metrics}.

TWO SCENARIOS, and the reason there are two (important):
    WaveformCalibration._fit_step_response really runs curve_fit
    (waveform.py:203) — it does NOT read the injected parameters back. But when
    the fit model CLASS equals the true model CLASS, fitting a single_exp curve
    with a single_exp model is exact to machine precision, so the designed
    inverse cancels the distortion algebraically and improvement_factor comes
    out ~3e12. That number is real but it is a TAUTOLOGY: it only proves the
    exact inverse of a known model cancels that model. Shown alone it would be
    a misleading figure.
    So P2 shows both:
      "matched"  single_exp truth, single_exp fit  -> algebraic ceiling (~3e12)
      "mismatch" 2-component multi_exp truth, 3-component multi_exp fit
                 -> finite, honest improvement (~1e2) with a visible residual
    The mismatch case is the one that describes real hardware: several tail
    time constants, a finite-order calibration model.

Produces:
    result_sqc/predistortion/P2_correction/predistortion_correction_sqc.png / .pdf
        (a) matched  : target / on-chip uncorrected / on-chip corrected
        (b) mismatch : same three (+ flat-top zoom insets)
        (c) residual |on-chip - target|, log-y, all four traces
        (d) AWG waveform actually played (over/undershoot that buys the fix)
    result_sqc/predistortion/P2_correction/predistortion_metrics_sqc.npz

Cost: LOW — analytical measurement path (no quantum simulation), seconds.
Cached to npz; pass --recompute to force a rebuild.

Usage:
    python result_sqc/predistortion/P2_correction/plot_correction_before_after_sqc.py
    python result_sqc/predistortion/P2_correction/plot_correction_before_after_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

from sqc.config import CONFIG  # noqa: E402

SUBDIR = "predistortion/P2_correction"
_CACHE = "predistortion_metrics_sqc"

# --- target waveform (R9: time axis derives from CONFIG.awg.dt) ----------
DT = CONFIG.awg.dt                 # 0.5 ns
T_END = 500.0                      # ns — long enough to expose the 400 ns tail
FLAT_START, FLAT_STOP = 50.0, 250.0   # flat-top edges (ns)
# A flat-top, not a bare step: the rising edge shows the slow tail creep and
# the falling edge shows the mirror-image undershoot. Both are what a flux
# pulse in a real sensing sequence actually suffers.

# Ground-truth distortion per scenario. "matched" reuses P1's canonical
# single_exp (A=0.05, tau=100 ns) so the gallery and this figure share physics.
SCENARIOS = ("matched", "mismatch")
MATCHED_AMP, MATCHED_TAU = 0.05, 100.0
MISMATCH_AMPS = (0.04, 0.02)
MISMATCH_TAUS = (80.0, 400.0)

TITLES = {
    "matched": ("(a) Matched model: single_exp truth, single_exp fit\n"
                "algebraic ceiling"),
    "mismatch": ("(b) Mismatched model: 2-component truth, 3-component fit\n"
                 "realistic device"),
}
CURVE = {          # trace -> style, shared by all panels
    "target": {"color": "k", "ls": "-", "lw": 1.6},
    "unc":    {"color": "#e41a1c", "ls": "-", "lw": 1.4},
    "cor":    {"color": "#377eb8", "ls": "--", "lw": 1.5},
}
SCEN_COLOR = {"matched": "#377eb8", "mismatch": "#ff7f00"}
ZOOM = (FLAT_START - 5.0, FLAT_START + 95.0)   # inset window on the rise


def make_target():
    """The desired on-chip flat-top waveform."""
    from sqc.control.waveform import Waveform

    t = np.arange(0.0, T_END, DT)
    samples = np.where((t > FLAT_START) & (t < FLAT_STOP), 1.0, 0.0)
    return Waveform(t_list=t, samples=samples)


def true_distortion(scenario):
    """Ground-truth distortion injected into the control line."""
    from sqc.hardware.distortion import (
        SingleExponentialDistortion, MultiExponentialDistortion,
    )
    if scenario == "matched":
        return SingleExponentialDistortion(amplitude=MATCHED_AMP,
                                           tau=MATCHED_TAU)
    return MultiExponentialDistortion(amplitudes=np.array(MISMATCH_AMPS),
                                      taus=np.array(MISMATCH_TAUS))


def run_validation(scenario):
    """Run the end-to-end predistortion validation. Returns the run() dict.

    The designer is left at its default (PredistortionDesigner(method="auto")),
    and fit_type is inferred by the workflow from the true-distortion class
    (_infer_fit_type): single_exp -> "single_exp", multi_exp -> "multi_exp"
    with n_exp_components=3. That 3-vs-2 component asymmetry is exactly the
    mismatch this figure is built to expose.
    """
    from sqc.workflows import PredistortionValidationWorkflow

    wf = PredistortionValidationWorkflow(
        target_waveform=make_target(),
        true_distortion=true_distortion(scenario),
    )
    return wf.run()


# metric keys copied straight out of workflow metrics (scalars only)
_SCALARS = ("rmse_uncorrected", "rmse_corrected", "improvement_factor",
            "settling_uncorrected_ns", "settling_corrected_ns")


def compute(recompute=False):
    """Run both scenarios, flatten waveforms + metrics into a npz dict."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    out = {"t": make_target().t_list, "dt": np.array([DT])}
    for s in SCENARIOS:
        r = run_validation(s)
        m = r["metrics"]
        out[f"{s}_target"] = r["target"].samples
        out[f"{s}_unc"] = r["on_chip_uncorrected"].samples
        out[f"{s}_cor"] = r["on_chip_corrected"].samples
        out[f"{s}_awg"] = r["awg_predistorted"].samples
        for k in _SCALARS:
            out[f"{s}_{k}"] = np.array([float(m[k])])
        out[f"{s}_inverse_type"] = np.array([m["inverse_model_type"]])
        # fit_params shape differs per scenario (single_exp: amplitude/tau;
        # multi_exp: amplitudes/taus arrays) -> store each field it actually has
        fp = m["calibration_fit_type"]
        for key in ("amplitude", "tau", "amplitudes", "taus"):
            if key in fp:
                out[f"{s}_fit_{key}"] = np.atleast_1d(np.asarray(fp[key],
                                                                 dtype=float))

    C.save_npz(SUBDIR, _CACHE, **out)
    return out


_RES_FLOOR = 1e-16   # log-y clip: matched residual reaches machine precision


def _panel_waveform(ax, d, s):
    """One scenario's target / uncorrected / corrected, with a rise zoom."""
    t = d["t"]
    ax.plot(t, d[f"{s}_target"], label="target (desired on-chip)",
            **CURVE["target"])
    ax.plot(t, d[f"{s}_unc"], label="on-chip, uncorrected", **CURVE["unc"])
    ax.plot(t, d[f"{s}_cor"], label="on-chip, corrected", **CURVE["cor"])
    ax.set_xlabel("Time (ns)", fontsize=11)
    ax.set_ylabel("Flux amplitude (a.u.)", fontsize=11)
    ax.set_title(TITLES[s], fontsize=10.5)
    ax.set_xlim(0, T_END)
    ax.set_ylim(-0.09, 1.14)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.92)
    ax.grid(True, ls=":", alpha=0.4)

    # Inset: the corrected trace sits on top of the target at this scale, so
    # zoom the post-edge flat-top where the uncorrected tail creeps upward.
    # Placed on the right (t>280 is empty): the lower-left would cover the
    # post-falling-edge residual tail, which is part of the story.
    axi = ax.inset_axes([0.56, 0.30, 0.41, 0.33])
    for key in ("target", "unc", "cor"):
        axi.plot(t, d[f"{s}_{key}"], **CURVE[key])
    axi.set_xlim(*ZOOM)
    axi.set_ylim(0.93, 1.008)
    axi.set_title("flat-top creep", fontsize=7.5, pad=2)
    axi.tick_params(labelsize=7)
    axi.grid(True, ls=":", alpha=0.4)


def plot(d, out):
    import matplotlib.pyplot as plt

    t = d["t"]
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.6))
    (ax_a, ax_b), (ax_c, ax_e) = axes

    _panel_waveform(ax_a, d, "matched")
    _panel_waveform(ax_b, d, "mismatch")

    # -- (c) residual, log-y: where correction quality actually shows --
    for s in SCENARIOS:
        col = SCEN_COLOR[s]
        for key, ls, tag in (("unc", "-", "uncorrected"),
                             ("cor", "--", "corrected")):
            res = np.maximum(np.abs(d[f"{s}_{key}"] - d[f"{s}_target"]),
                             _RES_FLOOR)
            ax_c.semilogy(t, res, ls=ls, lw=1.4, color=col,
                          alpha=1.0 if key == "cor" else 0.55,
                          label=f"{s}, {tag}")
    ax_c.axhline(1e-3, color="gray", lw=0.9, ls=":", alpha=0.9)
    ax_c.text(T_END * 0.03, 1.4e-3, "$10^{-3}$ settling tolerance",
              fontsize=8, color="gray")
    ax_c.set_xlabel("Time (ns)", fontsize=11)
    ax_c.set_ylabel(r"Residual $|$on-chip $-$ target$|$", fontsize=11)
    ax_c.set_title("(c) Residual (log-y): matched hits machine precision,\n"
                   "mismatch leaves a real floor", fontsize=10.5)
    ax_c.set_xlim(0, T_END)
    ax_c.set_ylim(1e-17, 3.0)
    ax_c.legend(fontsize=8, loc="center right", framealpha=0.92, ncol=1)
    ax_c.grid(True, which="both", ls=":", alpha=0.4)

    # -- (d) what the AWG actually plays --
    ax_e.plot(t, d["matched_target"], label="target", **CURVE["target"])
    for s in SCENARIOS:
        ax_e.plot(t, d[f"{s}_awg"], lw=1.4, color=SCEN_COLOR[s],
                  ls="-" if s == "matched" else "--",
                  label=f"AWG predistorted ({s})")
    ax_e.set_xlabel("Time (ns)", fontsize=11)
    ax_e.set_ylabel("AWG amplitude (a.u.)", fontsize=11)
    ax_e.set_title("(d) Waveform played by the AWG:\n"
                   "edge overshoot + tail droop is the price of a flat chip",
                   fontsize=10.5)
    ax_e.set_xlim(0, T_END)
    ax_e.legend(fontsize=8, loc="upper right", framealpha=0.92)
    ax_e.grid(True, ls=":", alpha=0.4)

    fig.suptitle("Waveform predistortion: before / after (sqc)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(os.path.join(out, "predistortion_correction_sqc.png"), dpi=300)
    fig.savefig(os.path.join(out, "predistortion_correction_sqc.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    d = compute(recompute="--recompute" in sys.argv)

    for s in SCENARIOS:
        res_c = np.abs(d[f"{s}_cor"] - d[f"{s}_target"]).max()
        res_u = np.abs(d[f"{s}_unc"] - d[f"{s}_target"]).max()
        print(f"  {s:9s} rmse {d[f'{s}_rmse_uncorrected'][0]:.3e} -> "
              f"{d[f'{s}_rmse_corrected'][0]:.3e}  "
              f"x{d[f'{s}_improvement_factor'][0]:.3g}   "
              f"max|e| {res_u:.2e} -> {res_c:.2e}   "
              f"settling {d[f'{s}_settling_uncorrected_ns'][0]:.0f} -> "
              f"{d[f'{s}_settling_corrected_ns'][0]:.0f} ns   "
              f"inv={d[f'{s}_inverse_type'][0]}")
        awg = d[f"{s}_awg"]
        print(f"            AWG range {awg.min():+.4f} .. {awg.max():+.4f}")

    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
