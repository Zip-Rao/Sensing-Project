#!/usr/bin/env python3
"""F4 ⭐ — Hybrid transient→Ramsey closed-loop speed-up (sqc).

This block has NO frozen `result/` counterpart — closed-loop frequency
calibration is a new sqc capability. Figures here demonstrate sqc's own
physical correctness, not a per-point overlay against src.

sqc source:
    sqc.workflows.FrequencyCalibrationWorkflow — default 2-stage hybrid
    (coarse transient+gradient, order=3  →  fine Ramsey+gradient), see run():
        history[i] = {phase, global_iter, cost, residual, V, ...}
        switch_iter, stage_boundaries, metrics{total_cost, coarse_iter,
        fine_iter, cost_transient_per_meas, cost_ramsey_per_meas}

Physics story (⭐ paper role): a cheap-but-approximate probe does the long
haul, then an expensive-but-exact probe finishes. The coarse transient phase
costs 2 mesolve per measurement and drives the residual down to the switch
threshold; the fine Ramsey phase costs 2·len(tau_list) = 800 per measurement
(400×) but converges to tight tolerance. Plotting residual against *cumulative
cost* — not iteration — is what exposes the speed-up: reaching the switch
point via Ramsey alone would have cost orders of magnitude more.

IMPORTANT (flux convention): V_seed / V_a / V_b are flux OFFSETS relative to
the qubit construction point (OPTIMAL_FLUX), not absolute flux — see
[[sqc-frequency-measurement-api]]. The original stub passed
V_seed=C.OPTIMAL_FLUX, which doubles the bias off-curve and stalls the loop
(the same bug fixed in F3).

Produces:
    result_sqc/frequency_calibration/F4_hybrid_convergence/freq_calibration_convergence_sqc.png / .pdf
        residual |f - f_target| vs iteration, coarse/fine phases shaded,
        switch_iter marked.
    result_sqc/frequency_calibration/F4_hybrid_convergence/freq_calibration_cost_sqc.png / .pdf
        residual vs cumulative measurement cost — the transient→Ramsey
        speed-up story.
    result_sqc/frequency_calibration/F4_hybrid_convergence/freq_calibration_sqc.npz  (raw history)

Cost: HIGH — coarse transient is cheap, but the fine Ramsey phase runs
double-sweep (~8.6 s per measurement). Cached to npz; --recompute to refresh.

Usage:
    python result_sqc/frequency_calibration/F4_hybrid_convergence/plot_hybrid_convergence_sqc.py
    python result_sqc/frequency_calibration/F4_hybrid_convergence/plot_hybrid_convergence_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "frequency_calibration/F4_hybrid_convergence"
_CACHE = "freq_calibration_sqc"
TWO_PI = 2.0 * np.pi

# Target set below the optimal-flux frequency so the loop has work to do.
F_TARGET_OFFSET_MHZ = -20.0  # rad·GHz target = f(OPTIMAL_FLUX) + 2π·(-20e-3)
# Flux OFFSETS (not absolute flux — see module docstring).
V_SEED = 0.0                 # start at the sweet spot
V_BOUNDS = (-0.05, 0.01)     # clamp the search to the monotone, on-curve region

PHASE_COLOR = {"coarse": "#4daf4a", "fine": "#377eb8"}

SWITCH_MHZ = 5.0            # coarse→fine handoff threshold (workflow default)
# Local flux sensitivity |df/dΦ| near the root, measured from the qubit model
# (≈1202 MHz per Φ₀ at offset −0.016). Used to size the fine stage's probe step.
SENSITIVITY_MHZ_PER_PHI0 = 1202.0
# The fine stage's blind first probe must match the scale it *arrives* at, i.e.
# the switch threshold — not the coarse stage's 0.01 Φ₀ default (= 12 MHz here),
# which would kick the search further from the root than it started and burn
# expensive Ramsey iterations recovering. This is what "fine" means.
FINE_FIRST_STEP = SWITCH_MHZ / SENSITIVITY_MHZ_PER_PHI0   # ≈0.0042 Φ₀

# Residual thresholds (MHz) for the cost-to-reach comparison. Total cost-to-
# converge is the WRONG metric (the two runs stop at different final accuracy);
# the fair one is "cost to first reach residual X".
THRESHOLDS_MHZ = np.array([10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.016])


def cost_to_reach(res_mhz, cost, thr):
    """Cumulative cost when |residual| first drops to ≤ thr (None if never)."""
    hit = np.asarray(res_mhz) <= thr
    return float(np.asarray(cost)[np.argmax(hit)]) if hit.any() else np.nan


def _target(q):
    return q.frequency + TWO_PI * (F_TARGET_OFFSET_MHZ * 1e-3)


class CountSolves:
    """Count TRUE solver calls (mesolve + sesolve) across the sqc call graph.

    The workflow's own cost model bills a transient measurement as a flat 2,
    which is badly wrong for the "fit" G₃ route (a moving omega_d invalidates the
    G₃ cache every iteration, forcing an ~80-call recalibration — measured 45×
    undercount). Both modules do `from qutip import mesolve` at import time, so
    the module attribute must be replaced; qutip.sesolve is called qualified, so
    it is patched on qutip itself. kernel_full runs on sesolve, so counting only
    mesolve would flatter the transient route instead.
    """

    def __enter__(self):
        import time
        import qutip
        import sqc.calibration.frequency as F
        import sqc.reconstruction.kernel as K
        self.n = 0
        self.marks = []          # true solver calls per measure() call
        self.time_marks = []     # wall-clock time after each measure() call
        self.t0 = time.perf_counter()
        self._targets = []

        def wrap(mod, name):
            orig = getattr(mod, name)
            self._targets.append((mod, name, orig))

            def counted(*a, **kw):
                self.n += 1
                return orig(*a, **kw)

            setattr(mod, name, counted)

        wrap(F, "mesolve")
        wrap(K, "mesolve")
        wrap(qutip, "sesolve")

        # Snapshot the counter after every measure() so each history row can be
        # placed on a TRUE cumulative-cost axis. Without this the plot would have
        # to use the workflow's billed cost, which undercounts transient routes
        # by 45-82x and would flatter exactly the strategy under test.
        FM = F.FrequencyMeasurement
        orig_measure = FM.measure
        self._targets.append((FM, "measure", orig_measure))
        marks = self.marks
        time_marks = self.time_marks
        t0 = self.t0

        def measure(inner_self, *a, **kw):
            out = orig_measure(inner_self, *a, **kw)
            marks.append(self.n)
            time_marks.append(time.perf_counter() - t0)
            return out

        FM.measure = measure
        return self

    def __exit__(self, *exc):
        for mod, name, orig in self._targets:
            setattr(mod, name, orig)
        return False


def true_residual_mhz(V):
    """Independent check: analytic |f(Φ_opt+V) − f_target| in MHz.

    Guards against a measurement method that converges to a residual it reports
    but has not actually achieved.
    """
    q = C.make_sqc_qubit()
    tgt = _target(q)
    q.change_flux(C.OPTIMAL_FLUX + float(V))
    return abs(q.frequency - tgt) / TWO_PI * 1e3


def run_hybrid():
    """Run the 2-stage transient→Ramsey hybrid. Returns the run() dict.

    Uses explicit stages rather than the flat preset for ONE reason: the preset
    does not forward step tuning to the fine stage, so it inherits
    first_bias_step=0.01 Φ₀ — a 12 MHz blind probe on a stage that arrives with
    <5 MHz residual. The coarse phase's head start is then thrown away and the
    hybrid ends up costing MORE than Ramsey-only. Everything else matches the
    preset (coarse: transient+gradient order=3; fine: Ramsey+gradient).
    """
    from sqc.workflows import FrequencyCalibrationWorkflow
    from sqc.workflows.frequency_calibration import CalibrationStage

    q = C.make_sqc_qubit()
    stages = [
        CalibrationStage(
            name="coarse", measure_method="transient", step_method="gradient",
            epsilon_f=TWO_PI * SWITCH_MHZ * 1e-3, max_iter=15, order=3,
            damping=0.8, first_bias_step=0.01, max_bias_step=0.02,
        ),
        CalibrationStage(
            name="fine", measure_method="ramsey", step_method="gradient",
            epsilon_f=1e-4, max_iter=10,
            damping=0.8, first_bias_step=FINE_FIRST_STEP,
            max_bias_step=0.02,
        ),
    ]
    wf = FrequencyCalibrationWorkflow(
        qubit=q, f_target=_target(q), stages=stages,
        V_seed=V_SEED, V_a=V_BOUNDS[0], V_b=V_BOUNDS[1],
    )
    return wf.run()


def run_ramsey_only():
    """Single-stage Ramsey-only baseline (same target/seed) for cost contrast."""
    from sqc.workflows import FrequencyCalibrationWorkflow
    from sqc.workflows.frequency_calibration import CalibrationStage

    q = C.make_sqc_qubit()
    stages = [CalibrationStage(
        name="ramsey_only", measure_method="ramsey", step_method="gradient",
        epsilon_f=1e-4, max_iter=15,
    )]
    wf = FrequencyCalibrationWorkflow(
        qubit=q, f_target=_target(q), stages=stages,
        V_seed=V_SEED, V_a=V_BOUNDS[0], V_b=V_BOUNDS[1],
    )
    return wf.run()


def run_transient_track():
    """Single-stage transient with drive TRACKING + kernel_full G₃.

    This is the configuration that actually wins (measured: 6 iters, 12 solver
    calls, true residual 0.00497 MHz — verified against the analytic dispersion,
    so it is not the transient being self-consistently overconfident).

    Two changes vs the coarse stage of the hybrid, both necessary:

    - ``drive_policy="track"`` predicts f_q at the point about to be measured, so
      the *measured* detuning stays <1 MHz instead of starting at 20 MHz — i.e.
      deep inside the transient linear window rather than past the 17.6 MHz fold
      that F2 characterises. This is what removes the coarse phase's wandering.
    - ``g3_source="kernel_full"`` is REQUIRED with track. The "fit" route
      (``_fit_g3_at_delta_max``) parks the qubit at the sweet spot, adds an
      artificial ±Δ, and fits an ODD polynomial — which assumes p_diff is
      symmetric about Δ=0. A moved omega_d offsets that scan window, breaking the
      assumption: measured G₁ swings through zero and flips sign, and the loop
      diverges (probed: 16 iters, 1432 solver calls, stalled at 4.07 MHz true
      residual). kernel_full integrates the kernel of the ACTUAL pulse at
      omega_d, so it has no symmetry assumption — and costs ~2 sesolve.
    """
    from sqc.workflows import FrequencyCalibrationWorkflow
    from sqc.workflows.frequency_calibration import CalibrationStage

    q = C.make_sqc_qubit()
    stages = [CalibrationStage(
        name="transient_track", measure_method="transient",
        step_method="gradient", epsilon_f=1e-4, max_iter=15,
        order=3, g3_source="kernel_full",
        drive_policy="track", sensitivity_source="secant",
        damping=0.8, first_bias_step=0.01, max_bias_step=0.02,
    )]
    wf = FrequencyCalibrationWorkflow(
        qubit=q, f_target=_target(q), stages=stages,
        V_seed=V_SEED, V_a=V_BOUNDS[0], V_b=V_BOUNDS[1],
    )
    return wf.run()


def _flatten(res):
    """Merged history dict → flat arrays (residual in MHz, cost, phase)."""
    h = res["history"]
    return {
        "iter": np.array([r["global_iter"] for r in h], dtype=float),
        "res_mhz": np.array([abs(r["residual"]) / TWO_PI * 1e3 for r in h]),
        "cost": np.array([r["cost"] for r in h], dtype=float),
        "V": np.array([r["V"] for r in h], dtype=float),
        "phase": np.array([r["phase"] for r in h]),
    }


def compute(recompute=False):
    """Run hybrid + Ramsey-only baseline, flatten histories. Cached to npz."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    def timed(fn, label):
        """Run one strategy under the true-solver counter; report + verify."""
        import time
        with CountSolves() as cnt:
            t0 = time.perf_counter()
            r = fn()
            secs = time.perf_counter() - t0
        n = cnt.n
        # TRUE cumulative cost aligned to history rows. Bracketing steppers make
        # extra unlogged measurements, so keep only the LAST len(history) marks
        # (all three strategies here use gradient, which logs every measurement).
        hist_n = len(r["history"])
        r["_true_cost"] = np.array(cnt.marks[-hist_n:], dtype=float) \
            if len(cnt.marks) >= hist_n else np.array(cnt.marks, dtype=float)
        r["_time_marks"] = np.array(cnt.time_marks[-hist_n:], dtype=float) \
            if len(cnt.time_marks) >= hist_n else np.array(cnt.time_marks, dtype=float)
        tr = true_residual_mhz(r["V_final"])
        print(f"{label:22s} iters={r['metrics']['stage_iters']} "
              f"conv={r['converged']}  reported|res|="
              f"{abs(r['residual'])/TWO_PI*1e3:.5f} MHz  TRUE|res|={tr:.5f} MHz")
        print(f"{'':22s} TRUE solver calls={n}  "
              f"(workflow-billed {r['metrics']['total_cost']})  {secs:.0f}s")
        return r, n, tr, secs

    hyb, n_hyb, tr_hyb, s_hyb = timed(run_hybrid, "hybrid sweet+fit")
    f = _flatten(hyb)
    m = hyb["metrics"]

    trk, n_trk, tr_trk, s_trk = timed(run_transient_track, "transient track+kf")
    t = _flatten(trk)

    base, n_base, tr_base, s_base = timed(run_ramsey_only, "ramsey-only")
    b = _flatten(base)

    data = {
        "iter": f["iter"], "res_mhz": f["res_mhz"], "cost": f["cost"],
        "V": f["V"], "phase": f["phase"],
        # Derived from stage_boundaries/metrics: the legacy switch_iter /
        # coarse_iter aliases are only emitted for the flat preset (stages=None),
        # and run_hybrid() passes explicit stages.
        "switch_iter": np.array([hyb["stage_boundaries"][0]]),
        "total_cost": np.array([m["total_cost"]]),
        "coarse_iter": np.array([m["stage_iters"][0]]),
        "fine_iter": np.array([m["stage_iters"][1]]),
        "cost_transient": np.array([m["stage_costs"][0]]),
        "cost_ramsey": np.array([m["stage_costs"][1]]),
        "base_iter": b["iter"], "base_res_mhz": b["res_mhz"],
        "base_cost": b["cost"], "base_V": b["V"],
        "base_total_cost": np.array([base["metrics"]["total_cost"]]),
        "f_target_offset_mhz": np.array([F_TARGET_OFFSET_MHZ]),
        # transient track+kernel_full (the winning strategy)
        "trk_iter": t["iter"], "trk_res_mhz": t["res_mhz"],
        "trk_cost": t["cost"], "trk_V": t["V"],
        # TRUE per-iteration cumulative cost (the honest x-axis for panels a/b)
        "cost_true": hyb["_true_cost"],
        "trk_cost_true": trk["_true_cost"],
        "base_cost_true": base["_true_cost"],
        # Measured cumulative wall-clock time after every logged iteration.
        "time_true": hyb["_time_marks"],
        "trk_time_true": trk["_time_marks"],
        "base_time_true": base["_time_marks"],
        # Analytic verification at every visited flux point. These arrays are
        # independent of each strategy's own frequency estimator.
        "verified_res_mhz": np.array([true_residual_mhz(v) for v in f["V"]]),
        "trk_verified_res_mhz": np.array([true_residual_mhz(v) for v in t["V"]]),
        "base_verified_res_mhz": np.array([true_residual_mhz(v) for v in b["V"]]),
        # TRUE solver counts (mesolve+sesolve) — the honest cost axis
        "n_true": np.array([n_hyb, n_trk, n_base], dtype=float),
        "true_res": np.array([tr_hyb, tr_trk, tr_base]),
        "secs": np.array([s_hyb, s_trk, s_base]),
        "strategy": np.array(["hybrid sweet+fit", "transient track+kernel_full",
                              "ramsey-only"]),
    }
    C.save_npz(SUBDIR, _CACHE, **data)
    return data


def _eps_mhz():
    return 1e-4 / TWO_PI * 1e3  # epsilon_f=1e-4 rad·GHz → MHz


def plot_convergence(d, out):
    """residual vs iteration, coarse/fine phases shaded, switch marked."""
    import matplotlib.pyplot as plt

    it, res, phase = d["iter"], d["res_mhz"], d["phase"]
    switch = float(d["switch_iter"][0])
    fig, ax = plt.subplots(figsize=(8, 5))

    # phase shading
    ax.axvspan(it.min() - 0.5, switch + 0.5, color=PHASE_COLOR["coarse"],
               alpha=0.12, zorder=0,
               label=f"coarse (transient, {int(d['coarse_iter'][0])} iters)")
    ax.axvspan(switch + 0.5, it.max() + 0.5, color=PHASE_COLOR["fine"],
               alpha=0.12, zorder=0,
               label=f"fine (Ramsey, {int(d['fine_iter'][0])} iters)")
    ax.axvline(switch + 0.5, color="k", ls="--", lw=1.1, zorder=2)
    ax.annotate("switch", xy=(switch + 0.5, res.max()), xytext=(6, -4),
                textcoords="offset points", fontsize=9, color="k")

    for ph in ("coarse", "fine"):
        sel = phase == ph
        if sel.any():
            ax.semilogy(it[sel], res[sel] + 1e-4, "o-", color=PHASE_COLOR[ph],
                        ms=6, lw=1.4, zorder=3)
    # the winning strategy + the baseline, for context against the hybrid
    if "trk_res_mhz" in d:
        ax.semilogy(d["trk_iter"], d["trk_res_mhz"] + 1e-4, "^-",
                    color="#984ea3", ms=7, lw=1.6, zorder=4,
                    label="transient track+kernel_full (single stage)")
    ax.semilogy(d["base_iter"], d["base_res_mhz"] + 1e-4, "s--", color="#e41a1c",
                ms=5, lw=1.2, mfc="none", zorder=3, label="Ramsey-only baseline")
    ax.axhline(_eps_mhz(), color="gray", ls="--", lw=0.9,
               label=f"tol ε={_eps_mhz():.3f} MHz")
    ax.set_xlabel("Global iteration", fontsize=11)
    ax.set_ylabel(r"|residual| $|f_q-f_{target}|$ (MHz)", fontsize=11)
    ax.set_title("F4 — Closed-loop convergence "
                 f"(target {F_TARGET_OFFSET_MHZ:+.0f} MHz)", fontsize=12)
    ax.legend(fontsize=9); ax.grid(True, which="both", ls=":", alpha=0.4)

    fig.tight_layout()
    fig.savefig(os.path.join(out, "freq_calibration_convergence_sqc.png"), dpi=300)
    fig.savefig(os.path.join(out, "freq_calibration_convergence_sqc.pdf"))
    plt.close(fig)


def plot_cost(d, out):
    """residual vs cumulative measurement cost + the speed-up crossover."""
    import matplotlib.pyplot as plt

    phase = d["phase"]
    fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(17.5, 5))

    # TRUE cost axis throughout — the workflow-billed cost undercounts the
    # transient routes by 45-82x, which would flatter them here.
    ck = "cost_true" if "cost_true" in d else "cost"
    bk = "base_cost_true" if "base_cost_true" in d else "base_cost"

    # hybrid trace, split by phase so the cheap coarse leg is visible
    for ph in ("coarse", "fine"):
        sel = phase == ph
        if sel.any():
            ax.loglog(d[ck][sel], d["res_mhz"][sel] + 1e-4, "o-",
                      color=PHASE_COLOR[ph], ms=6, lw=1.5,
                      label=f"hybrid — {ph}", zorder=3)
    # connect the two legs
    ax.loglog(d[ck], d["res_mhz"] + 1e-4, "-", color="gray", lw=0.8,
              alpha=0.6, zorder=2)
    # Ramsey-only baseline
    ax.loglog(d[bk], d["base_res_mhz"] + 1e-4, "s--", color="#e41a1c",
              ms=5, lw=1.3, mfc="none", label="Ramsey-only baseline", zorder=3)
    # the winning single-stage strategy, on the same TRUE cost axis
    tk = "trk_cost_true" if "trk_cost_true" in d else "trk_cost"
    if tk in d:
        ax.loglog(d[tk], d["trk_res_mhz"] + 1e-4, "^-", color="#984ea3",
                  ms=7, lw=1.6, zorder=4, label="transient track+kernel_full")

    ax.axhline(_eps_mhz(), color="gray", ls="--", lw=0.9,
               label=f"tol ε={_eps_mhz():.3f} MHz")
    # TRUE per-measurement cost, measured. The workflow's billed model (transient
    # = 2, Ramsey = 800, "400×") is wrong: kernel_full's sesolve calls dominate a
    # transient measurement, so the real ratio is ~5×, which is exactly why the
    # observed speed-up is ~4× and not ~400×.
    ct = d["n_true"][1] / max(len(d["trk_iter"]), 1) if "n_true" in d else np.nan
    cr = d["n_true"][2] / max(len(d["base_iter"]), 1) if "n_true" in d else np.nan
    ax.set_xlabel("Cumulative TRUE solver calls (mesolve + sesolve)",
                  fontsize=11)
    ax.set_ylabel(r"|residual| $|f_q-f_{target}|$ (MHz)", fontsize=11)
    ax.set_title("(a) Residual vs cumulative cost", fontsize=11)
    ax.annotate(f"TRUE per-measurement cost:\n"
                f"  transient (track+kf) $\\approx$ {ct:.0f}\n"
                f"  Ramsey double-sweep $\\approx$ {cr:.0f}   "
                f"({cr / max(ct, 1):.1f}$\\times$)\n"
                f"workflow bills transient as 2 (45-82$\\times$ undercount)",
                xy=(0.03, 0.05), xycoords="axes fraction", fontsize=8,
                color="dimgray",
                bbox=dict(boxstyle="round", fc="white", ec="lightgray", alpha=0.9))
    ax.legend(fontsize=9); ax.grid(True, which="both", ls=":", alpha=0.4)

    # -- (b) speed-up vs residual threshold: where each strategy actually pays --
    thr = THRESHOLDS_MHZ
    cb = np.array([cost_to_reach(d["base_res_mhz"], d[bk], t) for t in thr])
    ch = np.array([cost_to_reach(d["res_mhz"], d[ck], t) for t in thr])
    speedup = cb / ch
    ax2.loglog(thr, speedup, "o-", color=PHASE_COLOR["fine"], ms=7, lw=1.6,
               label="hybrid sweet+fit")
    if tk in d:
        ct = np.array([cost_to_reach(d["trk_res_mhz"], d[tk], t) for t in thr])
        ax2.loglog(thr, cb / ct, "^-", color="#984ea3", ms=7, lw=1.6,
                   label="transient track+kernel_full")
    ax2.axhline(1.0, color="k", ls="--", lw=1.0, label="break-even (1×)")
    ax2.axvline(SWITCH_MHZ, color=PHASE_COLOR["coarse"], ls=":", lw=1.2,
                label=f"switch threshold ({SWITCH_MHZ:.0f} MHz)")
    ax2.invert_xaxis()   # tightening tolerance to the right
    ax2.set_xlabel("Residual threshold reached (MHz)  →  tighter", fontsize=11)
    ax2.set_ylabel("Cost speed-up of hybrid (×)", fontsize=11)
    ax2.set_title("(b) Where each strategy pays off", fontsize=11)
    for t, s in zip(thr, speedup):
        if np.isfinite(s) and s >= 10:
            ax2.annotate(f"{s:.0f}×", xy=(t, s), xytext=(0, 8),
                         textcoords="offset points", ha="center", fontsize=9,
                         color="#984ea3")
    ax2.legend(fontsize=8.5, loc="upper right")
    ax2.grid(True, which="both", ls=":", alpha=0.4)

    # -- (c) TRUE solver cost to converge — the headline --
    if "n_true" in d:
        names = ["hybrid\nsweet+fit", "transient\ntrack+kernel_full",
                 "Ramsey\nonly"]
        cols = [PHASE_COLOR["fine"], "#984ea3", "#e41a1c"]
        n_true = d["n_true"]
        bars = ax3.bar(names, n_true, color=cols, alpha=0.85)
        ax3.set_yscale("log")
        ax3.set_ylabel("TRUE solver calls to converge\n(mesolve + sesolve)",
                       fontsize=10.5)
        ax3.set_title("(c) Honest cost to reach tolerance", fontsize=11)
        base_n = n_true[-1]
        for i, (bar, n, tr) in enumerate(zip(bars, n_true, d["true_res"])):
            if i == len(n_true) - 1:
                tag = "baseline"
            else:
                r = base_n / n
                tag = f"{r:.1f}× cheaper" if r >= 1 else f"{1/r:.1f}× costlier"
            ax3.annotate(f"{n:.0f}\n({tag})\ntrue res {tr:.4f} MHz",
                         xy=(bar.get_x() + bar.get_width() / 2, n),
                         xytext=(0, 4), textcoords="offset points",
                         ha="center", fontsize=8)
        # bottom well below the smallest bar (984) so it is not clipped by the
        # log axis, top with headroom for the annotations
        ax3.set_ylim(bottom=100, top=n_true.max() * 40)
        ax3.grid(True, axis="y", which="both", ls=":", alpha=0.4)
        ax3.annotate("sesolve counted 1:1 with mesolve (conservative — sesolve\n"
                     "is the cheaper call), so speed-ups are LOWER bounds",
                     xy=(0.5, 0.95), xycoords="axes fraction", ha="center",
                     va="top", fontsize=7.5, color="dimgray")

    fig.suptitle("F4 — Closed-loop frequency calibration: cost anatomy "
                 r"(drive tracking + kernel_full $G_3$ wins)", fontsize=12.5)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "freq_calibration_cost_sqc.png"), dpi=300)
    fig.savefig(os.path.join(out, "freq_calibration_cost_sqc.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    recompute = "--recompute" in sys.argv
    d = compute(recompute=recompute)

    tot, base_tot = float(d["total_cost"][0]), float(d["base_total_cost"][0])
    print(f"iters: coarse={int(d['coarse_iter'][0])} fine={int(d['fine_iter'][0])}"
          f"  switch@{int(d['switch_iter'][0])}")
    print(f"final |res|: hybrid={d['res_mhz'][-1]:.4f} MHz  "
          f"ramsey-only={d['base_res_mhz'][-1]:.4f} MHz")
    print(f"total cost:  hybrid={tot:.0f}  ramsey-only={base_tot:.0f}  "
          f"({base_tot / tot:.2f}× ) — total is NOT the fair metric "
          f"(different final accuracy)")
    if "n_true" in d:
        print("TRUE solver cost to converge (mesolve+sesolve), "
              "with independently verified residual:")
        base_n = d["n_true"][-1]
        for name, n, tr in zip(d["strategy"], d["n_true"], d["true_res"]):
            print(f"  {str(name):30s} {n:7.0f} calls  "
                  f"({base_n / n:6.1f}× vs Ramsey-only)  "
                  f"TRUE|res|={tr:.5f} MHz")
    # TRUE-cost keys where available (billed cost undercounts transient 45-82×)
    ck = "cost_true" if "cost_true" in d else "cost"
    bk = "base_cost_true" if "base_cost_true" in d else "base_cost"
    tk = "trk_cost_true" if "trk_cost_true" in d else None
    print("TRUE cost to first reach a residual threshold (speed-up vs Ramsey):")
    print(f"  {'thr (MHz)':>10} {'hybrid':>8} {'track+kf':>9} {'ramsey':>8}"
          f" {'hyb':>7} {'trk':>7}")
    for t in THRESHOLDS_MHZ:
        h = cost_to_reach(d["res_mhz"], d[ck], t)
        b = cost_to_reach(d["base_res_mhz"], d[bk], t)
        tr = cost_to_reach(d["trk_res_mhz"], d[tk], t) if tk else np.nan
        print(f"  {t:>10.3f} {h:>8.0f} {tr:>9.0f} {b:>8.0f}"
              f" {b / h:>6.1f}× {b / tr:>6.1f}×")

    out = C.out_dir(SUBDIR)
    plot_convergence(d, out)
    plot_cost(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
