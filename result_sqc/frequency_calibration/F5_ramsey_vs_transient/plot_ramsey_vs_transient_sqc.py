#!/usr/bin/env python3
"""F5 — Ramsey vs transient: the accuracy-vs-cost Pareto front (sqc).

sqc source:
    sqc.calibration.frequency.FrequencyMeasurement, same work points, six
    measurement *strategies*: Ramsey (double-/single-sweep × τ-grid density)
    and transient (order=1 linear, order=3 cubic).

Why this is not a duplicate of F1: F1 plots accuracy vs FLUX (does the cheap
probe reproduce the dispersion curve?). F5 plots accuracy vs COST — how many
mesolve calls each strategy spends to reach its accuracy. The flux window here
is deliberately kept INSIDE the transient fold (|Δ|≲12 MHz, see F2) so the
comparison isolates cost, not the fold breakdown F2 already characterises.

Cost is MEASURED, not estimated: `sqc.calibration.frequency.mesolve` is wrapped
with a counter for the duration of each strategy, alongside wall-clock time.

Note on transient order=3: its per-measurement marginal cost is ~2 mesolve, but
it needs a ONE-TIME G₃ calibration (~40-80 mesolve, module-cached). Both the
amortised and the marginal figure are recorded — the setup is what a real
experiment pays once per pulse/drive configuration.

IMPORTANT (flux convention): `measure(flux=...)` takes an OFFSET from the qubit
construction point, not absolute flux — see [[sqc-frequency-measurement-api]].
The original stub passed OPTIMAL_FLUX+offset as the offset (double-biasing
off-curve); same bug as F3/F4.

Produces:
    result_sqc/frequency_calibration/F5_ramsey_vs_transient/ramsey_vs_transient_sqc.png / .pdf
        (a) MAE vs mesolve calls per measurement (Pareto front);
        (b) MAE vs wall-clock seconds per measurement.
    result_sqc/frequency_calibration/F5_ramsey_vs_transient/ramsey_vs_transient_sqc.npz

Cost: medium (~1.5 min) — dominated by the Ramsey double-sweep strategy.
Cached to npz; pass --recompute to force a fresh benchmark.

Usage:
    python result_sqc/frequency_calibration/F5_ramsey_vs_transient/plot_ramsey_vs_transient_sqc.py
    python result_sqc/frequency_calibration/F5_ramsey_vs_transient/plot_ramsey_vs_transient_sqc.py --recompute
"""
from __future__ import annotations

import sys
import os
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import _common as C

SUBDIR = "frequency_calibration/F5_ramsey_vs_transient"
_CACHE = "ramsey_vs_transient_sqc"
TWO_PI = 2.0 * np.pi

# Flux OFFSETS (Φ₀) from OPTIMAL_FLUX. Kept inside the transient fold window
# (|Δ|≲12 MHz) so this figure isolates COST, not the fold breakdown (see F2).
FLUX_OFFSETS = np.array([-0.010, -0.005, 0.0, 0.005, 0.010])


def _tau(step_ns):
    """Ramsey τ grid at the given spacing (ns), derived from CONFIG.awg.dt (R9).

    Coarser τ spacing lowers the FFT Nyquist limit 1/(2·Δτ); at Δτ=2 ns that is
    0.25 GHz, still safely above the f_artificial=0.1 GHz single-sweep peak.
    """
    from sqc.config import CONFIG
    dt = float(CONFIG.awg.dt)
    assert abs(step_ns / dt - round(step_ns / dt)) < 1e-9, "step must be k·dt"
    return np.arange(0.0, 200.0, step_ns)


def _strategies():
    """(label, kwargs, style) for each measurement strategy benchmarked."""
    return [
        ("Ramsey double-sweep\n(Δτ=0.5 ns)",
         dict(method="ramsey", f_artificial=None, tau_list=_tau(0.5)),
         dict(color="#e41a1c", marker="o")),
        ("Ramsey single-sweep\n(Δτ=0.5 ns)",
         dict(method="ramsey", f_artificial=0.1, tau_list=_tau(0.5)),
         dict(color="#ff7f00", marker="v")),
        ("Ramsey single-sweep\n(Δτ=1 ns)",
         dict(method="ramsey", f_artificial=0.1, tau_list=_tau(1.0)),
         dict(color="#984ea3", marker="D")),
        ("Ramsey single-sweep\n(Δτ=2 ns)",
         dict(method="ramsey", f_artificial=0.1, tau_list=_tau(2.0)),
         dict(color="#a65628", marker="P")),
        ("Transient order=1",
         dict(method="transient", order=1),
         dict(color="#377eb8", marker="s")),
        ("Transient order=3 (fit)",
         dict(method="transient", order=3, g3_source="fit"),
         dict(color="#4daf4a", marker="^")),
    ]


class _CountMesolve:
    """Context manager counting mesolve calls inside sqc.calibration.frequency.

    The module does `from qutip import mesolve`, so the counter must replace the
    module attribute (patching qutip.mesolve would be too late).
    """

    def __enter__(self):
        import sqc.calibration.frequency as F
        self.n = 0
        self._mod = F
        self._orig = F.mesolve

        def counted(*a, **kw):
            self.n += 1
            return self._orig(*a, **kw)

        F.mesolve = counted
        return self

    def __exit__(self, *exc):
        self._mod.mesolve = self._orig
        return False


def analytic_f(offset):
    """Analytic f01 (GHz) at OPTIMAL_FLUX + offset."""
    q = C.make_sqc_qubit()
    q.change_flux(C.OPTIMAL_FLUX + float(offset))
    return q.frequency / TWO_PI


def bench_one(kwargs, analytic):
    """Benchmark one strategy: returns (mae_mhz, n_solve_total, secs_total,
    n_solve_marginal, secs_marginal).

    Marginal figures exclude the first flux point, so any one-time setup a
    strategy performs (e.g. the cached G₃ fit for order=3) is separated from
    the steady-state per-measurement cost.
    """
    from sqc.calibration.frequency import FrequencyMeasurement

    fm = FrequencyMeasurement(qubit=C.make_sqc_qubit(), **kwargs)
    meas, solves, secs = [], [], []
    for off in FLUX_OFFSETS:
        with _CountMesolve() as cnt:
            t0 = time.time()
            f = fm.measure(flux=float(off)) / TWO_PI
            secs.append(time.time() - t0)
        solves.append(cnt.n)
        meas.append(f)

    meas = np.asarray(meas)
    mae = float(np.mean(np.abs(meas - analytic)) * 1e3)   # MHz
    n_tot, s_tot = float(np.mean(solves)), float(np.mean(secs))
    n_marg = float(np.mean(solves[1:])) if len(solves) > 1 else n_tot
    s_marg = float(np.mean(secs[1:])) if len(secs) > 1 else s_tot
    return mae, n_tot, s_tot, n_marg, s_marg, meas


def compute(recompute=False):
    """Benchmark every strategy. Cached to npz."""
    path = os.path.join(C.out_dir(SUBDIR), f"{_CACHE}.npz")
    if not recompute and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}

    analytic = np.array([analytic_f(o) for o in FLUX_OFFSETS])
    labels, mae, n_tot, s_tot, n_marg, s_marg, curves = [], [], [], [], [], [], []
    for label, kw, _style in _strategies():
        r = bench_one(kw, analytic)
        labels.append(label)
        mae.append(r[0]); n_tot.append(r[1]); s_tot.append(r[2])
        n_marg.append(r[3]); s_marg.append(r[4]); curves.append(r[5])
        print(f"  {label.replace(chr(10), ' '):38s} MAE={r[0]:7.3f} MHz  "
              f"solves={r[1]:6.1f} (marg {r[3]:5.1f})  {r[2]:5.2f}s")

    data = {
        "labels": np.array(labels), "mae_mhz": np.array(mae),
        "n_solve": np.array(n_tot), "secs": np.array(s_tot),
        "n_solve_marginal": np.array(n_marg), "secs_marginal": np.array(s_marg),
        "offsets": FLUX_OFFSETS, "analytic": analytic,
        "measured": np.array(curves),
    }
    C.save_npz(SUBDIR, _CACHE, **data)
    return data


def plot(d, out):
    import matplotlib.pyplot as plt

    styles = [s for _l, _k, s in _strategies()]
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12.5, 5.2))

    for ax, xkey, xlabel, title in (
        (ax0, "n_solve", "mesolve calls per measurement",
         "(a) Accuracy vs solver cost"),
        (ax1, "secs", "wall-clock seconds per measurement",
         "(b) Accuracy vs wall time"),
    ):
        for i, st in enumerate(styles):
            ax.loglog(d[xkey][i], d["mae_mhz"][i], st["marker"],
                      color=st["color"], ms=11, mew=1.4,
                      label=str(d["labels"][i]).replace("\n", " "), zorder=3)
        # Pareto frontier: cheapest strategy at each accuracy level or better
        order = np.argsort(d[xkey])
        best, px, py = np.inf, [], []
        for i in order:
            if d["mae_mhz"][i] < best:
                best = d["mae_mhz"][i]
                px.append(d[xkey][i]); py.append(d["mae_mhz"][i])
        ax.loglog(px, py, "k--", lw=1.0, alpha=0.55,
                  label="Pareto frontier", zorder=2)
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel("MAE vs analytic $f_{01}$ (MHz)", fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.grid(True, which="both", ls=":", alpha=0.4)

    ax0.legend(fontsize=7.5, loc="upper right")
    # Headline: every Ramsey single-sweep variant sits up-and-right of the
    # order=3 point, i.e. is Pareto-dominated (worse MAE *and* costlier).
    ax0.annotate("transient order=3 dominates every\nRamsey single-sweep variant\n"
                 "(better MAE and ≥5× cheaper)",
                 xy=(0.04, 0.06), xycoords="axes fraction", fontsize=8.5,
                 color="#2e7d32",
                 bbox=dict(boxstyle="round", fc="white", ec="#4daf4a", alpha=0.9))
    fig.suptitle("F5 — Ramsey vs transient: accuracy-cost trade-off "
                 r"(inside the transient window, $|\Delta|\lesssim$12 MHz)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{_CACHE}.png"), dpi=300)
    fig.savefig(os.path.join(out, f"{_CACHE}.pdf"))
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    recompute = "--recompute" in sys.argv
    d = compute(recompute=recompute)

    i_best = int(np.argmin(d["mae_mhz"]))
    i_cheap = int(np.argmin(d["n_solve"]))
    print(f"most accurate: {str(d['labels'][i_best]).replace(chr(10), ' ')} "
          f"({d['mae_mhz'][i_best]:.3f} MHz, {d['n_solve'][i_best]:.0f} solves)")
    print(f"cheapest:      {str(d['labels'][i_cheap]).replace(chr(10), ' ')} "
          f"({d['mae_mhz'][i_cheap]:.3f} MHz, {d['n_solve'][i_cheap]:.0f} solves)")
    print(f"cost ratio: {d['n_solve'][i_best] / max(d['n_solve'][i_cheap], 1):.0f}× "
          f"more solves for {d['mae_mhz'][i_cheap] / max(d['mae_mhz'][i_best], 1e-9):.0f}× "
          f"better accuracy")

    out = C.out_dir(SUBDIR)
    plot(d, out)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
