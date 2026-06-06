"""Generate transient-method error-source figures.

Two figures (flux=0.0 sweet spot, flux=0.9 biased), each with:
  (left)  detuning Delta vs added flux fa, plus error curves
  (right) order-1 transient error vs analytic cubic law (G3/6G1)Delta^3,
          and order-3 corrected error.
Saves PNGs to the workspace.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from transient_error_analysis import (
    omega_q, p_diff, transient_measure, TWO_PI,
)


def robust_G1_G3(drange=0.3, n=41):
    """Wide-range odd-polynomial fit -> stable G1, G3."""
    d = np.linspace(-drange, drange, n)
    d = d[d != 0]
    pd = np.array([p_diff(x) for x in d])
    A = np.column_stack([d, d**3, d**5, d**7])
    coef, *_ = np.linalg.lstsq(A, pd, rcond=None)
    return coef[0], 6.0 * coef[1]


def build(flux_bias, G1, G3, n_scan=15, fa_max=0.01):
    omega_d = omega_q(flux_bias)
    fa = np.linspace(0, fa_max, n_scan)
    delta = np.array([omega_q(flux_bias + x) - omega_d for x in fa])
    f_true = omega_d + delta
    f1 = np.array([omega_d + transient_measure(d, G1, G3, 1) for d in delta])
    f3 = np.array([omega_d + transient_measure(d, G1, G3, 3) for d in delta])
    err1 = (f1 - f_true) / TWO_PI * 1e3          # MHz
    err3 = (f3 - f_true) / TWO_PI * 1e3
    err_pred = (G3 / (6.0 * G1)) * delta**3 / TWO_PI * 1e3
    return dict(fa=fa, total=flux_bias + fa, delta=delta / TWO_PI * 1e3,
                err1=err1, err3=err3, err_pred=err_pred, omega_d=omega_d)


def make_fig(flux_bias, G1, G3, fname):
    D = build(flux_bias, G1, G3)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

    # left: detuning that the qubit actually experiences
    axL.plot(D["total"], D["delta"], "k.-", lw=1.5, ms=7)
    axL.axvline(flux_bias, color="gray", ls="--", alpha=0.5, label="bias point")
    axL.set_xlabel(r"Total flux ($\Phi_0$)")
    axL.set_ylabel(r"True detuning $\Delta/2\pi$ (MHz)")
    axL.set_title(f"Detuning vs flux (bias={flux_bias})\n"
                  rf"$\kappa$ sets how fast $\Delta$ grows")
    axL.grid(True, alpha=0.3)
    axL.legend(fontsize=9)

    # right: error decomposition
    axR.axhline(0, color="gray", lw=0.5)
    axR.plot(D["total"], D["err1"], "r^-", ms=6,
             label="order-1 transient error (actual)")
    axR.plot(D["total"], D["err_pred"], "b--", lw=2,
             label=r"linear-trunc. prediction $\frac{G_3}{6G_1}\Delta^3$")
    axR.plot(D["total"], D["err3"], "gs-", ms=5,
             label="order-3 (cubic Newton) error")
    axR.set_xlabel(r"Total flux ($\Phi_0$)")
    axR.set_ylabel("Measurement error (MHz)")
    axR.set_title(f"Transient error source (bias={flux_bias})")
    axR.grid(True, alpha=0.3)
    axR.legend(fontsize=9)

    plt.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)

    # agreement metric
    m = np.abs(D["err_pred"]) > 1e-6
    ratio = (D["err1"][m] / D["err_pred"][m]) if m.any() else np.array([np.nan])
    rms1 = np.sqrt(np.mean(D["err1"]**2))
    rms3 = np.sqrt(np.mean(D["err3"]**2))
    return dict(ratio_mean=float(np.nanmean(ratio)),
                ratio_std=float(np.nanstd(ratio)),
                rms1=rms1, rms3=rms3,
                delta_max=float(D["delta"][-1]))


if __name__ == "__main__":
    G1, G3 = robust_G1_G3()
    print(f"Robust kernel integrals: G1={G1:.5f}, G3={G3:.4f}, "
          f"G3/G1={G3/G1:.3f} ns^2")
    dstar = np.sqrt(max(-2 * G1 / G3, 0)) / TWO_PI * 1e3
    print(f"Cubic turning point (order-3 safe-zone edge): Delta* = "
          f"{dstar:.2f} MHz\n")
    for bias, fn in ((0.0, "transient_err_source_flux00.png"),
                     (0.9, "transient_err_source_flux09.png")):
        stats = make_fig(bias, G1, G3, fn)
        print(f"flux={bias}: Delta_max={stats['delta_max']:.2f} MHz | "
              f"err_o1/err_pred = {stats['ratio_mean']:.4f}"
              f"+-{stats['ratio_std']:.4f} | "
              f"RMS o1={stats['rms1']:.4f} MHz o3={stats['rms3']:.5f} MHz | "
              f"-> {fn}")
