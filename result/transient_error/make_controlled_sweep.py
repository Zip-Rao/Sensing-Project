"""Controlled detuning sweep: order-1 vs order-3 vs cubic law, with the
order-3 breakdown (cubic turning point) marked explicitly."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from transient_error_analysis import p_diff, transient_measure, TWO_PI
from make_figures import robust_G1_G3

G1, G3 = robust_G1_G3()
dstar = np.sqrt(max(-2 * G1 / G3, 0))     # rad/ns

# sweep detuning directly (flux-independent view of the pure cubic effect)
delta = np.linspace(-0.20, 0.20, 81)      # rad/ns  (~ +-32 MHz)
err1, err3, pred = [], [], []
for d in delta:
    err1.append(transient_measure(d, G1, G3, 1) - d)
    err3.append(transient_measure(d, G1, G3, 3) - d)
    pred.append((G3 / (6.0 * G1)) * d**3)
err1 = np.array(err1) / TWO_PI * 1e3      # MHz
err3 = np.array(err3) / TWO_PI * 1e3
pred = np.array(pred) / TWO_PI * 1e3
dmhz = delta / TWO_PI * 1e3
dstar_mhz = dstar / TWO_PI * 1e3

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5))

axA.axhline(0, color="gray", lw=0.5)
axA.plot(dmhz, err1, "r-", lw=2, label="order-1 error (actual)")
axA.plot(dmhz, pred, "b--", lw=2, label=r"$\frac{G_3}{6G_1}\Delta^3$ (linear-trunc.)")
axA.set_xlabel(r"True detuning $\Delta/2\pi$ (MHz)")
axA.set_ylabel("Error (MHz)")
axA.set_title("Order-1 error IS the cubic linear-truncation term")
axA.grid(True, alpha=0.3)
axA.legend(fontsize=9)

axB.axhline(0, color="gray", lw=0.5)
axB.plot(dmhz, err1, "r-", lw=1.5, alpha=0.6, label="order-1 error")
axB.plot(dmhz, err3, "g-", lw=2, label="order-3 (cubic Newton) error")
for s in (+1, -1):
    axB.axvline(s * dstar_mhz, color="purple", ls=":", alpha=0.7)
axB.text(dstar_mhz, axB.get_ylim()[1] * 0.5,
         rf"  $\Delta^*\approx{dstar_mhz:.1f}$ MHz" + "\n  (Newton breakdown)",
         color="purple", fontsize=8, va="top")
axB.set_ylim(-8, 8)
axB.set_xlabel(r"True detuning $\Delta/2\pi$ (MHz)")
axB.set_ylabel("Error (MHz)")
axB.set_title("Order-3 suppresses the cubic error (within safe zone)")
axB.grid(True, alpha=0.3)
axB.legend(fontsize=9)

plt.tight_layout()
fig.savefig("transient_err_controlled_sweep.png", dpi=130)
plt.close(fig)

# quantify reduction strictly inside the safe zone |Delta| < 0.8 Delta*
safe = np.abs(delta) < 0.8 * dstar
r1 = np.sqrt(np.mean(err1[safe]**2))
r3 = np.sqrt(np.mean(err3[safe]**2))
print(f"G1={G1:.4f} G3={G3:.2f} G3/G1={G3/G1:.2f} ns^2  Delta*={dstar_mhz:.2f} MHz")
print(f"Inside safe zone |Delta|<{0.8*dstar_mhz:.1f} MHz: "
      f"RMS order1={r1:.4f} MHz, order3={r3:.5f} MHz, reduction x{r1/max(r3,1e-9):.0f}")
print("saved transient_err_controlled_sweep.png")
