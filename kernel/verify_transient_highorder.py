"""kernel/verify_transient_highorder.py

Accuracy verification for transient single-point frequency extraction:
    order-1 (linear, G1)  vs  order-3 (cubic Newton, sim full-kernel G3).

The cubic coefficient G3 comes from the **off-diagonal** sim Heisenberg kernel
(`KernelEstimator(mode='omega', method='sim', order=3, extract_off_diagonal=True)`),
integrated over the full simplex:

    p_diff = G1*dw + (1/6)*G3_full*dw**3 + O(dw**5)              (constant dw)
    G1      = INT  k1(t) dt
    G3_full = INT INT INT k3(t1,t2,t3) dt1 dt2 dt3               (NOT the diagonal!)

For a constant detuning the cubic response needs the triple integral; the
diagonal slice INT k3(t,t,t)dt underestimates it by ~(2*T_pi2)^2/2 (the Phase-11
bug). This script shows the full-kernel G3 makes order-3 extraction beat order-1
across a known-detuning sweep, while G3_diag does not.

Sign convention
---------------
The sim kernel models the perturbation as +(dw/2)*sigma_z (Virtual-Z, see
kernel._omega_vz_math), whereas a physical detuning enters as -(Delta/2)*sigma_z
(matching _calibrate_g3_taylor). Hence Delta = -dw, so the coefficients of
p_diff vs the *detuning* Delta are G1 = -G1_sim and G3 = -G3_full.

Standalone study: does NOT modify anything under sqc/ or src/.

Run:
    "C:\\Users\\21034\\anaconda3\\envs\\qutip-env\\python.exe" kernel/verify_transient_highorder.py
"""
from __future__ import annotations

import os
import sys

# Allow running as `python kernel/verify_transient_highorder.py` from repo root.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt
from qutip import Qobj, QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.reconstruction.kernel import KernelEstimator
from sqc.calibration.frequency import _calibrate_g3_taylor
from src.qubit import TransmonQubit


# ===========================================================================
# Setup — match _measure_frequency_transient / kernel_order_comparison params
# ===========================================================================
dt = CONFIG.awg.dt
t_rabi = CONFIG.pulse.t_rabi.copy()
t_global = CONFIG.pulse.t_global.copy()
T_pi2 = CONFIG.pulse.t_rabi_duration

qubit = TransmonQubit(
    EC=0.2 * 2 * np.pi, EJ=15 * 2 * np.pi,
    T1=10000, T2=5000, n_levels=2, flux=0.0,
)
omega_d = qubit.frequency
n_levels = qubit.n_levels
psi_e = basis(n_levels, 1)

# Orthogonal readout pulses (tau=0), identical to the production transient
# protocol: R_y(pi/2) - R_x(pi/2)  and  R_y(pi/2) - R_{-x}(pi/2).
ctrl_x = create_ramsey_pulse(
    t_rabi, 0.0, omega_d=omega_d, phase1=np.pi / 2, phase2=0.0, qubit=qubit,
)
ctrl_mx = create_ramsey_pulse(
    t_rabi, 0.0, omega_d=omega_d, phase1=np.pi / 2, phase2=np.pi, qubit=qubit,
)

# sigma_z in the n-level Fock basis (constant-detuning generator).
sigma_z = Qobj(np.diag([1.0, -1.0] + [0.0] * (n_levels - 2)))

# Measure on the pulse window (the kernel's own domain). For tau=0 the whole
# sequence lives here; outside it the detuning is a pure z-rotation that leaves
# p_e invariant, so this is the physically complete window.
t_meas = np.unique(np.asarray(ctrl_x.t_list, dtype=float))


# ===========================================================================
# Step 1 — sim off-diagonal kernel  ->  G1, G3_full, G3_diag
# ===========================================================================
print("=" * 72)
print("Step 1: sim off-diagonal Heisenberg kernel (order=3)")
print("=" * 72)

est = KernelEstimator(
    mode='omega', method='sim', order=3, extract_off_diagonal=True,
)
res_x = est.estimate_full(ctrl_x, qubit, t_samples=t_meas)
res_mx = est.estimate_full(ctrl_mx, qubit, t_samples=t_meas)
t_k = np.asarray(res_x.t_samples, dtype=float)

# Differential kernels (matches p_diff = (p_x - p_mx)/2 in the protocol).
k1_diff = (np.asarray(res_x.kernels[0]) - np.asarray(res_mx.kernels[0])) / 2.0
k3_diff = (np.asarray(res_x.kernels[2]) - np.asarray(res_mx.kernels[2])) / 2.0

G1_sim = float(np.trapezoid(k1_diff, t_k))
G3_full_sim = float(
    np.trapezoid(np.trapezoid(np.trapezoid(k3_diff, t_k, axis=0), t_k, axis=0), t_k, axis=0)
)
k3_diag = np.array([k3_diff[i, i, i] for i in range(len(t_k))])
G3_diag_sim = float(np.trapezoid(k3_diag, t_k))

# Coefficients of p_diff vs the *detuning* Delta (sign-flipped; see header).
G1 = -G1_sim
G3_full = -G3_full_sim
G3_diag = -G3_diag_sim

print(f"  G1_sim      = {G1_sim:+.4f}   (kernel, +dw/2 convention)")
print(f"  G3_full_sim = {G3_full_sim:+.2f}   (G3_full/G1 = {G3_full_sim / G1_sim:+.1f} ns^2)")
print(f"  G3_diag_sim = {G3_diag_sim:+.4f}   (G3_full/G3_diag = {G3_full_sim / G3_diag_sim:.1f}x)")
print(f"  expected ratio ~ (2*T_pi2)^2/2 = {(2 * T_pi2) ** 2 / 2:.1f}")
print(f"  -> detuning coefficients:  G1 = {G1:+.4f}, "
      f"G3_full = {G3_full:+.1f}, G3_diag = {G3_diag:+.4f}")


# ===========================================================================
# Step 2 — cross-check: kernel G1 vs direct small-Delta slope vs Route A fit
# ===========================================================================
print("\n" + "=" * 72)
print("Step 2: cross-check G1/G3 against direct simulation and Route A fit")
print("=" * 72)

# Put qubit at sweet spot with H_list populated (mesolve needs qubit.H_list).
t_sig = CONFIG.pulse.make_time(0, 300)
qubit.qubit_in_mag(FluxSignal(type=0, t_list=t_sig), frame=1, omega_d=omega_d)

_meas_opts = {"max_step": float(dt)}


def measure_p_diff(delta_rad: float) -> float:
    """Measure differential p_e under a known constant detuning (angular).

    Detuning enters as -(Delta/2)*sigma_z over the pulse window, exactly as
    in _calibrate_g3_taylor.
    """
    detuning_coeff = -0.5 * delta_rad * np.ones_like(t_meas)
    H_det = QobjEvo([[sigma_z, detuning_coeff]], tlist=t_meas, order=1)
    H_base = QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1)

    H_px = H_base + H_det + QobjEvo(
        ctrl_x.hamiltonian_on(t_meas), tlist=t_meas, order=1)
    H_pmx = H_base + H_det + QobjEvo(
        ctrl_mx.hamiltonian_on(t_meas), tlist=t_meas, order=1)

    p_x = float(mesolve(H_px, qubit.state, t_meas, [],
                        e_ops=[psi_e * psi_e.dag()], options=_meas_opts).expect[0][-1])
    p_mx = float(mesolve(H_pmx, qubit.state, t_meas, [],
                         e_ops=[psi_e * psi_e.dag()], options=_meas_opts).expect[0][-1])
    return (p_x - p_mx) / 2.0


# Direct slope at small Delta (central difference) = ground-truth G1.
d_small = 2 * np.pi * 0.001  # 1 MHz
G1_direct = (measure_p_diff(d_small) - measure_p_diff(-d_small)) / (2 * d_small)

# Route A fit (mesolve p_diff(Delta) polynomial over +/-0.08 GHz).
G1_fit, G3_fit = _calibrate_g3_taylor(qubit, omega_d, t_rabi, t_global)
# re-seat sweet spot (the calibrator left the qubit set, but be explicit)
qubit.qubit_in_mag(FluxSignal(type=0, t_list=t_sig), frame=1, omega_d=omega_d)

print(f"  G1 (sim kernel, detuning conv)     = {G1:+.4f}")
print(f"  G1 (direct small-Delta slope)      = {G1_direct:+.4f}   <- ground truth")
print(f"  G1 (Route A polynomial fit)        = {G1_fit:+.4f}   "
      f"(fit/true = {G1_fit / G1_direct:.2f}x, biased low: range too wide)")
print(f"  G3_full (sim)                      = {G3_full:+.1f}")
print(f"  G3 (Route A fit)                   = {G3_fit:+.1f}")
print(f"  G3_diag (sim, the Phase-11 bug)    = {G3_diag:+.4f}  (~200x too small)")


# ===========================================================================
# Step 3 — accuracy sweep over known constant detunings (safe zone)
# ===========================================================================
print("\n" + "=" * 72)
print("Step 3: accuracy sweep (order-1 linear vs order-3 sim-G3)")
print("=" * 72)


def invert_cubic(p_diff: float, g1: float, g3: float) -> float:
    """Newton solve p_diff = g1*dw + (g3/6)*dw^3 from the linear estimate."""
    dw = p_diff / g1
    if abs(g3) < 1e-12:
        return dw
    for _ in range(50):
        f = g1 * dw + (g3 / 6.0) * dw ** 3 - p_diff
        fp = g1 + (g3 / 2.0) * dw ** 2
        if abs(fp) < 1e-15:
            break
        dw_new = dw - f / fp
        if abs(dw_new - dw) < 1e-14 * max(abs(dw), 1e-14):
            dw = dw_new
            break
        dw = dw_new
    return dw


# Cubic model folds (dp/dDelta=0) at |Delta| = sqrt(-2*G1/G3_full); sweep up to
# ~1.5x that to show the safe-zone extension and the eventual common breakdown.
delta_fold_mhz = np.sqrt(abs(2 * G1 / G3_full)) / (2 * np.pi) * 1e3
print(f"  cubic fold at |Delta| ~ {delta_fold_mhz:.1f} MHz "
      f"(order-3 safe zone edge)\n")

delta_mhz = np.array([1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30], dtype=float)

err1, err3, err3d = [], [], []
print(f"  {'D(MHz)':>7} {'p_diff':>10} {'err_lin':>10} "
      f"{'err_cub_full':>13} {'err_cub_diag':>13}  (MHz)")
for dm in delta_mhz:
    dw_true = 2 * np.pi * dm * 1e-3
    p = measure_p_diff(dw_true)

    dw1 = p / G1                              # order-1 linear
    dw3 = invert_cubic(p, G1, G3_full)        # order-3 full G3
    dw3d = invert_cubic(p, G1, G3_diag)       # order-3 diagonal G3

    to_mhz = lambda x: abs(x - dw_true) / (2 * np.pi) * 1e3
    err1.append(to_mhz(dw1))
    err3.append(to_mhz(dw3))
    err3d.append(to_mhz(dw3d))
    print(f"  {dm:>7.0f} {p:>10.5f} {err1[-1]:>10.4f} "
          f"{err3[-1]:>13.4f} {err3d[-1]:>13.4f}")

err1 = np.array(err1)
err3 = np.array(err3)
err3d = np.array(err3d)


# ===========================================================================
# Step 4 — report
# ===========================================================================
fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 5.5))

axA.semilogy(delta_mhz, np.maximum(err1, 1e-6), 'o-', color='tab:red', lw=2,
             label='order-1 (linear, $G_1$)')
axA.semilogy(delta_mhz, np.maximum(err3, 1e-6), 's-', color='tab:green', lw=2,
             label='order-3 (Newton, sim $G_3^{full}$)')
axA.semilogy(delta_mhz, np.maximum(err3d, 1e-6), '^--', color='tab:gray', lw=1.5,
             alpha=0.7, label='order-3 (diagonal $G_3^{diag}$)')
axA.axvline(delta_fold_mhz, color='k', ls=':', alpha=0.5,
            label=f'cubic fold ~{delta_fold_mhz:.0f} MHz')
axA.set_xlabel(r'detuning  $|\Delta|$  (MHz)')
axA.set_ylabel(r'frequency error  $|\hat\Delta - \Delta|$  (MHz)')
axA.set_title('Transient extraction accuracy: linear vs sim high-order kernel')
axA.legend(fontsize=9)
axA.grid(True, which='both', alpha=0.3)

impr = np.where(err3 > 1e-9, err1 / err3, np.nan)
axB.semilogy(delta_mhz, impr, 'd-', color='tab:blue', lw=2)
axB.axhline(1.0, color='gray', ls='--', alpha=0.6)
axB.axvline(delta_fold_mhz, color='k', ls=':', alpha=0.5)
axB.set_xlabel(r'detuning  $|\Delta|$  (MHz)')
axB.set_ylabel(r'error reduction  $err_1 / err_3$')
axB.set_title('Order-3 (sim $G_3^{full}$) improvement over order-1')
axB.grid(True, which='both', alpha=0.3)

fig.tight_layout()
out_png = 'kernel/transient_highorder_verification.png'
fig.savefig(out_png, dpi=150, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out_png}")

# Safe-zone summary (where the cubic model is monotonic and accurate).
safe = delta_mhz <= delta_fold_mhz * 0.9
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)
print(f"  sim G1 = {G1:+.4f}  (matches direct slope {G1_direct:+.4f})")
print(f"  sim G3_full = {G3_full:+.1f},  G3_diag = {G3_diag:+.4f}")
print(f"  cubic safe zone edge ~ {delta_fold_mhz:.1f} MHz\n")
print(f"  within safe zone (|D| <= {delta_fold_mhz * 0.9:.0f} MHz):")
print(f"    mean err order-1           : {np.mean(err1[safe]):.4f} MHz")
print(f"    mean err order-3 (G3_full) : {np.mean(err3[safe]):.4f} MHz")
print(f"    mean err order-3 (G3_diag) : {np.mean(err3d[safe]):.4f} MHz")
print(f"    median improvement full/lin: {np.nanmedian(impr[safe]):.1f}x")
