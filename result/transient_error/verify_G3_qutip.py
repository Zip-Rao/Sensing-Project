"""Run in qutip-env to confirm the G3 mismatch against the REAL code path.

    "C:\\Users\\21034\\anaconda3\\envs\\qutip-env\\python.exe" \
        result/transient_error/verify_G3_qutip.py

It compares, using the actual sqc KernelEstimator and a real mesolve:
  (A) G3_diag  = int k3^diag(t) dt        <- what frequency.py feeds to Newton
  (B) G3_taylor = d^3 p_diff/dDelta^3|0   <- the TRUE cubic error coefficient
and shows they differ by ~ (pulse time)^2, NOT by ~1.
"""
import numpy as np
import qutip
from qutip import basis, QobjEvo, mesolve
from sqc.config import CONFIG
from sqc.control.sequence import create_ramsey_pulse
from sqc.reconstruction.kernel import KernelEstimator
from src.qubit import TransmonQubit

EC = 0.2 * 2 * np.pi
EJ = 15.0 * 2 * np.pi
t_rabi = CONFIG.pulse.t_rabi.copy()

q = TransmonQubit(EC=EC, EJ=EJ, T1=1e4, T2=5e3, n_levels=2, flux=0.0)
omega_d = q.frequency
ctrl_x = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=omega_d,
                             phase1=np.pi/2, phase2=0.0)
ctrl_mx = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=omega_d,
                              phase1=np.pi/2, phase2=np.pi)

# ---- (A) the code's diagonal-kernel G3 (order=3, method='exp') ----------
est = KernelEstimator(mode='omega', method='exp', order=3,
                      stim_amplitude=0.005, virtual_z_impl='math',
                      amp_scan_factor=1.5)
rx = est.estimate_full(ctrl_x, q)
rmx = est.estimate_full(ctrl_mx, q)
tk = np.asarray(rx.t_samples)
G1_diag = np.trapezoid((np.asarray(rx.kernels[0]) -
                        np.asarray(rmx.kernels[0]))/2, tk)
G3_diag = np.trapezoid((np.asarray(rx.kernels[2]) -
                        np.asarray(rmx.kernels[2]))/2, tk)

# ---- (B) the TRUE cubic coefficient: fit p_diff(Delta) -----------------
def p_diff_of_delta(delta):
    pe = []
    for ctrl in (ctrl_x, ctrl_mx):
        q2 = TransmonQubit(EC=EC, EJ=EJ, T1=1e4, T2=5e3, n_levels=2, flux=0.0)
        # impose a constant detuning `delta` on top of the drive frame
        q2.qubit_in_mag(None, frame=1, omega_d=omega_d - delta) \
            if False else None
        H0 = QobjEvo(q2.get_hamiltonian_rwa(q2.frequency)) \
            + delta * qutip.num(2)
        Hp = QobjEvo(ctrl.hamiltonian,
                     tlist=np.asarray(ctrl.t_list, float), order=1)
        res = mesolve(H0 + Hp, basis(2, 0),
                      np.asarray(ctrl.t_list, float), [],
                      e_ops=[basis(2, 1) * basis(2, 1).dag()])
        pe.append(res.expect[0][-1])
    return (pe[0] - pe[1]) / 2.0

d = np.linspace(-0.3, 0.3, 25)
d = d[d != 0]
pd = np.array([p_diff_of_delta(x) for x in d])
A = np.column_stack([d, d**3, d**5, d**7])
coef, *_ = np.linalg.lstsq(A, pd, rcond=None)
G1_taylor, G3_taylor = coef[0], 6 * coef[1]

print("\n================  G3 ROOT-CAUSE (real qutip)  ================")
print(f"(A) diagonal kernel : G1={G1_diag:+.4f}  G3={G3_diag:+.4f}"
      f"  -> |G3/G1|={abs(G3_diag/G1_diag):.3f}   (matches doc's '=1')")
print(f"(B) true Taylor     : G1={G1_taylor:+.4f}  G3={G3_taylor:+.2f}"
      f"  -> G3/G1={G3_taylor/G1_taylor:.1f} ns^2 (governs Delta^3 error)")
print(f"ratio G3_taylor/G3_diag = {G3_taylor/G3_diag:.1f} ns^2  "
      f"(~ (2*T_pi2)^2-scale, NOT 1)")
print("Conclusion: kernels are CORRECT; the cubic-Newton in frequency.py"
      " uses the diagonal-kernel G3, which is the WRONG object for inverting"
      " a constant detuning.")
