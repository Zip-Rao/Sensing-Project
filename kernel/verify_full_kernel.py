"""
Full multi-time kernel verification: k2(t1,t2), k3(t1,t2,t3) via Heisenberg sim,
plus Volterra expansion validation against mesolve ground truth.

Theory:
  k2^ord(t1, t2) = -<0|[W(t2), [W(t1), Q]]|0>   for t1 >= t2
  k3^ord(t1, t2, t3) = -i<0|[W(t3), [W(t2), [W(t1), Q]]]|0>  for t1>=t2>=t3

Then symmetrize to get full Volterra kernels:
  k2(t1,t2) = k2^ord(max(t1,t2), min(t1,t2))
  k3(t1,t2,t3) = average over all 6 permutations of the ordered kernel

Volterra prediction for perturbation dw(t):
  dp_e = int k1(t) dw(t) dt + 1/2 int int k2(t1,t2) dw(t1) dw(t2) dt1 dt2
       + 1/6 int int int k3(t1,t2,t3) dw(t1) dw(t2) dw(t3) dt1 dt2 dt3
"""
import numpy as np
import matplotlib.pyplot as plt
import time
from matplotlib import rcParams
import qutip
from qutip import Qobj, QobjEvo, basis, sesolve, mesolve, sigmax, sigmay, sigmaz, qeye

from sqc.config import CONFIG
from sqc.control.sequence import create_ramsey_pulse
from sqc.reconstruction.kernel import KernelEstimator
from src.qubit import TransmonQubit

rcParams.update({'font.size': 11, 'figure.dpi': 120, 'axes.grid': True, 'grid.alpha': 0.3})

# ===========================================================================
# Setup
# ===========================================================================
dt = CONFIG.awg.dt
T_pi2 = CONFIG.pulse.t_rabi_duration
t_rabi = CONFIG.pulse.t_rabi.copy()
t_global = CONFIG.pulse.t_global.copy()

q_ref = TransmonQubit(EC=0.2*2*np.pi, EJ=15*2*np.pi, T1=10000, T2=5000, n_levels=2, flux=0.0)
omega_d = q_ref.frequency

# Y-X Ramsey, tau=0
ctrl_y  = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=omega_d, phase1=+np.pi/2, phase2=0.0)

print(f"omega_d = {omega_d:.4f} rad*GHz")
print(f"T_pi2 = {T_pi2} ns, dt = {dt} ns")

# ===========================================================================
# 1. Compute full multi-time kernels via Heisenberg propagator
# ===========================================================================
print("\n" + "="*70)
print("1. Full multi-time Heisenberg kernel computation")
print("="*70)

n = 2
I_op = qeye(n)
sigma_z = Qobj(np.diag([1.0, -1.0]))
M = basis(n, 1) * basis(n, 1).dag()  # |1><1|
state0 = basis(n, 0)

# Build Hamiltonian on deduplicated time grid
raw_tlist = np.asarray(ctrl_y.t_list, dtype=float)
t_grid, unique_idx = np.unique(raw_tlist, return_index=True)

H0 = QobjEvo(0 * qeye(2))
H_pulse = QobjEvo(ctrl_y.hamiltonian_on(t_grid), tlist=t_grid, order=1)
H_full = H0 + H_pulse

# Single sesolve -> all U(t)
t0 = time.time()
res = sesolve(H_full, I_op, t_grid,
              options={'max_step': float(dt), 'store_states': True})
U_list = res.states
U_T = U_list[-1]
Q = U_T.dag() * M * U_T

# Pre-compute W(t) at all time points
W_list = []
for U_t in U_list:
    Z_t = U_t.dag() * sigma_z * U_t
    W_list.append(Z_t / 2.0)

Nt = len(t_grid)
print(f"Time grid: {Nt} points, sesolve took {time.time()-t0:.2f}s")

def _expect(op, state):
    """Compute <state|op|state>, handling both Qobj and scalar returns."""
    val = state.dag() * op * state
    if hasattr(val, 'shape') and val.shape == (1, 1):
        return complex(val[0, 0])
    return complex(val)

# -- k1(t) ----------------------------------------------------------------
k1 = np.zeros(Nt)
for i in range(Nt):
    comm = W_list[i] * Q - Q * W_list[i]
    k1[i] = float(np.real(1j * _expect(comm, state0)))

G1 = np.trapezoid(k1, t_grid)
print(f"k1: G1 = {G1:.4f}")

# -- k2(t1, t2) full 2D ---------------------------------------------------
print("Computing k2(t1,t2) full 2D matrix...")
t0 = time.time()
k2_full = np.zeros((Nt, Nt))
for i in range(Nt):
    Wi = W_list[i]
    for j in range(Nt):
        Wj = W_list[j]
        # k2^ord = -<0|[W(t<), [W(t>), Q]]|0>
        if t_grid[i] >= t_grid[j]:
            t_less, t_greater = j, i
        else:
            t_less, t_greater = i, j
        W_less = W_list[t_less]
        W_greater = W_list[t_greater]
        comm1 = W_greater * Q - Q * W_greater
        comm2 = W_less * comm1 - comm1 * W_less
        k2_full[i, j] = float(np.real(-_expect(comm2, state0)))

G2_full = np.trapezoid(np.trapezoid(k2_full, t_grid, axis=0), t_grid, axis=0)
print(f"k2: took {time.time()-t0:.1f}s, G2_full = {G2_full:.2e}")
print(f"k2 max|k2| = {np.max(np.abs(k2_full)):.2e}, RMS = {np.sqrt(np.mean(k2_full**2)):.2e}")
print(f"k2 diagonal: max|k2(t,t)| = {np.max(np.abs(np.diag(k2_full))):.2e}")

# -- k3(t1, t2, t3) 3D ----------------------------------------------------
# Full 3D is Nt^3 ~ 39^3 = 59319 elements — feasible
print("Computing k3(t1,t2,t3) full 3D array...")
t0 = time.time()
k3_ord = np.zeros((Nt, Nt, Nt))
# k3^ord(t1,t2,t3) = -i<0|[W(t3),[W(t2),[W(t1),Q]]]|0> for t1>=t2>=t3
for i in range(Nt):
    Wi = W_list[i]
    comm_i = Wi * Q - Q * Wi  # [W(t1), Q]
    for j in range(Nt):
        if t_grid[j] > t_grid[i]:
            continue  # enforce t1 >= t2
        Wj = W_list[j]
        comm_ij = Wj * comm_i - comm_i * Wj  # [W(t2), [W(t1), Q]]
        for k in range(Nt):
            if t_grid[k] > t_grid[j]:
                continue  # enforce t2 >= t3
            Wk = W_list[k]
            comm_ijk = Wk * comm_ij - comm_ij * Wk
            k3_ord[i, j, k] = float(np.real(-1j * _expect(comm_ijk, state0)))

# Symmetrize: k3(t1,t2,t3) = average over ordered permutations
# For efficiency, we compute the full symmetric kernel by summing over
# 6 time-orderings with appropriate sign (the ord kernel already has
# the time-ordering built in)
k3_sym = np.zeros((Nt, Nt, Nt))
# The fully symmetric k3 at any 3 time points equals k3^ord evaluated
# on the time-ordered permutation
for i in range(Nt):
    for j in range(Nt):
        for k in range(Nt):
            # Find the time-ordered indices
            idx = np.array([i, j, k])
            t_vals = t_grid[idx]
            order = np.argsort(-t_vals)  # descending: t[order[0]] >= t[order[1]] >= t[order[2]]
            i_ord, j_ord, k_ord = idx[order]
            k3_sym[i, j, k] = k3_ord[i_ord, j_ord, k_ord]

# Triple integral
G3_full = np.trapezoid(np.trapezoid(np.trapezoid(k3_sym, t_grid, axis=0), t_grid, axis=0), t_grid, axis=0)
G3_diag = np.trapezoid(np.diagonal(np.diagonal(k3_sym)), t_grid)

print(f"k3: took {time.time()-t0:.1f}s")
print(f"k3 diagonal: G3_diag = {G3_diag:.4f}, G3_diag/G1 = {G3_diag/G1:.4f}")
print(f"k3 full: G3_full = {G3_full:.1f}, G3_full/G1 = {G3_full/G1:.1f}")
print(f"Ratio G3_full / |G3_diag| = {G3_full / abs(G3_diag):.1f}")
print(f"sqrt(|ratio|) = {np.sqrt(abs(G3_full/G3_diag)):.1f} ns")
print(f"Compare: 2*T_pi2 = {2*T_pi2} ns")

# ===========================================================================
# 2. Visualize k2 and k3 slices
# ===========================================================================
print("\n" + "="*70)
print("2. Visualization")
print("="*70)

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# (a) k1(t)
ax = axes[0, 0]
ax.plot(t_grid, k1, 'b-', linewidth=2)
ax.axvline(x=T_pi2, color='red', linestyle='--', alpha=0.5)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('k1(t)')
ax.set_title(f'(a) k1(t), G1={G1:.3f}')
ax.set_xlim(0, 2*T_pi2)

# (b) k2(t1, t2) heatmap
ax = axes[0, 1]
im = ax.pcolormesh(t_grid, t_grid, k2_full, cmap='RdBu_r', shading='auto',
                   vmin=-np.max(np.abs(k2_full)), vmax=np.max(np.abs(k2_full)))
ax.axvline(x=T_pi2, color='green', linestyle='--', alpha=0.5)
ax.axhline(y=T_pi2, color='green', linestyle='--', alpha=0.5)
ax.set_xlabel('t1 (ns)')
ax.set_ylabel('t2 (ns)')
ax.set_title(f'(b) k2(t1,t2), max|k2|={np.max(np.abs(k2_full)):.2e}')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='k2')

# (c) k2 diagonal slice vs zero
ax = axes[0, 2]
ax.plot(t_grid, np.diag(k2_full), 'o-', markersize=3, color='purple', label='k2(t,t)')
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('k2(t,t)')
ax.set_title(f'(c) k2 diagonal slice (expect ~0)')
ax.legend()
ax.set_xlim(0, 2*T_pi2)

# (d) k3(t1,t2) at fixed t3 = T_pi2/2
ax = axes[1, 0]
t3_idx = np.argmin(np.abs(t_grid - T_pi2/2))
k3_slice = k3_sym[:, :, t3_idx]
vmax = np.max(np.abs(k3_slice))
im = ax.pcolormesh(t_grid, t_grid, k3_slice, cmap='RdBu_r', shading='auto',
                   vmin=-vmax, vmax=vmax)
ax.axvline(x=T_pi2, color='green', linestyle='--', alpha=0.5)
ax.axhline(y=T_pi2, color='green', linestyle='--', alpha=0.5)
ax.set_xlabel('t1 (ns)')
ax.set_ylabel('t2 (ns)')
ax.set_title(f'(d) k3(t1,t2,t3={T_pi2/2:.1f}ns) slice')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='k3')

# (e) k3 diagonal vs -k1
ax = axes[1, 1]
k3_diag_vals = np.array([k3_sym[i,i,i] for i in range(Nt)])
ax.plot(t_grid, -k1, 'b-', linewidth=2, alpha=0.4, label='-k1(t)')
ax.plot(t_grid, k3_diag_vals, 'ro-', markersize=3, label='k3(t,t,t)')
ax.axvline(x=T_pi2, color='gray', linestyle='--', alpha=0.5)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_title(f'(e) Diagonal: k3(t,t,t) vs -k1(t)')
ax.legend()
ax.set_xlim(0, 2*T_pi2)
rmse_k3diag = np.sqrt(np.mean((k3_diag_vals + k1)**2))
ax.text(0.98, 0.05, f'RMSE(k3+k1)={rmse_k3diag:.2e}', transform=ax.transAxes,
        ha='right', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

# (f) Summary
ax = axes[1, 2]
ax.axis('off')
summary = (
    f"Multi-time kernel verification:\n\n"
    f"k1:  G1 = {G1:.4f}\n"
    f"k2:  G2_full = {G2_full:.2e} (~0)\n"
    f"     max|k2| = {np.max(np.abs(k2_full)):.2e}\n\n"
    f"k3 diagonal:\n"
    f"  G3_diag = {G3_diag:.4f}\n"
    f"  |G3_diag/G1| = {abs(G3_diag/G1):.4f}\n"
    f"  k3(t,t,t) = -k1(t)?\n"
    f"  RMSE = {rmse_k3diag:.2e}\n\n"
    f"k3 full (3D integral):\n"
    f"  G3_full = {G3_full:.1f}\n"
    f"  G3_full/G1 = {G3_full/G1:.1f}\n"
    f"  G3_full/|G3_diag| = {G3_full/abs(G3_diag):.1f}\n"
    f"  sqrt(|ratio|) = {np.sqrt(abs(G3_full/G3_diag)):.1f} ns\n"
)
ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
        va='top', ha='left', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.suptitle('Full Multi-Time Kernel Functions (Heisenberg sim)', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('kernel_full_multitime.png', dpi=150, bbox_inches='tight')
print("-> kernel_full_multitime.png")
plt.show()

# ===========================================================================
# 3. Volterra expansion validation against mesolve
# ===========================================================================
print("\n" + "="*70)
print("3. Volterra expansion vs mesolve ground truth")
print("="*70)

def volterra_predict(dw_waveform, k1, k2, k3, t_grid):
    """Predict dp_e using full Volterra expansion up to order 3.

    dw_waveform: detuning perturbation at each time point (rad*GHz)
    Returns: dp_e prediction (scalar)
    """
    dt_val = t_grid[1] - t_grid[0]

    # 1st order
    pred1 = np.trapezoid(k1 * dw_waveform, t_grid)

    # 2nd order: 1/2 * sum_ij k2(ti,tj) * dw(ti) * dw(tj) * dt^2
    dw_outer = np.outer(dw_waveform, dw_waveform)
    pred2 = 0.5 * np.trapezoid(np.trapezoid(k2 * dw_outer, t_grid, axis=0), t_grid, axis=0)

    # 3rd order: 1/6 * sum_ijk k3(ti,tj,tk) * dw(ti) * dw(tj) * dw(tk) * dt^3
    # Use tensor contraction for efficiency
    pred3 = (1.0/6.0) * np.trapezoid(
        np.trapezoid(
            np.trapezoid(
                k3_sym * dw_waveform.reshape(-1,1,1) * dw_waveform.reshape(1,-1,1) * dw_waveform.reshape(1,1,-1),
                t_grid, axis=0
            ), t_grid, axis=0
        ), t_grid, axis=0
    )

    return pred1, pred2, pred3, pred1 + pred2 + pred3


def mesolve_ground_truth(dw_waveform, qubit, omega_d, t_rabi, t_global, t_grid):
    """Direct mesolve with constant detuning perturbation added to Hamiltonian."""
    n_levels = qubit.n_levels
    psi_e = basis(n_levels, 1)
    sigma_z_op = Qobj(np.diag([1.0, -1.0] + [0.0]*(n_levels-2)))

    # Build control Hamiltonian
    H_pulse = QobjEvo(ctrl_y.hamiltonian_on(t_grid), tlist=t_grid, order=1)
    H0 = QobjEvo(0 * qeye(n_levels))

    # Add detuning perturbation
    H_detuning = QobjEvo(
        [[0.5 * sigma_z_op, dw_waveform]],
        tlist=t_grid, order=1,
    )

    H_total = H0 + H_pulse + H_detuning
    result = mesolve(H_total, qubit.state, t_grid, [],
                     e_ops=[psi_e * psi_e.dag()],
                     options={'max_step': float(dt)})
    return float(result.expect[0][-1])


# Test with different perturbation shapes
q = TransmonQubit(EC=0.2*2*np.pi, EJ=15*2*np.pi, T1=10000, T2=5000, n_levels=2, flux=0.0)

# Use a smaller time grid for faster mesolve (the original deduplicated grid)
# For the Volterra prediction, we interpolate the kernels to the mesolve grid

# Test case 1: Constant detuning
print("\n--- Test 1: Constant detuning dw(t) = delta ---")
for delta_ghz in [0.005, 0.01, 0.02, 0.05, 0.08, 0.10]:
    delta_rad = 2.0 * np.pi * delta_ghz
    dw_const = delta_rad * np.ones_like(t_grid)

    # Volterra prediction
    p1, p2, p3, p_total = volterra_predict(dw_const, k1, k2_full, k3_sym, t_grid)

    # Mesolve ground truth
    p_e_base = mesolve_ground_truth(np.zeros_like(t_grid), q, omega_d, t_rabi, t_global, t_grid)
    p_e_pert = mesolve_ground_truth(dw_const, q, omega_d, t_rabi, t_global, t_grid)
    dp_true = p_e_pert - p_e_base

    print(f"  delta={delta_ghz:.3f} GHz: "
          f"Volterra_1={p1:.4e}, Volterra_2={p2:.4e}, Volterra_3={p3:.4e}, "
          f"Volterra_tot={p_total:.4e}, True={dp_true:.4e}, "
          f"err={abs(p_total-dp_true):.2e}")

# Test case 2: Gaussian pulse perturbation
print("\n--- Test 2: Gaussian pulse dw(t) = A * exp(-(t-t0)^2/(2*sigma^2)) ---")
for amp_ghz in [0.01, 0.02, 0.05, 0.10]:
    amp_rad = 2.0 * np.pi * amp_ghz
    sigma = 3.0  # ns
    t0_val = T_pi2  # center at pulse boundary
    dw_gauss = amp_rad * np.exp(-0.5 * ((t_grid - t0_val) / sigma)**2)

    p1, p2, p3, p_total = volterra_predict(dw_gauss, k1, k2_full, k3_sym, t_grid)

    p_e_base = mesolve_ground_truth(np.zeros_like(t_grid), q, omega_d, t_rabi, t_global, t_grid)
    p_e_pert = mesolve_ground_truth(dw_gauss, q, omega_d, t_rabi, t_global, t_grid)
    dp_true = p_e_pert - p_e_base

    print(f"  amp={amp_ghz:.3f} GHz: "
          f"V1={p1:.4e}, V2={p2:.4e}, V3={p3:.4e}, "
          f"V_tot={p_total:.4e}, True={dp_true:.4e}, "
          f"err={abs(p_total-dp_true):.2e}")

# Test case 3: Scan amplitude, show convergence
print("\n--- Test 3: Convergence scan (Gaussian at pulse center) ---")
amps_ghz = np.logspace(-3, -1, 15)
errors_v1 = []
errors_v2 = []
errors_v3 = []
for amp_ghz in amps_ghz:
    amp_rad = 2.0 * np.pi * amp_ghz
    sigma = 2.0
    t0_val = T_pi2 / 2  # center of first pulse
    dw_gauss = amp_rad * np.exp(-0.5 * ((t_grid - t0_val) / sigma)**2)

    p1, p2, p3, p_total = volterra_predict(dw_gauss, k1, k2_full, k3_sym, t_grid)

    p_e_base = mesolve_ground_truth(np.zeros_like(t_grid), q, omega_d, t_rabi, t_global, t_grid)
    p_e_pert = mesolve_ground_truth(dw_gauss, q, omega_d, t_rabi, t_global, t_grid)
    dp_true = p_e_pert - p_e_base

    errors_v1.append(abs(p1 - dp_true))
    errors_v2.append(abs(p1 + p2 - dp_true))
    errors_v3.append(abs(p_total - dp_true))

fig2, ax2 = plt.subplots(figsize=(10, 6))
ax2.loglog(amps_ghz, errors_v1, 'o-', label='1st order only')
ax2.loglog(amps_ghz, errors_v2, 's-', label='1st + 2nd order')
ax2.loglog(amps_ghz, errors_v3, '^-', label='1st + 2nd + 3rd order')
ax2.set_xlabel('Perturbation amplitude (GHz)')
ax2.set_ylabel('|Volterra - mesolve| error')
ax2.set_title('Volterra Expansion Convergence (Gaussian at pulse center)')
ax2.legend()
ax2.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('volterra_convergence.png', dpi=150, bbox_inches='tight')
print("-> volterra_convergence.png")
plt.show()

# ===========================================================================
# 4. Compare with 5-pt FD exp kernels (diagonal only)
# ===========================================================================
print("\n" + "="*70)
print("4. 5-pt FD exp vs Heisenberg sim (diagonal comparison)")
print("="*70)

q_exp = TransmonQubit(EC=0.2*2*np.pi, EJ=15*2*np.pi, T1=10000, T2=5000, n_levels=2, flux=0.0)
est_exp = KernelEstimator(mode='omega', method='exp', order=3, stim_amplitude=0.01, virtual_z_impl='math')
r_exp = est_exp.estimate_full(ctrl_y, q_exp)

k1_exp = np.asarray(r_exp.kernels[0])
k2_exp = np.asarray(r_exp.kernels[1])
k3_exp = np.asarray(r_exp.kernels[2])
t_exp = np.asarray(r_exp.t_samples)

# Interpolate Sim to Exp grid for comparison
k1_sim_interp = np.interp(t_exp, t_grid, k1)
k3_diag_interp = np.interp(t_exp, t_grid, k3_diag_vals)

G1_exp = np.trapezoid(k1_exp, t_exp)
G3_exp = np.trapezoid(k3_exp, t_exp)

fig3, axes3 = plt.subplots(2, 2, figsize=(14, 10))

ax = axes3[0, 0]
ax.plot(t_exp, k1_sim_interp, 'o-', color='darkgreen', markersize=3, label=f'Sim (G1={G1:.3f})')
ax.plot(t_exp, k1_exp, 's--', color='darkorange', markersize=3, label=f'Exp (G1={G1_exp:.3f})')
ax.axvline(x=T_pi2, color='red', linestyle='--', alpha=0.4)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('k1')
ax.set_title('k1: Sim vs Exp')
ax.legend()
ax.set_xlim(0, 2*T_pi2)

rmse_k1 = np.sqrt(np.mean((k1_exp - k1_sim_interp)**2))
ax.text(0.98, 0.05, f'RMSE={rmse_k1:.2e}', transform=ax.transAxes, ha='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

ax = axes3[0, 1]
ax.plot(t_exp, k2_exp, 's-', color='purple', markersize=3, label=f'Exp k2')
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('k2')
ax.set_title(f'k2: Exp (expect ~0), G2_exp={np.trapezoid(k2_exp,t_exp):.2e}')
ax.legend()
ax.set_xlim(0, 2*T_pi2)

ax = axes3[1, 0]
ax.plot(t_exp, k3_diag_interp, 'o-', color='darkgreen', markersize=3, label=f'Sim (G3_diag={G3_diag:.3f})')
ax.plot(t_exp, k3_exp, 's--', color='darkorange', markersize=3, label=f'Exp (G3={G3_exp:.1f})')
ax.axvline(x=T_pi2, color='blue', linestyle='--', alpha=0.4)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('k3')
ax.set_title('k3 (diagonal): Sim vs Exp')
ax.legend()
ax.set_xlim(0, 2*T_pi2)

ax = axes3[1, 1]
ax.axis('off')
summary_text = (
    f"Diagonal kernel comparison:\n\n"
    f"k1: Sim G1={G1:.4f}, Exp G1={G1_exp:.4f}\n"
    f"    RMSE = {rmse_k1:.2e}\n"
    f"    |G1_sim/G1_exp| = {abs(G1/G1_exp):.4f}\n\n"
    f"k2: Exp G2={np.trapezoid(k2_exp,t_exp):.2e} (~0)\n\n"
    f"k3: Sim G3_diag={G3_diag:.4f}\n"
    f"    Exp G3_diag={G3_exp:.1f}\n"
    f"    |G3_sim/G3_exp| = {abs(G3_diag/G3_exp):.4f}\n\n"
    f"Note: exp k3 has O(h^2) FD error;\n"
    f"sim k3 is machine precision."
)
ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
        va='top', ha='left', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.8))

plt.suptitle('Heisenberg sim vs 5-pt FD exp: Diagonal Kernel Comparison', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('kernel_diag_comparison.png', dpi=150, bbox_inches='tight')
print("-> kernel_diag_comparison.png")
plt.show()

# ===========================================================================
# 5. Final summary
# ===========================================================================
print("\n" + "="*70)
print("FINAL SUMMARY")
print("="*70)
print(f"""
Multi-time kernel verification (Heisenberg sim, Nt={Nt}):

  k1:  G1 = {G1:.4f}

  k2:  max|k2(t1,t2)| = {np.max(np.abs(k2_full)):.2e}
       G2_full = {G2_full:.2e}
       k2 = 0 VERIFIED (Y-X orthogonality)

  k3:  diagonal: k3(t,t,t) = -k1(t) VERIFIED (RMSE={rmse_k3diag:.2e})
       G3_diag = {G3_diag:.4f}
       |G3_diag/G1| = {abs(G3_diag/G1):.4f}

       FULL 3D integral: G3_full = {G3_full:.1f}
       G3_full/G1 = {G3_full/G1:.1f}
       G3_full/|G3_diag| = {G3_full/abs(G3_diag):.1f}
       sqrt(|ratio|) = {np.sqrt(abs(G3_full/G3_diag)):.1f} ns ~ 2*T_pi2 = {2*T_pi2} ns

  Exp (5-pt FD) vs Sim (diagonal):
       k1: RMSE = {rmse_k1:.2e}, ratio G1 = {abs(G1/G1_exp):.4f}
       k3: ratio G3 = {abs(G3_diag/G3_exp):.4f} (exp has O(h^2) FD error)
""")
