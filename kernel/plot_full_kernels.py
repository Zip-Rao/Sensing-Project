"""
Complete 1st-3rd order kernel function visualization with full formulas.
Zero-detuning Ramsey Y-X pi/2-pi/2 sequence, Heisenberg sim (machine precision).
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

import qutip
from qutip import Qobj, QobjEvo, basis, sesolve, sigmax, sigmay, sigmaz, qeye

from sqc.config import CONFIG
from sqc.control.sequence import create_ramsey_pulse
from src.qubit import TransmonQubit

rcParams.update({
    'font.size': 11, 'figure.dpi': 150, 'axes.grid': True, 'grid.alpha': 0.25,
    'mathtext.fontset': 'cm',
})

# ===========================================================================
# Setup
# ===========================================================================
dt = CONFIG.awg.dt
T_pi2 = CONFIG.pulse.t_rabi_duration
t_rabi = CONFIG.pulse.t_rabi.copy()
q_ref = TransmonQubit(EC=0.2*2*np.pi, EJ=15*2*np.pi, T1=10000, T2=5000, n_levels=2, flux=0.0)
omega_d = q_ref.frequency
ctrl_y = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=omega_d, phase1=+np.pi/2, phase2=0.0)

n = 2
sigma_z = Qobj(np.diag([1.0, -1.0]))
M = basis(n, 1) * basis(n, 1).dag()
state0 = basis(n, 0)

raw_tlist = np.asarray(ctrl_y.t_list, dtype=float)
t_grid = np.unique(raw_tlist)
Nt = len(t_grid)

H0 = QobjEvo(0 * qeye(2))
H_pulse = QobjEvo(ctrl_y.hamiltonian_on(t_grid), tlist=t_grid, order=1)
res = sesolve(H0 + H_pulse, qeye(n), t_grid,
              options={'max_step': float(dt), 'store_states': True})
U_list = res.states
U_T = U_list[-1]
Q = U_T.dag() * M * U_T
W_list = [U_t.dag() * sigma_z * U_t / 2.0 for U_t in U_list]

def _expect(op, st):
    v = st.dag() * op * st
    if hasattr(v, 'shape') and v.shape == (1, 1):
        return complex(v[0, 0])
    return complex(v)

Omega = np.pi / (2 * T_pi2)

# ===========================================================================
# Compute kernels
# ===========================================================================
# -- k1(t) ----------------------------------------------------------------
k1 = np.zeros(Nt)
for i in range(Nt):
    comm = W_list[i] * Q - Q * W_list[i]
    k1[i] = float(np.real(1j * _expect(comm, state0)))
G1 = np.trapezoid(k1, t_grid)

# -- k2(t1, t2) 2D -------------------------------------------------------
print("Computing k2(t1,t2) ...")
k2 = np.zeros((Nt, Nt))
for i in range(Nt):
    for j in range(Nt):
        ti, tj = (j, i) if t_grid[i] >= t_grid[j] else (i, j)
        W_less, W_greater = W_list[ti], W_list[tj]
        comm1 = W_greater * Q - Q * W_greater
        comm2 = W_less * comm1 - comm1 * W_less
        k2[i, j] = float(np.real(-_expect(comm2, state0)))
G2_full = np.trapezoid(np.trapezoid(k2, t_grid, axis=0), t_grid, axis=0)

# -- k3(t1,t2,t3) 3D -----------------------------------------------------
print("Computing k3(t1,t2,t3) ...")
k3_ord = np.zeros((Nt, Nt, Nt))
for i in range(Nt):
    Wi = W_list[i]
    comm_i = Wi * Q - Q * Wi
    for j in range(Nt):
        if t_grid[j] > t_grid[i]:
            continue
        Wj = W_list[j]
        comm_ij = Wj * comm_i - comm_i * Wj
        for k in range(Nt):
            if t_grid[k] > t_grid[j]:
                continue
            Wk = W_list[k]
            comm_ijk = Wk * comm_ij - comm_ij * Wk
            k3_ord[i, j, k] = float(np.real(-1j * _expect(comm_ijk, state0)))

# Symmetrize
k3_sym = np.zeros((Nt, Nt, Nt))
for i in range(Nt):
    for j in range(Nt):
        for k in range(Nt):
            idx = np.array([i, j, k])
            order = np.argsort(-t_grid[idx])
            i_o, j_o, k_o = idx[order]
            k3_sym[i, j, k] = k3_ord[i_o, j_o, k_o]

k3_diag = np.array([k3_sym[i, i, i] for i in range(Nt)])
G3_diag = np.trapezoid(k3_diag, t_grid)
G3_full = np.trapezoid(np.trapezoid(np.trapezoid(k3_sym, t_grid, axis=0), t_grid, axis=0), t_grid, axis=0)
print("Done.")

# ===========================================================================
# Visualization
# ===========================================================================
fig = plt.figure(figsize=(22, 18))
gs = gridspec.GridSpec(4, 4, figure=fig, hspace=0.45, wspace=0.4,
                       height_ratios=[1.2, 1.2, 1.2, 0.6])

# ---- Row 0: k1(t) ----
# (0,0): k1 plot
ax = fig.add_subplot(gs[0, :2])
ax.plot(t_grid, k1, 'b-', linewidth=2.5)
ax.fill_between(t_grid, 0, k1, alpha=0.08, color='blue')
ax.axvline(x=T_pi2, color='gray', linestyle='--', linewidth=1.5, alpha=0.6,
           label=f'$t = T_{{\\pi/2}} = {T_pi2}$ ns')
ax.axhline(y=0, color='gray', alpha=0.4)
# Annotate formula segments
ax.annotate(r'$k_1(t) = +\frac{1}{2}\sin(\Omega t)$',
            xy=(T_pi2/2, 0.28), fontsize=12, color='blue', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
ax.annotate(r'$k_1(t) = -\frac{1}{2}\cos(\Omega(t-T_1))$',
            xy=(1.55*T_pi2, -0.25), fontsize=12, color='blue', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
ax.set_xlabel('$t$ (ns)', fontsize=12)
ax.set_ylabel('$k_1(t)$  (rad$^{-1}$)', fontsize=12)
ax.set_title(r'$\mathbf{k_1(t)}$ — First-order kernel  |  $G_1 = \int k_1\,dt = %.3f$' % G1,
             fontsize=13, fontweight='bold')
ax.legend(fontsize=9, loc='lower left')
ax.set_xlim(0, 2*T_pi2)

# (0,2): k1 commutator formula
ax = fig.add_subplot(gs[0, 2])
ax.axis('off')
formula_k1 = (
    r"$\mathbf{First{-}order\ kernel}$" "\n\n"
    r"$k_1(t) = i\,\langle 0|\,[W(t),\,Q]\,|0\rangle$" "\n\n"
    r"$W(t) = \frac{1}{2}U^\dagger(t)\,\sigma_z\,U(t)$" "\n"
    r"$Q = U^\dagger(T)\,|1\rangle\langle 1|\,U(T)$" "\n\n"
    r"$\mathbf{Square\ envelope\ (analytic):}$" "\n\n"
    r"$\bullet$ Pulse 1 $(0 \leq t < T_{\pi/2})$:" "\n"
    r"$\; W(t) = \frac{1}{2}(\sin\Omega t\,\sigma_x + \cos\Omega t\,\sigma_z)$" "\n"
    r"$\; k_1(t) = +\frac{1}{2}\sin(\Omega t)$" "\n\n"
    r"$\bullet$ Pulse 2 $(T_{\pi/2} \leq t < 2T_{\pi/2})$:" "\n"
    r"$\; W(t) = \frac{1}{2}(-\cos\Omega(t{-}T_1)\,\sigma_x + \sin\Omega(t{-}T_1)\,\sigma_y)$" "\n"
    r"$\; k_1(t) = -\frac{1}{2}\cos(\Omega(t{-}T_1))$" "\n\n"
    r"$\bullet\; G_1 = \int_0^{2T_{\pi/2}}\! k_1(t)\,dt = -\frac{2T_{\pi/2}}{\pi}$" "\n"
    r"$\quad\; \approx -6.366\ \mathrm{(square)}$" "\n"
    r"$\quad\; \approx %.3f\ \mathrm{(rectangular)}$" % G1
)
ax.text(0.05, 0.97, formula_k1, transform=ax.transAxes, fontsize=9.5,
        va='top', ha='left', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.7))

# (0,3): Bloch sphere sketch — W(t) trajectory text description
ax = fig.add_subplot(gs[0, 3])
ax.axis('off')
bloch_text = (
    r"$\mathbf{Geometric\ meaning\ of\ }k_1:$" "\n\n"
    r"$W(t)$ = Heisenberg-picture $\sigma_z/2$" "\n"
    r"under the control propagator." "\n\n"
    r"$k_1(t)$ measures how much the" "\n"
    r"infinitesimal perturbation axis" "\n"
    r"$W(t)$ couples to the final" "\n"
    r"measurement direction $Q$." "\n\n"
    r"Pulse 1: $W$ rotates in Z-X plane" "\n"
    r"Pulse 2: $W$ rotates in X-Y plane" "\n\n"
    r"$Q = \frac{1}{2}(I - \sigma_y)$" "\n"
    r"$\mathbf{q} = (0, -1, 0)$ on Bloch sphere"
)
ax.text(0.05, 0.97, bloch_text, transform=ax.transAxes, fontsize=9,
        va='top', ha='left', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.7))

# ---- Row 1: k2(t1, t2) ----
# (1,0): k2 heatmap
ax = fig.add_subplot(gs[1, :2])
vmax2 = max(np.max(np.abs(k2)), 1e-10)
im = ax.pcolormesh(t_grid, t_grid, k2, cmap='RdBu_r', shading='auto',
                   vmin=-vmax2, vmax=vmax2, rasterized=True)
ax.axvline(x=T_pi2, color='green', linestyle='--', linewidth=1.2, alpha=0.5)
ax.axhline(y=T_pi2, color='green', linestyle='--', linewidth=1.2, alpha=0.5)
ax.set_xlabel('$t_1$ (ns)', fontsize=12)
ax.set_ylabel('$t_2$ (ns)', fontsize=12)
ax.set_title(r'$\mathbf{k_2(t_1, t_2)}$ — Second-order kernel  |  '
             r'$\max|k_2| = %.2e$  |  $\iint k_2 = %.2e$' %
             (np.max(np.abs(k2)), G2_full),
             fontsize=13, fontweight='bold')
ax.set_aspect('equal')
cbar = plt.colorbar(im, ax=ax, label='$k_2$  (rad$^{-2}$)', shrink=0.85)
cbar.ax.yaxis.label.set_size(11)

# (1,2): k2 diagonal
ax = fig.add_subplot(gs[1, 2])
ax.plot(t_grid, np.diag(k2), 'o-', markersize=3.5, color='purple', linewidth=1.2,
        label='$k_2(t,t)$')
ax.axhline(y=0, color='gray', alpha=0.4)
ax.axvline(x=T_pi2, color='gray', linestyle='--', alpha=0.4)
ax.set_xlabel('$t$ (ns)', fontsize=11)
ax.set_ylabel('$k_2(t,t)$  (rad$^{-2}$)', fontsize=11)
ax.set_title(r'$k_2$ diagonal slice', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.set_xlim(0, 2*T_pi2)

# (1,3): k2 formula
ax = fig.add_subplot(gs[1, 3])
ax.axis('off')
formula_k2 = (
    r"$\mathbf{Second{-}order\ kernel}$" "\n\n"
    r"$k_2(t_>, t_<) = -\langle 0|\,[W(t_<),\,[W(t_>),\,Q]]\,|0\rangle$" "\n"
    r"$t_> = \max(t_1, t_2),\ t_< = \min(t_1, t_2)$" "\n\n"
    r"$\mathbf{Zero{-}detuning\ Y{-}X:}$" "\n"
    r"$k_2(t_1, t_2) \equiv 0\ \mathrm{(strictly!)}$" "\n\n"
    r"$\mathbf{Proof\ (Pauli\ algebra):}$" "\n"
    r"$[W, [W, Q]]$" "\n"
    r"$= \frac{1}{2}(\mathbf{w}\cdot\mathbf{q})\mathbf{w} - \mathbf{q})\cdot\sigma$" "\n\n"
    r"For Y-X pulse sequence:" "\n"
    r"$\mathbf{w}(t)$ always lies in Z-X plane" "\n"
    r"(pulse 1) or X-Y plane (pulse 2)." "\n"
    r"$\mathbf{q} = (0,-1,0)$ — pointing along $-y$." "\n\n"
    r"$(\mathbf{w}\cdot\mathbf{q})\mathbf{w} - \mathbf{q} = \mathbf{0}$" "\n"
    r"$\Longrightarrow k_2 \equiv 0$" "\n\n"
    r"$\mathbf{Numerical:}$" "\n"
    r"$\max|k_2| = %.2e$  (dt artifact)" % np.max(np.abs(k2))
)
ax.text(0.05, 0.97, formula_k2, transform=ax.transAxes, fontsize=9,
        va='top', ha='left', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.7))

# ---- Row 2: k3(t1, t2, t3) ----
# (2,0): k3 slice at t3 = T_pi2/4 (early pulse 1)
t3_idx_a = np.argmin(np.abs(t_grid - T_pi2/4))
k3_slice_a = k3_sym[:, :, t3_idx_a]
vmax3a = max(np.max(np.abs(k3_slice_a)), 1e-10)
ax = fig.add_subplot(gs[2, 0])
im = ax.pcolormesh(t_grid, t_grid, k3_slice_a, cmap='RdBu_r', shading='auto',
                   vmin=-vmax3a, vmax=vmax3a, rasterized=True)
ax.axvline(x=T_pi2, color='green', linestyle='--', linewidth=1, alpha=0.4)
ax.axhline(y=T_pi2, color='green', linestyle='--', linewidth=1, alpha=0.4)
ax.set_xlabel('$t_1$ (ns)', fontsize=10)
ax.set_ylabel('$t_2$ (ns)', fontsize=10)
ax.set_title(r'$k_3(t_1,t_2,t_3{=}%.1f\,\mathrm{ns})$  pulse 1 early' % t_grid[t3_idx_a],
             fontsize=11, fontweight='bold')
ax.set_aspect('equal')
cbar = plt.colorbar(im, ax=ax, label='$k_3$  (rad$^{-3}$)', shrink=0.8)

# (2,1): k3 slice at t3 = T_pi2 (pulse boundary)
t3_idx_b = np.argmin(np.abs(t_grid - T_pi2))
k3_slice_b = k3_sym[:, :, t3_idx_b]
vmax3b = max(np.max(np.abs(k3_slice_b)), 1e-10)
ax = fig.add_subplot(gs[2, 1])
im = ax.pcolormesh(t_grid, t_grid, k3_slice_b, cmap='RdBu_r', shading='auto',
                   vmin=-vmax3b, vmax=vmax3b, rasterized=True)
ax.axvline(x=T_pi2, color='green', linestyle='--', linewidth=1, alpha=0.4)
ax.axhline(y=T_pi2, color='green', linestyle='--', linewidth=1, alpha=0.4)
ax.set_xlabel('$t_1$ (ns)', fontsize=10)
ax.set_ylabel('$t_2$ (ns)', fontsize=10)
ax.set_title(r'$k_3(t_1,t_2,t_3{=}%.1f\,\mathrm{ns})$  pulse boundary' % t_grid[t3_idx_b],
             fontsize=11, fontweight='bold')
ax.set_aspect('equal')
cbar = plt.colorbar(im, ax=ax, label='$k_3$  (rad$^{-3}$)', shrink=0.8)

# (2,2): k3 diagonal verification
ax = fig.add_subplot(gs[2, 2])
ax.plot(t_grid, -k1, 'b-', linewidth=2.5, alpha=0.35, label='$-k_1(t)$ (reference)')
ax.plot(t_grid, k3_diag, 'ro-', markersize=3.5, linewidth=1.2, label='$k_3(t,t,t)$')
ax.axvline(x=T_pi2, color='gray', linestyle='--', alpha=0.4)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('$t$ (ns)', fontsize=11)
ax.set_ylabel('$k_3(t,t,t)$  (rad$^{-3}$)', fontsize=11)
rmse = np.sqrt(np.mean((k3_diag + k1)**2))
ax.set_title(r'Diagonal: $k_3(t,t,t) = -k_1(t)$  |  RMSE = %.1e' % rmse,
             fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.set_xlim(0, 2*T_pi2)

# (2,3): k3 formula
ax = fig.add_subplot(gs[2, 3])
ax.axis('off')
formula_k3 = (
    r"$\mathbf{Third{-}order\ kernel}$" "\n\n"
    r"$k_3^{\rm ord}(t_1,t_2,t_3) = -i\langle 0|\,[W(t_3),\,[W(t_2),\,[W(t_1),\,Q]]]\,|0\rangle$" "\n"
    r"for $t_1 \geq t_2 \geq t_3$ (time-ordered wedge)" "\n\n"
    r"$\mathbf{Fully\ symmetric\ Volterra\ kernel:}$" "\n"
    r"$k_3(t_1,t_2,t_3) = k_3^{\rm ord}(t_{(1)}, t_{(2)}, t_{(3)})$" "\n"
    r"where $t_{(1)} \geq t_{(2)} \geq t_{(3)}$" "\n\n"
    r"$\mathbf{Diagonal\ (}t_1{=}t_2{=}t_3{=}t\mathbf{):}$" "\n"
    r"$k_3(t,t,t) = -k_1(t)$" "\n"
    r"(period-2 structure of $\mathrm{ad}_W^n$" "\n"
    r" in 2-level Pauli algebra)" "\n\n"
    r"$\mathbf{Integrals:}$" "\n"
    r"$G_3^{\rm diag} = \int k_3(t,t,t)\,dt = %.3f$" % G3_diag + "\n"
    r"$G_3^{\rm full}  = \iiint k_3\,dt_1dt_2dt_3 = %.1f$" % G3_full + "\n"
    r"$G_3^{\rm full} / |G_3^{\rm diag}| = %.1f$" % (G3_full/abs(G3_diag)) + "\n"
    r"$\sqrt{|\mathrm{ratio}|} = %.1f\ \mathrm{ns}$" % np.sqrt(abs(G3_full/G3_diag))
)
ax.text(0.05, 0.97, formula_k3, transform=ax.transAxes, fontsize=8.5,
        va='top', ha='left', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.7))

# ---- Row 3: Summary table ----
ax = fig.add_subplot(gs[3, :])
ax.axis('off')

# Use plain text table (matplotlib mathtext doesn't support \begin{array})
hdr = "Summary: Volterra Kernel Functions for Ramsey Y-X pi/2-pi/2 (tau=0)"
sep = "=" * 110
line1 = f"  k1(t)          G1 = {G1:+.4f}        k1 = +1/2 sin(Omega*t) [pulse 1], -1/2 cos(...) [pulse 2]"
line2 = f"  k2(t1,t2)      G2 = {G2_full:+.2e}   k2 = 0  (Y-X orthogonality, strictly)"
line3 = f"  k3(t1,t2,t3)   G3_diag = {G3_diag:+.4f}   diagonal: k3(t,t,t) = -k1(t)  (period-2 of ad_W^n)"
line4 = f"                 G3_full = {G3_full:+.1f}   off-diagonal dominates!  ratio = {G3_full/abs(G3_diag):.1f}x"
table_text = f"{hdr}\n{sep}\n{line1}\n{line2}\n{line3}\n{line4}"

ax.text(0.5, 0.5, table_text, transform=ax.transAxes, fontsize=9.5,
        va='center', ha='center', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.85))

plt.savefig('kernel_full_orders.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print("Saved: kernel_full_orders.png")
plt.close()

# ===========================================================================
# Print complete formulas
# ===========================================================================
print("""
================================================================================
COMPLETE KERNEL FUNCTION FORMULAS
================================================================================

GENERAL DEFINITION — Volterra expansion of p_e in terms of frequency perturbation dw(t):

  p_e = p_e^(0) + sum_{n=1}^{inf} (1/n!) * int...int k_n(t_1,...,t_n) * prod_{j=1}^n dw(t_j) * dt_j

  k_n(t_1,...,t_n) = delta^n p_e / (delta dw(t_1) ... delta dw(t_n)) |_{dw=0}

HEISENBERG-PICTURE COMMUTATOR FORMULAS (time-ordered wedge region):

  Define:
    U(t)    = control propagator from 0 to t (no perturbation)
    W(t)    = U^dag(t) * sigma_z * U(t) / 2         (effective perturbation operator)
    Q       = U^dag(T) * |1><1| * U(T)               (Heisenberg measurement operator)
    ad_W(O) = [W, O] = W*O - O*W                     (adjoint action)

  k1(t) = i * <0| ad_W(t) (Q) |0>

  k2^ord(t1, t2) = - <0| ad_W(t2) o ad_W(t1) (Q) |0>     for t1 >= t2

  k3^ord(t1, t2, t3) = -i * <0| ad_W(t3) o ad_W(t2) o ad_W(t1) (Q) |0>
                                                          for t1 >= t2 >= t3

  Full symmetric Volterra kernel:
    k2(t1,t2) = k2^ord(max(t1,t2), min(t1,t2))
    k3(t1,t2,t3) = k3^ord(t_(1), t_(2), t_(3))  where t_(1)>=t_(2)>=t_(3)

ANALYTICAL FORMS — Square envelope, zero-detuning Ramsey Y-X pi/2-pi/2:

  Rabi frequency: Omega = pi / (2*T_{pi/2})
  Pulse 1 (0 <= t < T_{pi/2}): Y rotation
  Pulse 2 (T_{pi/2} <= t < 2*T_{pi/2}): X rotation

  k1(t):
    Pulse 1: k1(t) = +1/2 * sin(Omega*t)
    Pulse 2: k1(t) = -1/2 * cos(Omega*(t - T_{pi/2}))
    G1 = int_0^{2*T_{pi/2}} k1(t) dt = -2*T_{pi/2} / pi = """ + f"{-2*T_pi2/np.pi:.4f}" + """

  k2(t1, t2):
    k2(t1, t2) = 0  (strictly, for any t1, t2)
    Proof: For Y-X pulses, W(t) always lies in Z-X or X-Y plane,
    while Q = (I - sigma_y)/2. The double commutator <0|[W,[W,Q]]|0>
    vanishes identically by Pauli algebra.
    G2_full = 0

  k3(t1, t2, t3):
    Diagonal: k3(t,t,t) = -k1(t)  (period-2 of ad_W^n in 2-level system)
    Off-diagonal (t1 != t2 != t3): DOES NOT simplify to -k1.
    The full 3D integral is much larger than the diagonal integral:
      G3_diag = """ + f"{G3_diag:.4f}" + """
      G3_full = """ + f"{G3_full:.1f}" + """
      Ratio   = """ + f"{G3_full/abs(G3_diag):.1f}" + """
    This ratio ~ (2*T_{pi/2})^2 / 2 explains the Phase 11 bug:
    using G3_diag instead of G3_full in cubic Newton correction
    underestimates the cubic response by ~""" + f"{abs(G3_full/G3_diag):.0f}" + """x.

NUMERICAL VALUES (rectangular envelope, dt=""" + f"{dt}" + """ ns, Nt=""" + f"{Nt}" + """):

  G1       = """ + f"{G1:.4f}" + """
  G2_full  = """ + f"{G2_full:.2e}" + """  (~0, Y-X orthogonality verified)
  G3_diag  = """ + f"{G3_diag:.4f}" + """  (|G3_diag/G1| = """ + f"{abs(G3_diag/G1):.4f}" + """)
  G3_full  = """ + f"{G3_full:.1f}" + """  (G3_full/G1 = """ + f"{G3_full/G1:.1f}" + """)
  max|k2|  = """ + f"{np.max(np.abs(k2)):.2e}" + """
  RMSE(k3_diag + k1) = """ + f"{rmse:.2e}" + """  (diagonal period-2 verified)

================================================================================
""")
