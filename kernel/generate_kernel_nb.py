"""Generate kernel_analytical_vs_numerical.ipynb using nbformat."""
import nbformat as nbf

nb = nbf.v4.new_notebook(
    metadata={
        "kernelspec": {"display_name": "qutip-env", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.13.0"},
    }
)

def md(source):
    nb.cells.append(nbf.v4.new_markdown_cell(source))

def code(source):
    nb.cells.append(nbf.v4.new_code_cell(source))

# ===== Cell 0: Title =====
md("""# First-Order and Third-Order Kernel Functions: Analytical vs Numerical

## Physical setup

**Zero-detuning orthogonal Ramsey Y-X $\\pi/2$-$\\pi/2$ sequence:**
- First pulse: $R_y(\\pi/2)$, phase $= +\\pi/2$
- Second pulse: $R_x(\\pi/2)$, phase $= 0$ (after $R_y(\\pi/2)$)

**Heisenberg-picture commutator formula** (from Dyson expansion):

$$k_n(t) = i^n \\langle 0| \\mathrm{ad}_W^n(Q) |0\\rangle$$

where:
- $W(t) = U_\\mathrm{ctrl}^\\dagger(t) \\sigma_z U_\\mathrm{ctrl}(t) / 2$ (effective perturbation operator)
- $Q = U_\\mathrm{ctrl}^\\dagger(T) |1\\rangle\\langle 1| U_\\mathrm{ctrl}(T)$ (Heisenberg measurement operator)
- $\\mathrm{ad}_W^n(Q) = [W, [W, \\ldots, [W, Q]\\ldots]]$ ($n$ nested commutators)

**Key results:**
- $k_1(t)$: single-peaked sinusoidal, $G_1 \\approx -6.37$ (square)
- $k_2(t)$: strictly zero at zero detuning (Y-X orthogonality)
- $k_3(t) = -k_1(t)$ in 2-level system, $|G_3/G_1| = 1$
""")

# ===== Cell 1: Imports =====
code("""import numpy as np
import matplotlib.pyplot as plt
import time, math
from matplotlib import rcParams
import qutip
from qutip import Qobj, QobjEvo, basis, sesolve, mesolve, sigmax, sigmay, sigmaz, qeye

from sqc.config import CONFIG
from sqc.control.sequence import create_ramsey_pulse
from sqc.reconstruction.kernel import KernelEstimator
from src.qubit import TransmonQubit

rcParams.update({'font.size': 12, 'figure.dpi': 120, 'axes.grid': True, 'grid.alpha': 0.3})

# Physical parameters
dt = CONFIG.awg.dt                    # 0.5 ns
T_pi2 = CONFIG.pulse.t_rabi_duration  # 10 ns
Omega_rabi = np.pi / (2 * T_pi2)      # Rabi frequency (rad/ns)

# Build Ramsey Y-X pulses
q_ref = TransmonQubit(EC=0.2*2*np.pi, EJ=15*2*np.pi, T1=10000, T2=5000, n_levels=2, flux=0.0)
omega_d = q_ref.frequency
t_rabi = CONFIG.pulse.t_rabi.copy()
t_global = CONFIG.pulse.t_global.copy()

ctrl_y  = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=omega_d, phase1=+np.pi/2, phase2=0.0)
ctrl_my = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=omega_d, phase1=-np.pi/2, phase2=0.0)

print(f'dt={dt} ns, T_pi2={T_pi2} ns, Omega_rabi={Omega_rabi:.4f} rad/ns')
print(f'omega_d={omega_d:.4f} rad*GHz (sweet spot)')
print(f'Pulse duration = {2*T_pi2} ns')
print(f'ctrl_y t_list has {len(ctrl_y.t_list)} time points')
""")

# ===== Cell 2: Analytical formulas =====
md(r"""## 1. Analytical Formulas (Square Envelope)

### 1a. Setup

Square Rabi envelope $\Omega(t) = \Omega = \pi/(2T_{\pi/2})$, duration $T_{\pi/2}$.

**Pulse 1** ($0 \le t < T_{\pi/2}$, Y rotation):
$$U_{\rm ctrl}(t,0) = R_y(\Omega t)$$
$$W(t) = \frac{1}{2}\bigl(\sin(\Omega t)\,\sigma_x + \cos(\Omega t)\,\sigma_z\bigr)$$

**Pulse 2** ($T_{\pi/2} \le t < 2T_{\pi/2}$, X rotation):
$$U_{\rm ctrl}(t,0) = R_x(\Omega(t-T_{\pi/2}))\,R_y(\pi/2)$$
$$W(t) = \frac{1}{2}\bigl(-\cos(\Omega(t-T_{\pi/2}))\,\sigma_x + \sin(\Omega(t-T_{\pi/2}))\,\sigma_y\bigr)$$

**Final measurement operator:**
$$U_{\rm ctrl}(T,0) = R_x(\pi/2)\,R_y(\pi/2)$$
$$Q = \frac{1}{2}(I + \sigma_x - \sigma_y)$$

### 1b. First-order kernel $k_1(t)$

$$k_1(t) = i\langle 0|[W(t), Q]|0\rangle = \begin{cases}
-\frac{1}{2}\sin(\Omega t) & t < T_{\pi/2} \\[6pt]
-\frac{1}{2}\cos(\Omega(t-T_{\pi/2})) & t > T_{\pi/2}
\end{cases}$$

$$G_1 = \int_0^{2T_{\pi/2}} k_1(t)\,dt = -\frac{2T_{\pi/2}}{\pi} \approx -6.366$$

### 1c. Second-order kernel $k_2(t)$

$$k_2(t) = -\langle 0|[W(t), [W(t), Q]]|0\rangle \equiv 0$$

**Algebraic proof:** For Y-X orthogonal pulses, $W(t)$ lies in Z-X plane (pulse 1)
or X-Y plane (pulse 2). The double commutator expectation in $|0\rangle$
vanishes by Pauli algebra. See `phase_10_kernel_theory_and_fix.md` for details.

### 1d. Third-order kernel $k_3(t)$

In a 2-level system, nested commutators have **period-2 structure**:
$$\mathrm{ad}_W^3(Q) = \mathrm{ad}_W^1(Q) = [W, Q]$$

Therefore:
$$k_3(t) = i^3\langle 0|\mathrm{ad}_W^3(Q)|0\rangle = -i\langle 0|[W,Q]|0\rangle = -k_1(t)$$

$$k_3(t) = \begin{cases}
+\frac{1}{2}\sin(\Omega t) & t < T_{\pi/2} \\[6pt]
+\frac{1}{2}\cos(\Omega(t-T_{\pi/2})) & t > T_{\pi/2}
\end{cases}$$

$$G_3 = -G_1 = +\frac{2T_{\pi/2}}{\pi} \approx +6.366, \quad |G_3/G_1| = 1$$
""")

# ===== Cell 3: Analytical evaluation =====
code(r"""# Evaluate analytical k1, k3 on a dense time grid
t_dense = np.linspace(0, 2*T_pi2, 500)

def analytical_k1_square(t, T1=T_pi2, Omega=Omega_rabi):
    '''Square-envelope k1(t) analytical solution.'''
    t = np.asarray(t)
    k1 = np.zeros_like(t)
    mask1 = (t >= 0) & (t < T1)
    k1[mask1] = -0.5 * np.sin(Omega * t[mask1])
    mask2 = (t >= T1) & (t < 2*T1)
    k1[mask2] = -0.5 * np.cos(Omega * (t[mask2] - T1))
    return k1

def analytical_k3_square(t, T1=T_pi2, Omega=Omega_rabi):
    '''Square-envelope k3(t) analytical solution: k3 = -k1 (2-level).'''
    return -analytical_k1_square(t, T1, Omega)

k1_analytic = analytical_k1_square(t_dense)
k3_analytic = analytical_k3_square(t_dense)

G1_analytic = np.trapezoid(k1_analytic, t_dense)
G3_analytic = np.trapezoid(k3_analytic, t_dense)

print(f'Analytical (square envelope, T_pi2={T_pi2} ns):')
print(f'  G1 = {G1_analytic:.4f}  (expected: -6.366)')
print(f'  G3 = {G3_analytic:.4f}  (expected: +6.366)')
print(f'  |G3/G1| = {abs(G3_analytic/G1_analytic):.4f}')
print(f'  max|k3 + k1| = {np.max(np.abs(k3_analytic + k1_analytic)):.2e}')

# Plot analytical kernels
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(t_dense, k1_analytic, 'b-', linewidth=2)
ax1.fill_between(t_dense, 0, k1_analytic, alpha=0.1, color='blue')
ax1.axvline(x=T_pi2, color='red', linestyle='--', alpha=0.5, label=f't = T_pi2 = {T_pi2} ns')
ax1.axhline(y=0, color='gray', alpha=0.4)
ax1.set_xlabel('t (ns)')
ax1.set_ylabel(r'$k_1(t)$ (rad$^{-1}$)')
ax1.set_title('Analytical $k_1(t)$ | Square Envelope')
ax1.legend()
ax1.set_xlim(0, 2*T_pi2)
ax1.annotate(r'$k_1=-\frac{1}{2}\sin(\Omega t)$', xy=(T_pi2/2, -0.25),
            fontsize=11, color='blue',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7))
ax1.annotate(r'$k_1=-\frac{1}{2}\cos(\Omega(t-T_1))$', xy=(T_pi2*1.5, -0.35),
            fontsize=11, color='blue',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7))

ax2.plot(t_dense, k3_analytic, 'r-', linewidth=2)
ax2.fill_between(t_dense, 0, k3_analytic, alpha=0.1, color='red')
ax2.axvline(x=T_pi2, color='blue', linestyle='--', alpha=0.5, label=f't = T_pi2 = {T_pi2} ns')
ax2.axhline(y=0, color='gray', alpha=0.4)
ax2.set_xlabel('t (ns)')
ax2.set_ylabel(r'$k_3(t)$ (rad$^{-3}$)')
ax2.set_title(r'Analytical $k_3(t) = -k_1(t)$ | Square Envelope')
ax2.legend()
ax2.set_xlim(0, 2*T_pi2)
ax2.annotate(r'$k_3=+\frac{1}{2}\sin(\Omega t)$', xy=(T_pi2/2, 0.25),
            fontsize=11, color='red',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.7))
ax2.annotate(r'$k_3=+\frac{1}{2}\cos(\Omega(t-T_1))$', xy=(T_pi2*1.5, 0.35),
            fontsize=11, color='red',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.7))

plt.suptitle('Analytical Kernel Functions (Square Envelope, Zero-Detuning Ramsey Y-X)',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('kernel_analytical_square.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved: kernel_analytical_square.png')
""")

# ===== Cell 4: Heisenberg sim =====
md(r"""## 2. Heisenberg Propagator (method='sim', order=3)

**Algorithm:**
1. `sesolve(qeye(n), t_grid)` $\to$ propagator $U(t)$ at all times
2. $W(t_j) = U^\dagger(t_j) \sigma_z U(t_j) / 2$
3. Evaluate nested commutators iteratively:
   - $k_1 = i \langle 0| [W, Q] |0\rangle$
   - $k_2 = -\langle 0| [W, [W, Q]] |0\rangle$
   - $k_3 = -i \langle 0| [W, [W, [W, Q]]] |0\rangle$

**Cost:** 1 sesolve + O($N_t \cdot d^3$) matrix multiplications.
**Precision:** Machine precision (no perturbation, no fitting, no FD).
**Any order** at the same cost.
""")

code("""print("Computing Heisenberg sim kernels (order=3)...")
q_sim = TransmonQubit(EC=0.2*2*np.pi, EJ=15*2*np.pi, T1=10000, T2=5000, n_levels=2, flux=0.0)
est_sim = KernelEstimator(mode='omega', method='sim', order=3, n_levels=2)

t0 = time.time()
result_sim_y = est_sim.estimate_full(ctrl_y, q_sim)
result_sim_my = est_sim.estimate_full(ctrl_my, q_sim)
t_sim = time.time() - t0

t_k = result_sim_y.t_samples
k1_sim = np.asarray(result_sim_y.kernels[0])
k2_sim = np.asarray(result_sim_y.kernels[1])
k3_sim = np.asarray(result_sim_y.kernels[2])

G1_sim = np.trapezoid(k1_sim, t_k)
G2_sim = np.trapezoid(k2_sim, t_k)
G3_sim = np.trapezoid(k3_sim, t_k)

print(f'Heisenberg sim completed in {t_sim:.2f}s')
print(f'  G1 = {G1_sim:.4f}')
print(f'  G2 = {G2_sim:.2e}  (expected ~0)')
print(f'  G3 = {G3_sim:.4f}')
print(f'  |G3/G1| = {abs(G3_sim/G1_sim):.4f}')
print(f'  |G2/G1| = {abs(G2_sim/G1_sim):.2e}')
""")

# ===== Cell 5: Analytical vs Heisenberg comparison =====
code("""# Compare analytical (square) vs Heisenberg sim (rectangular envelope)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

# k1 comparison
ax1.plot(t_dense, k1_analytic, 'b-', linewidth=2, alpha=0.4, label='Analytical (square)')
ax1.plot(t_k, k1_sim, 'o-', color='darkgreen', markersize=4, linewidth=1.5,
         label=f'Heisenberg sim (G1={G1_sim:.3f})')
ax1.axvline(x=T_pi2, color='red', linestyle='--', alpha=0.4)
ax1.axhline(y=0, color='gray', alpha=0.4)
ax1.set_xlabel('t (ns)')
ax1.set_ylabel('k1 (rad^-1)')
ax1.set_title('k1: Analytical vs Heisenberg sim')
ax1.legend()
ax1.set_xlim(0, 2*T_pi2)

# k3 comparison
ax2.plot(t_dense, k3_analytic, 'r-', linewidth=2, alpha=0.4, label='Analytical (square)')
ax2.plot(t_k, k3_sim, 'o-', color='darkorange', markersize=4, linewidth=1.5,
         label=f'Heisenberg sim (G3={G3_sim:.3f})')
ax2.axvline(x=T_pi2, color='blue', linestyle='--', alpha=0.4)
ax2.axhline(y=0, color='gray', alpha=0.4)
ax2.set_xlabel('t (ns)')
ax2.set_ylabel('k3 (rad^-3)')
ax2.set_title('k3: Analytical vs Heisenberg sim')
ax2.legend()
ax2.set_xlim(0, 2*T_pi2)

# k3 = -k1 verification (scatter)
ax3.scatter(k1_sim, -k3_sim, c=t_k, cmap='viridis', s=30, alpha=0.8)
ax3.plot([-0.65, 0.15], [-0.65, 0.15], 'k--', alpha=0.3, label='y = x (ideal)')
ax3.set_xlabel('k1 (sim)')
ax3.set_ylabel('-k3 (sim)')
ax3.set_title('Verification: k3 = -k1 (2-level)')
ax3.legend()
ax3.set_aspect('equal')
cbar = plt.colorbar(ax3.collections[0], ax=ax3, label='t (ns)')

plt.suptitle('Analytical (Square) vs Heisenberg sim (Rectangular Envelope)',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('kernel_analytic_vs_sim.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved: kernel_analytic_vs_sim.png')

# Compute RMSE between analytical and sim
k1_analytic_at_sim = analytical_k1_square(t_k)
k3_analytic_at_sim = analytical_k3_square(t_k)
rmse_k1 = np.sqrt(np.mean((k1_sim - k1_analytic_at_sim)**2))
rmse_k3 = np.sqrt(np.mean((k3_sim - k3_analytic_at_sim)**2))
print(f'RMSE(k1, analytic vs sim): {rmse_k1:.2e}')
print(f'RMSE(k3, analytic vs sim): {rmse_k3:.2e}')
print('(Small differences due to dt=0.5ns discretization of the rectangular envelope)')
""")

# ===== Cell 6: 5-pt FD explanation =====
md(r"""## 3. 5-Point Finite Difference (method='exp', order=3)

Samples $p_e$ at 5 Virtual Z amplitudes $\phi_z \in \{-2h, -h, 0, +h, +2h\}$
and applies fixed-coefficient stencils:

| Order | Formula | Accuracy |
|-------|---------|----------|
| $k_1$ | $(f_{-2} - 8f_{-1} + 8f_{+1} - f_{+2}) / 12h$ | $O(h^4)$ |
| $k_2$ | $(-f_{-2} + 16f_{-1} - 30f_0 + 16f_{+1} - f_{+2}) / 12h^2$ | $O(h^4)$ |
| $k_3$ | $(-f_{-2} + 2f_{-1} - 2f_{+1} + f_{+2}) / 2h^3$ | $O(h^2)$ |

where $f(\phi_z) \equiv p_e(\phi_z) - p_e(0)$.

**Cost:** $5 \times N_t$ mesolve calls.
**Precision:** $k_1$ $O(h^4)$, $k_2$ $O(h^4)$, $k_3$ $O(h^2)$.
**No ill-conditioning:** fixed-coefficient linear combination replaces Vandermonde
polynomial fitting (see `phase_10_kernel_theory_and_fix.md` for failure analysis).
**Recommended:** $h = 0.01$ rad for k3 (needs larger h than k1 due to $O(h^2)$ stencil).
""")

code("""print("Computing 5-pt FD exp kernels (order=3, h=0.01 rad)...")
q_exp = TransmonQubit(EC=0.2*2*np.pi, EJ=15*2*np.pi, T1=10000, T2=5000, n_levels=2, flux=0.0)
est_exp = KernelEstimator(
    mode='omega', method='exp', order=3,
    stim_amplitude=0.01,  # h = 0.01 rad (larger for better k3 precision)
    virtual_z_impl='math',
)

t0 = time.time()
result_exp_y = est_exp.estimate_full(ctrl_y, q_exp)
t_exp = time.time() - t0

t_k_exp = result_exp_y.t_samples
k1_exp = np.asarray(result_exp_y.kernels[0])
k2_exp = np.asarray(result_exp_y.kernels[1])
k3_exp = np.asarray(result_exp_y.kernels[2])

G1_exp = np.trapezoid(k1_exp, t_k_exp)
G2_exp = np.trapezoid(k2_exp, t_k_exp)
G3_exp = np.trapezoid(k3_exp, t_k_exp)

print(f'5-pt FD exp completed in {t_exp:.2f}s')
print(f'  h = 0.01 rad (stim_amplitude)')
print(f'  G1 = {G1_exp:.4f}')
print(f'  G2 = {G2_exp:.2e}')
print(f'  G3 = {G3_exp:.4f}')
print(f'  |G3/G1| = {abs(G3_exp/G1_exp):.4f}')
print(f'  |G2/G1| = {abs(G2_exp/G1_exp):.2e}')
""")

# ===== Cell 7: Sim vs Exp comparison =====
code("""# Full comparison: Heisenberg sim vs 5-pt FD exp
# Interpolate exp to sim time grid (sim uses deduplicated grid, exp uses full)
k1_exp_i = np.interp(t_k, t_k_exp, k1_exp)
k2_exp_i = np.interp(t_k, t_k_exp, k2_exp)
k3_exp_i = np.interp(t_k, t_k_exp, k3_exp)

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# (a) k1: sim vs exp
ax = axes[0, 0]
ax.plot(t_k, k1_sim, 'o-', color='darkgreen', markersize=4, linewidth=1.5,
        label=f'Sim (G1={G1_sim:.3f})')
ax.plot(t_k, k1_exp_i, 's--', color='darkorange', markersize=4, linewidth=1.5,
        label=f'Exp (G1={G1_exp:.3f})')
ax.axvline(x=T_pi2, color='red', linestyle='--', alpha=0.4)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('k1 (rad^-1)')
ax.set_title('(a) k1: Heisenberg sim vs 5-pt FD exp')
ax.legend()
ax.set_xlim(0, 2*T_pi2)

# (b) k1 residual
ax = axes[0, 1]
res_k1 = k1_sim - k1_exp_i
ax.plot(t_k, res_k1, 'o-', color='crimson', markersize=4, linewidth=1)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('delta k1 (rad^-1)')
rmse_k1_se = np.sqrt(np.mean(res_k1**2))
ax.set_title(f'(b) k1 residual: sim - exp, RMSE={rmse_k1_se:.2e}')
ax.set_xlim(0, 2*T_pi2)

# (c) G-integral bar chart
ax = axes[0, 2]
labels = ['G1', 'G2', 'G3']
vals_sim = [G1_sim, G2_sim, G3_sim]
vals_exp = [G1_exp, G2_exp, G3_exp]
x = np.arange(len(labels))
w = 0.35
b1 = ax.bar(x - w/2, vals_sim, w, color=['#2196F3','#9C27B0','#F44336'],
            alpha=0.8, label='Heisenberg sim')
b2 = ax.bar(x + w/2, vals_exp, w, color=['#2196F3','#9C27B0','#F44336'],
            alpha=0.4, label='5-pt FD exp', edgecolor='black', linewidth=1, hatch='//')
ax.axhline(y=0, color='gray', alpha=0.5)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('Integral value')
ax.set_title('(c) G_n = integral of k_n(t) dt')
ax.legend()
for bar, val in zip(b1, vals_sim):
    ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
            f'{val:.3f}', ha='center', va='bottom' if val > 0 else 'top', fontsize=8)
for bar, val in zip(b2, vals_exp):
    ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
            f'{val:.3f}', ha='center', va='bottom' if val > 0 else 'top', fontsize=8)

# (d) k2: sim vs exp
ax = axes[1, 0]
ax.plot(t_k, k2_sim, 'o-', color='purple', markersize=4, linewidth=1,
        label=f'Sim (G2={G2_sim:.2e})')
ax.plot(t_k_exp, k2_exp, 's--', color='magenta', markersize=4, linewidth=1,
        label=f'Exp (G2={G2_exp:.2e})')
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('k2 (rad^-2)')
ax.set_title('(d) k2: Heisenberg sim vs 5-pt FD exp (~0)')
ax.legend()
ax.set_xlim(0, 2*T_pi2)

# (e) k3: sim vs exp
ax = axes[1, 1]
ax.plot(t_k, k3_sim, 'o-', color='darkgreen', markersize=4, linewidth=1.5,
        label=f'Sim (G3={G3_sim:.3f})')
ax.plot(t_k, k3_exp_i, 's--', color='darkorange', markersize=4, linewidth=1.5,
        label=f'Exp (G3={G3_exp:.3f})')
ax.axvline(x=T_pi2, color='blue', linestyle='--', alpha=0.4)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('k3 (rad^-3)')
ax.set_title('(e) k3: Heisenberg sim vs 5-pt FD exp')
ax.legend()
ax.set_xlim(0, 2*T_pi2)

# (f) k3 residual
ax = axes[1, 2]
res_k3 = k3_sim - k3_exp_i
ax.plot(t_k, res_k3, 'o-', color='crimson', markersize=4, linewidth=1)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('delta k3 (rad^-3)')
rmse_k3_se = np.sqrt(np.mean(res_k3**2))
ax.set_title(f'(f) k3 residual: sim - exp, RMSE={rmse_k3_se:.2e}')
ax.set_xlim(0, 2*T_pi2)

plt.suptitle('Heisenberg sim vs 5-pt FD exp: Full Comparison',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('kernel_sim_vs_exp.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved: kernel_sim_vs_exp.png')
""")

# ===== Cell 8: W(t) Bloch sphere =====
md(r"""## 4. Geometric Picture: W(t) on the Bloch Sphere

$W(t) = U^\dagger(t)\sigma_z U(t)/2$ describes how the effective perturbation
axis evolves during the control pulse. The kernel $k_n(t)$ measures how $W(t)$
couples to the final measurement direction $Q$.

**Key insight:**
- Pulse 1 (Y): $W(t)$ rotates in the **Z-X plane**
- Pulse 2 (X): $W(t)$ rotates in the **X-Y plane**
- $Q$ points in the $(1,-1,0)$ direction
- $k_2=0$ because $W$ and $Q$ are always orthogonal in the right way
""")

code("""# Compute W(t) Bloch vector from sesolve
n = 2
sigma_z = qutip.Qobj(np.diag([1.0, -1.0]))
H0 = QobjEvo(0 * qeye(2))
H_pulse = QobjEvo(ctrl_y.hamiltonian, tlist=np.asarray(ctrl_y.t_list, dtype=float), order=1)
H_full = H0 + H_pulse

raw_tlist = np.asarray(ctrl_y.t_list, dtype=float)
t_grid = np.unique(raw_tlist)

res = sesolve(H_full, qeye(n), t_grid,
              options={'max_step': float(CONFIG.awg.dt), 'store_states': True})
U_list = res.states

# Extract W Bloch vector: W = (wx*sx + wy*sy + wz*sz)/2
W_bloch = np.zeros((len(t_grid), 3))
for i, U_t in enumerate(U_list):
    Z_t = U_t.dag() * sigma_z * U_t
    W_bloch[i, 0] = float((Z_t * sigmax()).tr().real) / 2
    W_bloch[i, 1] = float((Z_t * sigmay()).tr().real) / 2
    W_bloch[i, 2] = float((Z_t * sigmaz()).tr().real) / 2

# Visualization
fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# (a) W components vs time
ax = axes[0, 0]
ax.plot(t_grid, W_bloch[:, 0], 'b-', linewidth=1.5, label='wx')
ax.plot(t_grid, W_bloch[:, 1], 'r-', linewidth=1.5, label='wy')
ax.plot(t_grid, W_bloch[:, 2], 'g-', linewidth=1.5, label='wz')
ax.axvline(x=T_pi2, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('t (ns)')
ax.set_ylabel('Bloch component')
ax.set_title('(a) W(t) Bloch vector components')
ax.legend(fontsize=8)
ax.set_xlim(0, 2*T_pi2)

# (b) wx-wz plane
ax = axes[0, 1]
sc = ax.scatter(W_bloch[:, 0], W_bloch[:, 2], c=t_grid, cmap='viridis', s=10, alpha=0.8)
ax.plot(W_bloch[:, 0], W_bloch[:, 2], 'k-', alpha=0.15, linewidth=0.5)
ax.set_xlabel('wx')
ax.set_ylabel('wz')
ax.set_title('(b) W(t) in wx-wz plane (Z-X plane, pulse 1)')
ax.set_aspect('equal')
ax.axhline(y=0, color='gray', alpha=0.3)
ax.axvline(x=0, color='gray', alpha=0.3)
plt.colorbar(sc, ax=ax, label='t (ns)')

# (c) wx-wy plane
ax = axes[0, 2]
sc = ax.scatter(W_bloch[:, 0], W_bloch[:, 1], c=t_grid, cmap='viridis', s=10, alpha=0.8)
ax.plot(W_bloch[:, 0], W_bloch[:, 1], 'k-', alpha=0.15, linewidth=0.5)
ax.set_xlabel('wx')
ax.set_ylabel('wy')
ax.set_title('(c) W(t) in wx-wy plane (X-Y plane, pulse 2)')
ax.set_aspect('equal')
ax.axhline(y=0, color='gray', alpha=0.3)
ax.axvline(x=0, color='gray', alpha=0.3)
plt.colorbar(sc, ax=ax, label='t (ns)')

# (d) k1, k2, k3 overview
ax = axes[1, 0]
ax.plot(t_k, k1_sim, 'b-o', markersize=3, linewidth=1.2, label=f'k1 (G1={G1_sim:.3f})')
ax.plot(t_k, k2_sim, 's-', color='purple', markersize=3, linewidth=1.2, label=f'k2 (G2={G2_sim:.2e})')
ax.plot(t_k, k3_sim, 'r-^', markersize=3, linewidth=1.2, label=f'k3 (G3={G3_sim:.3f})')
ax.axvline(x=T_pi2, color='gray', linestyle='--', alpha=0.5)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_ylabel('kernel value')
ax.set_title('(d) k1, k2, k3 overview (Heisenberg sim)')
ax.legend(fontsize=8)
ax.set_xlim(0, 2*T_pi2)

# (e) k1 vs k3 verification
ax = axes[1, 1]
ax.plot(t_k, k1_sim, 'b-o', markersize=3, linewidth=1.5, label='k1 (sim)')
ax.plot(t_k, -k3_sim, 'r--s', markersize=3, linewidth=1.5, label='-k3 (sim)')
ax.axvline(x=T_pi2, color='gray', linestyle='--', alpha=0.5)
ax.axhline(y=0, color='gray', alpha=0.4)
ax.set_xlabel('t (ns)')
ax.set_title('(e) Verification: k3 = -k1')
ax.legend()
ax.set_xlim(0, 2*T_pi2)

# (f) Summary formulas
ax = axes[1, 2]
ax.axis('off')
summary_text = f'''Kernel commutator formulas:

k1 = i<0|[W, Q]|0>
k2 = -<0|[W, [W, Q]]|0>
k3 = -i<0|[W, [W, [W, Q]]]|0>

2-level simplification:
ad_W^3(Q) = ad_W^1(Q) = [W, Q]
=> k3 = -k1 (exact)

G1 = {G1_sim:.4f}
G2 = {G2_sim:.2e} (~0)
G3 = {G3_sim:.4f}
|G3/G1| = {abs(G3_sim/G1_sim):.4f}'''
ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
        va='top', ha='left', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.suptitle('W(t) Bloch Sphere Evolution and Kernel Geometry',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('kernel_bloch_geometry.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved: kernel_bloch_geometry.png')
""")

# ===== Cell 9: Final Summary =====
code("""# Final summary table
print("=" * 80)
print("FINAL SUMMARY: First-Order and Third-Order Kernel Functions")
print("=" * 80)

hdr = "+---------------------+------------------+------------------+------------------+"
print(hdr)
print("|                     | Analytical (sq)  | Heisenberg sim   | 5-pt FD exp      |")
print(hdr)
print(f"| G1                  | {G1_analytic:+.4f}         | {G1_sim:+.4f}         | {G1_exp:+.4f}         |")
print(f"| G2                  | 0 (strict)       | {G2_sim:+.2e}      | {G2_exp:+.2e}      |")
print(f"| G3                  | {G3_analytic:+.4f}         | {G3_sim:+.4f}         | {G3_exp:+.4f}         |")
print(f"| |G3/G1|             | {abs(G3_analytic/G1_analytic):.4f}            | {abs(G3_sim/G1_sim):.4f}            | {abs(G3_exp/G1_exp):.4f}            |")
print(f"| Time                | Instant          | {t_sim:.2f}s           | {t_exp:.2f}s           |")
print("| # QuTiP solves      | N/A              | 1 sesolve        | 5*N_t mesolve    |")
print("| Precision           | Exact            | Machine prec.    | k1 O(h4), k3 O(h2) |")
print(hdr)

print()
print("Key conclusions:")
print("  1. k1(t) = -1/2 sin(Omega*t) (pulse 1), -1/2 cos(Omega*(t-T1)) (pulse 2)")
print("     Single-peaked sinusoidal shape, peak at pulse center.")
print("  2. k2(t) = 0 (exactly at zero detuning)")
print("     Algebraic consequence of Y-X pulse orthogonality.")
print("  3. k3(t) = -k1(t) in 2-level system")
print("     Period-2 structure of nested commutators: ad_W^3 = ad_W^1.")
print("     |G3/G1| = 1 exactly.")
print("  4. Heisenberg sim: 1 sesolve -> all orders at machine precision.")
print("     No perturbation, no fitting, no finite differences.")
print("  5. 5-pt FD exp: matches sim, no Vandermonde ill-conditioning.")
print("     Fixed-coefficient stencils. h = 0.005 rad recommended.")
print()
print("Generated figures:")
print("  - kernel_analytical_square.png")
print("  - kernel_analytic_vs_sim.png")
print("  - kernel_sim_vs_exp.png")
print("  - kernel_bloch_geometry.png")
""")

# ===== Write notebook =====
with open("kernel_analytical_vs_numerical.ipynb", "w", encoding="utf-8") as f:
    nbf.write(nb, f)

print("Notebook written: kernel_analytical_vs_numerical.ipynb")
