"""
绘制 Transmon qubit 频率-磁通响应函数，展示不同信号幅度下的非线性效应。
用于解释 Wiener 反卷积在大幅度信号下失效的物理原因。
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

# ── Transmon 参数（与仿真一致） ──
EC = 0.2 * 2 * np.pi   # GHz
EJ0 = 10.0 * 2 * np.pi  # GHz

# 频率-磁通关系
def omega_01(Phi):
    EJ = EJ0 * np.abs(np.cos(np.pi * Phi))
    return np.sqrt(8 * EJ * EC) - EC

# 工作点
Phi_w = np.arctan(np.sqrt(2)) / np.pi  # ≈ 0.3927 Phi_0

# ── 图1：全局响应曲线 + 工作点标注 ──
Phi_full = np.linspace(0, 0.5, 1000)
omega_full = omega_01(Phi_full)

# 在工作点处做泰勒展开
omega_w = omega_01(Phi_w)
dPhi = 1e-6
kappa1 = (omega_01(Phi_w + dPhi) - omega_01(Phi_w - dPhi)) / (2 * dPhi)  # 一阶导
kappa2 = (omega_01(Phi_w + dPhi) - 2 * omega_w + omega_01(Phi_w - dPhi)) / (dPhi**2)  # 二阶导

# 线性近似和二阶近似
delta_Phi = Phi_full - Phi_w
omega_linear = omega_w + kappa1 * delta_Phi
omega_quadratic = omega_w + kappa1 * delta_Phi + 0.5 * kappa2 * delta_Phi**2

# 仿真中测试的幅度列表
amplitudes = [0.01, 0.02, 0.04, 0.06, 0.08, 0.10]
colors_amp = plt.cm.Reds(np.linspace(0.3, 0.9, len(amplitudes)))

fig, axes = plt.subplots(1, 3, figsize=(16, 5), gridspec_kw={'width_ratios': [1.2, 1, 1]})

# ── Panel (a): omega(Phi) 全局 + 线性近似 ──
ax = axes[0]
ax.plot(Phi_full, omega_full / (2 * np.pi), 'k-', lw=2, label=r'$\omega_{01}(\Phi)$')
ax.plot(Phi_full, omega_linear / (2 * np.pi), 'b--', lw=1.5, label='Linear approx.')
ax.plot(Phi_full, omega_quadratic / (2 * np.pi), 'g-.', lw=1.5, label='Quadratic approx.')

# 标注工作点
ax.axvline(Phi_w, color='gray', ls=':', lw=1)
ax.plot(Phi_w, omega_w / (2 * np.pi), 'ro', ms=8, zorder=5)
ax.annotate(r'$\Phi_w$', xy=(Phi_w, omega_w / (2 * np.pi)),
            xytext=(Phi_w + 0.03, omega_w / (2 * np.pi) + 0.3),
            fontsize=12, color='r',
            arrowprops=dict(arrowstyle='->', color='r', lw=1.5))

# 用色带标注各幅度范围
for i, A in enumerate(amplitudes):
    ax.axvspan(Phi_w - A, Phi_w + A, alpha=0.08, color=colors_amp[i])

ax.set_xlabel(r'Flux $\Phi / \Phi_0$', fontsize=13)
ax.set_ylabel(r'Frequency $\omega_{01}/2\pi$ (GHz)', fontsize=13)
ax.set_title('(a) Transmon frequency response', fontsize=13)
ax.legend(fontsize=10, loc='upper right')
ax.set_xlim(0.15, 0.65)
ax.set_ylim(omega_01(0.6) / (2 * np.pi) - 0.2, omega_01(0.2) / (2 * np.pi) + 0.2)

# ── Panel (b): 工作点附近放大，标注各幅度区间 ──
ax = axes[1]
Phi_zoom = np.linspace(Phi_w - 0.12, Phi_w + 0.12, 500)
omega_zoom = omega_01(Phi_zoom)
delta_zoom = Phi_zoom - Phi_w
omega_lin_zoom = omega_w + kappa1 * delta_zoom
omega_quad_zoom = omega_w + kappa1 * delta_zoom + 0.5 * kappa2 * delta_zoom**2

ax.plot(Phi_zoom, omega_zoom / (2 * np.pi), 'k-', lw=2, label=r'$\omega_{01}(\Phi)$')
ax.plot(Phi_zoom, omega_lin_zoom / (2 * np.pi), 'b--', lw=1.5, label='Linear')
ax.plot(Phi_zoom, omega_quad_zoom / (2 * np.pi), 'g-.', lw=1.5, label='Quadratic')

# 标注每个幅度
for i, A in enumerate(amplitudes):
    ax.axvspan(Phi_w - A, Phi_w + A, alpha=0.12, color=colors_amp[i])
    # 标注幅度值
    y_pos = omega_01(Phi_w + A) / (2 * np.pi)
    ax.annotate(f'{A}', xy=(Phi_w + A, y_pos),
                xytext=(Phi_w + A + 0.005, y_pos),
                fontsize=9, color=colors_amp[i], fontweight='bold',
                va='center')

ax.plot(Phi_w, omega_w / (2 * np.pi), 'ro', ms=8, zorder=5)
ax.set_xlabel(r'Flux $\Phi / \Phi_0$', fontsize=13)
ax.set_ylabel(r'Frequency $\omega_{01}/2\pi$ (GHz)', fontsize=13)
ax.set_title('(b) Zoom at working point', fontsize=13)
ax.legend(fontsize=10)

# ── Panel (c): 相对误差 δω_linear / δω_exact - 1 vs 幅度 ──
ax = axes[2]
A_scan = np.linspace(0.001, 0.12, 200)
# 对每个幅度，计算在整个 [-A, +A] 区间内线性近似的最大相对误差
max_rel_error = np.zeros_like(A_scan)
rms_rel_error = np.zeros_like(A_scan)
for i, A in enumerate(A_scan):
    Phi_test = np.linspace(Phi_w - A, Phi_w + A, 200)
    delta_omega_exact = omega_01(Phi_test) - omega_w
    delta_omega_linear = kappa1 * (Phi_test - Phi_w)
    # 避免除零
    mask = np.abs(delta_omega_exact) > 1e-10
    if mask.any():
        rel_err = np.abs((delta_omega_linear[mask] - delta_omega_exact[mask]) / delta_omega_exact[mask])
        max_rel_error[i] = np.max(rel_err)
        rms_rel_error[i] = np.sqrt(np.mean(rel_err**2))

ax.plot(A_scan, max_rel_error * 100, 'k-', lw=2, label='Max relative error')
ax.plot(A_scan, rms_rel_error * 100, 'b--', lw=1.5, label='RMS relative error')

# 标注仿真中的幅度点
rmse_wiener = [0.000154, 0.000172, 0.000439, 0.000881, 0.001479, 0.003017, 0.004856]
rmse_lm =     [0.000097, 0.000137, 0.000193, 0.000289, 0.000386, 0.000510, 0.001104]
amp_data =    [0.01,     0.02,     0.04,     0.05,     0.06,     0.08,     0.10]

for i, A in enumerate(amplitudes):
    idx = np.argmin(np.abs(A_scan - A))
    ax.plot(A, max_rel_error[idx] * 100, 'o', color=colors_amp[i], ms=8, zorder=5)
    ax.annotate(f'{A}', xy=(A, max_rel_error[idx] * 100),
                xytext=(A + 0.003, max_rel_error[idx] * 100 + 1),
                fontsize=9, color=colors_amp[i], fontweight='bold')

# 标注 ~5% 误差线作为参考
ax.axhline(5, color='gray', ls=':', lw=1)
ax.text(0.005, 5.5, '5% threshold', fontsize=9, color='gray')

ax.set_xlabel(r'Signal amplitude ($\Phi_0$)', fontsize=13)
ax.set_ylabel('Linearization error (%)', fontsize=13)
ax.set_title('(c) Linear approx. error vs amplitude', fontsize=13)
ax.legend(fontsize=10)
ax.set_xlim(0, 0.12)

plt.tight_layout()
plt.savefig('result/signal_amp/transmon_nonlinearity.png', dpi=200, bbox_inches='tight')
plt.savefig('result/signal_amp/transmon_nonlinearity.pdf', bbox_inches='tight')
print("Saved to result/signal_amp/transmon_nonlinearity.png")
plt.show()
