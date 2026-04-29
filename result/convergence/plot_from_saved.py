#!/usr/bin/env python3
"""
从已保存的 lm_parameters.pkl 中读取数据，绘制完整的LM收敛分析图。
无需重新运行仿真。

用法:
    python result/convergence/plot_from_saved.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pickle

# 设置科研制图规范
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'mathtext.fontset': 'stix',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

# ---------- 加载数据 ----------
script_dir = os.path.dirname(os.path.abspath(__file__))
pkl_path = os.path.join(script_dir, 'lm_parameters.pkl')

with open(pkl_path, 'rb') as f:
    results = pickle.load(f)

history = results['history']
residuals = history['res']          # list of ndarray, 每次迭代的残差向量
mu_history = history.get('mu', [])  # list of float, 阻尼参数
b_history = history.get('b', [])    # list of ndarray, 基函数系数

n_iter = len(residuals)
iter_indices = np.arange(1, n_iter + 1)
res_norms = [np.linalg.norm(r) for r in residuals]

print(f"已加载数据: {n_iter} 次迭代, "
      f"{len(mu_history)} 个mu值, "
      f"{len(b_history)} 组系数 (n_basis={len(b_history[0]) if b_history else '?'})")

# ============================================================
# 图1: 2x2 完整收敛分析面板
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# --- (a) 残差范数 ||r|| vs 迭代次数 ---
ax = axes[0, 0]
ax.plot(iter_indices, res_norms, 'bo-', linewidth=2, markersize=5,
        markerfacecolor='white', markeredgewidth=1.5)
ax.set_xlabel('Iteration Number')
ax.set_ylabel(r'Residual Norm $\|r\|$')
ax.set_title('(a) Residual Convergence')
ax.set_yscale('log')
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(0, n_iter + 1)
# 标注最终值
ax.annotate(f'Final: {res_norms[-1]:.2e}',
            xy=(iter_indices[-1], res_norms[-1]),
            xytext=(iter_indices[-1] * 0.6, res_norms[-1] * 5),
            arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
            fontsize=11, color='red')

# --- (b) 阻尼参数 mu ---
ax = axes[0, 1]
if len(mu_history) > 0:
    mu_indices = np.arange(1, len(mu_history) + 1)
    ax.plot(mu_indices, mu_history, 'rs-', linewidth=2, markersize=5,
            markerfacecolor='white', markeredgewidth=1.5)
    ax.set_yscale('log')
ax.set_xlabel('Iteration Number')
ax.set_ylabel(r'Damping Parameter $\mu$')
ax.set_title(r'(b) Damping Parameter $\mu$ Evolution')
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(0, n_iter + 1)

# --- (c) 步长 ||Δb|| ---
ax = axes[1, 0]
if len(b_history) > 1:
    step_norms = [np.linalg.norm(np.array(b_history[i+1]) - np.array(b_history[i]))
                  for i in range(len(b_history) - 1)]
    step_indices = np.arange(1, len(step_norms) + 1)
    ax.plot(step_indices, step_norms, 'g^-', linewidth=2, markersize=5,
            markerfacecolor='white', markeredgewidth=1.5)
    ax.set_yscale('log')
ax.set_xlabel('Iteration Number')
ax.set_ylabel(r'Step Size $\|\Delta b\|$')
ax.set_title(r'(c) Parameter Step Size $\|\Delta b\|$')
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(0, n_iter + 1)

# --- (d) 基函数系数演化（变化量最大的几个） ---
ax = axes[1, 1]
top_k = 0
if len(b_history) > 1:
    b_array = np.array(b_history)  # shape (n_iter+1, n_basis)
    n_basis = b_array.shape[1]
    final_change = np.abs(b_array[-1] - b_array[0])
    top_k = min(8, n_basis)
    top_indices = np.argsort(final_change)[-top_k:][::-1]
    cmap = plt.cm.tab10
    for rank, idx in enumerate(top_indices):
        ax.plot(np.arange(b_array.shape[0]), b_array[:, idx],
                '-', color=cmap(rank), linewidth=1.5,
                label=f'$b_{{{idx}}}$')
    ax.legend(fontsize=9, ncol=2, loc='best')
ax.set_xlabel('Iteration Number')
ax.set_ylabel('Coefficient Value')
ax.set_title(f'(d) Top-{top_k} Basis Coefficients Evolution')
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
save1 = os.path.join(script_dir, 'lm_convergence_full.png')
plt.savefig(save1, dpi=300)
print(f"完整收敛分析图已保存至: {save1}")
plt.show()

# ============================================================
# 图2: 残差向量的逐分量热图 — 展示哪些测量点先收敛
# ============================================================
if n_iter > 1:
    res_matrix = np.array(residuals)  # shape (n_iter, n_meas)
    fig2, ax2 = plt.subplots(figsize=(12, 5))
    im = ax2.imshow(np.abs(res_matrix).T, aspect='auto',
                    cmap='hot_r', origin='lower',
                    extent=[1, n_iter, 0, res_matrix.shape[1]])
    ax2.set_xlabel('Iteration Number')
    ax2.set_ylabel('Measurement Point Index')
    ax2.set_title('Absolute Residual per Measurement Point')
    plt.colorbar(im, ax=ax2, label=r'$|r_i|$')
    plt.tight_layout()
    save2 = os.path.join(script_dir, 'lm_residual_heatmap.png')
    plt.savefig(save2, dpi=300)
    print(f"残差热图已保存至: {save2}")
    plt.show()

# ============================================================
# 图3: 系数范数 ||b|| 随迭代的变化 — 检测过拟合
# ============================================================
if len(b_history) > 1:
    b_norms = [np.linalg.norm(b) for b in b_history]
    fig3, ax3 = plt.subplots(figsize=(8, 5))
    ax3.plot(np.arange(len(b_norms)), b_norms, 'mD-', linewidth=2,
             markersize=5, markerfacecolor='white', markeredgewidth=1.5)
    ax3.set_xlabel('Iteration Number')
    ax3.set_ylabel(r'$\|b\|$')
    ax3.set_title(r'Coefficient Norm $\|b\|$ Evolution')
    ax3.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    save3 = os.path.join(script_dir, 'lm_coefficient_norm.png')
    plt.savefig(save3, dpi=300)
    print(f"系数范数图已保存至: {save3}")
    plt.show()

print("\n所有图像绘制完成。")
