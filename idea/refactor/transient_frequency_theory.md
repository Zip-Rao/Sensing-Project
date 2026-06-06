# 瞬态频率标定：完整理论-代码对应

> 日期: 2026-06-06 | 代码: `sqc/calibration/frequency.py:_measure_frequency_transient`

---

## 0. 总览：从 Hamiltonian 到频率估计

瞬态测频通过**正交 Ramsey 脉冲 + 核函数灵敏度**，将 $p_e$ 差分直接转换为 detuning，不依赖 Ramsey τ 扫描。

**管线**：`制备 → τ=0 Ramsey 正交对 → p_diff → G_freq → Δω`

---

## 1. 系统设定

### 物理 Hamiltonian

Transmon qubit 在磁通偏置 $\Phi$ 下，旋转坐标系（RWA, 驱动频率 $\omega_d$）：

$$H(t) = \frac{1}{2}\underbrace{(\omega_q(\Phi) - \omega_d)}_{\textstyle \equiv\,\Delta}\,\sigma_z + H_{\text{ctrl}}(t;\omega_d)$$

其中 $H_{\text{ctrl}}$ 是微波控制脉冲（Ramsey 序列），$\Delta$ 是**待测的 detuning**。

```
代码 L220-228:
    Phi = FluxSignal(type=1, amplitude=flux)    # DC 磁通偏置
    qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)
    # → H(t) = (ω_q(Φ)-ω_d)·σ_z/2 + H_ctrl
```

### 正交 Ramsey 脉冲对

两个 τ=0 的 Ramsey 序列，仅在第二个 π/2 脉冲的相位上相差 π：

$$H_{\text{ctrl}}^X: \quad R_y(\pi/2) - R_x(\pi/2) \quad (\text{phase}_1=\pi/2, \text{phase}_2=0)$$
$$H_{\text{ctrl}}^{-X}: \quad R_y(\pi/2) - R_{-x}(\pi/2) \quad (\text{phase}_1=\pi/2, \text{phase}_2=\pi)$$

```
代码 L230-243:
    ctrl_x  = create_ramsey_pulse(phase1=π/2, phase2=0)    # Y/2 → X/2
    ctrl_mx = create_ramsey_pulse(phase1=π/2, phase2=π)    # Y/2 → -X/2
```

---

## 2. p_diff 的物理含义

### 演化末态

以基态 $|0\rangle$ 为初态，经过完整脉冲序列（含 detuning $\Delta$），末态为：

$$|\psi_X\rangle = U_{\text{total}}^X(\Delta)\,|0\rangle$$
$$|\psi_{-X}\rangle = U_{\text{total}}^{-X}(\Delta)\,|0\rangle$$

激发态概率：

$$p_X = |\langle 1|\psi_X\rangle|^2,\quad p_{-X} = |\langle 1|\psi_{-X}\rangle|^2$$

### p_diff 的一阶微扰展开

将 $p_X$ 和 $p_{-X}$ 对 $\Delta$ 做 Taylor 展开（$\Delta$ 是小量）：

$$p_X(\Delta) = p_0 + \left.\frac{dp_X}{d\Delta}\right|_0 \Delta + O(\Delta^2)$$
$$p_{-X}(\Delta) = p_0 + \left.\frac{dp_{-X}}{d\Delta}\right|_0 \Delta + O(\Delta^2)$$

由脉冲序列的对称性：$dp_{-X}/d\Delta = -dp_X/d\Delta$。定义差分信号：

$$p_{\text{diff}} \equiv \frac{p_X - p_{-X}}{2} = \left.\frac{dp_X}{d\Delta}\right|_0 \Delta + O(\Delta^3)$$

**关键**：$p_{\text{diff}}$ 对 $\Delta$ 是**奇函数**——偶次项（$\Delta^2$ 等）在差分中抵消。这是一阶模型精确到三阶（$O(\Delta^3)$ 误差）的根本原因。

```
代码 L246-271:
    res_x  = mesolve(H_x,  qubit.state, t_global)  → p_x
    res_mx = mesolve(H_mx, qubit.state, t_global)  → p_mx
    p_diff = (p_x - p_mx) / 2.0
```

---

## 3. 核函数灵敏度 $G_{\text{freq}}$

### 定义

$p_X$ 对 detuning $\Delta$ 的导数可以写成脉冲序列的一阶核函数的积分：

$$\left.\frac{dp_X}{d\Delta}\right|_0 = \int_0^T k_1^X(t)\,dt \equiv G_{\text{freq}}^X$$

其中 $k_1(t) = i\langle 0|[W(t), Q]|0\rangle$ 是一阶核函数，$W(t) = U_{\text{ctrl}}^\dagger(t)\,\sigma_z/2\,U_{\text{ctrl}}(t)$，$Q = U_{\text{ctrl}}^\dagger(T)\,|1\rangle\langle 1|\,U_{\text{ctrl}}(T)$。

同理 $dp_{-X}/d\Delta|_0 = G_{\text{freq}}^{-X}$，由对称性 $G_{\text{freq}}^{-X} = -G_{\text{freq}}^X$。

定义**差分核函数积分**：

$$G_{\text{freq}} \equiv \frac{G_{\text{freq}}^X - G_{\text{freq}}^{-X}}{2} = G_{\text{freq}}^X = \int_0^T \frac{k_1^X(t) - k_1^{-X}(t)}{2}\,dt$$

### 计算

`KernelEstimator(mode='omega', method='exp')` 通过 Virtual Z 双边差分测量 $k_1(t)$：

$$k_1(t_j) = \frac{p_e(+\phi_z) - p_e(-\phi_z)}{2\phi_z}\Bigg|_{t=t_j}$$

扫描 $t_j$ 覆盖整个脉冲序列，数值积分得到 $G_{\text{freq}}$。

```
代码 L277-296:
    estimator = KernelEstimator(mode='omega', method='exp', order=1)
    result_x  = estimator.estimate_full(ctrl_x, qubit)   → k_omega_x(t)
    result_mx = estimator.estimate_full(ctrl_mx, qubit)  → k_omega_mx(t)
    
    k_omega_diff = (k_omega_x - k_omega_mx) / 2.0       # 差分核函数
    G_freq = ∫ k_omega_diff(t) dt                        # 总灵敏度
```

---

## 4. 频率反演

### 线性关系

将 §2 和 §3 结合：

$$p_{\text{diff}} = G_{\text{freq}} \cdot \Delta$$

因此：

$$\Delta = \frac{p_{\text{diff}}}{G_{\text{freq}}}$$

$$\boxed{f_{\text{meas}} = \omega_d - \frac{p_{\text{diff}}}{G_{\text{freq}}}}$$

**无需 $\kappa = d\omega/d\Phi$ 换算**——$G_{\text{freq}}$ 直接是频率灵敏度（单位：无量纲，因为 $p_e$ 无量纲，$\Delta$ 单位 GHz）。

### 与旧方法的对比

| | 旧方法（Phase 10.5 前） | 新方法（Phase 10.5 后） |
|---|---|---|
| 核函数 | flux kernel ($k_1^{(\Phi)}$) | omega kernel ($k_1^{(\omega)}$) |
| $G$ | $G_{\text{flux}} = \int k_1^{(\Phi)} dt$ | $G_{\text{freq}} = \int k_1^{(\omega)} dt$ |
| 反演 | $\Delta = p_{\text{diff}} \cdot \kappa / G_{\text{flux}}$ | $\Delta = p_{\text{diff}} / G_{\text{freq}}$ |
| 需要 $\kappa$？ | **是**（额外一步，引入非线性误差） | **否**（直接） |

```
代码 L297-303:
    if abs(G_freq) < 1e-5: return omega_d
    delta_omega = p_diff / G_freq         # ← 直接！无 κ
    return float(omega_d - delta_omega)
```

---

## 5. 完整数据流

```
┌─────────────────────────────────────────────────────────────┐
│ Step 1: qubit 状态设置                                      │
│   Phi = FluxSignal(amplitude=flux)                          │
│   qubit.qubit_in_mag(Phi, omega_d)                          │
│                                                             │
│   物理: H(t) = Δ·σ_z/2 + H_ctrl(t)                          │
├─────────────────────────────────────────────────────────────┤
│ Step 2: 正交 Ramsey 脉冲对                                   │
│   ctrl_x  = R_y(π/2)-R_x(π/2), τ=0                         │
│   ctrl_mx = R_y(π/2)-R_{-x}(π/2), τ=0                      │
│                                                             │
│   物理: 两个仅在末态读出轴相差 π 的脉冲序列                    │
├─────────────────────────────────────────────────────────────┤
│ Step 3: 仿真 p_diff                                         │
│   p_x  = mesolve(H_base + ctrl_x,  |0>) → |⟨1|ψ_X⟩|²       │
│   p_mx = mesolve(H_base + ctrl_mx, |0>) → |⟨1|ψ_{-X}⟩|²    │
│   p_diff = (p_x - p_mx)/2                                   │
│                                                             │
│   物理: p_diff = G_freq · Δ + O(Δ³)                         │
│   偶次项在差分中严格抵消 → O(Δ³) 而非 O(Δ²)                   │
├─────────────────────────────────────────────────────────────┤
│ Step 4: 核函数估计                                           │
│   k_omega_x(t)  = KernelEstimator(omega, exp, order=1)      │
│   k_omega_mx(t) = KernelEstimator(omega, exp, order=1)      │
│   k_diff(t) = (k_x - k_mx)/2                                │
│                                                             │
│   物理: k₁(t) = i⟨0|[W(t), Q]|0⟩                            │
│   G_freq = ∫ k_diff(t) dt = dp_X/dΔ|₀                        │
├─────────────────────────────────────────────────────────────┤
│ Step 5: 频率反演                                             │
│   Δ = p_diff / G_freq                                       │
│   f_meas = ω_d - Δ                                          │
│                                                             │
│   物理: 一阶 Taylor 逆: Δ = (p - p₀) / (dp/dΔ|₀)             │
└─────────────────────────────────────────────────────────────┘
```

---

## 6. 有效范围与误差分析

### 线性区边界

$p_{\text{diff}}$ 的严格表达式（考虑三阶项）：

$$p_{\text{diff}}(\Delta) = G_{\text{freq}}\,\Delta + \frac{1}{6}G_3\,\Delta^3 + O(\Delta^5)$$

其中 $G_3 = \int k_3^{\text{diag}}(t)\,dt$（由 Heisenberg 传播子验证：$|G_3/G_1| = 1$ 在零 detuning）。

线性近似 $p_{\text{diff}} \approx G_{\text{freq}}\Delta$ 的**相对误差**：

$$\epsilon_{\text{rel}} \approx \frac{|G_3|}{6|G_{\text{freq}}|} \cdot \Delta^2 \approx \frac{1}{6}\Delta^2$$

对于 $\epsilon_{\text{rel}} < 5\%$：$|\Delta| < \sqrt{0.3} \approx 0.055$ GHz。这正是实验中观察到的交叉点（~0.07 GHz）。

### 三阶修正（可选）

当 $|\Delta|$ 超出线性区时，用 Newton 迭代解三阶方程：

$$\Delta^{(0)} = p_{\text{diff}} / G_{\text{freq}}$$
$$\Delta^{(k+1)} = \Delta^{(k)} - \frac{G_{\text{freq}}\Delta + G_3\Delta^3/6 - p_{\text{diff}}}{G_{\text{freq}} + G_3\Delta^2/2}$$

这在 $|\Delta| \approx 0.12$ GHz 处提供约 7.6% 改善，但无法根本解决大 detuning 下的非线性饱和。

---

## 7. 与标准 Ramsey τ-sweep 的关系

标准 Ramsey（$\tau$ 变化）测量的是**相位积累**：

$$\phi(\tau) = \Delta \cdot \tau \quad \Rightarrow \quad f_{\text{osc}} = \Delta / 2\pi$$

FTT 提取振荡频率 → $\Delta$。精度 $\sim 1/\tau_{\max} \approx 5$ MHz（$\tau_{\max}=200$ ns）。

瞬态方法（$\tau=0$）测量的是**瞬时灵敏度**：

$$p_{\text{diff}} = \underbrace{\left(\int k_1(t)\,dt\right)}_{G_{\text{freq}}} \cdot \Delta$$

$G_{\text{freq}}$ 由脉冲序列的核函数决定，通常 $G_{\text{freq}} \sim 6$（无量纲）。

**对比**：
- Ramsey: $p_e$ 的**频率** $\propto \Delta$ → 需要扫描 $\tau$ → 慢，但全范围有效
- Transient: $p_e$ 的**幅度** $\propto \Delta$ → 无需扫描 → 快，但限于小 $\Delta$

两者在 $|\Delta| \approx 0.07$ GHz 处形成互补交叉。
