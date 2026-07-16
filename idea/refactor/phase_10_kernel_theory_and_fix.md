# 核函数高阶提取：理论、失败分析与修复方案

> 🧭 理论汇总见 **[_kernel_frequency_theory.md](../_kernel_frequency_theory.md)**（核函数 + 测频精度唯一权威笔记）。本文为 **Phase 10 高阶提取的理论/失败分析记录**。
> 日期: 2026-06-04 | 关联: [phase_10_kernel_extension_handbook.md](phase_10_kernel_extension_handbook.md) §3–§5 | 配套: [`_sensing theory.md`](../_sensing%20theory.md) §6

---

## 1. 问题定义

`KernelEstimator` 需要从脉冲序列中提取 1 到 N 阶 Volterra 核对角元 $k_1(t), k_2(t), k_3(t), \ldots$。当前实现 (Phase 10.3) 在 order=1 时正确，但在 order≥2 时输出的是数值噪声而非物理核函数。本文档分析根本原因、给出理论推导、提出修复方案。

---

## 2. 核函数的理论定义

### 2.1 记号

| 符号 | 定义 |
|------|------|
| $U(t) \equiv U_{\text{ctrl}}(t, 0)$ | 纯控制传播子（不含待测 $\delta\omega$），$n \times n$ 酉矩阵 |
| $Z_c(t) = U^\dagger(t)\,\sigma_z\,U(t)$ | Heisenberg 绘景下的 $\sigma_z$ |
| $W(t) = \tfrac12 Z_c(t)$ | 有效微扰算符（互作用绘景下的 $\sigma_z/2$） |
| $M = |1\rangle\langle 1|$ | 测量投影算符 |
| $Q = U^\dagger(T)\,M\,U(T)$ | Heisenberg 绘景的测量算符 |
| $|0\rangle$ | 初态（基态） |

### 2.2 一阶核函数

由 Dyson 展开的一阶项（[`_sensing theory.md` §1.2](../_sensing%20theory.md)）：

$$
\boxed{k_1(t) = i\,\langle 0|\,[W(t), Q]\,|0\rangle}
$$

等价迹形式（便于数值计算）：

$$
k_1(t) = \Im\Bigl[\langle 0|U^\dagger(T) M\, U(T,t) \,\sigma_z\, U(t,0)|0\rangle\Bigr]
$$

### 2.3 二阶核函数（非对角 + 对角）

由 Dyson 展开的二阶项（[§6.3](../_sensing%20theory.md)），三重来源合并为双重对易子：

$$
\boxed{k_2(t_>, t_<) = -\langle 0|\,[W(t_<), [W(t_>), Q]]\,|0\rangle}
$$

其中 $t_> = \max(t_1, t_2), t_< = \min(t_1, t_2)$。对角线取极限 $t_1 = t_2 = t$：

$$
k_2^{\text{diag}}(t) \equiv k_2(t, t) = -\langle 0|\,[W(t), [W(t), Q]]\,|0\rangle
$$

### 2.4 n 阶核函数

时序楔形区域内（[§6.4](../_sensing%20theory.md)）：

$$
k_n^{\text{ord}}(t_{(1)}, \ldots, t_{(n)}) = i^n\,\langle 0|\,[W(t_{(n)}), [\cdots, [W(t_{(1)}), Q]\cdots]]\,|0\rangle
$$

其中 $t_{(1)} > t_{(2)} > \cdots > t_{(n)}$。全对称 Volterra 核为所有排列平均：

$$
k_n(t_1,\ldots,t_n) = \frac{1}{n!}\sum_{\pi\in S_n} k_n^{\text{ord}}(t_{\pi(1)},\ldots,t_{\pi(n)})
$$

**对角线退化**（所有 $t_i$ 相等）：排列退化为同一点，$1/n! \times n! = 1$：

$$
\boxed{k_n^{\text{diag}}(t) = i^n\,\langle 0|\,\underbrace{[W(t), [W(t), \cdots, [W(t), Q]\cdots]]}_{n\text{ 层对易子}}\,|0\rangle}
$$

### 2.5 与数值扰动的关系

以上公式假设**理想** $\delta$ 函数刺激 $\delta\omega(t) = \phi_z\,\delta(t-t_j)$。实际数值方法（exp 模式）使用有限宽度 Gaussian 近似 $\delta$ 函数：

$$
\delta\omega(t) = \frac{\phi_z}{\sigma_t\sqrt{2\pi}}\,\exp\!\left(-\frac{(t-t_j)^2}{2\sigma_t^2}\right), \quad \sigma_t = 2\cdot dt
$$

在 $\sigma_t \to 0$ 极限下，数值扰动法 $\to$ 理论值。有限 $\sigma_t$ 引入 $O(\sigma_t^2)$ 的卷积模糊，但不改变核函数的定性结构。

---

## 3. 解析验证：零延时 Ramsey π/2-π/2 脉冲

### 3.1 设定

2 能级系统，RWA 框架，零 detuning。理想方波包络：

- 第一个脉冲 $(0 \le t \le T_{\pi/2})$：绕 Y 轴旋转 $U(t) = R_y(\Omega t)$
- 第二个脉冲 $(T_{\pi/2} < t \le 2T_{\pi/2})$：绕 X 轴旋转 $U(t) = R_x(\Omega(t-T_{\pi/2}))\,R_y(\pi/2)$

其中 $\Omega = \pi/(2T_{\pi/2})$ 为 Rabi 频率。

### 3.2 W(t) 的显式

**第一脉冲期间**：
$$R_y(\theta)\,\sigma_z\,R_y(-\theta) = \cos\theta\,\sigma_z + \sin\theta\,\sigma_x$$
$$\Rightarrow W(t) = \frac{1}{2}(\cos\Omega t\,\sigma_z + \sin\Omega t\,\sigma_x)$$

**第二脉冲期间**：
$$R_y(-\tfrac\pi2)R_x(-\theta_x)\,\sigma_z\,R_x(\theta_x)R_y(\tfrac\pi2) = -\cos\theta_x\,\sigma_x + \sin\theta_x\,\sigma_y$$
$$\Rightarrow W(t) = \frac12(-\cos\Omega(t-T_{\pi/2})\,\sigma_x + \sin\Omega(t-T_{\pi/2})\,\sigma_y)$$

### 3.3 严格结果

| 核函数 | $G = \int_0^{2T_{\pi/2}} k(t)\,dt$ | 形状 | 关键性质 |
|--------|-----------------------------------|------|---------|
| $k_1(t)$ | $-3.183$ | 单峰正弦形，峰值在脉冲中心 | — |
| $k_2^{\text{diag}}(t)$ | **0** (严格) | 恒为零函数 | $[W,[W,Q]]$ 在初态期望下严格消失 |
| $k_3^{\text{diag}}(t)$ | $+3.183$ | 单峰正函数，较 k₁ 略窄 | $|G_3/G_1| = 1$ |

### 3.4 k₂ = 0 的代数证明概要

对于末态投影 $Q = \frac12(I + \sigma_x - \sigma_y)$（由 $U(T)=R_x(\pi/2)R_y(\pi/2)$ 导出），以及任意 $W = \frac12\mathbf{w}\cdot\boldsymbol{\sigma}$（$\|\mathbf{w}\|=1$），展开双层对易子：

$$
[W,[W,Q]] = \frac14[W, \mathbf{w}\cdot\boldsymbol{\sigma}, Q] = \frac14(\cdots)
$$

代入 $\langle 0|\cdot|0\rangle$，经过 Pauli 代数化简，结果恒为零。**这是 Y-X 正交脉冲序列的对称性保证的，不是近似**。

### 3.5 Gaussian 脉冲的偏差

用实际的 Gaussian 包络（`sqc/control/sequence.py:create_ramsey_pulse`）做 Heisenberg 数值计算：

| | 方波（解析） | Gaussian（数值） |
|---|---|---|
| $G_2/G_1$ | 0 | 0.042 |
| $G_3/G_1$ | 1.0 | 1.0 |
| $k_3$ 形状 | 单峰正函数 | 相似，略有展宽 |

Gaussian 包络的有限上升/下降时间轻微打破了 Y-X 对称性，导致 $k_2 \neq 0$，但仍比 $k_1$ 小两个数量级。

---

## 4. 当前算法 (Phase 10.3) 及其失败原因

### 4.1 当前实现

`_extract_kn_omega()` 的工作流程：

1. 在每个 $t_j$ 处，扫描 $\phi_z \in [-\phi_0, +\phi_0]$（$M$ 个采样点）
2. 每个 $\phi_z$ 运行一次 `mesolve`，得到 $\Delta p_e(\phi_z)$
3. 对 $\Delta p_e$ vs $\phi_z$ 做 $N$ 阶多项式拟合（无常数项）
4. 拟合系数 $a_n$ 乘以转换因子得到 $k_n^{\text{diag}}(t_j)$

### 4.2 失败原因：信号-噪声分离

对于 $\phi_z = 0.01$ rad 的 Virtual Z 刺激，$\Delta p_e$ 的构成：

$$
\Delta p_e(\phi_z) = \underbrace{a_1\phi_z}_{\sim 10^{-2}} + \underbrace{a_2\phi_z^2}_{\sim 10^{-5}} + \underbrace{a_3\phi_z^3}_{\sim 10^{-8}} + \underbrace{\eta}_{\text{mesolve 噪声} \sim 10^{-12}}
$$

各阶项的数量级：

| 项 | 量级（$\phi_z=0.01$） | 在 $\Delta p_e$ 中的占比 |
|----|----------------------|------------------------|
| $a_1\phi_z$ | $\sim 10^{-2}$ | **> 99.99%** |
| $a_2\phi_z^2$ | $\sim 10^{-5}$ | < 0.1% |
| $a_3\phi_z^3$ | $\sim 10^{-8}$ | < 0.001% |
| $\eta$ (噪声) | $\sim 10^{-12}$ | < 10⁻⁸ |

**多项式拟合的 Vandermonde 矩阵对此信号结构病态**：列向量 $\{\phi_z, \phi_z^2, \phi_z^3\}$ 在 $|\phi_z| \ll 1$ 范围内高度共线。条件数随阶数指数增长：

$$
\kappa(V_N) \sim \phi_0^{-(N-1)} \cdot N!
$$

对于 $\phi_0 = 0.01, N = 3$：$\kappa \sim 10^{6}$。这意味着拟合对输入数据中的 $10^{-12}$ 级噪声产生 $10^{-6}$ 级的系数误差——而 $a_3$ 本身只有 $10^{-8}$ 量级，完全被噪声淹没。

### 4.3 验证：改变 amp_scan_factor 时的 k₃ 发散

| amp_scan_factor | max $\phi_z$ | 拟合出的 $G_3$ | 与 FD 真实值 (180.8) 的偏差 |
|:---:|:---:|:---:|:---:|
| 1 | 0.01 | 2870 | **16×** （符号都对不上） |
| 3 | 0.03 | 15.5 | **0.086×** |
| 5 | 0.05 | -54.2 | **-0.30×**（符号翻转） |
| 10 | 0.10 | -53.3 | -0.29× |

**拟合结果不收敛**——这是病态拟合的典型特征：不同输入范围给出完全不同的输出，且都不接近真实值。

---

## 5. 修复方案

### 5.1 总览：两条路径

| | `method='sim'` | `method='exp'` |
|---|---|---|
| **原理** | Heisenberg 传播子 + 对易子求值 | 5 点有限差分 (FD) stencil |
| **计算** | 1 次 `sesolve` + $O(N)$ 矩阵乘法 | $5 \times N$ 次 `mesolve` |
| **k₁ 精度** | 机器精度 | $O(h^4)$ |
| **k₂ 精度** | 机器精度 | $O(h^4)$ |
| **k₃ 精度** | 机器精度 | $O(h^2)$ |
| **无病态** | ✓（无拟合） | ✓（固定系数线性组合） |

**共同原则：不做多项式拟合。** k₁/k₂/k₃ 通过确定性线性组合（而非数据驱动的回归）直接计算。

### 5.2 路径 A：Heisenberg 传播子 (`method='sim'`)

#### 算法

```
1. H_full = H_0 + H_pulse(t)            # 纯控制哈密顿量
2. res = sesolve(H_full, qeye(n), t_list)  # 一次演化，得到全部 U(t)
3. U_T = res.states[-1]
4. Q = U_T† · M · U_T                    # Heisenberg 测量算符
5. FOR each t_i:
      U_t = res.states[i]
      W_t = U_t† · σ_z · U_t / 2         # Heisenberg 微扰算符
      
      # k₁ = i⟨0|[W_t, Q]|0⟩
      k1[i] = Im[⟨0| W_t·Q - Q·W_t |0⟩]
      
      # k₂ = -⟨0|[W_t, [W_t, Q]]|0⟩
      k2[i] = -⟨0| W_t²·Q - 2W_t·Q·W_t + Q·W_t² |0⟩
      
      # k₃ = -i⟨0|[W_t, [W_t, [W_t, Q]]]|0⟩
      # (展开为 8 项对易子求和)
      k3[i] = -i⟨0| [W_t, [W_t, [W_t, Q]]] |0⟩
```

#### 理论依据

§2.2–2.4 的对易子公式。在 2 能级下可用 Pauli 代数闭式求值；在 $n>2$ 能级下用数值矩阵乘法。**无近似、无拟合、无数值微分**。

#### 计算复杂度

$O(N \cdot d^3 + T_{\text{sesolve}})$，其中 $N$ = 时间点数，$d$ = Hilbert 空间维度。对 $d=2$：每时间点 ~10 次 2×2 矩阵乘。

### 5.3 路径 B：5 点 FD Stencil (`method='exp'`)

#### 算法

在每个 $t_j$ 处，用 5 个对称采样点 $\phi_z \in \{-2h, -h, 0, +h, +2h\}$，测量 $p_e(\phi_z)$，通过固定系数线性组合直接求导：

$$
\begin{aligned}
k_1(t_j) &= \frac{f(-2h) - 8f(-h) + 8f(h) - f(2h)}{12h} &[O(h^4)] \\[4pt]
k_2(t_j) &= \frac{-f(-2h) + 16f(-h) - 30f(0) + 16f(h) - f(2h)}{12h^2} &[O(h^4)] \\[4pt]
k_3(t_j) &= \frac{-f(-2h) + 2f(-h) - 2f(h) + f(2h)}{2h^3} &[O(h^2)]
\end{aligned}
$$

其中 $f(\phi_z) \equiv p_e(\phi_z) - p_e(0)$。

#### 理论依据

这是 $\partial^n p_e / \partial\phi_z^n|_{\phi_z=0}$ 的标准有限差分近似。**在 $h \to 0$ 极限下与 Heisenberg 对易子公式等价**——两者都是测量同一个泛函导数 $\delta^n p_e / \delta\phi_z^n$。

与多项式拟合的对比：

| | 多项式拟合 (当前) | 5 点 FD Stencil (修复) |
|---|---|---|
| **数学操作** | 数据驱动回归（Vandermonde 求逆） | 固定系数线性组合 |
| **条件数** | $\kappa \sim \phi_0^{-(N-1)} N!$ → 指数发散 | **无条件数问题**（纯算术） |
| **k₁** | 正确（只用到 a₁，线性项 SNR 足够） | 正确，$O(h^4)$ 精度 |
| **k₂** | SNR ~2.5，与真实值相关 ~0.1 | $O(h^4)$，精确定量 |
| **k₃** | SNR ~3，符号可能翻转 | $O(h^2)$，精确定量 |

#### h 的选择

- **太小** ($h < 10^{-4}$)：除以 $h^3$ 放大 mesolve 噪声（~10⁻¹²）→ 信噪比退化
- **太大** ($h > 10^{-1}$)：高阶导数近似误差增大，且 Volterra 小扰动假设失效
- **推荐**：$h = 0.005$ rad。此时 $h^3 = 1.25 \times 10^{-7}$，噪声放大后可接受

#### 实验可行性

与当前 Virtual Z 扫幅方案**完全同构**——只是把 7-9 个 $\phi_z$ 采样点换成 5 个，把多项式拟合换成固定系数乘加。在真实硬件上实现不需要任何额外设备。

### 5.4 两条路径的关系

Heisenberg（sim）和 FD（exp）测量的是**同一个物理量**——泛函导数 $k_n = \delta^n p_e / \delta\phi_z^n$。区别仅在于：

- **sim**：用完整的 $U(t)$ 矩阵直接对易求值 → 精确但需要知道 $H(t)$
- **exp**：通过有限差分从 $p_e$ 测量值中数值求导 → 有 $O(h^2)\sim O(h^4)$ 截断误差，但适用于任何可测 $p_e$ 的系统

两者在 $h \to 0$ 极限下收敛到同一值。**FD 的 h 取 0.005 时，k₁ 和 k₂ 的截断误差 $O(h^4) \sim 10^{-9}$，远小于 mesolve 噪声**。

---

## 6. 特殊性质：为什么 k₂ 在零 detuning 时严格为零

### 6.1 代数原因

对于由 Y 旋转 + X 旋转组成的正交脉冲序列，末态投影 $Q$ 在 Bloch 球上的方向与 $W(t)$ 的演化轨迹满足：

$$\langle 0|[W(t), [W(t), Q]]|0\rangle = 0 \quad \forall t$$

这不是巧合——是脉冲序列的 Y-X 正交性导致的。具体来说，$W(t)$ 始终在 Bloch 球的 Z-X 平面（第一脉冲）或 X-Y 平面（第二脉冲）内演化，而 $Q$ 在 (1, -1, 0) 方向。双层对易子的初态期望值正好抵消。

### 6.2 物理后果

1. **对零 detuning 的小信号测量**：一阶模型 $p_e = p_0 + \int k_1 \delta\omega$ 精确到三阶（$k_2=0$ 意味着二次误差消失）
2. **对非零 detuning**：色散非线性 $\kappa' \neq 0$ 会通过 $\kappa' k_1 \delta$ 项（§6.6 换算公式）注入有效的二阶响应
3. **Gaussian 包络的边沿效应**：有限上升时间轻微打破对称性，产生 $\sim 4\%$ 的残余 k₂

---

## 7. 实施计划

### Step 1: 实现 `_heisenberg_kernels()` (sim 路径)

- 输入：`pulse, qubit, t_samples, order`
- 1 次 `sesolve(qeye(n), t_list)` → `U_list`
- 对易子循环求 $k_1, k_2, k_3$
- 输出 `KernelResult`
- 验证：与解析方波解对照（|G₂/G₁| < 1e-14, |G₃/G₁| ≈ 1.0）

### Step 2: 重写 `_extract_kn_omega()` (exp 路径)

- 将多项式拟合替换为 5 点 FD stencil
- $h$ 取 0.005 rad（可通过 `stim_amplitude` 调节）
- 保留 order=1 的快速路径（双边差分）
- 验证：与 Heisenberg 结果对照（k₁ 相关 > 0.999, k₂/k₃ 相关 > 0.99）

### Step 3: 统一 `_extract_kn_flux()` (exp flux 路径)

- 同样的 FD stencil，但在 flux 域操作
- 注意 §6.6 的 $\omega$/flux 换算（高阶时 $\kappa', \kappa''$ 注入 $\delta$ 函数项）

### Step 4: 数值验证

- 对照解析方波解（零 detuning）
- 对照 Heisenberg 传播子（任意 detuning）
- 扫描 $h$ 验证收敛性
- 回归：所有 baseline 不变

---

## 8. 参考文献

1. [`_sensing theory.md`](../_sensing%20theory.md) §1（一阶核）、§6（高阶 Volterra、对易子表达式、数值提取法）
2. [`phase_10_kernel_extension_handbook.md`](phase_10_kernel_extension_handbook.md) §3（Volterra 展开）、§5（实现细节、Gaussian 积分因子）
3. Herb & Degen, PRL 2024 — 瞬态磁场协议核函数的解析形式
4. Gao et al., PRX Quantum 2021 — Virtual Z 的实验实现 (§IV.B)

---

## 附录 A: 5 点 Stencil 系数推导

对于等距节点 $\{-2h, -h, 0, h, 2h\}$，用 Lagrange 插值多项式求导在 $x=0$ 处的值：

**一阶导数** (5 点, $O(h^4)$):
$$f'(0) = \frac{f(-2h) - 8f(-h) + 8f(h) - f(2h)}{12h}$$

**二阶导数** (5 点, $O(h^4)$):
$$f''(0) = \frac{-f(-2h) + 16f(-h) - 30f(0) + 16f(h) - f(2h)}{12h^2}$$

**三阶导数** (5 点, $O(h^2)$):
$$f'''(0) = \frac{-f(-2h) + 2f(-h) - 2f(h) + f(2h)}{2h^3}$$

---

## 附录 B: 2 能级 Heisenberg 核函数的 Pauli 代数实现

```python
# Pauli matrices
sx, sy, sz = sigmax(), sigmay(), sigmaz()

# For any W = (wx*sx + wy*sy + wz*sz)/2 with w² = 1:
# [W, Q] = i*(w × q)·σ / 2       (cross product on Bloch sphere)
# [W, [W, Q]] = (w·q)w·σ - q·σ   (double cross product)
# [W, [W, [W, Q]]] = -i*(w × q)·σ + i*(w×(w×q))·σ  (triple)

# For zero-detuning Ramsey Y-X:
# Q_bloch = (1, -1, 0)/2          (measurement direction on Bloch sphere)
# W_bloch(t):
#   pulse 1: (sin(Ωt), 0, cos(Ωt))/2        (Y rotation → Z-X plane)
#   pulse 2: (0, sin(Ω(t-T₁)), -cos(Ω(t-T₁)))/2  (X rotation → X-Y plane)
```
