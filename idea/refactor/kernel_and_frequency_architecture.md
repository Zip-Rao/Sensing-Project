# 核函数体系与频率标定体系 — 架构总结

> 日期: 2026-06-06 | 关联: [phase_10_kernel_extension_handbook.md](phase_10_kernel_extension_handbook.md), [phase_10_kernel_theory_and_fix.md](phase_10_kernel_theory_and_fix.md), [`_sensing theory.md`](../_sensing%20theory.md)

---

## 一、核函数体系 (`sqc/reconstruction/kernel.py`)

### 1.1 核心类

| 类 | 职责 |
|----|------|
| `KernelEstimator` | 核函数估计器，三个正交维度 |
| `KernelResult` | 核函数结果容器，支持序列化 (`save`/`load`) |

### 1.2 三维设计空间

| 维度 | 可选值 | 默认 | 含义 |
|------|--------|------|------|
| `mode` | `'flux'`, `'omega'` | `'flux'` | 刺激物理量：磁通 $\Phi$（$\Phi_0$）或频率 $\omega$（rad·GHz） |
| `method` | `'sim'`, `'exp'` | `'exp'` | 模拟策略：纯理论传播子 或 实验协议仿真 |
| `order` | `1, 2, 3, \ldots` | `1` | Volterra 阶数（1=线性，≥2=非线性对角核） |

### 1.3 合法组合

| (mode, method) | qubit 需求 | order=1 | order≥2 |
|:---:|:---:|:---:|:---:|
| (omega, sim) | 可选 | ✓ Heisenberg 传播子 | ✓ Heisenberg 对易子 |
| (omega, exp) | 必需 | ✓ Virtual Z 双边差分 | ✓ 5-pt FD stencil |
| (flux, exp) | 必需 | ✓ 单边差分（legacy） | ✓ 5-pt FD stencil |
| (flux, sim) | — | ❌ 非法 | ❌ 非法 |

### 1.4 两套计算路径（Phase 10-fix 后）

```
method='sim' ──→ _heisenberg_kernels()
                 │  1 次 sesolve(qeye(n), t_grid)
                 │  ├─ U(t) 传播子 (n×n 矩阵，所有时间点)
                 │  ├─ W(t) = U†(t) σ_z U(t) / 2
                 │  └─ k_n = i^n ⟨0| ad_W^n(Q) |0⟩  (嵌套对易子)
                 │
                 │  代价: O(1 sesolve + N_t × d³ 矩阵乘)
                 │  精度: 机器精度 (无近似、无拟合)
                 │  阶数: 任意阶同代价

method='exp' ──→ order=1:
                 │  flux:  _estimate_flux()   单边高斯刺激 (legacy)
                 │  omega: _estimate_omega()  Virtual Z 双边差分
                 │
                 order≥2:
                    flux:  _extract_kn_flux()   5-pt FD stencil (FluxSignal)
                    omega: _extract_kn_omega()  5-pt FD stencil (Virtual Z)
                    │
                    代价: 5 × N_t 次 mesolve
                    精度: k₁ O(h⁴), k₂ O(h⁴), k₃ O(h²)
```

### 1.5 5 点 FD stencil 公式

在 $\phi_z \in \{-2h, -h, 0, +h, +2h\}$ 处采样 $p_e$：

| 阶数 | 公式 | 精度 |
|------|------|------|
| $k_1$ | $(f_{-2} - 8f_{-1} + 8f_{+1} - f_{+2}) / 12h$ | $O(h^4)$ |
| $k_2$ | $(-f_{-2} + 16f_{-1} - 30f_0 + 16f_{+1} - f_{+2}) / 12h^2$ | $O(h^4)$ |
| $k_3$ | $(-f_{-2} + 2f_{-1} - 2f_{+1} + f_{+2}) / 2h^3$ | $O(h^2)$ |

### 1.6 核函数的物理含义

$$p_e = p_e^{(0)} + \int k_1(t)\,\delta\omega(t)\,dt + \frac12 \iint k_2(t_1,t_2)\,\delta\omega(t_1)\delta\omega(t_2)\,dt_1dt_2 + \cdots$$

- **$k_1(t)$**: 线性灵敏度——在时刻 $t$ 注入单位频率微扰引起的 $p_e$ 变化
- **$k_2(t,t)$**: 二阶非线性——两个同时刻的频率微扰产生的二次响应
- **$k_3(t,t,t)$**: 三阶非线性

**特殊性质**（零 detuning Ramsey π/2-π/2）：
- $k_2(t) \equiv 0$（严格，Y-X 正交脉冲的对称性保证）
- $|G_3/G_1| = 1$（三阶与一阶同量级）

### 1.7 核函数单位

| mode | k₁ 单位 | k₂ 单位 | 积分 $G_n = \int k_n dt$ |
|------|---------|---------|--------------------------|
| `'omega'` | rad⁻¹（无量纲） | rad⁻² | $G_1$ 无量纲 |
| `'flux'` | $1/(\Phi_0\cdot\text{ns})$ | $1/(\Phi_0^2\cdot\text{ns})$ | $G_1$ 量纲 $1/\Phi_0$ |

$\omega$/flux 换算（线性色散）: $k_n^{(\Phi)} = \kappa^n \, k_n^{(\omega)}$，$\kappa = d\omega/d\Phi$。

### 1.8 序列化

```python
result = estimator.estimate_full(pulse, qubit)
result.save('kernel.npz')          # numpy.savez, 逐 kernel 存储
loaded = KernelResult.load('kernel.npz')
loaded.k1                           # → k₁ ndarray
loaded.kernels[1]                   # → k₂ ndarray (order≥2)
```

---

## 二、频率标定体系 (`sqc/calibration/frequency.py`)

### 2.1 核心类

| 类 | 职责 |
|----|------|
| `FrequencyMeasurement` | 单点频率测量（统一入口） |
| `FluxResponseCalibration` | $f(\Phi)$ 曲线的多点标定 |
| `SinglePointFrequencyCalibration` | 闭环调谐（组合 `FrequencyMeasurement`） |

### 2.2 两种测频方法

| 方法 | 原理 | 时间 | 适用范围 | 依赖 |
|------|------|------|---------|------|
| `"transient"` | τ=0 正交 Ramsey + 核函数 $G_1$ | ~0.2 s | $|\Delta\omega|<0.05$ GHz | `KernelEstimator(mode='omega', order=1)` |
| `"ramsey"` | τ 扫描 + FFT 峰值检测 | 4–50 s | 全范围 | 无需核函数 |

### 2.3 Transient 测频管线（`_measure_frequency_transient`）

```
1. 构造正交 Ramsey 脉冲
   ctrl_y  = R_y(π/2), τ=0
   ctrl_my = R_{-y}(π/2), τ=0

2. 仿真 p_diff = (p_y - p_{-y}) / 2
   (通过 mesolve 计算末态 p_e)

3. 估计 omega 核函数
   estimator = KernelEstimator(mode='omega', method='exp', order=1)
   result = estimator.estimate_full(ctrl, qubit)
   G_freq = ∫ k₁(t) dt            ← 总灵敏度

4. 频率反演
   Δω = p_diff / G_freq           ← 直接！无需 κ
   f_meas = ω_drive - Δω
```

**关键**：Phase 10.5 后使用 omega kernel 直接路径，消除了旧的 $\kappa$ workaround（`delta_omega = p_diff * kappa / G_diff`）。

### 2.4 Ramsey τ-sweep 测频管线（`_fit_ramsey_frequency`）

```
1. 设定人工 detuning f_art = 0.1 GHz
   (保证 f_art > |Δ|，使频率符号可辨)

2. 扫描 τ ∈ [0, 200] ns
   每个 τ: phase2 = 2π·f_art·τ → Ramsey 脉冲 → mesolve → p_e

3. FFT 峰值检测
   p_e(τ) → 去均值 → 汉宁窗 → FFT → 正频率峰值 f_meas_raw
   Δ = f_art - f_meas_raw          ← 单扫模式
   f_meas = ω_drive + Δ

   双扫模式 (f_art=None): 两次扫描 (±50 MHz) → 无符号歧义
```

### 2.5 精度-速度 Trade-off

| |Δω| (GHz) | Transient 误差 | Ramsey 误差 | Transient 时间 | Ramsey 时间 | 推荐 |
|-----------|---------------|--------------|-------------|-------------|------|
| 0.002 | 9×10⁻⁵ | 7×10⁻⁴ | 0.2 s | 4.2 s | Transient |
| 0.03 | 6×10⁻⁴ | 1×10⁻² | 0.2 s | 4.5 s | Transient |
| **0.07** | **5×10⁻³** | **2×10⁻²** | **0.2 s** | **5.6 s** | **交叉点** |
| 0.19 | 1×10⁻¹ | 7×10⁻² | 0.2 s | 8.3 s | Ramsey |
| 1.73 | 1.73 | 1.07 | 0.2 s | 50.7 s | Ramsey |

- **交叉点** $|\Delta\omega| \approx 0.07$ GHz
- **小 detuning**: Transient 快 20–50× 且更准（FFT 低频分辨率不足）
- **大 detuning**: Ramsey τ-sweep 不依赖线性近似，全范围鲁棒

### 2.6 推荐两级策略

```
Step 1: FrequencyMeasurement(method="transient")  ~0.2 s
    → 如果 p_diff < 0.3 (线性安全区): 完成
    → 如果 p_diff >= 0.3: 转到 Step 2

Step 2: FrequencyMeasurement(method="ramsey")     5–50 s
    → FFT 全范围可靠
    → 可闭环迭代细化
```

---

## 三、核函数 ↔ 频率标定的关系

```
              ┌─────────────────────────────┐
              │     FrequencyMeasurement     │
              │  method="transient"          │
              └─────────────┬───────────────┘
                            │
              ┌─────────────▼───────────────┐
              │      KernelEstimator         │
              │  mode='omega'               │
              │  method='exp'               │
              │  order=1                    │
              │  virtual_z_impl='math'      │
              └─────────────┬───────────────┘
                            │
              ┌─────────────▼───────────────┐
              │   _estimate_omega()          │
              │   双边差分: k₁ = Δp/Δφ_z    │
              │   (Virtual Z 高斯 σ_z 冲激) │
              └─────────────┬───────────────┘
                            │
              ┌─────────────▼───────────────┐
              │   G_freq = ∫ k₁(t) dt       │
              │   Δω = p_diff / G_freq      │
              └─────────────────────────────┘
```

**依赖关系**：
- `FrequencyMeasurement(method="transient")` **依赖** `KernelEstimator`
- `FrequencyMeasurement(method="ramsey")` **不依赖**核函数体系
- `KernelEstimator` **不依赖**任何 calibration 类（符合重构依赖矩阵）

---

## 四、与 Gao 2021 全栈架构的对应

| Gao 2021 层 | 本项目模块 | 核函数/频率角色 |
|-------------|-----------|---------------|
| **Control software** (§V) | `sqc/experiments/` + `sqc/calibration/` | 实验协议选择、频率标定 workflow |
| **Control electronics** (§IV.B) | `sqc/hardware/` | Virtual Z 的硬件实现（相位寄存器） |
| **Microwave signal processing** (§IV.B) | `sqc/reconstruction/kernel.py` | 核函数 = IQ 解调前的前端灵敏度模型 |
| **Device** (§II) | `sqc/devices/transmon.py` | $\omega(\Phi)$ 色散关系、$\kappa$ 灵敏度 |

---

## 五、已有能力 & 限制

### 已有

| 能力 | 实现 |
|------|------|
| 一阶 omega 核函数（sim + exp） | ✓ 双边差分 / Heisenberg 传播子 |
| 一阶 flux 核函数（exp） | ✓ 单边高斯差分 (legacy) |
| 二阶/三阶 omega 核函数（sim） | ✓ Heisenberg 对易子 |
| 二阶/三阶 omega 核函数（exp） | ✓ 5-pt FD stencil |
| 二阶/三阶 flux 核函数（exp） | ✓ 5-pt FD stencil |
| 核函数序列化 | ✓ KernelResult.save/load (.npz) |
| Transient 单点测频 | ✓ _measure_frequency_transient() |
| Ramsey τ-sweep 测频 | ✓ _fit_ramsey_frequency() |
| 闭环频率调谐 | ✓ SinglePointFrequencyCalibration |
| $f(\Phi)$ 曲线标定 | ✓ FluxResponseCalibration(method="ramsey") |

### 限制

| 限制 | 影响 |
|------|------|
| k₂ ≈ 0 在零 detuning（物理严格） | 二阶修正对小信号无帮助 |
| FD stencil 的 k₃ 精度 O(h²) | 需要合理选择 h (~0.005 rad) |
| Heisenberg 仅 sim（需知道 H(t)） | 不能用于实验数据 |
| Exp 高阶仍比 order-1 慢 2.5× | 日常测频建议 order=1 |
| extract_off_diagonal 未实现 | 完整的 k₂(t_i, t_j) 尚不可用 |

---

## 六、相关文件索引

| 文件 | 角色 |
|------|------|
| `sqc/reconstruction/kernel.py` | KernelEstimator + KernelResult（核函数体系核心） |
| `sqc/calibration/frequency.py` | FrequencyMeasurement + 两种测频管线 |
| `sqc/control/pulse.py` | Pulse.get_kernel() legacy shim |
| `sqc/reconstruction/transient.py` | TransientReconstruction（反卷积，含 Hammerstein-Volterra） |
| `idea/_sensing theory.md` | 物理理论：核函数定义、Volterra 展开、对易子公式 |
| `idea/refactor/phase_10_kernel_extension_handbook.md` | Phase 10 设计文档（API 设计、合法组合矩阵） |
| `idea/refactor/phase_10_kernel_theory_and_fix.md` | 失败分析 + 修复方案的数学推导 |
| `kernel_order_comparison.ipynb` | 实验 notebook：Transient vs Ramsey τ-sweep 三方对比 |
