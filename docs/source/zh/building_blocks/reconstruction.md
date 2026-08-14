# 重建(reconstruction)

## 这层提供什么

`reconstruction` 层是全栈的反问题层,把 {doc}`experiments` 层产出的原始测量
数据(布居 $p_e$、相位 $\varphi$、控制核等)反演回物理磁通波形 $\Phi(t)$ 或
磁场 $B(t)$。它是波形重建这条产品主线的算法核心:从比特测得的响应倒推外界施加
的信号。

一条贯穿本层的纯函数约束:`reconstruct()` 只读测量数据,不运行新仿真。
若某算法需要正向仿真(如 LM 全密度矩阵反演),所需的 `qubit` / `control_pulse`
在构造时注入,而非在 `reconstruct()` 里临时构造。输出统一是物理单位的
{py:class}`~sqc.control.FluxSignal`(或 `np.ndarray` 形式的 $B$)。

```{note}
本层含两个标定类 `CryoscopeCalibration` / `DelayRamseyCalibration`。它们在
`sqc.reconstruction` 命名空间下导出(与对应重建类配套),但概念上属于标定:
继承 {py:class}`~sqc.calibration.Calibration`,产出 `CalibrationTable` 供重建时
`inversion="calibration"` 查表反演。详见 {doc}`calibration`。
```

## 类总览

按功能分四组:

**扩展点**

| 类 | 角色 |
|---|---|
| `Reconstruction` | 抽象基类,所有重建器的公共契约,本层扩展点 |

**基函数工具**(供 LM 反演等参数化用)

| 函数 | 角色 |
|---|---|
| `generate_basis_functions` | 生成 B-spline / Fourier / Legendre 基函数列表 |
| `basis_function_decomposition` | 最小二乘把信号分解到基函数上 |
| `regularization_matrix` / `R` | 鼓励平滑的正则化矩阵(`R` 是别名) |

**控制核估计**(瞬态反卷积用)

| 类 | 角色 |
|---|---|
| `KernelEstimator` | 用窄高斯激励逐点探测,估计线性/高阶 Volterra 控制核 |
| `KernelResult` | 核估计结果容器(`k1` 线性核、高阶核、`save`/`load`) |

**各协议重建器**(各接收对应实验的 `ExperimentResult`)

| 类 | 角色 | 反演公式/方法 |
|---|---|---|
| `RamseyReconstruction` | Ramsey 重建 | IQ / unwrap 相位 → $B=\dot\varphi/\kappa$ |
| `EchoReconstruction` | 差分回波重建 | $B=-\varphi/(2k\kappa t_\mathrm{int})$ |
| `TransientReconstruction` | 瞬态场重建 | Wiener / Hammerstein / LM 反卷积 |
| `CryoscopeReconstruction` | Cryoscope 重建 | $\dot\varphi\to h(t)$,标定表或色散反演 |
| `DelayRamseyReconstruction` | 延迟 Ramsey 重建 | $\varphi(t_d)\to\Phi_\mathrm{tail}$ |
| `PiPulseCompReconstruction` | π-脉冲补偿重建 | $\Phi_\mathrm{tail}=-z^*(\tau)$ |
| `CryoscopeCalibration` | (标定) | 方波扫 $h$ 建 $\varphi(h)$ 查表 |
| `DelayRamseyCalibration` | (标定) | 扫 $z$ 拟合 $\varphi_\mathrm{cal}(z)=\tau_R\kappa z$ |

## Reconstruction：重建抽象基类

所有重建器的公共契约,本层扩展点。只规定一个抽象方法:

- `reconstruct(*args, **kwargs)`:从测量数据反演物理信号。

两条不变量(继承时必须遵守):

1. `reconstruct()` 必须是纯函数,只读测量、不运行新仿真;需要正向仿真的算法(如 LM)
   把 `qubit`/`control_pulse` 在构造时注入。
2. 输出永远是物理单位的 `FluxSignal`(或 `np.ndarray` 形式的场)。

要添加自定义重建算法,继承 `Reconstruction` 并实现 `reconstruct()`,详见
{doc}`../extending`。

## 基函数工具

一组把信号参数化到基函数上的纯函数,供 LM 反演等做低维参数化与正则化用。

- `generate_basis_functions(basis_type, n_basis, t_min, t_max) -> list[callable]`:
  在 $[t_\min, t_\max]$ 上生成 `n_basis` 个基函数。`basis_type` 取
  `"bspline"`(3 阶 B-spline)、`"fourier"`(正交正弦余弦)、`"legendre"`
  (Legendre 多项式)。返回一列可调用对象,各接受时间数组、返回基值。
- `basis_function_decomposition(sig, t_array, basis_functions) -> np.ndarray`:
  最小二乘把信号 `sig`(数组,或带 `.signal` 属性的 `FluxSignal`)分解到基上,
  返回 `(n_basis,)` 系数。
- `regularization_matrix(n, basis_type) -> np.ndarray`:鼓励平滑的 $n\times n$
  正则化矩阵,Fourier 基按频率阶数平方惩罚高频,其余基用二阶差分逼近二阶导。
  `R` 是它的向后兼容别名。

## 控制核估计

控制核是瞬态反卷积的核心输入：它刻画了"在控制脉冲的某一时刻注入一个微小
扰动，对最终比特布居数 $p_e$ 产生多大影响"。`KernelEstimator` 在这个"扰动
—响应"映射上支持一个三维设计空间，覆盖从快速一阶线性核到完整三阶
off-diagonal Volterra 张量的全部路径。

### KernelEstimator

**构造**

`KernelEstimator(mode="omega", method="exp", order=1, virtual_z_impl="math",
extract_off_diagonal=False, ...)`

**全部字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `mode` | `str` | 扰动通道：`"flux"`（磁通）或 `"omega"`（Virtual Z 频率） | `"flux"` |
| `method` | `str` | 估计策略：`"exp"`（实验式仿真，含 qubit 色散）或 `"sim"`（纯理论 Heisenberg，无 qubit 色散） | `"exp"` |
| `order` | `int` | 最高 Volterra 阶数（1=仅线性核，≥2=高阶核） | `1` |
| `extract_off_diagonal` | `bool` | 是否提取完整 n-D off-diagonal 张量（仅 LM 可消费） | `False` |
| `virtual_z_impl` | `str` | VZ 实现方式：`"math"`（σ_z 高斯脉冲）或 `"hardware"`（相位旋转子脉冲，仅 CompositePulse） | `"math"` |
| `stim_amplitude` | `float` | 探头激励幅度（flux 模式为 Φ₀，omega 模式为 rad） | `0.0215` |
| `stim_width` | `float` | 探头宽度（ns，flux 模式用） | `3.0` |
| `auto_calibrate` | `bool` | 按 qubit 非谐性自动调整幅度（flux 模式） | `False` |
| `n_levels` | `int` | Hilbert 空间截断（sim 模式用，exp 模式从 qubit 取） | `2` |
| `anharmonicity` | `float` | 非谐性 α（GHz），n_levels≥3 且 method="sim" 时用 | `0.0` |
| `kappa` | `float` | 线性色散 dω/dΦ（GHz/Φ₀），仅 sim 模式 flux→omega 转换时用 | `None` |
| `n_amp_samples` | `int` | 高阶振幅扫描点数（已弃用，现用固定 5 点 FD stencil） | `5` |
| `amp_scan_factor` | `float` | 振幅扫描范围乘子 | `1.0` |
| `probe_sigma_t` | `float` 或 `None` | VZ 高斯探头宽度 σ_t（ns）；`None`→`2·dt`。越小 off-diagonal 涂抹偏差越小，但须 ≥dt | `None` |
| `richardson` | `bool` | 是否启用 Richardson σ_t→0 外推消除探头宽度偏差（仅 omega/math 模式有效） | `False` |
| `richardson_sigmas` | `tuple` 或 `None` | Richardson 外推的 σ_t 采样点（×dt）；`None`→`(2.0, 1.5, 1.0)` | `None` |

**方法**

- `estimate(pulse, qubit=None, t_samples=None) -> (t_samples, k1)`：
  返回时间轴与一阶线性核 $k_1$（向后兼容）。
- `estimate_full(pulse, qubit=None, t_samples=None) -> KernelResult`：
  返回完整结果对象（含所有阶核）。

### 设计空间：完整调度表

三维设计空间由 `(mode, method, order)` 以及 `extract_off_diagonal` 和
`virtual_z_impl` 标志构成。合法组合及其内部实现路径：

| mode | method | order | off-diag | vz_impl | 内部方法 | 复杂度 | 说明 |
|---|---|---|---|---|---|---|---|
| flux | exp | 1 | — | — | `_estimate_flux` | M×1 mesolve | 窄高斯磁通，单边差分 |
| omega | exp | 1 | — | math | `_omega_vz_math` | M×2 mesolve | σ_z 高斯脉冲，双边差分 |
| omega | exp | 1 | — | hardware | `_omega_vz_hardware` | M×2 mesolve | 相位旋转子脉冲 |
| omega | sim | 1 | — | — | `_heisenberg_kernels` | 1 sesolve | 精确对易子公式 |
| flux | exp | ≥2 | 否 | — | `_extract_kn_flux` | M×5 mesolve | 5 点 FD stencil |
| omega | exp | ≥2 | 否 | math | `_extract_kn_omega` | M×5 (或 ×15 Richardson) mesolve | 5 点 FD + Richardson 可选 |
| omega | exp | ≥2 | 否 | hardware | `_extract_kn_omega`（hw 分支） | M×5 mesolve | 相位旋转子脉冲，5 点 FD |
| omega | exp | ≥2 | 是 | math | `_extract_kn_offdiag_exp` | O(M²)~O(M³) mesolve | 混合偏导 FD，多 σ_z 同时脉冲 |
| omega | exp | ≥2 | 是 | hardware | `_extract_kn_offdiag_exp`（hw 分支） | O(M²)~O(M³) mesolve | 多时刻同时相位踢（`with_phase_kicks`） |
| omega | sim | ≥2 | 否 | — | `_heisenberg_kernels` | 1 sesolve | 嵌套对易子，任意阶 |
| omega | sim | ≥2 | 是 | — | `_heisenberg_kernels_offdiag` | 1 sesolve + O(Mⁿ) 矩阵乘法 | 完整 n-D 张量，≤3 阶 |
| flux | sim | any | any | — | **非法** | — | 直接抛 `ValueError` |

### 扰动模式详解

**flux（磁通）**：窄高斯磁通信号（`FluxSignal(type=3)`）通过
`qubit.qubit_under_mag()` 注入 qubit 色散模型。核单位为 `1/(Φ₀·ns)`。
探头宽度由 `stim_width` 控制。

**omega / math VZ（σ_z 高斯脉冲）**：在哈密顿量中加上
`(φ_z/2)·σ_z·gaussian(t−tⱼ)` 项。物理上等价于瞬间频率偏移，但有限
σ_t 会涂抹 off-diagonal 结构（k₃ 在 σ_t=2·dt 时偏低 ~15%）。
Richardson 外推（`richardson=True`）在多个 σ_t 上采样后外推至 σ_t→0
可消除此偏差。核单位为 `rad⁻¹`（无量纲）。

**omega / hardware VZ（相位旋转子脉冲）**：不往 H 里加东西，而是重建
整个脉冲序列——在探测时刻 tⱼ 之后的所有子脉冲相位旋转 ±φ_z。这是
实验中最接近 AWG 相位寄存器操作的方式。对 `CompositePulse` 有效，单
`Pulse` 会回退到 math 实现。相位踢是瞬时的（无 σ_t），因此 Richardson
外推不适用。

### 一阶核（order=1）

快速路径：每个探测点仅需 1–2 次 `mesolve`（sim 只需 1 次总 `sesolve`）。
适合大多数线性 Wiener 反卷积场景。

```python
# 一阶 flux（向后兼容）
est = KernelEstimator(mode='flux', method='exp')
t, k1 = est.estimate(pulse, qubit)

# 一阶 hardware VZ（实验级真实度）
est = KernelEstimator(mode='omega', method='exp', virtual_z_impl='hardware')
t, k1 = est.estimate(pulse, qubit)

# 一阶 sim（精确理论核，极快）
est = KernelEstimator(mode='omega', method='sim')
t, k1 = est.estimate(pulse)  # 无需 qubit
```

### 高阶对角核（order≥2, off_diagonal=False）

5 点中心差分 stencil，在 $\{-2h, -h, 0, +h, +2h\}$ 五个幅值上采样：

- $k_1 = (f_{-2} - 8f_{-1} + 8f_{+1} - f_{+2}) / 12h$
- $k_2 = (-f_{-2} + 16f_{-1} - 30f_0 + 16f_{+1} - f_{+2}) / 12h^2$
- $k_3 = (-f_{-2} + 2f_{-1} - 2f_{+1} + f_{+2}) / 2h^3$

每个探测点 5 次 `mesolve`；Richardson 用 3 个 σ_t 则 15 次。输出为 1-D
对角核 $k_n(t, t, …, t)$。Wiener / Hammerstein 反卷积仅消费对角核。

```python
# 三阶 exp/math 对角 + Richardson
est = KernelEstimator(mode='omega', method='exp', order=3,
                      richardson=True)
res = est.estimate_full(pulse, qubit)
k1, k2, k3 = res.kernels  # 均为 1-D

# 三阶 exp/hardware 对角
est = KernelEstimator(mode='omega', method='exp', order=3,
                      virtual_z_impl='hardware')
res = est.estimate_full(pulse, qubit)
```

### 高阶 off-diagonal 核（extract_off_diagonal=True）

返回完整 n-D Volterra 张量：$k_2$ 为 `(M,M)`，$k_3$ 为 `(M,M,M)`。
**仅有 LM 全密度矩阵反演能消费 off-diagonal 核**；Wiener / Hammerstein
路径只接受对角核。

**Sim 路径**（Heisenberg）：精确解，一次 `sesolve`，零微扰。对 $k_3$
计算 $M^3$ 个对易子，$M$ 大时较慢。

**Exp 路径**（有限差分）：混合偏导 FD 推广了 5 点 stencil——$k_2(i,j)$
用四角混合差 $[p(+,+)-p(+,-)-p(-,+)+p(-,-)]/4h^2$；$k_3(i,j,l)$ 按
退化情况分全等 / 两等一异 / 全异三类。含 $p_e$ 缓存避免重复计算。

**Hardware VZ**（`virtual_z_impl='hardware'`）：同时打多个相位踢，通过
`pulse.with_phase_kicks([(t_i,+h), (t_j,-h)])` 一次重建含多个跳变的脉冲。
单 `Pulse` 回退到 math 实现。

```python
# sim off-diagonal（精确，推荐用于 G₃ 标定）
est = KernelEstimator(mode='omega', method='sim', order=3,
                      extract_off_diagonal=True)
res = est.estimate_full(pulse, qubit)
k1 = res.kernels[0]   # (M,)
k3 = res.kernels[2]   # (M, M, M) — 完整三阶张量

# exp/hardware off-diagonal（实验级真实度）
est = KernelEstimator(mode='omega', method='exp', order=2,
                      extract_off_diagonal=True,
                      virtual_z_impl='hardware')
res = est.estimate_full(pulse, qubit)
k1 = res.kernels[0]   # (M,)
k2 = res.kernels[1]   # (M, M)
```

### 下游消费约束

| 核类型 | Wiener | Hammerstein | Hammerstein-Volterra | LM |
|---|---|---|---|---|
| 一维对角（$k_1$） | ✓ | ✓ | ✓ | ✓ |
| 一维对角（$k_1,k_2,k_3$） | ✗ | ✗ | ✓ | ✓ |
| 多维 off-diagonal | ✗ | ✗ | ✗ | ✓ |

### KernelResult

核估计结果容器。`kernels[n-1]` 为第 n 阶核：`kernels[0]`=$k_1$
（1-D），若提取则 `kernels[1]`=$k_2$（对角 1-D 或 off-diag 2-D），
`kernels[2]`=$k_3$。

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `t_samples` | `np.ndarray` | 核采样时间轴 (ns) |
| `kernels` | `list[np.ndarray]` | 各阶核列表，`kernels[0]`=k₁ |
| `mode` | `str` | 扰动通道 |
| `method` | `str` | 估计方法 |
| `order` | `int` | 最高核阶数 |
| `stim_amplitude` | `float` | 激励幅度 |
| `units` | `str` | 核单位（flux: `1/(Φ₀·ns)`, omega: `rad⁻¹`） |
| `off_diagonal` | `bool` | `True` 表示高阶核为多维张量 |

**属性**：`k1`（=`kernels[0]`，线性核快捷访问）。

**方法**：`save(path)` 序列化到 `.npz` / `load(path)`（类方法）反序列化。

## RamseyReconstruction：Ramsey 重建

从 Ramsey 相位反演磁场 $B(\tau)$。

**构造**

`RamseyReconstruction(qubit, method="unwrap", k_span=3)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 比特对象(提供频率灵敏度 $\kappa$) | —— |
| `method` | str | 相位提取路径 | `"unwrap"`(或 `"iq"`) |
| `k_span` | int | unwrap 分支搜索窗半宽 | `3` |

**方法**

- `reconstruct(measurement) -> np.ndarray`：

  - `method="iq"`:双通道 $\arctan2$ 取相位(需 `data["p_e_I"]`/`p_e_Q`)。
  - `method="unwrap"`:单通道 $\arccos$ + $k$-span 分支解缠绕(需 `data["p_e"]`)。

  两条路径都对相位求梯度、除以频率灵敏度 $\kappa=\mathrm{d}\omega/\mathrm{d}\Phi$ 得
  $B=\dot\varphi/\kappa$。

**输出**

返回 `np.ndarray`,即磁场 $B(\tau)$ 数组。

## EchoReconstruction：差分回波重建

从差分回波 $p_e$ 直接解析反演磁场。

**构造**

`EchoReconstruction(qubit, t_int, k)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 比特对象 | —— |
| `t_int` | float | 每个回波块的相互作用时间(ns) | —— |
| `k` | int | $\pi$ 脉冲对数 | —— |

**方法**

- `reconstruct(measurement) -> np.ndarray`：从 `data["p_e"]` 取
  $\varphi=\arcsin(2p_e-1)$,$B=-\varphi/(2k\kappa t_\mathrm{int})$。

**输出**

返回 `np.ndarray`,即磁场 $B$ 数组。

```{note}
`arcsin`/`arccos` 前统一 `clip` 到 $[-1,1]$:`mesolve` 积分噪声可能把 $p_e$ 顶出
$[0,1]$ 微许,不裁会产生静默 NaN。Ramsey/Echo/Transient 各路都有同款保护。
```

## TransientReconstruction：瞬态场重建

波形重建主线的核心反卷积器,把瞬态实验的 $\Delta p$ 与控制核反卷积回 $\Phi(t)$。

**构造**

`TransientReconstruction(method="wiener", lambda_reg=1e-3, qubit=None, control_pulse=None, basis_type="bspline", n_basis=8, max_iter=50, tol=1e-6, mu_init=0.01, use_adjoint=False, max_volterra_iter=5, volterra_tol=1e-4)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `method` | str | 反演算法 | `"wiener"`(或 `"hammerstein"`/`"hammerstein_volterra"`/`"lm"`) |
| `lambda_reg` | float | Wiener 正则化参数 | `1e-3` |
| `qubit` | `TransmonQubit` | 比特对象(Hammerstein/LM 需要) | `None` |
| `control_pulse` | `CompositePulse` | 控制脉冲(LM 需要) | `None` |
| `basis_type` | str | LM 基函数类型 | `"bspline"` |
| `n_basis` | int | LM 基函数数量 | `8` |
| `max_iter` | int | LM 最大迭代次数 | `50` |
| `tol` | float | LM 收敛容差 | `1e-6` |
| `mu_init` | float | LM 初始阻尼 | `0.01` |
| `use_adjoint` | bool | LM 使用伴随法 Jacobian | `False` |
| `max_volterra_iter` | int | Hammerstein-Volterra 最大迭代 | `5` |
| `volterra_tol` | float | Volterra 收敛容差 | `1e-4` |

**方法**

- `reconstruct(measurement, kernel=None, **kwargs) -> FluxSignal`：

  - `"wiener"`:线性 Wiener 反卷积(FFT 域)。
  - `"hammerstein"`:Hammerstein-Wiener 非线性块模型(需 `qubit`)。
  - `"hammerstein_volterra"`:用高阶核的迭代反演(需 `order>=2` 的核)。
  - `"lm"`:Levenberg-Marquardt 全密度矩阵反演(需 `control_pulse`)。

  `kernel` 可传 `np.ndarray`(线性核)或 `KernelResult`(含高阶核)。非对角多维核
  只有 `method="lm"` 能消费,Wiener/Hammerstein 路径只接受对角(一维)核,否则抛
  `ValueError`。

**输出**

返回 `FluxSignal`,即重建的磁通波形 $\Phi(t)$。

## CryoscopeReconstruction：Cryoscope 重建

从相位-截断数据重建波形,核心物理 $\mathrm{d}\varphi/\mathrm{d}t=2\pi\Delta f(h(t))$。

**构造**

`CryoscopeReconstruction(tau=None, inversion="calibration", calibration=None, qubit=None, use_sg_filter=True, sg_window=11, sg_poly=3)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `tau` | float | 标定方波长度(ns) | —— |
| `inversion` | str | 反演路径 | `"calibration"`(或 `"response"`) |
| `calibration` | `CalibrationTable` | $\varphi(h)$ 标定表(`inversion="calibration"`时必传) | `None` |
| `qubit` | `TransmonQubit` | 比特对象(`inversion="response"`时必传) | `None` |
| `use_sg_filter` | bool | 是否用 Savitzky-Golay 预平滑 | `True` |
| `sg_window` | int | SG 窗口宽度 | `11` |
| `sg_poly` | int | SG 多项式阶数 | `3` |

**方法**

- `reconstruct(measurement, ...) -> FluxSignal`：读 `data["varphi"]`、`axes["trunc"]`。
  - `inversion="calibration"`:用 `CryoscopeCalibration` 产的 $\varphi(h)$ 查表反演。
  - `inversion="response"`:用解析 Transmon 频率-磁通色散反演。

**输出**

返回 `FluxSignal`,即重建的磁通波形。

## DelayRamseyReconstruction：延迟 Ramsey 重建

把延迟 Ramsey 相位 $\varphi(t_d)$ 反演成拖尾磁通波形(预畸变主线)。

**构造**

`DelayRamseyReconstruction(inversion="response", qubit=None, calibration=None, tau_R=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `inversion` | str | 反演路径 | `"response"`(或 `"calibration"`) |
| `qubit` | `TransmonQubit` | 比特对象(`inversion="response"`时必传) | `None` |
| `calibration` | `CalibrationTable` | $\varphi(z)$ 标定表(`inversion="calibration"`时必传) | `None` |
| `tau_R` | float | Ramsey 自由演化时间(ns) | —— |

**方法**

- `reconstruct(measurement, ...) -> FluxSignal`：读 `data["varphi"]`、`axes["t_d"]`。

**输出**

返回 `FluxSignal`,即拖尾磁通波形 $\Phi_\mathrm{tail}$。

## PiPulseCompReconstruction：π-脉冲补偿重建

直接的重建:2D 扫描给出最优补偿高度 $z^*(\tau)$ 后,拖尾即 $\Phi_\mathrm{tail}=-z^*(\tau)$。

**构造**

`PiPulseCompReconstruction()`

(无额外字段——仅需 `measurement` 即可反演。)

**方法**

- `reconstruct(measurement) -> FluxSignal`：读 `data["z_star"]`、`axes["tau"]`。

**输出**

返回 `FluxSignal`,即拖尾磁通波形。

```{note}
因 $\pi$ 脉冲有限宽 $T$,$z^*(\tau)$ 测的是拖尾在 $[\tau, \tau+T]$ 上的滑动平均,
而非瞬时值。对指数拖尾这只带来常数衰减因子,衰减时间常数不变。
```

## 标定类(CryoscopeCalibration / DelayRamseyCalibration)

两个在本命名空间导出、但概念上属于标定的类(继承
{py:class}`~sqc.calibration.Calibration`),为对应重建器的 `"calibration"` 反演路径
建查表:

### CryoscopeCalibration

**构造**

`CryoscopeCalibration(qubit, t_rabi=None, tau=20.0)`

**方法**

- `calibrate() -> CalibrationTable`：扫方波高度 $h$、IQ 读出 + 模型引导解缠绕建
  $\varphi(h)$ 查表,`kind="phi_h"`,产出供 `CryoscopeReconstruction`(inversion="calibration")用。

### DelayRamseyCalibration

**构造**

`DelayRamseyCalibration(qubit, t_rabi=None, tau_R=None)`

**方法**

- `calibrate() -> CalibrationTable`：扫已知磁通高度 $z$、拟合线性斜率
  $\varphi_\mathrm{cal}(z)=\tau_R\kappa z$,产出 `kind="phi_z"`,供
  `DelayRamseyReconstruction`(inversion="calibration")用。

完整标定契约见 {doc}`calibration`。

## 最小用例

```python
import numpy as np
from sqc.devices import TransmonQubit
from sqc.simulation import ExperimentResult
from sqc.reconstruction import (
    RamseyReconstruction, generate_basis_functions, basis_function_decomposition,
)

qubit = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10_000, T2=8_000, flux=0.1)

# 1) 合成一段 Ramsey 测量(实际来自 RamseyExperiment.run())
tau = np.linspace(0, 250, 60)
p_e = 0.5 * (1 - np.cos(0.02 * tau))
meas = ExperimentResult(data={"p_e": p_e}, axes={"tau": tau})

# 2) 纯后处理反演磁场 B(τ)——不跑新仿真
B = RamseyReconstruction(qubit=qubit, method="unwrap").reconstruct(meas)

# 3) 把信号分解到 B-spline 基上(LM 反演等的参数化)
basis = generate_basis_functions("bspline", 8, tau[0], tau[-1])
coeffs = basis_function_decomposition(p_e, tau, basis)   # (8,) 系数
```

## 物理角色 / 扩展

- 本层对应实验数据的离线分析:把测到的布居/相位翻译回真正施加的磁通/磁场。
- 各重建器与 {doc}`experiments` 的实验类一一配对:Ramsey↔Ramsey、Echo↔差分回波、
  Transient↔瞬态、Cryoscope/延迟 Ramsey/π-补偿↔对应拖尾测量。
- `KernelEstimator` 产的控制核是瞬态反卷积的必需输入;基函数工具支撑 LM 的低维
  参数化;两个标定类为查表反演路径供 `CalibrationTable`。
- 要添加自定义重建算法,继承 `Reconstruction` 抽象基类实现 `reconstruct()`,
  遵守纯函数、输出 `FluxSignal` 两条不变量,需要正向仿真则在构造时注入 `qubit`。
  完整扩展指南见 {doc}`../extending`。
