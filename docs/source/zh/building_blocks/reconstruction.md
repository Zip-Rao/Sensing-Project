# 重建(reconstruction)

## 这层提供什么

`reconstruction` 层是全栈的**反问题层**——把 {doc}`experiments` 层产出的原始测量
数据(布居 $p_e$、相位 $\varphi$、控制核等)**反演**回物理磁通波形 $\Phi(t)$ 或
磁场 $B(t)$。它是波形重建这条产品主线的算法核心:从"比特看到了什么"倒推"外界施加
了什么"。

一条贯穿本层的**纯函数约束**:`reconstruct()` **只读测量数据、绝不跑新仿真**。
若某算法需要正向仿真(如 LM 全密度矩阵反演),所需的 `qubit` / `control_pulse`
在**构造时**注入,而非在 `reconstruct()` 里临时造。输出统一是物理单位的
{py:class}`~sqc.control.FluxSignal`(或 `np.ndarray` 形式的 $B$)。

```{note}
本层含两个标定类 `CryoscopeCalibration` / `DelayRamseyCalibration`。它们在
`sqc.reconstruction` 命名空间下导出(与对应重建类配套),但概念上属**标定**——
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

| 类 | 对应实验 | 反演公式/方法 |
|---|---|---|
| `RamseyReconstruction` | Ramsey | IQ / unwrap 相位 → $B=\dot\varphi/\kappa$ |
| `EchoReconstruction` | 差分回波 | $B=-\varphi/(2k\kappa t_\mathrm{int})$ |
| `TransientReconstruction` | 瞬态场感知 | Wiener / Hammerstein / LM 反卷积 |
| `CryoscopeReconstruction` | Cryoscope | $\dot\varphi\to h(t)$,标定表或色散反演 |
| `DelayRamseyReconstruction` | 延迟 Ramsey | $\varphi(t_d)\to\Phi_\mathrm{tail}$ |
| `PiPulseCompReconstruction` | π-脉冲补偿 | $\Phi_\mathrm{tail}=-z^*(\tau)$ |
| `CryoscopeCalibration` | (标定)| 方波扫 $h$ 建 $\varphi(h)$ 查表 |
| `DelayRamseyCalibration` | (标定)| 扫 $z$ 拟合 $\varphi_\mathrm{cal}(z)=\tau_R\kappa z$ |

## Reconstruction —— 重建抽象基类

所有重建器的公共契约,本层**扩展点**。只规定一个抽象方法:

- `reconstruct(*args, **kwargs)` —— 从测量数据反演物理信号。

两条不变量(继承时必须守):**① `reconstruct()` 必须是纯函数**——只读测量、不跑新
仿真;需要正向仿真的算法(如 LM)把 `qubit`/`control_pulse` 在构造时注入。
**② 输出永远是物理单位的 `FluxSignal`**(或 `np.ndarray` 形式的场)。

要加自定义重建算法,继承 `Reconstruction` 实现 `reconstruct()`,详见 {doc}`../extending`。

## 基函数工具

一组把信号参数化到基函数上的纯函数,供 LM 反演等做低维参数化与正则化用。

- `generate_basis_functions(basis_type, n_basis, t_min, t_max) -> list[callable]`
  —— 在 $[t_\min, t_\max]$ 上生成 `n_basis` 个基函数。`basis_type` 取
  `"bspline"`(3 阶 B-spline)、`"fourier"`(正交正弦余弦)、`"legendre"`
  (Legendre 多项式)。返回一列可调用对象,各接受时间数组、返回基值。
- `basis_function_decomposition(sig, t_array, basis_functions) -> np.ndarray`
  —— 最小二乘把信号 `sig`(数组,或带 `.signal` 属性的 `FluxSignal`)分解到基上,
  返回 `(n_basis,)` 系数。
- `regularization_matrix(n, basis_type) -> np.ndarray` —— 鼓励平滑的 $n\times n$
  正则化矩阵:Fourier 基按频率阶数平方惩罚高频,其余基用二阶差分逼近二阶导。
  `R` 是它的向后兼容别名。

## 控制核估计

### KernelEstimator

估计瞬态反卷积所需的**控制核**:在每个时间点用窄高斯磁通/频率激励扰动控制脉冲、
测比特响应,从而估出线性核 $k_1$ 乃至高阶 Volterra 核。`@dataclass`,关键字段:

- `mode` —— 激励通道:`"flux"`(磁通)或 `"omega"`(频率)。
- `method` —— `"exp"`(实验式:跑仿真测响应)或 `"sim"`。
- `order` —— 最高核阶数(默认 1,只估线性核;$\ge 2$ 估 Volterra 高阶核)。
- `extract_off_diagonal` —— 是否提取非对角(多维)高阶核(仅 LM 反演能消费)。
- 其余:`n_levels`、`anharmonicity`、`kappa`、`stim_amplitude`、`stim_width` 等
  探针参数。

**方法**:`estimate(pulse, qubit=None, t_samples=None) -> (t_samples, k1)` —— 返回
时间轴与线性核;`estimate_full(...) -> KernelResult` —— 返回完整结果对象(含高阶核)。

### KernelResult

核估计结果容器。字段:`t_samples`(核采样时间轴)、`kernels`(各阶核列表)、
`mode`、`method`、`order`、`stim_amplitude`、`units`、`off_diagonal`。属性
`k1`(线性核 = `kernels[0]`);方法 `save(path)`/`load(path)`。可直接传给
`TransientReconstruction.reconstruct(measurement, kernel=...)`。

## RamseyReconstruction —— Ramsey 重建

从 Ramsey 相位反演磁场 $B(\tau)$。`@dataclass`,字段:`qubit`、
`method`(`"iq"` 或 `"unwrap"`,默认 `"unwrap"`)、`k_span`(unwrap 分支搜索窗,默认 3)。

- `method="iq"` —— 双通道 $\arctan2$ 取相位(需 `data["p_e_I"]`/`p_e_Q`)。
- `method="unwrap"` —— 单通道 $\arccos$ + $k$-span 分支解缠绕(需 `data["p_e"]`)。

两条路径都对相位求梯度、除以频率灵敏度 $\kappa=\mathrm{d}\omega/\mathrm{d}\Phi$ 得
$B=\dot\varphi/\kappa$。`reconstruct(measurement) -> np.ndarray`。

## EchoReconstruction —— 差分回波重建

从差分回波 $p_e$ 直接解析反演磁场:$\varphi=\arcsin(2p_e-1)$,
$B=-\varphi/(2k\kappa t_\mathrm{int})$。`@dataclass`,字段:`qubit`、
`t_int`(每个回波块的相互作用时间)、`k`($\pi$ 脉冲对数)。
`reconstruct(measurement) -> np.ndarray`(读 `data["p_e"]`)。

```{note}
`arcsin`/`arccos` 前统一 `clip` 到 $[-1,1]$:`mesolve` 积分噪声可能把 $p_e$ 顶出
$[0,1]$ 微许,不裁会产生静默 NaN。Ramsey/Echo/Transient 各路都有同款保护。
```

## TransientReconstruction —— 瞬态场重建

波形重建主线的**核心反卷积器**,把瞬态实验的 $\Delta p$ 与控制核反卷积回 $\Phi(t)$。
`@dataclass`,关键字段 `method` 选四种算法:

- `"wiener"` —— 线性 Wiener 反卷积(FFT 域),字段 `lambda_reg`(正则化)。
- `"hammerstein"` —— Hammerstein-Wiener 非线性块模型(需 `qubit`)。
- `"hammerstein_volterra"` —— 用高阶核的迭代 Hammerstein-Volterra 反演
  (需 `order>=2` 的核;字段 `max_volterra_iter`、`volterra_tol`)。
- `"lm"` —— Levenberg-Marquardt 全密度矩阵反演(需 `control_pulse`;字段
  `basis_type`、`n_basis`、`max_iter`、`tol`、`mu_init`、`use_adjoint`)。

`reconstruct(measurement, kernel=None, **kwargs)`:`kernel` 可传 `np.ndarray`
(线性核)或 `KernelResult`(含高阶核)。非对角多维核只有 `method="lm"` 能消费,
Wiener/Hammerstein 路径只接受对角(一维)核,否则抛 `ValueError`。

## CryoscopeReconstruction —— Cryoscope 重建

从相位-截断数据重建波形,核心物理 $\mathrm{d}\varphi/\mathrm{d}t=2\pi\Delta f(h(t))$。
`@dataclass`,字段:`tau`(标定方波长)、`inversion`、`calibration`、`qubit`、
`use_sg_filter`/`sg_window`/`sg_poly`(Savitzky-Golay 相位导数预平滑)。两种反演:

- `inversion="calibration"` —— 用 `CryoscopeCalibration` 产的 $\varphi(h)$ 查表反演
  (需 `calibration`)。
- `inversion="response"` —— 用解析 Transmon 频率-磁通色散反演(需 `qubit`)。

`reconstruct(measurement, ...) -> FluxSignal`(读 `data["varphi"]`、`axes["trunc"]`)。

## DelayRamseyReconstruction —— 延迟 Ramsey 重建

把延迟 Ramsey 相位 $\varphi(t_d)$ 反演成拖尾磁通波形(预畸变主线)。`@dataclass`,
字段:`inversion`(`"response"` 解析色散 / `"calibration"` 查表)、`qubit`、
`calibration`、`tau_R`。`reconstruct(measurement, ...) -> FluxSignal`
(读 `data["varphi"]`、`axes["t_d"]`)。

## PiPulseCompReconstruction —— π-脉冲补偿重建

最直接的重建:2D 扫描给出最优补偿高度 $z^*(\tau)$ 后,拖尾即 $\Phi_\mathrm{tail}=-z^*(\tau)$。
`reconstruct(measurement) -> FluxSignal`(读 `data["z_star"]`、`axes["tau"]`)。

```{note}
因 $\pi$ 脉冲有限宽 $T$,$z^*(\tau)$ 测的是拖尾在 $[\tau, \tau+T]$ 上的滑动平均,
而非瞬时值。对指数拖尾这只带来常数衰减因子,衰减时间常数不变。
```

## 标定类(CryoscopeCalibration / DelayRamseyCalibration)

两个在本命名空间导出、但概念属**标定**的类(继承
{py:class}`~sqc.calibration.Calibration`),为对应重建器的 `"calibration"` 反演路径
建查表:

- `CryoscopeCalibration` —— 扫方波高度 $h$、IQ 读出 + 模型引导解缠绕建 $\varphi(h)$
  查表(`calibrate() -> CalibrationTable`,`kind="phi_h"`)。
- `DelayRamseyCalibration` —— 扫已知磁通高度 $z$、拟合线性斜率
  $\varphi_\mathrm{cal}(z)=\tau_R\kappa z$(`calibrate() -> CalibrationTable`,
  `kind="phi_z"`)。

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

- 本层对应实验数据的**离线分析**:把测到的布居/相位翻译回真正施加的磁通/磁场。
- 各重建器与 {doc}`experiments` 的实验类一一配对:Ramsey↔Ramsey、Echo↔差分回波、
  Transient↔瞬态、Cryoscope/延迟 Ramsey/π-补偿↔对应拖尾测量。
- `KernelEstimator` 产的控制核是瞬态反卷积的必需输入;基函数工具支撑 LM 的低维
  参数化;两个标定类为查表反演路径供 `CalibrationTable`。
- **要加自定义重建算法**,继承 `Reconstruction` 抽象基类实现 `reconstruct()`,
  守住"纯函数、输出 `FluxSignal`"两条不变量,需要正向仿真则在构造时注入 `qubit`。
  完整扩展指南见 {doc}`../extending`。
