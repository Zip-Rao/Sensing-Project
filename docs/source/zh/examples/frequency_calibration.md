# 频率标定

## 概述

频率标定是本平台的第二条产品主线：表征超导 transmon 量子比特的频率-磁通响应
$f_{01}(\Phi)$，并据此回答两个实用问题——**在给定磁通偏置下 qubit 频率是多少**，
以及**要把频率调到目标值需要施加多少磁通偏置**。与波形重建不同，频率标定不反演
外部信号，而是**表征器件本身**：投片后确定工作点、建立查表、必要时闭环整定到目标
频率。

在物理上，transmon 频率由外部磁通经 SQUID 环路调制约瑟夫森能而定：

$$f_{01}(\Phi) \approx \frac{1}{2\pi}\left(\sqrt{8 E_J(\Phi)\, E_C} - E_C\right),
\qquad E_J(\Phi) = E_{J0}\,|\cos(\pi\Phi/\Phi_0)|.$$

由于 $E_J \propto |\cos(\pi\Phi/\Phi_0)|$，$f_{01}(\Phi)$ 是关于 $\Phi=0$ 的
**偶函数**，在整数磁通量子处取极大（**甜点**，$\mathrm{d}f/\mathrm{d}\Phi=0$，
对磁通噪声一阶不敏感）。偏置到甜点一侧则获得传感灵敏度
（$\kappa = \mathrm{d}\omega/\mathrm{d}\Phi$ 非零），这正是
{doc}`waveform_reconstruction` 中 `flux_bias` 的来源。

频率标定管线的执行顺序是：单点测频 → 扫磁通建 $f(\Phi)$ 曲线 → 查表 →
闭环整定 → 多阶段编排。以下按此逻辑顺序展开各层的职责。

## 管道架构

### 第一步：单点测频 — 标定层 (`FrequencyMeasurement`)

一切频率标定的基本操作是在**一个磁通工作点**上精确测量 $f_{01}$。
{py:class}`~sqc.calibration.frequency.FrequencyMeasurement` 是只读、不整定的
单点频率计，内部驱动 Ramsey 序列（或瞬态正交 Ramsey），跑 `mesolve` 后从数据
中提取频率。

核心方法是 `measure(flux=None, omega_d=None) -> float`：在指定磁通偏置 `flux`
处运行测量，返回有符号角频率（rad·GHz）。`omega_d` 参数指定参考驱动频率——
测得失谐 $\hat\delta = \hat f_q - f_d$ 后加回即得绝对频率。这在后续闭环搜索中
至关重要：将前次估计作为 `omega_d` 传入，保证当次测量的失谐落在鉴频器的线性窗内。

`method` 选两种测量协议：

| `method` | 原理 | 代价 | 适用 |
|---|---|---|---|
| `"ramsey"` | $\tau$ 扫描 + FFT 取峰 | 数十次 `mesolve` | 稳健，通用 |
| `"transient"` | $\tau=0$ 正交 Ramsey + 核灵敏度 $G=\int k_1\,dt$ | 2 次 `mesolve` | 快速，适合 $\Delta\omega \approx 0$ |

Ramsey 模式默认单扫（`f_artificial=0.1` GHz），假设 $|\Delta| < 0.1$ GHz；
设 `f_artificial=None` 则走双扫，对任意失谐稳健且返回符号，代价 2×。

### 第二步：扫磁通建 $f(\Phi)$ 曲线 — 标定层 (`FluxResponseCalibration`)

有了单点测频能力，下一步是沿磁通轴扫描，逐点调用 `FrequencyMeasurement`，
构建频率-磁通查表。{py:class}`~sqc.calibration.frequency.FluxResponseCalibration`
负责这一过程。

`@dataclass` 字段：`qubit`、`h_list`（磁通扫描点，$\Phi_0$，默认 51 点）、
`method`（当前仅 `"ramsey"`）。`calibrate()` 返回 `kind="f_phi"` 的
{py:class}`~sqc.calibration.CalibrationTable`：`inputs` 为磁通、`outputs` 为
角频率。此表有两个下游用途：为闭环整定提供磁通上下界 `[V_a, V_b]`；直接正向/
反向查表。

### 第三步：查表与闭环整定 — 标定层 (`CalibrationTable` + `SinglePointFrequencyCalibration`)

#### 正向/反向查表

{py:class}`~sqc.calibration.CalibrationTable` 提供两个基于三次样条的查表方法：

- `evaluate(x)` — **正向**：磁通 → 频率（如"在 $\Phi=0.015$ 处 $f_{01}$ 是多少"）。
- `inverse(y)` — **反向**：频率 → 磁通（如"要得到目标频率该偏置多少磁通"）。

由于 $f(\Phi)$ 是偶函数、整体不单调，`inverse` 自动在单调段上构造反函数；
查询目标频率必须落在甜点一侧，否则解不唯一。

#### 闭环整定

查表给出的是**开环估计**——用插值反推所需磁通，但未考虑实际测量噪声与模型偏差。
{py:class}`~sqc.calibration.frequency.SinglePointFrequencyCalibration` 将
$f_q(V)$ 闭环反馈整定到目标频率 $f_\text{target}$（Vepsalainen 2022），通过
迭代测量-调整逼近真值。

核心字段：`f_target`（目标角频率）、`V_a`/`V_b`（磁通框界，来自前序
`FluxResponseCalibration` 确定的单调支）、`step_method`（根搜索方法）、
`measure_method`（每轮测频协议，内部委托给 `FrequencyMeasurement`）。

三种步进方法：

| `step_method` | 收敛速度 | 需要框界 | 特点 |
|---|---|---|---|
| `"secant"` | 超线性，1–3 轮 | 需要 | 配 `bracket_tightening`（regula falsi）自动缩框 |
| `"bisection"` | $O(\log_2)$，10–15 轮 | 需要 | 框宽每步折半，自动处理偶函数过甜点 |
| `"gradient"` | 阻尼 Newton | 不需要 | 只需 `V_seed`；`damping`（默认 0.8）抑超调；`best_V` 跟踪最优 |

`calibrate()` 返回 `kind="f01"` 的表，`fit_params["history"]` 含完整迭代轨迹
（每轮的 $V$、$f$、残差），可用于画收敛曲线。

```{note}
闭环每轮迭代同时执行两种独立更新：**磁通电压更新**（根搜索，将 $f_q$ 推向目标）
和**驱动频率更新**（观测器，`drive_policy` 设定 $f_d$ 使测量落在线性窗内）。
驱动频率不进入误差定义——$e_k = \text{measure}(V_k) - f_\text{target}$ 始终
相对固定目标——因此只影响测量可信度，不改变收敛目标。三种驱动策略
（`"sweet"`/`"target"`/`"track"`）配合不同测频协议构成四阶段状态机
（COARSE_ACQUIRE → TRACKING → LOCKED → REACQUIRE），由下一步的
`FrequencyCalibrationWorkflow` 编排。各策略的完整说明见
{doc}`../building_blocks/calibration`。
```

### 第四步：多阶段编排 — 工作流层 (`FrequencyCalibrationWorkflow`)

单次 `SinglePointFrequencyCalibration` 的 `measure_method` 在构造时固定，
无法在搜索过程中切换——粗阶段需要宽范围的 Ramsey，细阶段受益于快速瞬态法。
{py:class}`~sqc.workflows.frequency_calibration.FrequencyCalibrationWorkflow`
将多个标定阶段串联为有序管线，每阶段以前一阶段的最优磁通和频率估计为初值。

两种构造方式：

- **默认 hybrid 预设**（不传 `stages`）：自动构建两阶段瞬态→Ramsey 管线。
  `switch_residual`（默认 $2\pi \cdot 5$ MHz）控制粗-细交棒阈值。
- **显式管线**（传 `stages=[...]`）：每阶段为一个 `CalibrationStage`，独立
  指定 `measure_method`、`step_method`、`epsilon_f`、`drive_policy` 等。

`run()` 返回合并的迭代历史（各行标 `phase`/`global_iter`/`cost`）、
`V_final`、`residual`、`converged` 等。

### 底层支撑：设备层 (`TransmonQubit`)

所有标定操作的物理基础是 {py:class}`~sqc.devices.transmon.TransmonQubit`——
持有 $E_C$、$E_J$、退相干时间及当前偏置磁通。标定类内部通过
`qubit.qubit_in_mag(FluxSignal)` 设置直流磁通偏置、更新哈密顿量后运行
`mesolve`。`qubit.frequency` 提供甜点频率作为 Ramsey 测量的默认驱动参考。

## 使用方式

### 端到端管道

```python
import numpy as np
from sqc.devices.transmon import TransmonQubit
from sqc.calibration import FluxResponseCalibration, FrequencyMeasurement
from sqc.calibration import SinglePointFrequencyCalibration

# ── 器件：EC/EJ 用角频率 (rad·GHz) 传入 ──────────────────────────────
qubit = TransmonQubit(
    EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
    T1=10_000, T2=8_000, flux=0.0, n_levels=3,
)

# ── 1. 单点测频：甜点处的 f01 ───────────────────────────────────────
fm = FrequencyMeasurement(qubit=qubit, flux=0.0, method="ramsey")
f01 = fm.measure()                       # 有符号角频率 (rad·GHz)
print(f"甜点 f01 = {f01 / (2*np.pi):.4f} GHz")

# ── 2. 扫磁通建 f(Φ) 查表 ──────────────────────────────────────────
cal = FluxResponseCalibration(
    qubit=qubit,
    method="ramsey",
    h_list=np.linspace(-0.03, 0.03, 5),  # 演示用粗网格
)
table = cal.calibrate()                  # CalibrationTable, kind="f_phi"

# ── 3a. 正向查表：Φ → f ────────────────────────────────────────────
f_at_bias = table.evaluate(np.array([0.015]))

# ── 3b. 反向查表 + 闭环整定 ────────────────────────────────────────
f_target = table.outputs.max() * 0.999   # 略低于甜点（留在单调支）
tuner = SinglePointFrequencyCalibration(
    qubit=qubit,
    f_target=f_target,
    V_a=0.0, V_b=0.03,                   # 框界，取自 f(Φ) 单调支
    step_method="secant",                 # 割线法，典型 1–3 轮收敛
)
result = tuner.calibrate()               # CalibrationTable, kind="f01"
print("整定后偏置:", result.fit_params["V_opt"],
      "收敛:", result.fit_params["converged"])
```

### 多阶段 hybrid 管线

```python
from sqc.workflows.frequency_calibration import FrequencyCalibrationWorkflow

# 默认两阶段 hybrid：瞬态粗搜（三次修正）→ Ramsey 精调
wf = FrequencyCalibrationWorkflow(
    qubit=qubit,
    f_target=f_target,
    V_a=0.0, V_b=0.03,
    switch_residual=2 * np.pi * 5e-3,    # 5 MHz 粗-细交棒
    epsilon_f=1e-4,                       # 最终收敛容差
)
hybrid_result = wf.run()
print(f"V_final={hybrid_result['V_final']:.6f}, "
      f"residual={hybrid_result['residual']/(2*np.pi)*1e3:.2f} MHz, "
      f"converged={hybrid_result['converged']}")
```

```{note}
`FluxResponseCalibration` 在每个磁通点跑一次完整 Ramsey $\tau$ 扫描（数十次
`mesolve`），默认 51 点，完整标定是分钟级任务。5 点粗网格仅为演示流程；实际
运行请按精度需求加密，或只在关心的频段附近细扫。
```

## 结果解读

- `fm.measure()` 返回单个工作点的有符号角频率（rad·GHz），除以 $2\pi$ 得 GHz。
  默认单扫假设 $|\Delta| < 0.1$ GHz；若可能远离甜点，设 `f_artificial=None`
  走双扫模式。
- `table` 是 `kind="f_phi"` 的 `CalibrationTable`：`inputs` 为磁通点，`outputs`
  为角频率。`evaluate` 做三次样条插值（磁通→频率）；`inverse` 做反插值（频率→
  磁通），自动取单调段。典型曲线在 $\Phi=0$ 取极大，两侧对称下降——甜点提供
  一阶抗磁通噪声，偏置到旁边则获得传感灵敏度。
- 闭环整定的 `result.fit_params` 含 `V_opt`（最优磁通）、`converged`（是否容差
  内收敛）、`history`（逐轮 $V$/$f$/residual）。`history` 可直接用于画收敛曲线。
- 多阶段 hybrid 的 `run()` 返回合并的 `history`（各行标 `phase`），以及
  `stage_boundaries`（各阶段结束时的全局迭代序号）和 `metrics`（累计 `mesolve`
  调用代价）。
- `FrequencyMeasurement` 支持 `order=3` 三次 Newton 修正，通过 `g3_source`
  选择三次系数来源：`"fit"`（奇多项式拟合 $p_\text{diff}(\Delta)$，自适应范围）
  或 `"kernel_full"`（非对角核 $\iiint k_3\,dt^3$）。
- 各标定类的完整字段与方法选项见 {doc}`../building_blocks/calibration`。
