# 波形重建

## 概述

波形重建是本平台的第一条产品主线：将超导 transmon 量子比特作为时域磁通传感器，
从量子态测量数据中反演作用在比特上、形状未知的外部磁通波形 $\Phi(t)$。

在物理上，外部磁通通过 SQUID 环路调制约瑟夫森能，进而移动比特频率
$\omega_T(\Phi)$。将一个 Ramsey 型控制脉冲沿时间轴**滑动**，在每一延迟
$t_d$ 处测量比特的激发态布居 $p_e$。比特在该延迟附近累积的相干相位正比于磁通
引起的瞬时频率偏移，因此 $p_e(t_d)$ 是 $\Phi(t)$ 经脉冲响应平滑后的近似副本。
这一平滑作用由**控制核** $k(t)$ 描述——它是脉冲对点磁通激励的响应函数。

测量结果与待恢复波形之间构成卷积关系：

$$\Delta p(t) \approx (k * \Phi)(t),$$

其中 $\Delta p$ 是"有信号"与"零磁通参考"两次滑动测量之差，已扣除磁通无关本底。
**重建即为该卷积的反问题**：已知 $\Delta p$ 与控制核 $k$，求解 $\Phi(t)$。

## 管道架构

波形重建管线的执行顺序是：传感器建模 → 信号与探针构造 → 滑动测量 → 实验编排 →
核估计 → 反问题求解 → 全局编排。以下按此逻辑顺序展开各层的职责。

### 第一步：传感器建模 — 设备层 (`TransmonQubit`)

一切重建的物理基础是量子比特本身。{py:class}`~sqc.devices.transmon.TransmonQubit`
持有电荷能 $E_C$、约瑟夫森能 $E_J$、退相干时间 $T_1$/$T_2$ 以及偏置磁通
`flux`。对波形重建最关键的两个方法是：

- `frequency()` — 返回当前偏置下的 $|0\rangle \to |1\rangle$ 跃迁频率
  $\omega_T(\Phi)$。
- `frequency_sensitivity(flux)` — 数值计算 $d\omega/d\Phi$（中心差商），
  即频率对磁通的响应灵敏度 $\kappa$。
- `qubit_in_mag(FluxSignal)` — 预计算 `freq_coeffs` 和 `H_list`，供
  `mesolve()` 使用。**每次磁通信号改变后必须先调用此方法**重建时间依赖
  哈密顿量。

偏置磁通的选择直接影响灵敏度：$\Phi=0$（甜点）处 $\kappa=0$，一阶不敏感；
$\Phi \approx 0.1\,\Phi_0$ 附近 $\kappa$ 最大，为瞬态感知的最优工作点。

### 第二步：构造磁通信号 — 控制层 (`FluxSignal`)

待恢复的未知磁通波形由 {py:class}`~sqc.control.flux_signal.FluxSignal` 表示。
它是一个带物理语义的 `Waveform` 子类，支持 8 种信号类型：零（`type=0`）、
常数（`1`）、正弦（`2`）、高斯（`3`）、非对称冲激（`4`）、双峰（`5`）、
基展开（`6`）、复杂波包（`7`）、自定义原始采样（`8`）。通过 `amplitude`、
`width`、`center`、`rise`、`fall` 等参数控制波形形态，`signal` 属性返回
`(N,)` 采样数组。

仿真中，`TransientSensingExperiment.run()` 会自动构造一个参考零信号（`type=1`,
`amplitude=0`）用于扣除本底。

### 第三步：构造控制脉冲 — 控制层 (`Pulse` + `CompositePulse`)

探针是 Ramsey 型控制脉冲：两个 $\pi/2$ 脉冲，中间自由演化延迟 $\tau=0$（对
瞬态感知，脉冲直接滑过信号，无额外等待）。控制层的三层抽象构建此探针：

- {py:class}`~sqc.control.pulse.Pulse`：单个脉冲。持有 I/Q 包络（`Omega`、
  `Omega_Q`），支持实验坐标系 (`frame=0`) 和旋转坐标系 (`frame=1`)，含或不含
  旋波近似 (RWA)。其 `hamiltonian` 属性返回 QuTiP 列表格式的 $[H_0, c_0(t)]$
  对，可直接传入 `mesolve`。
- {py:class}`~sqc.control.pulse.CompositePulse`：多个 `Pulse` 的串联序列。
  通过拼接各脉冲的时间轴和哈密顿量系数形成全局表示。
- {py:class}`~sqc.control.sequence.create_ramsey_pulse`：工厂函数，便捷构造
  标准 Ramsey 脉冲（两个 $\pi/2$ + 零自由演化延迟）。

### 第四步：执行滑动测量 — 仿真层 (`SlidingMeasurementRunner`)

控制脉冲构造完毕后，{py:class}`~sqc.simulation.runner.SlidingMeasurementRunner`
将其沿信号时间轴逐点滑动，在每个延迟 $t_d$ 处执行一次完整的量子演化：

1. 确定当前延迟下控制脉冲的时间窗口，构造该窗口上的时变哈密顿量。
2. 调用 QuTiP 的 `mesolve()` 求解 Lindblad 主方程，跟踪比特密度矩阵演化。
3. 从末态提取激发态布居 $p_e(t_d)$。
4. 汇总全部延迟点的 $p_e$ 为测量曲线。

这一层是计算量最大的部分——每个延迟点对应一次完整的 `mesolve` 积分。

### 第五步：编排实验 — 实验层 (`TransientSensingExperiment`)

{py:class}`~sqc.experiments.TransientSensingExperiment` 将上述三步串联成
一个完整的测量协议：

1. 调用 `create_ramsey_pulse` 构造控制脉冲。
2. 实例化 `SlidingMeasurementRunner`，对磁通信号执行滑动测量，得到
   $p_e^{\text{sig}}$。
3. 对零磁通参考信号重复，得到 $p_e^{\text{ref}}$，计算
   $\Delta p = p_e^{\text{sig}} - p_e^{\text{ref}}$。
4. 调用 {py:class}`~sqc.reconstruction.kernel.KernelEstimator` 估计控制核
   $k(t)$（详见下一步）。

返回的 {py:class}`~sqc.simulation.result.ExperimentResult` 包含
`data["p_e"]`、`data["delta_p"]`、`data["kernel"]`、`data["flux_samples"]`
（仿真真值）以及轴信息 `axes["scan"]`、`axes["t_flux"]`。

### 第六步：控制核估计与反问题求解 — 重建层

重建层是反问题的算法核心，包含两个紧密协作的组件：核估计器提供前向模型的
平滑核，重建器执行反演。

#### 控制核估计 (`KernelEstimator`)

{py:class}`~sqc.reconstruction.kernel.KernelEstimator` 估计脉冲的点磁通响应
函数 $k(t)$：在每个时间点注入窄高斯磁通激励（`stim_amplitude`、`stim_width`），
测量比特激发态布居的变化，从而逐点构建核。`order=1` 只估计线性核 $k_1$；
`order >= 2` 额外估计高阶 Volterra 核 $k_2, k_3, \ldots$，供非线性反演使用。
返回的 `KernelResult` 对象可直接传入重建器的 `reconstruct()`。

#### 反问题求解 (`TransientReconstruction`)

{py:class}`~sqc.reconstruction.transient.TransientReconstruction` 从
$\Delta p$ 和 $k$ 反演 $\Phi(t)$，提供四种方法，按精度递增、代价递增排列：

| `method` | 算法 | 适用场景 |
|---|---|---|
| `"wiener"` | 频域 Wiener 反卷积 $X_f = \bar{K}_f Y_f / (|K_f|^2 + \lambda^2)$ | 小信号、线性区，速度最快 |
| `"hammerstein"` | Wiener 反卷积 + $\omega(\Phi)$ 色散逆映射 | 中等幅值，$\Phi$ 偏离甜点 |
| `"hammerstein_volterra"` | 迭代 Volterra 级数逆，利用高阶核逐次扣除非线性贡献 | 大信号、非线性明显 |
| `"lm"` | Levenberg-Marquardt 全密度矩阵迭代反演，基函数参数化 + 伴随 Jacobian | 最高精度、代价最大 |

除 LM 外的三种方法输入为 $\Delta p$ 和 $k$，输出 `FluxSignal`。LM 方法额外
需要 `qubit` 和 `control_pulse` 在构造时注入，因每轮迭代需运行正向密度矩阵
仿真计算 $p_e^{\text{sim}}$ 及残差 Jacobian。

与此层配套的还有 {py:class}`~sqc.reconstruction.basis.generate_basis_functions`
等基函数工具，为 LM 提供 B-spline / Fourier / Legendre 基的低维参数化与平滑
正则化矩阵。

```{note}
非对角多维高阶核（`KernelEstimator(extract_off_diagonal=True)` 产出）只有
`method="lm"` 能消费。Wiener / Hammerstein 路径仅接受对角（一维）核，传错会
抛出 `ValueError`。
```

### 第七步：全局时间轴 — 配置层 (`CONFIG`)

{py:class}`sqc.config.CONFIG` 是全局参数的唯一定义源。所有时间轴均通过
`CONFIG.awg.dt`（采样时间步长，$=1/\text{sample\_rate}$）或
`CONFIG.pulse.make_time(start, end)` 派生，保证从脉冲构造到滑动测量到重建的
每一步使用一致的离散化。`SensingWorkflow.configure()` 修改的脉冲/硬件参数
通过 `_sync_config()` 即时同步到 `CONFIG`。

### 第八步：端到端编排 — 工作流层 (`SensingWorkflow`)

{py:class}`~sqc.workflows.SensingWorkflow` 将以上七步封装为一个 `configure()`
$\to$ `run()` $\to$ `plot()` 接口，是面向用户的单一入口。它本身不实现物理，
而是负责：

- **配置管理**：按功能组（qubit、signal、reconstruction、pulse、hardware）
  存储参数，`configure()` 只修改显式传入的字段。
- **管道执行**：`run()` 依次调用 `_build_qubit()` → `_build_flux_signal()` →
  `_make_experiment().run()` → `_run_reconstruction()`。
- **批处理**：`sweep("signal.amplitude", [...])` 沿单参数扫描，自动管理
  CONFIG 同步；`compare(methods=[...])` 在同一份 $\Delta p$ 上对比多种重建
  算法，按 RMSE 选优，不重跑 `mesolve`。
- **可视化**：`plot()` 依据最近执行的操作自动选择图类型（$p_e$/$\Delta p$
  测量、重建波形 vs 真值、扫参曲线、多方法 overlay + 残差）。

## 使用方式

### 端到端管道（推荐）

```python
from sqc.workflows import SensingWorkflow

# 配置
wf = SensingWorkflow()
wf.configure(
    protocol="transient",       # 瞬态场感知（滑动 Ramsey）
    flux_bias=0.1,              # 偏置至 κ 有限处
    signal_type=4,              # 待恢复波形：非对称冲激
    signal_amplitude=0.02,      # 峰值磁通 (Φ₀)
    signal_center=100,          # 冲激中心 (ns)
    reconstruction="wiener",    # 重建算法
    lambda_reg=5.0,             # Wiener 正则化参数
)

# 执行：测量 + 重建
result = wf.run()
print("恢复波形点数:", len(result.reconstructed_signal.signal))

# 同一份测量上比较三种算法
cmp = wf.compare(methods=["wiener", "hammerstein", "lm"])
print("最优方法:", cmp.best)
for m in cmp.methods:
    print(f"  {m}: rmse={cmp.metrics[m]['rmse']:.3e}")

# 出图
wf.plot()
```

```{note}
`compare()` 复用 `run()` 缓存的测量数据，仅切换重建算法，不重跑 `mesolve`，
三种方法在同一批 $\Delta p$ 上比较。`best` 按 RMSE 选优；仿真中真值波形
`flux_samples` 已知，RMSE 为有效指标。
```

### 分步执行管道

若需对单层精细控制，可按管道的逻辑顺序直接操作各层对象：

```python
import numpy as np
from sqc.devices import TransmonQubit
from sqc.control.flux_signal import FluxSignal
from sqc.experiments.transient import TransientSensingExperiment
from sqc.reconstruction.transient import TransientReconstruction

# 1. 构造比特（传感器）
qubit = TransmonQubit(
    EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
    T1=10_000, T2=8_000, flux=0.1,
)

# 2. 构造磁通信号（被测对象）
t_list = np.arange(0, 400, 0.5)  # dt = 0.5 ns
flux_signal = FluxSignal(
    type=4, t_list=t_list,
    amplitude=0.02, center=100, rise=10, fall=10,
)

# 3–5. 运行实验（脉冲构造 + 滑动测量 + Δp 与核）
exp = TransientSensingExperiment(qubit=qubit, flux_signal=flux_signal)
meas = exp.run()

# 6. 反问题求解
rec = TransientReconstruction(method="wiener", lambda_reg=5.0)
recovered = rec.reconstruct(meas, kernel=meas.data["kernel"])
print(recovered.signal[:5])
```

## 结果解读

- `result.measurement.data` 包含四个关键数组：`p_e`（激发态布居）、`delta_p`
  （扣除本底后的信号差异）、`kernel`（控制核 $k$）、`flux_samples`（仿真注入
  的真值波形，仅仿真可得）。轴信息在 `result.measurement.axes` 中：`scan`
  （滑动延迟轴）、`t_flux`（波形时间轴）。
- `result.reconstructed_signal` 是恢复出的
  {py:class}`~sqc.control.flux_signal.FluxSignal`。`wf.plot()` 将其叠加在真
  值上，直观比较重建质量。
- `compare()` 的 `metrics` 为每种方法提供 `rmse`、`snr` 和 `peak`。典型结论：
  Wiener 最快、对小信号足够；Hammerstein 补偿了 $\omega(\Phi)$ 的弱非线性，
  大信号下 RMSE 更低；LM 最准确但代价最大（每轮迭代均需正向密度矩阵仿真）。
- 三个重建算法的物理差异详见 {doc}`../building_blocks/reconstruction` §
  TransientReconstruction。
- 需扫参数时使用 `wf.sweep("signal.amplitude", [...])` 获得 SNR-幅值曲线。
  完整 API 见 {doc}`../building_blocks/workflows`。
