# 工作流(workflows)

## 这层提供什么

`workflows` 层是全栈的**顶层编排层**——把 {doc}`experiments`、{doc}`reconstruction`、
{doc}`calibration` 三层组合成完整的研究流程,对用户暴露一个高层 API。它不引入新物理,
只负责"配置 → 测量 → 重建 → 对比 → 可视化"的串联与批处理。

`SensingWorkflow` 是平台**面向用户的单一入口**:一个类里 `configure()` 配参、`run()`
跑管道、`sweep()` 扫参、`compare()` 比较重建算法、`plot()` 自动出图。
`PredistortionValidationWorkflow` 则是预畸变主线的端到端验证工作流。

## 类总览

**扩展点**

| 类 | 角色 |
|---|---|
| `Workflow` | 抽象基类,所有工作流的公共契约(`run() -> dict`),本层扩展点 |

**工作流**

| 类 | 角色 |
|---|---|
| `SensingWorkflow` | 统一研究入口:配置/执行/扫参/对比/可视化 |
| `PredistortionValidationWorkflow` | 预畸变端到端验证:注入畸变→标定→设计逆滤波→验证改善 |

**结果容器(活跃)**

| 类 | 角色 |
|---|---|
| `WorkflowResult` | 单次 `run()` 结果(测量、重建信号、标定表、配置快照)|
| `SweepResult` | `sweep()` 扫参结果(每个取值一份 `WorkflowResult` + 指标)|
| `CompareResult` | `compare()` 算法对比结果(各方法信号 + 指标 + 最优者)|

**结果容器(未来功能占位)**

| 类 | 对应(未实现)方法 |
|---|---|
| `DiffReport` / `BenchmarkResult` / `NoiseReport` / `CVResult` | 差异报告 / 基准 / 噪声表征 / 交叉验证(post-v1)|

## Workflow —— 工作流抽象基类

所有工作流的公共契约,本层**扩展点**。只规定一个抽象方法:

- `run() -> dict` —— 执行工作流,返回命名结果字典。

要加自定义工作流,继承 `Workflow` 实现 `run()`,详见 {doc}`../extending`。

## SensingWorkflow —— 统一研究入口

平台面向用户的单一 API,一个类串起配置、协议执行、扫参、A/B 对比与可视化,内部委托
experiment / reconstruction / calibration 三层。核心方法:

- `configure(**params) -> self` —— 按功能分组设参(只改显式传入的,其余保持),支持链式。
  分组:**Qubit**(`EC`/`EJ`/`T1`/`T2`/`flux_bias`/`n_levels`,`EC`/`EJ` 单位 GHz)、
  **Protocol & Signal**(`protocol`、`signal_*`)、**Reconstruction**(`reconstruction`
  方法名及 `lambda_reg`/`lm_*`)、**Pulse**(`t_rabi_duration`/`t_global_*`,单位 ns)、
  **Hardware**(`sample_rate`,GSa/s)。
- `run(measure=True, reconstruct=True, calibrate=False) -> WorkflowResult` —— 跑一条
  完整感知管道。`calibrate=True` 仅 `cryoscope`/`delay_ramsey` 协议支持(测量前先标定)。
- `sweep(param_path, values) -> SweepResult` —— 沿点分路径 `"组.字段"`(如
  `"signal.amplitude"`)扫一串取值,每点跑一次并算 SNR/RMSE/peak 指标。
- `compare(methods, ...) -> CompareResult` —— 在同一份测量上比较多种重建算法,按 RMSE
  选最优。
- `plot()` —— 依据最近一次 `run`/`sweep`/`compare` 自动派发可视化。

协议 → (实验, 重建) 映射:`"transient"`、`"ramsey"`、`"echo"`、`"cryoscope"`、
`"delay_ramsey"` 各对应一对实验/重建类。

```{note}
若干进阶方法(如批量基准、噪声表征、交叉验证、多比特、持久化等)为 post-v1 占位,
当前抛 `NotImplementedError`,不进 v1 公开面。v1 用 `run`/`sweep`/`compare` 已覆盖
单管道、扫参与算法对比。
```

## PredistortionValidationWorkflow —— 预畸变端到端验证

预畸变主线的完整验证工作流(`@dataclass`)。给定目标波形与已知畸变,跑:注入畸变到
控制线 → 测传函 → 拟合模型 → 设计逆滤波器 → 施加预畸变 → 复测 → 算改善指标。字段:
`target_waveform`(目标片上波形)、`true_distortion`(真值畸变模型)、`designer`
(逆滤波器设计器,默认 `PredistortionDesigner(method="auto")`)、`control_line_params`;
以及协议驱动测量的 `qubit`、`measurement_protocol`(为 `None` 时走解析路径)。

`run() -> dict` 返回 `target`、`on_chip_uncorrected`、`on_chip_corrected`、
`awg_predistorted`、`inverse_model`、`measured_model`、`metrics`。`metrics` 含
`rmse_uncorrected`、`rmse_corrected`、`improvement_factor`(改善倍数)与两个整定时间。

## 结果容器

- `WorkflowResult` —— 单次 `run()` 结果:`config_snapshot`、`measurement`
  (`ExperimentResult`)、`reconstructed_signal`(`FluxSignal`)、
  `reconstruction_details`、`calibration`(`CalibrationTable`)。
- `SweepResult` —— `sweep()` 结果:`param_path`、`values`、`results`
  (每点一份 `WorkflowResult`)、`metrics`(键 → 逐点指标列表)。
- `CompareResult` —— `compare()` 结果:`methods`、`signals`(方法→`FluxSignal`)、
  `metrics`(方法→指标字典)、`best`(最低 RMSE 者)。
- `DiffReport` / `BenchmarkResult` / `NoiseReport` / `CVResult` —— 分别是差异报告、
  批量基准、噪声表征、交叉验证的返回类型;对应方法为 post-v1 占位,当前仅作接口预留。

## 最小用例

```python
import numpy as np
from sqc.workflows import SensingWorkflow
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.workflows import PredistortionValidationWorkflow

# ── 1. SensingWorkflow：配置 + 跑单次 ─────────────────────────────────
wf = SensingWorkflow()
wf.configure(
    protocol="transient",
    signal_amplitude=0.01,
    reconstruction="wiener",
    flux_bias=0.1,
)
result = wf.run()          # WorkflowResult
print(result.reconstructed_signal)   # FluxSignal

# ── 2. 扫参(信号幅度) ─────────────────────────────────────────────
sweep = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02, 0.05])
print(sweep.metrics["snr"])          # 每个幅度对应的 SNR

# ── 3. A/B 比较重建算法 ──────────────────────────────────────────────
wf.run(measure=True, reconstruct=False)   # 只测量、缓存数据
cmp = wf.compare(methods=["wiener", "hammerstein", "lm"])
print(f"最优方法: {cmp.best}")

# ── 4. PredistortionValidationWorkflow ────────────────────────────────
t = np.linspace(0, 500, 1000)
target = Waveform(t_list=t, samples=np.ones_like(t) * 0.05)
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)
val_wf = PredistortionValidationWorkflow(
    target_waveform=target,
    true_distortion=distortion,
)
val_result = val_wf.run()
print(f"改善倍数: {val_result['metrics']['improvement_factor']:.1f}×")
```

## 物理角色 / 扩展

- 本层对应实验室的**上层实验脚本**:选协议、调参、跑一批数据、出结果。它把
  {doc}`devices`→{doc}`hardware`→{doc}`control`→{doc}`simulation`→{doc}`experiments`
  →{doc}`reconstruction`→{doc}`calibration` 七层全部串联,是用户无需直接操作底层的便捷入口。
- `SensingWorkflow.sweep()` 自动管理 CONFIG 同步(脉冲/硬件组参数改变时调
  `_sync_config()`),保证每次 `run()` 用的是当前参数对应的时间轴。
- `compare()` 把测量缓存复用给多个重建算法,避免重复跑 `mesolve`;`best` 字段按 RMSE
  自动选优,有真值(仿真 `flux_samples`)时有意义,无真值时是对 NaN 最稳健的方法。
- **要加自定义工作流**,继承 `Workflow` 抽象基类实现 `run() -> dict`,在内部按需编排
  实验/重建/标定三层。完整扩展指南见 {doc}`../extending`。
