# 工作流(workflows)

## 这层提供什么

`workflows` 层是全栈的顶层编排层,把 {doc}`experiments`、{doc}`reconstruction`、
{doc}`calibration` 三层组合成完整的研究流程,对用户暴露一个高层 API。它不引入新物理,
只负责配置、测量、重建、对比、可视化各步骤的串联与批处理。

`SensingWorkflow` 是平台面向用户的单一入口:一个类里 `configure()` 配参、`run()`
跑管道、`sweep()` 扫参、`compare()` 比较重建算法、`plot()` 自动出图。
`PredistortionValidationWorkflow` 是预畸变主线的端到端验证工作流。

### 频率标定的两种用法

频率标定有三个相关类，按**新旧两套接口**组织：

| | 旧接口（一次性） | 新接口（事件驱动） |
|---|---|---|
| **入口** | `FrequencyCalibrationWorkflow` | `FrequencyCalibrationRuntime` |
| **核心** | 内部串联多个 `SinglePointFrequencyCalibration` | 内部驱动 `FrequencyStateMachine` |
| **怎么跑** | `wf.run()` 一口气跑完 | `runtime.run()` 逐步事件循环 |
| **能暂停吗** | 不能 | 能（`save_run` / `load_run` 持久化） |
| **能回退吗** | 不能（阶段间只能前进） | 能（Track→Reacquire→Track, Lock→Verify→Lock） |
| **共享组件** | `DampedSecantTracker`（步进公式） | 同一个 `DampedSecantTracker` |

```text
用户
 │
 ├── (旧) FrequencyCalibrationWorkflow.run()
 │       内部: stage1 → stage2 → ...  一次性跑完
 │       每 stage 内委托 DampedSecantTracker
 │
 └── (新) FrequencyCalibrationRuntime.run()
             │
             ├── FrequencyStateMachine  ← 纯状态转移（无 QuTiP）
             │     六个状态: Acquire → Track → Verify → Lock
             │     回退分支: Reacquire, SafeStop
             │     DampedSecantTracker ← 共享步进公式
             │
             └── SQCExecutor ← 命令 → QuTiP 测量
```

**选型指南**：

- 想一键跑完、不需要中途干预 → 用 `FrequencyCalibrationWorkflow`
- 想逐步控制、可回退、可持久化、可模拟长期运行 → 用 `FrequencyCalibrationRuntime` + `FrequencyStateMachine`

两者底层共享同一个 `DampedSecantTracker` 控制律，数值行为一致。

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
| `FrequencyCalibrationWorkflow` | 多阶段闭环频率标定 (legacy staged workflow) |
| `FrequencyStateMachine` | 事件驱动频率标定状态机 V2: Acquire→Track→Verify→Lock |
| `FrequencyCalibrationRuntime` | 状态机运行时:事件循环 + tracker + 持久化 |

**结果容器(活跃)**

| 类 | 角色 |
|---|---|
| `WorkflowResult` | 单次 `run()` 结果(测量、重建信号、标定表、配置快照) |
| `SweepResult` | `sweep()` 扫参结果(每个取值一份 `WorkflowResult` + 指标) |
| `CompareResult` | `compare()` 算法对比结果(各方法信号 + 指标 + 最优者) |

**结果容器(未来功能占位)**

| 类 | 对应(未实现)方法 |
|---|---|
| `DiffReport` / `BenchmarkResult` / `NoiseReport` / `CVResult` | 差异报告 / 基准 / 噪声表征 / 交叉验证(post-v1) |

## Workflow —— 工作流抽象基类

所有工作流的公共契约,本层扩展点。只规定一个抽象方法:

- `run() -> dict`:执行工作流,返回命名结果字典。

要添加自定义工作流,继承 `Workflow` 实现 `run()`,详见 {doc}`../extending`。

## SensingWorkflow —— 统一研究入口

平台面向用户的单一 API,一个类串起配置、协议执行、扫参、A/B 对比与可视化,内部委托
experiment / reconstruction / calibration 三层。

**构造**

`SensingWorkflow()`

**方法**

- `configure(**params) -> self`：按功能分组设参(只改显式传入的,其余保持),支持链式。
  分组:Qubit(`EC`/`EJ`/`T1`/`T2`/`flux_bias`/`n_levels`,`EC`/`EJ` 单位 GHz)、
  Protocol & Signal(`protocol`、`signal_*`)、Reconstruction(`reconstruction`
  方法名及 `lambda_reg`/`lm_*`)、Pulse(`t_rabi_duration`/`t_global_*`、
  `rotation_angle`/`rabi_rate`、`pulse_envelope`/`envelope_sigma`及两段相位)、
  Hardware(`sample_rate`,GSa/s)。
- `run(measure=True, reconstruct=True, calibrate=False) -> WorkflowResult`：跑一条
  完整感知管道。`calibrate=True` 仅 `cryoscope`/`delay_ramsey` 协议支持(测量前先标定)。
- `sweep(param_path, values) -> SweepResult`：沿点分路径 `"组.字段"`(如
  `"signal.amplitude"`)扫一串取值,每点跑一次并算 SNR/RMSE/peak 指标。
- `compare(methods, ...) -> CompareResult`：在同一份测量上比较多种重建算法,按 RMSE
  选最优。
- `plot()`：依据最近一次 `run`/`sweep`/`compare` 自动派发可视化。

协议 → (实验, 重建) 映射:`"transient"`、`"ramsey"`、`"echo"`、`"cryoscope"`、
`"delay_ramsey"` 各对应一对实验/重建类。

**输出**

- `run()` → `WorkflowResult`(含 `measurement`、`reconstructed_signal`、`config_snapshot`)
- `sweep()` → `SweepResult`(含 `results`、`metrics`)
- `compare()` → `CompareResult`(含 `signals`、`metrics`、`best`)

```{note}
若干进阶方法(如批量基准、噪声表征、交叉验证、多比特、持久化等)为 post-v1 占位,
当前抛 `NotImplementedError`,不进 v1 公开面。v1 用 `run`/`sweep`/`compare` 已覆盖
单管道、扫参与算法对比。
```

## FrequencyCalibrationWorkflow —— 多阶段频率标定 (legacy)

把 `SinglePointFrequencyCalibration` 编排为多阶段管线（如 transient→Ramsey hybrid）,每阶段以前一阶段的最优磁通和频率估计作初值。内部已委托给 `DampedSecantTracker`。

**构造**

`FrequencyCalibrationWorkflow(qubit, f_target, stages=None, V_seed=None, ...)`

不传 `stages` 时自动构建默认两阶段 hybrid: coarse(transient+gradient) → fine(ramsey)。详见 {doc}`calibration`。

## FrequencyStateMachine —— 事件驱动频率标定 V2

v2.19 新增的**事件驱动**六状态协议: Acquire → Track → Verify → Lock,带 Reacquire 恢复分支和 SafeStop 安全保持。

**设计原则**:
- 纯状态机,无 QuTiP 依赖——可脱离仿真环境做转移表测试
- 只有 **Track** 可以提出新的工作磁通偏置
- **Verify** 从进入冻结候选偏置到退出,Lock 同样不改偏置
- **Lock** 监测到漂移不直接调偏置——小漂移去 Verify,大跳变去 Reacquire
- 解析 $f(\Phi)$ 只做仿真 oracle,不进转移决策

**使用方式**（通过 `FrequencyCalibrationRuntime` 编排）:

```python
from sqc.workflows.frequency_runtime import FrequencyCalibrationRuntime
from sqc.workflows.frequency_state_machine import FrequencyCalibrationConfig

config = FrequencyCalibrationConfig(
    epsilon_enter=2*np.pi*20e-3, epsilon_final=2*np.pi*2e-3,
    N_verify=2, max_commands=30,
)
runtime = FrequencyCalibrationRuntime(qubit=q, f_target=f_target, config=config)
result = runtime.run()
# result["state"] → "lock", result["run_status"] → "calibrated"
```

详见 {doc}`calibration` 中「频率标定状态机 V2」一节。

## PredistortionValidationWorkflow —— 预畸变端到端验证

预畸变主线的完整验证工作流。给定目标波形与已知畸变,跑:注入畸变到
控制线 → 测传函 → 拟合模型 → 设计逆滤波器 → 施加预畸变 → 复测 → 算改善指标。

**构造**

`PredistortionValidationWorkflow(target_waveform, true_distortion, designer=None, control_line_params=None, qubit=None, measurement_protocol=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `target_waveform` | `Waveform` | 目标片上波形 | —— |
| `true_distortion` | `DistortionModel` | 真值畸变模型 | —— |
| `designer` | `PredistortionDesigner` | 逆滤波器设计器 | `PredistortionDesigner(method="auto")` |
| `control_line_params` | dict | 控制线构造参数 | `None` |
| `qubit` | `TransmonQubit` | 比特(协议驱动测量用) | `None` |
| `measurement_protocol` | str | 测量协议(`None`=解析路径) | `None` |

**方法**

- `run() -> dict`：执行完整验证管道。

**输出**

返回 `dict`,包含:

| 键 | 类型 | 含义 |
|---|---|---|
| `target` | `Waveform` | 目标波形 |
| `on_chip_uncorrected` | `Waveform` | 未校正的片上波形 |
| `on_chip_corrected` | `Waveform` | 校正后的片上波形 |
| `awg_predistorted` | `Waveform` | 预畸变后的 AWG 波形 |
| `inverse_model` | `DistortionModel` | 逆滤波器模型 |
| `measured_model` | `DistortionModel` | 测量拟合的传函模型 |
| `metrics` | dict | `rmse_uncorrected`、`rmse_corrected`、`improvement_factor`、整定时间 |

## 结果容器

### WorkflowResult

单次 `run()` 结果。

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `config_snapshot` | dict | 配置快照 |
| `measurement` | `ExperimentResult` | 原始测量数据 |
| `reconstructed_signal` | `FluxSignal` | 重建的磁通信号 |
| `reconstruction_details` | dict | 重建算法细节 |
| `calibration` | `CalibrationTable` | 标定表(可选) |

### SweepResult

`sweep()` 结果。

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `param_path` | str | 扫描参数路径 |
| `values` | list | 扫描取值列表 |
| `results` | `list[WorkflowResult]` | 每点一份结果 |
| `metrics` | dict | 键 → 逐点指标列表(如 `snr`、`rmse`) |

### CompareResult

`compare()` 结果。

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `methods` | `list[str]` | 比较的方法名列表 |
| `signals` | dict | 方法名 → `FluxSignal` |
| `metrics` | dict | 方法名 → 指标字典 |
| `best` | str | 最低 RMSE 的方法名 |

### 占位容器

`DiffReport` / `BenchmarkResult` / `NoiseReport` / `CVResult`:分别是差异报告、
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

- 本层对应实验室的上层实验脚本:选协议、调参、跑一批数据、出结果。它把
  {doc}`devices`→{doc}`hardware`→{doc}`control`→{doc}`simulation`→{doc}`experiments`
  →{doc}`reconstruction`→{doc}`calibration` 七层全部串联,是用户无需直接操作底层的便捷入口。
- `FrequencyCalibrationWorkflow` 是旧的多阶段频率标定接口,内部委托给 `DampedSecantTracker`。
- `FrequencyStateMachine` + `FrequencyCalibrationRuntime` 是 V2 事件驱动接口:六状态协议
  含完整守卫/预算/监视器迟滞,支持 `save_run`/`load_run` 持久化与断点恢复。
- `SensingWorkflow.sweep()` 自动管理 CONFIG 同步(脉冲/硬件组参数改变时调
  `_sync_config()`),保证每次 `run()` 用的是当前参数对应的时间轴。
- 对瞬态协议，`pulse.rotation_angle`与`pulse.rabi_rate`是互斥模式；扫描其中一项时
  工作流会自动清除另一项。固定驱动强度的小角度扫描应先设置`rabi_rate`，再扫描
  `pulse.t_rabi_duration`，每个扫描点都会重新构造脉冲并估计响应核。
- `compare()` 把测量缓存复用给多个重建算法,避免重复跑 `mesolve`;`best` 字段按 RMSE
  自动选优,有真值(仿真 `flux_samples`)时有意义,无真值时回退到对 NaN 最稳健的方法。
- 要添加自定义工作流,继承 `Workflow` 抽象基类实现 `run() -> dict`,在内部按需编排
  实验/重建/标定三层。完整扩展指南见 {doc}`../extending`。
