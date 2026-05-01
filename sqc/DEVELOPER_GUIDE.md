# SQC 开发者指南

本文面向后续维护和扩展 `sqc` 的开发者。它说明接口约定、模块职责、扩展流程、测试规范和兼容边界。项目总览和物理功能介绍见 [README.md](README.md)。

## 1. 开发原则

### 1.1 `sqc` 是新代码的唯一主线

新增功能优先放入 `sqc/`：

- 器件参数和算符放在 `sqc.devices`。
- 波形、pulse、sequence 放在 `sqc.control`。
- QuTiP 仿真运行和结果容器放在 `sqc.simulation`。
- 实验协议放在 `sqc.experiments`。
- 重建算法放在 `sqc.reconstruction`。
- 标定算法放在 `sqc.calibration`。
- 控制线、读出、失真、传递矩阵放在 `sqc.hardware`。
- 跨模块端到端流程放在 `sqc.workflows`。

`src/` 只作为兼容 facade。除非旧 notebook 或 `web_demo.py` 必须调用，否则不要把新业务逻辑写回 `src/`。

### 1.2 物理行为先冻结再修改

修改可能影响物理输出时，必须：

1. 先运行已有 regression，确认当前行为。
2. 修改实现。
3. 运行 unit、integration、regression。
4. 若 baseline 变化是有意的，说明物理原因并重生成对应 baseline。
5. 不要为通过测试随意放宽物理容差。

### 1.3 接口优先于脚本

新增能力应先设计成可导入、可测试的 Python API，再考虑 notebook 或 web demo 展示。

推荐返回：

- 单次实验：`ExperimentResult`
- 标定结果：`CalibrationTable`
- 波形结果：`Waveform` 或 `FluxSignal`
- 端到端流程：`dict`，其中至少包含原始对象和 `metrics`

## 2. 模块依赖规范

建议依赖方向如下：

```text
devices      -> base, qutip, numpy
control      -> waveform, devices 仅在 pulse 兼容层中允许
simulation   -> devices, control
experiments  -> devices, control, simulation, reconstruction
reconstruction -> control, simulation.result
calibration  -> devices, experiments, reconstruction, hardware
hardware     -> control
workflows    -> devices, control, hardware, experiments, reconstruction, calibration
src          -> sqc
```

规则：

- `devices` 不应依赖 `experiments` 或 `workflows`。
- `control.waveform` 应保持轻量，不依赖 QuTiP。
- `hardware.distortion` 和 `hardware.transfer_matrix` 应面向 `Waveform`，不要依赖旧 `Signal` 私有实现。
- `workflows` 可以组合多个模块，但不要把可复用算法埋在 workflow 私有函数里。
- `src` 只能向 `sqc` 单向导入。

## 3. 核心数据结构

### 3.1 `Waveform`

位置：`sqc.control.waveform.Waveform`

用途：通用采样波形。

字段：

| 字段 | 类型 | 说明 |
|---|---|---|
| `t_list` | `np.ndarray` | 时间轴 |
| `samples` | `np.ndarray` | 采样值 |
| `metadata` | `dict` | 附加信息 |

约定：

- `t_list.shape == samples.shape`
- 需要滤波或 FFT 的模型应检查时间轴均匀性。
- 方法返回新对象，不要原地修改输入，除非接口文档明确说明。

### 3.2 `FluxSignal`

位置：`sqc.control.flux_signal.FluxSignal`

用途：带磁通物理语义的波形，单位为 `Phi_0`。

兼容性：

- `signal` 是旧字段。
- `samples` 是新字段。
- 二者是 alias。

新代码应优先使用 `samples`。

### 3.3 `ExperimentResult`

位置：`sqc.simulation.result.ExperimentResult`

字段约定：

| 字段 | 内容 |
|---|---|
| `data` | 测量或计算得到的数组，例如 `p_e`、`kernel`、`delta_p` |
| `axes` | 坐标轴，例如 `t`、`tau`、`scan`、`t_samples` |
| `metadata` | 实验名、qubit spec、来源 |
| `config` | 运行参数 |

不要把核心数据塞进 `metadata`。`metadata` 只保存 provenance。

### 3.4 `CalibrationTable`

位置：`sqc.calibration.base.CalibrationTable`

用途：保存一维标定表，并提供插值和反查。

字段：

| 字段 | 内容 |
|---|---|
| `qubit_name` | qubit 名称 |
| `kind` | 标定类型，例如 `phi_h`、`transfer_function` |
| `inputs` | 输入轴 |
| `outputs` | 输出值 |
| `fit_params` | 拟合参数 |
| `metadata` | 方法、来源、版本 |

### 3.5 `DistortionModel`

位置：`sqc.hardware.distortion.DistortionModel`

必须实现：

```python
apply(waveform: Waveform) -> Waveform
step_response(t: np.ndarray) -> np.ndarray
impulse_response(t: np.ndarray) -> np.ndarray
frequency_response(omega: np.ndarray) -> np.ndarray
```

已实现模型：

- `SingleExponentialDistortion`
- `MultiExponentialDistortion`
- `FIRDistortion`
- `IIRDistortion`
- `CustomTransferDistortion`

### 3.6 `TransferMatrix`

位置：`sqc.hardware.transfer_matrix.TransferMatrix`

约定：

```text
Phi_j(omega) = sum_i H_ji(omega) V_i(omega)
```

关键接口：

```python
TransferMatrix.from_dc_matrix(dc_matrix, source_names, target_names)
TransferMatrix.apply(source_voltages)
TransferMatrix.H_ji(target, source)
TransferMatrix.diagonal()
TransferMatrix.off_diagonal()
```

`apply()` 的输入是 `dict[str, Waveform]`，输出是 `dict[str, FluxSignal]`。

## 4. 主要接口参考

### 4.1 Devices

| 接口 | 位置 | 说明 |
|---|---|---|
| `Device` | `sqc.devices.base` | 器件抽象基类 |
| `QubitSpec` | `sqc.devices.transmon` | 不可变 transmon 参数 |
| `TransmonQubit` | `sqc.devices.transmon` | 兼容旧代码的可变 transmon |
| `Cavity` / `Resonator` | `sqc.devices.transmon` / `sqc.devices.resonator` | cavity 兼容对象 |
| `ChipTopology` | `sqc.devices.chip` | 多 qubit / resonator 拓扑 |
| `Coupled_System` | `sqc.devices.transmon` | 旧双 qubit cavity 系统 |

### 4.2 Experiments

| 接口 | 位置 | 旧协议 |
|---|---|---|
| `RabiExperiment` | `sqc.experiments.rabi` | `Protocal(type=0)` |
| `RamseyExperiment` | `sqc.experiments.ramsey` | `Protocal(type=1)` |
| `DiffEchoExperiment` | `sqc.experiments.echo` | `Protocal(type=2)` |
| `TransientSensingExperiment` | `sqc.experiments.transient` | `Protocal(type=4)` |
| `CryoscopeExperiment` | `sqc.experiments.cryoscope` | `Protocal(type=5)` |

新增实验应继承 `sqc.experiments.base.Experiment`，实现：

```python
def build_sequence(self):
    ...

def run(self) -> ExperimentResult:
    ...
```

### 4.3 Reconstruction

| 接口 | 位置 | 说明 |
|---|---|---|
| `KernelEstimator` | `sqc.reconstruction.kernel` | 控制 pulse 到 kernel |
| `WienerReconstruction` | `sqc.reconstruction.wiener` | Wiener 反卷积 |
| `RamseyIQReconstruction` | `sqc.reconstruction.wiener` | Ramsey IQ 重建 |
| `RamseyUnwrapReconstruction` | `sqc.reconstruction.wiener` | Ramsey unwrap 重建 |
| `DiffEchoReconstruction` | `sqc.reconstruction.wiener` | differential echo 解析重建 |
| `CryoscopeReconstruction` | `sqc.reconstruction.cryoscope` | cryoscope 重建 |
| `HammersteinWienerReconstruction` | `sqc.reconstruction.hammerstein` | 弱非线性重建入口 |
| `LMReconstruction` | `sqc.reconstruction.numerical_inverse` | 接口保留，当前未完成 |

### 4.4 Calibration

| 接口 | 位置 | 输出 |
|---|---|---|
| `FluxResponseCalibration` | `sqc.calibration.flux_response` | `CalibrationTable` |
| `QubitFrequencyCalibration` | `sqc.calibration.qubit_frequency` | `CalibrationTable` |
| `TransferFunctionCalibration` | `sqc.calibration.transfer_function` | `CalibrationTable` |
| `PredistortionDesigner` | `sqc.calibration.predistortion` | inverse `DistortionModel` 或 AWG `Waveform` |

### 4.5 Hardware

| 接口 | 位置 | 说明 |
|---|---|---|
| `IdealProjectiveReadout` | `sqc.hardware.readout` | 理想投影读出 |
| `IQReadoutModel` | `sqc.hardware.readout` | IQ 读出模型 |
| `ControlLine` | `sqc.hardware.control_line` | AWG 到芯片控制线 |
| `DistortionModel` 子类 | `sqc.hardware.distortion` | 线性失真模型 |
| `TransferMatrix` | `sqc.hardware.transfer_matrix` | 多线串扰矩阵 |

### 4.6 Workflows

| 接口 | 位置 | 说明 |
|---|---|---|
| `PredistortionValidationWorkflow` | `sqc.workflows.predistortion_validation` | 预失真闭环验证 |
| `ZCrosstalkWorkflow` | `sqc.workflows.z_crosstalk` | 双 qubit Z 串扰提取和补偿 |

Workflow 的返回值建议包含：

- 输入或目标对象，例如 `target`
- 未补偿结果，例如 `on_chip_uncorrected`
- 补偿对象，例如 `awg_predistorted` 或 `compensation_pulse`
- 补偿后结果
- `metrics`

## 5. 扩展流程

### 5.1 新增实验协议

1. 在 `sqc/experiments/<name>.py` 新建 dataclass。
2. 继承 `Experiment`。
3. 明确输入 qubit、pulse、flux signal、runner。
4. `run()` 返回 `ExperimentResult`。
5. 加 unit 或 integration test。
6. 若需要旧入口，最后再给 `src.protocal` 增加 facade。

推荐模板：

```python
from dataclasses import dataclass
from sqc.experiments.base import Experiment
from sqc.simulation.result import ExperimentResult

@dataclass
class NewExperiment(Experiment):
    qubit: object

    def build_sequence(self):
        ...

    def run(self) -> ExperimentResult:
        return ExperimentResult(
            data={},
            axes={},
            metadata={"experiment": "NewExperiment"},
            config={},
        )
```

### 5.2 新增重建算法

1. 在 `sqc/reconstruction/<name>.py` 新建类。
2. 继承 `Reconstruction`。
3. 输入 `ExperimentResult`，输出 `FluxSignal`、`Waveform` 或数组。
4. 把正则化、窗口、滤波参数作为 dataclass 字段。
5. 加数值稳定性测试。

### 5.3 新增标定

1. 在 `sqc/calibration/<name>.py` 新建类。
2. 继承 `Calibration`。
3. `calibrate()` 返回 `CalibrationTable` 或明确的 filter/model 对象。
4. `fit_params` 必须足够重建模型。
5. 无测量数据时不要伪造结果，应抛 `NotImplementedError` 或要求显式输入。

### 5.4 新增硬件失真模型

1. 继承 `DistortionModel`。
2. 实现时域响应和频域响应。
3. `apply()` 返回新的 `Waveform`。
4. 检查时间轴均匀性。
5. 添加 step、impulse、frequency response 单测。

### 5.5 新增 workflow

1. 放入 `sqc/workflows/<name>.py`。
2. 继承 `Workflow`。
3. 组合已有 experiment、reconstruction、calibration、hardware API。
4. 返回 `metrics`，指标必须可测试。
5. 如果 workflow 代表新的物理交付物，增加 regression baseline。

## 6. 测试规范

测试分三层：

| 层级 | 目录 | 目的 |
|---|---|---|
| unit | `tests/unit` | 单模块逻辑、数学接口、轻量对象 |
| integration | `tests/integration` | 多模块协作、旧 facade 兼容、workflow 快速验证 |
| regression | `tests/regression` | 冻结物理行为和 baseline |

常用命令：

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
$env:MPLBACKEND='Agg'

& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest tests/unit -v
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest tests/integration -v
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest tests/regression -m regression -v
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest -v
```

当前全量期望：

```text
77 passed, 3 warnings
```

已知 warning：

- QuTiP coefficient RuntimeWarning，来自旧 diff echo baseline 路径。
- QuTiP `e_ops` FutureWarning。
- Matplotlib 中文字体 glyph warning，仅影响 `web_demo.py` 图中文字显示。

## 7. Baseline 规范

Baseline 文件在 `tests/baselines/`。

新增 baseline 的步骤：

1. 在 `tests/regression/generate_baselines.py` 增加 `_baseline_<name>()`。
2. 在 `BASELINES` 注册。
3. 新建 `tests/regression/test_<name>_baseline.py`。
4. 只生成新增 baseline，不要无意刷新旧 baseline。
5. 更新 `tests/baselines/README.md`。

生成命令：

```powershell
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m tests.regression.generate_baselines
```

只允许在物理行为有意改变时刷新旧 baseline。

## 8. `src` 兼容边界

`src` 是旧用户入口，不是新架构主线。

允许修改 `src` 的情况：

- 旧 notebook 或 `web_demo.py` 必须继续运行。
- 新 `sqc` 能力需要暴露给旧 API，且用户明确同意。
- 修复 facade 的返回 tuple 形状或导入兼容问题。

不建议修改 `src` 的情况：

- 新增实验、重建、标定、workflow。
- 重构内部算法。
- 添加未来功能的真实实现。

修改 `src` 后必须跑：

```powershell
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest tests/integration -v
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest tests/regression -m regression -v
```

## 9. 文档规范

文档分工：

- `sqc/README.md`：面向项目使用者和验收者，说明功能、原理、运行方式和边界。
- `sqc/DEVELOPER_GUIDE.md`：面向开发者，说明接口、扩展和测试规范。
- `tests/baselines/README.md`：只说明 baseline 文件和生成方式。
- `idea/refactor/*.md`：作为历史重构计划和阶段手册，原则上不再作为最新 API 文档。

新增公共接口时，至少更新：

1. `sqc/README.md` 的功能说明。
2. 本文件的接口表或扩展流程。
3. 对应测试。

## 10. 代码风格约定

- 优先使用 dataclass 表达纯配置或结果容器。
- 函数输入输出尽量显式，不依赖全局状态。
- 面向数组的计算使用 NumPy，避免手写低效循环。
- 需要插值、滤波、优化时优先使用 SciPy 标准实现。
- QuTiP 对象只在 device、simulation、experiment 必要位置出现。
- 不要在核心库函数中 `plt.show()`。
- 不要让测试依赖随机数；必须使用随机数时设置 seed。
- 错误边界要明确抛出 `ValueError`、`KeyError` 或 `NotImplementedError`，不要静默返回假结果。

## 11. 当前未完成事项

以下不是 bug，而是明确的后续开发任务：

- 修复并启用 LM 数值反演全流程。
- 将真实 cryoscope/transient 测得的传递函数采集链路接入 `TransferFunctionCalibration`。
- 为 P4/P5 增加新版 notebook 或 web demo tab。
- 实现 cavity number splitting、Ramsey revival、Wigner tomography。
- 添加噪声模型、并行优化、CPMG 频域重建、Volterra 非线性响应等扩展。

## 12. 快速检查清单

提交或交付前检查：

```text
[ ] 当前分支不是 master
[ ] 未误改 src，或 src 改动已说明风险
[ ] 新功能在 sqc 中实现
[ ] 公共接口已写入 README 或开发者指南
[ ] unit 测试通过
[ ] integration 测试通过
[ ] regression 测试通过
[ ] 若新增物理交付物，已新增 baseline
[ ] web_demo 或 notebook 兼容入口未被破坏
```
