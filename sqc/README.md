# SQC 超导量子控制重构项目说明

`sqc` 是本项目重构后的主包，目标是把原先集中在 `src/` 中的量子传感仿真代码，整理为更接近超导量子计算机全栈软件的结构。重构后的系统保留旧 notebook、`web_demo.py` 和 `src.*` 导入方式，同时把新业务逻辑放入分层清晰、可测试、可扩展的 `sqc/` 包。

当前重构覆盖 Phase 0 到 Phase 5：

- Phase 0：冻结旧物理行为，建立 pytest 与 baseline 回归体系。
- Phase 1：建立 `sqc` 包结构、基础数据容器、设备模型和兼容镜像。
- Phase 2：迁移主要实验协议，形成原生 experiment / runner / readout 接口。
- Phase 3：迁移波形重建、cryoscope、频率和 flux 标定框架。
- Phase 4：实现控制线失真、传递函数标定、预失真设计与验证 workflow。
- Phase 5：实现多 Z 线传递矩阵、芯片拓扑和双 qubit Z-crosstalk 补偿 demo。

## 可行性结论

截至当前状态，重构项目具备可运行、可验证、可继续开发的基础：

- 旧入口仍可使用：`src.qubit`、`src.signal`、`src.pulse`、`src.protocal`、`src.analysis` 保留兼容 facade。
- 新入口已经成型：新代码应优先从 `sqc.devices`、`sqc.control`、`sqc.experiments`、`sqc.reconstruction`、`sqc.calibration`、`sqc.hardware`、`sqc.workflows` 导入。
- 物理回归已冻结：Ramsey、diff echo、transient、cryoscope、predistortion、Z-crosstalk 都有 baseline。
- 工程可验收：全量测试当前通过 `77 passed`，仅保留 QuTiP 与 Matplotlib 的非功能性 warning。
- 风险可控：未完成的 LM 数值反演、真实硬件 cryoscope/transient 传递函数采集、cavity tomography 被明确标记为后续扩展，没有伪装成已完成能力。

一个实际边界需要注意：`Simulation.ipynb` 是旧 notebook，部分单元可能仍按旧返回值习惯写图，例如把 `Protocal(type=1).evolve()` 的 tuple 直接当作 y 值。核心 `src` 功能和 `web_demo.py` 已由 smoke 与回归测试覆盖，notebook 的展示单元若要完全顺滑，建议后续单独清理为新版 demo notebook。

## 总体架构

`sqc` 按全栈控制链路分层：

```text
sqc/
  devices/          物理器件：transmon、cavity、chip topology、coupler 抽象
  control/          控制对象：waveform、flux signal、pulse、pulse sequence
  simulation/       Hamiltonian 构造、QuTiP runner、结果容器、噪声工具
  experiments/      实验协议：Rabi、Ramsey、Diff Echo、Transient、Cryoscope
  reconstruction/   波形重建：kernel、Wiener、cryoscope、basis、Hammerstein
  calibration/      标定：flux response、频率、transfer function、predistortion
  hardware/         硬件链路：readout、control line、distortion、transfer matrix
  workflows/        端到端流程：predistortion validation、Z-crosstalk demo
```

推荐的数据流是：

```text
Device/QubitSpec
  -> Waveform/FluxSignal/Pulse
  -> Experiment 或 Workflow
  -> Runner/HamiltonianBuilder
  -> ExperimentResult
  -> Reconstruction/Calibration
  -> 新的 Waveform、CalibrationTable 或补偿策略
```

## 物理和工程原则

本项目聚焦软件仿真侧的控制闭环，而不是芯片设计、加工、低温系统或真实仪器驱动。它对应的是：

```text
器件模型 -> 控制波形 -> 哈密顿量仿真 -> 测量/读出 -> 重建/标定 -> 补偿/预失真
```

关键建模约定：

- 时间单位沿用旧代码，通常为 ns。
- `EC`、`EJ`、qubit frequency 沿用旧代码的角频率约定。
- flux 以 `Phi_0` 为单位。
- 新接口优先使用 `Waveform.samples`，旧接口可继续使用 `Signal.signal`。
- `ExperimentResult` 统一承载 `data`、`axes`、`metadata`、`config`。
- 物理行为改变必须经过 regression baseline 审核。

## 核心功能

### 1. 设备层

`QubitSpec` 是不可变参数对象，适合新代码保存纯物理参数。

```python
import numpy as np
from sqc.devices.transmon import QubitSpec

spec = QubitSpec(
    name="Q0",
    EC=2 * np.pi * 0.2,
    EJ=2 * np.pi * 15,
    T1=10000,
    T2=8000,
    flux_bias=0.0,
    n_levels=2,
)

frequency = spec.frequency()
sensitivity = spec.sensitivity()
```

`TransmonQubit` 保留旧的可变对象行为，适合兼容旧 notebook 和 `Protocal`。

```python
from sqc.devices.transmon import TransmonQubit

qubit = TransmonQubit(
    EC=2 * np.pi * 0.2,
    EJ=2 * np.pi * 15,
    T1=10000,
    T2=8000,
    n_levels=2,
)
```

`ChipTopology` 用于多 qubit / resonator 拓扑：

```python
from sqc.devices.chip import ChipTopology

chip = ChipTopology(qubits=[qa, qb])
dim = chip.hilbert_dim()
H0 = chip.hamiltonian_static()
```

### 2. 控制层

`Waveform` 是通用采样波形：

```python
from sqc.control.waveform import Waveform

waveform = Waveform(t_list=t, samples=x)
```

`FluxSignal` 是带物理语义的磁通波形，保留旧 `Signal(type=...)` 生成方式：

```python
from sqc.control.flux_signal import FluxSignal

phi = FluxSignal(type=2, t_list=t, amplitude=0.001, frequency=0.01)
```

`Pulse`、`CompositePulse` 和 `create_ramsey_pulse` 等函数来自旧代码迁移，仍可用于 Rabi、Ramsey、diff echo、cryoscope 等协议。

### 3. 实验层

原生实验接口返回 `ExperimentResult`：

```python
from sqc.experiments.ramsey import RamseyExperiment

result = RamseyExperiment(qubit).run()
p_e = result.data["p_e"]
tau = result.axes["tau"]
```

当前主要实验映射：

| 旧协议 | 新类 | 说明 |
|---|---|---|
| `Protocal(type=0)` | `RabiExperiment` | Rabi 振荡 |
| `Protocal(type=1)` | `RamseyExperiment` | Ramsey 干涉 |
| `Protocal(type=2)` | `DiffEchoExperiment` | differential echo 重建 |
| `Protocal(type=4)` | `TransientSensingExperiment` | 滑动 Ramsey / transient sensing |
| `Protocal(type=5)` | `CryoscopeExperiment` | cryoscope truncation demo |

### 4. 重建与标定

稳定能力已经迁入 `sqc.reconstruction`：

- `KernelEstimator`：从控制 pulse 自动估计 transient kernel。
- `WienerReconstruction`：滑动测量的线性 Wiener 反卷积。
- `DiffEchoReconstruction`：differential echo 的解析反演。
- `CryoscopeReconstruction`：cryoscope 相位微分重建。
- `HammersteinWienerReconstruction`：面向弱非线性响应的扩展入口。
- `LMReconstruction`：保留接口，但因 LM 收敛问题尚未完成，当前显式抛出 `NotImplementedError`。

标定层提供：

- `CalibrationTable`：一维插值表和反查表。
- `FluxResponseCalibration`：cryoscope flux response 标定。
- `QubitFrequencyCalibration`：qubit frequency 曲线标定。
- `TransferFunctionCalibration`：从测得 step response 或仿真 transfer model 构建传递函数表。

### 5. 控制线失真与预失真

P4 建立了 AWG 到芯片的控制线模型：

```python
from sqc.calibration.predistortion import PredistortionDesigner
from sqc.hardware.control_line import ControlLine
from sqc.hardware.distortion import SingleExponentialDistortion

line = ControlLine(
    name="Z0",
    kind="z",
    source="AWG0:CH1",
    target="Q0",
    transfer_function=SingleExponentialDistortion(amplitude=0.05, tau=20.0),
)

on_chip = line.apply(target_waveform)
awg = line.predistort(target_waveform, PredistortionDesigner(method="iir_inverse"))
corrected = line.apply(awg)
```

可用失真模型：

- `SingleExponentialDistortion`
- `MultiExponentialDistortion`
- `FIRDistortion`
- `IIRDistortion`
- `CustomTransferDistortion`

`PredistortionDesigner` 支持：

- `iir_inverse`
- `frequency_inverse`
- `fir_inverse`

端到端验证入口：

```python
from sqc.workflows.predistortion_validation import PredistortionValidationWorkflow

result = PredistortionValidationWorkflow(
    qubit=qubit,
    target_waveform=target,
    true_distortion=distortion,
).run()

factor = result["metrics"]["improvement_factor"]
```

### 6. Z-crosstalk 与多线传递矩阵

P5 建立多 Z 线串扰模型：

```python
from sqc.hardware.transfer_matrix import TransferMatrix

transfer = TransferMatrix.from_dc_matrix(
    dc_matrix=np.array([[1.0, 0.0], [0.04, 1.0]]),
    source_names=["QA", "QB"],
    target_names=["QA", "QB"],
)
```

约定为：

```text
Phi_j(omega) = sum_i H_ji(omega) V_i(omega)
```

Z-crosstalk demo：

```python
from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

result = ZCrosstalkWorkflow(
    chip=chip,
    flux_pulse_on_A=pulse_a,
    true_transfer_matrix=transfer,
).run()

print(result["fit_error_dB"])
print(result["compensation_factor"])
```

当前 demo 的补偿策略是：A 线产生对 B 的串扰，B 自己的 Z 线输出补偿波形来抵消 B 上的寄生 flux。它要求 `H_BB` 已知；如果没有 B 线自响应，workflow 会明确报错。

## 兼容旧代码

旧导入仍然有效：

```python
from src.qubit import TransmonQubit, Cavity, Coupled_System
from src.signal import Signal, CompositeSignal
from src.pulse import Pulse, CompositePulse, create_ramsey_pulse
from src.protocal import Protocal
from src.analysis import Analysis
```

兼容策略：

- `src.qubit`、`src.signal`、`src.pulse` 是镜像 facade。
- `src.protocal.Protocal(type=0/1/2/4/5)` 委托到 `sqc.experiments`。
- `src.analysis.Analysis` 委托到 `sqc.reconstruction` 和 `sqc.simulation`。
- 后续新增能力原则上先进入 `sqc/`，只有旧 notebook 必须调用时才增加 `src` facade。

## 安装和运行

推荐使用已有 conda 环境：

```powershell
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest -v
```

如果需要设置环境变量：

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
$env:MPLBACKEND='Agg'
$env:PYTHONIOENCODING='utf-8'
```

验证命令：

```powershell
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest tests/unit -v
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest tests/integration -v
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest tests/regression -m regression -v
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m pytest -v
```

`web_demo.py` smoke：

```powershell
$env:PYTHONIOENCODING='utf-8'
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -c "import web_demo"
```

## Baseline

`tests/baselines/` 保存冻结物理行为：

| Baseline | 覆盖内容 |
|---|---|
| `qubit_static.pkl` | qubit 静态属性 |
| `ramsey_default.pkl` | legacy Ramsey |
| `diff_echo_default.pkl` | legacy differential echo |
| `transient_default.pkl` | transient sensing |
| `cryoscope_default.pkl` | cryoscope |
| `predistortion_default.pkl` | P4 预失真闭环 |
| `z_crosstalk_default.pkl` | P5 Z-crosstalk 闭环 |

只有在明确改变物理行为时才允许重生成 baseline：

```powershell
& 'C:\Users\21034\anaconda3\envs\qutip-env\python.exe' -m tests.regression.generate_baselines
```

## 已知边界

- `LMReconstruction` 和相关数值反演函数保留接口，但当前因 LM 收敛问题显式抛出 `NotImplementedError`。
- `TransferFunctionCalibration` 支持 measured step response 或 known simulation model；真实 cryoscope/transient 硬件采集链路仍是后续工作。
- Cavity number splitting、Ramsey revival、Wigner tomography 在 P5 手册中是可选扩展，当前未实现。
- `web_demo.py` 仍是旧演示入口，没有新增 P4/P5 tab。
- `Simulation.ipynb` 是旧 notebook，建议后续另建新版 notebook 来展示 `sqc` 原生 workflow。

## 开发者文档

面向后续开发的接口、扩展规则、测试规范见：

- [DEVELOPER_GUIDE.md](DEVELOPER_GUIDE.md)
