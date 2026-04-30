# Sensing-Project 全栈化重构方案 — 主方案

> 版本:v1.0  
> 撰写日期:2026-04-30  
> 目标读者:具备 Python/QuTiP/超导量子计算基础的工程师或 AI agent  
> 目标:任何读者按本方案 + 6 份 phase handbook 可独立完成重构,无需进一步澄清

---

## 0. 文档导航

```
idea/refactor/
├── _refactor_plan.md       ← 本文件:主方案 (single source of truth)
├── phase_0_handbook.md     ← 测试基线 + 工具基础设施
├── phase_1_handbook.md     ← sqc/ 骨架 + ABC + 数据结构 + src/ 镜像
├── phase_2_handbook.md     ← 已实现协议 (case 1/2/4) 反向 wrapper + 实验对象化
├── phase_3_handbook.md     ← Track B 成果 (LM/Cryoscope/瞬态) 内化到 sqc/
├── phase_4_handbook.md     ← ControlLine + DistortionModel + Predistortion
└── phase_5_handbook.md     ← TransferMatrix + 双 qubit Z-crosstalk demo
```

阅读顺序:**先读本文件 §1–§13**(全局认知)→ **再读 phase 0 handbook**(必经先决条件)→ **逐个执行 phase 1–5 handbook**。

---

## 1. 项目定位与三大主线

### 1.1 项目本质

本项目是一个**研究型仿真框架**,以 flux-tunable Transmon qubit 为对象,研究**时变磁通/频率响应**的传感、重建、标定与反演。重构目标**不是**做通用量子算法模拟器,而是把它升级为**接近真实超导量子计算机控制栈架构**的科研平台。

### 1.2 三大主线 (优先级从高到低)

| # | 主线 | 物理目标 | 当前实现 | 重构后归属 |
|---|---|---|---|---|
| 1 | **波形重建** | 从 Ramsey/Cryoscope/瞬态测量恢复片上 Φ(t) 或 Δω(t) | case 1/2/4 + Wiener/LM/diff_echo | `sqc/reconstruction/` + `sqc/workflows/waveform_reconstruction.py` |
| 2 | **qubit 标定** | 标定 f₀₁、f(Φ)、灵敏度 κ、非线性系数 | `Calibration` 雏形 + 部分 case | `sqc/calibration/` + `sqc/workflows/qubit_calibration.py` |
| 3 | **波形预失真** | 建模 AWG→芯片 传递函数,设计 FIR/IIR/频域补偿 | **未实现** | `sqc/hardware/distortion.py` + `sqc/calibration/predistortion.py` + `sqc/workflows/predistortion_validation.py` |

### 1.3 与全栈架构的对应

参考真实 cQED 全栈:

```
Quantum algorithms              ← 不实现 (留接口)
Control software                ← 部分实现:experiments/, workflows/
Control electronics             ← 新增:hardware/electronics.py
Microwave signal processing     ← 部分实现:control/, hardware/distortion.py
Cryogenics and interconnects    ← 部分实现:hardware/control_line.py, transfer_matrix.py
Device                          ← 已有:devices/transmon.py, resonator.py, chip.py
```

本项目**重点落在** Device → Control pulse → Control-line transfer function → Calibration → Waveform reconstruction → Predistortion 这条链路。

---

## 2. 现状盘点 (基于代码事实,2026-04-30)

### 2.1 当前 src/ 文件清单

| 文件 | 行数 | 关键导出 |
|---|---|---|
| [src/qubit.py](../../src/qubit.py) | 645 | `TransmonQubit`、`Cavity`、`Coupled_System`、`ideal_iSWAP`、`simulate_iSWAP`、`ideal_CZ`、`simulate_CZ` |
| [src/signal.py](../../src/signal.py) | 282 | `Signal`、`CompositeSignal` |
| [src/pulse.py](../../src/pulse.py) | 691 | `Pulse`、`CompositePulse`、`create_pulse`、`create_ramsey_pulse`、`create_diff_echo_pulse`、`create_echo_pulse`、`create_cpmg_pulse`、`create_cryoscope_pulse` |
| [src/protocal.py](../../src/protocal.py) | 402 | `Protocal`(拼写故意保留)、`Calibration`、`IQ_readout` |
| [src/analysis.py](../../src/analysis.py) | 825 | `Analysis`、`generate_basis_functions`、`basis_function_decomposition`、`R`、`forward_simulation`、`compute_jacobian`、`compute_jacobian_finite_difference`、`levenberg_marquardt` |
| [web_demo.py](../../web_demo.py) | 639 | Gradio 界面,使用 protocol 0/1/4 |
| [Simulation.ipynb](../../Simulation.ipynb) | — | 主实验 notebook |

### 2.2 当前已实现的 protocol case

| case | 名称 | `Protocal.evolve()` 状态 | 物理含义 |
|---|---|---|---|
| 0 | Rabi 振荡 | ✓ | 单脉冲扫描 |
| 1 | Ramsey | ✓ | 含磁通信号 + π/2-τ-π/2 |
| 2 | 差分回波 | ✓ | π/2-[τ-π-τ′-(τ+t_int)-π-τ″]^k-π/2 |
| 3 | CPMG | × (`pass`) | 留空 |
| 4 | 瞬态磁场滑动测量 | ✓ | 滑动 Ramsey + kernel 计算 + Δp |
| 5 | Cryoscope | ✓ (实现不完整,测试不充分) | 截断扫描 + IQ 测量 |
| 6 | Cryoscope 标定 | × (`Calibration.calibrate` case 3 部分实现) | 扫描方波高度 |
| 7 | Cryoscope 测量 | × | 见 _TODO_master.md 1.1.3 |
| 8 | 瞬态频率标定 | × | 见 _TODO_master.md 1.2 |

### 2.3 当前已实现的 reconstruction 算法

| 算法 | 函数 | 状态 |
|---|---|---|
| Ramsey IQ 反演 | `Analysis.get_signal_from_ramsey_by_iq` | ✓ |
| Ramsey 解缠绕反演 | `Analysis.get_signal_from_ramsey_by_unwrap` | ✓ |
| 差分回波反演 | `Analysis.get_signal_from_diff_echo` | ✓ |
| Wiener 反卷积 | `Analysis.wiener_deconvolution` | ✓ |
| Hammerstein-Wiener | `Analysis.hammerstein_wiener_deconvolution` | ✓ |
| LM 数值反演 | `Analysis.numerical_inverse` + `levenberg_marquardt` | △ (Jacobian 收敛问题,见 _TODO_master.md 0.3) |
| Cryoscope 反演 | `Analysis.get_signal_from_cryoscope` + `get_h_from_phi` | △ (依赖 case 6/7 未完成) |

### 2.4 当前的设计债务 (重构后必须消除)

| # | 债务 | 影响 |
|---|---|---|
| D1 | `TransmonQubit.qubit_in_mag()` 把 `H_list`、`freq_coeffs`、`mag_signal`、`isinmag` 写到 qubit 实例 | 物理对象与实验状态耦合,无法并行多实验 |
| D2 | kernel 计算在 `Pulse.get_kernel()`、`CompositePulse.get_kernel()`、`Analysis.get_kernel()` 三处重复 | 维护成本,修改易遗漏 |
| D3 | `Protocal.evolve()` 是 200+ 行 match-case 巨函数,case 之间不能复用 | 难以加新 protocol;改一个 case 影响其他 |
| D4 | `IQ_readout()` 是模块级函数而非方法,且依赖 `qubit.H_list` 这种隐藏状态 | 接口不清晰 |
| D5 | `Signal` 一个类承担"控制波形包络"和"物理磁通信号"两种语义 | 概念混淆 |
| D6 | kernel 计算的刺激幅度硬编码 `0.0215`、宽度 `3 ns` | 物理参数变化时需手改 |
| D7 | 项目 0 个测试 | 重构无法验证物理结果不变 |
| D8 | LM Jacobian 伴随法与有限差分不一致 | 反演收敛问题(见 _TODO_master.md 0.3) |

---

## 3. 双 Track 编织模型 (核心)

### 3.1 两条平行 Track

| Track | 推进位置 | 内容 | 跟踪文档 |
|---|---|---|---|
| **Track A — 重构** | `sqc/`(新建) | 骨架/ABC/数据结构/反向 wrapper/实验对象化/迁移内化 | 本方案 + 6 份 handbook |
| **Track B — 功能开发** | `src/`(旧项目内) | LM 收敛修复 → 失真模型 → Cryoscope → 瞬态频率标定 | 现有 `_TODO_master.md` |

### 3.2 编织原则

1. Track B 的某项功能在 `src/` 里**实现并通过物理回归测试**后,Track A 的某个 phase 把它从"反向 wrapper"内化为 `sqc/` 原生实现,同时 `src/` 的相应符号转为"镜像导出"(`from sqc.* import *`)。
2. Track B 不进入 `sqc/`,**所有新功能开发先在 `src/` 完成**;`sqc/` 仅做骨架预留(ABC + `raise NotImplementedError`)和已稳定功能的内化。
3. 任何阶段 **`Simulation.ipynb` 和 `web_demo.py` 必须可运行**,不能因为重构破坏现有实验。

### 3.3 编织时序图

```
时间 →

Track B (src/):  [LM修复]──[失真模型]──[Cryoscope]──[瞬态标定]──···
                     │           │           │           │
                     ▼           ▼           ▼           ▼
Track A (sqc/): [P0]─[P1]─[P2]──[P3]─────[P4]────────[P3]──[P5]
                 测  骨  实验   反演内      失真内       标定内  多qbit
                 试  架  对象   化(LM)     化           化      crosstalk
```

虚线箭头表示"内化触发":Track B 一项功能完成后,触发 Track A 的相应 phase 把它迁移到 `sqc/`。

### 3.4 Track A 的 6 个 Phase

| Phase | 名称 | 触发条件 | 主要交付物 |
|---|---|---|---|
| **P0** | 测试基线 + 工具基础设施 | 立即开始 | `tests/`、`pytest.ini`、物理回归 baseline pickle |
| **P1** | sqc/ 骨架 + ABC + 数据结构 + src/ 镜像 | P0 完成 | `sqc/` 完整目录、所有 ABC、所有数据结构、`src/` 镜像通过 P0 测试 |
| **P2** | 已实现 case (1/2/4) 反向 wrapper + 实验对象化 | P1 完成 | `RamseyExperiment`、`TransientSensingExperiment`、`DiffEchoExperiment`、`KernelEstimator`、`ReadoutModel` |
| **P3** | LM/Cryoscope/瞬态标定 内化 | Track B 对应功能完成 | `LMReconstruction` 内化、`CryoscopeExperiment` 内化、`TransientFrequencyCalibration` 内化 |
| **P4** | ControlLine + DistortionModel + Predistortion | Track B 失真模型完成 | `hardware/control_line.py`、`hardware/distortion.py`、`calibration/predistortion.py`、`workflows/predistortion_validation.py` |
| **P5** | TransferMatrix + 双 qubit Z-crosstalk | P4 完成 | `hardware/transfer_matrix.py`、`workflows/z_crosstalk.py`、双 qubit demo notebook |

---

## 4. 目标架构

### 4.1 完整目录树

```text
Sensing-Project/
├── sqc/                            ← 新建主包 (Superconducting Quantum Control)
│   ├── __init__.py
│   ├── devices/                    ← 物理器件层
│   │   ├── __init__.py
│   │   ├── base.py                 ← Device ABC
│   │   ├── transmon.py             ← TransmonQubit, QubitSpec
│   │   ├── resonator.py            ← Resonator (= 现 Cavity)
│   │   ├── coupler.py              ← TunableCoupler ABC + 等效实现
│   │   └── chip.py                 ← ChipTopology, CoupledSystem
│   │
│   ├── hardware/                   ← 控制电子学/链路层
│   │   ├── __init__.py
│   │   ├── control_line.py         ← ControlLine (xy/z/readout)
│   │   ├── distortion.py           ← DistortionModel ABC + 子类
│   │   ├── transfer_matrix.py      ← TransferMatrix (多 qubit)
│   │   ├── electronics.py          ← AWG, LO, 抽象信号源 (P5+)
│   │   └── readout.py              ← ReadoutModel (理想 + IQ + 误判矩阵)
│   │
│   ├── control/                    ← 控制脉冲层
│   │   ├── __init__.py
│   │   ├── waveform.py             ← Waveform (= 现 Signal 的波形部分)
│   │   ├── flux_signal.py          ← FluxSignal (磁通语义)
│   │   ├── pulse.py                ← Pulse, CompositePulse
│   │   ├── sequence.py             ← PulseSequence + Ramsey/Echo/CPMG/Cryoscope 工厂
│   │   ├── schedules.py            ← (P5) 多通道时序对齐
│   │   └── gates.py                ← 单/双比特门 (DRAG/iSWAP/CZ)
│   │
│   ├── simulation/                 ← QuTiP 调用层
│   │   ├── __init__.py
│   │   ├── hamiltonian.py          ← HamiltonianBuilder (无副作用)
│   │   ├── runner.py               ← MesolveRunner, SlidingMeasurementRunner
│   │   ├── noise.py                ← 1/f 噪声、Lindblad 算符
│   │   └── result.py               ← ExperimentResult, MeasurementTrace
│   │
│   ├── experiments/                ← 实验协议对象
│   │   ├── __init__.py
│   │   ├── base.py                 ← Experiment ABC
│   │   ├── rabi.py                 ← RabiExperiment
│   │   ├── ramsey.py               ← RamseyExperiment
│   │   ├── echo.py                 ← EchoExperiment, DiffEchoExperiment
│   │   ├── cpmg.py                 ← CpmgExperiment (P3+)
│   │   ├── cryoscope.py            ← CryoscopeExperiment (P3+)
│   │   └── transient.py            ← TransientSensingExperiment
│   │
│   ├── calibration/                ← 标定 workflow
│   │   ├── __init__.py
│   │   ├── base.py                 ← Calibration ABC + CalibrationTable
│   │   ├── qubit_frequency.py      ← QubitFrequencyCalibration (Ramsey)
│   │   ├── flux_response.py        ← FluxResponseCalibration (Cryoscope/瞬态)
│   │   ├── readout.py              ← ReadoutCalibration (P5+)
│   │   ├── drag.py                 ← DragCalibration (可选,P5+)
│   │   ├── transfer_function.py    ← TransferFunctionCalibration (P4)
│   │   └── predistortion.py        ← PredistortionDesigner (P4)
│   │
│   ├── reconstruction/             ← 波形重建算法
│   │   ├── __init__.py
│   │   ├── base.py                 ← Reconstruction ABC
│   │   ├── kernel.py               ← KernelEstimator
│   │   ├── wiener.py               ← WienerReconstruction
│   │   ├── hammerstein.py          ← HammersteinWienerReconstruction
│   │   ├── numerical_inverse.py    ← LMReconstruction (P3 内化)
│   │   ├── cryoscope.py            ← CryoscopeReconstruction (P3 内化)
│   │   └── basis.py                ← 基函数生成 + 分解 + 正则化矩阵
│   │
│   └── workflows/                  ← 顶层科研流程
│       ├── __init__.py
│       ├── base.py                 ← Workflow ABC
│       ├── waveform_reconstruction.py
│       ├── qubit_calibration.py
│       ├── predistortion_validation.py
│       └── z_crosstalk.py          ← P5
│
├── src/                            ← 旧目录,永久镜像导出 (见 §8)
│   ├── __init__.py
│   ├── qubit.py
│   ├── signal.py
│   ├── pulse.py
│   ├── protocal.py
│   └── analysis.py
│
├── tests/                          ← P0 新建
│   ├── __init__.py
│   ├── conftest.py
│   ├── baselines/                  ← 物理回归 baseline pickle
│   │   ├── ramsey_default.pkl
│   │   ├── diff_echo_default.pkl
│   │   └── transient_default.pkl
│   ├── unit/
│   │   ├── test_devices.py
│   │   ├── test_control.py
│   │   ├── test_simulation.py
│   │   └── test_reconstruction.py
│   ├── integration/
│   │   ├── test_experiments.py
│   │   └── test_workflows.py
│   └── regression/
│       └── test_physics_baseline.py
│
├── docs/                           ← (可选,P4 后) 用户文档
├── idea/                           ← 设计文档
│   ├── _TODO_master.md
│   ├── _output.md                  ← 草稿,保留
│   └── refactor/                   ← 本方案所在
├── result/                         ← 实验结果(不变)
├── Simulation.ipynb                ← 主 notebook (不变)
├── web_demo.py                     ← Gradio demo (不变)
└── requirements.txt
```

### 4.2 模块职责一句话

| 层 | 模块 | 职责 |
|---|---|---|
| **devices** | transmon.py | Transmon 物理参数与算符,**不存任何实验状态** |
| | resonator.py | 多模腔模型,张量积算符 |
| | coupler.py | 可调耦合器(等效或显式) |
| | chip.py | 多 qubit 芯片拓扑 |
| **hardware** | control_line.py | 一条物理控制线(xy/z/readout)的元数据与传递函数 |
| | distortion.py | AWG→片上的失真模型 |
| | transfer_matrix.py | 多 qubit 串扰矩阵 H_ji(ω) |
| | readout.py | 投影/IQ/误判矩阵读出模型 |
| **control** | waveform.py | 通用时域波形(无物理语义) |
| | flux_signal.py | 磁通语义的信号(继承 Waveform) |
| | pulse.py | 微波脉冲 + 哈密顿量 list 格式 |
| | sequence.py | Ramsey/Echo/CPMG/Cryoscope 等脉冲序列工厂 |
| | gates.py | DRAG/iSWAP/CZ 门 |
| **simulation** | hamiltonian.py | 由 device + flux_signal + pulse 构造 QuTiP H_list,**无副作用** |
| | runner.py | mesolve/sliding measurement 的执行器 |
| | noise.py | 1/f 噪声生成 + Lindblad 算符 |
| | result.py | 统一结果数据结构 |
| **experiments** | rabi/ramsey/echo/cpmg/cryoscope/transient.py | 一个实验 = device + flux + sequence + readout + runner |
| **calibration** | qubit_frequency/flux_response/transfer_function/predistortion.py | 标定 workflow,产出 CalibrationTable 或 filter 系数 |
| **reconstruction** | kernel/wiener/hammerstein/numerical_inverse/cryoscope.py | 仅消费测量数据 + kernel + calibration table,产出重建波形 |
| **workflows** | waveform_reconstruction/qubit_calibration/predistortion_validation/z_crosstalk.py | 串联多个 experiment + reconstruction + calibration |

### 4.3 依赖矩阵 (高层只能依赖低层)

```
workflows  → experiments, calibration, reconstruction, simulation, control, hardware, devices
calibration → experiments, reconstruction, simulation, hardware, devices  (不依赖 workflows)
reconstruction → simulation, control, devices                            (不依赖 calibration/experiments/workflows)
experiments → simulation, control, hardware, devices                     (不依赖 calibration/reconstruction/workflows)
simulation → control, hardware, devices                                  (不依赖 experiments/calibration/reconstruction/workflows)
control → devices                                                        (允许引用 hardware.distortion 仅作类型注解)
hardware → devices
devices → (基础库)
```

**禁止反向依赖**。检测方法见 §9.4。

---

## 5. 核心数据结构

所有数据结构使用 `@dataclass(frozen=True)`(不可变)或显式声明可变字段;不引入 pydantic/attrs。

### 5.1 QubitSpec

```python
# sqc/devices/transmon.py
from dataclasses import dataclass, field

@dataclass(frozen=True)
class QubitSpec:
    """Transmon qubit 的纯参数描述,无实验状态。"""
    name: str                   # 唯一标识,如 "Q0"
    EC: float                   # 电容能量 (rad·GHz, 含 2π)
    EJ: float                   # 零磁通约瑟夫森能量 (rad·GHz)
    T1: float                   # 弛豫时间 (ns)
    T2: float                   # 退相干时间 (ns)
    flux_bias: float = 0.0      # 静态磁通偏置 (Φ₀)
    n_levels: int = 3           # Fock 截断
    
    def frequency(self, flux: float | None = None) -> float:
        """f₀₁(Φ) = √(8 EJ(Φ) EC) - EC, in rad·GHz."""
    
    def anharmonicity(self) -> float:
        """α = -EC."""
    
    def sensitivity(self, flux: float | None = None, 
                    delta: float = 1e-6) -> float:
        """κ = df₀₁/dΦ via central difference."""
```

**不变量**:`QubitSpec` **不可变**。所有"qubit 在某磁通下的频率"通过 `frequency(flux=...)` 调用,**不修改对象**。

### 5.2 Waveform / FluxSignal / MeasurementTrace

```python
# sqc/control/waveform.py
@dataclass
class Waveform:
    """通用时域波形,语义中立。"""
    t_list: np.ndarray           # shape (N,), 单位 ns
    samples: np.ndarray          # shape (N,), 单位由具体使用者定义
    metadata: dict = field(default_factory=dict)
    
    def value_at(self, t: float) -> float: ...
    def truncate(self, t_start: float, t_end: float) -> "Waveform": ...
    def copy(self) -> "Waveform": ...

# sqc/control/flux_signal.py
@dataclass
class FluxSignal(Waveform):
    """磁通信号 (单位 Φ₀)。samples 表示 ΔΦ(t)。"""
    pass

# sqc/simulation/result.py
@dataclass(frozen=True)
class MeasurementTrace:
    """一次或多次测量的结果。"""
    axis: np.ndarray             # 扫描参数轴 (如 tau_list, scan_list)
    p_e: np.ndarray              # 激发态概率
    p_e_iq: tuple[np.ndarray, np.ndarray] | None = None   # (I, Q) 双通道
    metadata: dict = field(default_factory=dict)
```

### 5.3 ControlLine

```python
# sqc/hardware/control_line.py
from typing import Literal

@dataclass
class ControlLine:
    """一条物理控制线。"""
    name: str                                 # 如 "Z0", "XY0"
    kind: Literal["xy", "z", "readout"]
    source: str                               # 如 "AWG0:CH1"
    target: str                               # 如 "Q0", "C01"
    transfer_function: "DistortionModel | None" = None
    
    def apply(self, awg_waveform: Waveform) -> Waveform:
        """awg_waveform → on_chip_waveform via transfer_function."""
```

### 5.4 PulseSequence

```python
# sqc/control/sequence.py
@dataclass
class PulseSequence:
    """统一的脉冲序列描述,跨通道。"""
    channels: dict[str, list["Pulse"]]        # 键为 ControlLine.name
    t_list: np.ndarray                        # 全局时间轴
    frame: Literal["lab", "rotating"]
    metadata: dict = field(default_factory=dict)
    
    def hamiltonian(self) -> list:
        """返回 QuTiP list 格式 [[H0, c0], [H1, c1], ...]"""
    
    def total_duration(self) -> float: ...
```

### 5.5 ExperimentResult

```python
# sqc/simulation/result.py
@dataclass
class ExperimentResult:
    """统一的实验结果容器。"""
    data: dict[str, np.ndarray]               # 命名数据 (p_e, varphi, B, ...)
    axes: dict[str, np.ndarray]               # 命名轴 (tau, scan, h, ...)
    metadata: dict                            # qubit_spec、experiment_class、git_sha 等
    config: dict                              # 调用时的输入参数
    
    def save(self, path: str) -> None: ...
    @classmethod
    def load(cls, path: str) -> "ExperimentResult": ...
```

`save/load` 使用 pickle 协议 5,**包含 ExperimentResult class 完全限定名,以便版本检查**。

### 5.6 TransferMatrix

```python
# sqc/hardware/transfer_matrix.py
@dataclass
class TransferMatrix:
    """多 qubit Z 线传递矩阵 H_ji(ω)。Phi_j(ω) = Σ_i H_ji(ω) V_i(ω)。"""
    elements: dict[tuple[str, str], np.ndarray]    # 键 (target, source) → 频域响应
    frequency_axis: np.ndarray                     # rad/ns
    time_axis: np.ndarray | None = None            # 可选,时域核函数
    
    def apply(self, source_to_voltage: dict[str, Waveform]) -> dict[str, FluxSignal]:
        """对每条 source 加 V_i,返回每个 target 的 Φ_j(t)。"""
    
    def diagonal(self) -> dict[str, np.ndarray]: ...
    def off_diagonal(self) -> dict[tuple[str, str], np.ndarray]: ...
```

### 5.7 CalibrationTable

```python
# sqc/calibration/base.py
@dataclass(frozen=True)
class CalibrationTable:
    """标定结果表,可序列化。"""
    qubit_name: str
    kind: Literal["f01", "f_phi", "kappa", "phi_h", "transfer_function"]
    inputs: np.ndarray                # 自变量(磁通/方波高度/...)
    outputs: np.ndarray               # 因变量(频率/相位/...)
    fit_params: dict                  # 拟合参数(多项式系数等)
    metadata: dict
    
    def evaluate(self, x: np.ndarray) -> np.ndarray:
        """插值/外推 outputs(x)。"""
    
    def inverse(self, y: np.ndarray) -> np.ndarray:
        """反函数(若单调)。"""
```

---

## 6. ABC 契约清单

所有 ABC 使用 `abc.ABC` + `@abstractmethod`(不引入 `typing.Protocol`,以减少 Python 版本依赖)。每个 ABC 必须定义 `__repr__` 包含关键参数,便于调试。

### 6.1 Device 层

```python
# sqc/devices/base.py
class Device(ABC):
    """所有物理器件的基类。"""
    
    @property @abstractmethod
    def name(self) -> str: ...
    
    @abstractmethod
    def hilbert_dim(self) -> int: ...
    
    @abstractmethod
    def hamiltonian_static(self) -> Qobj:
        """不含驱动/外场的孤立哈密顿量。"""
    
    @abstractmethod
    def collapse_operators(self) -> list[Qobj]: ...
```

### 6.2 Hardware 层

```python
# sqc/hardware/distortion.py
class DistortionModel(ABC):
    """AWG → 片上 的传递函数模型。"""
    
    @abstractmethod
    def apply(self, awg_waveform: Waveform) -> Waveform:
        """前向:输入 AWG 波形 → 输出片上波形。"""
    
    @abstractmethod
    def step_response(self, t: np.ndarray) -> np.ndarray: ...
    
    @abstractmethod
    def impulse_response(self, t: np.ndarray) -> np.ndarray: ...
    
    @abstractmethod
    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """复数频域响应 H(ω)。"""

# sqc/hardware/readout.py
class ReadoutModel(ABC):
    @abstractmethod
    def measure(self, state: Qobj | np.ndarray, 
                qubit_spec: QubitSpec) -> dict[str, float]:
        """返回 {"p_e": ..., "I": ..., "Q": ...}。"""
```

### 6.3 Control 层

```python
# sqc/control/pulse.py (复用现 Pulse, 接口稳定)
class PulseBase(ABC):
    @abstractmethod
    def hamiltonian(self) -> list:
        """QuTiP list 格式 [[op, coeffs_array], ...]。"""
    
    @property @abstractmethod
    def t_list(self) -> np.ndarray: ...
    
    @property @abstractmethod
    def frame(self) -> Literal["lab", "rotating"]: ...
```

### 6.4 Simulation 层

```python
# sqc/simulation/hamiltonian.py
class HamiltonianBuilder:
    """从 device + flux_signal + pulse 构造时变 H_list,无副作用。
    
    取代 src/qubit.py:qubit_in_mag (D1 债务)。
    """
    
    @staticmethod
    def build(qubit: QubitSpec, flux_signal: FluxSignal | None, 
              pulse: PulseBase | None,
              frame: Literal["lab", "rotating"] = "rotating",
              omega_d: float | None = None) -> tuple[list, np.ndarray]:
        """
        返回 (H_list, t_global)。
        H_list 形如 [[H0, coeffs], ...] 兼容 QobjEvo。
        QubitSpec 不被修改。
        """

# sqc/simulation/runner.py
class RunnerBase(ABC):
    @abstractmethod
    def run(self, *args, **kwargs) -> ExperimentResult: ...

class MesolveRunner(RunnerBase):
    def run(self, H_list: list, psi0: Qobj, t_list: np.ndarray,
            c_ops: list[Qobj], e_ops: list[Qobj]) -> ExperimentResult: ...

class SlidingMeasurementRunner(RunnerBase):
    def run(self, qubit: QubitSpec, flux_signal: FluxSignal,
            control_pulse: "CompositePulse",
            scan_list: np.ndarray | None = None) -> ExperimentResult: ...
```

### 6.5 Experiment 层

```python
# sqc/experiments/base.py
class Experiment(ABC):
    """一次实验 = 物理对象 + 控制脉冲 + 读出 + Runner。"""
    
    @abstractmethod
    def __init__(self, qubit: QubitSpec, flux_signal: FluxSignal | None,
                 readout: ReadoutModel, **kwargs): ...
    
    @abstractmethod
    def build_sequence(self) -> PulseSequence: ...
    
    @abstractmethod
    def run(self) -> ExperimentResult: ...
```

每个具体 Experiment 类必须实现:
- `RamseyExperiment(qubit, flux_signal, readout, t_rabi, tau_list, omega_d, phase1=0, phase2=0)`
- `DiffEchoExperiment(qubit, flux_signal, readout, t_rabi, tau_list, t_int, t_rep, k, omega_d)`
- `TransientSensingExperiment(qubit, flux_signal, readout, t_rabi, omega_d, scan_list=None)`
- `CryoscopeExperiment(qubit, flux_signal, readout, t_rabi, tau, trunc_list, omega_d)` (P3)
- `CpmgExperiment(...)` (P3+)
- `RabiExperiment(qubit, readout, t_rabi, omega_d)`

### 6.6 Calibration 层

```python
# sqc/calibration/base.py
class Calibration(ABC):
    @abstractmethod
    def __init__(self, qubit: QubitSpec, **kwargs): ...
    
    @abstractmethod
    def calibrate(self) -> CalibrationTable: ...
```

具体类:
- `QubitFrequencyCalibration` (Ramsey)
- `FluxResponseCalibration` (Cryoscope/瞬态/Ramsey 三种 backend)
- `TransferFunctionCalibration` (P4)
- `PredistortionDesigner` (P4,产出 IIR/FIR 系数而非 CalibrationTable)

### 6.7 Reconstruction 层

```python
# sqc/reconstruction/base.py
class Reconstruction(ABC):
    @abstractmethod
    def reconstruct(self, measurement: ExperimentResult, 
                    kernel: np.ndarray | None,
                    calibration: CalibrationTable | None,
                    **kwargs) -> Waveform | FluxSignal:
        """从测量结果反演物理信号。"""
```

具体类:
- `WienerReconstruction(lambda_reg, dt)`
- `HammersteinWienerReconstruction(qubit_spec, lambda_reg, dt)`
- `LMReconstruction(qubit_spec, control_pulse, basis_type, n_basis, lambda_reg, max_iter, tol)` (P3)
- `CryoscopeReconstruction(calibration: CalibrationTable, tau, dt, method='SG'|'inverse'|'diff')` (P3)
- `RamseyIQReconstruction(qubit_spec)`
- `RamseyUnwrapReconstruction(qubit_spec, k_span=3)`
- `DiffEchoReconstruction(qubit_spec, t_int, k)`

### 6.8 Workflow 层

```python
# sqc/workflows/base.py
class Workflow(ABC):
    @abstractmethod
    def run(self) -> dict[str, ExperimentResult | CalibrationTable | Waveform]: ...
```

---

## 7. 函数级新旧映射表

下表按 src/ 文件分组,列出每个公共符号迁移到 sqc/ 的位置。**所有迁移在 P1–P3 完成,P4–P5 仅新增**。

### 7.1 src/qubit.py → sqc/

| 旧符号 | 新位置 | 备注 |
|---|---|---|
| `TransmonQubit.__init__` | `sqc/devices/transmon.py: TransmonQubit.__init__` | 接口保持兼容,但内部用 `QubitSpec` |
| `TransmonQubit.calculate_frequency` | `QubitSpec.frequency` | 静态方法,接受 flux 参数 |
| `TransmonQubit.calculate_anharmonicity` | `QubitSpec.anharmonicity` | |
| `TransmonQubit.frequency_sensitivity` | `QubitSpec.sensitivity` | |
| `TransmonQubit.get_hamiltonian` | `TransmonQubit.hamiltonian_static` | 不变 |
| `TransmonQubit.get_hamiltonian_rwa` | `TransmonQubit.hamiltonian_rwa(omega_d)` | |
| `TransmonQubit.get_collapse_operators` | `TransmonQubit.collapse_operators` | |
| `TransmonQubit.generate_1f_noise` | `sqc/simulation/noise.py: generate_1f_noise` | 模块级函数 |
| `TransmonQubit.qubit_in_mag` | **拆分**:`HamiltonianBuilder.build` (P1 起) | **不再修改 qubit 对象**;返回 (H_list, t_list) |
| `TransmonQubit.qubit_under_mag` | `sqc/simulation/runner.py: 内部使用` | 私有,不导出 |
| `TransmonQubit.qubit_under_mag_hamiltonian` | 同上 | 私有 |
| `TransmonQubit.change_flux` | `QubitSpec.with_flux(new_flux) -> QubitSpec` | 返回新对象,不修改自身 |
| `TransmonQubit.optimal_work_point` | `QubitSpec.optimal_work_point` (静态) | |
| `TransmonQubit.simulate_gate` | `sqc/control/gates.py: simulate_single_qubit_gate` | 模块级函数 |
| `TransmonQubit.ideal_gate` | `sqc/control/gates.py: ideal_single_qubit_gate` | |
| `Cavity` | `sqc/devices/resonator.py: Resonator` | 类名变更 |
| `Coupled_System` | `sqc/devices/chip.py: CoupledSystem` | 拼写改正,接口保持 |
| `Coupled_System.simulate_iSWAP` | `sqc/control/gates.py: simulate_iSWAP` | 模块级函数 |
| `Coupled_System.prepare_ket11` | `sqc/control/gates.py: prepare_ket11` | |
| `ideal_iSWAP` | `sqc/control/gates.py: ideal_iSWAP` | 不变 |
| `simulate_iSWAP` (模块级) | `sqc/control/gates.py: simulate_iSWAP_simple` | 区分 |
| `ideal_CZ` | `sqc/control/gates.py: ideal_CZ` | |
| `simulate_CZ` (模块级) | `sqc/control/gates.py: simulate_CZ` | |

### 7.2 src/signal.py → sqc/

| 旧符号 | 新位置 | 备注 |
|---|---|---|
| `Signal` | **拆分**:`sqc/control/waveform.py: Waveform` + `sqc/control/flux_signal.py: FluxSignal` | type 0–8 → 工厂函数 `make_waveform(kind=..., t_list=..., **params)` |
| `Signal.value_at` | `Waveform.value_at` | |
| `Signal.truncate` | `Waveform.truncate` | 返回新对象(不再原地修改) |
| `Signal.update_signal` | `Waveform.with_params(**kwargs)` | 返回新对象 |
| `Signal.generate_basis_signal` | `sqc/reconstruction/basis.py: generate_basis_functions` | 复用 Analysis 中的同名函数 |
| `Signal.copy` | `Waveform.copy` | |
| `Signal.plot` | `Waveform.plot` | 保留 |
| `CompositeSignal` | `sqc/control/waveform.py: CompositeWaveform` | |

工厂映射(type → kind 字符串):

| 旧 type | 新 kind |
|---|---|
| 0 | `"zero"` |
| 1 | `"constant"` |
| 2 | `"sinusoidal"` |
| 3 | `"gaussian"` |
| 4 | `"asymmetric_pulse"` |
| 5 | `"double_peak"` |
| 6 | `"basis_expansion"` |
| 7 | `"wavepacket"` |
| 8 | `"custom"` |

### 7.3 src/pulse.py → sqc/

| 旧符号 | 新位置 | 备注 |
|---|---|---|
| `Pulse` | `sqc/control/pulse.py: Pulse` | 接口保持 |
| `CompositePulse` | `sqc/control/pulse.py: CompositePulse` | |
| `Pulse.get_kernel` / `CompositePulse.get_kernel` | `sqc/reconstruction/kernel.py: KernelEstimator.estimate(pulse, qubit, stim_amp=None, stim_width=None)` | **去重**;参数自动校准(_TODO 0.2) |
| `create_pulse` | `sqc/control/sequence.py: create_pulse` | |
| `create_ramsey_pulse` | `sqc/control/sequence.py: create_ramsey_pulse` | |
| `create_diff_echo_pulse` | `sqc/control/sequence.py: create_diff_echo_pulse` | |
| `create_echo_pulse` | `sqc/control/sequence.py: create_echo_pulse` | |
| `create_cpmg_pulse` | `sqc/control/sequence.py: create_cpmg_pulse` | |
| `create_cryoscope_pulse` | `sqc/control/sequence.py: create_cryoscope_pulse` | |

### 7.4 src/protocal.py → sqc/

| 旧符号 | 新位置 | 备注 |
|---|---|---|
| `Protocal.__init__` | (无对应,被 Experiment 替代) | |
| `Protocal.evolve` case 0 | `sqc/experiments/rabi.py: RabiExperiment.run` | |
| `Protocal.evolve` case 1 | `sqc/experiments/ramsey.py: RamseyExperiment.run` | |
| `Protocal.evolve` case 2 | `sqc/experiments/echo.py: DiffEchoExperiment.run` | |
| `Protocal.evolve` case 3 | `sqc/experiments/cpmg.py: CpmgExperiment.run` (P3+) | |
| `Protocal.evolve` case 4 | `sqc/experiments/transient.py: TransientSensingExperiment.run` | |
| `Protocal.evolve` case 5 | `sqc/experiments/cryoscope.py: CryoscopeExperiment.run` (P3) | |
| `Protocal.single_measurement` | `sqc/simulation/runner.py: MesolveRunner._single_measurement` | 私有 |
| `Protocal.sliding_measrement` | `sqc/simulation/runner.py: SlidingMeasurementRunner.run` | **拼写修正** |
| `Calibration.calibrate` case 0 | `sqc/calibration/qubit_frequency.py: QubitFrequencyCalibration` | |
| `Calibration.calibrate` case 1 | `sqc/calibration/flux_response.py: FluxResponseCalibration(method='ramsey')` | |
| `Calibration.calibrate` case 2 | `sqc/calibration/flux_response.py: FluxResponseCalibration(method='transient')` (P3) | |
| `Calibration.calibrate` case 3 | `sqc/calibration/flux_response.py: FluxResponseCalibration(method='cryoscope')` (P3) | |
| `IQ_readout` | `sqc/hardware/readout.py: IQReadoutModel.measure` | 类化 |

**`Protocal` 类名保留**(D3 债务在新代码中通过 `Experiment` 解决,旧 src/ 中不改名)。

### 7.5 src/analysis.py → sqc/

| 旧符号 | 新位置 | 备注 |
|---|---|---|
| `Analysis.get_expectation_values` | `sqc/simulation/result.py: extract_expectation` | 模块级 |
| `Analysis.get_population` | `sqc/simulation/result.py: extract_population` | |
| `Analysis.get_signal_from_ramsey_by_iq` | `sqc/reconstruction/wiener.py: RamseyIQReconstruction.reconstruct` | |
| `Analysis.get_signal_from_ramsey_by_unwrap` | `sqc/reconstruction/wiener.py: RamseyUnwrapReconstruction.reconstruct` | |
| `Analysis.get_signal_from_diff_echo` | `sqc/reconstruction/wiener.py: DiffEchoReconstruction.reconstruct` | |
| `Analysis.get_kernel` | (废弃,使用 `KernelEstimator`) | |
| `Analysis.wiener_deconvolution` | `sqc/reconstruction/wiener.py: WienerReconstruction.reconstruct` | |
| `Analysis.hammerstein_wiener_deconvolution` | `sqc/reconstruction/hammerstein.py: HammersteinWienerReconstruction.reconstruct` | |
| `Analysis.numerical_inverse` | `sqc/reconstruction/numerical_inverse.py: LMReconstruction.reconstruct` (P3) | |
| `Analysis.get_h_from_phi` | `sqc/calibration/flux_response.py: build_inverse_calibration` (P3) | |
| `Analysis.get_signal_from_cryoscope` | `sqc/reconstruction/cryoscope.py: CryoscopeReconstruction.reconstruct` (P3) | |
| `generate_basis_functions` | `sqc/reconstruction/basis.py: generate_basis_functions` | |
| `basis_function_decomposition` | `sqc/reconstruction/basis.py: basis_function_decomposition` | |
| `R` | `sqc/reconstruction/basis.py: regularization_matrix` | |
| `forward_simulation` | `sqc/reconstruction/numerical_inverse.py: _forward_simulation` (P3) | 私有 |
| `compute_jacobian` | `sqc/reconstruction/numerical_inverse.py: _compute_jacobian_adjoint` (P3) | |
| `compute_jacobian_finite_difference` | `sqc/reconstruction/numerical_inverse.py: _compute_jacobian_fd` (P3) | |
| `levenberg_marquardt` | `sqc/reconstruction/numerical_inverse.py: levenberg_marquardt` (P3) | |

---

## 8. 兼容层契约

### 8.1 永久镜像策略

`src/` 目录**永久存在**,作为 `sqc/` 的**镜像导出**。Notebook 和 web_demo.py 可以无限期使用旧 import,不强制弃用。

### 8.2 实现模式

每个 src/*.py 文件遵循**同一种模板**:

```python
# src/qubit.py (P1 之后)
"""
src.qubit — sqc.devices.transmon 的兼容镜像。

新代码应直接使用 sqc.devices.transmon。本文件仅为保持
Simulation.ipynb 和 web_demo.py 等历史代码可运行而存在。

本文件不应包含任何业务逻辑,所有改动须在 sqc/ 中进行。
"""
from sqc.devices.transmon import (
    TransmonQubit,
    QubitSpec,
)
from sqc.devices.resonator import Resonator as Cavity
from sqc.devices.chip import CoupledSystem as Coupled_System
from sqc.control.gates import (
    ideal_iSWAP,
    simulate_iSWAP_simple as simulate_iSWAP,
    ideal_CZ,
    simulate_CZ,
)

__all__ = [
    "TransmonQubit", "QubitSpec",
    "Cavity", "Coupled_System",
    "ideal_iSWAP", "simulate_iSWAP",
    "ideal_CZ", "simulate_CZ",
]
```

### 8.3 import 等价矩阵

**P1 完成后的等价关系**(每行左右两边在新代码里产生同一对象):

| 旧 import | 新 import |
|---|---|
| `from src.qubit import TransmonQubit` | `from sqc.devices.transmon import TransmonQubit` |
| `from src.qubit import Cavity` | `from sqc.devices.resonator import Resonator as Cavity` |
| `from src.qubit import Coupled_System` | `from sqc.devices.chip import CoupledSystem as Coupled_System` |
| `from src.signal import Signal` | `from sqc.control.flux_signal import Signal` (兼容别名) |
| `from src.signal import CompositeSignal` | `from sqc.control.waveform import CompositeWaveform as CompositeSignal` |
| `from src.pulse import Pulse, CompositePulse` | `from sqc.control.pulse import Pulse, CompositePulse` |
| `from src.pulse import create_ramsey_pulse` | `from sqc.control.sequence import create_ramsey_pulse` |
| `from src.protocal import Protocal` | `from sqc.experiments.legacy import Protocal` (P2 后,facade) |
| `from src.protocal import Calibration` | `from sqc.calibration.legacy import Calibration` (P3 后,facade) |
| `from src.protocal import IQ_readout` | `from sqc.hardware.readout import IQ_readout_legacy as IQ_readout` |
| `from src.analysis import Analysis` | `from sqc.reconstruction.legacy import Analysis` (P3 后,facade) |

`legacy.py` 文件在每个相关包内提供 facade,以保持旧 API 调用方式(类带 case dispatch)能继续工作。

### 8.4 关于 `Signal` 类的特殊处理

`Signal` 在旧代码里是"控制波形 + 物理磁通信号"两用对象。重构后:
- `sqc/control/waveform.py: Waveform` 是基类,通用波形语义。
- `sqc/control/flux_signal.py: FluxSignal(Waveform)` 加磁通语义。
- `sqc/control/flux_signal.py: Signal = FluxSignal` (兼容别名,旧 import 可用)。
- `Signal.type=0..8` 通过 `make_waveform(kind=..., **kwargs)` 工厂提供;`Signal(type=N, **kwargs)` 兼容构造在 `legacy.py` 实现。

### 8.5 拼写保留:Protocal vs Experiment

- 旧符号 `Protocal` **永久保留**(在 `sqc/experiments/legacy.py` 实现,转发给具体 Experiment 类)。
- 新符号 `Experiment` 是 ABC,所有新代码使用 `RamseyExperiment` 等具体类。
- 不批量改 `src/` 中的 `Protocal` 拼写,以避免 Notebook/demo 失效。

---

## 9. 测试基础设施

### 9.1 选型

- **pytest** + **pytest-xdist** (并行) + **pytest-cov** (可选,覆盖率)
- 不引入 hypothesis 等 property-based 框架
- 不引入 mock 库;若需要,使用 `unittest.mock`

### 9.2 目录布局

```
tests/
├── __init__.py
├── conftest.py                ← 全局 fixture (qubit_default, t_list_default 等)
├── baselines/                 ← 物理回归 baseline pickle
│   ├── ramsey_default.pkl
│   ├── diff_echo_default.pkl
│   └── transient_default.pkl
├── unit/                      ← 单元测试,纯函数/类
│   ├── test_devices.py
│   ├── test_waveform.py
│   ├── test_pulse.py
│   ├── test_hamiltonian.py
│   ├── test_runner.py
│   ├── test_reconstruction.py
│   └── test_calibration.py
├── integration/               ← 跨模块集成测试
│   ├── test_ramsey_pipeline.py
│   ├── test_transient_pipeline.py
│   └── test_workflows.py
└── regression/                ← 物理回归,锁定数值结果
    └── test_physics_baseline.py
```

### 9.3 物理回归 baseline 协议

**目的**:重构前后的物理结果必须一致。

**机制**:在 `src/` 改动之前,跑一次"baseline 生成脚本"(Phase 0 提供)记录三个标准协议在固定参数下的输出 pickle。重构期间,每次提交前跑 `pytest tests/regression -x`,确保结果误差在容限内。

**容限**:相对误差 ≤ 1e-6(QuTiP mesolve 的浮点误差量级)。

**baseline 锚定的参数**:

```python
# tests/conftest.py
@pytest.fixture
def qubit_default():
    return TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, n_levels=3
    )
```

### 9.4 依赖检测(可选)

在 `tests/conftest.py` 加一条 import-graph 检查,确保 `sqc.devices.*` 不 import 任何 `sqc.experiments.*` / `sqc.calibration.*` / `sqc.workflows.*` 模块。

---

## 10. 阶段路线总览

| Phase | Track A 内容 | 触发条件 | 估计行数 (新增) | 详情 |
|---|---|---|---|---|
| **P0** | 测试基线 + pytest 配置 + baseline pickle | 立即开始 | ~500 | [phase_0_handbook.md](phase_0_handbook.md) |
| **P1** | sqc/ 全骨架 + 所有 ABC + 数据结构 + src/ 镜像 | P0 完成 | ~1500 | [phase_1_handbook.md](phase_1_handbook.md) |
| **P2** | RamseyExperiment / DiffEcho / TransientSensing 内化;KernelEstimator 去重;ReadoutModel | P1 完成 + Track B LM 修复 | ~1000 | [phase_2_handbook.md](phase_2_handbook.md) |
| **P3** | LMReconstruction 内化;CryoscopeExperiment 内化;FluxResponseCalibration | Track B LM/Cryoscope/瞬态完成 | ~1500 | [phase_3_handbook.md](phase_3_handbook.md) |
| **P4** | ControlLine + DistortionModel + PredistortionDesigner + workflow | Track B 失真模型完成 | ~1500 | [phase_4_handbook.md](phase_4_handbook.md) |
| **P5** | TransferMatrix + 双 qubit Z-crosstalk demo | P4 完成 | ~1500 | [phase_5_handbook.md](phase_5_handbook.md) |

每个 phase handbook 包含:
1. 目标(一句话)
2. 前置条件(具体到 commit/文件/Track B 状态)
3. 任务清单(每条:文件路径、函数签名、估时、依赖)
4. 验收标准(可执行检查命令 + 数值阈值)
5. 测试要求(新增/修改的测试用例)
6. 风险与回滚
7. 输出物(文件清单 + 接口快照)

---

## 11. 风险登记 + 回滚策略

### 11.1 风险清单

| # | 风险 | 概率 | 影响 | 缓解措施 |
|---|---|---|---|---|
| R1 | 物理结果偏移 (浮点 / 算法等价性) | 中 | 高 | P0 物理回归 baseline + 每次提交跑 regression |
| R2 | `qubit_in_mag` 副作用消除导致 case 行为变化 | 中 | 高 | 旧符号在 facade 层保留副作用模拟 |
| R3 | QuTiP 版本升级使 H_list 格式变化 | 低 | 高 | requirements.txt 钉死 qutip 版本;CI 跑固定环境 |
| R4 | LM 收敛问题在迁移后被掩盖或恶化 | 中 | 中 | P3 起 baseline 包括 LM 反演结果 |
| R5 | Notebook 因 import 路径变化失效 | 高 | 中 | P1 完成立即跑 notebook 全部 cell;web_demo.py 端到端 smoke |
| R6 | 多 qubit 拓扑下 tensor product 维度爆炸 | 中 | 中 | P5 限制双 qubit + 单 cavity mode |
| R7 | 文档过期 / 与代码不同步 | 高 | 中 | 每个 phase 完成后,在主方案附录"已完成 phase 锚定的代码 commit SHA" |

### 11.2 回滚策略

**Phase 级回滚**:每个 phase 在独立 git 分支(`refactor/phase_N`)完成,合并 master 前必须满足:
- 所有 unit + integration + regression 测试通过
- Notebook 至少跑成功一个端到端 cell
- web_demo.py 启动 + Ramsey 功能可用

**单 commit 回滚**:任何 commit 引入 baseline 偏移 > 1e-6,该 commit 立即 revert。

**全局回滚**:若整个重构停滞,`src/` 永远是可工作版本(因为它从未被改变除了添加 facade)。删除 `sqc/` 即可回到原状。

---

## 12. 关键设计原则

按重要性排序:

1. **物理结果不变**。重构不修改任何物理算法的等价行为。任何性能优化或代码精简,都必须通过 baseline 回归测试。
2. **device/qubit 只描述物理参数**,不保存某次实验的临时状态(消除 D1 债务)。
3. **HamiltonianBuilder 无副作用**。输入 device + flux + pulse,输出 H_list 和 t_list,**不修改任何输入对象**。
4. **pulse/sequence 只描述控制**,不负责 kernel 或 reconstruction(消除 D2 债务)。
5. **experiment 是组合**,不是巨函数。一个 experiment = device + flux + sequence + readout + runner(消除 D3 债务)。
6. **calibration 是 workflow**,不是 protocol case 分支。
7. **reconstruction 只消费数据**,不做仿真;若需仿真(如 LM),通过依赖注入 device + control_pulse。
8. **hardware 层显式建模**控制电子学 → 控制线 → 失真 → crosstalk,反映真实 cQED 栈。
9. **永久镜像**。`src/` 不删除,旧代码无限期可工作。
10. **新代码用类型注解**,`from __future__ import annotations` 启用前向引用。
11. **不引入新依赖**(除测试)。只用 numpy/scipy/qutip/matplotlib/dataclasses。
12. **不引入配置文件**(yaml/toml)。所有参数通过 Python 对象传递。
13. **文档先行**:新增模块必须附 docstring(英文,遵循 numpy style),关键设计决策注释 `WHY`。
14. **不批量改拼写**(如 Protocal)以避免破坏现有 import。

---

## 13. 假设与待解决项

### 13.1 已采纳的假设(P1 实施时按此办)

| # | 假设 | 选择理由 |
|---|---|---|
| A1 | 测试框架使用 pytest + pytest-xdist | 项目当前 0 测试,pytest 是 Python 事实标准 |
| A2 | 物理回归用 pickle 锚定;容限 1e-6 相对误差 | mesolve 浮点稳定性约 1e-8,1e-6 留足缓冲 |
| A3 | ABC 用 `abc.ABC` + `@abstractmethod` | 比 `typing.Protocol` 兼容性更好,适合长期维护 |
| A4 | 数据结构用 `@dataclass`,不引入 pydantic/attrs | 减少依赖,标准库够用 |
| A5 | Phase 2 仅对 case 1/2/4 做实验对象化 | 这三个是当前 ✓ 状态;case 5/6/7/8 由 P3 在 Track B 完成后处理 |
| A6 | 不批量改 Protocal 拼写;新代码用 Experiment | 改拼写会破坏 web_demo.py 和 notebook |
| A7 | 文档语言中文,代码注释/docstring 英文 | 与现有 idea/ 风格一致 |
| A8 | 不引入 yaml/toml 配置 | 当前所有参数都通过 Python 对象传递,不需要配置层 |
| A9 | 不为 Track B 单独写 handbook | Track B 由 _TODO_master.md 跟踪;主方案只描述"内化触发条件" |
| A10 | requirements.txt 钉死 qutip 版本(P0 时执行) | 防止 QuTiP API 变更引入二次故障 |

### 13.2 待解决项(本方案不做硬性规定,留给执行者)

| # | 待解决项 | 决策时机 |
|---|---|---|
| Q1 | `sqc/__init__.py` 是否暴露顶层快捷 import (如 `from sqc import TransmonQubit`)? | P1 |
| Q2 | `Waveform.t_list` 是否应支持非均匀采样? | P1(若 mesolve 不支持非均匀,直接禁止) |
| Q3 | `ExperimentResult.save/load` 的 pickle 协议版本? | P0(建议 protocol=5) |
| Q4 | `KernelEstimator` 的刺激幅度自动校准算法(_TODO 0.2)? | P2 |
| Q5 | `LMReconstruction` 的 Jacobian 修复策略(_TODO 0.3)? | Track B,在 src/ 内完成 |
| Q6 | 是否在 P5 引入 schedule (类似 Qiskit Pulse 的多通道时序)? | P5 |
| Q7 | `tests/baselines/` 是否纳入 git? | P0(建议:小于 1 MB 的纳入,大于则用 git-lfs) |
| Q8 | 是否给 sqc 加一个对外 README? | P5 |

---

## 14. 附录:关键代码引用

为了让执行者快速定位,这里列出主方案中引用的关键代码位置。

### 14.1 现有代码痛点位置

| 痛点 | 文件:行 |
|---|---|
| qubit_in_mag 副作用 | [src/qubit.py:203-225](../../src/qubit.py#L203-L225) |
| Protocal.evolve 巨函数 | [src/protocal.py:46-186](../../src/protocal.py#L46-L186) |
| sliding_measrement 拼写 | [src/protocal.py:252](../../src/protocal.py#L252) |
| 三处重复 kernel | [src/pulse.py:160](../../src/pulse.py#L160) / [src/pulse.py:292](../../src/pulse.py#L292) / [src/analysis.py:142](../../src/analysis.py#L142) |
| kernel 硬编码 0.0215 | [src/pulse.py:319](../../src/pulse.py#L319) / [src/analysis.py:195](../../src/analysis.py#L195) |
| LM 伴随 Jacobian | [src/analysis.py:448-597](../../src/analysis.py#L448-L597) |

### 14.2 关键现有功能

| 功能 | 入口 |
|---|---|
| Ramsey 标准协议 | [src/protocal.py:61-99](../../src/protocal.py#L61-L99) |
| 差分回波 | [src/protocal.py:117-141](../../src/protocal.py#L117-L141) |
| 滑动测量 | [src/protocal.py:144-163](../../src/protocal.py#L144-L163) |
| Cryoscope 部分实现 | [src/protocal.py:164-183](../../src/protocal.py#L164-L183) |
| Wiener 反卷积 | [src/analysis.py:221-251](../../src/analysis.py#L221-L251) |
| LM 数值反演 | [src/analysis.py:272-297](../../src/analysis.py#L272-L297) |

### 14.3 全栈架构参考(已在 §1.3)

按真实 cQED 系统六层栈映射(Quantum algorithms / Control software / Control electronics / Microwave signal processing / Cryogenics and interconnects / Device)。本项目重点实现中间四层。

---

**文档结束**。下一步阅读 [phase_0_handbook.md](phase_0_handbook.md)。
