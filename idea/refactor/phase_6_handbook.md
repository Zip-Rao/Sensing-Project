# Phase 6 — 用户可操作接口补完

> **状态**: 规划阶段，待实施。本文件是 Phase 5 完成后发现的用户侧可用性差距的集中记录。
> **前置**: P0–P5 全部完成。
> **配套文件**:
> - 全局配置规范: [CLAUDE.md §R9](../../CLAUDE.md)
> - 重构总方案: [_refactor_plan.md](_refactor_plan.md) §8 (兼容层契约)、§13 (子代理执行模式)
> - 当前交接状态: [_handoff_state.md](_handoff_state.md)

---

## §1 问题诊断: P5 交付后的"可操作感"缺口

### 1.1 CONFIG 层: 半完成 migration

commit `afe8d87` 引入 `ReconstructionConfig / SimulationConfig / TransmonDefaults / ControlLineDefaults` 四个 dataclass，但只在 **AWG + PulseConfig** 两层真正接了代码（12 个模块读 `CONFIG.pulse.* / CONFIG.awg.dt`）。其余四层字段全部为"声明但未接线"状态。

commit `747cc0d` (Phase 6 pre-work) 修复了 reconstruction 层的 5 个文件，将它们从字面默认值改为 `field(default_factory=lambda: CONFIG.reconstruction.*)`。但以下字段仍处于 **声明但无代码读取** 状态：

| 字段 | CONFIG 默认 | 当前状态 |
|---|---|---|
| `transmon.EC / EJ / T1 / T2 / flux_bias / n_levels` | 0.2 / 15 / 10000 / 8000 / 0 / 3 | `to_dict()` 从未被任何代码调用 (除 config.py 自身 docstring) |
| `transmon.F01_RANGE / EJ_EC_RANGE / ALPHA_RANGE` | Gao 推荐区间 | 无运行时 sanity check，CLAUDE.md R6 纯人工检查 |
| `simulation.store_states / atol / rtol` | False / 1e-8 / 1e-6 | `runner.py` 和 `numerical_inverse.py` 各自 hardcode |
| `control_line.impedance / attenuation_db / delay` | 50 / 20 / 0 | `control_line.py` 自己 hardcode 了一份 |
| `reconstruction.lm_n_basis / lm_max_iter / lm_tol / lm_mu_init / lm_basis_type / lm_lambda` | 已在 `747cc0d` 接线 | ✅ 已修复 |
| `reconstruction.cryoscope_tau` | 100.0 | ✅ 已在 `747cc0d` 接线 |
| `reconstruction.lambda_reg / stim_amplitude / stim_width` | 10.0 / 0.0215 / 3.0 | ✅ 已在 `747cc0d` 接线 |

### 1.2 `reconfigure()` 覆盖不足

当前 [sqc/config.py:234-258](sqc/config.py#L234) 只接受 AWG + Pulse 参数：

```python
def reconfigure(sample_rate=..., t_rabi_duration=..., t_global_start=..., t_global_end=...):
```

用户想改 `lambda_reg` 或 `stim_amplitude` 必须绕道：

```python
# 方法 A — 直接改 CONFIG（只影响后续构造，已存在的实例不动）
CONFIG.reconstruction.lambda_reg = 5.0

# 方法 B — 构造时传参（需要知道每个类的字段名）
WienerReconstruction(lambda_reg=5.0)

# 方法 C — 手动构造新 Config（过于冗长）
from sqc.config import Config, ReconstructionConfig
cfg = Config(reconstruction=ReconstructionConfig(lambda_reg=5.0))
```

**缺少一个统一的"单点调参入口"**。

### 1.3 demo notebook 里"可调把手"是隐式的

[Simulation_sqc.ipynb](../../Simulation_sqc.ipynb) 每个实验 cell 用默认参数：
- 改 flux 幅度 → 用户必须知道 `FluxSignal(amplitude=...)` → `TransientSensingExperiment(flux_signal=...)` 三级嵌套
- 改 λ → 用户必须知道 `WienerReconstruction(lambda_reg=...)`
- 改 qubit 工作点 → 用户必须知道 `TransmonQubit(flux=..., EC=..., EJ=...)`
- 没有任何 cell 演示**扫参数看效果**的工作流

### 1.4 Workflow 层缺少面向用户的编排接口

当前 `Workflow` ABC ([sqc/workflows/base.py](sqc/workflows/base.py)) 只有一个 `run()` 方法，返回 `dict`。实际只有一个实现 `ZCrosstalkWorkflow`。

缺失的能力：
- **任务选择**：用户无法在 workflow 里选"只做重建"还是"重建 + 标定 + 预失真"全套
- **qubit 级操作**：无接口支持"对 QA 施加 CZ gate、对 QB 施加 Ramsey 读取"
- **结果消费**：`run()` 返回裸 `dict`，无类型安全、无可视化快捷入口

---

## §2 目标状态: 用户视角的操作入口

### 2.1 总体设计

原设想分三层（Session → Workflow → Notebook），P6d 将 Session 并入 Workflow，**Workflow 即用户唯一可见入口**。

```
用户代码
  │
  └─ wf = SensingWorkflow()              ← 唯一入口 ★
       │
       ├─ wf.configure(...)              ← 5 组参数一次性设置
       ├─ wf.run(...)                    ← 执行管线 (测量 + 重建 + 标定)
       ├─ wf.sweep("param", [...])       ← 参数扫描
       ├─ wf.compare(methods=[...])      ← 算法 A/B 对比
       ├─ wf.plot()                      ← 一键出图
       │
       ├─ wf.pipeline([...])             ← 链式多阶段管线 (stub)
       ├─ wf.multi_qubit({...})          ← 多 qubit 编排 (stub)
       ├─ wf.crosstalk("QA", "QB")       ← Z-crosstalk 测量 (stub)
       │
       ├─ wf.find_optimal_work_point()   ← 搜索最佳偏置 (stub)
       ├─ wf.detectability_limit()       ← 最小可检测幅度 (stub)
       ├─ wf.noise_characterize()        ← 噪声表征 (stub)
       ├─ wf.cross_validate()            ← 交叉验证 (stub)
       ├─ wf.benchmark([...])            ← 批量基准测试 (stub)
       │
       ├─ wf.save() / wf.load()          ← 结果持久化 (stub)
       └─ wf.diff(other)                 ← 结果对比 (stub)

内部委托 (用户不可见):
  ├─ Experiment 层  (P2) — TransientSensingExperiment, RamseyExperiment, ...
  ├─ Reconstruction 层 (P3) — TransientReconstruction, CryoscopeReconstruction, ...
  ├─ Calibration 层 (P3c) — CryoscopeCalibration, DelayRamseyCalibration
  └─ Config 层 (P0/P6a) — reconfigure(), CONFIG
```

Session 的原始职责（持有 Config + qubits + chip + 参数）全部由 `SensingWorkflow` 承担。用户从头到尾只跟一个类交互。

### 2.2 关键约束

- **不动 `src/`** (R1)
- **不引入新依赖** (R10): 不用 pydantic、click、yaml
- **所有时间轴从 `CONFIG.awg.dt` 派生** (R9)
- **不破坏已有 Experiment/Reconstruction 的向后兼容**

---

## §3 实施计划

### P6a: 补完 `reconfigure()` — 单点调参入口

**文件**: [sqc/config.py](sqc/config.py) `reconfigure()` 函数

**扩展现有签名**：

```python
def reconfigure(
    # AWG (已有)
    sample_rate: float | None = None,
    # Pulse (已有)
    t_rabi_duration: float | None = None,
    t_global_start: float | None = None,
    t_global_end: float | None = None,
    # Reconstruction (新增)
    lambda_reg: float | None = None,
    stim_amplitude: float | None = None,
    stim_width: float | None = None,
    lm_n_basis: int | None = None,
    lm_max_iter: int | None = None,
    lm_lambda: float | None = None,
    # Simulation (新增)
    atol: float | None = None,
    rtol: float | None = None,
    # Transmon (新增)
    n_levels: int | None = None,
    # 通用回退
    **kwargs,
) -> Config:
```

**行为**：返回**新** `Config` 对象，不修改全局 `CONFIG`。只覆盖显式传入的参数，其余沿用当前 `CONFIG`。

**用途**：

```python
cfg = reconfigure(
    sample_rate=4.0,          # AWG dt 自动变为 0.25 ns
    lambda_reg=5.0,           # Wiener + Hammerstein 默认 λ
    n_levels=3,               # Transmon Fock 截断
    t_rabi_duration=20,       # 脉冲窗口
)
# 用 cfg.reconstruction.lambda_reg 构造 reconstruction
# 用 cfg.pulse.t_rabi 构造 experiment
# 用 cfg.transmon.to_dict() 构造 qubit
```

**验收**：
- `reconfigure(lambda_reg=5.0).reconstruction.lambda_reg == 5.0`
- `reconfigure()(不传参).reconstruction == CONFIG.reconstruction`
- `CONFIG` 本身不变 (singleton 不受影响)

### P6b: Workflow 级面向用户接口

#### 6b.1 任务编排 Workflow

**新文件**: `sqc/workflows/pipeline.py` — 可组合的实验管线

```python
from dataclasses import dataclass, field
from typing import Literal

@dataclass
class SensingPipeline(Workflow):
    """端到端感测管线: 测量 → 重建 → (可选)标定 → (可选)预失真.

    Parameters
    ----------
    qubit : TransmonQubit
    protocol : {"transient", "ramsey", "diff_echo", "cryoscope"}
        感测协议。
    reconstruction : {"wiener", "hammerstein", "lm", "none"}
        重建算法。"none" 只跑测量不重建。
    calibrate : bool
        是否在重建后追加频率标定。
    predistort : bool
        是否追加预失真设计。
    """

    qubit: object
    protocol: Literal["transient", "ramsey", "diff_echo", "cryoscope"] = "transient"
    reconstruction: Literal["wiener", "hammerstein", "lm", "none"] = "wiener"
    calibrate: bool = False
    predistort: bool = False

    # 各阶段的参数化
    experiment_kwargs: dict = field(default_factory=dict)
    reconstruction_kwargs: dict = field(default_factory=dict)

    def run(self) -> dict:
        """执行管线，返回每步的 ExperimentResult / ReconstructionResult."""
        ...
```

**关键设计决策**：
- `reconstruction="none"` 允许用户只跑测量拿原始 `p_e`/`Δp` 数据
- `calibrate=True` 在重建后跑 `QubitFrequencyCalibration` 得到标定表
- `predistort=True` 用标定表 + 失真模型设计预失真 FIR
- 每步结果存 `dict`，key 为阶段名

#### 6b.2 Qubit 级门操作接口

**新文件**: `sqc/control/gate_operations.py` — 面向用户的单/双 qubit gate

```python
@dataclass
class GateOperation:
    """对指定 qubit 执行指定门操作。

    Parameters
    ----------
    qubit : TransmonQubit
    gate : {"X_pi", "X_pi2", "Y_pi", "Y_pi2", "Z", "H", "identity"}
        门类型。X_pi = π 脉冲绕 X 轴 (DRAG)。
    frame : {"lab", "rotating"}
    omega_d : float or None
    """

    qubit: object
    gate: Literal["X_pi", "X_pi2", "Y_pi", "Y_pi2", "Z", "H", "identity"] = "X_pi2"
    frame: Literal["lab", "rotating"] = "rotating"
    omega_d: float | None = None

    def apply(self) -> Qobj:
        """返回该门的 unitary (Qobj)。不执行演化，只给算符。"""

    def simulate(self, t_list: np.ndarray | None = None) -> qutip.Result:
        """在 qubit 上执行该门并返回 mesolve 结果。"""


@dataclass
class TwoQubitGate:
    """双 qubit 门操作 (iSWAP / CZ)。

    Parameters
    ----------
    qubit_A, qubit_B : TransmonQubit
    gate : {"iSWAP", "CZ", "sqrt_iSWAP"}
    coupler : TunableCoupler or None
    """

    qubit_A: object
    qubit_B: object
    gate: Literal["iSWAP", "CZ", "sqrt_iSWAP"] = "iSWAP"
    coupler: object | None = None

    def apply(self) -> Qobj: ...
    def simulate(self, t_list=None) -> qutip.Result: ...
```

**验收**：
- `GateOperation(qubit=q, gate="X_pi2").simulate()` 生成 Ramsey 的第一个 π/2 脉冲 unitary
- `TwoQubitGate(qA, qB, gate="CZ").simulate()` 与 `src/qubit.py:Coupled_System` 行为一致

### P6c: Demo notebook 参数扫描示范

**文件**: [Simulation_sqc.ipynb](../../Simulation_sqc.ipynb) — 新增 cell(s)

**新增 cell "7. 参数扫描演示"**，覆盖用户最常调的四样：

```python
# 7.1 扫 flux 幅度 → 观察 Δp 幅度变化
# 7.2 扫 Wiener λ → 观察重建波形平滑度
# 7.3 扫 qubit flux 工作点 → 观察不同灵敏度
# 7.4 用 reconfigure() 单点改全局参数
```

每个子 cell 格式：
1. 一行 markdown 说明调什么
2. `for` 循环扫参数
3. `plt.subplots` 对比结果

**验收**：用户在 notebook 里复制粘贴任一子 cell 即可独立运行，不需要回头跑 setup cell。

### P6d: Workflow 层升级 — 用户可见的统一科研入口

**动机**：P6a/P6b/P6c 各自补一个缺口，但用户仍要知道 `FluxSignal` → `TransientSensingExperiment` → `TransientReconstruction` 三层类名才能跑通一个实验。P6d 将"Session"概念**并入 Workflow 层**，让 Workflow 成为唯一的用户可见入口：一个对象搞定参数配置、协议执行、参数扫描、A/B 对比，以及更复杂的科研端到端流程。

**核心理念**：Workflow 不是"一个实验的脚本"，而是**一次研究交互的会话**——你告诉它要研究什么，它帮你跑、帮你比、帮你出图。

**新文件**: `sqc/workflows/sensing.py` — `SensingWorkflow(Workflow)` 类

---

#### 6d.1 API 总览

```
SensingWorkflow
  │
  ├─ 🔥 configure(**kwargs) → self       # 分组设置所有参数 (full impl)
  ├─ 🔥 run(measure, reconstruct, calibrate) → WorkflowResult  # 执行 (full impl)
  ├─ 🔥 sweep(param, values) → SweepResult  # 参数扫描 (full impl)
  ├─ 🔥 compare(methods, measurement) → CompareResult  # 重建算法对比 (full impl)
  ├─ 🔥 plot()                            # 一键出图 (full impl)
  │
  ├─ ⚠️ pipeline(stages) → WorkflowResult  # 多阶段链式管线 (stub)
  ├─ ⚠️ multi_qubit(chip_spec) → self     # 多 qubit 编排 (stub)
  ├─ ⚠️ crosstalk(drive, sense) → WorkflowResult  # qubit 间 crosstalk 测量 (stub)
  │
  ├─ 💡 save(path) / load(path)           # 结果持久化 (stub)
  ├─ 💡 diff(other) → DiffReport          # 结果对比 (stub)
  ├─ 💡 benchmark(configs) → BenchmarkResult  # 协议/算法基准测试 (stub)
  │
  ├─ 💡 find_optimal_work_point() → float  # 搜索最佳 flux 偏置点 (stub)
  ├─ 💡 detectability_limit() → float      # 最小可检测信号幅度 (stub)
  ├─ 💡 noise_characterize() → NoiseReport # 零信号噪声表征 (stub)
  └─ 💡 cross_validate(k_folds) → CVResult # 重建交叉验证 (stub)
```

**🔥 = Phase 6 完整实现**　**⚠️ = 接口定义 + docstring，内部 raise NotImplementedError**　**💡 = 接口定义 + docstring，留给后续 phase**

---

#### 6d.2 Core 1: `configure()` + `run()` — 基础管线 (full impl)

```python
from sqc.workflows import SensingWorkflow

wf = SensingWorkflow()

# 一个方法改所有参数（5 组，按 section 分组暴露）
wf.configure(
    # ═══════════════════════════════════════════════════════════════
    # Qubit
    # ═══════════════════════════════════════════════════════════════
    EC=0.2, EJ=15.0, T1=10_000, T2=8_000,
    flux_bias=0.0, n_levels=3,

    # ═══════════════════════════════════════════════════════════════
    # Protocol & Signal — 选协议 + 定义要传感的波形
    # ═══════════════════════════════════════════════════════════════
    protocol="transient",       # "transient"|"ramsey"|"echo"|"cryoscope"|"delay_ramsey"
    signal_type=4,              # FluxSignal type 0-8
    signal_amplitude=0.01,      # Φ₀
    signal_width=10, signal_center=100,
    signal_rise=10, signal_fall=10,
    signal_custom=None,         # np.ndarray, for type=8

    # ═══════════════════════════════════════════════════════════════
    # Reconstruction — 选算法 + 调超参
    # ═══════════════════════════════════════════════════════════════
    reconstruction="wiener",    # "wiener"|"hammerstein"|"lm"|"none"
    lambda_reg=5.0,
    lm_n_basis=50, lm_max_iter=10, lm_basis_type="fourier",
    use_adjoint=True,

    # ═══════════════════════════════════════════════════════════════
    # Pulse
    # ═══════════════════════════════════════════════════════════════
    t_rabi_duration=20, t_global_start=-50, t_global_end=400,

    # ═══════════════════════════════════════════════════════════════
    # Hardware
    # ═══════════════════════════════════════════════════════════════
    sample_rate=2.0,
)

# 一键执行
result = wf.run(measure=True, reconstruct=True, calibrate=False)
wf.plot()
```

`run()` 三开关语义：

| 开关 | 行为 | 产出 |
|---|---|---|
| `measure=True` | 构造 Experiment → mesolve | `ExperimentResult` (p_e, Δp, kernel) |
| `reconstruct=True` | 构造 Reconstruction → reconstruct() | `FluxSignal` (B(t)) |
| `calibrate=True` | 先标定再测量 (仅 cryoscope/delay_ramsey) | `CalibrationTable` |

`measure=False, reconstruct=True` 时需传入 `previous_measurement=`。

`configure()` 部分更新：只传要改的参数，其余保持不变。返回 `self` 支持链式调用。

---

#### 6d.3 Core 2: `sweep()` — 参数扫描 (full impl)

研究最频繁的操作：扫一个参数，看结果怎么变。

```python
wf = SensingWorkflow()
wf.configure(
    protocol="transient", reconstruction="wiener",
    signal_type=4, signal_center=100,
)

# 扫信号幅度 → 自动跑 N 次，收集所有结果
sweep_result = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02, 0.05])

# sweep_result:
#   .param_path   → "signal.amplitude"
#   .values       → [0.005, 0.01, 0.02, 0.05]
#   .results      → [WorkflowResult, WorkflowResult, ...]
#   .metrics      → {"snr": [...], "rmse": [...], "peak": [...]}

# 自动出对比图
wf.sweep_plot(metric="snr")        # SNR vs amplitude 曲线
wf.sweep_plot(overlay=True)       # 所有 B(t) 叠加

# 支持扫任意参数路径
wf.sweep("reconstruction.lambda_reg", [0.1, 1, 10, 100])
wf.sweep("qubit.flux_bias", [-0.1, -0.05, 0, 0.05, 0.1])
wf.sweep("pulse.t_rabi_duration", [10, 20, 40, 80])
```

**参数路径语法**：`"group.field"`，与 `configure()` 的分组名对应。

| 路径前缀 | 对应 configure() 参数 |
|---|---|
| `signal.*` | amplitude, width, center, rise, fall, type |
| `reconstruction.*` | lambda_reg, lm_n_basis, lm_max_iter, lm_basis_type |
| `qubit.*` | EC, EJ, flux_bias, n_levels, T1, T2 |
| `pulse.*` | t_rabi_duration, t_global_start, t_global_end |
| `hardware.*` | sample_rate |

---

#### 6d.4 Core 3: `compare()` — A/B 对比 (full impl)

同一份测量数据，切换重建算法/超参，横向对比。

```python
wf = SensingWorkflow()
wf.configure(protocol="transient", signal=dict(amplitude=0.01))

# 方式 A：先测，再换算法对比
wf.run(measure=True, reconstruct=False)  # 只测不重建
comparison = wf.compare(
    methods=["wiener", "hammerstein", "lm"],
    # measurement 自动从上一次 run() 取
)

# 方式 B：显式传入测量结果
comparison = wf.compare(
    measurement=some_previous_result,
    methods=["wiener", "lm"],
    method_kwargs={
        "wiener": {"lambda_reg": 5.0},
        "lm": {"n_basis": 50, "use_adjoint": True},
    },
)

# comparison:
#   .methods      → ["wiener", "hammerstein", "lm"]
#   .signals      → {"wiener": FluxSignal, "hammerstein": FluxSignal, "lm": FluxSignal}
#   .metrics      → {"wiener": {"rmse": ..., "runtime": ...}, ...}
#   .best         → "lm"  (最小 RMSE)

wf.compare_plot()  # 多子图叠加对比 + 残差
```

---

#### 6d.5 协议 → 底层类自动映射

`configure(protocol=...)` 自动选择，用户不需要知道类名：

| `protocol=` | Experiment 类 | Reconstruction 类 | 支持 `calibrate` |
|---|---|---|---|
| `"transient"` | `TransientSensingExperiment` | `TransientReconstruction` | ❌ |
| `"ramsey"` | `RamseyExperiment` | `RamseyReconstruction` | ❌ |
| `"echo"` | `DiffEchoExperiment` | `EchoReconstruction` | ❌ |
| `"cryoscope"` | `CryoscopeExperiment` | `CryoscopeReconstruction` | ✅ `CryoscopeCalibration` |
| `"delay_ramsey"` | `DelayRamseyExperiment` | `DelayRamseyReconstruction` | ✅ `DelayRamseyCalibration` |

---

#### 6d.6 数据结构

```python
@dataclass
class WorkflowResult:
    config_snapshot: dict              # configure() 时的参数快照
    measurement: ExperimentResult | None
    reconstructed_signal: FluxSignal | None
    reconstruction_details: dict       # LM 时有迭代历史等
    calibration: object | None

@dataclass
class SweepResult:
    param_path: str                    # "signal.amplitude"
    values: list                      # [0.005, 0.01, ...]
    results: list[WorkflowResult]
    metrics: dict[str, list]          # {"snr": [...], "rmse": [...]}

@dataclass
class CompareResult:
    methods: list[str]
    signals: dict[str, FluxSignal]
    metrics: dict[str, dict]          # {"wiener": {"rmse": ..., "runtime": ...}}
    best: str
```

---

#### 6d.7 `plot()` 行为

根据 `run()` / `sweep()` / `compare()` 实际执行了什么自动生成合适的图：

- **run(measure=True)**: `p_e` + `Δp` vs scan
- **run(reconstruct=True)**: `Δp` + `B(t)` 重建波形 (有原始信号则虚线叠加)
- **run(calibrate=True)**: 标定曲线
- **sweep()**: 指标 vs 扫描参数曲线
- **compare()**: 多方法 B(t) 叠加 + 残差子图

---

#### 6d.8 接口预留 (stub only, raise NotImplementedError)

以下方法在 P6 只定义签名 + docstring，内部 `raise NotImplementedError("planned for P6.1+")`，供后续 phase 填充。定义它们是为了**确立 Workflow 作为统一科研入口的完整心智模型**——用户看一个类的 API 就知道"这个平台能做什么"。

##### ⚠️ `pipeline(stages: list[tuple]) → WorkflowResult`

链式多阶段管线：标定 → 测量 → 重建 → 预失真 → 验证。

```python
wf.pipeline([
    ("calibrate",   {"method": "qubit_frequency"}),
    ("measure",     {"protocol": "cryoscope"}),
    ("reconstruct", {"method": "wiener"}),
    ("predistort",  {"target": "step"}),
    ("verify",      {"protocol": "cryoscope"}),
])
# → WorkflowResult.stages["calibrate"], .stages["measure"], ...
```

##### ⚠️ `multi_qubit(chip_spec: dict) → self`

注册多 qubit chip，后续可按名字引用。

```python
wf.multi_qubit({
    "QA": {"EC": 0.20, "EJ": 15.0, "flux_bias": 0.0},
    "QB": {"EC": 0.22, "EJ": 14.5, "flux_bias": 0.0},
})
wf.gate(on="QA", gate="X_pi")           # (依赖 P6b.2 GateOperation)
wf.sense(on="QB", protocol="transient")  # 在 QB 上感测
```

##### ⚠️ `crosstalk(drive: str, sense: str, **params) → WorkflowResult`

测量 qubit 间 Z-crosstalk 传输函数。内部调用 P5 的 `ZCrosstalkWorkflow`。

```python
wf.multi_qubit({...})
wf.crosstalk(drive="QA", sense="QB", flux_amplitude=0.05)
# → WorkflowResult 含 H_BA(ω) 传输矩阵
```

##### 💡 `save(path: str)` / `load(path: str) → SensingWorkflow`

结果持久化：保存 `WorkflowResult` 到 `.pkl`，或从 `.pkl` 恢复。

##### 💡 `diff(other: WorkflowResult | str) → DiffReport`

两个 WorkflowResult 的差异对比：信号 RMSE、SNR 差异、参数差异。

##### 💡 `benchmark(configs: list[dict]) → BenchmarkResult`

对一组标准信号/协议/算法组合跑批量对比，出排名表。

```python
wf.benchmark([
    {"protocol": "transient", "reconstruction": "wiener", "signal": "gaussian"},
    {"protocol": "transient", "reconstruction": "lm",     "signal": "gaussian"},
    {"protocol": "ramsey",    "reconstruction": "wiener", "signal": "gaussian"},
    ...
])
# → BenchmarkResult.table() 打印对比矩阵
```

##### 💡 `find_optimal_work_point(metric="sensitivity", flux_range=(-0.5, 0.5)) → float`

在给定范围内扫描 flux_bias，找到使 dω/dΦ 最大 (灵敏度最高) 的工作点。返回最佳 `flux_bias` 值。

```python
optimal_flux = wf.find_optimal_work_point()
wf.configure(flux_bias=optimal_flux)
```

##### 💡 `detectability_limit(protocol, signal_type, amplitude_range, confidence=0.95) → float`

二分搜索最小可检测信号幅度。对给定协议 + 信号类型，找到使重建 SNR ≥ threshold 的最小幅度。

```python
min_amp = wf.detectability_limit(protocol="transient", signal_type=3)
# → "transient 协议对高斯脉冲的最小可检测幅度为 0.0003 Φ₀"
```

##### 💡 `noise_characterize(n_repeats=100, protocol="ramsey") → NoiseReport`

跑 N 次零信号测量，统计分析测量噪声的均值、方差、频谱。

```python
noise = wf.noise_characterize(n_repeats=50)
# → NoiseReport(mean=..., std=..., psd=..., normality_test=...)
```

##### 💡 `cross_validate(k_folds=5, reconstruction="lm") → CVResult`

K-fold 交叉验证重建算法的泛化能力：将测量数据分 k 折，轮流留一折做验证。

```python
cv = wf.cross_validate(k_folds=5)
# → CVResult(mean_rmse=..., std_rmse=..., per_fold=[...])
```

---

#### 6d.9 实现策略

**文件**: `sqc/workflows/sensing.py`

**核心逻辑** (~250 行 full impl + ~120 行 stubs)：

1. `__init__`: 初始化空 `_params` dict + 默认 Config
2. `configure(**kwargs)`: 按分组 key 合并到 `_params`，返回 `self`
3. `run(measure, reconstruct, calibrate)`: 根据 `_params` 构造 qubit/signal/experiment/reconstruction 并执行
4. `sweep(param_path, values)`: 解析 `"group.field"` 路径 → for 循环更新参数 + `run()` → 收集结果
5. `compare(methods, measurement, method_kwargs)`: 对每个 method 构造对应 Reconstruction → 收集结果
6. `plot()`: 读 `_last_result`，按结果类型分派到不同的画图逻辑
7. 其余方法 (pipeline, multi_qubit, crosstalk, save, load, diff, benchmark, find_optimal_work_point, detectability_limit, noise_characterize, cross_validate): **stub** — 完整 docstring + `raise NotImplementedError("planned for P6.1+")`

**测试**: `tests/unit/test_workflow.py` (~12 tests)

- `test_configure_stores_params`
- `test_configure_partial_update`
- `test_run_measure_only`
- `test_run_full_pipeline`
- `test_run_reconstruct_without_measure_raises`
- `test_protocol_mapping` — 5 种 protocol 全部映射正确
- `test_calibrate_skipped_for_transient`
- `test_sweep_returns_correct_count` — 扫 4 个值得到 4 个结果
- `test_sweep_metrics_computed` — SweepResult.metrics 非空
- `test_compare_multiple_methods` — compare 3 种方法各自重建
- `test_compare_without_measurement_raises`
- `test_stub_methods_raise_not_implemented` — 所有 stub 方法正确 raise

---

#### 6d.10 验收标准

1. 基础管线可运行：
```python
from sqc.workflows import SensingWorkflow
wf = SensingWorkflow()
wf.configure(protocol="transient", signal_type=4, signal_amplitude=0.01,
             reconstruction="wiener", t_rabi_duration=20)
result = wf.run()
wf.plot()
```

2. `sweep()` 正确执行并返回可用的 SweepResult：
```python
sweep = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02])
assert len(sweep.results) == 3
```

3. `compare()` 正确对比多种重建算法：
```python
cmp = wf.compare(methods=["wiener", "lm"])
assert "wiener" in cmp.signals and "lm" in cmp.signals
```

4. `configure()` 部分更新：改一个参数不影响其他。
5. 切换 `protocol=` 自动用正确的 Experiment + Reconstruction 类。
6. 所有 stub 方法有完整 docstring 且 `raise NotImplementedError`。
7. `pytest tests/unit/test_workflow.py` 全部通过。
8. `git diff --quiet master -- src/` 返回 0 (R1)。

---

## §4 优先级与依赖

| 任务 | 优先级 | 预计改动量 | 依赖 | 风险 |
|---|---|---|---|---|
| P6a `reconfigure()` 补完 | 🔥🔥 最高 | ~30 行 | 无 | 低 |
| **P6d `SensingWorkflow` (full impl)** | 🔥🔥 最高 | ~250+120 行 | P6a | 中 — 新 API，5 协议映射 |
| P6d sweep + compare (full impl) | 🔥🔥 最高 | 含在上行 | P6d core | 中 — 参数路径解析 |
| P6d 科研接口 stubs (9 个) | 🔥 高 | 含在上行 | P6d core | 低 — 只定义签名 |
| P6c notebook 参数扫描 | 🔥 高 | ~60 行 | P6d | 低 — 纯 notebook |
| P6b.2 GateOperation | ⚠️ 中 | ~80 行 | 无 | 中 |
| P6b.1 SensingPipeline | 💡 低 | ~120 行 | P6a | 低 — 被 P6d 覆盖 |

**执行顺序**: P6a → P6d → P6c。P6b 降级（P6d 覆盖其主要用例）。

---

## §5 与已有代码的关系

### 不动 `src/` (R1)

所有新代码在 `sqc/workflows/` 内。`src_mirror/` 按需更新 facade。

### 与 P0–P5 的关系

- P6a: P0 (CONFIG) 的**完成态**
- **P6d: P2+P3+P6a 的聚合层** — 将 Experiment (P2) + Reconstruction (P3) + Config (P6a) + 参数扫描 + A/B 对比全部封装为 Workflow 统一入口。**这是面向研究用户的"前端"**，底层类成为"后端"
- P6b.1/6b.2: P5 (Workflow) 的泛化，P6d 出现后降为可选中层
- P6c: 用 P6d API 改写 notebook demo cells

### 测试策略

- P6a: `test_config.py` 新增
- P6d: `tests/unit/test_workflow.py` 新增 (~12 tests, 覆盖 configure/run/sweep/compare + stubs)
- P6b: `test_gate_operations.py` (unit) + `test_sensing_pipeline.py` (integration)
- P6c: 手动 notebook 验证

---

## §6 不做的事

- **不引入 CLI / YAML 配置** (R10)：保持纯 Python 构造
- **不引入 ORM / 数据库**：ExperimentResult 的序列化保持 pickle/npz
- **不做 GUI**：Gradio demo (`web_demo_v2.py`) 已存在，本次不扩展
- **不把 `src/` 的 `Protocal` 改成正确拼写** (R4)

---

> **更新记录**:
> - 2026-05-12 — 初次编写 (基于 P5 完成后的 CONFIG 审计 + 用户操作入口诊断)
> - 2026-05-15 — P6d 从 `SensingSession` 改为 `SensingWorkflow`，并入 Workflow 层作为统一用户入口；新增 `sweep()` / `compare()` full impl + 9 个科研端到端接口 stub (pipeline, multi_qubit, crosstalk, save/load, diff, benchmark, find_optimal_work_point, detectability_limit, noise_characterize, cross_validate)
