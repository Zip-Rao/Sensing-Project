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

```
用户代码层
  │
  ├─ session = Session(config)           ← 单点入口，持有 Config + qubits + chip
  │    ├─ session.reconfigure(...)       ← 一键覆盖任意层参数
  │    ├─ session.qubit("QA").apply_gate("X_pi")   ← qubit 级操作
  │    ├─ session.qubit("QA").measure(protocol="ramsey")  ← 单 qubit 测量
  │    └─ session.run(workflow)          ← 编排执行
  │
  ├─ Workflow 级 ──────────────────────────────────────────
  │   ├─ CalibrateWorkflow      校验 qubit freq / anharmonicity / T1 / T2
  │   ├─ SenseWorkflow          感测 + 波形重建 (可选中重建算法)
  │   ├─ PredistortWorkflow     失真测量 + FIR 设计 + 验证
  │   └─ GateWorkflow           对指定 qubit(s) 做指定门操作
  │
  └─ Notebook 级 ─────────────────────────────────────────
      ├─ %%cell magic 风格的单行 demo
      └─ 参数扫描 cell 模板
```

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

---

## §4 优先级与依赖

| 任务 | 优先级 | 预计改动量 | 依赖 | 风险 |
|---|---|---|---|---|
| P6a `reconfigure()` 补完 | 🔥 高 | ~30 行 | 无 | 低 — 纯新增参数，不改默认行为 |
| P6c notebook 参数扫描 | 🔥 高 | ~60 行 | P6a (推荐，非强制) | 低 — 纯 notebook |
| P6b.2 GateOperation | ⚠️ 中 | ~80 行 | 无 | 中 — 需与 `src/qubit.py` gate 语义对齐 |
| P6b.1 SensingPipeline | ⚠️ 中 | ~120 行 | P6a (构造器读 CONFIG) | 中 — 新 Workflow 需等价性测试 |
| P6b.1 TwoQubitGate | 💡 低 | ~100 行 | GateOperation | 中 — 需 ChipTopology + coupler |

---

## §5 与已有代码的关系

### 不动 `src/` (R1)

所有新代码在 `sqc/` 或 notebook 内。`src_mirror/` 按需更新 facade (保持 `src/` 的 `Protocal(type=N)` 入口依旧可用)。

### 与 P0–P5 交付的关系

- P6a 是 P0 (CONFIG) 的**完成态**——把 migration 做完
- P6b 是 P5 (Workflow) 的**泛化**——从硬编 ZCrosstalkDemo 到可组合管线
- P6c 是 P2 (Experiment) 的**用户侧展示**——把 API 能力暴露为可交互 demo

### 测试策略

- P6a: `test_config.py` 新增 `test_reconfigure_covers_all_layers`
- P6b: `test_gate_operations.py` 新增 (unit) + `test_sensing_pipeline.py` 新增 (integration)
- P6c: 手动在 notebook 里跑 (暂无 notebook 自动化测试)

---

## §6 不做的事

- **不引入 CLI / YAML 配置** (R10)：保持纯 Python 构造
- **不引入 ORM / 数据库**：ExperimentResult 的序列化保持 pickle/npz
- **不做 GUI**：Gradio demo (`web_demo_v2.py`) 已存在，本次不扩展
- **不把 `src/` 的 `Protocal` 改成正确拼写** (R4)

---

> **更新记录**: 2026-05-12 — 初次编写 (基于 P5 完成后的 CONFIG 审计 + 用户操作入口诊断)
