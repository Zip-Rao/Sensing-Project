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
闭环整定 → 编排运行。平台提供**两套编排接口**：旧的一次性管线（`FrequencyCalibrationWorkflow`）
和新增的事件驱动状态机（`FrequencyCalibrationRuntime` + `FrequencyStateMachine`）。
以下按逻辑顺序展开各层职责，并在最后详述两者的区别与选型。

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
（`"sweet"`/`"target"`/`"track"`）配合不同测频协议构成六状态事件驱动协议
（Acquire → Track → Verify → Lock + Reacquire），由 `FrequencyCalibrationRuntime`
编排。各策略的完整说明见 {doc}`../building_blocks/calibration`。旧的四阶段
pipeline 接口（`FrequencyCalibrationWorkflow`）仍可用，内部已委托给共享的
`DampedSecantTracker` 控制器。
```

### 第四步：编排运行 — 工作流层

前序三步都是在**单个磁通点**或**固定策略**下工作。实际标定需要根据搜索阶段切换测量协议
和驱动策略——粗阶段需要宽范围的 Ramsey，细阶段受益于快速瞬态法；正常追踪时需要驱动频率跟随
qubit 移动，验证阶段则需要冻结偏置做独立判定。平台提供两套编排方式。

#### 旧接口：一次性管线 `FrequencyCalibrationWorkflow`

{py:class}`~sqc.workflows.frequency_calibration.FrequencyCalibrationWorkflow`
将多个 `SinglePointFrequencyCalibration` 串联为有序管线，每阶段以前一阶段的最优
磁通和频率估计为初值。因为每个阶段的 `measure_method` 在构造时固定，管线是**单向、
不可逆的**——阶段之间只能前进，不能回退或重试。

两种构造方式：

- **默认 hybrid 预设**（不传 `stages`）：自动构建两阶段瞬态→Ramsey 管线。
  `switch_residual`（默认 $2\pi \cdot 5$ MHz）控制粗-细交棒阈值。
- **显式管线**（传 `stages=[...]`）：每阶段为一个 `CalibrationStage`，独立
  指定 `measure_method`、`step_method`、`epsilon_f`、`drive_policy` 等。

`run()` 返回合并的迭代历史（各行标 `phase`/`global_iter`/`cost`）、
`V_final`、`residual`、`converged` 等。

#### 新接口：事件驱动状态机 `FrequencyStateMachine` + `FrequencyCalibrationRuntime`

旧管线的问题在于**只有"成功前进"一条路**：如果粗搜阶段偏置越界、瞬态鉴频器失锁、
或验证未通过，管线无法自动回到粗搜阶段重试——它会继续跑后面的阶段，产生无意义的结果。
对于需要长期运行的真实实验场景（锁定状态需持续监测漂移、按计划审计、在失锁时自动恢复），
一次性管线完全不够。

V2 接口用一个**事件驱动的六状态协议**解决这些问题。它把标定过程拆成离散的状态，
每个状态有明确的职责和进入/退出条件，状态之间按守卫规则转移：

```text
Acquire ──→ Track ──→ Verify ──→ Lock
              ↑          │          │
              │          │          │
            Reacquire ◄─────────────┘
```

**六个状态各司其职**：

| 状态 | 职责 | 能改磁通吗？ |
|---|---|---|
| **Acquire** | 宽范围 Ramsey 测量绝对频率，得到初始"种子"估计 | 否 |
| **Track** | 用瞬态鉴频器做局部测频 + 割线灵敏度估计 + 磁通步进——**唯一**允许改偏置的状态 | **是** |
| **Verify** | 冻结当前偏置，用独立的 Ramsey 双扫判断是否达到最终容差 | 否（进入即冻结） |
| **Lock** | 长期稳频：低成本监测漂移，按计划触发 Ramsey 审计 | 否 |
| **Reacquire** | 失锁/越界/置信度不足后，重新跑宽范围捕获 | 否 |
| **SafeStop** | 预算耗尽/联锁/不可恢复故障 → 安全保持 | 仅 safe bias |

**关键设计规则**：

1. **只有 Track 能改磁通偏置**——其他状态连"微调"都不行。这是为了防止 Verify 阶段
   "边验证边改偏置"的自欺行为。
2. **Verify 进入瞬间冻结候选偏置**，直到退出都不变。偏置的实际值与冻结值偏差超过
   `bias_freeze_tolerance` 即报错。
3. **Lock 监测到漂移不能直接调偏置**——小漂移走 Lock→Verify 独立确认，确认后才可能
   去 Track；大跳变直接走 Lock→Reacquire 重新捕获。这保证了"锁定"意味着偏置从未被
   悄悄改动。
4. **解析 $f(\Phi)$ 只做仿真 oracle**，不进转移决策——状态机不假设你知道频率-磁通
   关系，所有决策基于实际测量结果。

**阈值体系**——四个层次控制状态转移的松紧：

```
epsilon_hold  <  epsilon_final  <  epsilon_enter  <  Delta_val
(10 kHz)         (100 kHz)         (5 MHz)           (20 MHz)
物理带宽目标     Verify 通过       候选进入阈值      局部有效窗口
```

- **候选条件**（Track → Verify）：$|\text{error}| + \text{uncertainty} \le \epsilon_\text{enter}$
- **验证条件**（Verify → Lock）：$|\text{error}| + \text{uncertainty} \le \epsilon_\text{final}$，且连续 $N_\text{verify}$ 次
- **局部有效条件**（是否必须退出 Track）：$|\text{error}| + \text{uncertainty} \le \Delta_\text{val}$

**Lock 的监视器迟滞**——这是防止噪声导致状态抖动的关键设计。Lock 状态使用低成本
瞬态监测器，其噪声特性不同于 Verify 的独立 Ramsey，因此不能简单套用同一个阈值：

| 监测结果 | 动作 |
|---|---|
| $U_\text{mon} \le \epsilon_\text{mon\_clear}$ | 一切正常，清零可疑计数 |
| $\epsilon_\text{mon\_clear} < U_\text{mon} \le \epsilon_\text{mon\_suspect}$ | 灰区：累计可疑计数，连续 $N_\text{mon\_suspect}$ 次 → Verify |
| $\epsilon_\text{mon\_suspect} < U_\text{mon} < \Delta_\text{mon\_reacquire}$ | 明确漂移，立即 Lock → Verify |
| $U_\text{mon} \ge \Delta_\text{mon\_reacquire}$ 或 reference lost | 大跳变，直接 Lock → Reacquire |

这样做的好处是：刚通过 Verify 的点不会因为监测器的单次噪声波动就退出 Lock；
真正的缓慢漂移会在连续命中灰区后被捕获；大的突变直接触发重捕获，不浪费验证次数。

**命令--事件架构**——状态机不直接做 I/O，而是通过命令和事件与外部执行器解耦：

```text
┌──────────────┐     command      ┌────────────┐     execute     ┌─────────────┐
│ StateMachine │ ───────────────→ │  Runtime   │ ──────────────→ │ SQCExecutor │
│  (纯逻辑)     │                  │ (编排层)    │                 │ (QuTiP 后端) │
│              │ ←─────────────── │            │ ←────────────── │             │
└──────────────┘     event        └────────────┘    result       └─────────────┘
```

- **命令**（状态机 → 外部）：`AcquireFrequency`、`TrackFrequency(bias, drive)`、
  `VerifyFrequency(frozen_bias, drive)`、`MonitorFrequency(locked_bias)`、`SafeHold(bias)`
- **事件**（外部 → 状态机）：`MeasurementSucceeded`（携带频率、不确定度、有效性标志）、
  `MeasurementTechnicalFailure`、`MeasurementRejected`、`TimerElapsed`、`InterlockTriggered`、
  `BudgetExhausted`、`CancelRequested`

每条事件携带 `command_id` 以匹配命令，保证幂等——状态机反复调用 `next_command()` 返回
同一命令，直到收到匹配事件后才转移到下一状态。

**持久化与断点恢复**：`FrequencyCalibrationRuntime` 提供 `save_run(dir)` 和
`load_run(dir, qubit)` 方法，将运行状态写入五个文件：

| 文件 | 内容 |
|---|---|
| `config.json` | 完整协议配置 |
| `commands.jsonl` | 每条命令及其对应事件 |
| `transitions.jsonl` | 每次状态转移（from / to / reason） |
| `checkpoint.json` | 状态机全量快照（可从此恢复） |
| `result.json` | 最终结果摘要 |

这对于长时间运行的实验至关重要——断电或异常退出后可以从 checkpoint 继续，不会丢失
已经完成的标定进度。

**与旧接口的关系**：两者共享同一个 `DampedSecantTracker` 控制律（阻尼割线步进公式），
数值行为一致。旧 `FrequencyCalibrationWorkflow` 的每个梯度阶段内部已委托给 tracker。
区别在于编排层——旧的是固定管线，新的是事件循环 + 状态转移表。

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

### 事件驱动状态机（V2 新接口）

上节的 `FrequencyCalibrationWorkflow` 适合"设好参数、跑完看结果"的场景。但真实实验中，
标定过程可能持续数小时，期间可能发生偏置漂移、测量失效、联锁触发等意外。V2 接口用一个
**可暂停、可恢复、可回退**的六状态协议来处理这些情况。

使用时只需构造 `FrequencyCalibrationConfig` 和 `FrequencyCalibrationRuntime`，
然后调用 `run()`。Runtime 内部自动完成：

1. 创建 `FrequencyStateMachine` 并启动
2. 循环：`next_command()` → `SQCExecutor` 执行测量 → `handle(event)` 转移状态
3. Track 状态下自动管理 `DampedSecantTracker` 生命周期（`initialize` → `propose` → `accept`）
4. 每 10 条命令自动打 checkpoint

```python
from sqc.workflows.frequency_runtime import FrequencyCalibrationRuntime
from sqc.workflows.frequency_state_machine import FrequencyCalibrationConfig

# ── 协议配置：所有阈值和预算集中管理 ──────────────────────────────────
config = FrequencyCalibrationConfig(
    # -- 阈值（均为角频率，rad·GHz）--
    epsilon_enter=2 * np.pi * 20e-3,      # 20 MHz  — 候选进入 Verify 的条件
    epsilon_final=2 * np.pi * 2e-3,       # 2 MHz   — Verify 通过的条件
    Delta_val=2 * np.pi * 50e-3,          # 50 MHz  — 局部有效窗口

    # -- 验证策略 --
    N_verify=2,                            # 需要连续 2 次通过才算验证成功

    # -- Lock 监视器迟滞 --
    epsilon_mon_clear=2 * np.pi * 5e-3,   # 5 MHz   — 监视器认为"干净"
    epsilon_mon_suspect=2 * np.pi * 10e-3,# 10 MHz  — 监视器认为"可疑"
    Delta_mon_reacquire=2 * np.pi * 30e-3,# 30 MHz  — 大跳变，直接重捕获
    N_mon_suspect=3,                       # 连续可疑 3 次 → Verify

    # -- 预算 --
    max_commands=200,                      # 总命令数上限
    max_reacquire_attempts=5,             # 重捕获次数上限
    max_wall_time=3600.0,                 # 1 小时墙上时间

    # -- Track 控制器参数 --
    damping=0.8,                           # 阻尼因子（<1 抑制过冲）
    first_bias_step=0.01,                  # 首次探测步长 (Φ₀)
    max_bias_step=0.02,                    # 每步最大偏置变化 (Φ₀)
)

# ── 运行 ──────────────────────────────────────────────────────────────
runtime = FrequencyCalibrationRuntime(qubit=qubit, f_target=f_target, config=config)
result = runtime.run()

# ── 结果解读 ──────────────────────────────────────────────────────────
print(f"终态: {result['state']}")               # "lock" / "safe_stop"
print(f"运行状态: {result['run_status']}")       # "calibrated" / "completed" / "failed"
print(f"最终频率: {result['f_final']/(2*np.pi):.6f} GHz")
print(f"候选偏置: {result['candidate_bias']:.6f} Φ₀")
print(f"命令数: {result['n_commands']}, 耗时: {result['elapsed']:.1f}s")

# ── 转移轨迹 ──────────────────────────────────────────────────────────
for t in result['transition_log']:
    print(f"  v{t['version']}: {t['from']} → {t['to']}  [{t['reason']}]")

# ── 持久化 ────────────────────────────────────────────────────────────
runtime.save_run("calibration_run_001")
# 写入 config.json / commands.jsonl / transitions.jsonl / checkpoint.json / result.json

# ... 数小时后，从断点继续 ...
# runtime2 = FrequencyCalibrationRuntime.load_run("calibration_run_001", qubit=qubit)
# runtime2.run()  # 状态机从最后一个 checkpoint 恢复
```

这段代码展示了 V2 接口的核心流程：配置阈值 → 启动运行 → 查看轨迹 → 持久化保存。
`transition_log` 中每条记录包含状态转移的起止状态和触发原因码（如 `TARGET_CANDIDATE`、
`DRIFT_SUSPECTED`、`REACQUIRE_LIMIT` 等），方便事后审计——为什么进入了 Reacquire？
连续几次 suspect 触发了 Verify？所有决策都有记录可查。

```{note}
**旧接口 vs 新接口**：`FrequencyCalibrationWorkflow` 是一次性管线（stage1→stage2→...），
适合快速原型和批量仿真。`FrequencyCalibrationRuntime` + `FrequencyStateMachine`
是事件驱动架构，支持状态回退（Track→Reacquire→Track）、长期 Lock 监视与自动审计、
持久化与断点恢复，适合长时间标定实验和多 scenario 可靠性测试。
两者底层共享同一个 `DampedSecantTracker` 控制律（阻尼割线步进公式），数值行为一致。
详见 {doc}`../building_blocks/workflows`。
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
- **V2 状态机结果**：`runtime.run()` 返回 `state`（终态，如 `"lock"`）、
  `run_status`（运行状态，如 `"calibrated"` 表示成功进入 Lock）、`f_final`（最终
  频率估计）、`candidate_bias`（锁定时的偏置值）、`transition_log`（完整转移轨迹，
  每条含 `from`/`to`/`reason`/`version`）、`n_commands`（消耗的命令数）、`elapsed`
  （墙上时间）。`transition_log` 是事后审计的关键——全部状态转移和触发原因一目了然。
- 各标定类的完整字段与方法选项见 {doc}`../building_blocks/calibration`；V2 状态机
  的阈值与配置见 {doc}`../building_blocks/workflows`。
