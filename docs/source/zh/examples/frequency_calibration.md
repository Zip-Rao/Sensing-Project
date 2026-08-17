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

```{admonition} 源码扩展阅读
:class: tip
如果你希望从实现层理解事件驱动状态机，可阅读
{doc}`频率标定状态机源码导读 <frequency_state_machine_reading_guide>`。导读采用多遍阅读法，
先解释 State、Command 和 Event，再沿一次完整运行深入测量后端、Track 控制器与 Runtime。
```

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

##### 把它理解为持续运行的闭环系统

理解这套架构时，重点不是记住六个状态，而是先区分两个相互配合的循环。外层是
**协议循环**：它判断当前证据是否足以继续局部控制、需要独立验证，还是必须重新捕获；
内层是 **Track 控制循环**：它只根据相邻工作点的残差和偏置计算下一步控制量。

外层由 `FrequencyStateMachine` 管理，内层由 `DampedSecantTracker` 管理，
`FrequencyCalibrationRuntime` 在两者之间调度命令、测量和事件。这样，割线算法只负责
“下一步走多远”，不能自行宣布测量有效或标定成功；状态机只负责协议判断，不直接运行
QuTiP 或访问硬件。

状态图描述的是协议拓扑；要理解代码实现，更直接的方法是沿着一轮测量观察信息如何在
各层之间往返：

```text
用户 / UI
    │  配置、run、request_cancel
    ▼
FrequencyCalibrationRuntime ── next_command / handle ── FrequencyStateMachine
    │                                  纯协议决策
    ├── propose / accept ── DampedSecantTracker
    │                       纯 Track 控制律
    ├── execute ── SQCExecutor ── FrequencyMeasurement
    │                                │
    │                                ▼
    │                         TransmonQubit / QuTiP / hardware
    └── journal / checkpoint
```

用户从 runtime 进入系统，而不是直接操作状态机。runtime 创建状态机，向它询问下一条
命令，再协调能够完成该命令的组件。状态机内部没有 QuTiP 或硬件 I/O；它只回答两个
问题：当前允许执行什么操作，以及返回的证据应该把协议带到哪个状态。

启动后，状态机首先发出 `AcquireFrequency`。`SQCExecutor` 把这个协议级请求翻译成
宽范围 Ramsey `FrequencyMeasurement`，测量最终落到 `TransmonQubit`、QuTiP 仿真或
硬件适配器。executor 再把结果标准化为 `MeasurementSucceeded` 或失败事件。事件经
runtime 返回状态机 reducer 后，状态机才根据证据选择 Track 或 Verify。

###### Acquire：先建立绝对参考

Acquire 的目标不是调节偏置，而是用宽范围 Ramsey 建立当前频率、偏置和不确定度组成的
绝对参考。种子可靠且已满足候选条件时，协议直接进入 Verify；可靠但仍需修正时，才把
这个参考交给 Track。无效或歧义测量会有限重试，而不是让局部控制器在未知分支上工作。

###### Track：在可信局部范围内改变偏置

Track 比其他状态多一层控制计算。执行 `TrackFrequency` 前，runtime 先让
`DampedSecantTracker` 提出下一偏置和预测 drive。tracker 是纯数值控制器：它知道当前
与上一次残差，但不知道协议应该 Verify 还是 Reacquire。局部测量返回后，状态机先执行
有效性守卫：可信且达到候选条件则进入 Verify；可信但未达到则继续 Track；局部模型
失效则进入 Reacquire。这样数值控制器不能自行判断它依赖的测量模型是否仍然可信。

###### Verify：把控制与成功判定隔开

Verify 改变的是证据来源，而不是偏置。runtime 保持候选偏置冻结，executor 执行独立
double-sweep Ramsey。这避免让 Track 使用的同一个局部估计器既控制系统又证明自己已经
收敛。通过后进入长期 Lock；结果可靠但未达标时返回 Track；证据歧义时返回 Reacquire。

###### Lock：切换到长期监测时间尺度

进入 Lock 表示候选工作点已经被独立验证，并不表示 runtime 结束。Lock 沿同一分层路径
运行低成本 monitor，并按 `audit_interval` 安排独立 Ramsey 审计。它只能请求 Verify 或
Reacquire，不能绕过协议直接修改偏置：疑似小漂移先由 Verify 确认，明显跳变或参考丢失
则重新捕获。

###### Reacquire：恢复已经失效的局部知识

Track、Verify 和 Lock 都依赖当前分支、局部灵敏度、drive 位置或锁定参考。任一参考不再
可信时，Reacquire 清除不应继续沿用的局部状态并重新执行宽范围捕获。新种子恢复后，可按
距离目标的情况回到 Track 或直接进入 Verify；只有连续重捕获失败或超过次数限制，才进入
SafeStop。因此 Reacquire 是恢复闭环，不是失败终点。

##### Runtime、安全包络与中断恢复

从调度角度看，每一轮都遵循同一条路径：runtime 先检查取消信号，向状态机索取 pending
command；若处于 Track，再让 tracker 完善偏置和 drive；随后在 I/O 前估算并预留成本，
调用 executor，把返回事件交给 reducer，最后记录 journal 并按策略保存 checkpoint。

预算、中断、journal 和 checkpoint 包围整个循环，而不属于某个科学状态。它们是
runtime 的横切职责：可以在 I/O 前阻止新命令，保留已经返回的测量证据，并通过
`SafeStop → SafeHold` 安全退出，而不把持久化或线程协调逻辑塞进纯 reducer。

用户调用 `request_cancel()` 时只设置线程安全的协作式取消信号。runtime 在命令边界和
monitor 等待期间将它转换为 `CancelRequested`；同步 `executor.execute()` 已经开始时，
当前调用仍须返回后才能停止。`KeyboardInterrupt` 默认复用同一路径。进入 SafeStop 后，
系统还会实际下发 `SafeHold`，收到确认事件后才把安全保持记为完成。

checkpoint 则处理进程在运行中间退出的情况。尚未收到匹配事件的 pending command 会保留
原 `command_id`，恢复后以同一 ID 重放。因此协议提供 at-least-once 恢复语义；真实硬件
若要求避免重复副作用，executor 仍需按 `command_id` 去重。

##### 配置与接口速查

前面的叙事说明各层为什么这样协作；下面保留完整状态职责、守卫、阈值、消息和持久化字段，
供配置实验与排查运行记录时直接查阅。

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

Track 按固定顺序应用有效性守卫：后端 `out_of_range`、可选 detuning 捕获窗
`linear_range - guard_margin`、可选实验 `min_confidence`、割线灵敏度范围
`[S_min, S_max]`，以及探针失谐有效性界 `Delta_val`。目标残差只参与候选判断和
偏置更新，不能代替 probe detuning。确定性后端用
`uncertainty_source="deterministic_zero"` 标记零统计误差；这是仿真假设，不是测得的
实验置信度。

**阈值体系**——四个层次控制状态转移的松紧：

```
epsilon_hold  <  epsilon_final  <  epsilon_enter  <  Delta_val
(10 kHz)         (100 kHz)         (5 MHz)           (20 MHz)
物理带宽目标     Verify 通过       候选进入阈值      局部有效窗口
```

- **候选条件**（Track → Verify）：$|\widehat r| + z\sigma_r \le \epsilon_\text{enter}$
- **验证条件**（Verify → Lock）：$|\widehat r| + z\sigma_r \le \epsilon_\text{final}$，且连续 $N_\text{verify}$ 次
- **局部有效条件**（是否必须退出 Track）：$|\widehat\Delta| + z\sigma_\Delta \le \Delta_\text{val}$

Track 还可通过 `expected_sensitivity_sign` 显式给出当前单调分支的
$\chi_\Phi=\operatorname{sign}(\partial f/\partial\Phi)$。首次探索步据此选择方向，后续
割线若与该方向不一致，会在发出下一次测量前进入 Reacquire。设为 `None` 时保留旧的
隐式正方向约定。

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
在解释上述目标残差之前，Lock monitor 还会独立检查探针失谐
$|\widehat\Delta_\text{mon}|+z\sigma_{\Delta,\text{mon}}\le\Delta_\text{val}$，并核对实际
偏置和驱动。越出局部可信范围时直接 Reacquire，执行值与命令不一致时进入 SafeStop；
因此“小残差”不会掩盖一个已经失效的局部探针。

**命令--事件契约**——上述跨层路径通过以下稳定消息保持解耦：

- **命令**（状态机 → 外部）：`AcquireFrequency`、`TrackFrequency(bias, drive)`、
  `VerifyFrequency(frozen_bias, drive)`、`MonitorFrequency(locked_bias, drive)`、`SafeHold(bias)`
- **事件**（外部 → 状态机）：`MeasurementSucceeded`（携带频率、不确定度、有效性标志）、
  `MeasurementTechnicalFailure`、`MeasurementRejected`、`TimerElapsed`、`InterlockTriggered`、
  `BudgetExhausted`、`CancelRequested`

每条事件携带 `command_id` 以匹配 pending command。反复调用 `next_command()` 会返回
同一命令，直到收到匹配事件后才转移。checkpoint 恢复时 pending command 会使用相同 ID
重放，因此必须由 executor 按 `command_id` 去重硬件副作用，不能只依赖 reducer。

**持久化与断点恢复**：`FrequencyCalibrationRuntime` 提供 `save_run(dir)` 和
`load_run(dir, qubit, executor=...)` 方法，将运行状态写入五个文件：

| 文件 | 内容 |
|---|---|
| `config.json` | 完整协议配置 |
| `commands.jsonl` | 每周期 journal：命令/事件、状态、测量、成本、diagnostics 和原因 |
| `transitions.jsonl` | 每次状态转移（from / to / reason） |
| `checkpoint.json` | 状态机全量快照（可从此恢复） |
| `result.json` | 最终结果摘要 |

checkpoint 会恢复配置、budget、tracker、Verify/Lock 计数、retry 状态、转移历史和
pending command。恢复语义为 **at-least-once**：尚未确认的硬件动作可能使用原 ID 再次
提交；这不构成跨硬件或进程故障的严格 exactly-once 保证。

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
    Delta_val=2 * np.pi * 50e-3,          # 50 MHz  — 探针失谐有效窗口
    confidence_multiplier=1.0,             # 确定性仿真；实验须预先指定 z_(1-beta)
    linear_range=2 * np.pi * 80e-3,       # 后端可信的局部范围
    min_confidence=None,                   # 仅在后端真实报告置信度时设置

    # -- 验证策略 --
    N_verify=2,                            # 需要连续 2 次通过才算验证成功
    max_verify_attempts_per_episode=5,
    max_verify_shots=10_000,
    verify_track_max_residual=2 * np.pi * 20e-3,

    # -- Lock 监视器迟滞 --
    epsilon_mon_clear=2 * np.pi * 5e-3,   # 5 MHz   — 监视器认为"干净"
    epsilon_mon_suspect=2 * np.pi * 10e-3,# 10 MHz  — 监视器认为"可疑"
    Delta_mon_reacquire=2 * np.pi * 30e-3,# 30 MHz  — 大跳变，直接重捕获
    N_mon_suspect=3,                       # 连续可疑 3 次 → Verify
    monitor_interval=1.0,                  # 低成本监测间隔（秒）
    audit_interval=60.0,                   # Ramsey 审计间隔（秒）
    require_periodic_audit=True,

    # -- 预算 --
    max_commands=200,                      # 总命令数上限
    max_reacquire_attempts=5,             # 重捕获次数上限
    max_wall_time=3600.0,                 # 1 小时墙上时间
    stop_after_lock_cycles=20,             # 有界示例；0 表示长期运行

    # -- Track 控制器参数 --
    damping=0.8,                           # 阻尼因子（<1 抑制过冲）
    first_bias_step=0.01,                  # 首次探测步长 (Φ₀)
    expected_sensitivity_sign=-1,          # 当前单调分支 χ_Φ；未知时可设 None
    max_bias_step=0.02,                    # 每步最大偏置变化 (Φ₀)
    require_monitor_local_validity=True,   # Lock 也必须报告可信探针失谐
)

# ── 运行 ──────────────────────────────────────────────────────────────
runtime = FrequencyCalibrationRuntime(
    qubit=qubit,
    f_target=f_target,
    config=config,
    checkpoint_directory="calibration_run_001",
)
result = runtime.run()

# ── 结果解读 ──────────────────────────────────────────────────────────
print(f"终态: {result['state']}")               # 本有界运行中为 "safe_stop"
print(f"运行状态: {result['run_status']}")       # SafeStop 后为 "completed"
print(f"安全保持: {result['safe_hold_confirmed']}")
print(f"保持目标: {result['hold_target_met']} "
      f"({result['hold_passes']}/{result['hold_samples']})")
print(f"是否中断: {result['interrupted']} {result['interrupt']}")
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
# runtime2 = FrequencyCalibrationRuntime.load_run(
#     "calibration_run_001", qubit=qubit, executor=idempotent_executor
# )
# runtime2.run()  # 状态机从最后一个 checkpoint 恢复
```

这段代码展示了 V2 接口的核心流程：配置阈值 → 启动运行 → 查看轨迹 → 持久化保存。
`transition_log` 中每条记录包含状态转移的起止状态和触发原因码（如 `TARGET_CANDIDATE`、
`DRIFT_SUSPECTED`、`REACQUIRE_LIMIT` 等），方便事后审计——为什么进入了 Reacquire？
连续几次 suspect 触发了 Verify？journal 每行还包含状态前后、偏置、drive、频率、残差、
不确定度、valid/ambiguity、shots、耗时、diagnostics 和转移原因。

进入 Lock 会把 lifecycle 标为 `CALIBRATED`，但不会结束 `run()`。Lock 会继续监测，
计划 Ramsey 审计也必须经过 Lock → Verify。本有界示例最终由
`stop_after_lock_cycles` 触发 SafeStop，`safe_hold_confirmed` 记录 SafeHold 是否确认。

UI/API 需要取消时，另一个线程只设置共享信号：

```python
runtime.request_cancel(reason="operator_stop", source="ui")
```

runtime 会在命令边界和 monitor 等待期间检查信号，记录 `CancelRequested`，写入配置的
紧急 checkpoint，再执行 SafeHold。`KeyboardInterrupt` 默认也转换到同一路径。这是
协作式取消：如果 `executor.execute()` 已经开始，当前同步测量必须返回后才能处理中断。

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
- **V2 状态机结果**：`runtime.run()` 返回 `state`（有界完成时为 `"safe_stop"`）、
  `run_status`（确认 SafeStop 后为 `"completed"`；`CALIBRATED` 只是首次进入 Lock 的
  非终止里程碑）、`safe_hold_confirmed`、`hold_target_met`（最近一次 Lock monitor
  是否满足 `epsilon_hold`）、`hold_passes`/`hold_samples`、`f_final`（最终
  频率估计）、`candidate_bias`（锁定时的偏置值）、`transition_log`（完整转移轨迹，
  每条含 `from`/`to`/`reason`/`version`）、`n_commands`（消耗的命令数）、`elapsed`
  （墙上时间）和逐周期 `journal`。`transition_log` 与 `journal` 共同构成事后审计轨迹。
- 各标定类的完整字段与方法选项见 {doc}`../building_blocks/calibration`；V2 状态机
  的阈值与配置见 {doc}`../building_blocks/workflows`。
