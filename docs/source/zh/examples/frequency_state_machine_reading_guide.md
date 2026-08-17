---
orphan: true
---

# 频率标定状态机源码导读

这篇导读不把状态机拆成类名清单，而是沿一次闭环运行逐层打开六个黑箱。每一遍只回答一个新的问题，并把上一遍暂时接受的黑箱变成这一遍的主角。

```text
第一遍：状态机用什么语言描述协议？
             ↓
第二遍：这些语言怎样形成正常状态路径？
             ↓
第三遍：一条命令怎样变成物理测量证据？
             ↓
第四遍：多轮 Track 怎样形成数值闭环？
             ↓
第五遍：谁把决策、控制和测量组织成持续运行？
             ↓
第六遍：失效、中断和进程恢复时，系统凭什么仍可信？
```

贯穿全文的最小循环是：

```text
当前状态
   ↓ next_command()
Command
   ↓ Runtime / Executor 执行
Event
   ↓ handle(event)
下一状态与 ReasonCode
```

```{figure} frequency_state_machine_datapath.svg
:alt: 频率标定状态机的命令—事件数据通路
:width: 100%
:align: center

频率标定状态机的 CPU 式数据通路。黑线传递命令或测量证据，蓝线传递控制字，红线表示异常路径，绿色虚线表示持久化与恢复。该图描述架构职责，不要求每个框对应一个 Python 类；六遍阅读只是在同一拓扑上依次激活不同路径。
```

## 如何沿这张图完成六遍阅读

第一遍先建立接口边界。只看 **Register File**、**Control Unit**、**CMD REG**、**EVT REG** 和 **Commit Unit**：寄存器保存已提交事实，控制器发出命令，执行侧返回事件，提交单元负责把接受的结果写回。此时不进入测量、割线或异常细节。

第二遍沿正常提交环走一圈：**Register File → Control Unit → Operand MUX → CMD REG → Issue Unit → Measurement Unit → EVT REG → Guard Comparator → Control Unit → Commit Unit → write-back**。Acquire、Track、Verify 和 Lock 是 `SR` 中的状态值，而不是四套串联执行器。

第三遍只打开 **CMD REG → Issue Unit → Measurement Unit → EVT REG**。从已锁存且获准执行的命令出发，追踪不同 role 怎样选择 Ramsey 或局部 probe，并怎样把频率、不确定度、实际偏置、实际驱动、shots 和耗时封装成 Event。

第四遍转向 `TR/ER → Secant ALU → Operand MUX` 反馈环。已接受的 Track 事件先提交到 Track/Estimate 寄存器，割线单元再从跨轮记忆计算下一轮 bias/drive proposal。这样可以同时看清“证据提交”和“下一步控制量计算”是两件事。

第五遍把主环放回持续运行环境。重点看 `BR`、**Issue Unit** 和 **Commit Unit** 如何处理预算、中断、调度与 identity，再沿绿色虚线理解 **Journal Memory** 如何保存历史、成本和 checkpoint；底部 write-back bus 表示一次接受事件的原子提交边界。

第六遍最后沿 **Guard Comparator → Control Unit / Exception Unit** 检查异常语义。事件无效、硬件失败、预算耗尽或用户取消怎样成为 fault class，何时 retry、Reacquire 或 SafeStop，以及 checkpoint 怎样恢复寄存器上下文。随后再用测试验证这些路径是否具有代码证据。

第一遍先认识这几个词；第二遍让它们跑起来；第三、四遍分别打开 `Executor` 和 Track 控制器；第五遍再把整个循环放回 Runtime；第六遍最后检查正常路径之外的恢复语义和测试证据。

## 阅读地图

| 遍次 | 这一遍打开的黑箱 | 核心问题 | 读完后的接力点 |
|---|---|---|---|
| 第一遍 | 状态机词汇与对象 | 状态、命令、事件和原因分别表达什么？ | 已能读懂一条协议记录，但还不知道状态怎样前进 |
| 第二遍 | 正常状态路径 | `Acquire → Track → Verify → Lock` 怎样逐轮发生？ | 已能追踪转换，但暂不解释测量数据怎样产生 |
| 第三遍 | 测量 Backend | Ramsey 与短脉冲 probe 怎样生成事件？ | 已理解单轮证据，但 Track 仍缺少跨轮控制记忆 |
| 第四遍 | Track 控制器 | 割线、偏置和 drive tracking 怎样跨轮协同？ | 已形成科学闭环，但还缺持续调度与资源管理 |
| 第五遍 | Runtime | 预算、中断、日志和 checkpoint 怎样包住闭环？ | 已理解正常运行生命周期，最后转向异常路径 |
| 第六遍 | 恢复与测试 | 失败怎样分类、恢复语义由什么测试支撑？ | 能从日志或测试反向定位到协议、控制、测量或运行层 |

建议顺序阅读。若已有状态机基础，可以快速浏览第一遍，但不要跳过 `command_id`、`state_version`、目标残差与 probe detuning 的区分；这些概念会贯穿后五遍。

## 第一遍：认识状态机的语言

```{figure} frequency_state_machine_pass1.svg
:alt: 第一遍高亮寄存器、控制器及命令事件边界
:width: 100%
:align: center

第一遍的激活路径：先识别寄存器中的协议词汇，再认识 Control Unit、Command、Event 与 Commit 的责任边界；其余功能单元暂时作为灰色黑箱。
```

在这张抽象图里，`CalibrationState` 和当前科学估计属于 **Register File**；`next_command()` / `handle(event)` 属于 **Control Unit** 的协议职责；Command 与 Event 分别跨过 `CMD REG` 和 `EVT REG` 边界；`ReasonCode` 最终进入 `LR · Last result`。这些是阅读类和字段时的定位坐标，而不是新的运行时类。

第一遍我们暂时不追踪具体状态转换，也不看 QuTiP、割线算法或 Runtime，只读 frequency_state_machine.py 的“词汇表”。

读完第一遍，你应该能够回答：

> 状态机现在处于什么状态、希望外界做什么、外界返回了什么、为什么发生转移，以及整个任务是否仍在运行。

---

### 1.1 阅读边界：先学语言，不追路径

#### 一、先建立整体心智模型

这不是一个主动操作硬件的对象，而是一个纯决策器。可以把它理解成一位控制室调度员：

```text
状态 State
   ↓ 决定
命令 Command
   ↓ 交给外部执行
事件 Event
   ↓ 返回状态机
转换原因 ReasonCode
   ↓
新状态 State
```

状态机自己不知道怎样运行 Ramsey，也不直接设置磁通偏置。它只会说：

```text
“请在这个偏置上做一次 Acquire”
“请在这个候选偏置上做一次 Track”
“请冻结偏置并做 Verify”
```

外部测量完成后，再回答它：

```text
“测量成功，频率是……”
“测量成功，但结果有歧义”
“硬件超时”
“用户请求取消”
```

这就是文件顶部示例的含义：

```python
command = machine.next_command()
event = executor.execute(command)
machine.handle(event)
```

第一遍要认识的就是这套交流语言。

---

### 1.2 两条坐标轴：协议阶段与运行生命周期

#### 二、State：协议现在处在哪一章

先看 FrequencyState：

```python
class FrequencyState(str, Enum):
    ACQUIRE = "acquire"
    TRACK = "track"
    VERIFY = "verify"
    LOCK = "lock"
    REACQUIRE = "reacquire"
    SAFE_STOP = "safe_stop"
```

这里的状态不是“动作”，而是协议当前所处的阶段。

##### Acquire

含义是：

> 当前还没有足够可信的绝对频率参考，需要进行宽范围捕获。

它通常通过 Ramsey 获得：

$$
\hat f_{\mathrm{acq}},\qquad \sigma_{\mathrm{acq}}.
$$

Acquire 本身不负责调偏置，只负责建立种子。

##### Track

含义是：

> 已经知道频率大致在哪里，可以使用局部、低成本的测量和控制算法逐步靠近目标。

Track 是唯一允许提出新磁通偏置的状态。它同时做：

- 短脉冲局部测频；
- 检查探针失谐是否有效；
- 估计频率—磁通灵敏度；
- 提出下一步偏置；
- 更新下一次测量的驱动频率。

##### Verify

含义是：

> Track 认为已经接近目标，但这个结论还没有得到独立测量确认。

进入 Verify 后冻结候选偏置，用独立 Ramsey 测量判断是否真正达到最终容差。Verify 只能测量，不能调整偏置。

##### Lock

含义是：

> 候选偏置已经通过独立 Verify，现在进入长期保持和漂移监测阶段。

Lock 不直接纠偏。检测到疑似小漂移时进入 Verify，大跳变或参考失效时进入 Reacquire。

##### Reacquire

含义是：

> 现有局部频率参考已经不可信，不能继续使用局部反演，需要重新进行宽范围捕获。

Reacquire 与 Acquire 使用相似的宽范围测量，但语义不同：

- Acquire 是首次建立种子；
- Reacquire 是运行过程中丢失种子后的恢复。

##### SafeStop

含义是：

> 本次协议不能或不应继续执行，需要转入安全偏置并停止科学测量。

典型原因包括：

- 用户取消；
- 硬件联锁；
- 预算耗尽；
- 实际偏置与命令不一致；
- 多次硬件失败。

因此，主状态拓扑可以先记成：

```text
Acquire → Track → Verify → Lock
            ↑        |
            └────────┘
              小修正

Track / Verify / Lock → Reacquire
所有状态              → SafeStop
```

---

#### 三、RunStatus：整个任务处于什么生命周期

接着看 RunStatus。

它和 `FrequencyState` 是两个不同维度。

例如：

```text
FrequencyState = LOCK
RunStatus      = CALIBRATED
```

表示协议位于 Lock，并且已经至少完成过一次独立验证。

##### 为什么不能只用 State

`FrequencyState` 回答：

> 协议下一步应该做什么？

`RunStatus` 回答：

> 整个运行任务处于什么生命周期？

具体包括：

- `READY`：对象已经创建，但还没启动；
- `RUNNING`：状态机正在执行；
- `CALIBRATED`：第一次成功完成 `Verify → Lock`；
- `COMPLETED`：有界仿真按策略正常结束；
- `SAFE_STOPPED`：因取消或联锁安全停止；
- `FAILED`：尚未成功标定便因不可恢复错误退出。

这里有一个容易忽略的点：进入 Lock 后，状态机仍然可以继续长期运行。因此：

```text
CALIBRATED ≠ 程序已经终止
```

它只表示“已经证明过一次标定成功”。

例如：

```text
state      = Lock
run_status = Calibrated
```

之后仍然会不断执行 monitor。

若配置要求测试只运行有限个 Lock 周期，则结束时可能变成：

```text
state      = SafeStop 或终止边界
run_status = Completed
```

---

### 1.3 决策政策：原因码与配置共同定义协议

#### 四、ReasonCode：不是只记录去了哪里，还记录为什么

看 ReasonCode 时，不建议逐项死记。应按“原因属于哪一类”理解。

##### 正常前进原因

```python
SEED_RELIABLE
CANDIDATE_REACHED
VERIFY_PASSED
```

分别表示：

```text
Acquire → Track：种子可靠
Track   → Verify：候选条件成立
Verify  → Lock：独立验证通过
```

##### 仍需继续修正

```python
CORRECTION_REQUIRED
```

可用于：

```text
Track  → Track
Verify → Track
```

虽然两个转换目标不同，但共同含义是“尚未达到最终目标，仍需局部控制”。

##### 局部测量失效

```python
LOCAL_RANGE_LOST
INVERSE_FAILED
BRANCH_DISCONTINUITY
SENSITIVITY_INVALID
CONFIDENCE_INSUFFICIENT
```

它们分别表达：

- 探针失谐超出有效区；
- 局部反演无法求解；
- 三阶反演跳到错误分支；
- 割线灵敏度异常；
- 测量置信度不足。

这些原因通常导致 Reacquire，因为继续使用局部控制已经不安全。

##### Lock 漂移原因

```python
MONITOR_CLEAR
MONITOR_SUSPECT
DRIFT_SUSPECTED
LARGE_FREQUENCY_JUMP
REFERENCE_LOST
AUDIT_DUE
```

其中：

- `MONITOR_CLEAR`：监测正常，保持 Lock；
- `MONITOR_SUSPECT`：本次可疑，但尚未达到连续次数；
- `DRIFT_SUSPECTED`：连续可疑，进入 Verify；
- `LARGE_FREQUENCY_JUMP`：大跳变，直接 Reacquire；
- `REFERENCE_LOST`：参考丢失，直接 Reacquire；
- `AUDIT_DUE`：周期性独立复核到期，进入 Verify。

这里 `AUDIT_DUE` 仍是代码中的旧命名。我们此前讨论过，更合适的未来命名是 `PERIODIC_VERIFY_DUE`，但当前代码尚未完成这个兼容重命名。

##### 运行安全原因

```python
TECHNICAL_TIMEOUT
INTERLOCK
BUDGET_EXHAUSTED
REACQUIRE_LIMIT
CANCELLED
```

这些原因描述的不是物理标定结果，而是运行系统无法继续。

##### 为什么必须使用 Enum

如果只写：

```python
reason = "measurement bad"
```

后续很难统计究竟是越界、低置信度、硬件失败还是参考丢失。

使用稳定的 `ReasonCode` 后，journal 可以可靠分析：

```text
有多少次 Track 因 LOCAL_RANGE_LOST 失败？
有多少次 Lock 因 DRIFT_SUSPECTED 进入 Verify？
是否经常因为 SENSITIVITY_INVALID 触发 Reacquire？
```

因此 `ReasonCode` 既服务于控制逻辑，也服务于复现实验和故障分析。

---

#### 五、Config：这台状态机采用什么政策

FrequencyCalibrationConfig 不是简单的参数集合，而是“一次运行的协议政策”。

同一套状态机代码，通过不同配置可以表达：

- 确定性仿真；
- 有限-shot 仿真；
- 短时间测试；
- 长期硬件运行；
- 更保守或更激进的控制策略。

##### 第一组：三个目标容差和一个局部范围

最核心的是：

```python
epsilon_enter
epsilon_final
epsilon_hold
Delta_val
```

满足：

$$
\epsilon_{\mathrm{hold}}
<
\epsilon_{\mathrm{final}}
<
\epsilon_{\mathrm{enter}}
<
\Delta_{\mathrm{val}}.
$$

它们分别回答四个不同问题。

###### `epsilon_enter`

Track 判断：

> 当前结果是否已经足够接近，可以停止调节并交给 Verify？

对应候选进入条件：

$$
|r|+z\sigma_r\leq\epsilon_{\mathrm{enter}}.
$$

它不是最终成功条件，只是进入 Verify 的门槛。

###### `epsilon_final`

Verify 判断：

> 独立 Ramsey 是否确认候选工作点达到最终容差？

$$
|r_{\mathrm{ver}}|+z\sigma_{\mathrm{ver}}
\leq\epsilon_{\mathrm{final}}.
$$

###### `epsilon_hold`

Lock 记录：

> 长期监测结果是否满足最严格的物理保持目标？

它用于长期性能记录，不应直接与 monitor 路由阈值混为一谈。

###### `Delta_val`

Track 判断：

> 当前短脉冲 probe 的实际失谐是否还处于经过验证的局部范围？

$$
|\widehat\Delta|+z\sigma_\Delta
\leq\Delta_{\mathrm{val}}.
$$

这里最重要的是区分：

$$
r=f-f_{\mathrm{target}},
\qquad
\Delta=f-f_{\mathrm{drive}}.
$$

前三个 `epsilon` 主要约束目标残差 \(r\)，而 `Delta_val` 约束探针失谐 \(\Delta\)。

##### `confidence_multiplier`

```python
confidence_multiplier = 1.0
```

对应公式中的 \(z\)：

$$
U=|\hat x|+z\sigma_x.
$$

代码默认 `1.0` 是为了兼容旧行为，并不意味着所有实验都应该采用一倍标准差。真机配置需要根据覆盖率和多重检测政策确定 \(z\)。

##### Acquire 配置

```python
sigma_acquire_max
max_acq_retries
acquisition_bias
```

分别表示：

- Acquire 允许的最大频率不确定度；
- 宽范围捕获最多重试多少次；
- 在哪个磁通偏置执行 Acquire。

##### Track 配置

```python
damping
first_bias_step
max_bias_step
guard_margin
linear_range
S_min
S_max
```

这些参数约束局部控制器：

- `damping`：割线步长阻尼；
- `first_bias_step`：没有割线信息时的首次探索步；
- `max_bias_step`：单次最大偏置变化；
- `guard_margin`：局部有效范围的安全余量；
- `linear_range`：额外的短脉冲捕获区约束；
- `S_min/S_max`：允许的频率—磁通灵敏度范围。

第一遍不需要理解割线公式，只要知道：

> 这些参数限制 Track 可以迈多大步，以及什么时候必须停止相信局部控制。

##### Verify 配置

```python
N_verify
max_verify_attempts_per_episode
max_verify_shots
bias_freeze_tolerance
verify_track_max_residual
```

其中：

- `N_verify`：需要连续通过几次独立 Verify；
- `max_verify_attempts_per_episode`：一次 Verify 阶段最多尝试多少次；
- `max_verify_shots`：一次 Verify episode 的最大 shots；
- `bias_freeze_tolerance`：实际偏置偏离冻结偏置多少算联锁错误；
- `verify_track_max_residual`：Verify 未通过时，偏差多小才允许返回 Track。

最后一个参数避免把所有 Verify 失败都送回 Track：

```text
小偏差、结果可靠 → Track
大偏差或不可信   → Reacquire
```

##### Lock 配置

```python
epsilon_mon_clear
epsilon_mon_suspect
Delta_mon_reacquire
N_mon_suspect
monitor_interval
audit_interval
```

形成三段路由：

$$
U_{\mathrm{mon}}\leq\epsilon_{\mathrm{mon,clear}}
$$

保持 Lock；

$$
\epsilon_{\mathrm{mon,clear}}
<
U_{\mathrm{mon}}
\leq
\epsilon_{\mathrm{mon,suspect}}
$$

累计可疑次数；

$$
U_{\mathrm{mon}}\geq\Delta_{\mathrm{mon,reacquire}}
$$

直接 Reacquire。

`monitor_interval` 控制低成本监测间隔，`audit_interval` 控制周期性独立 Ramsey 复核间隔。

##### 预算配置

```python
max_commands
max_wall_time
max_shots
max_solver_calls
max_reacquire_attempts
stop_after_lock_cycles
```

这些参数保证状态机不会无限消耗资源。

尤其要理解：

```python
stop_after_lock_cycles = 0
```

表示长期运行，不自动结束；测试中通常设置正整数，避免无限循环。

##### `__post_init__()` 为什么重要

配置类后半部分的 `__post_init__()` 不是业务流程，而是“协议合法性证明”。

例如它拒绝：

```text
epsilon_final 大于 epsilon_enter
负数的 monitor interval
N_verify 小于 1
max_verify_attempts 少于 N_verify
linear_range 小于 guard_margin
强制周期复核却把 audit_interval 设为 0
```

所以读配置时要形成这样的认识：

> Config 不仅储存参数，也防止用户构造逻辑上自相矛盾的状态机。

---

### 1.4 交互契约：Command、Event 与身份字段

#### 六、Command：状态机向外部提出的请求

命令定义从 AcquireFrequency 开始。

命令只描述“要做什么”，不包含执行结果。

##### `AcquireFrequency`

```python
AcquireFrequency(
    role="acquire",
    bias=...,
)
```

含义是：

> 请在指定偏置执行一次宽范围频率捕获。

`role` 可以是：

- `"acquire"`：首次捕获；
- `"reacquire"`：恢复捕获。

因此 Acquire 和 Reacquire 可以复用同一个测量接口。

##### `TrackFrequency`

```python
TrackFrequency(
    candidate_bias=...,
    predicted_drive=...,
    predicted_detuning=...,
    predicted_detuning_uncertainty=...,
    prediction_guard_complete=...,
    tracker_spec=...,
)
```

这是信息量最大的命令。

- `candidate_bias`：这一轮实际要测量的新偏置；
- `predicted_drive`：drive tracking 预测的驱动频率；
- `predicted_detuning`：执行测量前预测的探针失谐；
- `predicted_detuning_uncertainty`：预测失谐的不确定度；
- `prediction_guard_complete`：测前范围守卫是否具有完整证据；
- `tracker_spec`：割线灵敏度、步长等控制器诊断信息。

注意 `candidate_bias` 不是“已经验证的候选点”，而是 Track 下一步要尝试的工作点。这里 `candidate` 的命名容易造成误解。

##### `VerifyFrequency`

```python
VerifyFrequency(
    frozen_bias=...,
    drive=f_target,
)
```

含义很明确：

> 保持该偏置不变，使用目标频率作为参考执行独立 Ramsey 验证。

字段叫 `frozen_bias` 而不是 `candidate_bias`，是在接口层明确表达 Verify 不允许修改偏置。

##### `MonitorFrequency`

```python
MonitorFrequency(
    locked_bias=...,
)
```

表示：

> 在已经验证的锁定偏置上执行一次低成本漂移监测。

同样，`locked_bias` 表明这个命令没有“提出新偏置”的权力。

##### `SafeHold`

```python
SafeHold(
    bias=...,
    reason=...,
)
```

表示：

> 将设备置于安全偏置，并停止后续科学测量。

##### `command_id`

所有命令都有唯一的 `command_id`：

```python
command_id = uuid4()
```

它用于把返回事件和原始命令配对。

例如命令 A 已经超时，系统重试后发出命令 B。如果命令 A 的迟到结果随后返回，状态机可以通过 `command_id` 识别它是旧结果并忽略。

因此 `command_id` 是幂等、重试和 checkpoint 恢复的基础。

---

#### 七、Event：外界向状态机提交的证据

事件定义从 MeasurementSucceeded 开始。

##### `MeasurementSucceeded`

它不只是“测量没有报错”，而是一份结构化测量证据：

```python
MeasurementSucceeded(
    command_id=...,
    frequency=...,
    uncertainty=...,
    valid=...,
    ambiguous=...,
    method=...,
    shots=...,
    elapsed_time=...,
    applied_bias=...,
    applied_drive=...,
    probe_detuning=...,
    probe_detuning_uncertainty=...,
    diagnostics=...,
)
```

这些字段可以分成四层。

###### 频率估计

```python
frequency
uncertainty
```

即：

$$
\hat f,\qquad \sigma_f.
$$

状态机可以据此计算目标残差：

$$
r=\hat f-f_{\mathrm{target}}.
$$

###### 估计质量

```python
valid
ambiguous
```

它们不是同一个含义：

- `valid=False`：拟合、反演或质量检查失败；
- `ambiguous=True`：存在多个不能可靠区分的可能解。

即使频率字段存在，只要结果有歧义，也不能用于控制。

###### 实际执行证据

```python
applied_bias
applied_drive
```

它们表示硬件实际执行了什么，而不是命令要求了什么。

Verify 可以检查：

$$
|\Phi_{\mathrm{applied}}-\Phi_{\mathrm{frozen}}|
\leq
\delta\Phi_{\mathrm{freeze}}.
$$

Track 也可以检查实际偏置、驱动是否与命令一致。若不一致，说明控制系统与测量证据脱节，应触发 interlock。

###### 探针局部有效性

```python
probe_detuning
probe_detuning_uncertainty
```

对应：

$$
\widehat\Delta
=
\hat f-f_{\mathrm{drive}},
\qquad
\sigma_\Delta.
$$

它们只用于判断局部短脉冲 probe 是否仍然可信：

$$
|\widehat\Delta|+z\sigma_\Delta
\leq
\Delta_{\mathrm{val}}.
$$

这是不能用目标残差代替的字段。

###### 成本与追溯

```python
shots
elapsed_time
method
diagnostics
```

用于记录：

- 本次测量用了多少 shots；
- 花了多长时间；
- 使用 Ramsey 还是 transient；
- solver calls、灵敏度、置信度和原始数据引用。

`diagnostics` 是扩展区；状态转移所必需的核心物理量应尽量使用明确字段，而不是全部塞进 `diagnostics`。

##### `MeasurementTechnicalFailure`

MeasurementTechnicalFailure 表示：

> 测量没有产生可用于判断的科学结果。

例如：

- 硬件超时；
- 仪器通信失败；
- solver 抛出异常；
- 数据文件损坏。

它与 `valid=False` 的区别是：

```text
TechnicalFailure：测量过程没有正常完成
valid=False：测量完成了，但科学结果不可信
```

这种区别会影响重试策略。技术故障通常可以原状态重试；科学无效往往需要 Reacquire 或切换方法。

##### `MeasurementRejected`

MeasurementRejected 表示：

> 测量在技术上完成，但在进入状态机之前已被质量守卫拒绝。

例如 Backend 或 Runtime 已经判断：

```text
LOCAL_RANGE_LOST
SENSITIVITY_INVALID
CONFIDENCE_INSUFFICIENT
```

它与 `MeasurementSucceeded(valid=False)` 有部分语义重叠，但 `MeasurementRejected` 能直接携带稳定的 `ReasonCode`。

##### 自发事件

并非所有事件都对应一次测量命令。

###### `TimerElapsed`

表示定时器到期，例如周期性独立复核到期。

###### `InterlockTriggered`

表示外部安全联锁触发，不需要等待当前状态的普通测量逻辑。

###### `BudgetExhausted`

表示 command、shots、wall time 或 solver-call 预算耗尽。

###### `CancelRequested`

表示用户或上层程序请求停止。

这些事件之所以叫“自发事件”，是因为它们可以在没有匹配科学命令的情况下改变状态机。

---

#### 八、身份字段：如何避免旧事件污染新状态

事件里还有两个容易被忽略的字段：

```python
command_id
state_version
```

##### `command_id`

回答：

> 这个结果属于哪一条命令？

只接受当前 pending command 的返回事件。

##### `state_version`

回答：

> 这个结果属于状态机的哪一个版本？

每次真正发生状态转换时，`state_version` 会递增。

例如：

```text
version 3：Track
version 4：Verify
```

如果一个属于 Track version 3 的迟到事件在系统已经进入 Verify version 4 后返回，即使内容看起来正常，也不能再用于当前决策。

二者共同防止异步或重试环境中的竞态：

```text
command_id 解决“是哪次操作”
state_version 解决“属于哪个状态时期”
```

---

### 1.5 从词汇表回到一个真实对象

#### 九、初始化对象时，状态机记住了什么

最后看 FrequencyStateMachine.__init__()。

创建后初始状态是：

```python
_state = ACQUIRE
_run_status = READY
```

这意味着：

> 协议将从 Acquire 开始，但尚未正式运行。

其他内部字段可以按用途理解。

##### 当前最佳物理认识

```python
_f_hat
_uncertainty
_candidate_bias
```

表示状态机当前认为：

- 频率是多少；
- 该估计的不确定度是多少；
- 当前候选偏置是什么。

##### 连续判定计数

```python
_verify_streak
_monitor_streak
```

分别记录：

- Verify 已连续通过多少次；
- Lock 已连续出现多少次可疑 monitor。

##### Lock 长期记录

```python
_lock_entry_time
_hold_target_met
_hold_samples
_hold_passes
```

用于记录进入 Lock 的时间，以及多少个监测样本真正满足 `epsilon_hold`。

##### 运行资源

```python
_budget
```

记录 shots、solver calls、命令数量、Verify 尝试次数和 Reacquire 次数等。

##### 可复现性

```python
_transition_log
_last_event
_pending_command
_tracker_snapshot
```

分别用于：

- 保存状态转换历史；
- 保存最后一次事件；
- 保存已发出但尚未完成的命令；
- 保存 Track 控制器的可序列化状态。

##### `start()`

start() 只做生命周期启动：

```text
READY → RUNNING
启动 wall-time 预算
记录“运行开始”
```

它不会直接执行 Acquire。真正生成首条 Acquire 命令的是后面的 `next_command()`，那属于第二遍阅读内容。

---

### 1.6 第一遍收束：能解释一条记录

#### 第一遍读完后应该形成的完整图景

现在可以把所有概念放进一句话：

> `FrequencyStateMachine` 在某个 `FrequencyState` 中，根据 `FrequencyCalibrationConfig` 生成带唯一 `command_id` 的 `Command`；外部执行后返回携带测量证据的 `Event`；状态机用 `ReasonCode` 记录为什么进入下一个状态，同时用 `RunStatus` 表示整个任务是否仍在运行、已经标定、正常结束或安全停止。

压缩成代码语言就是：

```text
FrequencyState       当前协议阶段
RunStatus            整个任务生命周期
Config               运行政策
Command              状态机对外部的请求
Event                外部返回的证据或信号
ReasonCode           决策原因
command_id/version   防止迟到和重复事件
```

第一遍建议先停在这里，不要立即钻进 `_handle_track()`。你可以自己做一个很小的阅读练习：手工构造一条 `TrackFrequency` 和对应的 `MeasurementSucceeded`，逐字段解释“命令要求了什么”和“事件证明了什么”。能准确区分 `predicted_drive`、`applied_drive`、`frequency`、`probe_detuning` 与目标残差，就说明第一遍已经真正读懂了。


## 第二遍：追踪 Acquire → Track → Verify → Lock

```{figure} frequency_state_machine_pass2.svg
:alt: 第二遍高亮正常状态提交环
:width: 100%
:align: center

第二遍的激活路径：从已提交上下文生成一条 Command，经测量得到 Event，再由 Guard、Control 和 Commit 决定下一状态并原子写回。
```

这条环路就是软件中的一轮协议推进。`SR` 保存 Acquire、Track、Verify 或 Lock；Control Unit 根据当前 `SR` 选择命令；Guard Comparator 把证据变成条件码；Control Unit 决定 `next state + reason`；Commit Unit 才拥有写回 `SR/ER/LR` 的权限。

第二遍我们追踪一次完整的正常运行：

$$
\text{Acquire}\rightarrow\text{Track}\rightarrow\text{Verify}\rightarrow\text{Lock}.
$$

这一遍重点不是测量怎样实现，而是观察状态机的两扇门：

- `next_command()`：状态机向外输出命令；
- `handle(event)`：状态机接收结果并转移。

可以把每轮运行记成：

```text
当前状态
  ↓ next_command()
pending command
  ↓ 外部执行
event
  ↓ handle()
新状态
```

---

### 2.1 先看一轮协议怎样闭合

#### 1. 起点：`start()` 只启动任务

创建状态机后：

```python
machine = FrequencyStateMachine(config, f_target)
```

内部状态为：

```text
state      = Acquire
run_status = Ready
```

此时还不能直接调用 `next_command()`，需要先执行：

```python
machine.start()
```

start() 做三件事：

```text
RunStatus: Ready → Running
启动 wall-time 预算计时
记录一条“run started”转换日志
```

注意，`start()` 不执行 Ramsey，也不生成测量结果。它只是宣布：

> 这次运行正式开始，当前协议阶段是 Acquire。

---

#### 2. `next_command()`：状态机的输出侧

正常路径首先进入 next_command()。

它的完整逻辑可以压缩成：

```python
def next_command():
    检查任务是否已经结束
    检查预算是否耗尽
    如果已有 pending command:
        原样返回
    根据当前状态生成命令
    保存为 pending command
    消耗一次 command 预算
    返回命令
```

##### 2.1 为什么命令要保存为 pending

假设第一次调用：

```python
cmd1 = machine.next_command()
```

状态机产生：

```python
AcquireFrequency(
    command_id="abc-123",
    role="acquire",
    bias=0.0,
)
```

随后再次调用：

```python
cmd2 = machine.next_command()
```

在没有收到结果前：

```python
cmd2 is cmd1
cmd2.command_id == "abc-123"
```

状态机不会生成另一条新命令。

这是 `next_command()` 的幂等性：

> 一条命令在收到对应事件之前，始终是同一条 pending command。

它解决了 checkpoint、重试和通信中断问题。即使上层不确定命令是否已经发出，再调用一次也不会创建第二个科学操作。

##### 2.2 `_dispatch_command()` 负责状态到命令的映射

真正决定命令类型的是 _dispatch_command()：

```text
Acquire   → AcquireFrequency
Track     → TrackFrequency
Verify    → VerifyFrequency
Lock      → MonitorFrequency 或 VerifyFrequency
Reacquire → AcquireFrequency(role="reacquire")
SafeStop  → SafeHold
```

`next_command()` 负责通用管理，`_dispatch_command()` 负责协议语义。

这两个函数的分工是：

```text
next_command()
    管生命周期、预算、幂等和 pending

_dispatch_command()
    管当前状态应该请求什么操作
```

---

#### 3. 第一轮：Acquire 发出宽范围捕获命令

状态为 Acquire 时，`_dispatch_command()` 返回：

```python
AcquireFrequency(
    role="acquire",
    bias=config.acquisition_bias,
)
```

默认：

```python
acquisition_bias = 0.0
```

此时状态机内部变为：

```text
state              = Acquire
pending_command     = AcquireFrequency(...)
pending_command_id  = 当前 UUID
```

注意状态仍是 Acquire。发出命令本身不会引起状态转移，因为测量结果尚未回来。

完整顺序是：

```text
Acquire
   ↓ 生成命令
Acquire + pending AcquireFrequency
   ↓ 收到测量结果
Track 或 Verify
```

---

#### 4. `handle()`：所有事件的统一入口

测量后端执行命令并返回，例如：

```python
MeasurementSucceeded(
    command_id=cmd.command_id,
    frequency=f_acq,
    uncertainty=sigma_acq,
    valid=True,
    ambiguous=False,
    applied_bias=0.0,
)
```

然后调用：

```python
machine.handle(event)
```

统一入口是 handle()。

它先不关心当前是 Acquire 还是 Track，而是做通用检查。

##### 4.1 先处理自发事件

`handle()` 首先检查：

```text
InterlockTriggered
CancelRequested
BudgetExhausted
TimerElapsed
SafeHoldApplied
```

这些事件可以绕过普通测量分发，直接改变运行状态。

正常路径中的 `MeasurementSucceeded` 不匹配这些分支，因此继续向下。

##### 4.2 检查 `command_id`

随后比较：

```python
event.command_id == self._pending_command_id
```

假设当前 pending ID 是 `"abc-123"`，却收到：

```python
MeasurementSucceeded(command_id="old-456")
```

状态机直接忽略，不记录、不转移。

这是为了防止旧命令的迟到结果污染当前状态。

##### 4.3 检查 `state_version`

接着检查：

```python
event.state_version in (0, self._state_version)
```

其中 `0` 表示兼容同步执行器，没有显式填写版本。否则事件必须属于当前状态版本。

所以事件要同时回答：

```text
它属于哪条命令？
它属于哪个状态时期？
```

##### 4.4 记录测量成本

如果事件是 `MeasurementSucceeded`，状态机记录：

```python
event.shots
event.diagnostics["solver_calls"]
```

这发生在具体状态判断之前，因此即使随后结果被判定为无效，实际发生的测量成本仍会被计入。

##### 4.5 按当前状态分派

最后才执行：

```python
match self._state:
    case ACQUIRE:
        self._handle_acquire(event)
    case TRACK:
        self._handle_track(event)
    ...
```

也就是说，事件的解释依赖当前状态。

同样一条 `MeasurementSucceeded`：

- 在 Acquire 中被解释为频率种子；
- 在 Track 中被解释为局部控制测量；
- 在 Verify 中被解释为独立确认；
- 在 Lock 中被解释为漂移监测。

---

### 2.2 Acquire：先建立可靠的绝对种子

#### 5. Acquire 如何处理成功事件

Acquire 的处理器是 _handle_acquire()。

它首先匹配可靠测量：

```python
MeasurementSucceeded(
    valid=True,
    ambiguous=False,
    frequency=f,
    uncertainty=u,
) if u <= sigma_acquire_max
```

因此一个可靠 Acquire 结果必须同时满足：

```text
测量技术上成功
valid=True
ambiguous=False
不确定度不超过 sigma_acquire_max
```

满足后，状态机保存：

```python
self._f_hat = f
self._uncertainty = u
self._candidate_bias = event.applied_bias
self._candidate_source_event = "acquire"
self._tracker_snapshot = None
```

这里建立了 Track 或 Verify 后续需要的初始知识：

```text
当前频率估计      f_hat
当前不确定度      uncertainty
测量发生的偏置    candidate_bias
候选来源          acquire
```

然后计算目标残差：

$$
r_{\mathrm{acq}}
=
f_{\mathrm{acq}}-f_{\mathrm{target}}.
$$

候选条件由 `_check_candidate()` 判断：

$$
U_{\mathrm{acq}}
=
|r_{\mathrm{acq}}|
+
z\sigma_{\mathrm{acq}}
\leq
\epsilon_{\mathrm{enter}}.
$$

出现两个正常分支。

##### 5.1 已经接近目标：Acquire → Verify

如果：

$$
U_{\mathrm{acq}}\leq\epsilon_{\mathrm{enter}},
$$

则：

```python
_transition_to(
    FrequencyState.VERIFY,
    ReasonCode.CANDIDATE_REACHED,
)
```

因为 Acquire 已经提供了可靠且接近目标的偏置，没有必要先做一次 Track。

路径是：

```text
Acquire
  └─ 可靠且候选成立 → Verify
```

##### 5.2 可靠但尚未接近目标：Acquire → Track

否则：

```python
_transition_to(
    FrequencyState.TRACK,
    ReasonCode.SEED_RELIABLE,
)
```

含义是：

> 绝对频率参考已经可靠，但还需要局部控制把工作点移动到目标附近。

路径是：

```text
Acquire
  └─ 可靠但未达候选 → Track
```

最后 `_handle_acquire()` 调用：

```python
self._clear_pending()
```

清除刚才的 Acquire 命令。下一次 `next_command()` 才会根据新状态创建一条新命令。

---

#### 6. 一个具体的 Acquire 例子

假设：

$$
f_{\mathrm{target}}/(2\pi)=5.000\ \mathrm{GHz},
$$

Acquire 测得：

$$
\hat f_{\mathrm{acq}}/(2\pi)=5.030\ \mathrm{GHz},
$$

即目标残差为：

$$
r_{\mathrm{acq}}/(2\pi)=30\ \mathrm{MHz}.
$$

假设：

$$
\sigma_{\mathrm{acq}}/(2\pi)=0.5\ \mathrm{MHz},
\qquad
z=1,
\qquad
\epsilon_{\mathrm{enter}}/(2\pi)=5\ \mathrm{MHz}.
$$

则：

$$
U_{\mathrm{acq}}/(2\pi)
=
30+0.5
=
30.5\ \mathrm{MHz}.
$$

由于：

$$
30.5>5,
$$

因此：

```text
Acquire → Track
reason = SEED_RELIABLE
```

状态机现在知道频率约为 5.030 GHz，但还没有达到 5.000 GHz 的目标。

---

### 2.3 Track：局部证据与目标候选分开判断

#### 7. Track：发出局部测量命令

进入 Track 后，再次调用：

```python
cmd = machine.next_command()
```

`_dispatch_command()` 会先构造一个基础 `TrackFrequency`：

```python
TrackFrequency(
    candidate_bias=self._candidate_bias,
    predicted_drive=self._f_hat,
    predicted_detuning=0.0,
    predicted_detuning_uncertainty=self._uncertainty,
    prediction_guard_complete=True,
    tracker_spec=...,
)
```

这里的初始逻辑是：

```text
在 Acquire 测量得到的偏置开始
驱动频率设为刚得到的频率估计
```

所以即使目标还差 30 MHz，驱动可以跟随到当前频率附近：

$$
f_{\mathrm{drive}}\approx\hat f_{\mathrm{acq}}.
$$

这会使探针失谐：

$$
\Delta=f-f_{\mathrm{drive}}
$$

保持较小。

需要注意：状态机这里产生的只是基础 Track 命令。实际 Runtime 会结合 `DampedSecantTracker` 的快照，把它补全为真正的下一步偏置和预测驱动。这个过程属于第三、四遍内容。

---

#### 8. Track 测量返回什么

Backend 返回的正常 Track 事件大致是：

```python
MeasurementSucceeded(
    command_id=cmd.command_id,
    frequency=f_measured,
    uncertainty=sigma,
    valid=True,
    ambiguous=False,
    applied_bias=actual_bias,
    applied_drive=actual_drive,
    probe_detuning=f_measured - actual_drive,
    probe_detuning_uncertainty=sigma_delta,
    diagnostics={
        "s_hat": ...,
        "confidence": ...,
        "out_of_range": False,
    },
)
```

这里有两条并行信息：

```text
目标残差：
r = frequency - f_target

探针失谐：
Δ = frequency - applied_drive
```

例如：

```text
目标还差 30 MHz
驱动跟随当前频率
实际探针失谐只有 0.1 MHz
```

这在 Track 中完全合法：

$$
|r|=30\ \mathrm{MHz},
\qquad
|\Delta|=0.1\ \mathrm{MHz}.
$$

目标残差大，说明还需继续调节；探针失谐小，说明短脉冲测量仍然有效。

---

#### 9. Track 处理器的判断顺序

Track 事件进入 _handle_track()。

首先保存测量结果：

```python
error = f - self.f_target
self._f_hat = f
self._uncertainty = u
self._candidate_bias = event.applied_bias
```

这里的 `error` 是目标残差 \(r\)，不是探针失谐。

随后取得本轮 pending command：

```python
pending = self._pending_command
```

这样状态机可以比较：

```text
命令要求的偏置/驱动
和
事件报告的实际偏置/驱动
```

然后进入 `_track_guard_failure()`。

##### 9.1 Track guard 先判断测量能不能信

Track guard 的思想是：

> 在判断是否接近目标之前，先确认这次局部测量本身是否合法。

它依次检查：

1. `frequency` 和 `uncertainty` 是否有限；
2. 实际偏置是否等于命令偏置；
3. 实际驱动是否等于命令驱动；
4. Backend 是否报告 `out_of_range`；
5. 是否具有明确的 `probe_detuning`；
6. 探针失谐是否满足局部范围；
7. confidence 是否足够；
8. 灵敏度 \(s_{\mathrm{hat}}\) 是否正常。

核心范围条件是：

$$
|\widehat\Delta|
+
z\sigma_\Delta
\leq
\Delta_{\mathrm{val}}.
$$

只有 guard 全部通过，才有资格继续判断目标候选条件。

这种顺序很重要：

```text
先问：这个局部测量有效吗？
再问：测量结果接近目标了吗？
```

不能反过来。

---

#### 10. Track 的三个正常出口

##### 10.1 执行不一致：Track → SafeStop

如果实际偏置或驱动与命令不一致：

```python
guard_failure == ReasonCode.INTERLOCK
```

则：

```text
Track → SafeStop
run_status = Failed
```

这是安全失败，不是普通 Reacquire。因为此时无法确定测量对应哪个实际控制条件。

##### 10.2 局部证据失效：Track → Reacquire

如果出现：

```text
LOCAL_RANGE_LOST
CONFIDENCE_INSUFFICIENT
SENSITIVITY_INVALID
```

则：

```text
Track → Reacquire
```

状态机不允许继续外推局部反演。

##### 10.3 guard 有效：判断目标候选

如果 guard 全部通过，再计算：

$$
U_{\mathrm{tra}}
=
|r_{\mathrm{tra}}|
+
z\sigma_{\mathrm{tra}}.
$$

如果：

$$
U_{\mathrm{tra}}\leq\epsilon_{\mathrm{enter}},
$$

则：

```text
Track → Verify
reason = CANDIDATE_REACHED
```

同时记录：

```python
_candidate_source_event = "track"
_candidate_bias_version += 1
```

表示当前候选偏置来自 Track。

如果尚未满足候选条件：

```python
self._stay(
    FrequencyState.TRACK,
    ReasonCode.CORRECTION_REQUIRED,
)
```

即：

```text
Track → Track
```

下一轮 Runtime 会基于控制器快照提出新的偏置。

---

#### 11. 自循环 `_stay()` 并非“什么都没发生”

Track 未达到目标时调用：

```python
_stay(TRACK, CORRECTION_REQUIRED)
```

虽然状态名称没有变化，但 `_stay()` 仍然：

```python
self._state_version += 1
self._log_transition(TRACK, TRACK, reason)
```

因此：

```text
Track version 3 → Track version 4
```

这表示：

> 完成了一轮有效 Track 测量，虽然仍处于 Track，但协议已经进入下一轮。

随后 `_clear_pending()` 清除旧命令。下一次 `next_command()` 会创建具有新 `command_id` 的新 Track 命令。

完整的 Track 自循环是：

```text
Track version 3
  ↓ next_command
TrackFrequency command A
  ↓ measurement
MeasurementSucceeded A
  ↓ handle
Track version 4，CORRECTION_REQUIRED
  ↓ clear pending
下一条 TrackFrequency command B
```

---

#### 12. 一个 Track 数值例子

延续前面的 30 MHz 初始残差。

第一轮 Track 得到：

$$
r_1/(2\pi)=12\ \mathrm{MHz},
$$

$$
\Delta_1/(2\pi)=0.2\ \mathrm{MHz}.
$$

假设：

$$
\Delta_{\mathrm{val}}/(2\pi)=20\ \mathrm{MHz},
\qquad
\epsilon_{\mathrm{enter}}/(2\pi)=5\ \mathrm{MHz}.
$$

则：

```text
探针失谐 0.2 MHz < 20 MHz：局部测量有效
目标残差 12 MHz > 5 MHz：尚未成为候选
```

所以：

```text
Track → Track
reason = CORRECTION_REQUIRED
```

第二轮得到：

$$
r_2/(2\pi)=3\ \mathrm{MHz},
\qquad
\Delta_2/(2\pi)=0.1\ \mathrm{MHz}.
$$

此时：

```text
局部测量有效
目标残差进入 5 MHz 候选区
```

所以：

```text
Track → Verify
reason = CANDIDATE_REACHED
```

---

### 2.4 Verify：冻结控制量，取得独立证据

#### 13. Verify：命令为什么使用 frozen bias

进入 Verify 后：

```python
cmd = machine.next_command()
```

生成：

```python
VerifyFrequency(
    frozen_bias=self._candidate_bias,
    drive=self.f_target,
)
```

这里有两个关键语义。

##### 13.1 偏置被冻结

`frozen_bias` 是刚才 Track 达到候选条件时的实际偏置。

在 Verify 阶段，不再调用 tracker，也不提出新偏置。连续多次 Verify 都应使用同一个冻结偏置。

##### 13.2 驱动设为目标频率

```python
drive = self.f_target
```

Verify 直接测量候选频率相对于最终目标的偏差，而不是继续沿用 Track 的局部 drive tracking 参考。

Backend 会在该冻结偏置上执行独立双扫 Ramsey。

---

#### 14. Verify 首先检查“真的冻结了吗”

Verify 的处理器是 _handle_verify()。

正常事件到来后，它首先检查：

```python
abs(event.applied_bias - self._candidate_bias)
<= bias_freeze_tolerance
```

这一步发生在频率容差判断之前。

原因是：如果命令要求在偏置 \(\Phi_c\) 验证，但硬件实际在另一个偏置测量，那么无论频率结果多接近目标，都不能证明候选点有效。

若冻结偏置不一致：

```text
Verify → SafeStop
reason = INTERLOCK
run_status = Failed
```

---

#### 15. Verify 如何判断通过

偏置检查通过后，计算：

$$
r_{\mathrm{ver}}
=
f_{\mathrm{ver}}-f_{\mathrm{target}}.
$$

然后判断：

$$
U_{\mathrm{ver}}
=
|r_{\mathrm{ver}}|
+
z\sigma_{\mathrm{ver}}
\leq
\epsilon_{\mathrm{final}}.
$$

如果通过：

```python
self._verify_streak += 1
```

但不一定立即进入 Lock。

假设：

```python
N_verify = 2
```

第一次通过时：

```text
Verify → Verify
reason = VERIFY_PASSED
verify_streak = 1
```

状态仍是 Verify，但状态版本增加，旧 pending command 被清除。下一轮会生成一条新的 `VerifyFrequency`，再次在同一冻结偏置上测量。

第二次连续通过时：

```text
verify_streak = 2
Verify → Lock
reason = VERIFY_PASSED
```

同时：

```python
run_status = CALIBRATED
```

---

#### 16. 为什么第一次通过仍然 `_stay(VERIFY)`

这与 Track 自循环类似：

```text
Verify version 5
  ↓ 第一次独立 Ramsey 通过
Verify version 6
  ↓ 第二次独立 Ramsey 通过
Lock version 7
```

第一次 `_stay()` 不是重复使用同一份数据，而是要求下一轮生成新的 Verify 命令、重新测量。

所以连续命中的代码含义是：

```text
不同 command_id
不同测量事件
相同 frozen_bias
连续满足 epsilon_final
```

而不是把一次双扫 Ramsey 的正负分支分别算作两次命中。

---

#### 17. Verify 未通过时为什么可能返回 Track

若：

$$
U_{\mathrm{ver}}>\epsilon_{\mathrm{final}},
$$

Verify 不一定直接 Reacquire。

它还检查：

$$
U_{\mathrm{ver}}
\leq
\texttt{verify\_track\_max\_residual}.
$$

如果结果可靠、偏差不大，说明：

> 候选没有达到最终容差，但仍处于局部控制可以修正的范围。

于是：

```text
Verify → Track
reason = CORRECTION_REQUIRED
```

此时：

```python
_verify_streak = 0
```

并更新 tracker 快照中的：

```python
f
e
uncertainty
streak
```

这里有一处很重要的设计：Verify 返回 Track 时没有清空整个割线历史，而是用独立 Ramsey 结果更新当前点。

否则如果完全清空 tracker，下一次 Track 可能只在同一偏置重新测量，再次满足 `epsilon_enter`，立即返回 Verify，形成：

```text
Track ↔ Verify
```

的无效振荡。

若 Verify 偏差已经超过局部恢复范围，则：

```text
Verify → Reacquire
reason = CONFIDENCE_INSUFFICIENT
```

---

#### 18. 一个 Verify 数值例子

候选进入条件是：

$$
\epsilon_{\mathrm{enter}}/(2\pi)=5\ \mathrm{MHz}.
$$

Track 测得残差 3 MHz，因此进入 Verify。

但 Verify 使用更严格的：

$$
\epsilon_{\mathrm{final}}/(2\pi)=0.1\ \mathrm{MHz}.
$$

第一次独立 Ramsey：

$$
r_{\mathrm{ver},1}/(2\pi)=0.04\ \mathrm{MHz},
$$

$$
\sigma_{\mathrm{ver},1}/(2\pi)=0.02\ \mathrm{MHz}.
$$

所以：

$$
U_{\mathrm{ver},1}/(2\pi)
=
0.04+0.02
=
0.06\ \mathrm{MHz}.
$$

满足 0.1 MHz：

```text
Verify → Verify
verify_streak = 1
```

第二次：

$$
r_{\mathrm{ver},2}/(2\pi)=0.05\ \mathrm{MHz},
\qquad
\sigma_{\mathrm{ver},2}/(2\pi)=0.02\ \mathrm{MHz},
$$

得到：

$$
U_{\mathrm{ver},2}/(2\pi)=0.07\ \mathrm{MHz}.
$$

再次满足条件：

```text
Verify → Lock
verify_streak 达到 2
run_status = CALIBRATED
```

---

### 2.5 Lock：通过验证后继续长期监测

#### 19. 进入 Lock 时初始化了什么

Verify 连续通过后，代码不只是改变状态，还初始化长期运行信息：

```python
self._lock_entry_time = time.time()
self._monitor_streak = 0
self._last_audit_time = time.time()
self._hold_target_met = None
self._hold_samples = 0
self._hold_passes = 0
```

含义分别是：

```text
从什么时候开始 Lock
当前连续可疑次数归零
从现在开始计算周期性独立复核间隔
尚未获得 hold 监测样本
hold 样本总数归零
满足 epsilon_hold 的样本数归零
```

此时：

```text
state      = Lock
run_status = Calibrated
```

任务并没有终止。

---

#### 20. Lock 的下一条命令是什么

再次调用 `next_command()` 时，Lock 分支先检查周期性独立复核是否到期：

```text
if now - last_audit_time >= audit_interval:
    Lock → Verify
    return VerifyFrequency(...)
```

如果未到期，则返回：

```python
MonitorFrequency(
    locked_bias=self._candidate_bias,
)
```

因此 Lock 有两种输出：

```text
普通周期       → 低成本 MonitorFrequency
独立复核到期   → VerifyFrequency
```

无论哪一种，都不会生成新的磁通偏置。

---

#### 21. Lock monitor 如何决定下一步

Lock 的事件处理器是 _handle_lock()。

它不会一收到频率值就解释目标残差，而是先验证这次局部测量本身。`MonitorFrequency`
同时指定 `locked_bias` 和 `drive`，事件必须回报一致的实际偏置、实际驱动，以及：

$$
\widehat\Delta_{\mathrm{mon}}
=
\widehat\omega_{01}-\omega_d,
$$

$$
|\widehat\Delta_{\mathrm{mon}}|
+z\sigma_{\Delta,\mathrm{mon}}
\leq
\Delta_{\mathrm{val}}.
$$

执行值不一致是联锁错误，进入 SafeStop；探针失谐越界、诊断缺失或置信度不足说明局部
参考不再可信，进入 Reacquire。只有该守卫通过，下面的目标残差才有资格参与 Lock 路由。
内置确定性和有限采样后端都会提供这些字段；旧自定义后端可暂时设置
`require_monitor_local_validity=False` 迁移，但这会恢复较弱的旧契约。

成功测量后计算：

$$
r_{\mathrm{mon}}
=
f_{\mathrm{mon}}-f_{\mathrm{target}},
$$

$$
U_{\mathrm{mon}}
=
|r_{\mathrm{mon}}|
+
z\sigma_{\mathrm{mon}}.
$$

同时单独记录是否满足长期目标：

$$
U_{\mathrm{mon}}
\leq
\epsilon_{\mathrm{hold}}.
$$

然后根据 monitor 阈值路由。

##### 21.1 清晰区：Lock → Lock

如果：

$$
U_{\mathrm{mon}}
\leq
\epsilon_{\mathrm{mon,clear}},
$$

则：

```text
Lock → Lock
reason = MONITOR_CLEAR
monitor_streak = 0
```

这是正常长期自循环。

##### 21.2 可疑区：先累计次数

如果：

$$
\epsilon_{\mathrm{mon,clear}}
<
U_{\mathrm{mon}}
\leq
\epsilon_{\mathrm{mon,suspect}},
$$

则：

```python
monitor_streak += 1
```

在尚未达到 `N_mon_suspect` 时：

```text
Lock → Lock
reason = MONITOR_SUSPECT
```

达到连续次数时：

```text
Lock → Verify
reason = DRIFT_SUSPECTED
```

注意，不会直接进入 Track。

##### 21.3 大跳变：Lock → Reacquire

如果达到大跳变阈值：

```text
Lock → Reacquire
reason = LARGE_FREQUENCY_JUMP
```

此时局部参考可能已经失效，不值得先尝试局部修正。

---

### 2.6 把四个状态重新连成一条时间线

#### 22. 第二遍的完整正常时间线

把前面的代码连起来，就是：

```text
machine.start()
state = Acquire
run_status = Running

next_command()
→ AcquireFrequency A

handle(Acquire result)
→ 保存频率种子
→ 目标尚远
→ Acquire → Track
→ 清除 command A

next_command()
→ TrackFrequency B

handle(Track result B)
→ 局部 guard 有效
→ 尚未达到候选条件
→ Track → Track
→ 清除 command B

next_command()
→ TrackFrequency C

handle(Track result C)
→ 局部 guard 有效
→ 候选条件成立
→ Track → Verify
→ 冻结 candidate bias
→ 清除 command C

next_command()
→ VerifyFrequency D

handle(Verify result D)
→ 偏置冻结正确
→ 最终容差通过
→ Verify → Verify
→ streak = 1
→ 清除 command D

next_command()
→ VerifyFrequency E

handle(Verify result E)
→ 再次通过
→ Verify → Lock
→ run_status = Calibrated
→ 初始化 Lock 记录
→ 清除 command E

next_command()
→ MonitorFrequency F

handle(Monitor result F)
→ monitor clear
→ Lock → Lock
```

对应的状态序列是：

```text
Acquire(v0)
  ↓ SEED_RELIABLE
Track(v1)
  ↓ CORRECTION_REQUIRED
Track(v2)
  ↓ CANDIDATE_REACHED
Verify(v3)
  ↓ VERIFY_PASSED
Verify(v4)
  ↓ VERIFY_PASSED
Lock(v5)
  ↓ MONITOR_CLEAR
Lock(v6)
```

---

### 2.7 第二遍收束：从记录反推一轮决策

#### 23. 第二遍最应该掌握的三个原则

第一，**发出命令不等于状态转换**：

```text
next_command() 只创建并冻结 pending command
handle(event) 才根据证据改变状态
```

第二，**每轮结果都要清除 pending command**：

```text
命令 A
→ 事件 A
→ 状态转移
→ clear pending
→ 才能生成命令 B
```

第三，**正常路径上存在三种不同自循环**：

```text
Track  → Track   仍需局部修正
Verify → Verify  还差连续验证次数
Lock   → Lock    长期监测正常
```

它们虽然都是状态名称不变，但物理含义完全不同，而且每次都会增加 `state_version`、记录转换原因并生成新的测量命令。

第二遍读懂的标志，是你能拿任意一条 journal 记录，从：

```text
state_before
command
event
transition_reason
state_after
```

反推出这一轮为什么发生。下一遍再进入 Backend，追踪这些命令究竟怎样变成 Ramsey、短脉冲测量和 `MeasurementSucceeded`。


## 第三遍：沿命令进入测量后端

```{figure} frequency_state_machine_pass3.svg
:alt: 第三遍高亮测量后端数据通路
:width: 100%
:align: center

第三遍的激活路径：Command 从执行边界进入 Issue 与 Measurement Unit，测量结果被封装成 Event；其他模块仅说明输入来自哪里、结果将去哪里。
```

图中的 **Measurement Unit** 对应 Executor 与 Backend 的职责集合，而不是单个类。`CMD REG` 强调 Backend 接收的是已经确定 role、bias、drive 和 identity 的命令；`EVT REG` 强调 Backend 返回证据，不直接写状态寄存器，也不自行决定 Acquire、Track、Verify 或 Lock 的转移。

第三遍的目标是打开第二遍里暂时视为黑箱的部分：

```python
event = executor.execute(command)
```

我们要回答：

1. `Command` 如何被翻译成具体测量？
2. Acquire、Track、Verify、Lock 分别使用什么 probe？
3. 原始布居怎样变成频率？
4. 确定性仿真和有限-shot 仿真有什么区别？
5. Backend 应负责什么，不应负责什么？

主要阅读两个文件：

- frequency_backends.py（`sqc/workflows/frequency_backends.py` 源码约第 139 行）：状态机命令与测量实现之间的适配层；
- frequency.py（`sqc/calibration/frequency.py` 源码约第 777 行）：真正的单点频率测量算法。

整体调用链是：

```text
Command
  ↓
SQCExecutor.execute()
  ↓
_execute_acquire / track / verify / monitor
  ↓
FrequencyMeasurement.measure()
  ↓
Ramsey 或 transient 底层算法
  ↓
frequency / populations / estimator diagnostics
  ↓
MeasurementSucceeded
  ↓
FrequencyStateMachine.handle()
```

### 3.1 先确定 Backend 的责任边界

#### 3.1 为什么需要 Backend 层

状态机只认识抽象协议：

```text
AcquireFrequency
TrackFrequency
VerifyFrequency
MonitorFrequency
```

底层测量只认识物理参数：

```text
flux
omega_d
tau_list
pulse shape
measurement method
```

Backend 的作用就是把两种语言互相翻译。

例如状态机发出：

```python
TrackFrequency(
    candidate_bias=0.94,
    predicted_drive=2 * np.pi * 5.03,
)
```

Backend 将其翻译为：

```python
measurement.measure(
    flux=0.94,
    omega_d=2 * np.pi * 5.03,
)
```

测量得到频率后，再包装成：

```python
MeasurementSucceeded(
    frequency=...,
    applied_bias=0.94,
    applied_drive=2 * np.pi * 5.03,
    probe_detuning=frequency - applied_drive,
)
```

所以：

```text
状态机决定“为什么测、测完以后去哪”
Backend 决定“如何执行这次测量”
FrequencyMeasurement 决定“如何从物理响应估计频率”
```

这种分层保证状态转移不会散落进测量代码。

---

#### 3.2 从 `SQCExecutor.execute()` 开始

入口是 SQCExecutor.execute()（`sqc/workflows/frequency_backends.py` 源码约第 168 行）。

它本质上是命令分派器：

```python
match cmd:
    case AcquireFrequency():
        event = self._execute_acquire(cmd)
    case TrackFrequency():
        event = self._execute_track(cmd)
    case VerifyFrequency():
        event = self._execute_verify(cmd)
    case MonitorFrequency():
        event = self._execute_monitor(cmd)
    case SafeHold():
        event = self._execute_safe_hold(cmd)
```

映射关系为：

| Command | Backend 方法 | 默认 probe |
|---|---|---|
| `AcquireFrequency` | `_execute_acquire()` | 扫 \(\tau\) Ramsey，单扫 |
| `TrackFrequency` | `_execute_track()` | \(\tau=0\) 正交短脉冲 |
| `VerifyFrequency` | `_execute_verify()` | 扫 \(\tau\) Ramsey，双扫 |
| `MonitorFrequency` | `_execute_monitor()` | \(\tau=0\) 正交短脉冲 |
| `SafeHold` | `_execute_safe_hold()` | 不做科学测量 |

`execute()` 还有一个重要职责：统计实际 solver calls。

---

#### 3.3 `_SolverCallCounter`：测量成本怎样被记录

每条命令都在：

```python
with _SolverCallCounter() as counter:
    ...
```

内部执行。

它临时包装：

```text
qutip.mesolve
qutip.sesolve
frequency 模块中的 mesolve
kernel 模块中的 mesolve
```

每次调用都累计一次。执行结束后，将实际计数写入事件：

```python
event.diagnostics["solver_calls"] = counter.total
event.diagnostics["solver_call_breakdown"] = {
    "mesolve": ...,
    "sesolve": ...,
}
```

需要区分两个成本概念。

##### 测前成本估计

estimate_cost()（`sqc/workflows/frequency_backends.py` 源码约第 207 行） 在执行前给 Runtime 一个保守估计：

```text
“这条命令大约需要多少 shots 和 solver calls？”
```

Runtime 用它判断是否还有足够预算。

##### 测后实际成本

`_SolverCallCounter` 记录真正发生的 solver 调用。

因此：

```text
estimate_cost()  → 测前预算预留
counter.total    → 测后实际记账
```

二者不必完全相同。测前估计应该保守，测后日志用于真实成本分析。

---

#### 3.4 测量对象为什么采用 lazy creation

`SQCExecutor` 内部保存三个测量对象：

```python
_ramsey_meas
_transient_meas
_verify_meas
```

它们初始都是 `None`，只有第一次需要时才创建：

```python
def _get_ramsey():
    if self._ramsey_meas is None:
        self._ramsey_meas = FrequencyMeasurement(...)
    return self._ramsey_meas
```

这样做有两个作用：

1. 不使用某种测量时，不创建对应对象；
2. 同一种测量在多轮状态机中复用相同脉冲配置和内部缓存。

三个对象的默认配置不同。

##### Acquire Ramsey

```python
FrequencyMeasurement(
    method="ramsey",
    f_artificial=0.1,
)
```

即单扫 Ramsey，人工失谐为 \(0.1\ \mathrm{GHz}=100\ \mathrm{MHz}\)。

##### Track/Monitor transient

```python
FrequencyMeasurement(
    method="transient",
    order=3,
    g3_source="kernel_full",
)
```

默认使用三阶修正和完整三阶响应核。

##### Verify Ramsey

```python
FrequencyMeasurement(
    method="ramsey",
    f_artificial=None,
)
```

`None` 表示正负人工失谐双扫 Ramsey。

---

#### 3.5 `FrequencyMeasurement`：统一的单点频率计

FrequencyMeasurement（`sqc/calibration/frequency.py` 源码约第 777 行） 是 Backend 下面的统一测量接口。

最重要的方法是：

```python
measure(
    flux=...,
    omega_d=...,
) -> float
```

它解决的是：

> 在一个指定磁通偏置和驱动参考下，估计量子比特绝对频率。

内部根据 `method` 分派：

```python
match self.method:
    case "ramsey":
        return _fit_ramsey_frequency(...)
    case "transient":
        return _measure_frequency_transient(...)
```

它是只读测量对象，不负责：

- 修改状态机状态；
- 决定是否进入 Verify；
- 计算下一磁通偏置；
- 触发 Reacquire。

这些决策都在更上层。

---

#### 3.6 `omega_d` 是整个测量链的关键参考

`FrequencyMeasurement.measure()` 中，驱动频率按以下优先级确定：

```text
调用参数 omega_d
    ↓ 若没有
对象配置 self.omega_d
    ↓ 若没有
qubit.frequency
```

它返回的是：

$$
\hat\omega_{01}
=
\omega_d+\widehat\Delta.
$$

所以测量首先估计量子比特相对当前驱动的失谐，再把驱动频率加回来形成绝对频率。

在真实实验中，真正的 \(\omega_{01}\) 本来未知，因此只能设置一个参考驱动 \(\omega_d\)，然后测量：

$$
\Delta=\omega_{01}-\omega_d.
$$

在 Track 中，把前一轮预测作为下一轮 `omega_d`，就是 drive tracking：

```text
前一轮估计频率
       ↓
下一轮驱动参考
       ↓
使实际 probe detuning 保持较小
```

`qubit.frequency` 作为默认参考在仿真中方便，但它可能包含仿真对象已知的构造频率。因此闭环关键路径都应显式传入 `omega_d`，避免无意使用仿真 oracle。

---

### 3.2 Acquire：宽范围种子怎样产生

#### 3.7 Acquire 后端：建立宽范围频率种子

入口是 _execute_acquire()（`sqc/workflows/frequency_backends.py` 源码约第 234 行）。

核心调用：

```python
meas = self._get_ramsey()
frequency = meas.measure(flux=cmd.bias)
```

默认 Acquire 测量器为：

```python
method="ramsey"
f_artificial=0.1
```

即对一组 `tau_list` 执行一次 Ramsey 扫描，并加入已知人工失谐：

$$
f_a=0.1\ \mathrm{GHz}.
$$

测得 Ramsey 条纹频率 \(f_{\mathrm{meas}}\) 后，按当前符号约定恢复：

$$
\Delta=f_a-f_{\mathrm{meas}},
$$

$$
\hat\omega_{01}
=
\omega_d+2\pi\Delta.
$$

##### 为什么加入人工失谐

普通 Ramsey FFT 主要得到振荡频率绝对值，符号可能不明确。若人工失谐满足：

$$
f_a>|\Delta|_{\max},
$$

则可以从：

$$
f_{\mathrm{meas}}=|f_a-\Delta|=f_a-\Delta
$$

恢复带符号失谐。

##### 当前 Acquire 的边界

虽然注释称其为 wide-range acquisition，但当前默认单扫仍要求：

$$
|\Delta|<100\ \mathrm{MHz}.
$$

它比 Track 的局部窗口宽得多，但不是数学意义上的无限捕获范围。若硬件启动误差可能超过该范围，Acquire 需要采用更宽范围 spectroscopy、多尺度 Ramsey，或双扫方案。

这是读代码时必须看到的实现边界。

##### Acquire 返回的事件

确定性后端返回：

```python
MeasurementSucceeded(
    frequency=frequency,
    uncertainty=0.0,
    valid=True,
    ambiguous=False,
    method="ramsey",
    shots=0,
    applied_bias=cmd.bias,
    diagnostics={
        "role": cmd.role,
        "circuits": len(tau_list),
        "uncertainty_source": "deterministic_zero",
    },
)
```

这里：

- `shots=0` 表示这是确定性期望值仿真，不是没有测量成本；
- `circuits` 和实际 solver calls 仍被记录；
- `uncertainty=0` 明确来自 `deterministic_zero`，不能解释成实验上具有无限精度。

---

#### 3.8 Ramsey 底层怎样得到频率

Ramsey 测量位于 frequency.py（`sqc/calibration/frequency.py` 源码约第 110 行）。

对每个自由演化时间 \(\tau\)，执行：

$$
\frac{\pi}{2}
\rightarrow
\tau
\rightarrow
\frac{\pi}{2}
\rightarrow
\text{readout}.
$$

得到：

$$
p_e(\tau).
$$

然后 `_fft_peak()`：

1. 检查条纹对比度；
2. 去除直流分量；
3. 计算 FFT；
4. 找主峰；
5. 使用二次插值提高频率分辨率。

##### 单扫模式

配置：

```python
f_artificial = 0.1
```

只生成一条 Ramsey 扫描：

$$
p_+(\tau).
$$

成本约为：

$$
N_\tau
$$

条电路。

优点是成本较低，缺点是需要已知捕获范围。

##### 双扫模式

配置：

```python
f_artificial = None
```

内部使用：

$$
+f_a,\qquad -f_a,
$$

当前：

$$
f_a=0.05\ \mathrm{GHz}=50\ \mathrm{MHz}.
$$

得到两条峰频：

$$
f_p=|f_a-\Delta|,
\qquad
f_n=|-f_a-\Delta|.
$$

利用：

$$
\Delta
=
\frac{f_n^2-f_p^2}{4f_a}
$$

恢复带符号失谐。

成本约为：

$$
2N_\tau.
$$

这就是 Verify 使用双扫 Ramsey 的原因：它成本较高，但避免让单扫的符号假设成为最终确认依据。

---

### 3.3 Track：局部 probe 怎样变成频率

#### 3.9 Track 后端：低成本局部短脉冲测频

入口是 _execute_track()（`sqc/workflows/frequency_backends.py` 源码约第 267 行）。

核心调用：

```python
frequency = measurement.measure(
    flux=cmd.candidate_bias,
    omega_d=cmd.predicted_drive,
)
```

这里必须使用命令给出的两个参数：

```text
candidate_bias    本轮实际测量工作点
predicted_drive   drive tracking 给出的参考频率
```

测量完成后，Backend 显式计算：

$$
\widehat\Delta
=
\hat f-\omega_{d,\mathrm{applied}}.
$$

代码是：

```python
detuning = frequency - cmd.predicted_drive
```

并写入：

```python
probe_detuning=detuning
probe_detuning_uncertainty=0.0
```

这个字段是状态机 Track guard 的唯一正确物理输入。目标残差：

$$
r=\hat f-f_{\mathrm{target}}
$$

由状态机另行计算，Backend 不应把二者混在一起。

---

#### 3.10 transient probe 具体测什么

底层是 _measure_frequency_transient()（`sqc/calibration/frequency.py` 源码约第 502 行）。

它不是扫 \(\tau\)，而是执行两条 \(\tau=0\) 的正交序列：

```text
Ry → Rx
Ry → R−x
```

更准确地说，两条 Ramsey 型控制分支使用不同的第二脉冲相位。

分别得到末态激发概率：

$$
p_{+X},\qquad p_{-X}.
$$

构造差分：

$$
p_d
=
\frac{p_{+X}-p_{-X}}{2}.
$$

差分能够：

- 抑制共同偏置；
- 保留失谐的一阶奇响应；
- 在共振附近提供带符号局部测量。

因此一次 transient 测量只需要两条电路，而不是 \(N_\tau\) 或 \(2N_\tau\) 条 Ramsey 扫描。

---

#### 3.11 怎样从差分布居反演频率

响应模型为：

$$
p_d
=
G_1\delta\omega
+
\frac{G_3}{6}\delta\omega^3.
$$

代码内部的 `delta_omega` 符号约定与论文中的 \(\Delta\) 有一次转换，因此最终使用：

$$
\hat\omega_{01}
=
\omega_d-\widehat{\delta\omega}.
$$

从用户角度更直观地看，就是：

$$
\hat\omega_{01}
=
\omega_d+\widehat\Delta.
$$

##### 一阶估计

先通过一阶响应核积分获得：

$$
G_1=\int k_1(t)\,dt.
$$

然后：

$$
\widehat{\delta\omega}^{(1)}
=
\frac{p_d}{G_1}.
$$

##### 三阶修正

默认：

```python
order=3
g3_source="kernel_full"
```

求解：

$$
p_d
=
G_1\widehat{\delta\omega}
+
\frac{G_3}{6}\widehat{\delta\omega}^{\,3}.
$$

`_solve_cubic_detuning()` 使用线性结果作为初始点，并检查 cubic fold，避免跳到远处的伪根。

##### \(G_3\) 的两种来源

##### `g3_source="fit"`

人为生成若干已知失谐下的响应，拟合奇多项式得到 \(G_1,G_3\)。

##### `g3_source="kernel_full"`

计算完整非对角三阶响应核：

$$
G_3
=
\iiint
k_3(t_1,t_2,t_3)
\,dt_1dt_2dt_3.
$$

当前状态机默认使用这一条。

这两个来源用于响应系数构建或数值交叉检查，不代表 Track 每轮都重新采集完整 response-calibration data。测量对象和相关系数具有缓存行为。

---

#### 3.12 Track 怎样标记局部越界

得到频率后，Backend 计算：

```python
out_of_range = (
    linear_range is not None
    and abs(detuning) + guard_margin > linear_range
)
```

并写入：

```python
diagnostics["out_of_range"]
diagnostics["detuning"]
```

但 Backend 只报告测量事实，不执行状态转移。

真正决定：

```text
Track → Reacquire
```

的是状态机 `_track_guard_failure()`。

因此责任边界是：

```text
Backend：
    计算 detuning
    报告 out_of_range
    返回测量证据

StateMachine：
    根据 Config 判断证据是否可接受
    选择 Track / Verify / Reacquire / SafeStop
```

---

#### 3.13 Track 也允许 Ramsey 模式

`SQCExecutor.track_method` 默认：

```python
track_method = "transient"
```

但也允许：

```python
track_method = "ramsey"
```

此时 `_get_track_measurement()` 返回 Verify 使用的双扫 Ramsey 测量器。

这主要用于：

- 对照实验；
- transient probe 不可用的环境；
- 测试不同测量成本和收敛行为。

即使 Track 使用 Ramsey，Backend 仍会填充：

```python
probe_detuning
probe_detuning_uncertainty
```

状态机并不依赖具体测量方法，只依赖事件协议。

这体现了 Command/Event 边界的价值：

```text
状态机判断局部证据
但不必知道证据来自 transient 还是 Ramsey
```

---

### 3.4 Verify、Lock 与 SafeHold：同是测量，不同权限

#### 3.14 Verify 后端：冻结偏置下的独立 Ramsey

入口是 _execute_verify()（`sqc/workflows/frequency_backends.py` 源码约第 314 行）。

核心代码：

```python
meas = self._get_verify()

frequency = meas.measure(
    flux=cmd.frozen_bias,
    omega_d=cmd.drive,
)
```

Verify 测量器为：

```python
FrequencyMeasurement(
    method="ramsey",
    f_artificial=None,
)
```

所以它执行的是扫 \(\tau\) 双扫 Ramsey。

##### 为什么与 Track 独立

Track 默认使用：

```text
τ=0 短脉冲差分 probe
局部 G1/G3 反演
```

Verify 使用：

```text
完整 τ 扫描
FFT/条纹拟合
正负人工失谐双扫
```

两者的控制序列和估计器都不同，因此 Track 连续命中不能冒充 Verify。

##### Verify 返回什么

确定性后端返回：

```python
MeasurementSucceeded(
    frequency=frequency,
    uncertainty=0.0,
    method="ramsey",
    applied_bias=cmd.frozen_bias,
    applied_drive=cmd.drive,
    diagnostics={
        "circuits": 2 * len(tau_list),
        "verifier_spec": ...,
        "uncertainty_source": "deterministic_zero",
    },
)
```

状态机随后检查：

1. `applied_bias` 是否仍等于冻结候选；
2. 目标残差是否满足 `epsilon_final`；
3. 是否达到连续 `N_verify` 次。

Backend 不决定是否进入 Lock。

---

#### 3.15 Lock monitor 后端：复用短脉冲，但不复用控制权

入口是 _execute_monitor()（`sqc/workflows/frequency_backends.py` 源码约第 348 行）。

它调用：

```python
frequency = transient.measure(
    flux=cmd.locked_bias,
    omega_d=self.f_target,
)
```

关键区别不在测量器本身，而在控制语义：

```text
Track transient：
    使用 candidate_bias 和 predicted_drive
    测完后控制器可能提出下一偏置

Lock transient：
    使用 locked_bias 和目标驱动
    测完后只能保持 Lock、进入 Verify 或 Reacquire
```

因为 Lock 中：

$$
\omega_d=\omega_{\mathrm{target}},
$$

所以：

$$
\Delta_{\mathrm{mon}}
=
\omega_{01}-\omega_d
=
\omega_{01}-\omega_{\mathrm{target}}
=
r_{\mathrm{mon}}.
$$

在这个特定场景中，probe detuning 与目标残差数值相同；但这是因为驱动恰好设置为目标，不代表两个概念一般可以混用。

一次 monitor 默认只需两个 transient 电路，远少于 Verify 的 \(2N_\tau\)。

---

#### 3.16 SafeHold 后端

`SafeHold` 不执行频率测量：

```python
return SafeHoldApplied(
    command_id=cmd.command_id,
    applied_bias=cmd.bias,
)
```

确定性仿真中它只是记录事件；真实硬件适配器中应在这里：

1. 设置安全 DAC 偏置；
2. 停止或关闭科学测量序列；
3. 返回实际应用的安全偏置；
4. 保证相同 `command_id` 重放时具有幂等性。

当前 `SQCExecutor` 仍是 QuTiP/器件对象参考后端，不是完整硬件驱动实现。

---

### 3.5 证据模型、成本与异常边界

#### 3.17 确定性 Backend 和有限-shot Backend

这是第三遍最重要的边界之一。

##### `SQCExecutor`

默认后端直接使用 QuTiP 求出的期望值：

```python
uncertainty = 0.0
shots = 0
diagnostics["uncertainty_source"] = "deterministic_zero"
```

它适合：

- 算法逻辑测试；
- 无采样噪声的数值收敛；
- solver-call 成本比较；
- 状态转移的确定性验证。

但这里的零不确定度表示：

> 当前数值后端没有加入统计采样误差。

它不表示真实实验测量具有零误差。

##### `FiniteShotSQCExecutor`

FiniteShotSQCExecutor（`sqc/workflows/frequency_backends.py` 源码约第 531 行） 在确定性布居之上加入：

- 每条电路的有限 shots；
- 二项采样；
- 可配置 readout assignment error；
- 按 Acquire、Track、Verify、Monitor 分开的随机数流；
- 可 checkpoint 的 RNG state。

它才会产生非零统计不确定度。

---

#### 3.18 有限-shot transient 如何工作

首先取得确定性末态布居：

```python
details = meas.measure_details(...)
```

其中包含：

```text
plus_x population
minus_x population
G1
G3
```

然后对每个分支执行二项采样：

$$
\hat p_\pm
=
\frac{\mathrm{Binomial}(N,p_\pm)}{N}.
$$

若考虑 assignment error：

$$
p_{\mathrm{obs}}
=
p_{01}
+
(1-p_{01}-p_{10})p.
$$

构造：

$$
\hat p_d
=
\frac{\hat p_+-\hat p_-}{2}.
$$

再使用已经固定的 \(G_1,G_3\) 反演频率。

统计误差采用 delta method：

$$
\sigma_{p_d}
=
\sqrt{
\frac{
p_+(1-p_+)+p_-(1-p_-)
}{
4N
}
},
$$

$$
\sigma_\Delta
\approx
\left|
\frac{\sigma_{p_d}}{G_1}
\right|.
$$

事件中记录：

```python
uncertainty=uncertainty
probe_detuning_uncertainty=uncertainty
shots=2 * shots_per_circuit
uncertainty_source="binomial_delta_method"
```

需要注意：当前误差传播主要按局部一阶灵敏度 \(G_1\) 计算，即使中心估计使用三阶反演。这是清晰的近似，不是完整的非线性置信区间。

---

#### 3.19 有限-shot Ramsey 如何工作

有限-shot Ramsey 先获得每个 \(\tau\) 上的确定性布居，再分别做二项采样。

然后 _ramsey_peak_with_uncertainty()（`sqc/workflows/frequency_backends.py` 源码约第 45 行）：

1. 用 FFT 找初始峰；
2. 以正弦模型进行加权非线性最小二乘；
3. 从 Jacobian 构造近似信息矩阵；
4. 取拟合协方差中的频率方差；
5. 拟合失败或矩阵奇异时使用由 \(\tau\) 扫描跨度决定的 fallback 分辨率。

对于双扫 Ramsey，两条峰的不确定度传播到：

$$
\Delta
=
\frac{f_n^2-f_p^2}{4f_a},
$$

近似得到：

$$
\sigma_\Delta
=
\frac{
\sqrt{(f_n\sigma_n)^2+(f_p\sigma_p)^2}
}{
2|f_a|
}.
$$

最后转换到角频率：

$$
\sigma_\omega=2\pi\sigma_\Delta.
$$

事件标记：

```python
method="ramsey_finite_shot"
uncertainty_source="weighted_ramsey_fit_covariance"
```

这比确定性 `sigma=0` 更接近实验数据结构，但仍是仿真中的二项噪声模型，不等于已经获得真机置信度。

---

#### 3.20 `measure()` 与 `measure_details()` 的区别

`FrequencyMeasurement` 有两个入口。

##### `measure()`

只返回：

```python
frequency: float
```

适合确定性 Backend，因为它只需要最终频率，再自行包装事件。

##### `measure_details()`

返回：

```python
FrequencyMeasurementDetails(
    frequency=...,
    method=...,
    populations=...,
    estimator=...,
    circuits=...,
)
```

有限-shot 后端需要原始分支布居，因此必须使用 `measure_details()`：

```text
先得到理想分支布居
再进行有限-shot 采样
再重新估计频率与不确定度
```

这也说明有限-shot 后端不是简单地给确定性频率加高斯噪声，而是在更接近测量数据的位置注入二项采样。

---

#### 3.21 异常如何穿过 Backend

每个 `_execute_*()` 通常使用：

```python
try:
    ...
    return MeasurementSucceeded(...)
except Exception as exc:
    return MeasurementTechnicalFailure(
        command_id=cmd.command_id,
        reason=str(exc),
        elapsed_time=...,
    )
```

所以：

```text
底层 solver 异常
脉冲构造异常
拟合代码异常
```

不会直接让状态机主循环崩溃，而会被转换成结构化技术失败事件。

状态机再根据：

```text
当前状态
技术重试次数
最大重试预算
```

决定原状态重试还是 SafeStop。

Backend 应报告失败事实，不应自己决定状态转移。

---

### 3.6 第三遍收束：沿一条命令穿过测量栈

#### 3.22 一次 Track 命令的完整跨层路径

把第三遍全部串起来：

```text
1. StateMachine 当前处于 Track

2. next_command()
   → TrackFrequency(
       candidate_bias=Φn+1,
       predicted_drive=ωd,n+1
     )

3. Runtime 将命令交给 SQCExecutor

4. SQCExecutor.execute()
   → _execute_track()

5. _get_track_measurement()
   → FrequencyMeasurement(method="transient")

6. FrequencyMeasurement.measure()
   → _measure_frequency_transient()

7. 在 Φn+1、ωd,n+1 下执行两条 τ=0 正交序列
   → p+x, p−x

8. 形成差分
   → pd = (p+x − p−x)/2

9. 使用 G1/G3 反演
   → Δ̂
   → f̂ = ωd,n+1 + Δ̂

10. Backend 计算
    → probe_detuning = f̂ − ωd,n+1

11. 包装事件
    → MeasurementSucceeded(
        frequency=f̂,
        applied_bias=Φn+1,
        applied_drive=ωd,n+1,
        probe_detuning=Δ̂,
        ...
      )

12. StateMachine.handle(event)
    → 检查执行一致性和局部范围
    → 再判断目标残差
    → Track / Verify / Reacquire / SafeStop
```

---

#### 3.23 第三遍应形成的架构认识

读完后，应该能清楚区分四层：

```text
FrequencyStateMachine
    决定协议状态与转移

FrequencyCalibrationRuntime
    调度命令、预算、中断和恢复

SQCExecutor / FiniteShotSQCExecutor
    将命令适配为具体测量并包装事件

FrequencyMeasurement
    从 Ramsey 或 transient 响应估计频率
```

也应能准确说出四个状态的默认测量方法：

```text
Acquire：
    扫 τ 单扫 Ramsey
    约 Nτ 条电路
    用于较宽范围种子

Track：
    τ=0 正交短脉冲
    2 条电路
    用于局部快速控制

Verify：
    扫 τ 双扫 Ramsey
    约 2Nτ 条电路
    用于冻结候选的独立确认

Lock monitor：
    τ=0 正交短脉冲
    2 条电路
    只监测，不更新偏置
```

第三遍最重要的三个结论是：

1. Backend 产生测量证据，但不决定状态转移；
2. Track 的 `probe_detuning` 必须由实际频率与实际驱动形成，不能用目标残差代替；
3. 确定性 `sigma=0` 和有限-shot 统计不确定度必须明确区分，二者都不能冒充真机实验置信度。

## 第四遍：Track 控制器如何跨轮形成闭环

```{figure} frequency_state_machine_pass4.svg
:alt: 第四遍高亮 Track 割线反馈环
:width: 100%
:align: center

第四遍的激活路径：接受的 Track 证据写入 `TR/ER`，Secant ALU 读取跨轮记忆生成 proposal，再由 Operand MUX 组成下一条 Track 命令。
```

这里的 `TR · Track` 对应 `TrackerSnapshot` 的已接受历史，`ER · Estimate` 提供当前频率与控制点，**Secant ALU** 对应 `DampedSecantTracker.propose()` 的纯计算，Operand MUX 对应 Runtime 把 proposal、drive tracking 和命令字段组合起来的过程。ALU 输出只是 proposal，仍须经过测量、守卫与提交。

第四遍聚焦 Track 控制器。前三遍已经知道：

```text
状态机决定是否继续 Track
Backend 在指定偏置上测量频率
```

现在要打开中间的问题：

> 如果还需要 Track，下一次磁通偏置和驱动频率应该取多少？

这部分主要由两处共同完成：

- frequency_control.py（`sqc/calibration/frequency_control.py` 源码约第 1 行）：纯控制算法；
- frequency_runtime.py（`sqc/workflows/frequency_runtime.py` 源码约第 433 行）：把控制算法接入逐轮测量。

完整链路是：

```text
上一轮 TrackSnapshot
        ↓
DampedSecantTracker.propose()
        ↓
TrackProposal：下一偏置和局部灵敏度
        ↓
Runtime 预测下一驱动和探针失谐
        ↓
TrackFrequency
        ↓
Backend 测量
        ↓
MeasurementSucceeded
        ↓
StateMachine.handle() 决定状态转移
        ↓
Runtime _post_process()
        ↓
DampedSecantTracker.accept()
        ↓
新的 TrackSnapshot
```

---

### 4.1 先建立跨轮控制器的记忆模型

#### 4.1 为什么控制器要单独拆出来

旧实现中的梯度闭环一次调用就会从开始运行到停止，不适合：

- 每轮检查用户中断；
- 每轮检查预算；
- 每轮保存 checkpoint；
- Track 中途进入 Verify；
- Verify 失败后返回 Track；
- 测量失败时进入 Reacquire。

因此控制器被拆成逐步接口：

```python
proposal = tracker.propose(snapshot)
estimate = measure(proposal.V_next)
result = tracker.accept(snapshot, estimate, proposal)
snapshot = result.snapshot
```

控制器自身是纯计算对象：

```text
不运行 Ramsey
不调用 QuTiP
不设置磁通偏置
不决定状态转换
不写 journal
```

它只回答：

> 根据已有两个测量点，建议下一偏置在哪里？

这种拆分也让旧 `SinglePointFrequencyCalibration` 和新状态机能够复用同一套偏置更新算法。

---

#### 4.2 先认识四个数据对象

##### 1. `FrequencyEstimate`

FrequencyEstimate（`sqc/calibration/frequency_control.py` 源码约第 26 行） 是控制器能够理解的测量结果：

```python
FrequencyEstimate(
    frequency=...,
    uncertainty=...,
    valid=...,
    ambiguous=...,
    method=...,
    shots=...,
    elapsed_time=...,
    diagnostics=...,
)
```

它和状态机的 `MeasurementSucceeded` 很相似，但层次不同：

```text
MeasurementSucceeded
    是 Command/Event 协议中的外部事件

FrequencyEstimate
    是纯控制器使用的测量值
```

Runtime 负责把前者转换为后者。

目前 tracker 主要使用：

```python
frequency
uncertainty
```

`valid` 与 `ambiguous` 已经由 Runtime 和状态机提前检查，不会把无效事件传给 `accept()`。

---

##### 2. `TrackSnapshot`

TrackSnapshot（`sqc/calibration/frequency_control.py` 源码约第 64 行） 是控制器的完整记忆。

它不是只有当前点，而是同时保存当前点和前一点：

```python
V
V_prev
f
f_prev
e
e_prev
```

其中：

$$
e=f-f_{\mathrm{target}}.
$$

这些字段组成割线所需的两个点：

$$
(V_{n-1},e_{n-1}),
\qquad
(V_n,e_n).
$$

另外还保存：

```python
best_V
best_abs_e
streak
n_iter
s_hat
uncertainty
uncertainty_prev
```

含义分别是：

- `best_V`：历史上目标残差最小的偏置；
- `best_abs_e`：对应的最小绝对残差；
- `streak`：连续进入 tracker 内部容差的次数；
- `n_iter`：种子以后已经接受的测量次数；
- `s_hat`：最近一次割线灵敏度；
- 当前和前一次频率不确定度。

整个 Snapshot 都能序列化，所以 checkpoint 恢复后不必重新建立割线历史。

---

##### 3. `TrackProposal`

TrackProposal（`sqc/calibration/frequency_control.py` 源码约第 89 行） 表示控制器提出的下一步建议：

```python
TrackProposal(
    V_next=...,
    step=...,
    s_hat=...,
    converged=...,
    diagnostics=...,
)
```

它只决定偏置，不决定驱动频率。类注释明确规定：

> 调用者必须根据 Proposal 解析驱动频率，并负责执行测量。

所以：

```text
Tracker 决定下一偏置
Runtime 决定下一驱动
Backend 执行测量
```

---

##### 4. `TrackStepResult`

TrackStepResult（`sqc/calibration/frequency_control.py` 源码约第 104 行） 是接受测量后的结果：

```python
TrackStepResult(
    snapshot=new_snapshot,
    is_best=...,
    stopped_by=...,
)
```

`stopped_by` 可能是：

```text
tolerance
max_iter
predicate
running
```

但在新状态机架构里，是否离开 Track 最终由状态机决定，而不是单靠这个字段。

---

#### 4.3 `initialize()`：第一份 Snapshot 从哪里来

一个容易误解的地方是：

> Tracker 的初始 Snapshot 默认不是直接由 Acquire 创建，而是由第一次成功的 Track 测量创建。

进入 Track 后，状态机先生成一条基础 `TrackFrequency`。Backend 完成这次测量，状态机处理事件；随后 Runtime 的 `_post_process()` 发现还没有 tracker snapshot，于是执行：

```python
seed = FrequencyEstimate(
    frequency=event.frequency,
    uncertainty=event.uncertainty,
    ...
)

snapshot = tracker.initialize(
    seed,
    V_seed=event.applied_bias,
)
```

initialize()（`sqc/calibration/frequency_control.py` 源码约第 190 行） 计算：

$$
e_0=f_0-f_{\mathrm{target}},
$$

并创建：

```python
TrackSnapshot(
    V=V_seed,
    V_prev=V_seed,
    f=f_0,
    f_prev=None,
    e=e_0,
    e_prev=None,
    best_V=V_seed,
    best_abs_e=abs(e_0),
    n_iter=0,
    s_hat=None,
)
```

此时只有一个测量点，因此还不能计算割线灵敏度。

---

#### 4.4 为什么第一次只能走探索步

下一轮 Runtime 调用：

```python
proposal = tracker.propose(snapshot)
```

由于：

```python
snapshot.n_iter == 0
```

控制器没有两个不同偏置上的结果，只能使用固定探索步：

```python
step = first_bias_step * chi_phi * sign(e)
V_next = V - step
```

其中

$$
\chi_\Phi
=
\operatorname{sign}\!\left(\frac{\partial f}{\partial V}\right)
\in\{-1,+1\}
$$

由 `expected_sensitivity_sign` 给出。`sign(e)` 回答“频率位于目标哪一侧”，
$\chi_\Phi$ 回答“改变偏置后频率往哪一侧走”；只有两者结合，才能确定偏置动作方向。

假设：

```text
当前偏置 V = 0.05 Φ0
目标残差 e > 0
first_bias_step = 0.01 Φ0
expected_sensitivity_sign = +1
```

则：

```text
step   = +0.01 Φ0
V_next = 0.05 - 0.01 = 0.04 Φ0
```

如果残差为负：

```text
step   = -0.01 Φ0
V_next = 0.05 - (-0.01) = 0.06 Φ0
```

如果当前工作点位于负斜率分支，即 `expected_sensitivity_sign=-1`，上述偏置动作会反转。
参数为 `None` 时保留旧的隐式正方向约定，用于兼容既有调用；论文协议和硬件运行应显式
配置已知分支方向。

如果这次探索没有产生可用梯度，控制器会继续复用带 $\chi_\Phi$ 的探索步，而不会凭空
编造灵敏度。

---

### 4.2 从两个工作点得到下一偏置

#### 4.5 割线灵敏度怎样得到

有两个不同偏置上的测量后：

$$
(V_{n-1},e_{n-1}),
\qquad
(V_n,e_n),
$$

计算：

$$
\Delta e=e_n-e_{n-1},
$$

$$
\Delta V=V_n-V_{n-1}.
$$

局部频率—偏置灵敏度为：

$$
\hat s_n
=
\frac{\Delta e}{\Delta V}.
$$

因为目标频率是常数：

$$
e_n-e_{n-1}
=
f_n-f_{n-1},
$$

所以它也等价于：

$$
\hat s_n
=
\frac{f_n-f_{n-1}}{V_n-V_{n-1}}
\approx
\frac{\partial f}{\partial V}.
$$

代码中还计算逆梯度：

$$
\frac{\Delta V}{\Delta e},
$$

用于直接计算偏置修正。

---

#### 4.6 阻尼割线步怎样计算

若局部线性近似为：

$$
e(V)\approx e_n+\hat s_n(V-V_n),
$$

令下一点残差为零，可得普通割线根：

$$
V_\star
=
V_n-\frac{e_n}{\hat s_n}.
$$

代码没有一步跳到根，而是乘阻尼系数：

$$
\Delta V_{\mathrm{control}}
=
\lambda\frac{e_n}{\hat s_n},
$$

$$
V_{n+1}
=
V_n-\Delta V_{\mathrm{control}},
$$

其中：

$$
0<\lambda\leq1.
$$

对应实现：

```python
grad_inv = dV / de
step = damping * e * grad_inv
V_next = snapshot.V - step
```

阻尼的作用是降低：

- 局部非线性造成的过冲；
- 有噪声割线造成的大步跳变；
- 接近 fold 时灵敏度失真带来的风险。

默认：

```python
damping = 0.8
```

即只走预计根位置的 80%。

在真正发出下一条测量命令前，受约束模式还要求：

$$
\operatorname{sign}(\hat s_n)=\chi_\Phi.
$$

若割线幅值虽然位于 `[S_min, S_max]`，符号却与已知分支方向不一致，proposal 保持
当前偏置并标记 `sensitivity_direction_mismatch`。Runtime 把它写入
`tracker_spec.sensitivity_direction_valid=False`，状态机的测前 guard 随即以
`SENSITIVITY_INVALID` 进入 Reacquire。错误割线因此不会先驱动硬件、再等待测后发现。

---

#### 4.7 为什么还要限制最大偏置步长

计算完成后：

```python
step = np.clip(
    step,
    -max_bias_step,
    max_bias_step,
)
```

即：

$$
|\Delta V_{\mathrm{control}}|
\leq
\Delta V_{\max}.
$$

即使割线给出很大的修正，也不会单步超过配置上限。

这在以下情况下尤其重要：

- 两次频率差很小，使逆梯度异常大；
- 测量噪声改变割线符号；
- 工作点靠近甜点，灵敏度很小；
- 局部线性近似已经开始失效。

此外，tracker 还支持：

```python
V_lo
V_hi
```

对最终 `V_next` 施加硬偏置边界。新 Runtime 当前主要使用 `max_bias_step`，没有默认给 tracker 传入全局 `V_lo/V_hi`。

---

#### 4.8 零梯度时为什么退回探索步

如果：

$$
|\Delta e|\approx0,
$$

则：

$$
\frac{\Delta V}{\Delta e}
$$

不可用。

代码不会除以接近零的数，而是重新使用：

```python
first_bias_step * chi_phi * sign(e)
```

这相当于说：

> 当前两点没有提供可靠的局部斜率信息，再做一个受限探索点，而不是用不稳定的割线外推。

但要注意：连续出现零梯度可能说明：

- 测量分辨率不足；
- 工作点位于低灵敏度区；
- 实际执行偏置没有变化；
- 响应进入平台区；
- 短脉冲反演异常。

状态机的 `S_min/S_max`、灵敏度符号和其他 Track guard 会进一步决定是否应该 Reacquire。

---

#### 4.9 一个完整的割线数值例子

为了直观，下面用 MHz 和 \(\Phi_0\) 表示。

目标频率：

$$
f_{\mathrm{target}}=5.000\ \mathrm{GHz}.
$$

第一次 Track 测量：

```text
V0 = 0.050 Φ0
f0 = 5.100 GHz
e0 = +100 MHz
```

初始化 Snapshot：

```text
V = V_prev = 0.050
e = +100 MHz
e_prev = None
n_iter = 0
```

由于没有梯度，第一次 Proposal 使用探索步：

```text
first_bias_step = 0.010 Φ0
V1 = 0.050 - 0.010 = 0.040 Φ0
```

在 \(V_1\) 测得：

```text
f1 = 5.050 GHz
e1 = +50 MHz
```

现在：

$$
\Delta e
=
50-100
=
-50\ \mathrm{MHz},
$$

$$
\Delta V
=
0.040-0.050
=
-0.010\Phi_0.
$$

因此：

$$
\hat s
=
\frac{-50}{-0.010}
=
5000\ \mathrm{MHz}/\Phi_0.
$$

若阻尼：

$$
\lambda=0.8,
$$

下一修正步为：

$$
\Delta V_{\mathrm{control}}
=
0.8\frac{50}{5000}
=
0.008\Phi_0.
$$

所以：

$$
V_2
=
0.040-0.008
=
0.032\Phi_0.
$$

这就是 `propose()` 的完整计算。

---

### 4.3 偏置更新与 drive tracking 怎样协同

#### 4.10 drive tracking 为什么与偏置建议分开

`TrackProposal` 只给出：

```text
V_next
s_hat
```

下一驱动由 Runtime 调用：

```python
DampedSecantTracker.resolve_track_drive(
    V=proposal.V_next,
    V_prev=snapshot.V,
    f_prev=snapshot.f,
    s_hat=proposal.s_hat,
)
```

计算：

$$
\omega_{d,n+1}
=
\hat f_n
+
\hat s_n(V_{n+1}-V_n).
$$

它不是直接把驱动设为最终目标，而是预测新偏置上的量子比特频率。

延续上例：

```text
当前频率 f1 = 5.050 GHz
当前偏置 V1 = 0.040 Φ0
下一偏置 V2 = 0.032 Φ0
灵敏度 s = 5000 MHz/Φ0
```

则：

$$
f_{\mathrm{pred},2}
=
5.050\ \mathrm{GHz}
+
5000\frac{\mathrm{MHz}}{\Phi_0}
(0.032-0.040)\Phi_0,
$$

$$
f_{\mathrm{pred},2}
=
5.010\ \mathrm{GHz}.
$$

所以新驱动设为：

$$
\omega_{d,2}/(2\pi)=5.010\ \mathrm{GHz}.
$$

若真实新频率是 5.012 GHz，则短脉冲实际只需测量：

$$
\Delta_2/(2\pi)
=
5.012-5.010
=
2\ \mathrm{MHz}.
$$

虽然相对于最终目标仍差：

$$
r_2/(2\pi)=12\ \mathrm{MHz},
$$

probe 只看见 2 MHz 的局部失谐。

这正是：

$$
r\neq\Delta
$$

以及 drive tracking 能让局部测量跨越较大目标区间的原因。

---

#### 4.11 没有割线时 drive 怎样处理

如果：

```python
s_hat is None
```

`resolve_track_drive()` 返回上一频率：

```python
omega_d = f_prev
```

这意味着首次探索步假设新偏置频率不会偏离上一频率太远。

但这只是操作策略，不是完整的测前安全证明。因为没有灵敏度时，无法严格预测偏置步会造成多大失谐。

所以 Runtime 还需要单独处理 blind step 的预测守卫。

---

#### 4.12 测前 prediction guard

Runtime 在真正执行下一 Track 测量之前，要估计：

$$
|\Delta_{\mathrm{pred}}|
+
z\sigma_{\Delta,\mathrm{pred}}
\leq
\Delta_{\mathrm{val}}.
$$

有两种情况。

##### 情况一：还没有割线灵敏度

如果 `proposal.s_hat is None`，默认无法完成严格预测：

```python
predicted_detuning = None
prediction_guard_complete = False
prediction_source = "unavailable"
```

如果用户配置：

```python
blind_step_detuning_bound
```

则可使用外部给定的保守上界：

```python
predicted_detuning = blind_step_detuning_bound
predicted_detuning_uncertainty = 0.0
prediction_guard_complete = True
prediction_source = "configured_blind_step_bound"
```

这里的零不确定度不表示实验无误差，而表示这个字段本身被解释为已包含风险的保守上界。

如果：

```python
require_complete_prediction_guard=True
```

但又没有 `blind_step_detuning_bound`，首次 blind step 会在测前被拒绝。严格硬件配置必须同时提供可辩护的首次步失谐上界，或者采用其他安全 bootstrap 策略。

##### 情况二：已经有割线

有两个测量点后，Runtime 传播频率不确定度。

首先估计灵敏度不确定度：

$$
\sigma_s
\approx
\frac{
\sqrt{\sigma_{f,n}^2+\sigma_{f,n-1}^2}
}{
|V_n-V_{n-1}|
}.
$$

因为 drive 设置为预测的新频率，预测失谐均值记为：

$$
\widehat\Delta_{\mathrm{pred}}=0.
$$

预测不确定度近似为：

$$
\sigma_{\Delta,\mathrm{pred}}
=
\sqrt{
\sigma_{f,n}^2
+
\left[
(V_{n+1}-V_n)\sigma_s
\right]^2
}.
$$

然后状态机 `validate_command()` 检查：

$$
z\sigma_{\Delta,\mathrm{pred}}
\leq
\Delta_{\mathrm{val}},
$$

以及额外 `linear_range` 和 `guard_margin`。

这个传播是清晰的一阶近似，当前没有建模相邻频率估计之间的相关性，也不是完整的非线性置信区间。

---

#### 4.13 Runtime 为什么要替换 pending command

状态机的 `_dispatch_command()` 只产生基础命令：

```python
TrackFrequency(
    command_id=固定ID,
    candidate_bias=当前偏置,
    predicted_drive=当前频率估计,
)
```

Runtime 调用 tracker 后得到更完整的命令：

```python
TrackFrequency(
    command_id=同一个固定ID,
    candidate_bias=proposal.V_next,
    predicted_drive=预测驱动,
    predicted_detuning=...,
    tracker_spec={
        "s_hat": ...,
        "step": ...,
        "prediction_source": ...,
    },
)
```

随后执行：

```python
machine.replace_pending_command(enriched_command)
```

为什么必须替换状态机内部的 pending command？

因为后续需要：

- 状态机比较实际执行值与真正命令值；
- checkpoint 保存实际准备执行的命令；
- 恢复后重放完全相同的偏置和驱动；
- journal 精确记录控制建议。

替换时必须保留原 `command_id`，否则同一轮会被错误解释成另一条命令。

---

### 4.4 新证据怎样更新控制器，而不越权决定状态

#### 4.14 测量后 `accept()` 做什么

Backend 完成测量后，调用顺序是：

```text
machine.handle(event)
journal 记录
runtime._post_process(cmd, event)
```

只有以下条件满足时才更新 tracker：

```text
cmd 是 TrackFrequency
event 是 MeasurementSucceeded
event.valid=True
event.ambiguous=False
状态机没有因该结果进入 Reacquire
```

这意味着越界或失效测量不会污染割线历史。

已有 Snapshot 时，Runtime 从实际命令重建 Proposal：

```python
TrackProposal(
    V_next=cmd.candidate_bias,
    step=cmd.tracker_spec["step"],
    s_hat=cmd.tracker_spec["s_hat"],
)
```

再将事件转换为：

```python
FrequencyEstimate(
    frequency=event.frequency,
    uncertainty=event.uncertainty,
)
```

然后调用：

```python
result = tracker.accept(
    old_snapshot,
    estimate,
    proposal,
)
```

---

#### 4.15 `accept()` 如何滚动历史

假设旧 Snapshot 是：

```text
当前点     (Vn, fn, en, σn)
前一点     (Vn−1, fn−1, en−1, σn−1)
```

接受新测量后变为：

```text
当前点     (Vn+1, fn+1, en+1, σn+1)
前一点     (Vn, fn, en, σn)
```

对应代码：

```python
V=V_next
V_prev=snapshot.V
f=f_next
f_prev=snapshot.f
e=e_next
e_prev=snapshot.e
uncertainty=estimate.uncertainty
uncertainty_prev=snapshot.uncertainty
```

它同时更新：

$$
e_{n+1}=f_{n+1}-f_{\mathrm{target}}.
$$

`n_iter` 增加一，表示又接受了一次由 Proposal 发起的测量。

---

#### 4.16 best point 有什么用

每次 `accept()` 比较：

$$
|e_{n+1}|
<
\texttt{best\_abs\_e}.
$$

若更好，则更新：

```python
best_V = V_next
best_abs_e = abs(e_next)
is_best = True
```

它为旧 batch 流程、诊断和失败回退保留最佳历史工作点。

但当前状态机不会因为 `best_V` 更好就自动回写物理偏置。状态机的 `_candidate_bias` 仍代表最近实际测量偏置。

所以：

```text
best_V 是控制器历史信息
candidate_bias 是状态机当前物理候选
```

不能把记录中的最佳点理解成设备已经自动回到该偏置。

---

#### 4.17 tracker 的 streak 与 Verify streak 不同

`accept()` 会计算：

```python
streak = (
    snapshot.streak + 1
    if abs(e_next) <= epsilon_f
    else 0
)
```

这是 tracker 内部的目标残差 streak。

它不能替代状态机 Verify 的 `verify_streak`，因为：

- tracker streak 来自 Track 的同一局部估计器；
- Verify streak 来自冻结偏置下的独立 Ramsey；
- tracker 使用纯残差 `abs(e)`；
-状态机候选和 Verify 还会加入不确定度上界。

因此：

```text
TrackSnapshot.streak
    控制算法内部的连续容差记录

FrequencyStateMachine._verify_streak
    独立验证连续通过记录
```

两者物理证据不同。

---

#### 4.18 为什么 Runtime 把 tracker 容差设为 `epsilon_final`

Runtime 初始化 tracker 时：

```python
DampedSecantTracker(
    converge_streak=1,
    max_iter=999,
    epsilon_f=config.epsilon_final,
)
```

看起来奇怪，因为 Track 进入 Verify 的门槛是：

```python
epsilon_enter
```

为什么 tracker 不用 `epsilon_enter`？

因为候选路由属于状态机：

$$
|r|+z\sigma_r
\leq
\epsilon_{\mathrm{enter}}
\quad\Rightarrow\quad
\text{Track → Verify}.
$$

tracker 不应自己宣布候选成功。

更重要的是，Verify 可能测得：

$$
\epsilon_{\mathrm{final}}
<
|r_{\mathrm{ver}}|
<
\epsilon_{\mathrm{enter}}.
$$

这时 Verify 会返回 Track，要求进一步修正。

如果 tracker 的内部容差也是 `epsilon_enter`，它会认为当前点已经收敛，提出零步长；Track 再次测量同一点后又进入 Verify，形成：

```text
Verify → Track → Verify → Track
```

所以 tracker 使用更严格的 `epsilon_final`，保证可靠 Verify miss 返回 Track 后仍会提出真实修正。

---

#### 4.19 为什么 tracker 的 `max_iter=999`

`DampedSecantTracker` 原本支持自己的：

```python
max_iter
stop_predicate
converge_streak
```

这些能力是为了兼容旧 batch 闭环。

新状态机中，全局停止和预算由以下部分管理：

```text
FrequencyCalibrationConfig.max_commands
Budget
Track guard
Verify 路由
Reacquire
SafeStop
```

因此 Runtime 将 tracker 的 `max_iter` 设得很大，避免纯控制器提前抢走状态机的生命周期决策权。

这体现了一个架构原则：

> Tracker 决定下一步怎么走；状态机决定还允不允许继续走。

---

#### 4.20 `proposal.converged` 在新架构中的含义

`propose()` 在以下情况下返回：

```python
converged=True
```

- tracker streak 达到内部要求；
- tracker `n_iter` 达到 `max_iter`。

但在状态机正常路径中：

- 达到 `epsilon_enter` 的 Track 结果已经被状态机送入 Verify；
- Verify miss 返回 Track 时 tracker 使用更严格的 `epsilon_final`；
- `max_iter=999`，实际预算由状态机控制。

因此正常情况下，状态机不会依赖 tracker 的 `converged=True` 宣布标定成功。

Runtime 只有在 `proposal.converged=False` 时才用 Proposal 丰富命令；否则保留状态机基础 Track 命令。最终是否进入 Verify 或 Lock仍由测量事件和状态机判断。

---

### 4.5 第四遍收束：把控制器放回完整 Track 循环

#### 4.21 Verify 返回 Track 时怎样保留控制连续性

在第二遍已经看到，Verify 可靠但未达最终容差时：

```text
Verify → Track
```

状态机会更新现有 tracker snapshot：

```python
tracker_snapshot["f"] = f_verify
tracker_snapshot["e"] = f_verify - f_target
tracker_snapshot["uncertainty"] = sigma_verify
tracker_snapshot["streak"] = 0
```

但保留：

```text
V_prev
f_prev
e_prev
```

以及原割线历史。

这样下一轮 `propose()` 能把独立 Verify 的频率作为当前点，同时利用之前 Track 的偏置差与频率差继续计算修正。

这里有一个细节：Verify 测量发生在同一个冻结候选偏置上，因此当前 `V` 不变，只更新了该点的频率证据。它相当于用更可信的 Ramsey 结果替换当前局部估计，而不是制造一个新的偏置点。

---

#### 4.22 一次完整 Track 控制循环

把第四遍串起来：

```text
1. 第一次成功 Track 测量
   → Runtime initialize()
   → TrackSnapshot(V0, f0, e0, n_iter=0)

2. 下一轮 Runtime _enrich_command()
   → tracker.propose(snapshot)
   → 无割线，使用 first_bias_step
   → TrackProposal(V1, s_hat=None)

3. Runtime 解析 drive
   → 因 s_hat=None，先使用 f0
   → 构造/检查 blind-step prediction guard

4. Runtime 用相同 command_id 替换 pending command

5. Backend 在 V1、ωd,1 上测量
   → 返回 f1、probe_detuning、uncertainty

6. StateMachine.handle(event)
   → 先判断 probe 是否有效
   → 再决定 Track / Verify / Reacquire

7. 若仍可保留 Track 历史
   → Runtime _post_process()
   → tracker.accept()
   → Snapshot 滚动到 (V1, f1, e1)

8. 再下一轮 propose()
   → 用 (V0,e0)、(V1,e1) 计算 s_hat
   → 生成阻尼割线偏置 V2

9. Runtime resolve_track_drive()
   → 预测 V2 上的频率
   → 使下一次 probe detuning 接近零

10. 重复，直到状态机进入 Verify 或恢复/停止状态
```

---

#### 4.23 控制器与其他层的责任边界

读完后应能画出：

```text
DampedSecantTracker
    输入：历史偏置、频率、残差
    输出：下一偏置和割线灵敏度
    不做 I/O，不决定状态

Runtime
    调用 propose/accept
    解析 drive
    传播预测不确定度
    丰富并保存 pending command

Backend
    实际执行偏置和驱动
    返回频率、detuning 与测量不确定度

FrequencyStateMachine
    检查局部有效性
    判断候选进入
    决定 Track / Verify / Reacquire / SafeStop
```

第四遍最重要的结论是：

1. `TrackSnapshot` 是可 checkpoint 的控制器完整记忆；
2. 偏置更新使用阻尼割线，drive tracking 使用同一割线预测下一频率；
3. tracker 内部收敛不等于独立 Verify 成功；
4. 测前 prediction guard 与测后 probe-detuning guard 是两道不同安全防线；
5. 状态机拥有协议决策权，tracker 只拥有“下一步建议权”。

## 第五遍：Runtime 如何组织持续运行

```{figure} frequency_state_machine_pass5.svg
:alt: 第五遍高亮运行时预算调度和持久化路径
:width: 100%
:align: center

第五遍的激活路径：预算和计时进入 Issue Unit，运行结果经 Commit 写回，同时追加到 Journal Memory；checkpoint 允许恢复已提交上下文。
```

在图中，`BR · Budget` 是动态计数器，Config 中的上限属于 `CR · Config ROM`；**Issue Unit** 对应 Runtime 的测前取消、预算预留、调度和 identity 检查；**Journal Memory** 对应增长型历史、成本账本与 checkpoint，而不是普通“日志寄存器”。这一区分正好解释 Config budget 与运行时 `Budget` 对象为何不是同一个东西。

第五遍打开的是整个系统的“运行中枢”：

```python
result = runtime.run()
```

前四遍分别解释了：

```text
状态机：
    决定当前需要什么证据，以及收到证据后转到哪里

Backend：
    把 Command 变成 Ramsey 或 transient 测量

Track 控制器：
    根据历史测量提出下一偏置和 drive
```

但这些部件本身不会自动形成一个持续运行的系统。真正负责把它们串起来的是：

`FrequencyCalibrationRuntime`

主要阅读文件是：

- `sqc/workflows/frequency_runtime.py`
- `sqc/workflows/frequency_state_machine.py` 中的 `Budget`、`snapshot()`、`restore()`
- `tests/unit/test_frequency_runtime.py`

完整运行关系是：

```text
创建或恢复状态机
        ↓
取得下一条 Command
        ↓
补全 Track 控制参数
        ↓
测前有效性与预算检查
        ↓
检查取消信号
        ↓
必要时等待监测周期
        ↓
Executor 执行测量
        ↓
状态机处理 Event
        ↓
记录 journal 和成本
        ↓
更新 Track 控制器
        ↓
保存 checkpoint
        ↓
进入下一轮
```

### 5.1 Runtime 先把各层组织成一轮可执行生命周期

#### 5.1 Runtime 解决的不是测量问题，而是生命周期问题

状态机能够回答：

> 当前处于什么状态，下一条命令是什么？

Backend 能够回答：

> 给我一条命令，我怎样完成测量？

Track 控制器能够回答：

> 根据前两个工作点，下一偏置应该是多少？

但持续运行还需要解决另一组问题：

```text
什么时候正式开始？
一轮操作按什么顺序执行？
预算不足时，如何保证昂贵测量尚未开始？
用户在等待期间请求停止怎么办？
用户在测量期间请求停止怎么办？
进入 SafeStop 后，安全偏置是否真的已经施加？
运行中断后，如何从同一条 pending command 恢复？
每次测量和状态转换怎样保存？
有限-shot随机状态怎样保持可复现？
```

这些都不是某个物理测量算法内部应该处理的事情，因此集中在 Runtime。

可以把 Runtime 理解成控制系统的操作系统：

```text
状态机       提供协议决策
控制器       提供数值建议
Backend      提供测量能力
Runtime      提供调度、资源管理、中断和持久化
```

#### 5.2 Runtime 初始化时建立了什么

`FrequencyCalibrationRuntime` 的用户输入主要是：

```python
FrequencyCalibrationRuntime(
    qubit=qubit,
    f_target=f_target,
    config=config,
    executor=executor,
    cancellation_token=token,
    checkpoint_directory=directory,
)
```

其中：

- `qubit`：被测对象；
- `f_target`：目标角频率；
- `config`：协议阈值、时间间隔和预算上限；
- `executor`：实际测量后端；
- `cancellation_token`：线程安全的取消信号；
- `checkpoint_directory`：中断时持久化运行状态的位置。

如果没有传入 `config`，Runtime 创建默认 `FrequencyCalibrationConfig`。

如果没有传入 `executor`，Runtime 创建：

```python
SQCExecutor(
    qubit=qubit,
    config=config,
    f_target=f_target,
)
```

如果没有传入取消信号，则创建自己的 `CancellationToken`。

但是初始化 Runtime 时并不立即创建并启动状态机。状态机的创建发生在第一次 `run()` 中。这使“配置 Runtime”和“真正开始占用运行预算”成为两个阶段。

Runtime 还保存一些不属于状态机本体的运行数据：

```text
_machine
_tracker
_journal
_cost_ledger
_checkpoints
_last_monitor_time
_interrupt_handled
_interrupt_request
_safe_hold_confirmed
```

这里已经体现出一条边界：

```text
_machine：
    科学协议的状态

Runtime 自身字段：
    执行环境、日志、取消与持久化状态
```

#### 5.3 run() 怎样启动或恢复一次运行

`run()` 首先检查 `_machine`。

如果还没有状态机：

```python
self._machine = FrequencyStateMachine(
    config=self.config,
    f_target=self.f_target,
)
self._machine.start()
```

同时清空本次运行的 journal、cost ledger 和 monitor 计时。

如果状态机已经存在，但仍处于 `READY`，只调用 `start()`。

如果状态机来自 `load_run()`，通常已经是 `RUNNING`、`CALIBRATED` 或某个终止状态，Runtime 不会重新构造 Acquire 起点，而是直接从恢复的状态继续。

因此有三种入口：

```text
新运行：
    没有 machine
    → 创建
    → start()
    → Acquire

已构造但未启动：
    machine.run_status == READY
    → start()

磁盘恢复：
    machine 已由 snapshot 恢复
    → 保留状态、预算、pending command 和历史
```

之后进入：

```python
self._run_event_loop()
```

循环退出后，`run()` 汇总结果：

```text
最终状态与 RunStatus
最终频率和候选偏置
transition log
journal
cost ledger
成本汇总
命令数和耗时
SafeHold 是否确认
Lock hold 统计
中断信息
Backend provenance
是否真正发生 Verify → Lock
```

特别值得注意的是：

```python
completed_verified_lock
```

并不是简单检查最终 `run_status`，而是搜索转换历史中是否真实发生过：

```text
Verify → Lock
```

这能区分“完成了某种运行”与“确实通过独立验证进入 Lock”。

#### 5.4 主循环为什么允许 RUNNING 和 CALIBRATED

循环条件不是只有：

```python
run_status == RUNNING
```

而是：

```text
run_status 为 RUNNING 或 CALIBRATED

或者：

已经进入 SafeStop，
但 SafeHold 尚未获得确认，
且还允许重试
```

`CALIBRATED` 不表示 Runtime 立即退出。它表示：

> 已经通过 Verify 进入 Lock，但长期监测仍在运行。

因此正常生命周期是：

```text
READY
  ↓ start
RUNNING
  ↓ Verify → Lock
CALIBRATED
  ↓ Lock monitor 持续运行
COMPLETED / SAFE_STOPPED / FAILED
```

如果 `stop_after_lock_cycles=0`，Lock 可以长期持续，直到：

- 用户取消；
- 时间、命令、shots 或 solver 预算耗尽；
- 参考失效；
- 外部 interlock；
- 进程被外部终止。

测试必须配置有限的 Lock cycle、命令数或时间预算，否则同步 `run()` 可以按设计长期不返回。

#### 5.5 一轮 Runtime 的真实执行顺序

主循环的核心不是简单的：

```python
command = machine.next_command()
event = executor.execute(command)
machine.handle(event)
```

而是：

```text
1. 检查取消信号
2. next_command()
3. enrich Track command
4. 替换 pending command
5. 测前 guard
6. 估计成本并检查预算
7. 再检查取消信号
8. 必要时等待 monitor interval
9. 等待后再检查预算
10. 再检查取消信号
11. executor.execute()
12. 补充 state_version
13. machine.handle()
14. 写 journal 与 cost ledger
15. 检查测量期间到达的取消请求
16. 处理 SafeHold 确认
17. Track 后处理
18. 周期性内存 checkpoint
```

这个顺序不是任意的。每一步都在保护不同的边界。

#### 5.6 第一道边界：取消检查先于取科学命令

每轮最先执行：

```python
if self._process_cancel_request():
    continue
```

如果取消已经到达，Runtime 不再生成或执行新的科学测量，而是把取消信号转换为：

```python
CancelRequested(...)
```

交给状态机。

状态机收到该自发事件后：

```text
任意科学状态
    ↓ CancelRequested
SafeStop
run_status = SAFE_STOPPED
pending command 被清除
```

循环随后 `continue`，下一轮 `next_command()` 返回 `SafeHold`。

因此在运行开始前就调用：

```python
runtime.request_cancel()
runtime.run()
```

不会先执行 Acquire。测试明确验证：

```text
科学命令数 = 0
最后执行 SafeHold
safe_hold_confirmed = True
```

#### 5.7 Track 命令为何要在 Runtime 中补全

状态机生成的初始 `TrackFrequency` 只是一条基础命令。Runtime 在 `_enrich_command()` 中读取 tracker snapshot，调用：

```python
proposal = tracker.propose(snapshot)
```

再补全：

```text
candidate_bias
predicted_drive
predicted_detuning
predicted_detuning_uncertainty
prediction_guard_complete
tracker_spec
```

补全后必须执行：

```python
machine.replace_pending_command(cmd)
```

这样 checkpoint 保存的不是状态机最初生成的粗略命令，而是即将真正交给 Backend 的完整命令。

这一顺序形成：

```text
next_command()
    → 创建 command_id
    → 状态机保存 baseline pending command

_enrich_command()
    → 保留 command_id
    → 补全真实偏置和 drive

replace_pending_command()
    → checkpoint 保存真实执行参数
```

如果进程在补全后、测量前崩溃，恢复时能够重放同一个偏置、drive 和 command ID。

#### 5.8 第二道边界：测前命令有效性检查

补全 Track 命令后，Runtime 调用：

```python
preflight_failure = machine.validate_command(cmd)
```

当前主要检查 Track prediction guard。

如果命令预测已越出局部范围，或者配置要求完整预测但证据不足，Runtime 不调用 Backend，而是构造：

```python
MeasurementRejected(
    reason_code=preflight_failure,
    diagnostics={
        "stage": "premeasurement",
        ...
    },
)
```

再正常交给状态机处理并写入 journal。

这与测量后失败不同：

```text
测前失败：
    尚未执行物理测量
    来源是 command prediction guard

测后失败：
    Backend 已经返回证据
    来源是实际 probe detuning、置信度或硬件不一致
```

两者都通过 Event 进入状态机，而不是 Runtime 直接修改状态。

### 5.2 预算与时间：昂贵操作开始前必须先获准

#### 5.9 Config 中的预算与 Budget 对象有什么区别

这是之前讨论过但在 Runtime 中才能完整理解的问题。

`FrequencyCalibrationConfig` 保存的是政策上限：

```text
max_commands
max_wall_time
max_shots
max_solver_calls
max_reacquire_attempts
max_verify_attempts_per_episode
max_verify_shots
stop_after_lock_cycles
```

它回答：

> 这次运行最多允许使用多少资源？

`Budget` 保存的是动态计数：

```text
commands_issued
wall_time_start
total_shots
solver_calls
reacquire_attempts
verify_attempts_this_episode
verify_shots_this_episode
lock_cycles_completed
```

它回答：

> 这次运行已经使用了多少资源？

关系是：

```text
FrequencyCalibrationConfig
        上限和政策
             ↓ 比较
Budget
        当前累计值
```

Config 通常在一次运行中保持不变；Budget 会在每个事件后更新，并被写入 checkpoint。

#### 5.10 命令预算何时消耗

`next_command()` 在创建新的 pending command 时执行：

```python
self._budget.consume_command()
```

但如果已经存在 pending command，`next_command()` 原样返回，不重复增加命令数。

因此命令预算统计的是：

> 状态机新发出了多少个协议命令身份？

而不是：

> Runtime 调用了多少次 `next_command()`？

重放同一 pending command 不会再次扣除命令预算。

还有一个细节：结果中的 `n_commands` 根据 journal 中实际记录的 command 行计算，并且会包括执行过的 `SafeHold`；状态机的 `commands_issued` 则是在协议命令生成时累计。两者目的不同，不保证在所有中断和测前拒绝路径上完全相等。

#### 5.11 为什么预算检查分为“现有预算检查”和“测前成本检查”

状态机 `next_command()` 首先调用：

```python
budget.check(config)
```

它检查已经累计的：

```text
commands
shots
solver calls
Verify attempts/shots
wall time
Lock cycles
```

如果已经达到上限，直接进入 SafeStop，不再发出科学命令。

但假设当前只用了 900 shots，上限是 1000，而下一次 Verify 预计需要 200 shots。仅检查：

```text
900 < 1000
```

仍会错误地允许一条必然越预算的测量。

因此 Runtime 取得命令后，还调用 Backend：

```python
cost = executor.estimate_cost(cmd)
```

再执行：

```python
machine.reserve_execution_budget(
    estimated_shots=cost["shots"],
    estimated_solver_calls=cost["solver_calls"],
)
```

其判断类似：

```text
当前实际消耗 + 下一测量预计成本
    是否仍不超过上限
```

若不能容纳：

```text
清除 pending science command
进入 SafeStop
不执行 Backend 科学测量
下一条命令变成 SafeHold
```

测试验证了 `max_shots=1`、而测量预估成本为 2 时，Executor 收到的全部命令都只能是 `SafeHold`。

#### 5.12 预算预留为什么不立即扣除成本

`reserve_execution_budget()` 名字中有 reserve，但当前实现并不把预计 shots 和 solver calls直接加入 Budget。

它只回答：

> 这次测量预计能否放进剩余预算？

真实成本在事件返回后由状态机记录：

```python
budget.consume_shots(event.shots)
budget.consume_solver_calls(actual_solver_calls)
```

因此：

```text
测前：
    使用 estimate_cost 做 admission control
    不扣真实账

测后：
    使用 Event 中的实际 shots 和 solver calls
    更新 Budget
```

这样不会因为保守估计高于实际成本而永久多扣预算。

代价是：当前 Runtime 是同步单命令执行模型。若未来允许多个命令并发，仅做检查而不占用 reservation 会产生竞态；届时需要真正的未决预留账本。当前单线程顺序执行下不存在该问题。

#### 5.13 Verify 还有独立的 episode 预算

除全局 shots 预算外，Verify 还限制：

```text
max_verify_attempts_per_episode
max_verify_shots
```

一次 Verify episode 指从进入 Verify 开始，直到：

- 通过进入 Lock；
- 可靠未达标返回 Track；
- 歧义进入 Reacquire；
- 失败进入其他恢复路径。

每次 Verify 事件后：

```python
record_verify_attempt(event.shots)
```

离开本次 Verify episode 时：

```python
reset_verify_episode()
```

因此 Verify 不会因为连续低质量测量无限消耗 shots，也不会把上一轮 Verify episode 的次数错误带入下一次独立候选验证。

测前预算检查还专门判断：

```text
本 episode 已使用 Verify shots
+ 下一条 Verify 预计 shots
≤ max_verify_shots
```

#### 5.14 Lock 的 monitor interval 怎样调度

状态机在 Lock 中决定下一条命令是：

```text
VerifyFrequency  若周期性独立复核到期
MonitorFrequency 否则执行低成本监测
```

Runtime 负责 `MonitorFrequency` 的实际时间间隔。

若：

```python
config.monitor_interval > 0
```

且之前已经完成过一次 monitor，则计算：

```python
delay = monitor_interval - (
    current_time - last_monitor_time
)
```

如果仍需等待，不使用不可中断的 `sleep()`，而是：

```python
cancellation_token.wait(delay)
```

这很关键。普通 `sleep(60)` 会让用户取消后仍必须等待一分钟；Event wait 可以在取消信号到达时立即返回。

测试验证：

```text
monitor_interval = 60 s
等待期间请求取消
→ 不再执行第二次 Monitor
→ 转入 SafeStop
→ 执行 SafeHold
```

第一次进入 Lock 时 `_last_monitor_time=0`，因此第一条 monitor 通常立即执行。第二条开始才按与上一条 monitor 的时间差等待。

`monitor_interval=0` 表示每次循环都可立即监测，不代表关闭 monitor。

#### 5.15 周期性独立复核与 monitor interval 不属于同一个时钟

两种时间参数作用不同：

```text
monitor_interval：
    相邻低成本 MonitorFrequency 之间的间隔
    Runtime 负责等待

audit_interval：
    Lock 中相邻独立 Ramsey 复核之间的间隔
    状态机在 dispatch 时判断是否进入 Verify
```

代码仍使用兼容旧名 `audit_interval`，概念上更适合称为“周期性独立复核间隔”。

如果周期到期：

```text
Lock
  ↓ AUDIT_DUE
Verify
  ↓ VerifyFrequency
独立 Ramsey
```

它不会让 Monitor 直接获得 Verify 权限。

配置：

```python
require_periodic_audit=True
```

时，`audit_interval` 必须大于零，否则 Config 构造直接失败。

默认：

```python
audit_interval=0
```

表示关闭周期性复核。因此若论文方法要求定期独立 Ramsey，运行配置必须显式启用，不能仅依赖默认值。

### 5.3 中断与 SafeHold：停止协议和设备安全是两件事

#### 5.16 CancellationToken 是什么

`CancellationToken` 内部使用：

```python
threading.Event
threading.Lock
```

任何 UI、API 或信号处理线程都可以调用：

```python
runtime.request_cancel(
    reason="user_requested",
    source="api",
)
```

请求内容包括：

```text
reason
source
requested_at
```

它采用 first-request-wins：

```text
第一次 request() → True，保存请求
后续 request()   → False，不覆盖第一次原因
```

这样能保留真正首先触发中断的来源。例如 UI 和键盘几乎同时请求停止时，不会让后到请求覆盖先到的审计信息。

Token 本身不修改状态机。它只设置线程安全信号；Runtime 线程负责把它转换成 `CancelRequested` Event。

#### 5.17 Runtime 一轮检查几次取消信号

取消不是只在每轮末尾检查一次。

一轮中至少有这些检查点：

```text
A. 取得下一命令前

B. 命令补全和预算预留后

C. monitor 等待返回后

D. 正式进入 executor.execute() 前

E. executor.execute() 返回并处理事件后
```

因此在以下阶段请求取消，都能阻止下一条科学命令：

```text
两个命令之间
Track 命令补全后
预算检查后
Lock monitor 等待期间
当前同步测量完成后
```

但是它是协作式取消，不是抢占式取消。

#### 5.18 为什么同步 execute() 不能被立即中断

一旦进入：

```python
event = executor.execute(cmd)
```

Runtime 线程就被同步测量占用。Token 可以被其他线程设置，但 Runtime 无法在 Python 调用内部自动跳出。

因此如果用户在测量期间取消：

```text
1. 取消信号立即被设置
2. 当前 executor.execute() 继续运行
3. 测量返回 Event
4. 状态机先处理并记录该 Event
5. Runtime 检查取消信号
6. 转入 SafeStop
7. 不再执行下一条科学命令
```

测试中的 Executor 在第一次 Acquire 执行期间调用：

```python
runtime.request_cancel("during_measurement")
```

结果是：

```text
Acquire 完成一次
不会开始第二条科学命令
随后执行 SafeHold
```

所以当前中断延迟上界大致是：

> 当前同步命令的最长执行时间。

若硬件测量可能持续很久，要实现更细粒度中断，需要 Backend 自身支持：

- 分块执行；
- 硬件任务取消接口；
- 超时；
- 在内部扫描点之间检查 token；
- 异步 future/cancel。

这不是 Runtime 外层多检查几次标志就能解决的。

#### 5.19 取消请求怎样变成可审计事件

`_process_cancel_request()` 只处理一次取消信号：

```python
if interrupt_handled:
    return False
```

第一次发现请求后，它保存：

```text
取消前状态
取消前 state_version
当时的 pending command
取消原因、来源和时间
```

然后向状态机提交：

```python
CancelRequested(**request)
```

并在 journal 中增加一条没有 command 的自发事件：

```text
command = None
event = CancelRequested
state_before
state_after = SafeStop
interrupt_reason
interrupt_source
requested_at
interrupted_command
interrupted_command_id
transition_reason = CANCELLED
```

这里 `command=None` 很重要，因为取消不是某条科学命令的测量结果，而是外部自发事件。

若当时有 pending Acquire 或 Verify，它的类型和 ID仍被写入：

```text
interrupted_command
interrupted_command_id
```

方便以后判断取消发生在什么操作附近。

#### 5.20 KeyboardInterrupt 怎样接入同一体系

`run()` 外层捕获 `KeyboardInterrupt`。

默认：

```python
handle_keyboard_interrupt=True
```

此时不会直接让 Python 堆栈退出，而是转换为：

```python
request_cancel(
    reason="keyboard_interrupt",
    source="keyboard",
)
```

然后复用相同的：

```text
CancelRequested
→ SafeStop
→ SafeHold
```

流程。

这样 API 中断和 Ctrl+C 最终具有统一日志、状态转移和安全收尾语义。

若设置：

```python
handle_keyboard_interrupt=False
```

Runtime 会原样重新抛出 `KeyboardInterrupt`。这适合上层应用自己拥有更高层的异常处理策略，但也意味着该 Runtime 不保证自动完成 SafeHold。

#### 5.21 进入 SafeStop 为什么还不能立即退出循环

状态机进入 `SafeStop` 只表示协议决策：

> 不再执行科学测量，应进入安全保持。

它不证明硬件已经设置到安全偏置。

真正动作是下一条：

```python
SafeHold(
    bias=0.0,
    reason=...
)
```

Backend 执行后必须返回：

```python
SafeHoldApplied(
    applied_bias=...
)
```

Runtime 才设置：

```python
safe_hold_confirmed = True
```

因此：

```text
SafeStop：
    协议状态

SafeHold：
    请求硬件采取安全动作

SafeHoldApplied：
    硬件或仿真后端确认动作完成
```

终止状态不能替代设备确认。

#### 5.22 SafeHold 失败时怎样处理

循环专门允许在进入 SafeStop 后继续执行 SafeHold，直到：

- 收到 `SafeHoldApplied`；
- 或超过允许的技术重试次数。

当前循环条件使用 `safe_hold_attempts <= max_technical_retries`，语义相当于：

```text
一次初始尝试
+ 配置允许的重试机会
```

测试构造第一次 SafeHold 不确认、第二次确认的 Executor，验证 Runtime 确实再次发送 SafeHold，并以确认事件作为 journal 最后一条记录。

如果始终不确认，Runtime 最终可以退出，但：

```python
safe_hold_confirmed == False
```

调用方必须把它视为重要安全状态，而不是普通完成。

#### 5.23 中断 checkpoint 为什么保存两次

发现取消时，Runtime 立即：

```text
记录 CancelRequested
保存内存 snapshot
尝试 save_run(checkpoint_directory)
```

这第一次保存保证即使 SafeHold 本身失败或进程随后崩溃，也至少记录：

```text
取消已经被接受
状态机已经进入 SafeStop
当时中断了哪条命令
```

当收到 `SafeHoldApplied` 后，再次调用持久化。

第二次保存把：

```text
safe_hold_confirmed = True
```

写入 `result.json`。

因此磁盘记录可以区分：

```text
已经收到停止请求，但安全动作尚未确认

与

停止请求和安全动作都已确认
```

持久化采用 best effort。保存失败不会阻止 SafeHold 尝试；错误记录在：

```text
interrupt_checkpoint_saved
interrupt_checkpoint_error
```

安全动作优先级高于日志写盘。

### 5.4 运行证据：日志、成本、checkpoint 与 provenance

#### 5.24 journal 记录的是什么

每完成一条 command→event 周期，Runtime 写一条 journal。

基本身份字段包括：

```text
command
command_id
event
state_before
state_after
state_version
run_status
```

测量事件还记录：

```text
commanded_bias
commanded_drive

frequency
uncertainty
residual
valid
ambiguous
method
shots
elapsed_time

applied_bias
applied_drive
probe_detuning
probe_detuning_uncertainty

confidence_multiplier
diagnostics
transition_reason
```

对于 Monitor，还记录：

```text
hold_target_met
```

这使日志能够回答：

```text
想执行什么？
实际执行了什么？
测到了什么？
不确定度是多少？
probe 相对 drive 偏了多少？
状态从哪里转到哪里？
为什么转移？
```

`_json_safe()` 会把：

- NumPy 数组；
- NumPy 标量；
- Enum；
- tuple；

转换成 JSON 可序列化形式。测试直接执行：

```python
json.dumps(result["journal"])
```

验证完整 journal 可以序列化。

#### 5.25 transition log 与 journal 有什么区别

两者相关但不相同。

`transition_log` 是状态机的协议历史：

```text
version
from
to
reason
candidate_bias
f_hat
uncertainty
```

它只关心：

> 状态如何变化，为什么变化？

`journal` 是 Runtime 的执行历史：

```text
Command
Event
命令参数
测量参数
成本
状态前后
诊断信息
```

它关心：

> 这轮实际做了什么，得到了什么？

例如一条测量可能是 Track 自循环：

```text
transition log：
    Track → Track
    reason = TARGET_NOT_REACHED

journal：
    TrackFrequency
    candidate_bias = ...
    predicted_drive = ...
    actual frequency = ...
    probe_detuning = ...
    shots = ...
```

反过来，取消事件没有科学 Command，但仍会有状态转换和一条自发 journal 记录。

#### 5.26 cost ledger 为什么不直接等于 journal

Journal 已包含 shots 和 elapsed time，但 Runtime 仍维护独立 `cost_ledger`。

它把每条执行压缩为统一成本结构：

```text
command_id
command
state_before
state_after
event
success
shots
circuits
mesolve_calls
sesolve_calls
solver_calls
elapsed_time
```

这样生成总成本时不必重新解释每种 Event：

```python
cost_summary = {
    "shots": ...,
    "circuits": ...,
    "mesolve_calls": ...,
    "sesolve_calls": ...,
    "solver_calls": ...,
}
```

两者用途不同：

```text
journal：
    科学复现和协议审计

cost ledger：
    资源分析和性能统计
```

Budget 又与二者不同：

```text
Budget：
    在线决定还能否继续

cost ledger：
    离线说明实际花了多少

journal：
    说明为什么花、得到什么证据
```

#### 5.27 周期性内存 checkpoint 保存了什么

每当：

```python
len(journal) % 10 == 0
```

Runtime 将：

```python
machine.snapshot()
```

加入 `_checkpoints`。

`StateMachineSnapshot` 保存：

```text
协议状态和 RunStatus
state_version

candidate_bias 及其版本/来源
最新频率和不确定度

完整 tracker_snapshot

Verify streak
Lock monitor streak
Lock 进入时间与累计时间
last_audit_time
hold 统计

完整 Budget

transition log
last_event

pending command
pending command_id
technical retries
```

它是状态机可恢复状态的完整记录，而不只是“当前 State 枚举”。

但是 `_checkpoints` 默认只是 Runtime 内存中的列表。除非调用 `save_run()`，进程退出后这些周期 checkpoint 不会自动成为磁盘文件。

这一点必须明确区分：

```text
self._checkpoints.append(snapshot)
    → 内存 checkpoint

save_run(directory)
    → 磁盘持久化
```

#### 5.28 save_run() 在磁盘写什么

`save_run(directory)` 当前写出：

```text
config.json
commands.jsonl
cost-ledger.jsonl
transitions.jsonl
checkpoint.json
result.json
```

`config.json` 包含：

```text
FrequencyCalibrationConfig 全部字段
f_target
run_id
backend_provenance
可选 RNG state
```

`commands.jsonl` 是 journal。

`cost-ledger.jsonl` 是成本账本。

`transitions.jsonl` 是状态转移历史。

`checkpoint.json` 是最新即时生成的完整状态机 snapshot，不要求 `_checkpoints` 列表里恰好已有一份。

`result.json` 保存当前结果摘要和中断/SafeHold状态。

因此一个运行目录同时包含：

```text
政策
执行历史
状态历史
恢复状态
结果摘要
随机性来源
```

#### 5.29 Backend provenance 有什么作用

Runtime 调用：

```python
executor.provenance()
```

若 Backend 提供该接口，就把信息写入结果和配置。

有限-shot Executor 可以报告：

```text
backend 类型
master_seed
默认 shots_per_circuit
不同角色的 shots
assignment error
独立 RNG 角色
```

过程仿真 Backend 还可以报告：

```text
drift_per_command
jump_schedule
reference_loss_commands
inner backend
```

这样同一条状态历史可以追溯到具体证据生成模型。

如果自定义 Executor 没有 `provenance()`，至少记录类名。

#### 5.30 有限-shot RNG 状态怎样恢复

`FiniteShotSQCExecutor` 为不同角色维护独立 RNG。

如果 Executor 提供：

```python
rng_state()
restore_rng_state(state)
```

`save_run()` 会把当前 RNG 状态写入 `config.json`，`load_run()` 会在恢复 Runtime 后把状态交还给 Executor。

这使恢复后的下一次随机采样从保存点继续，而不是重新从初始 seed 开始。

仅保存 `master_seed` 不够，因为运行到中途时各 RNG 已消费不同数量的随机数。真正可重复恢复需要保存 bit-generator 当前状态。

但这种可复现仍针对仿真随机流。真实硬件测量无法通过恢复伪随机数状态重现同一物理噪声。

### 5.5 跨进程恢复：保存命令身份，也保存算法记忆

#### 5.31 pending command 为什么是恢复语义的核心

假设状态机已经生成：

```python
VerifyFrequency(
    command_id="abc",
    frozen_bias=0.03,
)
```

并写入 checkpoint，但进程在结果处理前退出。

恢复后如果直接重新从 Acquire 开始，会丢失：

- 已完成的 Track；
- 当前冻结候选；
- Verify streak；
- 预算；
- 命令身份。

正确恢复是：

```text
checkpoint 保存 pending command
        ↓
load_run() 反序列化
        ↓
machine.next_command()
        ↓
发现已有 pending command
        ↓
原样返回同一 command_id
```

测试明确构造 Verify pending command，保存并恢复，然后验证：

```text
恢复后 Executor 第一条收到 VerifyFrequency
command_id 与保存前相同
没有重新执行 Acquire
```

#### 5.32 这种恢复为什么属于 at-least-once，而不是 exactly-once

checkpoint 能证明：

> 这条命令已经发出或准备执行，但没有已处理的匹配 Event。

它不能总能证明：

> 外部硬件是否已经执行完成，只是结果尚未写回。

例如进程可能在以下时刻崩溃：

```text
硬件已经完成 Verify
    ↓
Event 尚未写入 checkpoint
    ↓
进程退出
```

恢复后会再次执行同一 pending command。

因此当前语义更接近：

```text
at-least-once execution
```

而不是：

```text
exactly-once execution
```

`command_id` 保证身份一致和事件去重，但不能让不支持幂等的物理硬件动作自动变成 exactly-once。

对于只读 Ramsey 测量，重复通常可以接受；对于带偏置执行的 Track，硬件适配层最好按绝对目标值设置偏置和 drive，而不是执行“再增加一次”的相对动作。绝对设定更容易实现幂等重放。

若未来需要 exactly-once，需要硬件任务系统支持：

- 查询 command ID 是否已执行；
- 持久化设备侧任务状态；
- 幂等命令；
- 或事务式执行/确认协议。

#### 5.33 state_version 在恢复和异步事件中做什么

Runtime 在执行命令前保存：

```python
version_before = machine.state_version
```

如果 Backend 返回 Event 的 `state_version` 仍为默认 0，Runtime 将其补为执行前版本。

状态机随后同时检查：

```text
command_id 是否匹配 pending command
state_version 是否对应预期状态版本
```

`command_id` 防止旧命令的结果被新命令接收。

`state_version` 防止即使命令身份看似合理，事件却属于已经变化的协议状态。

每次真实转换和自循环都会增加 `state_version`，所以：

```text
Track → Track
```

虽然枚举状态没变，也形成新的协议版本。

#### 5.34 Tracker 为什么在磁盘恢复后可以继续

`StateMachineSnapshot` 保存的是序列化后的 `_tracker_snapshot`，其中包含：

```text
V, V_prev
f, f_prev
e, e_prev
uncertainty, uncertainty_prev
s_hat
best_V, best_abs_e
streak
n_iter
```

恢复状态机时，这个字典重新放回 `_tracker_snapshot`。

Runtime 的 `_tracker` 控制器对象本身可能还是 `None`，但下一次进入 `_enrich_command()` 时会按 Config 重新创建一个等价的 `DampedSecantTracker`，然后用恢复的 snapshot 继续 `propose()`。

因此恢复分为：

```text
算法参数：
    从 Config 重建 tracker 对象

算法历史：
    从 tracker_snapshot 恢复
```

不需要序列化整个 Python 控制器实例。

#### 5.35 恢复后 wall-time 预算的语义

`Budget` 保存：

```python
wall_time_start
```

它是原运行开始时的绝对时间戳。

恢复后不会自动把 wall-time 重新从零开始。因此长时间暂停后再次恢复，原 wall-time 预算可能已经耗尽。

这代表当前政策是：

> wall-time 约束覆盖从首次启动开始的实际日历时间，包括停机间隔。

这未必适合所有实验。如果希望预算只统计 Runtime 真正活跃的时间，需要改成累计 active duration，并在 save/load 时记录暂停边界。

当前读代码时不能假设恢复会重置时间预算。

### 5.6 第五遍收束：用两条时间线检查 Runtime 边界

#### 5.36 一次正常 Track 轮次的完整 Runtime 时间线

把第五遍所有部件连起来：

```text
1. Runtime 检查 cancellation token
   → 没有取消

2. machine.next_command()
   → 生成 baseline TrackFrequency
   → 保存 pending command
   → commands_issued += 1

3. Runtime _enrich_command()
   → tracker.propose(snapshot)
   → 得到 V_next、step、s_hat
   → 预测 drive 与 detuning uncertainty

4. replace_pending_command()
   → 同一 command_id
   → pending payload 变成真实执行参数

5. validate_command()
   → prediction guard 通过

6. executor.estimate_cost()
   → 估计 shots / solver calls

7. reserve_execution_budget()
   → 剩余预算足够
   → 不预扣实际成本

8. Runtime 再检查取消信号

9. executor.execute(command)
   → 在 candidate_bias 和 predicted_drive 上测量
   → 返回 frequency、uncertainty、probe_detuning
   → 返回实际 solver calls

10. Runtime 补充 event.state_version

11. machine.handle(event)
    → 检查 applied bias/drive
    → 检查实际 probe detuning
    → 计算目标 residual
    → Track 自循环、进入 Verify 或 Reacquire
    → 实际 shots/solver calls 加入 Budget
    → 清除 pending command

12. Runtime _record_journal()
    → 保存命令、测量和转换原因

13. Runtime _record_cost()
    → 保存成本账本

14. Runtime 再检查取消信号
    → 防止测量期间的取消被遗漏

15. 若仍处于可接受 Track 路径
    → _post_process()
    → tracker.accept()
    → 滚动 TrackSnapshot

16. 每十条 journal
    → 保存内存 snapshot

17. 下一轮
```

#### 5.37 一次取消路径的完整时间线

如果取消发生在命令之间：

```text
1. UI 调用 runtime.request_cancel()
2. CancellationToken 保存第一条请求
3. Runtime 在下一检查点发现 requested=True
4. 创建 CancelRequested
5. 状态机清除 pending command并进入 SafeStop
6. journal 记录取消原因、来源和被中断命令
7. 保存取消 checkpoint
8. 下一轮 next_command() 返回 SafeHold
9. Backend 执行安全偏置
10. 返回 SafeHoldApplied
11. Runtime 设置 safe_hold_confirmed=True
12. 再保存一次 checkpoint
13. 退出
```

如果取消发生在同步测量期间：

```text
1. 当前 execute() 继续完成
2. 返回 Event
3. 状态机处理 Event并记录成本
4. Runtime 发现 cancellation token
5. 进入上述 CancelRequested → SafeHold 路径
6. 不再开始下一条科学命令
```

#### 5.38 Runtime 不负责什么

Runtime 虽然是调度中枢，但不应吸收所有逻辑。

它不负责：

```text
定义 Acquire/Track/Verify/Lock 的转换条件
    → 状态机负责

计算候选是否满足 epsilon_enter
    → 状态机 guard 负责

从布居反演频率
    → FrequencyMeasurement 负责

决定割线偏置步
    → DampedSecantTracker 负责

伪造 uncertainty 或 confidence
    → Backend 必须报告真实证据模型

直接把 monitor 漂移变成偏置更新
    → 必须先由状态机转入 Verify
```

Runtime 只负责：

```text
什么时候调用谁
以什么顺序调用
调用前是否还允许
调用后记录什么
如何响应外部停止
如何恢复运行
```

#### 5.39 第五遍读完后应形成的总体认识

前四遍的最小模型是：

```text
Command → Measurement → Event → Transition
```

第五遍以后，应扩展为：

```text
                    Config
                      ↓
Cancellation → Runtime event loop ← Budget
                      ↓
           next_command / enrich
                      ↓
        validate / reserve / wait
                      ↓
                  Executor
                      ↓
                    Event
                      ↓
             StateMachine.handle
                 ↙          ↘
        journal/cost       tracker.accept
                 \          /
                  checkpoint
                      ↓
                    resume
```

第五遍最重要的结论是：

1. Runtime 是生命周期协调器，不是测量算法或状态决策器。
2. 取消信号会在一轮中的多个安全边界检查，但无法强制抢占已经进入的同步 `execute()`。
3. Config 保存预算上限，`Budget` 保存动态计数；测前检查预计成本，测后记录实际成本。
4. 进入 SafeStop 不等于硬件已安全，必须执行并确认 `SafeHoldApplied`。
5. Journal、transition log、cost ledger 和 checkpoint 分别服务执行复现、协议历史、成本分析和恢复，不能互相替代。
6. checkpoint 保存 pending command 和 command ID，恢复不会重新从 Acquire 开始。
7. 当前命令恢复更接近 at-least-once；exactly-once 需要硬件侧幂等或任务查询支持。
8. 周期内存 checkpoint 不等于磁盘持久化；真正跨进程恢复需要 `save_run()`。
9. 有限-shot恢复还会保存 RNG 状态，从而延续仿真随机流，但这不等于重现实验噪声。
10. Lock 是持续运行状态，因此 Runtime 必须有取消、时间预算或有限 Lock cycle 才能保证同步调用最终返回。

## 第六遍：异常恢复与测试证据

```{figure} frequency_state_machine_pass6.svg
:alt: 第六遍高亮异常恢复与证据路径
:width: 100%
:align: center

第六遍的激活路径：Guard Comparator 产生条件码或故障类别，Control Unit 选择协议回退，Exception Unit 执行技术恢复或安全动作，Journal/checkpoint 保存可审计证据。
```

抽象图刻意把三件事分开：Guard Comparator 只分类证据，Control Unit 决定 retry、Reacquire 或 SafeStop，Exception Unit 才承担 SafeHold 等外部动作。陈旧或重复 Event 由 `CMD/EVT REG` 中的 identity 与 `LR` 中的最近结果共同识别；恢复后仍必须回到 Control Unit，不能绕过提交边界直接修改状态。

第六遍不再沿正常的：

```text
Acquire → Track → Verify → Lock
```

向前读，而是从一个问题出发：

> 当证据失效、测量失败、事件重复、进程中断或硬件状态不可信时，系统怎样避免继续使用错误历史？

这一遍主要阅读：

- `sqc/workflows/frequency_state_machine.py` 中各状态 handler；
- `_handle_technical_failure()`；
- `command_id`、`state_version` 和 pending command；
- `snapshot()` / `restore()`；
- `tests/unit/test_frequency_state_machine.py`；
- `tests/unit/test_frequency_runtime.py`；
- `tests/integration/test_frequency_state_machine_smoke.py`。

整体恢复结构是：

```text
科学证据仍可信但尚未达标
        → 留在 Track 或 Verify

局部参考失效，但系统仍可测量
        → Reacquire

暂时性技术失败
        → 原状态重试

硬件执行不一致或安全联锁
        → SafeStop

重复、陈旧或错误身份的事件
        → 忽略

进程中断
        → snapshot / save_run
        → 恢复 pending command
```

### 6.1 先按语义分类失败，再决定回退方向

#### 6.1 先区分三类“失败”

理解异常路径最重要的第一步，是不要把所有失败都叫作“测量失败”。

代码实际上处理三种不同问题。

##### 第一类：科学有效性失败

测量过程可能正常完成，但结果不能支持当前局部推断，例如：

```text
probe detuning 越出局部有效范围
灵敏度异常
置信度不足
结果有歧义
参考丢失
发生大幅频率跳变
```

此时硬件不一定坏了，测量程序也不一定抛异常；问题是：

> 当前局部模型和参考不能再安全使用。

这类问题通常进入：

```text
Reacquire
```

##### 第二类：技术执行失败

例如：

```text
求解器异常
硬件通信超时
仪器未返回数据
Backend 抛出异常
测量任务被技术性拒绝
```

它们被包装为：

```python
MeasurementTechnicalFailure
```

或某些：

```python
MeasurementRejected
```

这类失败通常先在当前状态重试。

##### 第三类：安全一致性失败

例如：

```text
实际偏置与命令偏置不一致
实际 drive 与命令 drive 不一致
Verify 时偏置没有真正冻结
外部 interlock 触发
```

这不是“再测一次可能就好”的普通失败，而是：

> 软件认为设备执行了 A，但设备报告实际执行了 B。

因此直接进入：

```text
SafeStop
```

三者的基本关系是：

```text
局部科学模型不再可信
    → Reacquire

测量暂时没执行成功
    → retry

设备执行状态与命令不一致
    → SafeStop
```

#### 6.2 Reacquire 不是技术重试

Reacquire 是正式协议状态，不是一个异常捕获器。

进入 Reacquire 表示：

> 当前局部反演参考已经不能继续使用，需要重新获得一个较宽范围的绝对频率种子。

在 Reacquire 中，状态机发出的仍是：

```python
AcquireFrequency
```

但带有：

```python
role="reacquire"
```

而且偏置优先使用最后一个候选偏置：

```python
bias = candidate_bias
```

如果没有候选偏置，才使用：

```python
config.acquisition_bias
```

因此 Reacquire 的物理含义是：

```text
保持对当前工作区附近的关注
        ↓
重新做较宽范围绝对测频
        ↓
重新建立频率种子
        ↓
决定去 Track 还是直接 Verify
```

它不是简单地回到初始 Acquire 状态，因为系统仍保留“最后工作点在哪里”的信息。

#### 6.3 哪些 Track 结果进入 Reacquire

Track 成功返回 `MeasurementSucceeded` 后，代码按固定顺序检查 guard。

如果发生以下问题：

```text
probe detuning 越界
缺少可靠 probe detuning
测量 uncertainty 非法
confidence 不足
sensitivity 超出允许范围
diagnostics 明确报告 out_of_range
```

通常转移为：

```text
Track → Reacquire
```

并记录具体 `ReasonCode`：

```text
LOCAL_RANGE_LOST
CONFIDENCE_INSUFFICIENT
SENSITIVITY_INVALID
```

如果：

```python
valid=False
```

或：

```python
ambiguous=True
```

也进入：

```text
Track → Reacquire
reason = CONFIDENCE_INSUFFICIENT
```

这里的判断不是：

> 频率离目标太远，所以 Reacquire。

而是：

> 当前局部 probe 无法可靠告诉我们频率在哪里，所以 Reacquire。

只要 drive tracking 让实际 probe detuning 仍在有效窗口内，即使目标残差很大，也可以继续 Track。

#### 6.4 为什么 Track 执行不一致不去 Reacquire

Track guard 还比较：

```text
event.applied_bias 与 command.candidate_bias
event.applied_drive 与 command.predicted_drive
```

允许误差由：

```text
actuation_bias_tolerance
actuation_drive_tolerance
```

定义。

如果不一致，返回：

```text
INTERLOCK
```

状态机执行：

```text
Track → SafeStop
run_status = FAILED
```

原因是 Reacquire 只能修复“我不知道频率在哪里”，不能修复：

> 我不知道硬件到底执行了哪个偏置和 drive。

如果实际执行参数都不可信，即使重新测到一个频率，也无法把频率与控制工作点正确对应。因此这里必须停止，而不是继续自动恢复。

#### 6.5 Verify 有三种失败方向

Verify 的恢复逻辑比 Track 更细。

##### 可靠且达到最终容差

如果：

```text
|r_ver| + z σ_ver ≤ epsilon_final
```

增加 Verify streak。

达到 `N_verify` 后：

```text
Verify → Lock
```

##### 可靠但需要局部修正

如果没有达到 `epsilon_final`，但仍满足：

```text
|r_ver| + z σ_ver ≤ verify_track_max_residual
```

则：

```text
Verify → Track
reason = CORRECTION_REQUIRED
```

这表示独立 Ramsey 可信，只是确认候选仍有可由局部控制器修正的小残差。

##### 不适合局部恢复

如果结果虽以成功事件返回，但：

```text
valid=False
ambiguous=True
```

或残差已经超过 Verify 的局部恢复范围，则：

```text
Verify → Reacquire
```

这表示不能安全地把 Verify 结果作为局部 Track 的新起点。

#### 6.6 Verify 的偏置冻结失败为什么直接 SafeStop

Verify 命令携带：

```python
frozen_bias
```

测量事件必须报告：

```python
applied_bias
```

状态机检查：

```text
|applied_bias - candidate_bias|
≤ bias_freeze_tolerance
```

如果不满足：

```text
Verify → SafeStop
reason = INTERLOCK
run_status = FAILED
```

因为 Verify 的科学意义建立在：

> 磁通偏置冻结后，用独立 Ramsey 重新确认同一个候选点。

如果测量时偏置变了，那么即使 Ramsey 结果很好，它验证的也不是原候选点。此时不能返回 Track，也不能 Reacquire 后假装这次验证有效，而应停止并检查执行链。

#### 6.7 一个需要按当前代码精确理解的细节

当前 Verify 对两种低质量结果的处理并不完全相同。

如果 Backend 返回：

```python
MeasurementSucceeded(
    valid=False
)
```

或：

```python
MeasurementSucceeded(
    ambiguous=True
)
```

代码立即进入 Reacquire。

但如果 Backend 返回：

```python
MeasurementRejected(
    reason_code=CONFIDENCE_INSUFFICIENT
)
```

当前 `_handle_verify()` 会把它送入通用技术失败处理，先在 Verify 重试，而不是立即 Reacquire。

所以当前实际语义是：

```text
成功执行但证据无效/歧义
    → Reacquire

Backend 拒绝产生 Verify 结果
    → 作为技术失败重试
```

这可以解释为：

- 前者已经证明科学证据不可用；
- 后者可能只是一次暂时无法完成的测量。

但硬件 Backend 在设计 Event 时必须保持一致，否则同一种“低置信度”可能因事件类型不同走向不同路径。

#### 6.8 Lock 的恢复逻辑为什么与 Track 不同

Lock 中的低成本 monitor 没有偏置更新权。

如果 monitor 清晰：

```text
Lock → Lock
reason = MONITOR_CLEAR
```

如果落入疑似漂移区，累计：

```python
monitor_streak
```

达到 `N_mon_suspect` 后：

```text
Lock → Verify
reason = DRIFT_SUSPECTED
```

如果漂移介于 suspect 上界和大跳变阈值之间，也直接进入 Verify。

如果达到大跳变条件：

```text
Lock → Reacquire
reason = LARGE_FREQUENCY_JUMP
```

如果：

```python
valid=False
```

或：

```python
ambiguous=True
```

则：

```text
Lock → Reacquire
reason = REFERENCE_LOST
```

因此 Lock 恢复策略是：

```text
小且清晰
    → 继续 Lock

可能漂移
    → Verify

大跳变或参考丢失
    → Reacquire
```

不存在：

```text
Lock monitor
    → 直接调整偏置
```

#### 6.9 Lock 技术失败不会立即解释成参考丢失

如果 monitor 以成功事件返回，但 `valid=False`，代码认为参考丢失，进入 Reacquire。

如果 Backend 根本没有成功完成测量，而是返回：

```python
MeasurementTechnicalFailure
```

当前代码先在 Lock 中技术重试。

区别是：

```text
完成了测量，结果证明参考无效
    → Reacquire

测量本身没有正常完成
    → retry
```

只有技术重试耗尽后才进入 SafeStop。

这样避免把一次临时仪器超时误解释成量子比特发生了频率跳变。

### 6.2 Reacquire 与技术重试使用不同的恢复预算

#### 6.10 Reacquire 成功后清除什么

Reacquire 获得可靠种子后会更新：

```text
f_hat
uncertainty
candidate_bias
candidate_bias_version
candidate_source_event = "reacquire"
```

并清除：

```text
tracker_snapshot
verify_streak
monitor_streak
hold_target_met
hold_samples
hold_passes
technical_retries
```

这表示：

> 旧局部控制历史、旧验证连续命中和旧 Lock 保持统计都不再属于新参考。

尤其是 `_tracker_snapshot=None` 很重要。旧割线灵敏度是在已经失效的局部参考附近建立的，Reacquire 后不能继续沿用。

但系统不会清空整个运行历史：

```text
transition log 仍保留
Budget 仍累计
journal 仍保留
run_id 不变
```

所以 Reacquire 是同一次运行中的恢复，不是重新创建一场全新运行。

#### 6.11 Reacquire 成功后去 Track 还是 Verify

新种子同样使用候选进入条件：

```text
|r_acq| + z σ_acq ≤ epsilon_enter
```

如果满足：

```text
Reacquire → Verify
```

如果不满足，但种子可靠：

```text
Reacquire → Track
```

因此 Reacquire 不是固定返回 Track。

完整恢复路径是：

```text
Reacquire
   ↓ 宽范围测频
可靠且已接近目标
   → Verify

可靠但仍需修正
   → Track
```

#### 6.12 Reacquire attempt 在什么时候记账

`_handle_reacquire()` 一进入就执行：

```python
budget.record_reacquire()
```

也就是说，每处理一次 Reacquire 命令返回的 Event，就累计一次 attempt，包括：

- 成功测量；
- 技术失败；
- 被拒绝；
- 无效或歧义结果。

它统计的是：

> 实际完成了多少次 Reacquire 尝试。

不是：

> 进入 Reacquire 状态多少次。

一次 Reacquire 状态内的多个技术重试会分别消耗 attempt 预算。

#### 6.13 为什么“最后一次允许的 Reacquire 成功”仍能继续

假设：

```python
max_reacquire_attempts = 1
```

第一次 Reacquire 成功后：

```text
reacquire_attempts = 1
state = Track 或 Verify
```

下一次 `next_command()` 不会因为计数已经等于上限而立即停止，因为 Reacquire 上限检查只在当前状态仍为 Reacquire 时阻止下一次宽范围捕获。

所以：

```text
最后一次允许的尝试成功
    → 可以继续 Track/Verify

最后一次允许的尝试失败
    → Reacquire limit
    → SafeStop
```

测试专门验证了这个边界，避免把“最多一次 Reacquire”错误实现为“第一次 Reacquire 即使成功也必须停止”。

#### 6.14 技术重试与 Reacquire attempt 是两套计数

状态机还有：

```python
_technical_retries
```

技术重试计数回答：

> 当前连续发生了多少次技术执行失败？

Reacquire attempt 回答：

> 宽范围恢复测量一共尝试了多少次？

两者用途不同。

例如一次 Reacquire 技术失败会同时影响：

```text
technical_retries += 1
reacquire_attempts += 1
```

前者限制连续技术异常，后者限制宽范围恢复资源。

任何成功的正常测量都会把：

```python
_technical_retries = 0
```

所以技术失败计数是连续失败 streak，不是整次运行的技术错误总数。

#### 6.15 max_technical_retries=N 的准确含义

`_handle_technical_failure()` 先执行：

```python
technical_retries += 1
```

然后判断：

```python
if technical_retries > max_retries:
    SafeStop
else:
    stay and retry
```

因此：

```python
max_technical_retries = 2
```

表示：

```text
第 1 次失败：允许重试
第 2 次失败：允许重试
第 3 次失败：重试机会耗尽，SafeStop
```

也就是允许两次失败后的重试，而不是总共只允许两次测量尝试。

Acquire 使用独立配置：

```python
max_acq_retries
```

其他状态通常使用：

```python
max_technical_retries
```

测试验证 `max_acq_retries=1` 时：

```text
第一次 Acquire 失败 → 留在 Acquire
第二次失败         → SafeStop
```

#### 6.16 技术重试会不会复用同一 command ID

不会。

一条技术失败事件被当前 handler 接受后，handler 最后执行：

```python
_clear_pending()
```

下一轮 `next_command()` 创建新的 Command 和新的 ID。

所以技术重试是：

```text
第一次 TrackFrequency，command_id=A
    ↓ TechnicalFailure
Track → Track
清除 A

第二次 TrackFrequency，command_id=B
```

而不是反复执行同一个 pending command。

同一 ID 的重复返回只发生在：

```text
pending command 尚未被匹配 Event 消费
```

例如：

- Runtime 重复调用 `next_command()`；
- checkpoint 恢复未完成命令；
- 调度层重新读取当前 pending command。

#### 6.17 MeasurementTechnicalFailure 也要记录部分成本

技术失败不等于零成本。

Backend 可能已经执行了部分电路或 solver 调用后才失败，因此事件包含：

```text
shots
solver_calls
elapsed_time
diagnostics
```

状态机在分派失败路径之前先执行：

```python
budget.consume_shots(event.shots)
budget.consume_solver_calls(event.solver_calls)
```

测试构造：

```text
shots = 3
solver_calls = 11
```

即使事件最终进入技术重试，Budget 仍累计这部分成本。

否则一个不断执行一半后失败的硬件任务可能绕过预算限制，无限消耗资源。

#### 6.18 MeasurementRejected 不一定等于 TechnicalFailure

`MeasurementRejected` 表示：

> 某一层决定不接受或不执行这次测量。

它可能来自：

```text
测前 prediction guard
Backend 置信度门限
参考失效
资源或硬件策略拒绝
```

状态机要结合：

```text
当前状态
reason_code
```

解释它。

在 Track 中，以下拒绝原因直接进入 Reacquire：

```text
LOCAL_RANGE_LOST
CONFIDENCE_INSUFFICIENT
SENSITIVITY_INVALID
```

其他拒绝通常走技术重试。

在 Verify 和 Lock 中，当前通用 `MeasurementRejected` 主要走技术重试，除非 Backend 把问题包装为成功但无效/歧义的科学结果。

因此 Event 类型和 `ReasonCode` 共同定义语义，不能只看类名。

### 6.3 SafeStop：终止原因与安全动作必须分别记录

#### 6.19 SafeStop 有三种不同 RunStatus

进入 `FrequencyState.SAFE_STOP` 时，`run_status` 不一定相同。

##### SAFE_STOPPED

用于外部主动安全停止：

```text
CancelRequested
InterlockTriggered
```

例如用户取消后：

```text
state = SafeStop
run_status = SAFE_STOPPED
```

##### FAILED

用于尚未成功标定前的失败，例如：

```text
Track interlock
Verify freeze interlock
技术重试耗尽
Reacquire limit
预算在标定完成前耗尽
```

##### COMPLETED

用于已经进入 Lock 或已经被标记为 `CALIBRATED` 后，由有界运行条件正常结束：

```text
stop_after_lock_cycles 达到上限
长期运行预算到期
```

此时：

```text
state = SafeStop
run_status = COMPLETED
```

所以不能仅看到 `SafeStop` 就断言运行失败。

应同时读取：

```text
state
run_status
transition reason
safe_hold_confirmed
completed_verified_lock
```

#### 6.20 SafeStop 仍然需要 SafeHold

状态机进入 SafeStop 后，下一条命令是：

```python
SafeHold(
    bias=0.0,
    reason=...
)
```

这一点在第五遍已经讲过，但在恢复路径中尤其重要：

```text
协议不再继续
```

不等于：

```text
设备已经进入安全偏置
```

Runtime 必须看到：

```python
SafeHoldApplied
```

才设置：

```python
safe_hold_confirmed=True
```

SafeHold 失败可以重试。最终若：

```text
state = SafeStop
safe_hold_confirmed = False
```

说明协议已经停止，但硬件安全动作没有获得确认，必须由上层系统报警或人工处置。

### 6.4 身份、幂等与恢复共同抵御迟到或重复事件

#### 6.21 自发事件为什么不要求 command ID

以下事件不是某条测量命令的响应：

```text
InterlockTriggered
CancelRequested
BudgetExhausted
TimerElapsed
```

它们可以在命令之间甚至存在 pending command 时到达。

因此 `handle()` 先处理这些自发事件，再检查 command ID。

例如 Cancel 到达时：

```text
当前有 pending VerifyFrequency
        ↓
CancelRequested
        ↓
清除 pending Verify
进入 SafeStop
```

如果先要求 Cancel 携带 Verify 的 command ID，就无法表达真正的外部中断。

#### 6.22 command_id 怎样阻止重复事件

对命令响应事件，状态机比较：

```python
event.command_id
```

与：

```python
_pending_command_id
```

如果不一致：

```python
return
```

不转换状态，不清除 pending command。

例如：

```text
pending command ID = B
旧事件 ID = A
```

旧事件被视为：

```text
stale or mismatched
```

并忽略。

测试先处理一次 command A，使状态进入 Verify，再重复提交 A 的成功事件。第二次事件不能增加 Verify streak，也不能进入 Lock。

#### 6.23 state_version 为什么还要再检查一次

command ID 识别具体命令，state version 识别协议世代。

状态机要求：

```python
event.state_version in (0, current_state_version)
```

其中 0 是兼容同步 Backend 的默认值；Runtime 通常在执行后把它补成命令执行前的版本。

如果事件携带其他版本：

```python
return
```

例如：

```text
当前 state_version = 7
事件声称来自 version = 8
```

即使命令 ID 看似正确，也不接受。

两道检查分别解决：

```text
command_id：
    这是哪一条命令的结果？

state_version：
    它属于哪一代协议状态？
```

#### 6.24 陈旧事件为什么连成本都不计

身份检查发生在 Budget 记账之前：

```text
先检查 command_id
再检查 state_version
然后才 consume shots / solver calls
```

因此错误身份事件不会：

- 改变状态；
- 清除 pending command；
- 增加 shots；
- 增加 solver calls；
- 改写 last event。

这是正确的状态机语义，因为它没有证据证明该事件属于本次运行的当前命令。

但 Runtime 仍可能把 Executor 返回的事件写入执行 journal。正常内置 Backend 总是回传匹配 ID；自定义或远程 Backend 若可能错配，调用方应把 journal 中“收到过该事件”和状态机中“接受了该事件”区分开。

#### 6.25 next_command() 的幂等性解决什么

只要 pending command 尚未收到匹配事件：

```python
cmd1 = machine.next_command()
cmd2 = machine.next_command()
```

两次返回同一对象语义和同一：

```text
command_id
```

而且不会重复增加：

```text
commands_issued
```

这允许 Runtime 在以下情况下安全重取命令：

```text
预算检查后需要重新读取
checkpoint 恢复
调度器重试获取
网络层重复拉取
```

一旦匹配 Event 被处理，handler 清除 pending command，下一次才生成新 ID。

#### 6.26 状态自循环为什么也增加 state_version

技术重试或正常 Track 自循环调用：

```python
_stay(state, reason)
```

虽然：

```text
from = Track
to = Track
```

枚举状态不变，但：

```python
state_version += 1
```

因为已经完成了一轮新的协议判断。

例如：

```text
Track version 4
    ↓ 测量仍需修正
Track version 5
```

这两个 Track 不是同一协议时刻。旧 version 4 的迟到事件不能再污染 version 5。

#### 6.27 snapshot 恢复的不只是 State 枚举

`StateMachineSnapshot` 保存：

```text
state
run_status
state_version

candidate_bias
candidate_bias_version
candidate_source_event

f_hat
uncertainty

tracker_snapshot

verify_streak
monitor_streak
Lock 时间和 hold 统计

Budget

transition_log
last_event

pending_command
pending_command_id
technical_retries
```

因此恢复能够继续：

- 同一个候选点；
- 同一个 Verify streak；
- 同一段割线历史；
- 同一份资源预算；
- 同一个技术失败 streak；
- 同一条尚未完成的命令。

如果只保存：

```json
{"state": "verify"}
```

则无法知道 Verify 在验证哪个偏置，也不知道已经连续通过几次，因此不具备真实恢复能力。

#### 6.28 Command 怎样被序列化

pending command 保存为：

```python
{
    "type": "VerifyFrequency",
    "payload": {
        "command_id": "...",
        "frozen_bias": ...,
        "drive": ...
    }
}
```

恢复时 `_deserialise_command()` 根据 `type` 重新构造具体 dataclass。

支持：

```text
AcquireFrequency
TrackFrequency
VerifyFrequency
MonitorFrequency
SafeHold
```

SafeHold 的 `ReasonCode` 会从字符串重新转换成 Enum。

未知命令类型会抛出错误，而不是静默降级。这比“忽略无法识别的命令并重新 Acquire”安全，因为后者可能掩盖 checkpoint 版本不兼容。

#### 6.29 恢复 pending command 为什么是 at-least-once

如果 checkpoint 保存了 pending Verify：

```text
VerifyFrequency ID=A
```

恢复后：

```python
next_command()
```

会返回同一个 A，而不是重新创建 Acquire。

这保证：

```text
不会遗漏未确认完成的命令
```

但不能保证 exactly-once。

崩溃可能发生在：

```text
硬件已经执行 A
Event 尚未持久化
进程退出
```

恢复后 A 会再次执行。

因此当前恢复语义是：

```text
at least once
```

要安全使用这一语义，硬件命令最好是绝对设定：

```text
把偏置设为 0.03 Φ0
```

而不是非幂等增量：

```text
在当前偏置上再增加 0.03 Φ0
```

#### 6.30 哪些恢复信息当前没有完全持久化

`save_run()` 已保存状态机、journal、cost ledger、Config、Backend provenance 和可选 RNG state。

但仍有一些执行环境信息需要注意。

##### monitor 的 Runtime 本地计时

Runtime 的：

```python
_last_monitor_time
```

不是 `StateMachineSnapshot` 字段。

状态机保存了：

```python
_last_audit_time
_lock_entry_time
```

但恢复后 Runtime 的低成本 monitor 间隔计时默认重新从 0 开始，因此第一条恢复后的 monitor 可能立即执行。

##### 内存 checkpoint 列表

`_checkpoints` 的历史列表不会完整写入磁盘；`save_run()` 写的是当前最新状态快照。

##### 外部硬件真实状态

checkpoint 保存命令和软件认知，但不会自动查询 DAC、AWG、LO 或量子比特是否仍处在保存时的物理状态。

真实硬件恢复前通常还需要：

```text
读取当前硬件状态
核对偏置和 drive
检查参考源
必要时重新 Reacquire
```

### 6.5 测试证据能证明什么，也不能证明什么

#### 6.31 Reacquire 是否真的恢复了实验参考

在确定性仿真中，Reacquire 可以通过新的 Ramsey 结果恢复频率种子。

但在真实实验中，“参考恢复”通常还需要证明：

```text
读出链可用
时钟/LO参考正常
脉冲标定仍有效
频率峰没有混叠
磁通执行可信
```

当前状态机通过：

```text
valid
ambiguous
uncertainty
applied_bias
diagnostics
```

接收这些证据，但具体怎样从硬件原始数据生成它们，属于实验 Backend。

因此现有代码证明的是：

> 只要 Backend 能可靠报告这些字段，状态机能够按协议恢复。

它不等于已经在真实仪器上验证了所有参考丢失情形。

#### 6.32 故障注入 Backend 测试了什么

`FaultInjectionExecutor` 可以按状态注入：

```text
MeasurementTechnicalFailure
MeasurementRejected
ambiguous=True
```

并配置：

```python
fail_count
```

测试验证两次故障注入后，第三次调用会恢复到真实 Executor。

它证明：

- 状态机和 Runtime 可以接收可控故障；
- 重试路径能够被确定性复现；
- 不需要依赖随机偶发错误测试恢复逻辑。

但它不证明：

- 真机故障的时间分布；
- 网络故障持续时间；
- 仪器取消能否及时响应；
- 硬件恢复后数据一定无偏。

故障注入是协议测试工具，不是真机可靠性证据。

#### 6.33 ProcessSimulationExecutor 测试了什么

过程仿真包装器可以加入：

```text
drift_per_command
jump_schedule
reference_loss_commands
```

例如：

```python
jump_schedule={0: 2 * np.pi * 1e-3}
reference_loss_commands={1}
```

第一次测量频率加入受控跳变；第二次返回参考丢失拒绝。

它使以下状态路径可测试：

```text
Lock suspected drift → Verify
large jump → Reacquire
reference loss → Reacquire
```

但当前过程模型是按命令序号叠加的确定性偏移或预设事件，不等于真实的：

- 连续时间随机漂移；
- 频率噪声 PSD；
- 非平稳跳变过程；
- Allan deviation；
- 长期闭环稳定性实验。

#### 6.34 单元测试主要证明什么

`test_frequency_state_machine.py` 主要把 Backend 替换为人工 Event，因此能够精确覆盖：

```text
每条状态转换
每个 ReasonCode
阈值边界
重试次数
Reacquire limit
预算
幂等性
陈旧事件
snapshot
命令类型
transition log
```

它证明的是：

> 给定某种 Event，协议决策符合定义。

它不执行真实 Ramsey 或 transient 物理模拟。

#### 6.35 Runtime 测试主要证明什么

`test_frequency_runtime.py` 使用轻量 RecordingExecutor 等测试替身，重点覆盖：

```text
取消信号 first-request-wins
测量前取消
测量中取消
KeyboardInterrupt
SafeHold 确认与重试
monitor interval 可中断等待
预算阻止额外测量
journal 可序列化
cost ledger
save/load
RNG 状态恢复
pending command 重放
```

它证明的是：

> 调度、中断和恢复时序正确。

它通常不证明频率估计器的物理准确性。

#### 6.36 Integration smoke 测试主要证明什么

`test_frequency_state_machine_smoke.py` 使用真实：

```text
TransmonQubit
SQCExecutor
FrequencyMeasurement
Ramsey/transient 模拟
```

覆盖：

```text
Acquire 能产生测量
Acquire 可直接进入 Verify
Verify 双次通过可进入 Lock
Track 确实在建议偏置测量
有限-shot可重复
solver calls 被记录
故障注入可以穿过 Backend
```

它证明跨层接口真正能接起来。

但 smoke 测试的参数规模和运行时间有限，不能代替：

- 全范围数值扫描；
- 长时间漂移仿真；
- 真机压力测试；
- 硬件恢复演练。

### 6.6 第六遍收束：从故障时间线反向定位责任层

#### 6.37 怎样从测试反向读源码

读异常体系时，推荐按测试路径而不是按文件行号阅读。

##### 路径一：局部范围丢失

```text
test_probe_detuning_outside_validated_range_goes_to_reacquire
    ↓
_track_guard_failure()
    ↓
_handle_track()
    ↓
Track → Reacquire
```

##### 路径二：Verify 偏置不冻结

```text
test_verify_bias_mismatch_safe_stops
    ↓
_handle_verify()
    ↓
bias_freeze_tolerance
    ↓
Verify → SafeStop
```

##### 路径三：技术重试

```text
test_technical_failure_retries_then_safe_stop
    ↓
_handle_technical_failure()
    ↓
technical_retries > max_retries
```

##### 路径四：陈旧事件

```text
test_stale_event_ignored
    ↓
command_id mismatch
    ↓
handle() early return
```

##### 路径五：恢复 pending Verify

```text
test_resume_reuses_pending_command_and_does_not_restart_acquire
    ↓
save_run()
    ↓
pending command serialization
    ↓
load_run()
    ↓
next_command() idempotent return
```

这种读法能把“测试名称”变成源码导航。

#### 6.38 完整的 Track 失效恢复时间线

假设 Track 实际 probe detuning 越界：

```text
1. Runtime 发出 TrackFrequency A
2. Backend 在 candidate_bias 和 predicted_drive 上测量
3. 返回 MeasurementSucceeded
4. probe_detuning 超过 Delta_val
5. machine.handle(event A)
6. _track_guard_failure() 返回 LOCAL_RANGE_LOST
7. Track → Reacquire
8. 清除 pending A
9. Runtime 不更新 tracker snapshot
10. 下一轮发出 AcquireFrequency B
    role = reacquire
    bias = last candidate_bias
11. 宽范围 Ramsey 获得新种子
12. 清除旧 tracker snapshot
13. 若已接近目标 → Verify
14. 否则 → Track
```

关键点是第 9 步：

> 导致 Reacquire 的 Track 测量不能再被控制器接受为新割线点。

Runtime 的 `_post_process()` 会检查状态已经进入 Reacquire，因此拒绝更新 tracker。

#### 6.39 完整的技术失败恢复时间线

假设 Track Backend 暂时超时：

```text
1. TrackFrequency A
2. Backend 返回 MeasurementTechnicalFailure A
3. 记录部分 shots / solver calls
4. technical_retries = 1
5. Track → Track
   reason = TECHNICAL_TIMEOUT
6. 清除 pending A
7. 下一轮生成 TrackFrequency B
8. Backend 成功
9. technical_retries = 0
10. 正常继续 Track/Verify
```

若连续失败超过上限：

```text
Track → SafeStop
run_status = FAILED
→ SafeHold
→ SafeHoldApplied
```

#### 6.40 完整的进程恢复时间线

假设运行在 Verify 中保存：

```text
state = Verify
candidate_bias = 0.03
verify_streak = 1
pending VerifyFrequency ID=A
```

恢复过程：

```text
1. load_run() 读取 config.json
2. 重建 FrequencyCalibrationConfig
3. 恢复 run_id
4. 创建 Runtime 和 Executor
5. 恢复 Backend RNG state
6. 读取 journal 与 cost ledger
7. 读取 checkpoint.json
8. 重建 Budget
9. 重建 StateMachineSnapshot
10. 创建新的 FrequencyStateMachine
11. restore(snapshot)
12. pending command A 被反序列化
13. 调用 restored.run()
14. next_command() 返回 A
15. 不执行 Acquire
16. Executor 重新执行 Verify A
17. 状态机继续 verify_streak
```

这里恢复的是协议连续性，不是自动证明硬件物理状态仍连续。

#### 6.41 第六遍应形成的整体认识

到这一遍，状态机不应再被理解为一个“成功时往右走”的流程图，而应被理解为一个证据分级系统：

```text
证据可信、目标未到
    → 继续局部控制

证据可信、需要独立确认
    → Verify

局部参考失效，但设备仍可测
    → Reacquire

测量暂时未完成
    → 技术重试

执行参数不一致或安全联锁
    → SafeStop

事件身份不属于当前协议世代
    → 忽略

进程中断但状态可恢复
    → checkpoint + pending command replay
```

第六遍最重要的结论是：

1. Reacquire 修复局部科学参考，不负责修复硬件执行不一致。
2. 技术失败先在原状态重试；连续失败超过配置上限才 SafeStop。
3. `max_technical_retries=N` 表示允许 N 次失败后重试，第 N+1 次才耗尽。
4. Reacquire attempt 与 technical retry 是两套不同预算。
5. 最后一次允许的 Reacquire 如果成功，可以继续 Track 或 Verify。
6. Verify 偏置未冻结和 Track 实际执行参数不一致属于 interlock，不能自动恢复。
7. command ID 与 state version 共同防止重复、陈旧或跨状态事件污染当前协议。
8. snapshot 必须包含 pending command、Budget、tracker 和 streak，只有 State 枚举远远不够。
9. 当前恢复属于 at-least-once；真实硬件需要幂等绝对命令或设备侧任务查询。
10. 单元测试证明协议逻辑，Runtime 测试证明调度恢复，integration smoke 才证明真实模拟层能够接通。
11. 故障注入和过程仿真是确定性测试证据，不是真机长期可靠性或稳定性证据。
12. SafeStop 是协议终止状态，`SafeHoldApplied` 才是安全动作确认。


