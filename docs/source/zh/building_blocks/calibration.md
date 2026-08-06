# 标定(calibration)

## 这层提供什么

`calibration` 层是全栈的标定与整定层,测量比特与控制线当前的状态,
并在需要时把它调到目标状态。它服务两条产品主线:频率标定(测/整 $f_{01}(\Phi)$)
与预畸变(测控制线传函、设计补偿滤波器)。产出统一封装成
{py:class}`~sqc.calibration.CalibrationTable`,供下游查表反演或整定使用。

与 {doc}`reconstruction` 层的纯函数约束不同,本层的标定类会主动运行仿真去测量:
每个 `Calibration` 子类实现 `calibrate()`,内部驱动比特、跑 `mesolve`、拟合,最终返回一张
`CalibrationTable`。

```{note}
{doc}`reconstruction` 层还导出两个概念属标定的类 `CryoscopeCalibration` /
`DelayRamseyCalibration`(与其重建器配套导出)。它们同样继承本层的
{py:class}`~sqc.calibration.Calibration`,产出 `CalibrationTable` 供 `inversion="calibration"`
查表反演。分述见 {doc}`reconstruction`。
```

## 类总览

按功能分五组:

**扩展点 + 结果容器**

| 类 | 角色 |
|---|---|
| `Calibration` | 抽象基类,所有标定工作流的公共契约,本层扩展点 |
| `CalibrationTable` | 标定结果容器,带 `evaluate`(插值)/ `inverse`(反插值) |

**频率标定 / 测量 / 整定**

| 类 | 角色 |
|---|---|
| `FluxResponseCalibration` | 扫磁通建 $f(\Phi)$ 查表(Ramsey 逐点测频) |
| `FrequencyMeasurement` | 单点 $f_{01}$ 测量(Ramsey,或进阶的瞬态法) |
| `SinglePointFrequencyCalibration` | 单点频率**整定**:闭环反馈把 $f_q(V)$ 调到目标 |
| `FrequencyCalibrationWorkflow` | 频率标定四阶段状态机的编排实现 |
| `CalibrationStage` | 状态机单阶段的参数容器(`FrequencyCalibrationWorkflow` 使用) |

**波形 / 控制线标定**

| 类 | 角色 |
|---|---|
| `WaveformCalibration` | 传函测量 + 预畸变设计的统一入口 |
| `PredistortionDesigner` | 独立的逆滤波器设计器(FIR/IIR/频域反演) |

**调度**

| 类 | 角色 |
|---|---|
| `CalibrationScheduler` | 标定任务调度器:注册表 + 依赖图,按序执行 |

## Calibration：标定抽象基类

所有标定工作流的公共契约,本层扩展点。只规定一个抽象方法:

- `calibrate() -> CalibrationTable`:运行标定工作流,返回结果表。

要添加自定义标定,继承 `Calibration` 实现 `calibrate()`,详见 {doc}`../extending`。

## CalibrationTable：标定结果容器

统一的标定结果载体。

**构造**

`CalibrationTable(name, qubit_name=None, kind=None, inputs=None, outputs=None, fit_params=None, metadata=None)`

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `name` | str | 标定项名称 |
| `qubit_name` | str | 所属比特名 |
| `kind` | str | 标定类型(`"f_phi"`、`"f01"`、`"transfer_function"`、`"predistortion"` 等) |
| `inputs` | `np.ndarray` | 自变量(如磁通 $\Phi$) |
| `outputs` | `np.ndarray` | 因变量(如频率 $\omega$) |
| `fit_params` | dict | 拟合/迭代细节 |
| `metadata` | dict | 附加元数据 |

**方法**

- `evaluate(x) -> np.ndarray`：在查询点 `x` 上插值 `outputs`(如给磁通取频率)。
  基于 SciPy 三次样条,点数不足时自动降到二次/线性,并外插。
- `inverse(y) -> np.ndarray`：反插值,找使 `outputs≈y` 的 `inputs`(需 `outputs` 单调,
  内部取单调段构造反函数)。

## FluxResponseCalibration：f(Φ) 曲线标定

扫一串直流磁通点,在每点用 Ramsey 测频,建 $f(\Phi)$ 查表。

**构造**

`FluxResponseCalibration(qubit, method="ramsey", h_list=None, tau=None, t_rabi=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特 | —— |
| `method` | str | 测量方法 | `"ramsey"` |
| `h_list` | `np.ndarray` | 磁通扫描点($\Phi_0$) | `linspace(-0.03, 0.03, 51)` |
| `tau` | float | 每点自由进动时间(ns) | —— |
| `t_rabi` | `np.ndarray` | Rabi 脉冲时间轴(ns) | `CONFIG.pulse.t_rabi` |

**方法**

- `calibrate() -> CalibrationTable`：逐磁通点跑 Ramsey 测频,返回 `kind="f_phi"` 的表,
  `inputs=磁通`、`outputs=角频率`。

**输出**

`CalibrationTable`,`kind="f_phi"`,`inputs` 为磁通扫描点、`outputs` 为对应角频率。

```{note}
`method="transient"`(未知瞬态信号 → $\Delta\omega(\Phi)$ 多项式拟合)属未来功能,
当前抛 `NotImplementedError`,不进 v1 公开面。
```

## FrequencyMeasurement：单点 f01 测量

在单个磁通工作点测 $f_{01}$(只读、不整定)。

**构造**

`FrequencyMeasurement(qubit, flux=0.0, method="ramsey", tau_list=None, t_rabi=None, t_global=None, f_artificial=0.1)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特 | —— |
| `flux` | float | 测量磁通点($\Phi_0$) | `0.0`(甜点) |
| `method` | str | 测量方法 | `"ramsey"`(或 `"transient"`) |
| `tau_list` | `np.ndarray` | 自由演化时间扫描(ns) | `CONFIG.pulse.tau_list` |
| `t_rabi` | `np.ndarray` | Rabi 时间轴(ns) | `CONFIG.pulse.t_rabi` |
| `t_global` | `np.ndarray` | 全局时间轴(ns) | `CONFIG.pulse.t_global` |
| `f_artificial` | float 或 None | 人工失谐(GHz),`None`=双扫模式 | `0.1` |

**方法**

- `measure(flux=None) -> float`：返回有符号角频率(rad·GHz)。
  - `method="ramsey"`(默认):Ramsey $\tau$ 扫描 + FFT 取峰。单扫模式(`f_artificial=0.1`,
    假设 $|\Delta|<0.1$ GHz)较快;`f_artificial=None` 走双扫模式,对任意 $|\Delta|$ 稳健、
    给符号,代价 2×。
  - `method="transient"`(进阶):$\tau=0$ 正交 Ramsey(R_y–R_x 与 R_y–R_{-x})差分读出
    + 控制核灵敏度 $G=\int k_1\,\mathrm{d}t$ 直接反出 $\Delta\omega$。比 $\tau$ 扫描省,但依赖弱信号
    线性近似,适合 $|\Delta\omega|$ 近零。`order>=3` 时加三次 Newton 修正,三次系数 $G_3$
    由 `g3_source` 选:`"fit"`(奇多项式拟合 $p_\mathrm{diff}(\Delta)$,自适应扫描范围)或
    `"kernel_full"`(完整非对角核三重积分 $\iiint k_3$)。
- `calibrate() -> CalibrationTable`：把单点测量打包成 `kind="f01"` 的表。

**输出**

`measure()` 返回 `float`(有符号角频率);`calibrate()` 返回 `CalibrationTable`(`kind="f01"`)。

```{note}
瞬态法是单点**频率测量**,与 §8 边界排除的瞬态**频率标定**(`FluxResponseCalibration`
的 `method="transient"`,建 $f(\Phi)$ 曲线)不是一回事:前者已实现且公开,后者未实现。
瞬态**波形重建**则是 v1 核心特性,见 {doc}`reconstruction`。
```

## SinglePointFrequencyCalibration：单点频率整定

把 $f_q(V)$ 闭环反馈整定到目标频率 $f_\mathrm{target}$(Vepsalainen 2022)。

**构造**

`SinglePointFrequencyCalibration(qubit, method="closed_loop", f_target=None, V_a=None, V_b=None, epsilon_f=1e-4, max_iter=20, measure_method="ramsey", step_method="secant", bracket_tightening=True, drive_policy="sweet", ...)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特 | —— |
| `method` | str | 整定方法 | `"closed_loop"` |
| `f_target` | float | 目标频率(rad·GHz) | —— |
| `V_a` | float | 磁通电压括号下界 | —— |
| `V_b` | float | 磁通电压括号上界 | —— |
| `epsilon_f` | float | 收敛容差(GHz·2π) | `1e-4` |
| `max_iter` | int | 最大迭代次数 | `20` |
| `measure_method` | str | 测频方法 | `"ramsey"`(或 `"transient"`) |
| `step_method` | str | 根搜索步进器 | `"secant"`(或 `"bisection"`/`"gradient"`) |
| `bracket_tightening` | bool | regula falsi 缩窄区间 | `True` |
| `drive_policy` | str 或 callable | 驱动频率策略 | `"sweet"`(或 `"target"`/`"track"`/callable) |
| `sensitivity_source` | str | 灵敏度来源 | `"secant"`(或 `"model"`) |
| `linear_range` | float | 线性窗 $\Delta_\mathrm{lin}$ | —— |
| `rho` | float | 越界阈值比例 | `0.6` |
| `converge_streak` | int | 连续收敛轮数要求 | `1` |
| `omega_d_seed` | float | 首轮 track 驱动频率(rad·GHz) | `None` |
| `damping` | float | 梯度阻尼因子(仅 `step_method="gradient"`) | `0.8` |
| `V_seed` | float | 梯度初始磁通(仅 `step_method="gradient"`) | —— |

**方法**

- `calibrate() -> CalibrationTable`：执行闭环整定,返回 `kind="f01"` 的表,
  `fit_params["history"]` 含完整迭代轨迹。

**输出**

`CalibrationTable`,`kind="f01"`,`fit_params["history"]` 含每轮迭代的
$V_k$、$f_{q,k}$、$e_k$、$f_{d,k}$、$\delta_k$ 和 `out_of_range` 标志。

### 每轮的两种更新:磁通电压 与 驱动频率

闭环每一轮同时跑**两种相互独立的更新**,作用在两个不同的执行器上:

- **磁通电压更新**(根搜索,真正把 $f_q$ 推向目标):梯度步用归一化 Newton 步
  $V_{k+1}=V_k-\alpha\,e_k/\hat s_k$($\alpha=$ `damping`,
  $\hat s_k=\partial f_q/\partial V$),secant / bisection 步则用无模型割线 / 二分。
  局部线性模型下误差按 $e_{k+1}\approx(1-\alpha)e_k$ 单调收缩,收到**唯一**不动点
  $f_q=f_\mathrm{target}$。这是**唯一**能改变收敛目标的更新。
- **驱动频率更新**(观测器,让测量保持可信):`drive_policy` 每轮设定驱动频率 $f_d$,
  目的只是让**当轮的频率测量**落在线性区。它**不进入误差定义**,因此改变不了收敛
  目标,只保证 $e_k$ 可信。

**保证安全的不变式:** `measure()` 返回的是**绝对**频率 $f_d+\hat\delta$,闭环按
$e_k=\mathrm{measure}(V_k)-f_\mathrm{target}$ 算误差。所以无论 $f_d$ 在哪,控制误差
永远相对固定的绝对目标——驱动跟着 qubit 走**不会**把残余误差藏掉。

`drive_policy`(仅梯度步;其他步进器会报错)接受字符串或
`callable(state) -> omega_d`:

| 策略 | $f_d$ | 阶段 | 用途 |
|---|---|---|---|
| `"sweet"`(默认) | `qubit.frequency` | **COARSE_ACQUIRE** | 向后兼容;目标靠近甜点、或宽范围粗获取时够用 |
| `"target"` | `f_target` | **LOCKED** | 测得失谐**即**控制误差,终局最简形式 |
| `"track"` | $\hat f_{q,k}+\hat s_k\,(V_{k+1}-V_k)$ | **TRACKING** | 预测下一点的 $f_q$,搜索移动时让 $\|\delta\|$ 留在线性窗内 |
| `callable` | 用户自定义 | 任意 | 任意反馈/滤波/预测律 |

### 四阶段状态机与 FrequencyCalibrationWorkflow

把上面两种更新与不同的测频协议、`drive_policy` 组合起来,整个闭环标定对应一个四阶段
状态机。编排由 `FrequencyCalibrationWorkflow` 负责,每个阶段是一个
`SinglePointFrequencyCalibration`,由 workflow 按序串联(每段以前一段的最优磁通和频率估计作初值):

```text
COARSE_ACQUIRE
    宽范围协议(Ramsey)得到初始频率
    drive_policy="sweet",仅更新磁通电压
        ↓  (残余失谐进入线性窗)
TRACKING
    驱动频率跟随预测的 qubit 频率
    局部协议(瞬态)高精度测频
    梯度下降更新磁通电压
    drive_policy="track"
        ↓  (误差降到目标容差附近)
LOCKED
    驱动频率固定为目标频率
    局部协议直接测量目标误差
    仅更新磁通电压(去抖收敛)
    drive_policy="target"
        ↓  (失谐越界 |δ|>ρ·Δ_lin 或测量失效)
REACQUIRE
    重新捕获:回到宽范围阶段重估频率
```

- **COARSE_ACQUIRE → TRACKING → LOCKED** 是正常前进路径,依次由 `"sweet"` /
  `"track"` / `"target"` 三种 `drive_policy` 驱动;每阶段自身的 `epsilon_f` 达标即
  交棒下一段。
- **REACQUIRE** 不是自动跳转,而是**编排层的响应**:当 `linear_range`
  ($\Delta_\mathrm{lin}$)配 `rho`(默认 0.6)在 `out_of_range` 信号里报告
  $|\delta|>\rho\,\Delta_\mathrm{lin}$,或测量失效时,调用方监听该信号(经
  `stop_predicate`)插入一个新的宽范围阶段重估频率。

#### FrequencyCalibrationWorkflow

**构造**

`FrequencyCalibrationWorkflow(stages, seed_drive_from_prev=True)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `stages` | `list[CalibrationStage]` | 顺序执行的阶段列表 | —— |
| `seed_drive_from_prev` | bool | 把上一段的频率估计作为下一段的 `omega_d_seed` | `True` |

**方法**

- `run() -> CalibrationTable`：按序执行所有阶段,返回末段的 `CalibrationTable`,
  `fit_params` 含全状态机轨迹。

#### CalibrationStage

**构造**

`CalibrationStage(name, calibration, stop_predicate=None)`

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `name` | str | 阶段名(`"COARSE_ACQUIRE"`/`"TRACKING"`/`"LOCKED"`/`"REACQUIRE"`) |
| `calibration` | `SinglePointFrequencyCalibration` | 该阶段的标定实例(含 `drive_policy`、`measure_method` 等) |
| `stop_predicate` | callable | 阶段结束条件(默认:自身 `epsilon_f` 达标即停) |

## WaveformCalibration：波形/控制线标定

预畸变主线的统一入口。

**构造**

`WaveformCalibration(method="transfer_function", distortion=None, measurement_protocol=None, qubit=None, control_line=None, fit_type="single_exp", transfer_model=None, predistortion_method="auto", n_taps=72, regularization=1e-6)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `method` | str | 标定路径 | `"transfer_function"`(或 `"predistortion"`) |
| `distortion` | `DistortionModel` | 畸变模型(解析路径:直接取阶跃响应) | `None` |
| `measurement_protocol` | str | 量子仿真测量协议 | `None`(或 `"cryoscope"`/`"delay_ramsey"`/`"transient"`/`"pi_pulse"`) |
| `qubit` | `TransmonQubit` | 比特(测量路径需要) | `None` |
| `control_line` | `ControlLine` | 控制线(测量路径需要) | `None` |
| `fit_type` | str | 传函拟合类型 | `"single_exp"`(或 `"multi_exp"`/`"fir"`/`"iir"`) |
| `transfer_model` | `DistortionModel` | 已测传函(预畸变路径需要) | `None` |
| `predistortion_method` | str | 预畸变设计方法 | `"auto"` |
| `n_taps` | int | FIR 滤波器阶数 | `72` |
| `regularization` | float | 岭正则化参数 | `1e-6` |

**方法**

- `calibrate() -> CalibrationTable`：按 `method` 返回 `kind="transfer_function"` 或
  `"predistortion"` 的表。
- `to_distortion_model() -> DistortionModel`：把标定结果直接转成 `DistortionModel`。

**输出**

`CalibrationTable`,`kind` 取决于 `method`;`to_distortion_model()` 返回可注入控制线的
`DistortionModel`。

## PredistortionDesigner：逆滤波器设计器

独立的预畸变滤波器设计器,可单用或被 `WaveformCalibration` 调用。

**构造**

`PredistortionDesigner(method="auto", n_taps=72, regularization=1e-6)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `method` | str | 设计方法 | `"auto"`(或 `"fir_inverse"`/`"iir_inverse"`/`"frequency_inverse"`) |
| `n_taps` | int | FIR 滤波器阶数 | `72` |
| `regularization` | float | 岭正则化参数 | `1e-6` |

**方法**

- `design(transfer_model, dt) -> DistortionModel`：给定传函模型返回其逆模型。
  `"auto"` 对指数型模型走解析 IIR 逆,其余走频域反演。
- `predistort(target, transfer=..., 或 inverse_model=...) -> Waveform`：直接把预畸变
  施加到目标波形上。
- `check_pole_stability(b, a) -> bool`：检查 IIR 滤波器极点是否都在单位圆内(稳定)。

**输出**

`design()` 返回 `DistortionModel`(逆模型);`predistort()` 返回预畸变后的 `Waveform`。

## CalibrationScheduler：标定任务调度器

标定任务的调度器:维护带依赖的命名标定任务注册表,按序执行。

**构造**

`CalibrationScheduler()`

**方法**

- `register(name, cal_class, depends_on=...)`：注册任务。
- `register_defaults()`：装入标准注册表(如 `frequency_closed_loop` 依赖
  `flux_response_ramsey`、`waveform_predistortion` 依赖
  `waveform_transfer_function`)。
- `run(name, **kw)`：按名跑(自动注入依赖结果)。
- `run_next()`：跑下一个就绪任务。
- `get_result(name)` / `status()`：查询。

**输出**

`run()` 返回 `CalibrationTable`;`status()` 返回 `{"ready": [...], "running": ..., "completed": [...], "failed": [...]}`。

```{note}
Kelly 2018 的 DAG 自动化(`check_state → maintain → auto_calibrate`)目前是接口 stub,
留作未来实现。当前调度为手动控制。
```

## 最小用例

```python
import numpy as np
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.calibration import (
    CalibrationTable,
    WaveformCalibration,
    PredistortionDesigner,
    CalibrationScheduler,
)

# 1) CalibrationTable 插值查表
table = CalibrationTable(
    name="flux_response",
    kind="f_phi",
    inputs=np.linspace(-0.03, 0.03, 7),
    outputs=np.array([5.5, 5.7, 5.9, 6.0, 5.9, 5.7, 5.5]) * 2 * np.pi,
)
f_at_zero = table.evaluate(np.array([0.0]))           # 插值
phi_for_target = table.inverse(np.array([5.8 * 2 * np.pi]))  # 反插值

# 2) 测传函 → 设计预畸变
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)
wf_cal = WaveformCalibration(distortion=distortion, fit_type="single_exp")
tf_result = wf_cal.calibrate()                         # kind="transfer_function"
inv_model = wf_cal.to_distortion_model()               # 转成 DistortionModel

# 3) 独立设计逆滤波器
designer = PredistortionDesigner(method="auto")
inverse = designer.design(distortion, dt=1.0)          # dt=1 ns

# 4) 调度器执行默认注册表
scheduler = CalibrationScheduler()
scheduler.register_defaults()
print(scheduler.status())   # {"ready": ["flux_response_ramsey", ...], ...}
```

## 物理角色 / 扩展

- 本层对应真实控制室里的标定程序:投片后第一步测 $f(\Phi)$,然后整定工作频率,
  最后测控制线传函并烧入补偿滤波器。三步分别对应 `FluxResponseCalibration`→
  `SinglePointFrequencyCalibration`→`WaveformCalibration`。
- `FrequencyCalibrationWorkflow` 把 `SinglePointFrequencyCalibration` 编排为四阶段
  状态机(COARSE_ACQUIRE → TRACKING → LOCKED → REACQUIRE),自动在宽范围 Ramsey
  与高精度瞬态之间切换,收敛到 kHz 级。
- `CalibrationTable` 的 `evaluate`/`inverse` 是 {doc}`reconstruction` 层
  `inversion="calibration"` 路径的基础:重建器把标定表传入,在反演时调 `inverse(φ)`
  把测到的相位翻回磁通幅度。
- `CalibrationScheduler` 的依赖注入机制自动把 `FluxResponseCalibration` 的输出
  (磁通范围)传给 `SinglePointFrequencyCalibration`,省去手动传参。
- 要添加自定义标定工作流,继承 `Calibration` 抽象基类实现 `calibrate() -> CalibrationTable`,
  然后向 `CalibrationScheduler` 注册并声明依赖。完整扩展指南见 {doc}`../extending`。
