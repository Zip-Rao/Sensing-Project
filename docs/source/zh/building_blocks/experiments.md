# 实验(experiments)

## 这层提供什么

`experiments` 层是全栈的协议编排层,把下面几层的组件(器件、磁通信号、
控制脉冲、求解器、读出模型)按某个传感协议组装成一次完整、可运行的实验。
每个实验类回答同一个问题的不同版本:给定一组脉冲序列与参数扫描,比特末态
布居或相位如何随参数变化。它自己不做数值积分(交给 {doc}`simulation` 层的
`mesolve`),也不做场重建(交给 {doc}`reconstruction` 层),只负责编排并产出
原始测量数据。

一条贯穿本层的约定:每个实验都是一个 `@dataclass`,构造时传入器件与可选参数,
由 `run()` 一次执行并返回结果。绝大多数实验返回统一的
{py:class}`~sqc.simulation.ExperimentResult`(`data`/`axes`/`metadata`/`config`
四字典),可直接传给 reconstruction 层。所有默认参数都对齐旧 `src/protocal.py`
的对应 case,使物理行为可回归。

```{note}
唯一的例外:{py:class}`~sqc.experiments.RabiExperiment` 的 `run()` 返回原始
QuTiP `Result`(而非 `ExperimentResult`),以保持与旧 case 0 的返回类型一致。
其余 7 个实验都返回 `ExperimentResult`。
```

## 类总览

| 类 | 角色 | 产出关键量 |
|---|---|---|
| `Experiment` | 抽象基类(本层扩展点) | —— |
| `RabiExperiment` | Rabi 振荡:定幅驱动、扫时长 | `p_e(t)`(原始 QuTiP `Result`) |
| `RamseyExperiment` | Ramsey:$\pi/2-\tau-\pi/2$,测自由演化累积相位 | `p_e(\tau)` |
| `DiffEchoExperiment` | 差分回波:$k$ 次重复累积相位,增敏弱信号 | `p_e(\tau)` |
| `TransientSensingExperiment` | 瞬态场感知:滑动测量 + 控制核 | `delta_p`、`kernel` |
| `CryoscopeExperiment` | Cryoscope:扫截断延迟,IQ 读出测 $\varphi(t_d)$ | `varphi(t_d)` |
| `DelayRamseyExperiment` | 延迟 Ramsey:滑短 Ramsey 过下降沿测拖尾 | `varphi(t_d)` |
| `PiPulseCompensationExperiment` | π-脉冲补偿:2D 扫 $(\tau, z)$ 测拖尾 | `z_star(\tau)` |

这些实验对应三条产品主线:波形重建(Ramsey / DiffEcho / Transient /
Cryoscope 恢复 $\Phi(t)$)、频率标定(Ramsey `f(\Phi)`)、预畸变
(DelayRamsey / PiPulseComp 测 AWG→芯片拖尾)。各自的端到端管道见
{doc}`../examples/waveform_reconstruction`、{doc}`../examples/frequency_calibration`、
{doc}`../examples/predistortion`。

## Experiment：实验抽象基类

所有传感实验的公共契约,是本层的扩展点。它规定两个抽象方法:

- `build_sequence()`:构建本实验的控制脉冲序列(有些实验按参数逐点构建,
  此方法返回 `None`)。
- `run(*args, **kwargs)`:执行整个实验并返回结果。

另提供一个共享辅助方法 `_route_flux(signal)`:把磁通信号经可选的
{py:class}`~sqc.hardware.ControlLine` 畸变(`control_line` 字段)后再耦合到比特;
`control_line=None`(默认)时直通。这是预畸变主线注入传函失真的统一入口。

要添加自定义传感协议,继承 `Experiment` 并实现 `build_sequence()` 与 `run()`,
`run()` 建议返回 `ExperimentResult`,详见 {doc}`../extending`。

## RabiExperiment：Rabi 振荡

施加定幅驱动脉冲,测激发态布居随时间的变化,用于标定 $\pi$/$\pi/2$ 脉冲时长。

**构造**

`RabiExperiment(qubit, t_rabi=None, omega_d=None, omega_rabi=1.0)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特 | —— |
| `t_rabi` | `np.ndarray` | Rabi 时间轴(ns) | `make_time(0, 40)` |
| `omega_d` | float | 驱动频率(rad·GHz) | 比特甜点频率 |
| `omega_rabi` | float | 定值 Rabi 频率 $\Omega$(rad·GHz) | `1.0` |

由于包络是固定幅度 $\Omega$,激发态布居呈现完整的 Rabi 振荡
$p_e = \sin^2(\Omega t / 2)$,在 `t_rabi` 上约含 $\Omega\,t_\mathrm{rabi}[-1] / (2\pi)$
个周期(默认参数下约 6 个完整的 0↔1 振荡)。增大 `omega_rabi` 可加快驱动、
增加周期数;$\pi$ 脉冲位于第一个 $p_e=1$ 峰,即 $t_\pi = \pi/\Omega$。

**方法**

- `build_sequence() -> Pulse`：返回定幅 {py:class}`~sqc.control.Pulse`(旋转系、RWA),
  其 `amplitude = omega_rabi`。
- `run() -> qutip.Result`：在全局时间轴 `CONFIG.pulse.t_global` 上跑 `mesolve`,
  返回原始 QuTiP `Result`。

**输出**

返回 QuTiP `Result`(⚠️ **非** `ExperimentResult`,是本层唯一的例外)。
取 `result.expect[0]` 即得 $p_e(t)$。

## RamseyExperiment：Ramsey 干涉

$\pi/2-\tau-\pi/2$ 序列,通过测自由演化期间累积的相位反映外磁通 $\Phi(t)$ 引起的
比特频率漂移,是波形重建与频率标定两条主线的核心协议。

**构造**

`RamseyExperiment(qubit, flux_signal=None, omega_d=None, tau_list=None, t_global=None, phase1=0.0, phase2=0.0, rotation_angle=π/2, rabi_rate=None, envelope="square", envelope_sigma=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特 | —— |
| `flux_signal` | `FluxSignal` | 施加的磁通信号 | 常值信号 |
| `omega_d` | float | 驱动频率(rad·GHz) | 比特甜点频率 |
| `tau_list` | `np.ndarray` | 自由演化时间扫描(ns) | `CONFIG.pulse.tau_list` |
| `t_global` | `np.ndarray` | 全局时间轴(ns) | `CONFIG.pulse.t_global` |
| `phase1` | float | 第一个 $\pi/2$ 脉冲相位 | `0.0` |
| `phase2` | float | 第二个 $\pi/2$ 脉冲相位 | `0.0` |
| `rotation_angle` | float 或 None | 每个控制脉冲的积分转角(rad) | $\pi/2$ |
| `rabi_rate` | float 或 None | 固定峰值 Rabi 速率；与 `rotation_angle` 互斥 | `None` |
| `envelope` | str 或 array | 方波、高斯或自定义控制包络 | `"square"` |
| `envelope_sigma` | float 或 None | 高斯包络标准差(ns) | 脉冲时长的 $1/4$ |

**方法**

- `run() -> ExperimentResult`：先把磁通信号投影到全局时间轴、经 `qubit_in_mag`
  耦合进比特(旋转系),再逐 `tau` 构建 Ramsey 序列跑演化,收集末态 $p_e$。

设置 `envelope="gaussian"` 可令两个控制脉冲采用高斯包络；实验类会把包络与
转角参数直接转发给 `create_ramsey_pulse()`。

**输出**

返回 `ExperimentResult`,其中:

| 键 | 含义 |
|---|---|
| `data["p_e"]` | $p_e(\tau)$ 数组 |
| `data["flux_samples"]` | 磁通采样数组 |
| `axes["tau"]` | 自由演化时间轴 |
| `axes["t_flux"]` | 磁通时间轴 |

## DiffEchoExperiment：差分回波

差分回波协议:序列 $\pi/2-[\tau-\pi-\tau'-(\tau+t_\mathrm{int})-\pi-\tau'']^k-\pi/2$,
通过 $k$ 次回波重复累积相位来增强弱信号灵敏度。

**构造**

`DiffEchoExperiment(qubit, flux_signal=None, k=5, t_rabi=None, t_int=None, t_rep=None, tau_list=None, t_global=None, omega_d=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特 | —— |
| `flux_signal` | `FluxSignal` | 施加的磁通信号 | 高斯信号 |
| `k` | int | 回波重复次数 | `5` |
| `t_rabi` | `np.ndarray` | Rabi 时间轴(ns) | `CONFIG.pulse.t_rabi` |
| `t_int` | float | 相互作用时间(ns) | —— |
| `t_rep` | float | 重复周期(ns) | —— |
| `tau_list` | `np.ndarray` | 自由演化扫描轴(ns) | `CONFIG.pulse.tau_list` |
| `t_global` | `np.ndarray` | 全局时间轴(ns) | `CONFIG.pulse.t_global` |
| `omega_d` | float | 驱动频率(rad·GHz) | 比特甜点频率 |

**方法**

- `run() -> ExperimentResult`：把基础磁通信号复制 $2k$ 份拼接、投影到全局轴耦合
  进比特,逐 `tau` 构建差分回波脉冲跑演化。

**输出**

返回 `ExperimentResult`,其中:

| 键 | 含义 |
|---|---|
| `data["p_e"]` | $p_e(\tau)$ 数组 |
| `axes["tau"]` | 自由演化时间轴 |
| `metadata["k"]` | 回波重复次数 |
| `metadata["t_int"]` | 相互作用时间 |
| `metadata["omega_d"]` | 驱动频率 |

## TransientSensingExperiment：瞬态场感知

瞬态磁场感知,基于滑动测量:让 Ramsey 控制脉冲(`tau=0`)在磁通信号上逐延迟
滑动,测每个延迟的 $p_e$,并对零磁通参考做差得 $\Delta p$,同时算出用于反卷积的
控制核。

**构造**

`TransientSensingExperiment(qubit, flux_signal=None, flux_signal_zero=None, t_rabi=None, omega_d=None, rotation_angle=π/2, rabi_rate=None, envelope="square", envelope_sigma=None, phase1=π/2, phase2=0.0, scan_list=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特 | —— |
| `flux_signal` | `FluxSignal` | 施加的磁通信号 | 非对称脉冲 |
| `flux_signal_zero` | `FluxSignal` | 零磁通参考信号 | 常零信号 |
| `t_rabi` | `np.ndarray` | Rabi 时间轴(ns) | `CONFIG.pulse.t_rabi` |
| `omega_d` | float | 驱动频率(rad·GHz) | 比特甜点频率 |
| `rotation_angle` | float或None | 每段控制脉冲的目标转角(rad) | $\pi/2$ |
| `rabi_rate` | float或None | 固定峰值Rabi频率；启用时转角由包络积分决定 | `None` |
| `envelope` | str或数组 | 方波、Gaussian或自定义包络 | `"square"` |
| `envelope_sigma` | float或None | Gaussian包络标准差(ns) | 脉冲时长的$1/4$ |
| `phase1`,`phase2` | float | 两段控制脉冲的旋转轴相位(rad) | $\pi/2,0$ |
| `scan_list` | `np.ndarray` | 延迟扫描轴(ns) | 自动生成 |

**方法**

- `run() -> ExperimentResult`：调 {py:class}`~sqc.simulation.SlidingMeasurementRunner`
  对信号与零参考各跑一次滑动测量,再经 {doc}`reconstruction` 层的 `KernelEstimator`
  算控制核。

**输出**

返回 `ExperimentResult`,其中:

| 键 | 含义 |
|---|---|
| `data["p_e"]` | 各延迟的激发态布居 |
| `data["delta_p"]` | $\Delta p$(相对零参考的差值) |
| `data["kernel"]` | 控制核数组 |
| `axes["scan"]` | 延迟扫描轴 |
| `axes["t_samples"]` | 核采样时间轴 |

$\Delta p$ 与 kernel 是波形重建中反卷积的输入。
结果的`metadata["rotation_angle"]`记录离散包络积分得到的实际转角。实验会针对所构造
的脉冲重新计算控制核，因此改变转角、时长、相位或包络后不应复用旧核。

## CryoscopeExperiment：Cryoscope

Cryoscope 协议:扫截断延迟 $t_d$,在每个截断处把磁通信号截断、耦合进比特,
用 IQ Ramsey 读出测累积相位 $\varphi(t_d)$;$\varphi$ 对 $t_d$ 求导得瞬时频率,
再经标定曲线反演得 $\Phi(t)$。

**构造**

`CryoscopeExperiment(qubit, flux_signal=None, t_rabi=None, tau=20.0, trunc_list=None, omega_d=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特 | —— |
| `flux_signal` | `FluxSignal` | 施加的磁通信号 | 正弦信号 |
| `t_rabi` | `np.ndarray` | Rabi 时间轴(ns) | `CONFIG.pulse.t_rabi` |
| `tau` | float | IQ 读出的自由进动时间(ns) | `20.0` |
| `trunc_list` | `np.ndarray` | 截断时间列表(逆序) | —— |
| `omega_d` | float | 驱动频率(rad·GHz) | 比特甜点频率 |

```{note}
`trunc_list` 中任何值超出 `flux_signal.t_list[-1]` 会在构造时抛 `ValueError`
(防止静默截断失败)。加长 `trunc_list` 时记得同步加长 `flux_signal` 的时间窗。
```

**方法**

- `run() -> ExperimentResult`：用 {py:class}`~sqc.hardware.IQReadoutModel` 逐截断点
  测 I/Q 投影,算原始相位 $\varphi=\arctan2(0.5-p_{e,I},\,p_{e,Q}-0.5)$,再用模型引导
  的解缠绕(锚定到理论累积相位 $\int_0^{t_d}(\omega_q(\Phi(t))-\omega_d)\,\mathrm{d}t$)
  得到连续 $\varphi(t_d)$。

**输出**

返回 `ExperimentResult`,其中:

| 键 | 含义 |
|---|---|
| `data["varphi"]` | 解缠绕后的累积相位 $\varphi(t_d)$ |
| `data["p_e_I"]` | I 通道激发态布居 |
| `data["p_e_Q"]` | Q 通道激发态布居 |
| `axes["trunc"]` | 截断时间轴(时间升序) |

## DelayRamseyExperiment：延迟 Ramsey

延迟 Ramsey(Ramsey 层析),用于测量磁通脉冲拖尾(预畸变主线):在下降沿后延迟
$t_d$ 处滑一段短 Ramsey,用 IQ 解调读相位 $\varphi(t_d)$,反演出拖尾波形。

**构造**

`DelayRamseyExperiment(qubit, flux_signal=None, t_d_list=None, tau_R=None, t_rabi=None, omega_d=None, t_fall=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特(需偏置在磁通敏感点) | —— |
| `flux_signal` | `FluxSignal` | 含拖尾的磁通信号 | 方波+指数拖尾 |
| `t_d_list` | `np.ndarray` | 下降沿后延迟列表(ns) | —— |
| `tau_R` | float | Ramsey 自由演化时间(ns) | —— |
| `t_rabi` | `np.ndarray` | Rabi 时间轴(ns) | `CONFIG.pulse.t_rabi` |
| `omega_d` | float | 驱动频率(rad·GHz) | 比特甜点频率 |
| `t_fall` | float | 下降沿时刻(ns),$t_d$ 相对它 | —— |

```{note}
`run_baseline` 字段已弃用:绝对相位参考现由模型引导解缠绕(见
{doc}`reconstruction`)解析给出,不再需要基线测量;字段仅为 API 兼容保留。
```

**方法**

- `run() -> ExperimentResult`：逐 $t_d$ 把拖尾片段窗在自由演化区间(两个 $\pi/2$
  脉冲期间置零以免失谐),经 IQ 读出测相位,再用与标定共享的模型引导解缠绕得
  $\varphi(t_d)$。

**输出**

返回 `ExperimentResult`,其中:

| 键 | 含义 |
|---|---|
| `data["varphi"]` | 解缠绕后的相位 $\varphi(t_d)$ |
| `data["varphi_raw"]` | 原始(未解缠绕)相位 |
| `data["p_e_I"]` | I 通道激发态布居 |
| `data["p_e_Q"]` | Q 通道激发态布居 |
| `axes["t_d"]` | 延迟时间轴 |

## PiPulseCompensationExperiment：π-脉冲补偿

π-脉冲补偿协议,同样用于测磁通脉冲拖尾:在下降沿后延迟 $\tau$ 处,同时施加一个
高度为 $z$ 的补偿磁通脉冲和一个 $\pi$ 脉冲;当 $z$ 恰好抵消拖尾时比特回到共振并翻转到
$|1\rangle$($p_e=1$)。对 $(\tau, z)$ 做 2D 扫描,取 $z^*(\tau)=\arg\max_z p_e$
得拖尾波形。

**构造**

`PiPulseCompensationExperiment(qubit, flux_signal=None, tau_list=None, z_list=None, T_pi=None, omega_bias=None, t_rabi=None, t_fall=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `qubit` | `TransmonQubit` | 被测比特 | —— |
| `flux_signal` | `FluxSignal` | 磁通信号(含下降沿拖尾) | —— |
| `tau_list` | `np.ndarray` | 下降沿后延迟列表(ns) | —— |
| `z_list` | `np.ndarray` | 补偿磁通高度($\Phi_0$) | `linspace(-0.01, 0.01, 21)` |
| `T_pi` | float | $\pi$ 脉冲时长(ns) | —— |
| `omega_bias` | float | 偏置频率(rad·GHz) | —— |
| `t_rabi` | `np.ndarray` | Rabi 时间轴(ns) | `CONFIG.pulse.t_rabi` |
| `t_fall` | float | 下降沿时刻(ns) | —— |

**方法**

- `run() -> ExperimentResult`：在 $(\tau, z)$ 网格上逐点跑 `mesolve` 得
  $p_e(\tau, z)$,再对每个 $\tau$ 用抛物线插值取峰(避免粗 `z_list` 的阶梯效应)
  得亚分辨率 $z^*(\tau)$。

**输出**

返回 `ExperimentResult`,其中:

| 键 | 含义 |
|---|---|
| `data["p_e"]` | 2D 数组 `(n_tau, n_z)`,各网格点的激发态布居 |
| `data["z_star"]` | 1D 数组,最优补偿高度(拖尾波形) |
| `axes["tau"]` | 延迟时间轴 |
| `axes["z"]` | 补偿高度轴 |

## 最小用例

```python
import numpy as np
from sqc.devices import TransmonQubit
from sqc.experiments import RabiExperiment, RamseyExperiment

qubit = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10_000, T2=8_000)

# 1) Rabi:标定 π 脉冲——返回原始 QuTiP Result
rabi = RabiExperiment(qubit=qubit).run()
p_e_t = rabi.expect[0]                        # p_e(t)

# 2) Ramsey:测自由演化累积相位——返回 ExperimentResult
ramsey = RamseyExperiment(
    qubit=qubit,
    tau_list=np.array([0.0, 50.0, 100.0]),    # 缩短扫描以加速演示
).run()
print(ramsey.data["p_e"])                     # p_e(τ)
print(ramsey.axes["tau"])                     # 自由演化时间轴
```

## 物理角色 / 扩展

- 每个实验类对应实验室里一套完整的脉冲序列与参数扫描方案:Rabi 标定门时长、
  Ramsey/回波做相位干涉、Cryoscope/延迟 Ramsey/π-脉冲补偿测波形与拖尾。
- 实验层负责编排,而非物理或数值:它调用 {doc}`control` 出脉冲、
  {doc}`simulation` 跑演化、{doc}`hardware` 的读出模型取信号,自己只负责把这些按
  协议串联起来并收集原始数据。
- 产出的 `ExperimentResult` 是 {doc}`reconstruction` / {doc}`calibration` 层的
  统一输入。
- 要添加自定义传感协议,继承 `Experiment` 抽象基类,实现 `build_sequence()` 与
  `run()`(建议返回 `ExperimentResult`),需要注入传函畸变时声明 `control_line`
  字段并用 `_route_flux()`。完整扩展指南见 {doc}`../extending`。
