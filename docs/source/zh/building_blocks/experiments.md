# 实验(experiments)

## 这层提供什么

`experiments` 层是全栈的**协议编排层**——把下面几层的零件(器件、磁通信号、
控制脉冲、求解器、读出模型)按某个**传感协议**组装成一次完整的、可运行的实验。
每个实验类回答同一个问题的不同版本:"施加这样的脉冲序列、扫这样的参数,比特末态
布居/相位怎么随参数变化?"它自己不做数值积分(交给 {doc}`simulation` 层的
`mesolve`),也不做场重建(交给 {doc}`reconstruction` 层)——它只负责**编排与产出
原始测量数据**。

一条贯穿本层的约定:**每个实验都是一个 `@dataclass`,构造时给器件与可选参数,
`run()` 一次跑完、返回结果**。绝大多数实验返回统一的
{py:class}`~sqc.simulation.ExperimentResult`(`data`/`axes`/`metadata`/`config`
四字典),可直接传给 reconstruction 层。所有默认参数都对齐旧 `src/protocal.py`
的对应 case,保证物理行为可回归。

```{note}
唯一的例外:{py:class}`~sqc.experiments.RabiExperiment` 的 `run()` 返回原始
QuTiP `Result`(而非 `ExperimentResult`),以保持与旧 case 0 的返回类型一致。
其余 7 个实验都返回 `ExperimentResult`。
```

## 类总览

| 类 | 传感协议 | 产出关键量 |
|---|---|---|
| `Experiment` | 抽象基类(本层扩展点) | —— |
| `RabiExperiment` | Rabi 振荡:定幅驱动、扫时长 | `p_e(t)`(原始 QuTiP `Result`) |
| `RamseyExperiment` | Ramsey:$\pi/2-\tau-\pi/2$,测自由演化累积相位 | `p_e(\tau)` |
| `DiffEchoExperiment` | 差分回波:$k$ 次重复累积相位,增敏弱信号 | `p_e(\tau)` |
| `TransientSensingExperiment` | 瞬态场感知:滑动测量 + 控制核 | `delta_p`、`kernel` |
| `CryoscopeExperiment` | Cryoscope:扫截断延迟,IQ 读出测 $\varphi(t_d)$ | `varphi(t_d)` |
| `DelayRamseyExperiment` | 延迟 Ramsey:滑短 Ramsey 过下降沿测拖尾 | `varphi(t_d)` |
| `PiPulseCompensationExperiment` | π-脉冲补偿:2D 扫 $(\tau, z)$ 测拖尾 | `z_star(\tau)` |

这些实验对应三条产品主线:**波形重建**(Ramsey / DiffEcho / Transient /
Cryoscope 恢复 $\Phi(t)$)、**频率标定**(Ramsey `f(\Phi)`)、**预畸变**
(DelayRamsey / PiPulseComp 测 AWG→芯片拖尾)。各自的端到端管道见
{doc}`../examples/waveform_reconstruction`、{doc}`../examples/frequency_calibration`、
{doc}`../examples/predistortion`。

## Experiment —— 实验抽象基类

所有传感实验的公共契约,是本层的**扩展点**。它规定两个抽象方法:

- `build_sequence()` —— 构建本实验的控制脉冲序列(有些实验按参数逐点构建,
  此方法返回 `None`)。
- `run(*args, **kwargs)` —— 跑完整个实验并返回结果。

另提供一个共享辅助方法 `_route_flux(signal)`:把磁通信号经可选的
{py:class}`~sqc.hardware.ControlLine` 畸变(`control_line` 字段)后再耦合到比特;
`control_line=None`(默认)时直通。这是预畸变主线注入传函失真的统一入口。

要加自定义传感协议,继承 `Experiment` 实现 `build_sequence()` / `run()`,
`run()` 建议返回 `ExperimentResult`,详见 {doc}`../extending`。

## RabiExperiment —— Rabi 振荡

施加**定幅**驱动脉冲、测激发态布居随时间的变化,用于标定 $\pi$/$\pi/2$ 脉冲时长。
`@dataclass`,主要字段:`qubit`、`t_rabi`(Rabi 时间轴,默认 `make_time(0, 40)`)、
`omega_d`(驱动频率,默认取比特甜点频率)。

`build_sequence()` 返回一个定幅 {py:class}`~sqc.control.Pulse`(旋转系、RWA)。
`run()` 在全局时间轴 `CONFIG.pulse.t_global` 上跑 `mesolve`,**返回原始 QuTiP
`Result`**(见开篇 note),取 `result.expect[0]` 即得 $p_e(t)$。

## RamseyExperiment —— Ramsey 干涉

$\pi/2-\tau-\pi/2$ 序列,通过测自由演化期间累积的相位来反映外磁通 $\Phi(t)$ 引起的
比特频率漂移,是波形重建与频率标定两条主线的核心协议。`@dataclass`,主要字段:
`qubit`、`flux_signal`(默认常值信号)、`omega_d`、`tau_list`(自由演化时间扫描,
默认 `CONFIG.pulse.tau_list`)、`t_global`、`phase1`/`phase2`(两个 $\pi/2$ 脉冲相位)。

`run()` 先把磁通信号投影到全局时间轴、经 `qubit_in_mag` 耦合进比特(旋转系),
再逐 `tau` 构建 Ramsey 序列跑演化,收集末态 $p_e$。返回的 `ExperimentResult`:
`data["p_e"]` 为 $p_e(\tau)$、`data["flux_samples"]` 为磁通采样、`axes["tau"]`
为自由演化轴、`axes["t_flux"]` 为磁通时间轴。

## DiffEchoExperiment —— 差分回波

差分回波协议:序列 $\pi/2-[\tau-\pi-\tau'-(\tau+t_\mathrm{int})-\pi-\tau'']^k-\pi/2$,
通过 $k$ 次回波重复累积相位来增强弱信号灵敏度。`@dataclass`,主要字段:`qubit`、
`flux_signal`(默认高斯信号)、`k`(回波重复次数,默认 5)、`t_rabi`、
`t_int`(相互作用时间)、`t_rep`(重复周期)、`tau_list`、`t_global`、`omega_d`。

`run()` 把基础磁通信号复制 $2k$ 份拼接、投影到全局轴耦合进比特,逐 `tau` 构建
差分回波脉冲跑演化。返回的 `ExperimentResult`:`data["p_e"]` 为 $p_e(\tau)$、
`axes["tau"]`、`metadata` 记录 `k` / `t_int` / `omega_d`。

## TransientSensingExperiment —— 瞬态场感知

瞬态磁场感知,基于**滑动测量**:让 Ramsey 控制脉冲(`tau=0`)在磁通信号上逐延迟
滑动,测每个延迟的 $p_e$,并对零磁通参考做差得 $\Delta p$,同时算出用于反卷积的
控制核。`@dataclass`,主要字段:`qubit`、`flux_signal`(默认非对称脉冲)、
`flux_signal_zero`(零磁通参考,默认常零)、`t_rabi`、`omega_d`、`scan_list`。

`run()` 调 {py:class}`~sqc.simulation.SlidingMeasurementRunner` 对信号与零参考各跑
一次滑动测量,再经 {doc}`reconstruction` 层的 `KernelEstimator` 算控制核。返回的
`ExperimentResult`:`data["p_e"]`、`data["delta_p"]`、`data["kernel"]`、
`axes["scan"]`(延迟轴)、`axes["t_samples"]`(核采样轴)。$\Delta p$ 与 kernel
正是波形重建里反卷积的输入。

## CryoscopeExperiment —— Cryoscope

Cryoscope 协议:扫**截断延迟** $t_d$,在每个截断处把磁通信号截断、耦合进比特,
用 IQ Ramsey 读出测累积相位 $\varphi(t_d)$;$\varphi$ 对 $t_d$ 求导即得瞬时频率、
再经标定曲线反演得 $\Phi(t)$。`@dataclass`,主要字段:`qubit`、`flux_signal`
(默认正弦)、`t_rabi`、`tau`(IQ 读出的自由进动时间)、`trunc_list`(截断时间列表,
逆序)、`omega_d`。

`run()` 用 {py:class}`~sqc.hardware.IQReadoutModel` 逐截断点测 I/Q 投影,算原始相位
$\varphi=\arctan2(0.5-p_{e,I},\,p_{e,Q}-0.5)$,再用模型引导的解缠绕
(锚定到理论累积相位 $\int_0^{t_d}(\omega_q(\Phi(t))-\omega_d)\,\mathrm{d}t$)得到
连续 $\varphi(t_d)$。返回的 `ExperimentResult`:`data["varphi"]`、`data["p_e_I"]`、
`data["p_e_Q"]`、`axes["trunc"]`(时间升序)。

```{note}
`trunc_list` 中任何值超出 `flux_signal.t_list[-1]` 会在构造时抛 `ValueError`
(防止静默截断失败)。加长 `trunc_list` 时记得同步加长 `flux_signal` 的时间窗。
```

## DelayRamseyExperiment —— 延迟 Ramsey

延迟 Ramsey(Ramsey 层析),用于**测量磁通脉冲拖尾**(预畸变主线):在下降沿后延迟
$t_d$ 处滑一段短 Ramsey,用 IQ 解调读相位 $\varphi(t_d)$,反演出拖尾波形。
`@dataclass`,主要字段:`qubit`(需偏置在磁通敏感点)、`flux_signal`(默认方波+
指数拖尾)、`t_d_list`(下降沿后延迟)、`tau_R`(Ramsey 自由演化时间)、`t_rabi`、
`omega_d`、`t_fall`(下降沿时刻,$t_d$ 相对它)。

`run()` 逐 $t_d$ 把拖尾片段窗在自由演化区间(两个 $\pi/2$ 脉冲期间置零以免失谐),
经 IQ 读出测相位,再用与标定共享的模型引导解缠绕得 $\varphi(t_d)$。返回的
`ExperimentResult`:`data["varphi"]`、`data["varphi_raw"]`、`data["p_e_I"]`、
`data["p_e_Q"]`、`axes["t_d"]`。

```{note}
`run_baseline` 字段已弃用:绝对相位参考现由模型引导解缠绕(见
{doc}`reconstruction`)解析给出,不再需要基线测量;字段仅为 API 兼容保留。
```

## PiPulseCompensationExperiment —— π-脉冲补偿

π-脉冲补偿协议,同样用于**测磁通脉冲拖尾**:在下降沿后延迟 $\tau$ 处,同时施加一个
高度为 $z$ 的补偿磁通脉冲和一个 $\pi$ 脉冲;当 $z$ 恰好抵消拖尾时比特回到共振、翻转到
$|1\rangle$($p_e=1$)。对 $(\tau, z)$ 做 2D 扫描,取 $z^*(\tau)=\arg\max_z p_e$
即得拖尾波形。`@dataclass`,主要字段:`qubit`、`flux_signal`、`tau_list`、
`z_list`(补偿高度,默认 `linspace(-0.01, 0.01, 21)`)、`T_pi`、`omega_bias`、
`t_rabi`、`t_fall`。

`run()` 在 $(\tau, z)$ 网格上逐点跑 `mesolve` 得 $p_e(\tau, z)$,再对每个 $\tau$ 用
抛物线插值取峰(避免粗 `z_list` 的阶梯效应)得亚分辨率 $z^*(\tau)$。返回的
`ExperimentResult`:`data["p_e"]`(2D `(n_tau, n_z)`)、`data["z_star"]`(拖尾波形)、
`axes["tau"]`、`axes["z"]`。

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

- 每个实验类对应实验室里一套完整的**脉冲序列 + 参数扫描**方案:Rabi 标定门时长、
  Ramsey/回波做相位干涉、Cryoscope/延迟 Ramsey/π-脉冲补偿测波形与拖尾。
- 实验层是"编排"而非"物理"或"数值":它调用 {doc}`control` 出脉冲、
  {doc}`simulation` 跑演化、{doc}`hardware` 的读出模型取信号,自己只负责把这些按
  协议串起来并收集原始数据。
- 产出的 `ExperimentResult` 是 {doc}`reconstruction` / {doc}`calibration` 层的
  统一输入。
- **要加自定义传感协议**,继承 `Experiment` 抽象基类,实现 `build_sequence()` 与
  `run()`(建议返回 `ExperimentResult`),需要注入传函畸变时声明 `control_line`
  字段并用 `_route_flux()`。完整扩展指南见 {doc}`../extending`。
