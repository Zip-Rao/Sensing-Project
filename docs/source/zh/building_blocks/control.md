# 控制(control)

## 这层提供什么

`control` 层定义全栈共用的**时域信号类型**,并提供构建实验控制序列所需的全部组件:
磁通轨迹(`FluxSignal`)、微波脉冲包络(`Pulse`)、多脉冲实验序列(`CompositePulse`)
及内置的比特门函数。

本层是信号类型的**定义处**。`Waveform` 是最底层的时域容器,`FluxSignal` 是它的子类(
samples 以 $\Phi_0$ 计)。{doc}`hardware` 层的 `ControlLine` / `DistortionModel` /
`TransferMatrix` 都以这两个类型为输入输出——它们在 hardware 层的 `import` 里声明,不
在这里重复定义,原因是 hardware 依赖 control,而非反过来。

本层的职责止于**脉冲设计**:把参数化的信号组装成 QuTiP list 格式 Hamiltonian,准备好
交给 {doc}`simulation` 层的 `HamiltonianBuilder` / `MesolveRunner` 做时间演化。
AWG 输出到芯片的传函畸变不在本层处理,那是 {doc}`hardware` 层的职责。

一条贯穿本层的约定:**所有时间轴以 ns 计,且从 `sqc.config.CONFIG.awg.dt` 派生**
(用 `np.arange`,不用 `np.linspace`)。自由演化间隔等内部时间轴由全局 `dt` 统一生成,
保证脉冲序列能无缝拼到仿真的统一时间栅格上。

## 类总览

| 类 / 函数 | 角色 | 说明 |
|---|---|---|
| `Waveform` | 信号基元 | 语义中立的时域信号容器 $(t\_list, samples)$,本层最底层类型 |
| `CompositeWaveform` | 组合信号 | 多个 `Waveform` 按时间顺序拼接 |
| `FluxSignal` | 磁通信号 | `Waveform` 子类,samples 以 $\Phi_0$ 计,支持 9 种类型化构造 |
| `CompositeSignal` | 组合磁通 | 多个 `FluxSignal` 顺序拼接 |
| `PulseBase` | 抽象基类 | 任意控制脉冲的公共契约,本层扩展点 |
| `Pulse` | 单脉冲 | 单个控制脉冲(实验室系/旋转系),以一个 `FluxSignal` 为 Rabi 包络 |
| `CompositePulse` | 组合脉冲 | 多个 `Pulse` 按时间拼接成序列 |
| `PulseSequence` | 序列容器 | 有序脉冲集合的 dataclass |
| `create_*` 函数 | 序列构建函数 | 生成 Ramsey / 回波 / Cryoscope 等标准脉冲序列 |
| `make_drag_envelope` | 包络生成函数 | 生成 DRAG 脉冲的 I/Q 包络对 |
| `ideal_*` / `simulate_*` | 两比特门 | iSWAP / CZ 的理想矩阵与仿真 |

## Waveform —— 时域信号基元

本层最底层的信号类型,`@dataclass`,**语义中立**——只是 $(t\_list, samples)$ 这个
容器,不带物理单位含义。它是 `FluxSignal` 的基类、`Pulse` 包络的载体,也是
{doc}`hardware` 层 `ControlLine` / `DistortionModel` 输入输出的统一类型。

**字段**:`t_list`(时间点 ns)、`samples`(各时刻信号值)、`metadata`。
两数组形状不一致时构造报错。

**属性/方法**:`duration`、`n_points`、`value_at(t)`(采样保持,越界返 0)、
`samples_on(t_global)`(线性插值到全局时间轴,窗外为 0)、`truncate(t_start, t_end)`
(返回**新**波形,窗外置零)、`copy()`、`plot(ax=None)`。

`CompositeWaveform` 是其子类,持 `components` 列表;类方法
`CompositeWaveform.from_components(components)` 把多个波形按时间顺序拼接。

## FluxSignal —— 磁通信号

`Waveform` 的子类,**samples 以 $\Phi_0$ 为单位**,是磁通轨迹的物理化表达。它兼容
旧 `src.signal.Signal`,支持按 `type` 参数化构造(`Signal = FluxSignal` 别名保留)。

**构造**:`FluxSignal(type=0, t_list=None, trigger=0.0, **kwargs)`。`type` 取值:

| type | 波形 | type | 波形 |
|---|---|---|---|
| 0 | 零信号 | 5 | 双峰 |
| 1 | 常数 | 6 | 基函数展开(B-spline/Fourier/Legendre) |
| 2 | 正弦 | 7 | 复波包 |
| 3 | 高斯 | 8 | 用户自定义原始采样 |
| 4 | 非对称脉冲(双指数) | | |

**参数**(按 type 取用):`amplitude`、`frequency`、`phase`、`center`、`width`、
`offset`、`rise`/`fall`(type 4)、`noise_level`/`seed`;type 6 另有 `n_basis`、
`b`(基系数)、`basis_type`。type 8 必须传 `signal=<数组>`。

**兼容属性/方法**:`signal`(=`samples` 别名,可读写)、`type`、`params`、
`basis_functions`;`value_at(t)`、`samples_on(t_global)`(含 `trigger` 偏移)、
`truncate(t_start, t_end)`(**就地**置零,与基类返回新对象不同)、`update_signal(**kwargs)`
(改参数并就地重生成)、`copy()`(保留当前 samples)、`plot()`。

`CompositeSignal(signals)` 把多个 `FluxSignal` 顺序拼接(`CompositeWaveform` 子类)。

### make_drag_envelope —— DRAG 包络生成函数

`make_drag_envelope(t_list, amplitude, sigma, beta, phi=0.0) -> (FluxSignal, FluxSignal)`
—— 生成 DRAG 脉冲的 I/Q 包络对 $(\Omega_I, \Omega_Q)$。`beta` 取 $-1/\alpha$
(非谐性倒数),`phi=0` 给 X 旋转、`phi=\pi/2` 给 Y 旋转。返回两个 `type=8` 的
`FluxSignal`,可直接喂给 `Pulse` 的 `Omega` / `Omega_Q`。

## PulseBase —— 脉冲抽象基类

任意控制脉冲的公共契约,是本层的**扩展点**。子类须提供三个实例属性:
`hamiltonian`(QuTiP list 格式时变 Hamiltonian)、`t_list`(时间轴)、
`frame`(参考系,0=实验室系,1=旋转系)。要加自定义脉冲类型,继承 `PulseBase`
并在 `__init__` 里备好这三者,详见 {doc}`../extending`。

## Pulse —— 单控制脉冲

单个控制脉冲。它以一个 `FluxSignal`(或标量)为 Rabi 包络 `Omega`,构建
list 格式 Hamiltonian——**实验室系**(`frame=0`)或**旋转系**(`frame=1`,可选 RWA)。

**构造**:`Pulse(frame=0, omega_d=0.0, phase=0.0, Omega=None, Omega_Q=None,
is_rwa=True, qubit=None, trigger=0.0)`。`omega_d` 驱动频率(rad·GHz)、`phase` 旋转轴
相位、`Omega_Q` 是 DRAG 的 Q 分量(`None` 则单正交)、`trigger` 序列内绝对起始时刻。
构造后 `hamiltonian` / `t_list` 即备好。

**方法**

- `get_Rabi_frequency(t) -> float` / `get_hamiltonian_at(t) -> Qobj` —— 单时刻 Rabi 值 / Hamiltonian。
- `get_hamiltonian() -> (H_list, t_list)` —— 本地时间轴上的 list 格式 Hamiltonian。
- `hamiltonian_on(t_global) -> list` —— 把本地 Hamiltonian 投影到全局时间轴(含 `trigger` 偏移,窗外为 0),这是拼接到仿真统一栅格的规范通路。
- `get_angle(qubit=None) -> (theta, phi)` —— 数值演化求脉冲对比特的旋转角与相位;`get_angle_simple() -> float` —— 简单 Rabi 角 $\int\Omega\,\mathrm{d}t$。
- `with_phase_shift(phi_z) -> Pulse` —— 返回相位加 `phi_z` 的**新** `Pulse`(虚拟 Z 门)。

```{note}
`Pulse.get_kernel()` 已**弃用**,转发到 {doc}`reconstruction` 层的 `KernelEstimator`。
新代码请直接用 `KernelEstimator`。
```

## CompositePulse —— 组合脉冲序列

多个 `Pulse` 按时间拼接成一条序列。**构造**:`CompositePulse(pulses)`。
`frame` / `omega_d` 取首个子脉冲的值。

**方法**:`get_Omega(t)` / `get_hamiltonian_at(t)` —— 全序列在时刻 t 的 Rabi 值 /
Hamiltonian;`hamiltonian_on(t_global)` —— 汇总各子脉冲(各带自己的 `trigger`)在全局
轴上的贡献;`with_phase_shift(phi_z, from_time=None)` —— 对 `trigger >= from_time` 的
子脉冲加相位,返回**新**序列;`plot()`。`get_kernel()` 同样已弃用(转发 `KernelEstimator`)。

## PulseSequence —— 序列容器

`@dataclass`,有序脉冲集合:字段 `pulses`(脉冲列表)、`metadata`。

## 序列构建函数

一组用于构建标准实验序列的函数,返回 `CompositePulse`(`create_pulse` 返回
`QobjEvo`)。所有自由演化间隔都用全局 `dt` 生成;每个子脉冲带绝对 `trigger`,
可投影到仿真时间轴。

- `create_pulse(qubit, frame, type, t_list, omega_d, phase, angle=None, trigger=0.0, Omega_Q=None, **kwargs) -> QobjEvo`
  —— 单脉冲,给定目标 `angle` 时自动标定幅度(DRAG 时按同一因子缩放 Q 包络)。
- `create_ramsey_pulse(t_rabi, tau, omega_d=0.0, phase1=π/2, phase2=0.0, qubit=None, trigger=0.0)`
  —— Ramsey:$\pi/2 - \tau - \pi/2$。
- `create_echo_pulse(t_rabi, tau, omega_d=0.0, phase1=0, phase2=0, phase3=0, trigger=0.0)`
  —— 自旋回波:$\pi/2 - \tau - \pi - \tau - \pi/2$。
- `create_diff_echo_pulse(t_rabi, tau, t_int, t_rep, k, omega_d, ...)` —— 差分回波(k 次重复)。
- `create_cpmg_pulse(t_rabi, tau, n, omega_d, phase1=0.0, phase2=π/2, phase3=π/2, trigger=0.0)`
  —— CPMG 序列:$\pi/2 - [\tau/2 - \pi - \tau - \pi - \tau/2] - \pi/2$,含 `n` 个 π 脉冲。
- `create_cryoscope_pulse(t_rabi, tau, omega_d, phase1=π/2, phase2=0.0, trigger=0.0)`
  —— Cryoscope:$\pi/2 - \tau - \pi/2$(用于波形重建主线)。
- `create_pi_pulse_compensation_pulse(t_rabi, T_pi, omega_d, phase=0.0, qubit=None, trigger=0.0)`
  —— π-脉冲补偿协议的微波驱动(补偿磁通经 `qubit_in_mag` 单独施加,不在此序列内)。

## 两比特门

两比特门的理想矩阵与仿真函数($|00\rangle,|01\rangle,|10\rangle,|11\rangle$ 基):

- `ideal_iSWAP() -> Qobj` / `ideal_CZ() -> Qobj` —— 理想 $4\times4$ 幺正矩阵。
- `simulate_iSWAP(qubit1, qubit2, g) -> (U_eff, leakage, F)` —— 交换耦合 `g`(GHz)下
  iSWAP 仿真,返回有效幺正、泄漏、平均门保真度。
- `simulate_CZ(qubit1, qubit2, g) -> (U_eff_corr, leakage, F)` —— 磁通脉冲实现 CZ,
  含相位校正。

```{note}
两比特门的完整器件模型(可调耦合、多模腔)见 {doc}`devices` 的 `CoupledSystem`。
```

## 最小用例

```python
import numpy as np
from sqc.config import CONFIG
from sqc.devices import TransmonQubit
from sqc.control import FluxSignal, create_ramsey_pulse

# 1) 直接构造一个磁通信号:高斯磁通脉冲
t = CONFIG.pulse.t_signal.copy()
flux = FluxSignal(type=3, t_list=t, amplitude=0.02, center=t[-1] / 2, width=20.0)
print(flux.samples.max())              # 峰值磁通(Φ_0)

# 2) 搭一条 Ramsey 序列(π/2 - τ - π/2),返回 CompositePulse
qubit = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10_000, T2=8_000)
t_rabi = CONFIG.pulse.t_rabi.copy()
ramsey = create_ramsey_pulse(t_rabi, tau=50.0, omega_d=qubit.frequency, qubit=qubit)

# 3) 投影到全局时间轴 → 交给 simulation 层构 Hamiltonian
t_global = CONFIG.pulse.t_global.copy()
H_ctrl = ramsey.hamiltonian_on(t_global)   # QuTiP list 格式,喂给 mesolve
```

## 物理角色 / 扩展

- `FluxSignal` 对应真实施加到 Z 线上的**磁通轨迹**(理想、未畸变的那份);它经
  {doc}`hardware` 层 `ControlLine` / `DistortionModel` 后,才变成芯片上真实的波形。
- `Pulse` / `CompositePulse` 对应微波 XY 驱动:Rabi 包络、DRAG 校正、多脉冲实验序列。
- 序列构建函数对应实验中常用的标准脉冲序列(Ramsey、回波、Cryoscope 等)。
- **要加自定义脉冲类型**,继承 `PulseBase` 提供 `hamiltonian` / `t_list` / `frame`;
  **要加自定义波形类型**,继承 `Waveform`(如 `FluxSignal` 那样带物理单位)。完整扩展
  指南见 {doc}`../extending`。
