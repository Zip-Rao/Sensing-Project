# 控制(control)

## 这层提供什么

`control` 层定义全栈共用的时域信号类型,并提供构建实验控制序列所需的全部组件:
磁通轨迹(`FluxSignal`)、微波脉冲包络(`Pulse`)、多脉冲实验序列(`CompositePulse`)
及内置的比特门函数。

本层是信号类型的定义处。`Waveform` 是最底层的时域容器,`FluxSignal` 是它的子类(
samples 以 $\Phi_0$ 计)。{doc}`hardware` 层的 `ControlLine` / `DistortionModel` /
`TransferMatrix` 都以这两个类型为输入输出;它们在 hardware 层的 `import` 里声明,不
在这里重复定义,原因是 hardware 依赖 control,而非反过来。

本层的职责止于脉冲设计:把参数化的信号组装成 QuTiP list 格式 Hamiltonian,准备好
交给 {doc}`simulation` 层的 `HamiltonianBuilder` / `MesolveRunner` 做时间演化。
AWG 输出到芯片的传函畸变不在本层处理,那是 {doc}`hardware` 层的职责。

一条贯穿本层的约定:所有时间轴以 ns 计,且从 `sqc.config.CONFIG.awg.dt` 派生
(用 `np.arange`,不用 `np.linspace`)。自由演化间隔等内部时间轴由全局 `dt` 统一生成,
保证脉冲序列能拼到仿真的统一时间栅格上。

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

## Waveform：时域信号基元

本层最底层的信号类型,语义中立:只是 $(t\_list, samples)$ 这个容器,不带物理单位
含义。它是 `FluxSignal` 的基类、`Pulse` 包络的载体,也是 {doc}`hardware` 层
`ControlLine` / `DistortionModel` 输入输出的统一类型。

**构造**

`Waveform(t_list, samples, metadata=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `t_list` | `np.ndarray` | 时间点(ns) | —— |
| `samples` | `np.ndarray` | 各时刻信号值 | —— |
| `metadata` | dict | 附加元数据 | `None` |

两数组形状不一致时构造报错。

**只读属性**：`duration`(时长)、`n_points`(点数)。

**方法**

- `value_at(t) -> float`：采样保持取值,越界返 0。
- `samples_on(t_global) -> np.ndarray`：线性插值到全局时间轴,窗外为 0。
- `truncate(t_start, t_end) -> Waveform`：返回**新**波形,窗外置零。
- `copy() -> Waveform`：深拷贝。
- `plot(ax=None)`：快速可视化。

`CompositeWaveform` 是其子类,持 `components` 列表;类方法
`CompositeWaveform.from_components(components)` 把多个波形按时间顺序拼接。

## FluxSignal：磁通信号

`Waveform` 的子类,samples 以 $\Phi_0$ 为单位,是磁通轨迹的物理化表达。它兼容
旧 `src.signal.Signal`,支持按 `type` 参数化构造(`Signal = FluxSignal` 别名保留)。

**构造**

`FluxSignal(type=0, t_list=None, trigger=0.0, **kwargs)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `type` | int | 信号类型(0–8) | `0` |
| `t_list` | `np.ndarray` | 时间轴(ns) | —— |
| `trigger` | float | 触发偏移(ns) | `0.0` |
| `signal` | `np.ndarray` | 信号采样值(samples 别名,可读写) | —— |
| `params` | dict | 构造参数记录 | —— |

`type` 取值:

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

**方法**

- `value_at(t) -> float` / `samples_on(t_global) -> np.ndarray`：含 `trigger` 偏移。
- `truncate(t_start, t_end)`：**就地**置零(与基类返回新对象不同)。
- `update_signal(**kwargs)`：改参数并就地重生成。
- `copy() -> FluxSignal`：保留当前 samples 的深拷贝。
- `plot()`：快速可视化。

`CompositeSignal(signals)` 把多个 `FluxSignal` 顺序拼接(`CompositeWaveform` 子类)。

**输出**

本身即为信号对象;`samples_on()` 返回投影到全局轴上的 `np.ndarray`。

### make_drag_envelope：DRAG 包络生成函数

`make_drag_envelope(t_list, amplitude, sigma, beta, phi=0.0) -> (FluxSignal, FluxSignal)`
生成 DRAG 脉冲的 I/Q 包络对 $(\Omega_I, \Omega_Q)$。`beta` 取 $-1/\alpha$
(非谐性倒数),`phi=0` 给 X 旋转、`phi=\pi/2` 给 Y 旋转。返回两个 `type=8` 的
`FluxSignal`,可直接喂给 `Pulse` 的 `Omega` / `Omega_Q`。

## PulseBase：脉冲抽象基类

任意控制脉冲的公共契约,是本层的扩展点。子类须提供三个实例属性:
`hamiltonian`(QuTiP list 格式时变 Hamiltonian)、`t_list`(时间轴)、
`frame`(参考系,0=实验室系,1=旋转系)。要加自定义脉冲类型,继承 `PulseBase`
并在 `__init__` 里备好这三者,详见 {doc}`../extending`。

## Pulse：单控制脉冲

单个控制脉冲。它以一个 `FluxSignal`(或标量)为 Rabi 包络 `Omega`,构建 list
格式 Hamiltonian:实验室系(`frame=0`)或旋转系(`frame=1`,可选 RWA)。

**构造**

`Pulse(frame=0, omega_d=0.0, phase=0.0, Omega=None, Omega_Q=None, is_rwa=True, qubit=None, trigger=0.0)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `frame` | int | 参考系(0=实验室系,1=旋转系) | `0` |
| `omega_d` | float | 驱动频率(rad·GHz) | `0.0` |
| `phase` | float | 旋转轴相位 | `0.0` |
| `Omega` | `FluxSignal` 或 float | Rabi 包络(I 分量) | `None` |
| `Omega_Q` | `FluxSignal` 或 float | DRAG Q 分量(`None` 则单正交) | `None` |
| `is_rwa` | bool | 是否做旋波近似 | `True` |
| `qubit` | `TransmonQubit` | 比特(角度计算用) | `None` |
| `trigger` | float | 序列内绝对起始时刻(ns) | `0.0` |

构造后 `hamiltonian` / `t_list` 即备好。

**方法**

- `get_Rabi_frequency(t) -> float` / `get_hamiltonian_at(t) -> Qobj`：单时刻 Rabi 值 / Hamiltonian。
- `get_hamiltonian() -> (H_list, t_list)`：本地时间轴上的 list 格式 Hamiltonian。
- `hamiltonian_on(t_global) -> list`：把本地 Hamiltonian 投影到全局时间轴(含 `trigger` 偏移,
  窗外为 0),这是拼接到仿真统一栅格的规范通路。
- `get_angle(qubit=None) -> (theta, phi)`：数值演化求脉冲对比特的旋转角与相位;
  `get_angle_simple() -> float`：简单 Rabi 角 $\int\Omega\,\mathrm{d}t$。
- `with_phase_shift(phi_z) -> Pulse`：返回相位加 `phi_z` 的**新** `Pulse`(虚拟 Z 门)。

**输出**

`get_hamiltonian()` 返回 `(H_list, t_list)`;`hamiltonian_on()` 返回全局轴上的 list
格式 Hamiltonian(可直接喂给 `mesolve`);`with_phase_shift()` 返回新 `Pulse`。

```{note}
`Pulse.get_kernel()` 已**弃用**,转发到 {doc}`reconstruction` 层的 `KernelEstimator`。
新代码请直接用 `KernelEstimator`。
```

## CompositePulse：组合脉冲序列

多个 `Pulse` 按时间拼接成一条序列。

**构造**

`CompositePulse(pulses)`

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `pulses` | `list[Pulse]` | 子脉冲列表(按时间序) |
| `frame` | int | 参考系(取首个子脉冲的值) |
| `omega_d` | float | 驱动频率(取首个子脉冲的值) |

**方法**

- `get_Omega(t) -> float` / `get_hamiltonian_at(t) -> Qobj`：全序列在时刻 t 的 Rabi 值 /
  Hamiltonian。
- `hamiltonian_on(t_global) -> list`：汇总各子脉冲(各带自己的 `trigger`)在全局轴上的贡献。
- `with_phase_shift(phi_z, from_time=None) -> CompositePulse`：对 `trigger >= from_time`
  的子脉冲加相位,返回**新**序列。
- `plot()`：可视化各子脉冲包络。

`get_kernel()` 同样已弃用(转发 `KernelEstimator`)。

**输出**

`hamiltonian_on()` 返回拼接后的全局 list 格式 Hamiltonian;`with_phase_shift()` 返回
新 `CompositePulse`。

## PulseSequence：序列容器

有序脉冲集合的 dataclass。

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `pulses` | `list` | 脉冲列表 |
| `metadata` | dict | 附加元数据 |

## 序列构建函数

一组用于构建标准实验序列的函数,返回 `CompositePulse`(`create_pulse` 返回
`QobjEvo`)。所有自由演化间隔都用全局 `dt` 生成;每个子脉冲带绝对 `trigger`,
可投影到仿真时间轴。

- `create_pulse(qubit, frame, type, t_list, omega_d, phase, angle=None, trigger=0.0, Omega_Q=None, **kwargs) -> QobjEvo`
  单脉冲,给定目标 `angle` 时自动标定幅度(DRAG 时按同一因子缩放 Q 包络)。
- `create_ramsey_pulse(t_rabi, tau, omega_d=0.0, phase1=π/2, phase2=0.0, qubit=None, trigger=0.0, rotation_angle=None, rabi_rate=None, envelope="square", envelope_sigma=None) -> CompositePulse`
  可参数化 Ramsey 序列。默认仍为$\pi/2-\tau-\pi/2$方波；`rotation_angle`
  按目标转角自动归一化包络面积，`rabi_rate`固定包络峰值并由脉冲时长决定实际转角，
  二者互斥。`envelope`可取`"square"`、`"gaussian"`或与`t_rabi`同形的自定义数组。

```{note}
固定脉冲时长扫描`rotation_angle`会改变Rabi幅值；验证小转角缩短控制时间时，应固定
`rabi_rate`并同步改变`t_rabi`。每一种控制配置都需要重新估计响应核。
```
- `create_echo_pulse(t_rabi, tau, omega_d=0.0, phase1=0, phase2=0, phase3=0, trigger=0.0) -> CompositePulse`
  自旋回波:$\pi/2 - \tau - \pi - \tau - \pi/2$。
- `create_diff_echo_pulse(t_rabi, tau, t_int, t_rep, k, omega_d, ...) -> CompositePulse`
  差分回波(k 次重复)。
- `create_cpmg_pulse(t_rabi, tau, n, omega_d, phase1=0.0, phase2=π/2, phase3=π/2, trigger=0.0) -> CompositePulse`
  CPMG 序列:$\pi/2 - [\tau/2 - \pi - \tau - \pi - \tau/2] - \pi/2$,含 `n` 个 π 脉冲。
- `create_cryoscope_pulse(t_rabi, tau, omega_d, phase1=π/2, phase2=0.0, trigger=0.0) -> CompositePulse`
  Cryoscope:$\pi/2 - \tau - \pi/2$(用于波形重建主线)。
- `create_pi_pulse_compensation_pulse(t_rabi, T_pi, omega_d, phase=0.0, qubit=None, trigger=0.0) -> CompositePulse`
  π-脉冲补偿协议的微波驱动(补偿磁通经 `qubit_in_mag` 单独施加,不在此序列内)。

## 两比特门

两比特门的理想矩阵与仿真函数($|00\rangle,|01\rangle,|10\rangle,|11\rangle$ 基):

- `ideal_iSWAP() -> Qobj` / `ideal_CZ() -> Qobj`：理想 $4\times4$ 幺正矩阵。
- `simulate_iSWAP(qubit1, qubit2, g) -> (U_eff, leakage, F)`：交换耦合 `g`(GHz)下
  iSWAP 仿真,返回有效幺正、泄漏、平均门保真度。
- `simulate_CZ(qubit1, qubit2, g) -> (U_eff_corr, leakage, F)`：磁通脉冲实现 CZ,
  含相位校正。
- `simulate_cz_from_flux(t, phi_chip, qubit1, qubit2, g, ...) -> CZResult`
  (v2.6)：波形驱动 CZ —— 接收**显式片上磁通波形** ``phi_chip(t)`` (单位 Φ₀),
  经 Transmon 色散关系 $\omega(\Phi) = \sqrt{8 E_C E_J|\cos(\pi\Phi)|} - E_C$
  转换为失谐,构建完整时变两量子比特 Hamiltonian,从演化算符计算所有门指标。
  返回 ``CZResult``,包含 ``conditional_phase``、``phase_error``、``leakage``、
  ``infidelity`` 以及关键非计算态布居 (``pop_20``, ``pop_02``)。支持
  ``store_trajectories=True`` 记录相位和泄漏的时间轨迹。这使得**预失真→CZ 门**
  的端到端验证成为可能：将经控制线畸变(及校正后)的波形直接送入门仿真。

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
- 要加自定义脉冲类型,继承 `PulseBase` 提供 `hamiltonian` / `t_list` / `frame`;
  要加自定义波形类型,继承 `Waveform`(如 `FluxSignal` 那样带物理单位)。完整扩展
  指南见 {doc}`../extending`。
