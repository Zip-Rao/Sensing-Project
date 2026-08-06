# 仿真(simulation)

## 这层提供什么

`simulation` 层是全栈的求解器层,夹在 {doc}`control`(产出脉冲与 Hamiltonian
素材)和 {doc}`experiments`(编排完整协议)之间,只负责三件事:把器件与信号
组装成时变 Hamiltonian、驱动 QuTiP 求解器跑时间演化、把输出装进统一的结果
容器。它不含任何协议编排逻辑(那在 experiments 层),也不含任何重建/反演逻辑
(那在 {doc}`reconstruction` 层)。

一条贯穿本层的设计约束:Hamiltonian 构造是纯函数。`HamiltonianBuilder` 只读
`qubit` / `flux_signal` / `pulse`,输出 QuTiP list 格式的 `H_list` 与全局时间轴,
绝不就地修改任何输入对象。这样同一组器件与信号可被不同参考系、不同求解器反复
复用,不留副作用。

所有时间轴以 ns 计,且统一从 `sqc.config.CONFIG.awg.dt` 派生(见 {doc}`control`),
本层沿用同一栅格,保证脉冲、磁通、自由演化能无缝拼接。

## 类总览

| 类 / 函数 | 角色 | 说明 |
|---|---|---|
| `HamiltonianBuilder` | Hamiltonian 构造器 | 纯函数:器件+磁通+脉冲 → QuTiP list 格式 `H_list` 与全局时间轴 |
| `RunnerBase` | 抽象基类 | 所有求解器的公共契约,本层扩展点 |
| `MesolveRunner` | 主求解器 | 直接封装 `qutip.mesolve`,接收 `H_list` 跑主方程演化 |
| `SlidingMeasurementRunner` | 滑动测量求解器 | 让控制脉冲在磁通信号上滑动,逐延迟测激发态布居 $p_e$ |
| `ExperimentResult` | 结果容器 | 可序列化的 `data/axes/metadata/config`,带 `save`/`load` |
| `MeasurementTrace` | 单条测量迹 | frozen dataclass,一维扫描($axis, p_e$) |
| `extract_expectation` | 抽取函数 | 从 QuTiP result 取期望值数组 |
| `extract_population` | 抽取函数 | 从 `result.states` 取某 Fock 能级布居 |
| `generate_1f_noise` | 噪声生成 | FFT 法生成 $1/f$ 噪声时间序列 |

## HamiltonianBuilder：时变 Hamiltonian 构造器

本层的入口。把器件、磁通信号、控制脉冲组装成 `mesolve` 所需的 QuTiP list 格式
Hamiltonian。它是纯函数、无副作用:所有方法是 `@staticmethod`,不改任何输入。

**接口**

`HamiltonianBuilder.build(qubit, flux_signal=None, pulse=None, frame="rotating", omega_d=None) -> (H_list, t_global)`

| 参数 | 类型 | 含义 |
|---|---|---|
| `qubit` | `QubitSpec` 或 `TransmonQubit` | 器件参数(后者自动调 `.spec()`) |
| `flux_signal` | `FluxSignal` | 磁通信号;`None` 时退化为静态 Hamiltonian |
| `pulse` | `Pulse` 或 `CompositePulse` | 控制脉冲;非 `None` 时追加进 `H_list` |
| `frame` | str | `"lab"`(实验室系)或 `"rotating"`(旋转系) |
| `omega_d` | float | 驱动频率;默认取比特甜点频率 `spec.frequency()` |

**返回的 `H_list` 结构**(与 `src/qubit.py:qubit_in_mag` 逐位对齐):

$$H_0 = \tfrac{\alpha}{2}(n^2 - n), \qquad
H_1 = \big(n + \tfrac12 I\big)\cdot \omega_{01}(\Phi(t))$$

其中 $H_0$ 是常数系数项(非谐性),$H_1$ 的系数数组来自逐时刻代入磁通
$\omega_{01}(\Phi(t))=$ `spec.frequency(spec.flux_bias + s)`。

**输出**

返回 `(H_list, t_global)` 元组:`H_list` 可直接喂给 `MesolveRunner.run()`,
`t_global` 为全局时间轴。

```{note}
两种参考系都写成 `n + 0.5*I`。那个 `0.5*I` 是一个常数能量平移,不影响动力学,
保留它是为了与 `src/qubit.py:qubit_in_mag` 的输出逐位一致(回归基线依赖此)。
```

## RunnerBase：求解器抽象基类

所有仿真求解器的公共契约,是本层的扩展点。它只规定一个抽象方法:

- `run(*args, **kwargs)`:执行时间演化,返回一个 `ExperimentResult`。

要接入新的求解器(如 `sesolve`、蒙特卡洛轨迹、外部数值后端),继承 `RunnerBase`
并实现 `run()`、令其返回 `ExperimentResult` 即可,详见 {doc}`../extending`。

## MesolveRunner：主方程求解器

本层的主力求解器,直接封装 `qutip.mesolve`。

**构造**

`MesolveRunner(options=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `options` | dict 或 None | 传给 `mesolve` 的选项字典 | `None`(用默认) |

**方法**

- `run(H_list, psi0, t_list, c_ops, e_ops, store_states=False) -> ExperimentResult`：

| 参数 | 类型 | 含义 |
|---|---|---|
| `H_list` | list | `HamiltonianBuilder.build` 产出的 list 格式 Hamiltonian |
| `psi0` | `Qobj` | 初态 |
| `t_list` | `np.ndarray` | 时间轴(ns) |
| `c_ops` | `list[Qobj]` | 坍缩算符(常取 `qubit.get_collapse_operators()`) |
| `e_ops` | `list[Qobj]` | 期望算符(如投影算符) |
| `store_states` | bool | 是否保留每一步的态(供 `extract_population` 用) |

内部用 `QobjEvo(H_list, tlist=t_list, order=1)` 包裹 `H_list`。

**输出**

返回 `ExperimentResult`,其中:
- `data["expect"]`：期望值数组(形状 `(len(t_list),)`)
- `data["states"]`：态列表(未开 `store_states` 时为 `None`)
- `axes["t"]`：时间轴

## SlidingMeasurementRunner：滑动测量求解器

把控制脉冲在磁通信号上逐延迟滑动,在每个延迟点跑一次演化、测末态激发态布居
$p_e$,得到 $p_e$ 随延迟变化的一维迹。它内化了旧
`src/protocal.py:Protocal.sliding_measrement` / `single_measurement` 的逻辑
(旧拼写 `sliding_measrement` 仅在 `src_mirror` 兼容层保留)。

**构造**

`SlidingMeasurementRunner(options=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `options` | dict 或 None | 传给 `mesolve` 的选项字典 | `None`(用默认) |

**方法**

- `run(qubit, phi_signal, control_pulse, scan_list=None) -> ExperimentResult`：

| 参数 | 类型 | 含义 |
|---|---|---|
| `qubit` | `TransmonQubit` | 比特(需支持 `qubit_under_mag`) |
| `phi_signal` | `FluxSignal` | 磁通信号(需有 `.t_list` 与 `.value_at()`) |
| `control_pulse` | `CompositePulse` | 控制脉冲 |
| `scan_list` | `np.ndarray` | 延迟扫描轴;`None` 时自动生成 |

**输出**

返回 `ExperimentResult`,其中:
- `data["p_e"]`：各延迟的激发态布居
- `axes["scan"]`：延迟轴
- `metadata["elapsed"]`：耗时

```{note}
`SlidingMeasurementRunner` 内部对每个延迟点用收紧的容差
(`atol=1e-10, rtol=1e-8`)跑 `mesolve`;样本数多时耗时较长,`run()` 会打印总耗时。
```

## 结果容器与抽取函数

### ExperimentResult

可序列化的结果容器。

**构造**

`ExperimentResult(data=None, axes=None, metadata=None, config=None)`

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `data` | dict | 具名数据数组(`expect`、`p_e`、`states` 等) |
| `axes` | dict | 具名坐标轴(`t`、`scan`、`tau` 等) |
| `metadata` | dict | 自由元数据(求解器名、耗时、协议类型等) |
| `config` | dict | 仿真配置快照(EC、EJ、n_levels 等) |

**方法**

- `save(path)`：pickle 落盘(protocol 5)。
- `ExperimentResult.load(path) -> ExperimentResult`(classmethod)：读回,类型不符时抛
  `TypeError`。

### MeasurementTrace

`@dataclass(frozen=True)`,表示单条一维扫描迹。

**字段**

| 字段 | 类型 | 含义 |
|---|---|---|
| `axis` | `np.ndarray` | 扫描轴 |
| `p_e` | `np.ndarray` | 激发态布居 |
| `p_e_iq` | `np.ndarray` 或 None | I/Q 投影(可选) |
| `metadata` | dict | 附加元数据 |

### 抽取函数

- `extract_expectation(result, e_ops_index=0) -> np.ndarray`：从 QuTiP
  `mesolve`/`sesolve` result 取第 `e_ops_index` 个期望值数组;无 `.expect` 时抛
  `ValueError`。
- `extract_population(result, level) -> np.ndarray`：从 `result.states`
  逐步投影取第 `level` 个 Fock 能级的布居(需求解时开 `store_states=True`);无
  `.states` 时抛 `ValueError`。

## generate_1f_noise：1/f 噪声生成

`generate_1f_noise(t_list, amplitude, f_min, f_max, seed=None) -> np.ndarray`
FFT 法生成实值 $1/f$ 噪声时间序列:在 $[f_{\min}, f_{\max}]$ 频段内按
$1/\sqrt{f}$ 谱密度填复高斯谱、强制厄米对称后逆变换取实部。`seed` 固定后可复现。
用于给磁通/频率注入 $1/f$ 类退相干噪声。

## 最小用例

```python
import numpy as np
from qutip import basis, num
from sqc.config import CONFIG
from sqc.devices import TransmonQubit
from sqc.control import FluxSignal
from sqc.simulation import HamiltonianBuilder, MesolveRunner, extract_expectation

# 1) 器件 + 磁通信号(高斯磁通脉冲,峰值 0.02 Φ_0)
qubit = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10_000, T2=8_000)
t = CONFIG.pulse.t_signal.copy()
flux = FluxSignal(type=3, t_list=t, amplitude=0.02, center=t[-1]/2, width=20.0)

# 2) 纯函数构造时变 Hamiltonian(旋转系),不改 qubit / flux
H_list, t_global = HamiltonianBuilder.build(qubit, flux_signal=flux, frame="rotating")

# 3) 交给主求解器跑演化,取 ⟨n⟩ 期望
n_lev = qubit.n_levels
res = MesolveRunner().run(
    H_list, psi0=basis(n_lev, 1), t_list=t_global,
    c_ops=qubit.get_collapse_operators(), e_ops=[num(n_lev)],
)
n_expect = extract_expectation(res)      # ⟨n⟩(t),形状 (len(t_global),)

res.save("run.pkl")                       # 结果落盘(可 ExperimentResult.load 读回)
```

## 物理角色 / 扩展

- `HamiltonianBuilder` 对应把物理器件与施加的磁通/微波场写成含时 Hamiltonian
  这一步,是仿真的建模环节,不涉及任何数值积分。
- `MesolveRunner` 对应开放系统的主方程时间演化(含 $T_1$/$T_2$ 耗散);
  `SlidingMeasurementRunner` 对应瞬态场感知里滑动测量窗那类扫描式实验。
- `ExperimentResult` 对应一次实验/仿真的原始输出记录,是 {doc}`experiments`、
  {doc}`reconstruction` 各层的统一数据接口。
- 要接入自定义求解器,继承 `RunnerBase` 实现 `run()` 并返回 `ExperimentResult`;
  完整扩展指南见 {doc}`../extending`。
