# 器件(devices)

## 这层提供什么

`devices` 层是全栈 cQED 架构的**最底层**,只建模物理量子器件本身——超导
Transmon 比特、读出/耦合谐振腔、以及把它们组装起来的芯片拓扑与耦合系统。
这一层**只描述物理**(能量、频率、Hilbert 空间、Hamiltonian、耗散算符),
不涉及任何控制脉冲、时间演化或重建逻辑;那些职责在上面的 `control`、
`simulation` 等层。

一条贯穿本层的设计约束:**器件是"纯参数对象"**——它描述物理,但不存储实验状态。
要表示受扰动的器件(如处于非零磁通的比特),用返回**新对象**的方法
(`QubitSpec.with_flux`),而不是就地修改。`QubitSpec` 因此是 frozen dataclass。

## 类总览

| 类 | 角色 | 说明 |
|---|---|---|
| `Device` | 抽象基类 | 所有物理器件的公共契约 |
| `QubitSpec` | 参数内核 | 不可变的 Transmon 参数描述,`HamiltonianBuilder` 的输入 |
| `TransmonQubit` | 运行时器件 | 可调频 Transmon,持有 QuTiP 算符与 Hamiltonian,各上层默认接收的器件对象 |
| `Resonator` | 谐振腔 | 多模读出/耦合谐振腔 |
| `CoupledSystem` | 耦合系统 | 双比特经可调耦合器接多模腔,用于两比特门 |
| `ChipTopology` | 芯片拓扑 | 多比特+谐振腔容器,把单器件算符嵌入整芯片 Hilbert 空间 |

## 物理模型与单位

Transmon 物理(参照 Gao 2021 §II.B–C,Eq. 13–20):

$$H = 4 E_C n^2 - E_J(\Phi)\cos\varphi, \qquad
E_J(\Phi) = E_{J,0}\,\lvert\cos(\pi \Phi/\Phi_0)\rvert$$

$$\omega_{01} = \sqrt{8 E_J E_C} - E_C, \qquad \alpha = -E_C$$

**单位约定(重要)**:`EC`、`EJ` 以 **rad·GHz** 传入(已含 $2\pi$ 因子),
磁通以 $\Phi_0$ 为单位,时间以 ns 计,$\hbar = 1$。所以 $E_C/h = 200\,\mathrm{MHz}$
对应传参 `EC = 2*np.pi*0.2`;`frequency()` 返回的也是 rad·GHz(即**角频率**
$\omega_{01}$),除以 $2\pi$ 得到以 GHz 计的物理频率 $f_{01}$。全站默认器件参数
($E_C/h=200\,\mathrm{MHz}$、
$E_J/h=15\,\mathrm{GHz}$、$T_1=10\,\mu\mathrm{s}$、$T_2=8\,\mu\mathrm{s}$)集中在
`sqc.config.CONFIG`。

## Device —— 器件抽象基类

所有物理器件的公共契约,是本层的**扩展点**。子类必须是纯参数对象(只描述物理、
不存实验状态),并实现以下抽象成员:

- `name`(property)—— 人类可读的器件标识符。
- `hilbert_dim() -> int` —— 完整 Hilbert 空间维数。
- `hamiltonian_static() -> Qobj` —— 无驱动、无外场的**时间无关** Hamiltonian。
- `collapse_operators() -> list[Qobj]` —— 用于耗散/退相的 Lindblad 坍缩算符。

要加自定义器件类型(如 fluxonium),继承 `Device` 并实现上述 4 个成员即可,
详见 {doc}`../extending`。

## QubitSpec —— 不可变参数记录

一个可调频 Transmon 的**纯参数描述**,也是本层的**参数内核**:`@dataclass(frozen=True)`,
不持有任何实验状态,所有方法都是纯函数(取参数、返回值,不改 `self`)。`simulation`
层的 `HamiltonianBuilder` 直接以 `QubitSpec` 为输入构造仿真所需的 Hamiltonian
(传入 `TransmonQubit` 时,它会调 `.spec()` 取出内部的 `QubitSpec`)。

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `name` | str | 器件标识符 | —— |
| `EC` | float | 充电能(rad·GHz) | —— |
| `EJ` | float | 零磁通处约瑟夫森能(rad·GHz) | —— |
| `T1` | float | 弛豫时间(ns) | —— |
| `T2` | float | 退相时间(ns) | —— |
| `flux_bias` | float | 静态磁通偏置($\Phi_0$) | `0.0` |
| `n_levels` | int | Hilbert 空间截断维数 | `3` |

**方法**

- `EJ_at(flux=None) -> float` —— 给定磁通处的有效约瑟夫森能
  $E_J(\Phi)=E_{J,0}\lvert\cos(\pi\Phi/\Phi_0)\rvert$;`flux=None` 时用 `flux_bias`。
- `frequency(flux=None) -> float` —— 跃迁角频率 $\omega_{01}(\Phi)=\sqrt{8E_J(\Phi)E_C}-E_C$(rad·GHz)。
- `anharmonicity() -> float` —— 非谐性 $\alpha=-E_C$(rad·GHz)。
- `sensitivity(flux=None, delta=1e-6) -> float` —— 中心差分数值求 $\mathrm{d}\omega_{01}/\mathrm{d}\Phi$(rad·GHz/$\Phi_0$)。
- `with_flux(new_flux) -> QubitSpec` —— 返回 `flux_bias=new_flux` 的**新** `QubitSpec`,不改 `self`。
- `optimal_work_point() -> float`(staticmethod)—— $\lvert\mathrm{d}f/\mathrm{d}\Phi\rvert$ 最大处的磁通,单位 $\Phi_0$,值为 $\arctan(\sqrt2)/\pi$。

## TransmonQubit —— 可调频 Transmon 器件

进入实验流程的**运行时器件对象**。它内部持有一个 `QubitSpec`,并把 QuTiP 算符、
Hamiltonian 与坍缩算符都构造好,实现了 `Device` 接口(`name`/`hilbert_dim`/
`hamiltonian_static`/`collapse_operators`)。`experiments`、`control`、`calibration`、
`simulation` 各层**默认接收的正是 `TransmonQubit`**,`sqc.config.CONFIG` 产出的默认
比特也是此类型。相较之下,`QubitSpec` 只描述物理参数;要进入实验流程,用 `TransmonQubit`。

**构造参数**:`TransmonQubit(EC, EJ, T1, T2, flux=0.0, state=0, n_levels=3, name="Q")`。
注意这里参数名是 `flux`(不是 `flux_bias`),且可传初始 `state`(Fock 态索引或 `Qobj`)。

**关键属性**(构造后即可读):`EC`、`EJ`(当前磁通处的有效值)、`EJ_0`(零磁通值)、
`flux`、`n_levels`、`T1`、`T2`、`frequency`、`anharmonicity`、算符 `a`/`a_dag`/`n`/`I`、
`hamiltonian`、`c_ops`、`state`。

**访问器与 Device 接口**

- `spec() -> QubitSpec` —— 取内部 `QubitSpec`。
- `name`(property)、`hilbert_dim()`、`hamiltonian_static()`、`collapse_operators()` —— `Device` 契约。
- `calculate_frequency() -> float` / `calculate_anharmonicity() -> float` —— 由当前 `EJ, EC` 算 $\omega_{01}$ / $\alpha$(rad·GHz)。
- `frequency_sensitivity(flux, delta_flux=1e-6) -> float` —— 中心差分 $\mathrm{d}\omega_{01}/\mathrm{d}\Phi$(rad·GHz/$\Phi_0$)。

**Hamiltonian 与坍缩算符**

- `get_hamiltonian() -> Qobj` —— 实验室系时间无关 $H$。
- `get_hamiltonian_rwa(omega_d) -> Qobj` —— 旋转系(RWA)$H=\Delta n+\tfrac{\alpha}{2}(n^2-n)$,$\Delta=f_{01}-\omega_d$。
- `get_collapse_operators() -> list[Qobj]` —— 由 `T1, T2` 生成 $[\sqrt{\gamma_1}\,a,\ \sqrt{\gamma_\phi}\,n]$。

```{note}
把比特置入磁通信号、构造时变 Hamiltonian 是 `simulation` 层
`HamiltonianBuilder` 的职责,见 {doc}`simulation`。
```

**噪声、态投影与门仿真**

- `generate_1f_noise(t_lists, amplitude, f_min, f_max) -> np.ndarray` —— 生成 1/f 噪声时间序列(转调 `sqc.simulation.noise`)。
- `calculate_state_projection(target_state) -> float` —— 当前 `state` 落在 `target_state`(Fock 索引或 `Qobj`)上的概率。
- `ideal_gate(theta, phi) -> Qobj` —— $\{|0\rangle,|1\rangle\}$ 子空间的理想单比特旋转($2\times2$ 幺正)。
- `simulate_gate(theta, phi, T, sigma) -> (Qobj, float, float)` —— DRAG 单比特门仿真,返回 `(U_eff, leakage, fidelity)`。

## Resonator —— 多模谐振腔

多模读出/耦合谐振腔,继承 `Device`。每个模是一个 Fock 空间,整体用张量积构造算符。

**构造**:`Resonator(frequencies, n_levels)` —— `frequencies` 为各模频率(rad·GHz)列表,
`n_levels` 为各模 Fock 截断维数列表(两者等长)。

**属性**:`frequencies`、`M`(模数)、`n_levels`、模算符列表 `a`/`adag`/`n`、整体单位算符 `I`、`hamiltonian`。

**方法**

- `name`(property)、`hilbert_dim()`、`hamiltonian_static()`、`collapse_operators()`(当前返回空列表)—— `Device` 契约。
- `get_hamiltonian() -> Qobj` —— 线性 Hamiltonian $\sum_i f_i n_i$。
- `get_hamiltonian_rot(omega_d) -> Qobj` —— 旋转系 $\sum_i (f_i-\omega_d) n_i$。

## CoupledSystem —— 双比特耦合系统

两个 Transmon 经可调耦合接入多模谐振腔,继承 `Device`,用于两比特门仿真。

**构造**:`CoupledSystem(qubit1, qubit2, cavity)` —— 两个 `TransmonQubit` 加一个 `Resonator`。

**属性**:张量积算符 `a1`/`a2`/`acav`、耦合强度 `g`(2×M)、子 Hamiltonian `H_q1`/`H_q2`/`H_cav`/`H_0`、总 `H`。

**方法**

- `name`(property)、`hilbert_dim()`、`hamiltonian_static()`、`collapse_operators()` —— `Device` 契约。
- `initialize_g() -> list[list[float]]` —— 初始化默认比特-模耦合强度。
- `control_g(g, is_time_dependent=False, f=None) -> None` —— 设置(可含时变)比特-腔耦合并重建 `H`。
- `build_system(is_rwa=False, omega_d=None) -> None` —— 构建完整系统 Hamiltonian(可选 RWA)。
- `prepare_ket11(T, sigma) -> (Qobj, float)` —— 用 DRAG π 脉冲制备 $|11\rangle$,返回 `(final_state, fidelity)`。
- `simulate_iSWAP(g) -> (Qobj, float, float)` —— 可调耦合下 iSWAP 门仿真,返回 `(U_cal, leakage, F_loc)`。

## ChipTopology —— 多比特芯片拓扑

`@dataclass`,充当整芯片所有器件的容器,并负责把单器件算符**嵌入(lift)**到
完整 Hilbert 空间。继承 `Device`。

**字段**:`qubits`(比特列表)、`resonators`(谐振腔列表)、
`couplings`(`(qubit_name, resonator_name) -> g` 字典)、
`control_lines`(名字 → `ControlLine` 字典)、`transfer_matrix`(频率相关传递矩阵 `TransferMatrix`,可为 `None`)。

**方法**

- `name`(property)、`hilbert_dim()` —— 整芯片总维数(各子系统维数之积)。
- `lift_qubit_op(op, qubit_name) -> Qobj` —— 把单比特算符嵌入整芯片空间(目标槽放 `op`,其余放单位算符);找不到该比特则抛 `ValueError`。
- `hamiltonian_static() -> Qobj` —— 各器件静态 Hamiltonian 之和(注:基础实现**不含**耦合项,耦合系统请用 `CoupledSystem`)。
- `collapse_operators() -> list[Qobj]` —— 汇总并嵌入所有器件的坍缩算符。
- `get_qubit(name) -> TransmonQubit` —— 按名字取比特,找不到抛 `ValueError`。

## 最小用例

```python
import numpy as np
from sqc.devices import QubitSpec, TransmonQubit

# 首选:不可变参数记录(EC/EJ 以 rad·GHz 传入,含 2π 因子)
spec = QubitSpec(name="Q0", EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
                 T1=10_000, T2=8_000)

print(spec.frequency() / (2 * np.pi))          # 甜点(sweet spot)频率(GHz)
print(spec.anharmonicity() / (2 * np.pi))       # 非谐性 ≈ -0.2 GHz

# 移到最灵敏工作点(返回新对象,不改原 spec)
sweet = QubitSpec.optimal_work_point()           # Φ_0
biased = spec.with_flux(sweet)
print(biased.sensitivity() / (2 * np.pi))        # df/dΦ(GHz/Φ_0)

# 运行时器件对象:实验各层默认接收的类型,持有现成算符与 Hamiltonian
qubit = TransmonQubit(EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
                      T1=10_000, T2=8_000)
H = qubit.get_hamiltonian()                       # 时间无关 H(Qobj)
c_ops = qubit.get_collapse_operators()            # [√γ1·a, √γφ·n]
```

## 物理角色 / 扩展

- `TransmonQubit` / `QubitSpec` 对应真实芯片上的一个物理比特;其参数直接来自器件表征或流片参数(推荐范围见 `sqc.config`:$E_J/E_C\in[40,80]$、$f_{01}\in[4,8]\,\mathrm{GHz}$)。
- `Resonator` 对应读出/耦合谐振腔;`CoupledSystem` 对应双比特+耦合器+腔的两比特门实验单元;`ChipTopology` 对应整块多比特芯片的布局与连接。
- **要加自定义器件类型**(如 fluxonium、多结比特),继承 `Device` 抽象基类,实现
  `name`、`hilbert_dim()`、`hamiltonian_static()`、`collapse_operators()` 四个成员,
  并保持"纯参数对象"约定(不存实验状态)。完整扩展指南见 {doc}`../extending`。
