# 硬件(hardware)

## 这层提供什么

`hardware` 层建模从 AWG 输出到芯片比特之间的这段模拟链路:信号离开任意波形
发生器后,经同轴线、衰减器、滤波器、bias-tee 传到芯片,一路上被线缆与滤波网络
的传递函数所畸变。这一层把这段畸变刻画为线性时不变(LTI)系统,并提供两个方向
的操作:前向(AWG 波形 → 芯片上真实波形)与逆向/预畸变(想要的芯片波形 → 该发
的 AWG 波形)。

这一层只描述链路本身的传函畸变与读出转导,不涉及比特物理(见 {doc}`devices`)、
脉冲设计(见 {doc}`control`)或时间演化(见 {doc}`simulation`)。

本层的输入输出统一使用 {doc}`control` 层定义的 `Waveform`(语义中立的时域容器)和
`FluxSignal`(其子类,samples 以 $\Phi_0$ 计);两个类型在 control 层定义,hardware
层通过 `import` 引用。信号类型不在本层重新定义的原因是代码依赖方向:hardware 依赖
control,而 control 不依赖 hardware,把类型定义放在 control 可以避免循环导入。

一条贯穿本层的约定:畸变模型是纯 LTI 描述。给定采样步长 `dt`,同一个模型能一致
地给出时域 `apply()`、阶跃响应、冲激响应和频率响应四种视图,四者互相自洽,直流
增益恒为 1(即静态偏置不被畸变改变,只有暂态被扭曲)。

## 类总览

| 类 | 角色 | 说明 |
|---|---|---|
| `DistortionModel` | 抽象基类 | 所有控制线畸变的公共契约(LTI),本层主扩展点 |
| `SingleExponentialDistortion` | 畸变子类 | 单指数拖尾,最常见的 Z 线畸变模型 |
| `MultiExponentialDistortion` | 畸变子类 | K 个单指数拖尾之和 |
| `FIRDistortion` | 畸变子类 | 有限冲激响应滤波器 |
| `IIRDistortion` | 畸变子类 | 无限冲激响应滤波器 |
| `CustomTransferDistortion` | 畸变子类 | 用户在频率栅格上直接给定 $H(\omega)$ |
| `CascadeDistortion` | 畸变子类 | 多级畸变串联(Rol 2020 §IV) |
| `ControlLine` | 控制线 | 一条物理控制线(xy/z/readout),可挂一个畸变模型 |
| `TransferMatrix` | 传递矩阵 | 多条 Z 线的频域传递关系 $\Phi_j(\omega)=\sum_i H_{ji}(\omega)V_i(\omega)$ |
| `ReadoutModel` | 抽象基类 | 读出模型的公共契约,读出侧扩展点 |
| `IdealProjectiveReadout` | 读出子类 | 对 $|1\rangle$ 的投影测量,返回 $p_e$ |
| `IQReadoutModel` | 读出子类 | 两条 $\pi/2$ 相位差 Ramsey 序列做 I/Q 解调 |

## 物理模型与单位

**单位约定**:时间以 ns 计,采样步长 `dt` 来自 `sqc.config.CONFIG.awg.dt`;
角频率 $\omega$ 以 rad/ns 计;磁通以 $\Phi_0$ 为单位;畸变幅度 `amplitude` 为无量纲
分数(典型 $0.001$–$0.1$),拖尾时间常数 `tau` 以 ns 计(典型 $10$–$1000$ ns)。

控制线畸变的物理来源(参照 Gao 2021 §III.D 控制线传递函数、§V.E 用 Cryoscope
标定畸变拖尾):AWG 发出的理想阶跃,经线缆与低通网络后不会瞬间到位,而是带一条
缓慢弛豫的指数拖尾。单指数模型是最常见的刻画:

$$s(t) = 1 - A\,e^{-t/\tau}, \qquad H(s) = (1-A) + \frac{A}{1 + s\tau}$$

其中 $A$=`amplitude`、$\tau$=`tau`。DC 增益 $H(0)=1$。所有子类都提供四种自洽视图:

- `apply(waveform, dt)`：时域,把畸变作用到均匀采样波形(经双线性变换离散化)。
- `step_response(t)`：阶跃响应 $s(t)$。
- `impulse_response(t)`：冲激响应 $h(t)$。
- `frequency_response(omega)`：复频响 $H(\omega)$。

**`smooth` 标志(所有畸变子类通用)**:

- `smooth=False`(默认)：保留直通(delta 函数)路径,输入的间断得以保留,频谱较宽;
  阶跃输入在 $t=0$ 处有跳变。这是向后兼容的行为。
- `smooth=True`：纯低通,无直通路径;理想阶跃产生平滑爬升(如 $1-e^{-t/\tau}$),
  无跳变。此时 `amplitude` 被忽略(拖尾权重固定为 1 以保证 DC 增益 = 1)。

## DistortionModel：畸变抽象基类

所有控制线畸变的公共契约,是本层的主扩展点。低层接口用 numpy 数组,便捷方法
可直接处理 `control` 层对象。子类必须实现 4 个抽象方法(即上面的四种视图):

- `apply(waveform: np.ndarray, dt: float) -> np.ndarray`：时域畸变。
- `step_response(t) -> np.ndarray` / `impulse_response(t) -> np.ndarray`。
- `frequency_response(omega) -> np.ndarray`：复频响。

基类另提供两个便捷包装(子类无需重写):

- `apply_to_waveform(wf) -> Waveform`：接收一个 `Waveform`,返回畸变后的新 `Waveform`。
- `apply_to_signal(signal) -> FluxSignal`：接收一个 `FluxSignal`,返回 `type=8`
  (用户自定义原始采样)的新 `FluxSignal`,这是全栈中重建/畸变信号的规范通路。

要加自定义畸变类型,继承 `DistortionModel` 实现上述 4 个抽象方法即可,详见
{doc}`../extending`。

## 内建畸变子类

### SingleExponentialDistortion：单指数拖尾

最常见的 Z 线畸变模型,即上文 $s(t)=1-A e^{-t/\tau}$。

**构造**

`SingleExponentialDistortion(amplitude=0.01, tau=100.0, smooth=False)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `amplitude` | float | 拖尾幅度(无量纲分数) | `0.01` |
| `tau` | float | 时间常数(ns) | `100.0` |
| `smooth` | bool | 纯低通模式(无直通路径) | `False` |

**方法**(除继承的 4 种 LTI 视图外)

- `design_inverse(dt, formula="bilinear") -> IIRDistortion`：设计一个 IIR 逆滤波器,
  与前向畸变级联后抵消指数拖尾。`formula` 可取 `"bilinear"`(双线性变换)
  或 `"rol2020"`(直接 z 域极零点设计,Rol 2020 Eq. S22)。当 $|A|>0.5$ 时会告警
  (Rol 2020 验证范围为 $|A|\le 0.1$)。

**输出**

`apply()` 返回 `np.ndarray`(畸变后采样);`apply_to_waveform()` 返回 `Waveform`;
`apply_to_signal()` 返回 `FluxSignal`;`design_inverse()` 返回 `IIRDistortion`。

### MultiExponentialDistortion：多指数拖尾

K 个单指数拖尾之和,刻画多个时间尺度并存的畸变。

**构造**

`MultiExponentialDistortion(amplitudes, taus, smooth=False)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `amplitudes` | `np.ndarray` | 各分量拖尾幅度(长度 K) | —— |
| `taus` | `np.ndarray` | 各分量时间常数(ns,长度 K) | —— |
| `smooth` | bool | 纯低通模式 | `False` |

两数组需等长。只读属性 `n_components` 返回 K。

**方法**(除继承的 4 种 LTI 视图外)

- `design_inverse(dt, formula="bilinear", n_iir_stages=3, fir_taps=72,
  fir_threshold_ns=30.0) -> CascadeDistortion`：按 Rol 2020 §IV 的迭代分解设计级联逆
  滤波器:$\tau$ 大于 `fir_threshold_ns` 的慢分量各配一个 IIR 级(至多 `n_iir_stages`
  个),再用一个 FIR 级校正快速残差。

**输出**

`apply()` 返回 `np.ndarray`;`design_inverse()` 返回 `CascadeDistortion`。

### FIRDistortion：有限冲激响应

$y[n]=\sum_k b_k\,x[n-k]$。

**构造**

`FIRDistortion(taps, smooth=False, smooth_tau=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `taps` | `np.ndarray` | FIR 系数 $b_k$ | —— |
| `smooth` | bool | 纯低通模式 | `False` |
| `smooth_tau` | float | 平滑低通时间常数(仅 `smooth=True` 时用) | `None` |

**类方法**

- `FIRDistortion.design_from_residual(s_residual, dt, n_taps=72,
  ridge=1e-6) -> FIRDistortion`：用最小二乘(带岭正则)求一组 taps,把残差阶跃响应
  校正回单位阶跃,供 `MultiExponentialDistortion.design_inverse` 的 FIR 残差级使用。

**输出**

`apply()` 返回 `np.ndarray`;`design_from_residual()` 返回 `FIRDistortion` 实例。

### IIRDistortion：无限冲激响应

$y[n]=\sum_k b_k x[n-k]-\sum_{k>0} a_k y[n-k]$,$a_0=1$。

**构造**

`IIRDistortion(b_coeffs, a_coeffs, smooth=False, smooth_tau=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `b_coeffs` | `np.ndarray` | 前馈系数 | —— |
| `a_coeffs` | `np.ndarray` | 反馈系数(若 $a_0\neq 1$ 会自动归一化) | —— |
| `smooth` | bool | 纯低通模式 | `False` |
| `smooth_tau` | float | 平滑低通时间常数 | `None` |

`SingleExponentialDistortion.design_inverse` 的返回类型即此类。

**输出**

`apply()` 返回 `np.ndarray`。

### CustomTransferDistortion：自定义频域传函

用户在频率栅格上直接给定复 $H(\omega)$,内部走 FFT → 乘 $H(\omega)$ → IFFT。

**构造**

`CustomTransferDistortion(omega_grid, H_grid, smooth=False, smooth_tau=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `omega_grid` | `np.ndarray` | 频率栅格(rad/ns) | —— |
| `H_grid` | `np.ndarray` | 对应复频响 | —— |
| `smooth` | bool | 纯低通模式 | `False` |
| `smooth_tau` | float | 平滑低通时间常数 | `None` |

栅格外的频点按端点值外推。

**输出**

`apply()` 返回 `np.ndarray`。

### CascadeDistortion：多级串联

把若干畸变级顺序作用:$\text{out}=\text{stage}_N(\cdots\text{stage}_1(\text{in}))$;
复合频响是各级频响之积;空 `stages` 即恒等系统。

**构造**

`CascadeDistortion(stages)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `stages` | `list[DistortionModel]` | 畸变模型列表(按作用顺序) | —— |

上文各 `design_inverse` 方法返回的即是此类。

**输出**

`apply()` 返回 `np.ndarray`;`frequency_response()` 返回各级频响之积。

## ControlLine：单条物理控制线

建模一条物理控制线(xy / z / readout),可选挂一个畸变模型。
它把畸变模型包装成"信号源 → 目标"的链路对象,并额外持有链路的物理参数。

**构造**

`ControlLine(name, kind, source, target, transfer_function=None, impedance=50.0, attenuation_db=20.0, delay=0.0, ...)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `name` | str | 控制线标识符(如 `"Z0"`) | —— |
| `kind` | str | 线类型 | `"xy"` / `"z"` / `"readout"` |
| `source` | str | 信号源标识(如 `"AWG0:CH1"`) | —— |
| `target` | str | 目标器件(如 `"Q0"`) | —— |
| `transfer_function` | `DistortionModel` 或 None | 挂载的畸变模型 | `None` |
| `impedance` | float | 特征阻抗(Ω) | `50.0` |
| `attenuation_db` | float | 室温到芯片总衰减(dB) | `20.0` |
| `delay` | float | 传播延迟(ns) | `0.0` |
| `metadata` | dict | 附加元数据 | `{}` |

**方法**

- `apply(awg_waveform: Waveform) -> Waveform`：**前向**:AWG 波形 → 芯片上波形。
  无 `transfer_function` 时原样返回(若 `delay>0` 则做时间平移)。
- `predistort(target_waveform, designer) -> Waveform`：**逆向**:想要的芯片波形 →
  该发的 AWG 波形,委托给 `calibration` 层的 `PredistortionDesigner`(见 {doc}`calibration`)。

**输出**

`apply()` 返回 `Waveform`(片上波形);`predistort()` 返回 `Waveform`(预畸变后的 AWG 波形)。

## TransferMatrix：多线传递矩阵

在频域刻画多条 Z 线之间的传递关系。设某目标比特 $j$ 上的磁通由所有
源线电压共同决定:

$$\Phi_j(\omega) = \sum_i H_{ji}(\omega)\,V_i(\omega)$$

对角元 $H_{ii}$ 是每条线对自身目标的响应;非对角元 $H_{ji}\,(j\neq i)$ 描述一条线
的电压串扰到另一比特磁通上的部分。

```{note}
本页只把 `TransferMatrix` 当作**前向**器件(把 AWG 电压映射为芯片磁通)介绍。
如何从测量中**反解**串扰矩阵,不在 v1 公开范围内。
```

**构造**

`TransferMatrix(elements, frequency_axis, time_axis=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `elements` | dict | `(target_name, source_name) -> H_ji(ω)` 复数数组字典 | —— |
| `frequency_axis` | `np.ndarray` | 频率栅格(rad/ns) | —— |
| `time_axis` | `np.ndarray` | 时域轴(可选) | `None` |

**只读属性**：`sources` / `targets`(排序后的源/目标名列表)。

**方法**

- `H_ji(target, source) -> np.ndarray`：取单个传函元素。
- `diagonal() -> dict`：各目标的自响应 $H_{ii}$。
- `off_diagonal() -> dict`：串扰元素 $H_{ji}\,(j\neq i)$。
- `apply(source_voltages: dict[str, Waveform]) -> dict[str, FluxSignal]`：**前向**:
  对每个目标 $j$ 做 $\Phi_j=\text{IFFT}(\sum_i H_{ji}\cdot\text{FFT}(V_i))$,不等长
  波形按最长者补零对齐。

**类方法**

- `TransferMatrix.from_dc_matrix(dc_matrix, source_names, target_names, n_freq=256)
  -> TransferMatrix`：从一个直流串扰矩阵构造**频率平坦**的传递矩阵($H_{ji}(\omega)$
  对所有 $\omega$ 取常数 `dc_matrix[j, i]`)。

**输出**

`apply()` 返回 `dict[str, FluxSignal]`(目标名 → 芯片磁通);`H_ji()` 返回 `np.ndarray`。

## 读出模型(readout)

读出侧的扩展点与两个内建实现。测量结果统一以 `dict[str, float]` 返回。

### ReadoutModel：读出抽象基类

读出模型的公共契约,是**读出侧扩展点**。子类实现单一抽象方法
`measure(*args, **kwargs) -> dict[str, float]`。

### IdealProjectiveReadout：理想投影测量

对 $|1\rangle$ 的投影测量。

**构造**

`IdealProjectiveReadout()`

**方法**

- `measure(state, qubit=None) -> dict[str, float]`：返回 `{"p_e": float}`,
  即激发态布居 $p_e=\langle 1|\rho|1\rangle$。

**输出**

`dict` 含单键 `"p_e"`。

### IQReadoutModel：I/Q 读出

跑两条 $\pi/2$ 相位差的 Ramsey 序列做 I/Q 解调(取代旧 `src/protocal.py:IQ_readout`)。

**构造**

`IQReadoutModel(tau=20.0, t_rabi=None, omega_d=None)`

**字段**

| 字段 | 类型 | 含义 | 默认 |
|---|---|---|---|
| `tau` | float | 自由进动时间(ns) | `20.0` |
| `t_rabi` | `np.ndarray` | Rabi 时间轴(ns) | `CONFIG.pulse.t_rabi` |
| `omega_d` | float | 驱动频率(rad·GHz) | `None`(用比特频率) |

**方法**

- `measure(qubit, **extra) -> dict[str, float]`：返回 `{"p_e_I": float, "p_e_Q": float}`;
  要求 `qubit` 已先调用 `qubit_in_mag(...)` 把 `H_list` 备好。

**输出**

`dict` 含 `"p_e_I"` 和 `"p_e_Q"` 两键。

## 最小用例

```python
import numpy as np
from sqc.config import CONFIG
from sqc.control.waveform import Waveform
from sqc.hardware import SingleExponentialDistortion, ControlLine

# 采样步长统一从 CONFIG 派生;时间轴用 arange(不用 linspace)
dt = CONFIG.awg.dt
t = np.arange(0.0, 200.0, dt)

# 一个理想阶跃波形(AWG 想发的)
awg = Waveform(t_list=t, samples=np.where(t >= 20.0, 1.0, 0.0))

# 单指数拖尾畸变:A=2%,tau=50 ns
dist = SingleExponentialDistortion(amplitude=0.02, tau=50.0)

# 挂到一条 Z 线上,前向得到芯片上真实波形
line = ControlLine(name="Z0", kind="z", source="AWG0:CH1",
                   target="Q0", transfer_function=dist)
on_chip = line.apply(awg)          # 阶跃边沿被指数拖尾展宽

# 设计逆滤波器,级联后抵消拖尾(预畸变思路)
inv = dist.design_inverse(dt, formula="rol2020")
corrected = inv.apply_to_waveform(on_chip)   # 拖尾被补偿,恢复为理想阶跃
```

## 物理角色 / 扩展

- `DistortionModel` 家族对应真实控制线的传递函数:线缆、衰减器、低通滤波器、
  bias-tee 共同造成的 LTI 畸变。`SingleExponentialDistortion` 对应 Cryoscope 标定中
  最常见的慢弛豫拖尾(Gao 2021 §V.E)。
- `ControlLine` 对应芯片上一条物理布线(xy 驱动线 / z 磁通线 / 读出线);
  `TransferMatrix` 对应多条 Z 线并存时的频域响应网络。
- `ReadoutModel` 家族对应读出链路:`IdealProjectiveReadout` 是无噪投影测量的理想
  上限,`IQReadoutModel` 对应真实的 I/Q 解调读出。
- 要加自定义畸变类型(如带纹波的传函、非最小相位滤波器),继承 `DistortionModel`
  实现 `apply` / `step_response` / `impulse_response` / `frequency_response` 四个方法;
  要加自定义读出,继承 `ReadoutModel` 实现 `measure`。完整扩展指南见 {doc}`../extending`。
