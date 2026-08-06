# 预畸变

## 概述

预畸变（predistortion）是本平台的第三条产品主线：任意波形发生器 (AWG) 输出的
电压波形经过控制线（线缆、滤波器、bias-tee 等）传输到达超导芯片时会产生线性
畸变——理想阶跃信号出现指数拖尾。预畸变的核心思想是**在 AWG 端预先施加控制线
传递函数的逆变换**，使波形经过控制线的正向畸变后恰好还原为目标波形。

在物理上，控制线建模为一个线性时不变 (LTI) 系统，由传递函数 $H(\omega)$ 或
等价的阶跃响应 $s(t)$ 完整描述。片上实际磁通波形是 AWG 输出与控制线冲激响应
的卷积：

$$\Phi_\text{chip}(t) = (h * V_\text{AWG})(t)
\quad\Longleftrightarrow\quad
\Phi_\text{chip}(\omega) = H(\omega) \, V_\text{AWG}(\omega).$$

预畸变在 AWG 波形前端级联一个**逆滤波器** $H^{-1}(\omega)$，使得总传递函数
恢复为恒等映射：

$$V_\text{AWG} = H^{-1} * \Phi_\text{target}
\;\Longrightarrow\;
\Phi_\text{chip} = H * H^{-1} * \Phi_\text{target} = \Phi_\text{target}.$$

## 管道架构

预畸变管线的执行顺序是：定义目标波形 → 建模控制线畸变 → 测量传递函数 →
设计逆滤波器 → 验证补偿效果。以下按此逻辑顺序展开各层的职责。

### 第一步：定义目标波形 — 控制层 (`Waveform`)

一切预畸变的起点是**目标波形**——你希望在芯片上出现的磁通信号。
{py:class}`~sqc.control.waveform.Waveform` 是通用时域信号的数据结构，持有时
间轴 `t_list` 和采样值 `samples`，提供 `value_at(t)`（采样保持取值）、
`truncate(t_start, t_end)`（时间段截断）、`samples_on(t_global)`（投影到全局
时间轴）、`copy()`（深拷贝）等语义无关操作。

```{note}
预畸变管线中目标波形和 AWG 输出均以 `Waveform` 传递——此时信号尚未到达比特，
无需赋以磁通语义。物理磁通信号对应的专用子类是
{py:class}`~sqc.control.flux_signal.FluxSignal`，仅在量子仿真测量路径中涉及。
```

### 第二步：建模控制线畸变 — 硬件层 (`DistortionModel` + `ControlLine`)

目标波形从 AWG 输出后，需经过一条物理控制线才能到达芯片。这一传输过程的畸变
效应由两层抽象建模：畸变模型定义数学形式，控制线封装物理语境。

#### 畸变模型

{py:class}`~sqc.hardware.distortion.DistortionModel` 是所有畸变模型的抽象基类，
定义统一的 LTI 接口：

- `apply(waveform, dt)`：将畸变施加到离散采样波形。
- `step_response(t)` / `impulse_response(t)` / `frequency_response(omega)`：
  阶跃响应 $s(t)$、冲激响应 $h(t)$、复频率响应 $H(\omega)$。
- `apply_to_waveform(wf)`：对 `Waveform` 对象的便捷封装。

六种具体模型覆盖从简单到复杂的各类场景：

| 类 | 数学形式 | 典型应用 |
|---|---|---|
| `SingleExponentialDistortion` | $s(t) = 1 - A e^{-t/\tau}$ | 最常见的磁通线拖尾 |
| `MultiExponentialDistortion` | $s(t) = 1 - \sum_k A_k e^{-t/\tau_k}$ | 多级滤波、多路径反射 |
| `IIRDistortion` | $y[n] = \sum b_k x[n-k] - \sum_{k>0} a_k y[n-k]$ | 通用 IIR 滤波、逆滤波器 |
| `FIRDistortion` | $y[n] = \sum b_k x[n-k]$ | 通用 FIR 滤波、残差补偿 |
| `CustomTransferDistortion` | 用户定义 $H(\omega)$ 频率网格 | 任意频域测量 LTI 系统 |
| `CascadeDistortion` | $H(\omega) = \prod_k H_k(\omega)$ | 多级串联补偿 (IIR + FIR) |

所有模型支持 `smooth` 参数：`smooth=False`（默认）保留直流直通路径，阶跃输入
在 $t=0$ 处有跳变；`smooth=True` 为纯低通行为，无直通路径。

#### 控制线

{py:class}`~sqc.hardware.control_line.ControlLine` 将畸变模型置于物理语境中，
建模一条具体的控制线（xy/z/readout）。它持有 `transfer_function`（一个
`DistortionModel` 实例）以及名称、源/目标、阻抗、衰减等辅助属性。两个关键
方法对应管线的两个方向：

- `apply(awg_waveform)` — **正向**：AWG 波形 $\to$ 控制线 $\to$ 片上波形。
- `predistort(target, designer)` — **逆向**：目标片上波形 $\to$ 逆滤波 $\to$ 所需 AWG 波形。

### 第三步：测量传递函数 — 标定层 (`WaveformCalibration`)

有了畸变模型和控制线，下一步是**测量**控制线的实际传递函数——即向系统输入已
知信号、观测输出、拟合模型参数。这是标定层的 `WaveformCalibration` 的职责。

{py:class}`~sqc.calibration.waveform.WaveformCalibration` 提供统一的测量与拟
合接口，两条路径：

- **解析路径**（`distortion` + `fit_type`）：直接对已知 `DistortionModel` 的
  解析阶跃响应函数做拟合。`fit_type` 可选 `"single_exp"`（单指数）、
  `"multi_exp"`（多指数和）、`"fir"`（有限冲激响应）、`"iir"`（无限冲激响应）。
- **量子仿真路径**（`measurement_protocol` + `qubit` + `control_line`）：通过
  量子比特协议测量阶跃响应（详见下文）。

拟合完成后，`to_distortion_model()` 将拟合参数导出为对应的 `DistortionModel`
子类实例——这是**正向**畸变模型，代表控制线对信号的畸变作用。

### 第四步：设计逆滤波器 — 标定层 (`PredistortionDesigner`)

测得正向传递函数后，下一步是求其逆。{py:class}`~sqc.calibration.waveform.PredistortionDesigner`
是独立的逆滤波器设计器：

| `method` | 算法 | 适用条件 |
|---|---|---|
| `"auto"` | 根据模型类型自动选择 | 默认，推荐 |
| `"iir_inverse"` | 解析 IIR 逆（双线性变换或 z 域零极点） | 指数型畸变，有闭式解 |
| `"fir_inverse"` | 频域 Wiener 逆 + FIR 截断 | 一般 LTI 畸变 |
| `"frequency_inverse"` | $H^{-1} = \bar{H}/(|H|^2 + \varepsilon^2)$ | 任意频域传递函数 |

对于最常见的单指数畸变 $s(t) = 1 - A e^{-t/\tau}$，`"iir_inverse"` 通过双线
性变换解析得到 IIR 逆滤波器系数；`SingleExponentialDistortion` 自带
`design_inverse()` 方法封装此逻辑。对于多指数畸变，`design_inverse()` 按时间
常数降序排列，对慢分量各分配 IIR 逆级，再用 FIR 补偿快分量残差 (Rol 2020,
§IV)。

`regularization` 参数控制逆的激进程度：值越大补偿越保守、越稳定；值越小
$|H|$ 近零处的噪声放大风险越大。

```{note}
`to_distortion_model()` 返回的是拟合出的**正向**畸变模型；逆滤波器需再经
`PredistortionDesigner.design()` 得到。两者不可混淆。
```

### 第五步：端到端编排 — 工作流层 (`PredistortionValidationWorkflow`)

将上述四步串联成一条完整验证管道的是
{py:class}`~sqc.workflows.PredistortionValidationWorkflow`。它不引入新物理，
仅编排已有各层：

1. 用真值畸变构造 `ControlLine`。
2. **正向未补偿**：目标 $\to$ 控制线 $\to$ 片上波形（含拖尾）。
3. 通过 `WaveformCalibration` 测量传递函数、拟合正向畸变模型。
4. 通过 `PredistortionDesigner` 设计逆滤波器。
5. **施加预畸变**：目标 $\to$ 逆滤波 $\to$ AWG 波形。
6. **正向验证**：预畸变 AWG 波形 $\to$ 控制线 $\to$ 片上波形，与目标比较，
   计算 RMSE、改善倍数与整定时间。

仅需传入 `target_waveform` 和 `true_distortion`，调用 `run()` 即可获得完整的
补偿前后对比。`measurement_protocol` 字段用于切换标定的解析/量子仿真路径。

### 量子仿真测量路径

以上管道默认走解析路径——直接对已知畸变模型的阶跃响应公式做拟合。若要仿真
**实验室中通过 qubit 测量控制线传函**的真实流程，需在
`WaveformCalibration` 或 `PredistortionValidationWorkflow` 上传入
`measurement_protocol` 和 `qubit`。此时系统转而调用以下三层：

**设备层** — {py:class}`~sqc.devices.transmon.TransmonQubit` 提供量子比特
的物理模型。不同协议对偏置磁通要求不同：Cryoscope 在甜点 ($\Phi=0$) 工作
最佳；延迟 Ramsey 与瞬态协议在 $\kappa = d\omega/d\Phi$ 最大处灵敏度最高。

**实验层** — 对应的实验类运行量子仿真，向控制线注入测试信号（默认阶跃），
经 `mesolve()` 获取比特响应。**重建层** — 对应的重建类将比特响应反演回片上
磁通波形，归一化即得阶跃响应。

协议驱动的核心是内部的 `_ProtocolDrivenMeasurement`：构造实验 $\to$ 运行
`mesolve` $\to$ 重建片上波形 $\to$ 归一化。支持的四种协议：

| `measurement_protocol` | 实验类 | 重建类 | 原理 |
|---|---|---|---|
| `"cryoscope"` | `CryoscopeExperiment` | `CryoscopeReconstruction` | 方波磁通 + Ramsey 相位解调 |
| `"delay_ramsey"` | `DelayRamseyExperiment` | `DelayRamseyReconstruction` | 延迟 Ramsey 测拖尾 |
| `"transient"` | `TransientSensingExperiment` | `TransientReconstruction` | 滑动 Ramsey + Wiener 反卷积 |
| `"pi_pulse"` | `PiPulseCompensationExperiment` | `PiPulseCompReconstruction` | $\pi$ 脉冲补偿反推拖尾 |

## 使用方式

### 端到端验证（推荐）

```python
import numpy as np
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.workflows import PredistortionValidationWorkflow

# 目标片上波形：0.05 Φ₀ 的平台
dt = 1.0  # ns
t = np.arange(0, 500, dt)
target = Waveform(t_list=t, samples=np.ones_like(t) * 0.05)

# 真值畸变：单指数拖尾，幅度 0.04，时间常数 200 ns
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)

# 一键验证
wf = PredistortionValidationWorkflow(
    target_waveform=target,
    true_distortion=distortion,
)
result = wf.run()

m = result["metrics"]
print(f"补偿前 RMSE = {m['rmse_uncorrected']:.3e}")
print(f"补偿后 RMSE = {m['rmse_corrected']:.3e}")
print(f"改善倍数    = {m['improvement_factor']:.1f}×")
print(f"补偿前整定时间 = {m['settling_uncorrected_ns']:.2f} ns")
print(f"补偿后整定时间 = {m['settling_corrected_ns']:.2f} ns")
```

### 分步执行管道

以下按管道的逻辑顺序逐层操作：测传函 → 设计逆 → 施加预畸变。

```python
from sqc.calibration.waveform import WaveformCalibration, PredistortionDesigner
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion
import numpy as np

# 定义目标
dt = 1.0
t = np.arange(0, 500, dt)
target = Waveform(t_list=t, samples=np.ones_like(t) * 0.05)

# 已知畸变
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)

# 1. 测量传递函数并拟合畸变模型
wf_cal = WaveformCalibration(
    distortion=distortion, fit_type="single_exp",
)
fwd_model = wf_cal.to_distortion_model()  # 正向畸变模型

# 2. 设计逆滤波器
designer = PredistortionDesigner(method="auto", regularization=1e-4)
inverse = designer.design(fwd_model, dt=dt)

# 3. 施加预畸变：目标波形 → AWG 波形
awg_waveform = inverse.apply_to_waveform(target)

print(f"目标波形点数: {target.n_points}")
print(f"AWG 波形点数: {awg_waveform.n_points}")
print(f"逆模型类型: {type(inverse).__name__}")
```

### 量子仿真测量路径

```python
from sqc.devices import TransmonQubit
from sqc.hardware.control_line import ControlLine
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.calibration.waveform import WaveformCalibration

# 比特与畸变
qubit = TransmonQubit(
    EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
    T1=10_000, T2=8_000, flux=0.0,  # 甜点适合 Cryoscope
)
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)

# 控制线
line = ControlLine(
    name="Z0", kind="z", source="AWG0", target="Q0",
    transfer_function=distortion,
)

# 通过 Cryoscope 量子仿真测量阶跃响应
cal = WaveformCalibration(
    qubit=qubit,
    control_line=line,
    measurement_protocol="cryoscope",
    method="transfer_function",
    fit_type="single_exp",
)
fwd_model = cal.to_distortion_model()  # 从量子仿真拟合的正向模型
```

## 结果解读

- `result["target"]`：目标片上波形（`Waveform`）。
- `result["on_chip_uncorrected"]`：无预畸变时正向畸变后的片上波形（含拖尾）。
- `result["on_chip_corrected"]`：经预畸变补偿后的片上波形。
- `result["awg_predistorted"]`：施加了逆滤波器后应送入 AWG 的波形。
- `result["inverse_model"]`：设计出的逆滤波器对象。
- `result["measured_model"]`：标定测量并拟合得到的正向畸变模型。

`metrics` 提供四个关键指标：

- `rmse_uncorrected` / `rmse_corrected`：补偿前后片上波形与目标的均方根误差。
- `improvement_factor`：改善倍数 $=$ `rmse_uncorrected / rmse_corrected`。
- `settling_uncorrected_ns` / `settling_corrected_ns`：补偿前后波形整定至容差
  （默认 $10^{-3}$）所需的时间。

```{important}
解析路径的改善倍数可以非常大。原因在于 `WaveformCalibration` 直接对已知畸变
的解析阶跃响应函数做拟合——相当于已知 $H$ 精确形式后再求逆，补偿后 RMSE 可
低至机器精度（$\sim 10^{-15}$），改善倍数为机器精度倒数级别。这是解析上限，
**不代表真实系统性能**。

要获得有物理意义的改善倍数，应通过量子仿真测量路径
（传 `measurement_protocol="cryoscope"` 及 `qubit`），由 `mesolve` 仿真测
量阶跃响应。此时拟合受限于测量分辨率与噪声，改善倍数落在有限范围（典型
$10$–$10^2\times$）。
```

- `SingleExponentialDistortion` 自带 `design_inverse()` 方法，通过双线性变换
  解析计算 IIR 逆滤波器。这是最常用的解析路径。
- 更复杂的畸变（多指数和、任意频域传递函数）使用 `MultiExponentialDistortion`
  的级联设计或 `PredistortionDesigner` 的频域方法。
- 各畸变模型的完整字段与物理含义见 {doc}`../building_blocks/hardware`。
- 标定与逆滤波器设计的完整 API 见 {doc}`../building_blocks/calibration`。
