# 标定(calibration)

## 这层提供什么

`calibration` 层是全栈的**标定与整定层**——把"比特/控制线**当前**是什么状态"测出来,
并在需要时把它**调到**目标状态。它服务两条产品主线:频率标定(测/整 $f_{01}(\Phi)$)
与预畸变(测控制线传函、设计补偿滤波器)。产出统一封装成
{py:class}`~sqc.calibration.CalibrationTable`,供下游查表反演或整定使用。

与 {doc}`reconstruction` 层的**纯函数**约束不同,本层的标定类**会主动跑仿真**去测量:
每个 `Calibration` 子类实现 `calibrate()`,内部驱动比特、跑 `mesolve`、拟合,最终返回一张
`CalibrationTable`。

```{note}
{doc}`reconstruction` 层还导出两个概念属标定的类 `CryoscopeCalibration` /
`DelayRamseyCalibration`(与其重建器配套导出)。它们同样继承本层的
{py:class}`~sqc.calibration.Calibration`,产出 `CalibrationTable` 供 `inversion="calibration"`
查表反演。分述见 {doc}`reconstruction`。
```

## 类总览

按功能分四组:

**扩展点 + 结果容器**

| 类 | 角色 |
|---|---|
| `Calibration` | 抽象基类,所有标定工作流的公共契约,本层扩展点 |
| `CalibrationTable` | 标定结果容器,带 `evaluate`(插值)/ `inverse`(反插值)|

**频率标定 / 测量 / 整定**

| 类 | 角色 |
|---|---|
| `FluxResponseCalibration` | 扫磁通建 $f(\Phi)$ 查表(Ramsey 逐点测频)|
| `FrequencyMeasurement` | 单点 $f_{01}$ 测量(Ramsey,或进阶的瞬态法)|
| `SinglePointFrequencyCalibration` | 单点频率**整定**:闭环反馈把 $f_q(V)$ 调到目标 |

**波形 / 控制线标定**

| 类 | 角色 |
|---|---|
| `WaveformCalibration` | 传函测量 + 预畸变设计的统一入口 |
| `PredistortionDesigner` | 独立的逆滤波器设计器(FIR/IIR/频域反演)|

**调度**

| 类 | 角色 |
|---|---|
| `CalibrationScheduler` | 标定任务调度器:注册表 + 依赖图,按序执行 |

## Calibration —— 标定抽象基类

所有标定工作流的公共契约,本层**扩展点**。只规定一个抽象方法:

- `calibrate() -> CalibrationTable` —— 跑标定工作流,返回结果表。

要加自定义标定,继承 `Calibration` 实现 `calibrate()`,详见 {doc}`../extending`。

## CalibrationTable —— 标定结果容器

统一的标定结果载体(`@dataclass`)。核心字段:`name`、`qubit_name`、`kind`(如
`"f_phi"`、`"f01"`、`"transfer_function"`、`"predistortion"`)、`inputs`/`outputs`
(自变量/因变量,如磁通↔频率)、`fit_params`(拟合/迭代细节)、`metadata`。

两个查表方法(基于 SciPy 三次样条,点数不足时自动降到二次/线性,并外插):

- `evaluate(x) -> np.ndarray` —— 在查询点 `x` 上插值 `outputs`(如给磁通取频率)。
- `inverse(y) -> np.ndarray` —— 反插值:找使 `outputs≈y` 的 `inputs`(需 `outputs` 单调,
  内部取单调段构造反函数)。

## FluxResponseCalibration —— f(Φ) 曲线标定

扫一串直流磁通点,在每点用 Ramsey 测频,建 $f(\Phi)$ 查表。`@dataclass`,字段:
`qubit`、`method`(`"ramsey"`)、`h_list`(磁通扫描点,默认 `linspace(-0.03, 0.03, 51)`)、
`tau`(每点自由进动时间)、`t_rabi`(Rabi 脉冲时间轴)。`calibrate()` 返回
`kind="f_phi"` 的表,`inputs=磁通`、`outputs=角频率`。

```{note}
`method="transient"`(未知瞬态信号 → $\Delta\omega(\Phi)$ 多项式拟合)属未来功能,
当前抛 `NotImplementedError`,不进 v1 公开面。
```

## FrequencyMeasurement —— 单点 f01 测量

在单个磁通工作点测 $f_{01}$(只读、不整定)。`@dataclass`,字段:`qubit`、`flux`
(测量磁通点,默认 0 即甜点)、`method`、`tau_list`、`t_rabi`、`t_global`、
`f_artificial`。核心方法 `measure(flux=None) -> float` 返回有符号角频率;
`calibrate()` 把单点测量打包成 `kind="f01"` 的表。

- `method="ramsey"`(默认)—— Ramsey τ 扫描 + FFT 取峰。单扫模式(`f_artificial=0.1`,
  假设 $|\Delta|<0.1$ GHz)快;`f_artificial=None` 走双扫模式,对任意 $|\Delta|$ 稳健、
  给符号,代价 2×。
- `method="transient"`(**进阶**)—— τ=0 正交 Ramsey(R_y–R_x 与 R_y–R_{-x})差分读出
  + 控制核灵敏度 $G=\int k_1\,\mathrm{d}t$ 直接反出 $\Delta\omega$。比 τ 扫描省,但依赖弱信号
  线性近似,适合 $|\Delta\omega|$ 近零。`order>=3` 时加三次 Newton 修正,三次系数 $G_3$
  由 `g3_source` 选:`"fit"`(奇多项式拟合 $p_\mathrm{diff}(\Delta)$,自适应扫描范围)或
  `"kernel_full"`(完整非对角核三重积分 $\iiint k_3$)。

```{note}
瞬态法是单点**频率测量**,与 §8 边界排除的瞬态**频率标定**(`FluxResponseCalibration`
的 `method="transient"`,建 $f(\Phi)$ 曲线)不是一回事——前者已实现且公开,后者未实现。
瞬态**波形重建**则是 v1 核心特性,见 {doc}`reconstruction`。
```

## SinglePointFrequencyCalibration —— 单点频率整定

把 $f_q(V)$ 闭环反馈整定到目标频率 $f_\mathrm{target}$(Vepsalainen 2022)。`@dataclass`,
字段:`qubit`、`method`(`"closed_loop"`)、`f_target`、`V_a`/`V_b`(括号边界,
来自前序 `FluxResponseCalibration`)、`epsilon_f`(收敛容差,默认 $10^{-4}$ GHz·2π)、
`max_iter`、`measure_method`(`"ramsey"` 或 `"transient"`)、`bracket_tightening`
(regula falsi 缩窄括号,默认开)、`step_method` 三选一:

- `"secant"`(默认)—— 割线法,超线性收敛,典型 1–3 次迭代。
- `"bisection"` —— 二分法,$O(\log_2)$ 收敛,每步框宽折半,便于可视化。自动处理
  $f(\Phi)$ 偶函数越过甜点的情形(自动分裂括号)。
- `"gradient"` —— 阻尼割线(数值梯度 Newton 步),**不要求预先括号** $V_a$/$V_b$,
  只需 `V_seed`;阻尼因子 `damping`(默认 0.8)抑制超调;`best_V` 追踪历史最优点。

`calibrate()` 返回 `kind="f01"` 的表,`fit_params["history"]` 含完整迭代轨迹。

## WaveformCalibration —— 波形/控制线标定

预畸变主线的统一入口(`@dataclass`),`method` 选两条路:

- `"transfer_function"` —— 测控制线阶跃响应、拟合成 `DistortionModel`。字段
  `distortion`(解析路径:直接取模型阶跃响应),或 `measurement_protocol`
  (`"cryoscope"`/`"delay_ramsey"`/`"transient"`/`"pi_pulse"`,配 `qubit`+`control_line`,
  走真实量子仿真测量);`fit_type`(`"multi_exp"`/`"single_exp"`/`"fir"`/`"iir"`)。
  `to_distortion_model()` 直接把结果转成 `DistortionModel`。
- `"predistortion"` —— 给定已测传函 `transfer_model`,设计补偿滤波器(内部委托
  `PredistortionDesigner`);字段 `predistortion_method`、`n_taps`、`regularization`。

`calibrate()` 按 `method` 返回 `kind="transfer_function"` 或 `"predistortion"` 的表。

## PredistortionDesigner —— 逆滤波器设计器

独立的预畸变滤波器设计器(`@dataclass`),可单用或被 `WaveformCalibration` 调用。
字段:`method`(`"auto"`/`"fir_inverse"`/`"iir_inverse"`/`"frequency_inverse"`)、
`n_taps`、`regularization`。核心方法:

- `design(transfer_model, dt) -> DistortionModel` —— 给定传函模型返回其逆模型。
  `"auto"` 对指数型模型走解析 IIR 逆,其余走频域反演。
- `predistort(target, transfer=..., 或 inverse_model=...) -> Waveform` —— 直接把预畸变
  施加到目标波形上。
- `check_pole_stability(b, a) -> bool` —— 检查 IIR 滤波器极点是否都在单位圆内(稳定)。

## CalibrationScheduler —— 标定任务调度器

标定"控制室"(`@dataclass`):维护带依赖的命名标定任务注册表,按序执行。
`register(name, cal_class, depends_on=...)` 注册任务;`register_defaults()` 装入标准注册表
(如 `frequency_closed_loop` 依赖 `flux_response_ramsey`、`waveform_predistortion` 依赖
`waveform_transfer_function`)。执行:`run(name, **kw)` 按名跑(自动注入依赖结果)、
`run_next()` 跑下一个就绪任务;`get_result(name)`、`status()` 查询。

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

- 本层对应真实控制室里的**标定程序**:投片后第一步测 $f(\Phi)$,然后整定工作频率,
  最后测控制线传函并烧入补偿滤波器。三步分别对应 `FluxResponseCalibration`→
  `SinglePointFrequencyCalibration`→`WaveformCalibration`。
- `CalibrationTable` 的 `evaluate`/`inverse` 是 {doc}`reconstruction` 层
  `inversion="calibration"` 路径的基础:重建器把标定表传入,在反演时调 `inverse(φ)`
  把测到的相位翻回磁通幅度。
- `CalibrationScheduler` 的依赖注入机制自动把 `FluxResponseCalibration` 的输出
  (磁通范围)穿给 `SinglePointFrequencyCalibration`,省去手动传参。
- **要加自定义标定工作流**,继承 `Calibration` 抽象基类实现 `calibrate() -> CalibrationTable`,
  然后向 `CalibrationScheduler` 注册并声明依赖。完整扩展指南见 {doc}`../extending`。
