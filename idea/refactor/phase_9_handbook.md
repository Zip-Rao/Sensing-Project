# Phase 9 — 级联预失真滤波器 + 协议驱动的阶跃响应测量

> **状态**: 规划阶段，待实施。本 phase 同时落地两块互补的预失真能力：
>   - **Plan A**：基于 Rol 2020 的 IIR 级联 + FIR 残差修正架构，把"单滤波器求逆"升级为"多级联补偿"
>   - **Plan B**：把 `WaveformCalibration` 的阶跃响应测量从解析公式切换到调用真实重建协议（Cryoscope / delay Ramsey / 瞬态磁场 / π 脉冲补偿）
>
> **前置**: P0–P6 已完成；P7 / P8 非硬依赖（详见 §6 协调说明）。
>
> **理论配套**:
> - Rol 等 2020 *Time-domain characterization and correction of on-chip distortion of control pulses in a quantum processor*（[本地副本](../../../scholaraio/data/libraries/papers/Rol-2020-Time-domain-characterization-and-correction-of-on-chip-distortion-of-control-pulses-in-a-quantum-processor/paper.md)）
> - [`_sensing theory.md`](../_sensing%20theory.md) "预失真"节、"阶跃响应测量"节、"设计滤波器"节
>
> **配套文件**:
> - 重构总方案: [_refactor_plan.md](_refactor_plan.md) §5.3 (ControlLine)、§5.4 (Calibration)、§8 (兼容层契约)
> - 当前交接状态: [_handoff_state.md](_handoff_state.md)
> - 上游 phase: [phase_8_handbook.md](phase_8_handbook.md)（filter function 适配）

---

## §0 实施前 API 核对（必读）

下列事实已与代码交叉验证，作为 P9 实施时的硬约束：

| 事项 | 现状 | P9 应对 |
|------|------|---------|
| `Experiment`（[base.py:10](../../sqc/experiments/base.py#L10)） | **是 ABC，不是 dataclass**，无字段 | `control_line` 字段**逐子类**追加；基类只放 `_route_flux()` 方法 |
| `Waveform`/`FluxSignal` 数据字段 | 使用 `.samples`（Waveform 父类的 dataclass 字段） | **禁止**用 `.signal`（那是 `src/signal.py` 旧 API）；`apply_to_signal` 必须用 `.samples` 并保持原类型 |
| `FluxSignal.copy()` | 通过 `__init__` 重建（[flux_signal.py:385](../../sqc/control/flux_signal.py#L385)），用 `_params` 重算 `.samples` | 若 `copy()` 后再覆盖 `.samples`，params 与 samples 不一致；`apply_to_signal` 应走 `FluxSignal(type=8, signal=out)` 构造 |
| `FluxSignal(type=1, t_list=..., amplitude=A)` | `type=1` 是 constant；`amplitude` 是 kwarg | 构造单位阶跃用此模式 |
| `FluxSignal(type=8, t_list=..., signal=arr)` | `type=8` 是 user-defined raw samples（见 [flux_signal.py:212](../../sqc/control/flux_signal.py#L212)） | `apply_to_signal` 返回值用此构造 |
| `DistortionModel.apply_to_waveform`（[distortion.py:84](../../sqc/hardware/distortion.py#L84)） | 返回 `Waveform`（**会丢失** FluxSignal 类型） | 新增 `apply_to_signal` 必须返回 `FluxSignal`，**不要** 复用 `apply_to_waveform` |
| `PredistortionDesigner.method`（[waveform.py:320](../../sqc/calibration/waveform.py#L320)） | Literal: `"auto"\|"fir_inverse"\|"iir_inverse"\|"frequency_inverse"` | 扩展加 `"cascade"`；`_auto_method` 增加多指数判断 |
| `WaveformCalibration.method`（[waveform.py:68](../../sqc/calibration/waveform.py#L68)） | Literal: `"transfer_function"\|"predistortion"`；`calibrate()` 还接受 `"simulation"` 作 alias | **不动**；新增 `measurement` 字段独立于 `method` |
| `WaveformCalibration.predistortion_method` | 与 `PredistortionDesigner.method` 同名但不同对象 | 不要混淆；P9.A.4 改的是后者 |
| 实验类名（精确） | `CryoscopeExperiment`、`DelayRamseyExperiment`、`TransientSensingExperiment`、`PiPulseCompensationExperiment`、`RamseyExperiment`、`EchoExperiment` | 注册表与 import 必须用精确名 |
| 重建类名（精确） | `CryoscopeReconstruction`、`DelayRamseyReconstruction`、`TransientReconstruction`、`PiPulseCompReconstruction` | 注意 `PiPulseComp`（无 "ensation"） |
| `CryoscopeReconstruction` 构造 | `inversion: "calibration"\|"response"`；`response` 模式需 `qubit` | `StepResponseMeasurement` 走 `response` 路径，避免依赖标定步 |
| `ExperimentResult` | `exp.run()` 返回此类型 | 直接 `rec.reconstruct(result)` |
| `rec.reconstruct()` 返回 | `FluxSignal`（使用 `.samples`、`.t_list`） | `StepResponseMeasurement.measure()` 取 `.samples` |

**子代理执行时若发现实际 API 与本表不符，立即 abort 并在 `_handoff_state.md` 标 `DECISION_NEEDED:`，不要自行猜测。**

---

## §1 问题诊断

### 1.1 Plan A — 级联补偿能力缺失

当前 [`sqc/hardware/distortion.py`](../../sqc/hardware/distortion.py) 提供 5 个 `DistortionModel` 子类，但只支持单滤波器或并行叠加：

| 模型 | 数学结构 |
|------|---------|
| `SingleExponentialDistortion` | 一阶低通 |
| `MultiExponentialDistortion` | K 个一阶低通**并行**求和（$H(s) = (1-\sum a_k) + \sum a_k/(1+s\tau_k)$） |
| `FIRDistortion` | 卷积 |
| `IIRDistortion` | 二阶 IIR |
| `CustomTransferDistortion` | 用户给 $H(\omega)$ |

[`ControlLine.transfer_function`](../../sqc/hardware/control_line.py#L55) 只接受**单个** `DistortionModel`。[`PredistortionDesigner`](../../sqc/calibration/waveform.py#L304) 只设计**单个**逆滤波器：

```python
# 当前能力
inverse = designer.design(transfer_model, dt)  # 单滤波器 → 单逆滤波器
```

但 Rol 2020 §IV 描述的工业标准做法是**串联级联**：

```
input → IIR₁ → IIR₂ → ... → IIR_N → FIR → output
```

每个 IIR 各自补偿一个指数分量（$\tau$ > 30 ns），FIR 补偿残余 < 30 ns 的快变畸变。

**后果**：当前框架无法复现论文 Fig. 2(a) 的 0.1% 精度结果——多指数失真只能用频域整体求逆（`_multi_exp_to_iir_inverse` fallback 到 `frequency_inverse`），无法做迭代式逐分量补偿。

### 1.2 Plan B — 阶跃响应测量绕开量子层

[`WaveformCalibration._measure_step_response`](../../sqc/calibration/waveform.py#L121-L124)：

```python
def _measure_step_response(self):
    t = np.linspace(0, self.t_max, self.n_points)
    step = self.distortion.step_response(t)   # ← 直接调数学公式
    return t, step
```

**后果**：
1. 拟合精度只反映 `curve_fit` 性能，**不反映量子测量本身**的精度损失（unwrap 误差、SG 滤波噪声、截断瞬态、SNR 衰减）
2. 已有的 [`CryoscopeExperiment`](../../sqc/experiments/cryoscope.py)、[`DelayRamseyExperiment`](../../sqc/experiments/delay_ramsey.py)、[`TransientSensingExperiment`](../../sqc/experiments/transient.py)、[`PiPulseCompExperiment`](../../sqc/experiments/pi_pulse_comp.py) **没有被接入**预失真链路
3. [`PredistortionValidationWorkflow`](../../sqc/workflows/predistortion_validation.py#L84-L91) 用 `method="simulation"` 走解析路径，未真正闭环

### 1.3 Experiment 不感知 ControlLine

所有 `*Experiment.run()` 都直接 `qubit.qubit_in_mag(flux_signal)`——`ControlLine.transfer_function` 即便设置，也**不会**被注入到量子模拟里。这是 Plan B 的底层障碍。

| 实验类 | 当前行为 | 期望行为 |
|--------|---------|---------|
| `CryoscopeExperiment` | 直接用 `flux_signal` | 经 `control_line` 畸变后再喂给 qubit |
| `DelayRamseyExperiment` | 同上 | 同上 |
| `TransientSensingExperiment` | 同上 | 同上 |
| `PiPulseCompExperiment` | 同上 | 同上 |

---

## §2 目标状态

### 2.1 Plan A — 级联滤波器架构

#### 2.1.1 核心变更

| 组件 | 旧行为 | 新行为 | **如何回旧** |
|------|--------|--------|-------------|
| `DistortionModel` 族 | 5 个叶子类，单一模型 | 增加 `CascadeDistortion` | 不用即可 |
| 单指数逆设计 | `_single_exp_to_iir_inverse`（双线性变换） | 保留；新增 `ExponentialIIRDesigner` | 前者保留不动 |
| FIR 残差优化 | 无 | 新增 `FIRResidualDesigner` | 不用即可 |
| 预失真求逆 | `PredistortionDesigner.design()` 返回单 `DistortionModel` | 新增 `method="cascade"` | 不选该方法即旧行为 |
| `_auto_method()` | 对 Exponential→iir_inverse，否则→frequency_inverse | **不变**。cascade 必须显式指定 `method="cascade"`，**不加入 auto** | 无影响 |
| `ControlLine.transfer_function` | 接受 `DistortionModel` | 天然兼容 `CascadeDistortion`（is-a 关系） | 无影响 |

#### 2.1.2 用户 API（预期）

```python
# 旧路径（单滤波器，保留）
designer = PredistortionDesigner(method="frequency_inverse")
inverse = designer.design(transfer_model, dt=0.5)

# 新路径（级联）
designer = PredistortionDesigner(method="cascade",
                                 n_iir_stages=3,
                                 fir_taps=72,
                                 fir_threshold_ns=30.0)
inverse = designer.design(transfer_model, dt=0.5)
# inverse: CascadeDistortion([iir_1, iir_2, iir_3, fir])

# 应用：和单滤波器一致
predistorted = inverse.apply_to_waveform(target)
```

### 2.2 Plan B — 协议驱动的阶跃响应测量

#### 2.2.1 核心变更

| 组件 | 旧行为 | 新行为 | **如何回旧** |
|------|--------|--------|-------------|
| `Experiment` 基类 | 无 `control_line` | `_route_flux()` 助手 + 子类加 `control_line: None`（默认） | 不传 control_line = passthrough |
| 各 `*Experiment.run()` | 直接 `qubit.qubit_in_mag(flux_signal)` | 开头插 `flux = self._route_flux(flux)` | 无 control_line 时为 identity |
| 阶跃响应测量 | 仅 `_measure_step_response` 调解析公式 | 新增 `StepResponseMeasurement`，可选协议 | 不用即可 |
| `WaveformCalibration` | `method="simulation"` 走解析 | 新增 `measurement: None`（默认），不给则走解析 | 默认就是旧行为 |
| `PredistortionValidationWorkflow` | 用解析模型 | 新增 `measurement_protocol: None`（默认），不给则走解析 | 默认就是旧行为 |

#### 2.2.2 用户 API（预期）

```python
# 保留旧路径（解析）
cal = WaveformCalibration(distortion=true_dist, method="transfer_function",
                          fit_type="multi_exp")
table = cal.calibrate()

# 新路径（协议驱动）
meas = StepResponseMeasurement(
    qubit=qubit,
    control_line=ControlLine(name="Z0", kind="z",
                              transfer_function=true_dist, ...),
    protocol="cryoscope",
    step_amplitude=0.05,
    t_max=500.0,
)
cal = WaveformCalibration(measurement=meas, fit_type="multi_exp")
table = cal.calibrate()
```

### 2.3 两块结合后的端到端流程

```
true_distortion
    ↓ [Plan B] StepResponseMeasurement(protocol="cryoscope")
        ↓ ControlLine.apply(V_step) → qubit 经历畸变 flux
        ↓ CryoscopeExperiment.run() → IQ → φ(t)
        ↓ CryoscopeReconstruction.reconstruct() → s_R(t)
s_R(t) （含量子测量的精度损失）
    ↓ [现有] WaveformCalibration._fit_step_response → {a_k, τ_k}
multi_exp 模型
    ↓ [Plan A] PredistortionDesigner(method="cascade")
        ↓ 逐 (a_k, τ_k) 调 ExponentialIIRDesigner → IIR_k
        ↓ 级联模拟残差，剩余分量交给 FIRResidualDesigner
CascadeDistortion([iir_1, ..., iir_N, fir])
    ↓ apply_to_waveform(target_pulse)
预失真后的 AWG 波形
    ↓ [Plan B 闭环验证] 串到 control_line 前面，再次 StepResponseMeasurement
最终阶跃响应 → 评估 ||s_corrected - u(t)||₂
```

### 2.4 向后兼容一览

**所有新增能力均为 opt-in**。不传新参数等于旧行为：

| 你想做什么 | 用哪个路径 | 旧行为是否保留 |
|-----------|-----------|---------------|
| 用实验测阶跃响应（含 control_line 失真） | `Exp(control_line=cl).run()` | ✅ `control_line=None` 时 `_route_flux` 为 identity |
| 用解析模型测阶跃响应（不跑量子模拟） | `WaveformCalibration(distortion=D)` | ✅ 不给 `measurement` 就是旧路径 |
| 用重建协议测阶跃响应 | `WaveformCalibration(measurement=meas)` | ✅ 不给 measurement 不触发 |
| 工作流用解析路径 | `PredistortionValidationWorkflow(...)` | ✅ 不给 `measurement_protocol` 就是旧路径 |
| 工作流跑 Cryoscope 闭环 | `PredistortionValidationWorkflow(measurement_protocol="cryoscope")` | ✅ opt-in |
| 单滤波器预失真（现有一切方法） | `PredistortionDesigner(method="iir_inverse")` 等 | ✅ 全部保留，一个字不动 |
| 级联预失真 | `PredistortionDesigner(method="cascade")` | ✅ 新增，不影响旧方法 |
| auto 自动选方法 | `PredistortionDesigner(method="auto")` | ✅ 行为完全不变 |

---

## §3 受影响文件

### Tier A — Plan A（级联滤波器）

| 文件 | 改动内容 |
|------|---------|
| [sqc/hardware/distortion.py](../../sqc/hardware/distortion.py) | 新增 `CascadeDistortion(DistortionModel)`；实现 `apply` / `step_response` / `impulse_response` / `frequency_response` 四接口（前向逐级 apply，频响逐级相乘） |
| [sqc/calibration/waveform.py](../../sqc/calibration/waveform.py) | 新增 `ExponentialIIRDesigner`（Rol 2020 §IV.B 公式）；新增 `FIRResidualDesigner`（CMA-ES / `scipy.optimize.differential_evolution`）；`PredistortionDesigner` 增加 `method="cascade"` 分支、`n_iir_stages` / `fir_taps` / `fir_threshold_ns` 字段 |
| [sqc/hardware/__init__.py](../../sqc/hardware/__init__.py) | 导出 `CascadeDistortion` |
| [sqc/calibration/__init__.py](../../sqc/calibration/__init__.py) | 导出 `ExponentialIIRDesigner`、`FIRResidualDesigner` |

### Tier B — Plan B（协议驱动测量）

| 文件 | 改动内容 |
|------|---------|
| [sqc/experiments/base.py](../../sqc/experiments/base.py) | `Experiment`（ABC，非 dataclass）新增 `_route_flux(signal)` 助手方法（读 `self.control_line` 属性）；**字段由子类各自添加**（见 P9.B.1）|
| [sqc/hardware/distortion.py](../../sqc/hardware/distortion.py) | `DistortionModel` 新增 `apply_to_signal(FluxSignal) -> FluxSignal` 适配方法（返回 `FluxSignal(type=8, ...)`，**不要** 复用 `apply_to_waveform`——后者会丢失类型）|
| [sqc/experiments/cryoscope.py](../../sqc/experiments/cryoscope.py) | `CryoscopeExperiment` 增加 `control_line` 字段；`run()` 在调 `qubit.qubit_in_mag()` 前插入 `phi_truncated = self._route_flux(phi_truncated)` |
| [sqc/experiments/delay_ramsey.py](../../sqc/experiments/delay_ramsey.py) | `DelayRamseyExperiment` 增加 `control_line` 字段；同上 |
| [sqc/experiments/transient.py](../../sqc/experiments/transient.py) | `TransientSensingExperiment` 增加 `control_line` 字段；同上 |
| [sqc/experiments/pi_pulse_comp.py](../../sqc/experiments/pi_pulse_comp.py) | `PiPulseCompensationExperiment` 增加 `control_line` 字段；同上 |
| [sqc/experiments/ramsey.py](../../sqc/experiments/ramsey.py) | `RamseyExperiment` 增加 `control_line` 字段；同上（向后兼容；不影响纯标定用例） |
| [sqc/experiments/echo.py](../../sqc/experiments/echo.py) | `EchoExperiment` 增加 `control_line` 字段；同上 |
| [sqc/calibration/step_response.py](../../sqc/calibration/step_response.py) | **新建**。包含 `StepResponseMeasurement` 数据类 + `_PROTOCOL_REGISTRY` 字典（协议名 → (Experiment 类, Reconstruction 类) 工厂） |
| [sqc/calibration/waveform.py](../../sqc/calibration/waveform.py) | `WaveformCalibration` 新增 `measurement: StepResponseMeasurement \| None` 字段；`_measure_step_response()` 优先调用 `measurement.measure()` |
| [sqc/workflows/predistortion_validation.py](../../sqc/workflows/predistortion_validation.py) | `PredistortionValidationWorkflow` 新增 `measurement_protocol: str \| None`、`qubit: TransmonQubit \| None`；当 protocol 给定时构造 `StepResponseMeasurement` 走协议路径 |
| [sqc/calibration/__init__.py](../../sqc/calibration/__init__.py) | 导出 `StepResponseMeasurement` |

### Tier C — 测试

| 文件 | 改动内容 |
|------|---------|
| [tests/unit/test_cascade_distortion.py](../../tests/unit/test_cascade_distortion.py) | **新建**。验证 `CascadeDistortion` 的前向/逆向数学正确性（2 阶级联 vs `MultiExponentialDistortion` 并行 sanity check） |
| [tests/unit/test_exponential_iir_designer.py](../../tests/unit/test_exponential_iir_designer.py) | **新建**。验证 Rol 2020 公式：对已知 (A, τ) 的单指数失真，逆滤波后阶跃响应误差 < 0.1% |
| [tests/integration/test_cascade_predistortion.py](../../tests/integration/test_cascade_predistortion.py) | **新建**。3 阶 multi_exp 失真 → 级联预失真 → 残差 RMSE 优于单 `frequency_inverse` 至少 2× |
| [tests/integration/test_step_response_measurement.py](../../tests/integration/test_step_response_measurement.py) | **新建**。逐协议验证 Plan B：Cryoscope / delay Ramsey 协议测得的 s(t) 与解析 s(t) RMSE < 1%（参考 Rol 2020 Fig. S2） |
| [tests/integration/test_predistortion_protocol_driven.py](../../tests/integration/test_predistortion_protocol_driven.py) | **新建**。`PredistortionValidationWorkflow(measurement_protocol="cryoscope")` 端到端闭环验证：improvement_factor > 5× |

### Tier D — 文档

| 文件 | 改动内容 |
|------|---------|
| [docs/architecture.md](../../docs/architecture.md) | 新增 P9 章节：级联架构图 + 协议驱动测量章节；更新 `DistortionModel` 列表加入 `CascadeDistortion`；更新 `Calibration` 章节加入 `StepResponseMeasurement` |
| [idea/_sensing theory.md](../_sensing%20theory.md) | 已完成（"设计滤波器"节已有 IIR 级联 + FIR 描述） |

---

## §4 执行计划

### P9.0: 安全检查点

- [ ] `git add -A && git commit -m "safety checkpoint: before P9 cascade + protocol-driven measurement"`
- [ ] `pytest tests/regression -m regression` 确认全绿
- [ ] 在 `_handoff_state.md` 记录回溯 SHA

### Plan A：级联滤波器（可独立提交）

#### P9.A.1: CascadeDistortion 模型

- [ ] 在 `sqc/hardware/distortion.py` 实现 `CascadeDistortion(DistortionModel)`
  - 字段: `stages: list[DistortionModel] = field(default_factory=list)`
  - `apply(waveform, dt)`: `for s in self.stages: waveform = s.apply(waveform, dt); return waveform`
  - `step_response(t)`: 对单位阶跃信号调 `apply`（数值方式，不用解析合成）
  - `impulse_response(t)`: 对 δ-脉冲调 `apply`
  - `frequency_response(omega)`: `np.prod([s.frequency_response(omega) for s in stages], axis=0)`
- [ ] 单元测试 `test_cascade_distortion.py`：
  - 两个 `SingleExponentialDistortion` 级联 vs 串行手动 apply 一致
  - DC gain = ∏(各级 DC gain)
  - 空 `stages=[]` 等同于恒等系统

#### P9.A.2: ExponentialIIRDesigner

- [ ] 在 `sqc/calibration/waveform.py` 实现 `ExponentialIIRDesigner`
  - 输入: `amplitude: float, tau: float, fs: float`
  - 公式（论文 Eq. S22）：
    ```
    α = 1 - exp(1 / (fs · τ · (1 + A)))
    k = A / ((1 + A)(1 - α))            # if A < 0
    k = A / (1 - α)                      # if A ≥ 0（实际 tail 几乎总是 A > 0）
    b₀ = 1 - k + k·α,  b₁ = -(1 - k)(1 - α)
    a₀ = 1,            a₁ = -(1 - α)
    ```
  - 返回 `IIRDistortion(b_coeffs=[b₀, b₁], a_coeffs=[a₀, a₁])`
- [ ] 单元测试 `test_exponential_iir_designer.py`：
  - 对 `SingleExponentialDistortion(A=0.05, τ=50)` 设计 IIR → 级联 → 阶跃响应应 ≈ u(t)，RMSE < 1e-3
  - 边界：A → 0 时 b ≈ [1, 0], a ≈ [1, 0]（近似恒等）

#### P9.A.3: FIRResidualDesigner

- [ ] 在 `sqc/calibration/waveform.py` 实现 `FIRResidualDesigner`
  - 字段: `n_taps: int = 72`, `optimizer: Literal["cmaes", "differential_evolution"] = "differential_evolution"`
  - 输入: 残差阶跃响应 `s_residual(t)`、采样间隔 `dt`
  - 输出: `FIRDistortion(taps=...)` 使 `s_residual ⋆ h_FIR ≈ u(t)` 在前 30 ns 内
  - 优化目标: `||s_corrected[:n_short] - u(t)[:n_short]||₂`
  - 默认 `differential_evolution`（scipy 自带，无新依赖）；如用户装了 `cma` 包则可切换 cmaes
- [ ] 单元测试：
  - 已知短时标失真，FIR 设计后残差 RMSE 显著下降

#### P9.A.4: PredistortionDesigner 增加 cascade 分支

> **API 现状**：[`PredistortionDesigner.method`](../../sqc/calibration/waveform.py#L320-L322) 当前是 `Literal["auto", "fir_inverse", "iir_inverse", "frequency_inverse"]`。新增 `"cascade"` 需扩展该 Literal。

- [ ] 扩展 `PredistortionDesigner.method` Literal:
  ```python
  method: Literal[
      "auto", "fir_inverse", "iir_inverse",
      "frequency_inverse", "cascade",   # 新增
  ] = "auto"
  ```
- [ ] `PredistortionDesigner` 增加字段:
  - `n_iir_stages: int = 3`
  - `fir_taps: int = 72`
  - `fir_threshold_ns: float = 30.0`
- [ ] `_auto_method()` **保持不动**。cascade 是显式选择，不加入 auto（避免改变 `method="auto"` 的默认行为）
- [ ] `method="cascade"` 流程（新增 `_design_cascade(transfer_model, dt)` 方法）:
  1. **类型检查**：若不是 `MultiExponentialDistortion`，先尝试 `MultiExponentialDistortion` 拟合其阶跃响应；若仍不收敛则 fallback 到 `frequency_inverse`
  2. 按 τ 降序排序 `(amplitudes, taus)`
  3. 选 `τ > fir_threshold_ns` 且总数 ≤ `n_iir_stages` 的分量
  4. 对每个 (a_k, τ_k) 调 `ExponentialIIRDesigner(amplitude=-a_k, tau=τ_k, fs=1/dt)` → `iir_k`
     - **注意符号**：`MultiExponentialDistortion` 的 step response 是 `1 - Σ a_k·exp(-t/τ_k)`（amplitudes 正），而 `ExponentialIIRDesigner` 公式假设 `s(t) = g(1 + A·exp(-t/τ))`——传入时 A 取负号
  5. 用 `CascadeDistortion([iir_1, ..., iir_N])` 对原 transfer 做前向 apply（在单位阶跃信号上），求残差阶跃响应
  6. 若 fir_taps > 0，调 `FIRResidualDesigner` 设计 FIR 补 < 30 ns 残差；否则跳过 FIR
  7. 返回 `CascadeDistortion([iir_1, ..., iir_N] + ([fir] if fir else []))`
- [ ] 集成测试 `test_cascade_predistortion.py`：
  - 输入：3 分量 `MultiExponentialDistortion(amplitudes=[0.06, 0.03, 0.01], taus=[200.0, 50.0, 10.0])`
  - 设计级联 → 应用到方波 → 测得 corrected step response
  - 验证 `RMSE(corrected, ideal_step) < RMSE(frequency_inverse_corrected, ideal_step) / 2`

### Plan B：协议驱动的阶跃响应测量（可独立提交）

#### P9.B.1: Experiment 基类增加 control_line 字段

> **重要**：[`sqc/experiments/base.py`](../../sqc/experiments/base.py) 当前的 `Experiment` 是 `ABC`（**不是** `@dataclass`），不能直接添加字段。子类（如 `CryoscopeExperiment`）才是 dataclass。两种实现方案二选一：
>
> **方案 A（推荐）**：保持 `Experiment` 为 ABC，仅在基类提供 `_route_flux()` 助手方法（读 `self.control_line` 属性，子类各自添加字段）
> **方案 B**：把 `Experiment` 改造为 dataclass 基类（需所有子类协调改造，工作量大）
>
> 选 A 路径。

- [ ] `sqc/experiments/base.py` `Experiment` ABC 增加助手方法（**不是字段**）：
  ```python
  def _route_flux(self, signal):
      """Pass flux signal through optional control_line distortion.

      Subclasses should declare `control_line: object | None = None`
      as a dataclass field. None ⇒ passthrough.
      """
      cl = getattr(self, "control_line", None)
      if cl is None or getattr(cl, "transfer_function", None) is None:
          return signal.copy()
      return cl.transfer_function.apply_to_signal(signal)
  ```
- [ ] 在每个 `*Experiment` 子类的 dataclass 字段列表追加（与现有字段并列）：
  ```python
  control_line: object | None = None  # ControlLine, duck-typed 避免循环依赖
  ```
- [ ] **不删除** 任何现有字段；向后兼容（不传 `control_line` 即旧行为）

#### P9.B.2: DistortionModel.apply_to_signal 适配

> **重要**：`FluxSignal` 继承自 `Waveform`，使用 `.samples`（**不是** `.signal`，那是 `src/signal.py` 旧 API）。`apply_to_waveform` 当前返回 `Waveform`，会丢失 `FluxSignal` 类型——`apply_to_signal` 必须保持原类型。

- [ ] `sqc/hardware/distortion.py` 的 `DistortionModel` 基类增加:
  ```python
  def apply_to_signal(self, signal):
      """Apply distortion to a FluxSignal, preserving its type.

      Implementation note: `FluxSignal.copy()` rebuilds from `_params`,
      so mutating `.samples` afterwards leaves params/samples
      inconsistent. Instead, rebuild as `type=8` (user-defined) with
      the distorted samples — this is the canonical "raw samples"
      pathway used elsewhere in sqc (see delay_ramsey.py:138).

      Parameters
      ----------
      signal : FluxSignal

      Returns
      -------
      FluxSignal (type=8) with distorted samples.
      """
      from sqc.control.flux_signal import FluxSignal

      dt = float(signal.t_list[1] - signal.t_list[0])
      out = self.apply(np.asarray(signal.samples), dt)
      return FluxSignal(type=8, t_list=signal.t_list.copy(), signal=out)
  ```
- [ ] 单元测试：
  - 与 `apply_to_waveform` 数值一致（仅比较 `.samples`）
  - 输入 `FluxSignal` 时输出仍为 `FluxSignal`（`isinstance(out, FluxSignal) == True`）
  - 输出的 `.samples` 与输入 `apply()` 的结果按位相等
- [ ] **DECISION_NEEDED**：若上游协议依赖 `signal._type`（例如某些重建走 type-specific 分支），强制 `type=8` 可能破坏行为。建议在子代理执行时先 grep `self._type` / `signal._type` / `_params["amplitude"]` 等使用模式，确认无 type-specific 依赖再合入；若有则改回 `signal.copy()` + 直接覆盖 `.samples`。

#### P9.B.3: 修改各 Experiment.run() 调 _route_flux

按顺序逐个修改并测试，每改完一个跑该实验的 baseline test：

- [ ] `CryoscopeExperiment.run()`: 在 `self.qubit.qubit_in_mag(phi_truncated, ...)` 之前插 `phi_truncated = self._route_flux(phi_truncated)`
- [ ] `DelayRamseyExperiment.run()`: 同样位置
- [ ] `TransientSensingExperiment.run()`: 同样位置（注意 sliding measurement 内部循环，确保畸变只施加一次而非每次滑动重复）
- [ ] `PiPulseCompensationExperiment.run()`: 同样位置
- [ ] `RamseyExperiment.run()`: 同样位置（如果不影响纯频率标定 baseline，验证向后兼容）
- [ ] `EchoExperiment.run()`: 同样位置

每个 baseline 必须保持不变（前提：调用方未传 `control_line`）。

#### P9.B.4: StepResponseMeasurement

> **类名核对**（已确认）：
> - 实验类：`CryoscopeExperiment`、`DelayRamseyExperiment`、`TransientSensingExperiment`、`PiPulseCompensationExperiment`
> - 重建类：`CryoscopeReconstruction`、`DelayRamseyReconstruction`、`TransientReconstruction`、`PiPulseCompReconstruction`
> - 输出对象：`exp.run()` 返回 `ExperimentResult`；`rec.reconstruct(result)` 返回 `FluxSignal`（使用 `.samples`）

- [ ] 新建 `sqc/calibration/step_response.py`，实现 `StepResponseMeasurement` 数据类（`@dataclass`）
- [ ] `_PROTOCOL_REGISTRY`（模块级常量）:
  ```python
  from sqc.experiments.cryoscope import CryoscopeExperiment
  from sqc.experiments.delay_ramsey import DelayRamseyExperiment
  from sqc.experiments.transient import TransientSensingExperiment
  from sqc.experiments.pi_pulse_comp import PiPulseCompensationExperiment
  from sqc.reconstruction.cryoscope import CryoscopeReconstruction
  from sqc.reconstruction.delay_ramsey import DelayRamseyReconstruction
  from sqc.reconstruction.transient import TransientReconstruction
  from sqc.reconstruction.pi_pulse_comp import PiPulseCompReconstruction

  _PROTOCOL_REGISTRY = {
      "cryoscope":    (CryoscopeExperiment,           CryoscopeReconstruction),
      "delay_ramsey": (DelayRamseyExperiment,         DelayRamseyReconstruction),
      "transient":    (TransientSensingExperiment,    TransientReconstruction),
      "pi_pulse":     (PiPulseCompensationExperiment, PiPulseCompReconstruction),
  }
  ```
- [ ] `measure()` 方法:
  1. 构造 `FluxSignal(type=1, t_list=np.arange(0, t_max, dt), amplitude=step_amplitude)`
  2. 取出 `(Exp, Rec) = _PROTOCOL_REGISTRY[self.protocol]`
  3. 构造 `exp = Exp(qubit=self.qubit, control_line=self.control_line, flux_signal=v_in, ...)`
  4. `result = exp.run()`
  5. 按协议构造 `rec`（参数从 `self` 或 `result.metadata` 取，例：Cryoscope 需 `tau` 和 `inversion="response"` + `qubit`）
  6. `flux_R = rec.reconstruct(result)`
  7. 归一化: 返回 `(flux_R.t_list, flux_R.samples / self.step_amplitude)`
- [ ] 协议特定的默认参数（在 `__post_init__` 中按协议选择）:
  - cryoscope: `qubit.flux = 0`（甜点）
  - delay_ramsey: `qubit.flux = κ_max_flux`（最优灵敏度）
  - transient: `qubit.flux = κ_max_flux`
  - pi_pulse: `qubit.flux = κ_max_flux`
  - `κ_max_flux` 可由 `qubit.frequency_sensitivity(flux)` 在 `np.linspace(0, 0.25, 100)` 上扫描后取 argmax 得到
- [ ] **DECISION_NEEDED**：是否在 `StepResponseMeasurement` 内自动设置 `qubit.flux`，还是要求调用方提前设置？保守做法：调用方负责设置；`__post_init__` 只校验当前 flux 是否符合协议预期，不符合时打 warning。

#### P9.B.5: WaveformCalibration 接入测量回调

> **现状**：[`WaveformCalibration`](../../sqc/calibration/waveform.py#L26) 已是 `@dataclass`，字段顺序: `method, distortion, fit_type, n_exp_components, t_max, n_points, transfer_model, predistortion_method, n_taps, regularization, dt`。新增 `measurement` 字段需放在已有字段之后（保持向后兼容）。

- [ ] `WaveformCalibration` 增加字段（追加在 `dt: float | None = None` 之后）：
  ```python
  # -- protocol-driven measurement (P9.B) --
  measurement: object | None = None  # StepResponseMeasurement, duck-typed
  ```
- [ ] 重写 `_measure_step_response()`:
  ```python
  def _measure_step_response(self) -> tuple[np.ndarray, np.ndarray]:
      if self.measurement is not None:
          return self.measurement.measure()
      # 保留旧路径（解析）
      t = np.linspace(0, self.t_max, self.n_points)
      step = self.distortion.step_response(t)
      return t, step
  ```
- [ ] 兼容性：
  - `distortion=None` 且 `measurement is not None` 时合法（旧调用方不受影响）
  - `distortion is not None` 且 `measurement is not None` 时：优先用 `measurement`（隐式覆盖解析路径）；可选打 warning 提示二选一
  - 两者都 None 时：保持现有错误行为（`step_response(t)` 抛 AttributeError）
- [ ] 集成测试 `test_step_response_measurement.py`：
  - 各协议测得的 s(t) 与解析 s(t) RMSE < 1%

#### P9.B.6: PredistortionValidationWorkflow 端到端协议驱动

- [ ] 增加字段:
  ```python
  qubit: object | None = None
  measurement_protocol: str | None = None  # None ⇒ 解析路径
  ```
- [ ] `run()` 中 step 3 改为:
  ```python
  if self.measurement_protocol is None:
      cal = WaveformCalibration(distortion=self.true_distortion,
                                method="transfer_function", ...)
  else:
      meas = StepResponseMeasurement(qubit=self.qubit, control_line=line,
                                      protocol=self.measurement_protocol, ...)
      cal = WaveformCalibration(measurement=meas, ...)
  ```
- [ ] 集成测试 `test_predistortion_protocol_driven.py`：
  - `measurement_protocol="cryoscope"` → improvement_factor > 5×

### P9.C: 联合演示（可选，作为 Notebook 收尾）

- [ ] 在 `Simulation_sqc.ipynb` 新增 section：
  1. 真实 multi_exp 失真注入 `ControlLine`
  2. `StepResponseMeasurement(protocol="cryoscope")` 测 s(t)
  3. `PredistortionDesigner(method="cascade")` 设计级联逆
  4. 同一 protocol 测预失真后的 s(t)，画对比图（参考 Rol 2020 Fig. 2(a)）

---

## §5 验证标准

### Plan A
- [ ] `test_cascade_distortion.py`：级联 = 串行手动 apply，DC gain 正确
- [ ] `test_exponential_iir_designer.py`：单指数补偿残差 RMSE < 1e-3
- [ ] `test_cascade_predistortion.py`：3 阶 multi_exp 失真，级联预失真 RMSE 优于单频域逆至少 2×
- [ ] 现有 `test_transfer_function_cal.py` 全绿（向后兼容）

### Plan B
- [ ] 所有现有 `*Experiment` baseline test 全绿（不传 `control_line` 即旧行为）
- [ ] `test_step_response_measurement.py`：Cryoscope / delay Ramsey 测得 s(t) 与解析 RMSE < 1%
- [ ] `test_predistortion_protocol_driven.py`：`measurement_protocol="cryoscope"` 闭环 improvement_factor > 5×

### 联合
- [ ] `Simulation_sqc.ipynb` Notebook 跑通，对比图复现 Rol 2020 Fig. 2(a) 趋势
- [ ] `pytest tests/regression -m regression` 全绿
- [ ] `git diff --quiet master -- src/` 返回 0（R1）

---

## §6 依赖与协调

### 6.1 P9 内部：Plan A 与 Plan B 的关系

**两者完全独立，可并行实施**。

| Plan | 改的对象 | 是否依赖另一个 |
|------|---------|---------------|
| Plan A | 滤波器**设计**逻辑（distortion.py / waveform.py） | 否——可以用解析阶跃响应单独验证 |
| Plan B | 阶跃响应**测量**链路（experiments / step_response.py） | 否——可以用现有单滤波器 `PredistortionDesigner` 单独验证 |

实施建议：
- **小型团队**: 先 Plan A 完整提交，再 Plan B；最后跑 §4.C 联合演示
- **并行团队**: A 走 distortion.py / calibration/waveform.py，B 走 experiments/* + step_response.py，**仅共用 `WaveformCalibration` 一个文件**——B 加 `measurement` 字段，A 加 `method="cascade"` 字段，互不冲突

### 6.2 与 Phase 8 的关系

**P9 不强依赖 P8**。两者改的对象完全不同：

| Phase | 主要改动对象 |
|-------|------------|
| P8 | `delay_ramsey.py` / `ramsey.py` 的 filter function 路径 |
| P9 | `*Experiment` 注入 `control_line` 失真；新增 cascade 设计 |

**P8 与 P9.B 在 `delay_ramsey.py` 和 `ramsey.py` 上有文件级交集，但行级隔离**：
- P8 改 `run()` 内的 flux 调制函数（zero-out vs filter function）
- P9.B 在 `run()` 开头新增 `flux = self._route_flux(flux)`，位置在 P8 改动之前

#### 情形 1：P9 先实施，P8 后实施

P9 完成后，`delay_ramsey.py / ramsey.py` 已经有 `_route_flux()` 调用。P8 实施时：
- P8 的 `use_filter_function` 分支照常添加，**位置在 `_route_flux()` 之后**
- P8 的 kernel 计算用的 flux signal 是**经 control_line 畸变后**的版本——这是正确的物理（kernel 应当对真实到达 qubit 的信号求得）
- **P8 不需要做调整**

#### 情形 2：P8 先实施，P9 后实施

P8 完成后，`delay_ramsey.py / ramsey.py` 已经有 `use_filter_function` 分支。P9.B 实施时：
- P9.B 在 `run()` 开头插 `flux = self._route_flux(flux)`，位置**早于** P8 的 zero-out / 连续施加分歧点
- 后续两条路径都基于畸变后的 flux 工作——物理上正确
- **P9 不需要做调整**

#### 情形 3：同 phase 并行

两个修改点行级不重叠（P8 改 `run()` 中段的 flux modulation 逻辑，P9.B 改 `run()` 开头的 flux routing），merge 时不会冲突。但为减少代码 review 复杂度，建议 P8.1 / P8.4 完成后再做 P9.B.3。

### 6.3 与 P7（时间轴统一）的关系

非硬依赖。P9 使用 `CONFIG.awg.dt` 派生时间轴（遵循 R9），与 P7 是否完成无关。

### 6.4 风险

1. **CMA-ES 依赖**：`cma` 包是新依赖（违反 R10）。**解决**：默认用 scipy `differential_evolution`，`cma` 作可选优化器（用户安装后才启用）
2. **协议驱动测量性能**：每次测一条阶跃响应需要跑完整 Cryoscope / delay Ramsey 序列（数百次 mesolve），单次约 30 秒——比解析方案慢 100×。**应对**：默认仍用 `analytical`；`protocol` 模式只在集成测试和 Notebook 用
3. **kernel 缓存与 control_line 失真的交互**：P8 计算的 kernel 应基于 control_line 畸变后的 flux——若 P8 直接缓存原始信号的 kernel，会出错。**应对**：在 P9.B.3 修改 delay_ramsey 时，确保 kernel 计算用的是 `_route_flux()` 后的信号
4. **Rol 2020 公式适用范围**：公式假设 A > 0 时 `k = A/((1+A)(1-α))`，A < 0 时同公式。论文实际只验证了 |A| ≤ 0.1 区间。**应对**：在 `ExponentialIIRDesigner` 内对 A > 0.5 或 A < -0.5 打 warning

---

## §7 不做的事

- **不修改** `src/` 下任何文件（R1）——`Protocal` 拼写、`sliding_measrement` 拼写保持原样
- **不引入** `cma` / `pydantic` / `yaml` 等新依赖（R10）；CMA-ES 用 scipy `differential_evolution` 作默认
- **不删除** 任何现有 `PredistortionDesigner` 求逆方法（保留 `iir_inverse` / `fir_inverse` / `frequency_inverse`）
- **不强迫** 所有 `Experiment` 调用方传 `control_line`——默认 `None` 即旧行为，向后兼容
- **不在** P9 内推广到多 qubit 串扰矩阵（`_sensing theory.md` "多 qubit 推广"一节属于未来 phase）
- **不修改** 现有 baseline pickle——所有新增能力默认走可选路径，不影响 `pytest tests/regression -m regression`
- **不要求** Plan A 和 Plan B 必须同 commit 提交——它们彼此独立，可分批 PR

---

## §8 阶段交接

完成 P9 后必须更新 `_handoff_state.md`：

```markdown
## P9 完成 — YYYY-MM-DD

### 新增能力
- CascadeDistortion + IIR 级联 + FIR 残差修正（Rol 2020）
- 协议驱动的阶跃响应测量（Cryoscope / delay Ramsey / transient / pi_pulse）
- PredistortionValidationWorkflow 端到端闭环

### 关键决策
- CMA-ES 用 scipy differential_evolution 作默认优化器（避免新依赖）
- Experiment 的 control_line 字段默认 None（向后兼容）
- StepResponseMeasurement 与 WaveformCalibration 解耦（用户可单用测量器）

### 下一步建议
- P10 候选: 多 qubit 串扰矩阵标定（_sensing theory.md 的多 qubit 推广节）
- 或: cma 包可选支持，提高 FIR 优化收敛速度
```
