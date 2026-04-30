# Phase 4 Handbook — ControlLine + DistortionModel + 预失真

> 前置阅读:[_refactor_plan.md](_refactor_plan.md) §1.3、§4–§6.2、[phase_3_handbook.md](phase_3_handbook.md)  
> 估计工时:5–7 天  
> 触发条件:Phase 3 完成 + Track B `_TODO_master.md 1.3` (DistortionModel 在 src/ 内已实现)  
> 完成标志:`ControlLine`、`DistortionModel` 子类、`PredistortionDesigner`、predistortion validation workflow 全部实现并通过 baseline

---

## 1. 目标

把项目推进到真实超导量子计算机控制链路层。本 phase 是**新增功能**(不是从 src/ 迁移),因为 Track B 1.3 在 src/ 中已完成基础失真模型,但**预失真设计与验证 workflow 主要在 sqc/ 中新增**。

主要交付:
1. `ControlLine` 完整实现(P1 只是 dataclass 占位)。
2. `DistortionModel` ABC + 子类:`SingleExponentialDistortion`、`MultiExponentialDistortion`、`FIRDistortion`、`IIRDistortion`、`CustomTransferDistortion`。
3. `TransferFunctionCalibration` — 测量 step response、提取 H(ω)。
4. `PredistortionDesigner` — 根据测得 H(ω) 设计 FIR/IIR 补偿滤波器。
5. `PredistortionValidationWorkflow` — 端到端补偿前后对比。
6. 把 Track B 在 src/ 中实现的 `DistortionModel` 内化(若 Track B 把它放在 src/signal.py 或 src/distortion.py)。

---

## 2. 前置条件

| # | 条件 | 验证 |
|---|---|---|
| 2.1 | Phase 3 完成 | LM/Cryoscope 内化完成,baseline 通过 |
| 2.2 | Track B `_TODO_master.md 1.3.1–1.3.4` 完成 | src/ 中存在 DistortionModel 类(可能在 src/signal.py 或 src/distortion.py) |
| 2.3 | (推荐) Track B `_TODO_master.md 2.3` (预失真应用) 至少完成 2.3.1 | 模拟失真环境已有 |

---

## 3. 任务清单

### 3.1 完整实现 ControlLine

**文件**:`sqc/hardware/control_line.py`

```python
# sqc/hardware/control_line.py (P4 完整版)
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from sqc.control.waveform import Waveform


@dataclass
class ControlLine:
    """Models a physical control line (xy / z / readout) with
    optional transfer function distortion.
    """
    name: str
    kind: Literal["xy", "z", "readout"]
    source: str                       # e.g., "AWG0:CH1"
    target: str                       # e.g., "Q0", "C01"
    transfer_function: "DistortionModel | None" = None
    metadata: dict = field(default_factory=dict)
    
    def apply(self, awg_waveform: Waveform) -> Waveform:
        """Forward-model: AWG waveform → on-chip waveform.
        
        If no transfer_function is set, returns awg_waveform unchanged.
        """
        if self.transfer_function is None:
            return awg_waveform.copy()
        return self.transfer_function.apply(awg_waveform)
    
    def predistort(self, target_waveform: Waveform,
                   designer: "PredistortionDesigner") -> Waveform:
        """Inverse: target on-chip waveform → required AWG waveform."""
        return designer.predistort(target_waveform, self.transfer_function)
```

### 3.2 DistortionModel 子类

**文件**:`sqc/hardware/distortion.py`

P1 阶段只有 ABC,P4 完整实现五个子类。

```python
# sqc/hardware/distortion.py (P4 完整版)
from __future__ import annotations
from abc import abstractmethod
from dataclasses import dataclass, field

import numpy as np
from scipy.signal import lfilter, freqz

from sqc.control.waveform import Waveform
from sqc.hardware.distortion_base import DistortionModel  # ABC from P1


@dataclass
class SingleExponentialDistortion(DistortionModel):
    """h(t) = (1 - amp) * δ(t) + (amp/τ) * exp(-t/τ) * Θ(t)
    
    A common model for slow charge-redistribution tails.
    """
    amplitude: float                # tail amplitude (dimensionless)
    tau: float                      # tail time constant (ns)
    
    def apply(self, awg_waveform: Waveform) -> Waveform:
        t = awg_waveform.t_list
        dt = t[1] - t[0]
        # Convolve via single-exponential filter
        b, a = self._iir_coeffs(dt)
        out = lfilter(b, a, awg_waveform.samples)
        return Waveform(t_list=t.copy(), samples=out,
                         metadata={**awg_waveform.metadata,
                                   "distortion": "SingleExp"})
    
    def step_response(self, t: np.ndarray) -> np.ndarray:
        return 1 + self.amplitude * (np.exp(-t / self.tau) - 1)
    
    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        h = np.zeros_like(t)
        h[0] = 1.0 - self.amplitude
        h += (self.amplitude / self.tau) * np.exp(-t / self.tau)
        return h
    
    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        return (1 - self.amplitude) + self.amplitude / (1 + 1j * omega * self.tau)
    
    def _iir_coeffs(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        # Implementation: discretize via bilinear/impulse-invariant
        # ... per Track B's implementation ...
        ...


@dataclass
class MultiExponentialDistortion(DistortionModel):
    """h(t) = δ(t) - sum_k amp_k * (1 - exp(-t/τ_k)) for slow tails."""
    amplitudes: np.ndarray            # shape (K,)
    taus: np.ndarray                   # shape (K,)
    
    def apply(self, awg_waveform: Waveform) -> Waveform: ...
    def step_response(self, t: np.ndarray) -> np.ndarray: ...
    def impulse_response(self, t: np.ndarray) -> np.ndarray: ...
    def frequency_response(self, omega: np.ndarray) -> np.ndarray: ...


@dataclass
class FIRDistortion(DistortionModel):
    """Generic FIR filter h[n], n=0..len(taps)-1."""
    taps: np.ndarray
    dt: float
    
    def apply(self, awg_waveform: Waveform) -> Waveform: ...
    def step_response(self, t: np.ndarray) -> np.ndarray: ...
    def impulse_response(self, t: np.ndarray) -> np.ndarray: ...
    def frequency_response(self, omega: np.ndarray) -> np.ndarray: ...


@dataclass
class IIRDistortion(DistortionModel):
    """Generic IIR filter b/a coefficients."""
    b: np.ndarray
    a: np.ndarray
    dt: float
    
    def apply(self, awg_waveform: Waveform) -> Waveform: ...
    def step_response(self, t: np.ndarray) -> np.ndarray: ...
    def impulse_response(self, t: np.ndarray) -> np.ndarray: ...
    def frequency_response(self, omega: np.ndarray) -> np.ndarray: ...


@dataclass
class CustomTransferDistortion(DistortionModel):
    """User-supplied complex frequency-domain H(ω) on a fixed grid."""
    omega: np.ndarray                  # rad/ns
    H: np.ndarray                       # complex
    
    def apply(self, awg_waveform: Waveform) -> Waveform:
        # FFT, multiply by H, inverse FFT
        ...
    def step_response(self, t: np.ndarray) -> np.ndarray: ...
    def impulse_response(self, t: np.ndarray) -> np.ndarray: ...
    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        return np.interp(omega, self.omega, self.H.real) + \
               1j * np.interp(omega, self.omega, self.H.imag)
```

注:DistortionModel ABC 应该挪到独立文件 `sqc/hardware/distortion_base.py`(避免子类与 ABC 同文件循环引用),或让 ABC 与子类共享 `distortion.py`。本方案采用第二种(共享一个文件)。

### 3.3 TransferFunctionCalibration

**文件**:`sqc/calibration/transfer_function.py`

```python
# sqc/calibration/transfer_function.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from sqc.devices.transmon import TransmonQubit
from sqc.calibration.base import Calibration, CalibrationTable
from sqc.calibration.flux_response import FluxResponseCalibration
from sqc.experiments.cryoscope import CryoscopeExperiment
from sqc.experiments.transient import TransientSensingExperiment
from sqc.reconstruction.cryoscope import CryoscopeReconstruction
from sqc.reconstruction.wiener import WienerReconstruction


@dataclass
class TransferFunctionCalibration(Calibration):
    """Measure step response of a control line, fit transfer function H(ω).
    
    Strategy:
    1. Apply known step / impulse via ControlLine.
    2. Reconstruct on-chip waveform via Cryoscope (slow) or Transient (fast).
    3. Compare reconstructed vs target → measured h(t) or H(ω).
    4. Fit IIR / multi-exp / FIR model.
    """
    qubit: TransmonQubit
    method: Literal["cryoscope", "transient", "hybrid"] = "hybrid"
    fit_type: Literal["multi_exp", "iir", "fir", "custom"] = "multi_exp"
    n_exp_components: int = 3
    
    def calibrate(self) -> CalibrationTable:
        # 1. Measure step response
        step_response = self._measure_step_response()
        # 2. Fit chosen model
        params = self._fit(step_response)
        return CalibrationTable(
            qubit_name=self.qubit.spec().name,
            kind="transfer_function",
            inputs=step_response.t_list,
            outputs=step_response.samples,
            fit_params=params,
            metadata={"method": self.method, "fit_type": self.fit_type},
        )
    
    def _measure_step_response(self):
        if self.method == "cryoscope":
            ...
        elif self.method == "transient":
            ...
        elif self.method == "hybrid":
            # Cryoscope for slow tails, transient for fast edge
            ...
    
    def _fit(self, step_response) -> dict:
        if self.fit_type == "multi_exp":
            from scipy.optimize import curve_fit
            ...
        elif self.fit_type == "iir":
            ...
        elif self.fit_type == "fir":
            ...
        elif self.fit_type == "custom":
            return {"H_complex": np.fft.fft(step_response.samples)}
```

### 3.4 PredistortionDesigner

**文件**:`sqc/calibration/predistortion.py`

```python
# sqc/calibration/predistortion.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from sqc.control.waveform import Waveform
from sqc.hardware.distortion import (
    DistortionModel, FIRDistortion, IIRDistortion,
)


@dataclass
class PredistortionDesigner:
    """Designs FIR/IIR predistortion filter given measured transfer function.
    
    Predistortion: V_comp(ω) = H(ω)^(-1) * Phi_target(ω)
    """
    method: Literal["fir_inverse", "iir_inverse", "frequency_inverse"] = "fir_inverse"
    n_taps: int = 64                    # for FIR
    regularization: float = 1e-4        # avoid noise amplification at notches
    
    def design(self, transfer_model: DistortionModel,
               omega_grid: np.ndarray | None = None) -> DistortionModel:
        """Returns a new DistortionModel that approximates H(ω)^(-1).
        
        Apply this returned model to the target waveform to get the AWG
        waveform that will produce that target on-chip.
        """
        if omega_grid is None:
            omega_grid = np.linspace(-np.pi, np.pi, 4096)
        H = transfer_model.frequency_response(omega_grid)
        # Regularized inverse (Wiener-style): H_inv = conj(H) / (|H|^2 + λ²)
        H_inv = np.conj(H) / (np.abs(H)**2 + self.regularization**2)
        
        match self.method:
            case "frequency_inverse":
                return CustomTransferDistortion(
                    omega=omega_grid, H=H_inv,
                )
            case "fir_inverse":
                # IFFT to get time-domain inverse, truncate to n_taps
                h_inv = np.fft.fftshift(np.fft.ifft(H_inv).real)
                h_inv_trunc = h_inv[len(h_inv)//2 - self.n_taps//2 :
                                     len(h_inv)//2 + self.n_taps//2]
                return FIRDistortion(
                    taps=h_inv_trunc,
                    dt=2 * np.pi / (omega_grid[-1] - omega_grid[0]),
                )
            case "iir_inverse":
                # Fit IIR model to H_inv (e.g., via vector fitting or Prony)
                ...
    
    def predistort(self, target: Waveform, transfer: DistortionModel) -> Waveform:
        inverse_model = self.design(transfer)
        return inverse_model.apply(target)
```

### 3.5 PredistortionValidationWorkflow

**文件**:`sqc/workflows/predistortion_validation.py`

```python
# sqc/workflows/predistortion_validation.py
from __future__ import annotations
from dataclasses import dataclass, field

import numpy as np

from sqc.devices.transmon import TransmonQubit
from sqc.control.waveform import Waveform
from sqc.hardware.control_line import ControlLine
from sqc.hardware.distortion import DistortionModel
from sqc.calibration.predistortion import PredistortionDesigner
from sqc.calibration.transfer_function import TransferFunctionCalibration
from sqc.experiments.cryoscope import CryoscopeExperiment
from sqc.workflows.base import Workflow


@dataclass
class PredistortionValidationWorkflow(Workflow):
    """End-to-end predistortion validation:
    
    1. Inject known distortion into ControlLine.
    2. Measure step response via CryoscopeExperiment.
    3. Fit H(ω) via TransferFunctionCalibration.
    4. Design inverse via PredistortionDesigner.
    5. Re-measure with predistorted AWG waveform.
    6. Compute residual: RMSE, settling time, frequency error.
    """
    qubit: TransmonQubit
    target_waveform: Waveform
    true_distortion: DistortionModel    # ground truth (for simulation)
    designer: PredistortionDesigner = field(default_factory=PredistortionDesigner)
    
    def run(self) -> dict:
        # 1. Build control line with the true distortion
        line = ControlLine(name="Z0", kind="z", source="AWG0",
                            target="Q0",
                            transfer_function=self.true_distortion)
        
        # 2. Without predistortion: AWG → on-chip is distorted
        on_chip_uncorrected = line.apply(self.target_waveform)
        
        # 3. Measure transfer function
        cal = TransferFunctionCalibration(
            qubit=self.qubit, method="cryoscope", fit_type="multi_exp",
        )
        measured_table = cal.calibrate()
        # Convert table → DistortionModel
        measured_model = self._table_to_model(measured_table)
        
        # 4. Design inverse
        inverse_model = self.designer.design(measured_model)
        
        # 5. Predistort
        awg_predistorted = inverse_model.apply(self.target_waveform)
        on_chip_corrected = line.apply(awg_predistorted)
        
        # 6. Compute residuals
        rmse_uncorrected = self._rmse(on_chip_uncorrected, self.target_waveform)
        rmse_corrected = self._rmse(on_chip_corrected, self.target_waveform)
        settle_uncorrected = self._settling_time(on_chip_uncorrected,
                                                   self.target_waveform)
        settle_corrected = self._settling_time(on_chip_corrected,
                                                 self.target_waveform)
        
        return {
            "target": self.target_waveform,
            "on_chip_uncorrected": on_chip_uncorrected,
            "on_chip_corrected": on_chip_corrected,
            "awg_predistorted": awg_predistorted,
            "metrics": {
                "rmse_uncorrected": rmse_uncorrected,
                "rmse_corrected": rmse_corrected,
                "improvement_factor": rmse_uncorrected / max(rmse_corrected, 1e-30),
                "settling_uncorrected_ns": settle_uncorrected,
                "settling_corrected_ns": settle_corrected,
            },
        }
    
    @staticmethod
    def _rmse(measured: Waveform, target: Waveform) -> float: ...
    @staticmethod
    def _settling_time(w: Waveform, target: Waveform,
                        tolerance: float = 0.001) -> float: ...
    @staticmethod
    def _table_to_model(table) -> DistortionModel: ...
```

### 3.6 把 Track B 1.3 的 src/ 实现内化

如果 Track B 把 DistortionModel 放在 `src/signal.py` 或新建的 `src/distortion.py`:

1. 把那个文件的 `DistortionModel` 类逐字符 port 到 `sqc/hardware/distortion.py`(选最贴合的子类)。
2. 在 `src/` 中加镜像导出。
3. 跑 baseline 验证。

如果 Track B 1.3 中 `apply()` 方法依赖了 `Signal` 类:

- P4 内化时,新版 `apply` 输入是 `Waveform`,但 facade 中应支持 `Signal` 输入(其 type 8 + samples 等价于 Waveform)。

### 3.7 把 PredistortionDesigner 暴露给旧 API(可选)

若 Track B 在 src/analysis.py 中实现了 `design_iir_correction`(_TODO_master.md 1.3.4),P4 通过 `Analysis.design_iir_correction` facade 转发:

```python
# src/analysis.py (P4 增量)
def design_iir_correction(measured_step_response, ...):
    """Backward-compat facade."""
    from sqc.calibration.predistortion import PredistortionDesigner
    designer = PredistortionDesigner(method="iir_inverse")
    ...
```

---

## 4. 验收标准

### 4.1 失真模型基本功能

```bash
python -c "
import numpy as np
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion

t = np.linspace(0, 100, 1000)
target = Waveform(t_list=t, samples=np.where(t > 10, 1.0, 0.0))
dist = SingleExponentialDistortion(amplitude=0.05, tau=20.0)
distorted = dist.apply(target)
print('step settles to', distorted.samples[-1], '(should be 1.0)')
"
```

期望:`step settles to ~1.0`(±0.001)。

### 4.2 预失真闭环验证

```bash
python -c "
... 完整 PredistortionValidationWorkflow run ...
"
```

期望:`improvement_factor > 10`(补偿后 RMSE 比未补偿低至少 10 倍)。

### 4.3 ControlLine 端到端

```bash
python -c "
import numpy as np
from sqc.hardware.control_line import ControlLine
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.control.waveform import Waveform
line = ControlLine(name='Z0', kind='z', source='AWG0', target='Q0',
                    transfer_function=SingleExponentialDistortion(0.05, 20))
t = np.linspace(0,100,1000)
w = Waveform(t_list=t, samples=np.where(t>10, 1.0, 0.0))
out = line.apply(w)
print('control line OK, on-chip shape:', out.samples.shape)
"
```

### 4.4 物理回归(关键)

```bash
pytest tests/regression -m regression -v
```

P4 引入新 baseline:
- `predistortion_default.pkl`(对一个标准失真+目标信号的全流程结果)

### 4.5 web_demo / Notebook

- 增加 web_demo 中的 "Predistortion" tab(可选)。
- Simulation.ipynb 中能跑通预失真闭环 demo cell。

---

## 5. 测试要求

### 5.1 必须新增的测试

| 文件 | 用例 |
|---|---|
| `tests/unit/test_distortion_models.py` | 五个子类各自的 step_response/impulse_response/frequency_response 数学正确性 |
| `tests/unit/test_control_line.py` | apply / predistort 行为 |
| `tests/unit/test_predistortion_designer.py` | FIR/IIR/frequency 三种方法都能产生有效的逆模型;逆模型 × 原模型 ≈ 单位响应 |
| `tests/integration/test_transfer_function_cal.py` | 已知失真→cryoscope 测量→Calibration 表→重建模型,fit 误差 < 5% |
| `tests/integration/test_predistortion_validation.py` | 端到端 RMSE 改善 > 10x |
| `tests/regression/test_predistortion_baseline.py` | 与 baseline 比对 |

### 5.2 必须保持的回归

所有已有 baseline 仍然 pass。

---

## 6. 风险与回滚

| # | 风险 | 缓解 |
|---|---|---|
| 6.1 | FIR 逆滤波在 |H(ω)| ≈ 0 处放大噪声 | regularization 默认 1e-4;在 designer 中显式限频带 |
| 6.2 | IIR 逆滤波不稳定(极点出单位圆) | 设计后检查极点,失败时 fallback 到 FIR |
| 6.3 | TransferFunctionCalibration 的 hybrid 方法在不同时间尺度上拼接出错 | 显式定义 cross-over frequency,加单元测试覆盖 |
| 6.4 | DistortionModel.apply 在 Waveform 时间步不一致时出错 | 在 apply 内 assert dt 一致;不一致则抛 ValueError |
| 6.5 | Track B 1.3 的 src/ 实现与本方案接口不一致 | P4 开始前先盘点 Track B 实现的 API,在 sqc/ 中保持兼容 |

**回滚**:本 phase 全部新增,不影响 src/。任何模块出问题,删除对应 sqc/ 文件即可。

---

## 7. 输出物

### 7.1 文件清单(新增)

```
sqc/hardware/control_line.py            (P1 占位 → P4 完整)
sqc/hardware/distortion.py              (P1 ABC → P4 完整含 5 个子类)
sqc/calibration/transfer_function.py
sqc/calibration/predistortion.py
sqc/workflows/predistortion_validation.py
tests/unit/test_distortion_models.py
tests/unit/test_control_line.py
tests/unit/test_predistortion_designer.py
tests/integration/test_transfer_function_cal.py
tests/integration/test_predistortion_validation.py
tests/regression/test_predistortion_baseline.py
tests/baselines/predistortion_default.pkl
```

### 7.2 文件清单(修改)

```
src/analysis.py                          (添加 design_iir_correction facade,可选)
tests/regression/generate_baselines.py   (添加 predistortion_default 生成)
```

### 7.3 接口快照

| 符号 | 来源 |
|---|---|
| `ControlLine(name, kind, source, target, transfer_function, metadata)` | sqc.hardware.control_line |
| `SingleExponentialDistortion(amplitude, tau)` | sqc.hardware.distortion |
| `MultiExponentialDistortion(amplitudes, taus)` | 同上 |
| `FIRDistortion(taps, dt)` | 同上 |
| `IIRDistortion(b, a, dt)` | 同上 |
| `CustomTransferDistortion(omega, H)` | 同上 |
| `TransferFunctionCalibration(qubit, method, fit_type, n_exp_components)` | sqc.calibration.transfer_function |
| `PredistortionDesigner(method, n_taps, regularization)` | sqc.calibration.predistortion |
| `PredistortionDesigner.design(transfer_model, omega_grid) -> DistortionModel` | 同上 |
| `PredistortionDesigner.predistort(target, transfer) -> Waveform` | 同上 |
| `PredistortionValidationWorkflow(qubit, target_waveform, true_distortion, designer)` | sqc.workflows.predistortion_validation |
| `PredistortionValidationWorkflow.run() -> dict` | 同上 |

---

## 8. 完成确认清单

```
- [ ] sqc/hardware/control_line.py 完整实现
- [ ] sqc/hardware/distortion.py 5 个子类完整实现 (Single/MultiExp, FIR, IIR, Custom)
- [ ] sqc/calibration/transfer_function.py 完整实现
- [ ] sqc/calibration/predistortion.py 完整实现
- [ ] sqc/workflows/predistortion_validation.py 完整实现
- [ ] Track B 1.3 在 src/ 中的 DistortionModel 已迁移到 sqc/hardware/distortion.py
- [ ] tests/baselines/predistortion_default.pkl 生成
- [ ] tests/unit 新增 3 个文件,各 ≥ 3 个用例
- [ ] tests/integration 新增 2 个文件
- [ ] tests/regression 新增 1 个文件
- [ ] pytest tests/unit -v 全部 pass
- [ ] pytest tests/integration -v 全部 pass
- [ ] pytest tests/regression -m regression -v 全部 pass
- [ ] PredistortionValidationWorkflow 改善 factor > 10
- [ ] python web_demo.py 仍跑通(若加 Predistortion tab 则也跑通)
- [ ] Simulation.ipynb 加预失真 demo cell 并可运行
```

---

**下一步**:[phase_5_handbook.md](phase_5_handbook.md) — 多 qubit Z-crosstalk demo。
