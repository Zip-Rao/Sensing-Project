# Phase 3 Handbook — Track B 成果内化:LM/Cryoscope/瞬态频率标定

> 前置阅读:[_refactor_plan.md](_refactor_plan.md) §6.5–§6.7、§7.4–§7.5、[phase_2_handbook.md](phase_2_handbook.md)  
> 估计工时:5–7 天(取决于 Track B 完成度)  
> 触发条件:Track B 已完成 _TODO_master.md 阶段 0+1(LM 修复 + Cryoscope + 瞬态标定 + 失真模型)中**至少**:0.3 (LM 修复)、1.1 (Cryoscope)、1.2 (瞬态标定)
> 完成标志:`LMReconstruction` / `CryoscopeExperiment` / `TransientFrequencyCalibration` / `FluxResponseCalibration` / `CryoscopeReconstruction` 全部实现并通过对应 baseline

---

## 1. 目标

把 Track B 在 src/ 内已完成的功能,**内化**到 sqc/ 的对应位置。本 phase 是"功能性重构"的核心:从 src/ 把代码搬过来,清理设计债务,加测试。

涉及的功能模块:
1. **CryoscopeExperiment**(case 5/6/7 内化)
2. **TransientFrequencyCalibration**(case 8 内化)
3. **FluxResponseCalibration**(整合 Calibration case 1/2/3)
4. **LMReconstruction**(`Analysis.numerical_inverse` + `levenberg_marquardt` + Jacobian 全部内化)
5. **CryoscopeReconstruction**(`Analysis.get_signal_from_cryoscope` + `get_h_from_phi` 内化)
6. **RamseyIQReconstruction / RamseyUnwrapReconstruction / DiffEchoReconstruction**(`Analysis.get_signal_from_*` 类化)
7. **WienerReconstruction / HammersteinWienerReconstruction**(`Analysis.wiener_deconvolution` 等类化)
8. **basis 模块**(`generate_basis_functions` / `basis_function_decomposition` / `R` 内化)

---

## 2. 前置条件

| # | 条件 | 验证 |
|---|---|---|
| 2.1 | Phase 2 完成 | `pytest tests/regression` 全部 pass,Experiment 类已在用 |
| 2.2 | Track B `_TODO_master.md 0.3` (LM 收敛修复) 完成 | `Analysis.numerical_inverse` 在 src/ 内能稳定收敛;有相应 baseline |
| 2.3 | Track B `_TODO_master.md 1.1` (Cryoscope case 6/7) 完成 | `Protocal(type=5).evolve()` 行为稳定;`Calibration(type=3).calibrate()` 行为稳定 |
| 2.4 | Track B `_TODO_master.md 1.2` (瞬态标定 case 8) 完成 | 在 src/protocal.py 中有 case 8 实现;`Analysis.calibrate_frequency_response` 之类的方法已实现 |
| 2.5 | Track B 已重生成相关 baseline | `tests/baselines/` 包含 cryoscope_default.pkl、transient_calib_default.pkl、lm_default.pkl |

如果 Track B 任一前置项尚未完成,**该项内化推迟**。例如:LM 修复未完成,本 phase 仅做 Cryoscope/瞬态标定/Ramsey-IQ 等内化,LM 部分留到下一次 P3。

---

## 3. 任务清单

### 3.1 内化 basis 模块

**文件**:`sqc/reconstruction/basis.py`

把 `src/analysis.py:325-400` 的 `generate_basis_functions`、`basis_function_decomposition`、`R` 直接 port 过来。

```python
# sqc/reconstruction/basis.py
from __future__ import annotations
from typing import Literal

import numpy as np


BasisType = Literal["bspline", "fourier", "legendre"]


def generate_basis_functions(basis_type: BasisType, n_basis: int,
                              t_min: float, t_max: float) -> list:
    """Generate a list of basis functions on [t_min, t_max].
    
    Verbatim port from src/analysis.py:327.
    """
    # ... port from src/analysis.py ...


def basis_function_decomposition(signal_or_obj, t_array: np.ndarray,
                                  basis_functions: list) -> np.ndarray:
    """Decompose a signal onto basis functions via least squares.
    
    Verbatim port from src/analysis.py:361.
    """
    # ... port ...


def regularization_matrix(n: int, basis_type: BasisType) -> np.ndarray:
    """Smoothness regularization matrix for basis_type.
    
    Replaces R(n, basis_type) in src/analysis.py:386.
    """
    # ... port ...


# Backward-compat alias
R = regularization_matrix
```

`src/analysis.py` 镜像层重新导出:
```python
# (P3 之后 src/analysis.py 顶部新增)
from sqc.reconstruction.basis import (
    generate_basis_functions, basis_function_decomposition, R,
)
```

### 3.2 内化 Wiener / RamseyIQ / RamseyUnwrap / DiffEcho reconstruction

**文件**:`sqc/reconstruction/wiener.py`

把 `Analysis.wiener_deconvolution`、`get_signal_from_ramsey_by_iq`、`get_signal_from_ramsey_by_unwrap`、`get_signal_from_diff_echo` 类化。

```python
# sqc/reconstruction/wiener.py
from __future__ import annotations
from dataclasses import dataclass

import numpy as np

from sqc.devices.transmon import QubitSpec
from sqc.reconstruction.base import Reconstruction
from sqc.simulation.result import ExperimentResult
from sqc.control.flux_signal import FluxSignal


@dataclass
class WienerReconstruction(Reconstruction):
    lambda_reg: float = 1.0
    
    def reconstruct(self, measurement: ExperimentResult,
                    kernel: np.ndarray,
                    calibration=None,
                    dt: float | None = None) -> FluxSignal:
        """Linear Wiener deconvolution. Verbatim from Analysis.wiener_deconvolution."""
        delta_p = measurement.data["delta_p"]
        if dt is None:
            dt = measurement.axes["scan"][1] - measurement.axes["scan"][0]
        
        N_del = len(delta_p)
        N_ker = len(kernel)
        N = N_del - N_ker + 1
        N_fft = N_del
        p_pad = np.zeros(N_fft); p_pad[:N_del] = delta_p
        k_pad = np.zeros(N_fft); k_pad[:N_ker] = kernel
        Y = np.fft.fft(p_pad)
        H = np.fft.fft(k_pad)
        H_conj = np.conj(H)
        G = H_conj / (np.abs(H)**2 + self.lambda_reg**2)
        X_w = Y * G / dt
        x_rec = np.real(np.fft.ifft(X_w))
        x_final = x_rec[:N]
        x_lists = np.linspace(0, (N-1)*dt, N)
        return FluxSignal(type=8, t_list=x_lists, signal=x_final)


@dataclass
class RamseyIQReconstruction(Reconstruction):
    qubit: QubitSpec
    
    def reconstruct(self, measurement, kernel=None, calibration=None) -> np.ndarray:
        """B(τ) from Ramsey IQ data. Replaces Analysis.get_signal_from_ramsey_by_iq."""
        # ... port src/analysis.py:50-74 ...


@dataclass
class RamseyUnwrapReconstruction(Reconstruction):
    qubit: QubitSpec
    k_span: int = 3
    
    def reconstruct(self, measurement, kernel=None, calibration=None) -> np.ndarray:
        """B(τ) from Ramsey p_e via phase unwrapping. Replaces Analysis.get_signal_from_ramsey_by_unwrap."""
        # ... port src/analysis.py:75-129 ...


@dataclass
class DiffEchoReconstruction(Reconstruction):
    qubit: QubitSpec
    t_int: float
    k: int
    
    def reconstruct(self, measurement, kernel=None, calibration=None) -> np.ndarray:
        """B from differential echo p_e. Replaces Analysis.get_signal_from_diff_echo."""
        p_e_list = measurement.data["p_e"]
        varphi = np.arcsin(2 * np.array(p_e_list) - 1)
        kappa = self.qubit.sensitivity()
        return -varphi / (2 * self.k * kappa * self.t_int)
```

### 3.3 内化 HammersteinWienerReconstruction

**文件**:`sqc/reconstruction/hammerstein.py`

```python
# sqc/reconstruction/hammerstein.py
from __future__ import annotations
from dataclasses import dataclass

import numpy as np
from sqc.devices.transmon import QubitSpec
from sqc.reconstruction.base import Reconstruction
from sqc.reconstruction.wiener import WienerReconstruction


@dataclass
class HammersteinWienerReconstruction(Reconstruction):
    qubit: QubitSpec
    lambda_reg: float = 1.0
    
    def reconstruct(self, measurement, kernel, calibration=None,
                    dt=None) -> np.ndarray:
        """Block-structured B(t) reconstruction.
        
        Replaces Analysis.hammerstein_wiener_deconvolution.
        """
        # 1. Linear Wiener gives ω(t) = Δω(t)
        wiener = WienerReconstruction(lambda_reg=self.lambda_reg)
        omega_signal = wiener.reconstruct(measurement, kernel, dt=dt)
        omega = omega_signal.samples
        omega_lists = omega_signal.t_list
        
        # 2. Inverse Transmon dispersion: B = (1/π) arccos(...)
        EC, EJ = self.qubit.EC, self.qubit.EJ
        B = (1 / np.pi) * np.arccos(
            (omega + self.qubit.frequency() + EC) ** 2 / (8 * EC * EJ)
        )
        return B
```

### 3.4 内化 LMReconstruction (核心)

**文件**:`sqc/reconstruction/numerical_inverse.py`

把 `src/analysis.py:272-825` 全部相关代码 port 过来:
- `numerical_inverse` 方法 → `LMReconstruction.reconstruct`
- `forward_simulation` → 私有 `_forward_simulation`
- `compute_jacobian` (adjoint) → 私有 `_compute_jacobian_adjoint`
- `compute_jacobian_finite_difference` → 私有 `_compute_jacobian_fd`
- `levenberg_marquardt` → 私有 `_levenberg_marquardt`

```python
# sqc/reconstruction/numerical_inverse.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from sqc.devices.transmon import TransmonQubit
from sqc.control.pulse import CompositePulse
from sqc.control.flux_signal import FluxSignal
from sqc.reconstruction.base import Reconstruction
from sqc.reconstruction.basis import (
    generate_basis_functions, basis_function_decomposition,
    regularization_matrix,
)


@dataclass
class LMReconstruction(Reconstruction):
    """Levenberg-Marquardt full-density-matrix waveform reconstruction.
    
    Replaces Analysis.numerical_inverse + module-level helpers.
    """
    qubit: TransmonQubit
    control_pulse: CompositePulse
    basis_type: Literal["bspline", "fourier", "legendre"] = "fourier"
    n_basis: int = 100
    lambda_reg: float = 100.0
    max_iter: int = 10
    tol: float = 1e-6
    mu_init: float = 1e-3
    use_adjoint: bool = True            # if False, fall back to finite difference
    
    def reconstruct(self, measurement, kernel=None, calibration=None,
                    initial_guess: np.ndarray | None = None) -> tuple[FluxSignal, dict]:
        """Returns (B_opt: FluxSignal, history: dict).
        
        measurement.data must contain "p_meas" (1-D array of p_e).
        """
        p_meas = measurement.data["p_meas"]
        t_list = measurement.axes["t_signal"]
        
        B_init = FluxSignal(
            type=6, t_list=t_list,
            n_basis=self.n_basis, basis_type=self.basis_type,
        )
        if initial_guess is None:
            initial_guess = np.zeros_like(t_list)
        b_init = basis_function_decomposition(initial_guess, t_list,
                                               B_init.basis_functions)
        B_init.update_signal(b=b_init)
        
        b_opt, history = self._levenberg_marquardt(
            p_meas, t_list, b_init, B_init,
        )
        B_init.update_signal(b=b_opt)
        return B_init, history
    
    # --- private helpers (port from src/analysis.py:402-825) ---
    
    def _forward_simulation(self, B, t_meas, H_list, t_evolve_list): ...
    def _compute_jacobian_adjoint(self, B, t_meas, results, H_list,
                                   t_evolve_list): ...
    def _compute_jacobian_fd(self, B, t_meas, p_sim, H_list,
                              t_evolve_list, eps=1e-6): ...
    def _levenberg_marquardt(self, p_meas, t_list, b_init, B_init): ...
    def _build_h_for_qubit(self, t_list, t_meas): ...
```

**关键**:Track B 0.3 应该已经修复了 Jacobian 不一致问题,因此 P3 内化时应使用**修复后的版本**。在 commit message 中显式说明:"内化的是 Track B 0.3 修复后的版本"。

### 3.5 实现 CryoscopeExperiment

**文件**:`sqc/experiments/cryoscope.py`

把 `src/protocal.py:164-183` (case 5) port 过来,但**修复**:
- 当前 case 5 的 truncate 是 in-place 修改 Phi(因为 Signal.truncate 是 in-place);新版用 Waveform.truncate 返回新对象。
- 现 case 5 trunc_list 用 `[140:20:-1]` 倒序,需保留以匹配 baseline。

```python
# sqc/experiments/cryoscope.py
from __future__ import annotations
from dataclasses import dataclass, field

import numpy as np
from qutip import basis

from sqc.devices.transmon import TransmonQubit
from sqc.control.flux_signal import FluxSignal
from sqc.hardware.readout import IQReadoutModel
from sqc.simulation.result import ExperimentResult
from sqc.experiments.base import Experiment


@dataclass
class CryoscopeExperiment(Experiment):
    """Cryoscope: scan truncation delay, measure φ(t_d) via IQ Ramsey.
    
    Replaces Protocal.evolve case 5.
    """
    qubit: TransmonQubit
    flux_signal: FluxSignal | None = None
    t_rabi: np.ndarray = field(default_factory=lambda: np.linspace(0, 10, 20))
    tau: float = 100.0
    trunc_list: np.ndarray | None = None
    
    def __post_init__(self):
        if self.flux_signal is None:
            self.flux_signal = FluxSignal(
                type=2, t_list=np.linspace(0, 80, 160), amplitude=0.01,
            )
        if self.trunc_list is None:
            self.trunc_list = self.flux_signal.t_list[140:20:-1]
    
    def build_sequence(self):
        return None
    
    def run(self) -> ExperimentResult:
        readout = IQReadoutModel(tau=self.tau, t_rabi=self.t_rabi)
        p_e_I_list = []
        p_e_Q_list = []
        
        for trunc in self.trunc_list:
            phi_truncated = self.flux_signal.copy()
            # Use legacy truncate (in-place mutation; same semantics as src)
            phi_truncated.truncate(0, trunc)
            self.qubit.qubit_in_mag(phi_truncated, frame=1,
                                     omega_d=self.qubit.frequency)
            measured = readout.measure(self.qubit)
            p_e_I_list.append(measured["p_e_I"])
            p_e_Q_list.append(measured["p_e_Q"])
        
        p_e_I = np.asarray(p_e_I_list)
        p_e_Q = np.asarray(p_e_Q_list)
        # Reverse to match legacy time axis ordering
        varphi = np.arctan2(p_e_Q - 0.5, p_e_I - 0.5)[::-1]
        
        return ExperimentResult(
            data={"varphi": varphi, "p_e_I": p_e_I, "p_e_Q": p_e_Q},
            axes={"trunc": np.asarray(self.trunc_list)},
            metadata={"experiment": "CryoscopeExperiment"},
            config={"tau": self.tau},
        )
```

### 3.6 实现 CryoscopeReconstruction

**文件**:`sqc/reconstruction/cryoscope.py`

```python
# sqc/reconstruction/cryoscope.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.interpolate import interp1d

from sqc.calibration.base import CalibrationTable
from sqc.reconstruction.base import Reconstruction


@dataclass
class CryoscopeReconstruction(Reconstruction):
    """Cryoscope waveform reconstruction.
    
    Replaces Analysis.get_signal_from_cryoscope + get_h_from_phi.
    """
    calibration: CalibrationTable     # phi(h) calibration
    tau: float                         # calibration square-pulse length
    method: Literal["calib_inverse", "SG_diff", "diff"] = "calib_inverse"
    
    def reconstruct(self, measurement, kernel=None, calibration=None,
                    dt: float | None = None) -> np.ndarray:
        varphi = measurement.data["varphi"]
        trunc_list = measurement.axes["trunc"]
        if dt is None:
            dt = abs(trunc_list[0] - trunc_list[1])
        
        if self.method == "calib_inverse":
            delta_phi = np.diff(varphi, prepend=0)
            delta_phi_norm = delta_phi * self.tau / dt
            h_meas = self.calibration.inverse(delta_phi_norm)
            return h_meas
        elif self.method == "SG_diff":
            # Savitzky-Golay differentiation followed by f^{-1}
            from scipy.signal import savgol_filter
            d_varphi = savgol_filter(varphi, window_length=11,
                                      polyorder=3, deriv=1, delta=dt)
            h_meas = self.calibration.inverse(d_varphi)
            return h_meas
        elif self.method == "diff":
            d_varphi = np.gradient(varphi, dt)
            h_meas = self.calibration.inverse(d_varphi)
            return h_meas
        else:
            raise ValueError(f"unknown method: {self.method}")
```

### 3.7 实现 FluxResponseCalibration (整合)

**文件**:`sqc/calibration/flux_response.py`

把现 `Calibration.calibrate` 的 case 1/2/3 整合成一个类,通过 `method` 参数选 backend。

```python
# sqc/calibration/flux_response.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.interpolate import interp1d

from sqc.devices.transmon import TransmonQubit
from sqc.calibration.base import Calibration, CalibrationTable
from sqc.control.flux_signal import FluxSignal
from sqc.experiments.cryoscope import CryoscopeExperiment
from sqc.experiments.transient import TransientSensingExperiment
from sqc.experiments.ramsey import RamseyExperiment


@dataclass
class FluxResponseCalibration(Calibration):
    """Calibrate f(Φ), κ, etc. via Ramsey, Cryoscope, or transient sensing.
    
    Replaces Calibration.calibrate case 1/2/3 + the case-5/6 cryoscope flow.
    """
    qubit: TransmonQubit
    method: Literal["ramsey", "cryoscope", "transient"] = "cryoscope"
    h_list: np.ndarray | None = None
    tau: float = 100.0
    
    def __post_init__(self):
        if self.h_list is None:
            self.h_list = np.linspace(-0.03, 0.03, 21)
    
    def calibrate(self) -> CalibrationTable:
        match self.method:
            case "cryoscope":
                return self._calibrate_cryoscope()
            case "transient":
                return self._calibrate_transient()
            case "ramsey":
                return self._calibrate_ramsey()
    
    def _calibrate_cryoscope(self) -> CalibrationTable:
        """Scan square-pulse height h, build φ(h) lookup table.
        
        Direct port of Calibration.calibrate case 3.
        """
        from sqc.hardware.readout import IQReadoutModel
        readout = IQReadoutModel(tau=self.tau)
        varphi_list = []
        
        def _make_signal(h):
            t_list = np.linspace(0, self.tau + 20, 240)
            samples = np.zeros_like(t_list)
            samples[(t_list >= 10) & (t_list <= self.tau + 10)] = h
            return FluxSignal(type=8, t_list=t_list, signal=samples)
        
        Phi = _make_signal(0.0)
        for h in self.h_list:
            Phi = _make_signal(h)
            self.qubit.qubit_in_mag(Phi, frame=1, omega_d=self.qubit.frequency)
            measured = readout.measure(self.qubit)
            varphi = np.arctan2(measured["p_e_Q"] - 0.5,
                                 measured["p_e_I"] - 0.5)
            varphi_list.append(varphi)
        
        varphi_list = np.unwrap(np.asarray(varphi_list), period=np.pi)
        
        return CalibrationTable(
            qubit_name=self.qubit.spec().name,
            kind="phi_h",
            inputs=np.asarray(self.h_list),
            outputs=varphi_list,
            fit_params={"tau": self.tau, "method": "cryoscope"},
            metadata={},
        )
    
    def _calibrate_transient(self) -> CalibrationTable:
        """Scan with known ramp/sine signal, fit polynomial coefficients
        of Δω(Φ).
        
        Direct port of (Track B) Calibration case 2 / Analysis.calibrate_frequency_response.
        """
        # Implementation per Track B's _frequency_calibration_via_transient.md
        ...
    
    def _calibrate_ramsey(self) -> CalibrationTable:
        """Scan flux Φ in static steps, fit Ramsey detuning at each step.
        
        Replaces Calibration case 0/1.
        """
        ...
```

#### 3.7.1 CalibrationTable 反函数实现

补全 `sqc/calibration/base.py` 中 `CalibrationTable.evaluate` 和 `inverse`:

```python
# sqc/calibration/base.py (补全)
@dataclass(frozen=True)
class CalibrationTable:
    qubit_name: str
    kind: Literal[...]
    inputs: np.ndarray
    outputs: np.ndarray
    fit_params: dict
    metadata: dict
    
    def evaluate(self, x: np.ndarray) -> np.ndarray:
        from scipy.interpolate import interp1d
        f = interp1d(self.inputs, self.outputs, kind="cubic",
                     fill_value="extrapolate")
        return f(x)
    
    def inverse(self, y: np.ndarray) -> np.ndarray:
        """Build inverse via sorting (assumes monotonic, else takes monotonic chunk)."""
        idx = np.argsort(self.outputs)
        x_sorted = self.inputs[idx]
        y_sorted = self.outputs[idx]
        mask = np.diff(y_sorted, prepend=-np.inf) > 1e-12
        from scipy.interpolate import interp1d
        f_inv = interp1d(y_sorted[mask], x_sorted[mask], kind="cubic",
                         fill_value="extrapolate")
        return f_inv(y)
```

### 3.8 实现 TransientFrequencyCalibration (case 8 内化)

**文件**:`sqc/calibration/qubit_frequency.py`

```python
# sqc/calibration/qubit_frequency.py
from __future__ import annotations
from dataclasses import dataclass

import numpy as np

from sqc.devices.transmon import TransmonQubit
from sqc.calibration.base import Calibration, CalibrationTable


@dataclass
class QubitFrequencyCalibration(Calibration):
    """Calibrate f₀₁ via Ramsey free precession.
    
    Replaces Calibration case 0.
    """
    qubit: TransmonQubit
    
    def calibrate(self) -> CalibrationTable: ...


@dataclass
class TransientFrequencyCalibration(Calibration):
    """Calibrate Δω(Φ) polynomial coefficients via transient sensing.
    
    Replaces Track B's case 8 implementation.
    """
    qubit: TransmonQubit
    polynomial_order: int = 3
    test_signal_kind: Literal["ramp", "sine", "gaussian"] = "ramp"
    
    def calibrate(self) -> CalibrationTable:
        # Port from Track B's analysis.calibrate_frequency_response
        ...
```

### 3.9 内化 RamseyIQ Reconstruction 与 Analysis 类的 facade

**文件修改**:`src/analysis.py`

把 `Analysis` 类改为 facade:

```python
# src/analysis.py (P3 之后)
"""src.analysis — facade over sqc.reconstruction.

Internally delegates to sqc/reconstruction/* classes. Legacy method names
preserved to keep Simulation.ipynb and other historical code working.
"""
from __future__ import annotations

import numpy as np

from sqc.reconstruction.basis import (
    generate_basis_functions, basis_function_decomposition,
    regularization_matrix as R,
)
from sqc.reconstruction.kernel import KernelEstimator
from sqc.reconstruction.wiener import (
    WienerReconstruction, RamseyIQReconstruction,
    RamseyUnwrapReconstruction, DiffEchoReconstruction,
)
from sqc.reconstruction.hammerstein import HammersteinWienerReconstruction
from sqc.reconstruction.numerical_inverse import LMReconstruction
from sqc.reconstruction.cryoscope import CryoscopeReconstruction


class Analysis:
    """Legacy facade over sqc.reconstruction.* classes."""
    
    def get_expectation_values(self, result, e_ops_index):
        from sqc.simulation.result import extract_expectation
        return extract_expectation(result, e_ops_index)
    
    def get_population(self, result, level):
        from sqc.simulation.result import extract_population
        return extract_population(result, level)
    
    def get_signal_from_ramsey_by_iq(self, qubit, tau_list,
                                      p_e_list_I, p_e_list_Q):
        from sqc.simulation.result import ExperimentResult
        spec = qubit.spec() if hasattr(qubit, "spec") else qubit
        recon = RamseyIQReconstruction(qubit=spec)
        meas = ExperimentResult(
            data={"p_e_I": np.asarray(p_e_list_I),
                  "p_e_Q": np.asarray(p_e_list_Q)},
            axes={"tau": np.asarray(tau_list)},
        )
        return recon.reconstruct(meas)
    
    def get_signal_from_ramsey_by_unwrap(self, qubit, tau_list, p_e_list,
                                          k_span=3):
        spec = qubit.spec() if hasattr(qubit, "spec") else qubit
        recon = RamseyUnwrapReconstruction(qubit=spec, k_span=k_span)
        # ... build ExperimentResult and delegate ...
        ...
    
    def get_signal_from_diff_echo(self, qubit, p_e_list, t_int, k):
        spec = qubit.spec() if hasattr(qubit, "spec") else qubit
        recon = DiffEchoReconstruction(qubit=spec, t_int=t_int, k=k)
        ...
    
    def get_kernel(self, control_pulse, qubit):
        return KernelEstimator().estimate(control_pulse, qubit)
    
    def wiener_deconvolution(self, delta_p, kernel, dt, lambdas):
        recon = WienerReconstruction(lambda_reg=lambdas)
        ...
    
    def hammerstein_wiener_deconvolution(self, qubit, delta_p, kernel,
                                          dt, lambdas):
        spec = qubit.spec() if hasattr(qubit, "spec") else qubit
        recon = HammersteinWienerReconstruction(qubit=spec, lambda_reg=lambdas)
        ...
    
    def numerical_inverse(self, qubit, control_pulse, p_meas, t_list,
                          B_guess, basis_type='fourier', n_basis=100,
                          lambdas=100.0, max_iter=10, tol=1e-6):
        recon = LMReconstruction(
            qubit=qubit, control_pulse=control_pulse,
            basis_type=basis_type, n_basis=n_basis,
            lambda_reg=lambdas, max_iter=max_iter, tol=tol,
        )
        ...
    
    def get_h_from_phi(self, h_list, phi_list):
        # Build a CalibrationTable and return interpolation functions
        ...
    
    def get_signal_from_cryoscope(self, qubit, trunc_list, varphi_meas,
                                   h_of_phi, tau, dt):
        # Reconstruct via CryoscopeReconstruction
        ...


# Module-level helpers (port verbatim, but they live in sqc.reconstruction.basis now)
# Here we just re-export them for backward compat:
__all__ = [
    "Analysis",
    "generate_basis_functions",
    "basis_function_decomposition",
    "R",
    # forward_simulation, compute_jacobian, etc. moved into LMReconstruction
    # but we expose them as module-level functions for legacy support:
    "forward_simulation",
    "compute_jacobian",
    "compute_jacobian_finite_difference",
    "levenberg_marquardt",
]


def forward_simulation(*args, **kwargs):
    from sqc.reconstruction.numerical_inverse import _forward_simulation
    return _forward_simulation(*args, **kwargs)


def compute_jacobian(*args, **kwargs):
    from sqc.reconstruction.numerical_inverse import _compute_jacobian_adjoint
    return _compute_jacobian_adjoint(*args, **kwargs)


def compute_jacobian_finite_difference(*args, **kwargs):
    from sqc.reconstruction.numerical_inverse import _compute_jacobian_fd
    return _compute_jacobian_fd(*args, **kwargs)


def levenberg_marquardt(*args, **kwargs):
    from sqc.reconstruction.numerical_inverse import _levenberg_marquardt
    return _levenberg_marquardt(*args, **kwargs)
```

### 3.10 更新 src/protocal.py:Calibration 为 facade

```python
# src/protocal.py 中的 Calibration 类 (P3)
class Calibration:
    """Legacy facade over sqc.calibration.* classes."""
    
    def __init__(self, qubit, type=0, **kwargs):
        self.qubit = qubit
        self.type = type
        self.params = kwargs
    
    def calibrate(self):
        from sqc.calibration.flux_response import FluxResponseCalibration
        from sqc.calibration.qubit_frequency import (
            QubitFrequencyCalibration, TransientFrequencyCalibration,
        )
        match self.type:
            case 0:
                cal = QubitFrequencyCalibration(qubit=self.qubit)
                return cal.calibrate()
            case 1:
                cal = FluxResponseCalibration(qubit=self.qubit, method="ramsey")
                return cal.calibrate()
            case 2:
                cal = TransientFrequencyCalibration(qubit=self.qubit, **self.params)
                return cal.calibrate()
            case 3:
                cal = FluxResponseCalibration(qubit=self.qubit, method="cryoscope",
                                               **self.params)
                table = cal.calibrate()
                return table.inputs, table.outputs, cal.tau   # legacy tuple
            case _:
                raise ValueError(f"Unknown calibration type {self.type}")
```

---

## 4. 验收标准

### 4.1 各内化项独立验证

```bash
# Cryoscope
python -c "
import numpy as np
from sqc.devices.transmon import TransmonQubit
from sqc.experiments.cryoscope import CryoscopeExperiment
q = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000)
exp = CryoscopeExperiment(qubit=q)
result = exp.run()
print('Cryoscope OK, varphi shape:', result.data['varphi'].shape)
"

# LM reconstruction
python -c "
... (类似的 smoke test)
"

# FluxResponseCalibration
python -c "
... (类似)
"
```

### 4.2 物理回归(关键)

```bash
pytest tests/regression -m regression -v
```

**全部 pass**,包括 P3 引入的新 baseline:
- `cryoscope_default.pkl`(case 5 输出)
- `lm_default.pkl`(LM 反演输出)
- `transient_calib_default.pkl`(case 8 输出)

### 4.3 facade 等价性

```bash
python -c "
import numpy as np
from src.qubit import TransmonQubit
from src.protocal import Protocal, Calibration
from src.analysis import Analysis

q = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000)

# Cryoscope via legacy
proto = Protocal(type=5)
proto.initialize(q)
trunc_list, varphi, Phi, p_e_list = proto.evolve(q)
print('legacy Cryoscope OK')

# Analysis via legacy
ana = Analysis()
# ... call various methods, verify outputs
"
```

---

## 5. 测试要求

### 5.1 必须新增的测试

| 文件 | 用例 |
|---|---|
| `tests/integration/test_cryoscope_experiment.py` | CryoscopeExperiment 与 Track B 完成的 src/protocal.py case 5/6/7 输出一致 |
| `tests/integration/test_lm_reconstruction.py` | LMReconstruction.reconstruct 与 Track B 修复后的 Analysis.numerical_inverse 输出一致 |
| `tests/integration/test_flux_response_cal.py` | 三个 method (ramsey/cryoscope/transient) 各自输出 CalibrationTable,数值与 Track B 实现一致 |
| `tests/unit/test_basis_module.py` | generate_basis_functions、basis_function_decomposition、regularization_matrix 与 src/analysis.py 模块级函数等价 |
| `tests/unit/test_calibration_table.py` | CalibrationTable.evaluate / inverse 行为正确;inverse 失败时给清晰错误 |
| `tests/integration/test_legacy_analysis_facade.py` | Analysis 类的所有公共方法返回与 P2 之前版本数值等价 |

### 5.2 必须保持的回归

- 所有 P0 baseline + P3 新增 baseline 全部 pass

---

## 6. 风险与回滚

| # | 风险 | 缓解 |
|---|---|---|
| 6.1 | LM 内化时把 Track B 0.3 修复未完整 port,导致 LM 重新出现收敛问题 | 内化前在 src/ 跑一次完整 LM 反演,记录 b_opt 与 history 作为 baseline;内化后比对 |
| 6.2 | Cryoscope truncate 在新 Waveform 体系下行为不同(返回新对象 vs in-place) | CryoscopeExperiment 内显式调用 `.copy()` 后 `.truncate()`,保留 in-place 语义 |
| 6.3 | CalibrationTable.inverse 在非单调段失败 | 加 mask 取单调子段,失败时抛 ValueError 而非返回错误结果 |
| 6.4 | facade 转发后旧 Notebook 因为 import 路径间接,启动慢 | 顶层 import 重构为 lazy import |
| 6.5 | LM 在新结构里依赖 qubit.qubit_in_mag 副作用 | LMReconstruction 显式调用 `qubit.qubit_in_mag(B_curr)` 在每个 trial 之前,与旧实现一致 |

**回滚**:
- 若 LM 内化破坏数值行为,**立即 revert**,LM 留在 src/ 内,FluxResponse/Cryoscope 内化继续。
- 若 Cryoscope 内化破坏 web_demo,**立即 revert**,在主方案 §13 假设 A5 上更新优先级。

---

## 7. 输出物

### 7.1 文件清单(新增)

```
sqc/reconstruction/basis.py
sqc/reconstruction/wiener.py                    (P3 完整实现)
sqc/reconstruction/hammerstein.py
sqc/reconstruction/numerical_inverse.py
sqc/reconstruction/cryoscope.py
sqc/calibration/qubit_frequency.py
sqc/calibration/flux_response.py
sqc/experiments/cryoscope.py
sqc/experiments/cpmg.py                          (若 Track B CPMG 完成)
tests/unit/test_basis_module.py
tests/unit/test_calibration_table.py
tests/integration/test_cryoscope_experiment.py
tests/integration/test_lm_reconstruction.py
tests/integration/test_flux_response_cal.py
tests/integration/test_legacy_analysis_facade.py
tests/baselines/cryoscope_default.pkl            (新增)
tests/baselines/lm_default.pkl                   (新增)
tests/baselines/transient_calib_default.pkl     (新增)
```

### 7.2 文件清单(修改)

```
src/analysis.py                                  (Analysis 改为 facade)
src/protocal.py                                  (Calibration 改为 facade,case 5 改为 facade)
sqc/calibration/base.py                          (CalibrationTable.evaluate / inverse 实现)
tests/regression/generate_baselines.py           (加 cryoscope/lm/transient_calib 三个 baseline)
```

### 7.3 接口快照

| 符号 | 来源 |
|---|---|
| `WienerReconstruction(lambda_reg)` | sqc.reconstruction.wiener |
| `RamseyIQReconstruction(qubit)` | 同上 |
| `RamseyUnwrapReconstruction(qubit, k_span=3)` | 同上 |
| `DiffEchoReconstruction(qubit, t_int, k)` | 同上 |
| `HammersteinWienerReconstruction(qubit, lambda_reg)` | sqc.reconstruction.hammerstein |
| `LMReconstruction(qubit, control_pulse, basis_type, n_basis, lambda_reg, max_iter, tol, mu_init, use_adjoint)` | sqc.reconstruction.numerical_inverse |
| `CryoscopeReconstruction(calibration, tau, method)` | sqc.reconstruction.cryoscope |
| `CryoscopeExperiment(qubit, flux_signal, t_rabi, tau, trunc_list)` | sqc.experiments.cryoscope |
| `FluxResponseCalibration(qubit, method, h_list, tau)` | sqc.calibration.flux_response |
| `TransientFrequencyCalibration(qubit, polynomial_order, test_signal_kind)` | sqc.calibration.qubit_frequency |
| `QubitFrequencyCalibration(qubit)` | 同上 |
| `CalibrationTable.evaluate(x) / inverse(y)` | sqc.calibration.base |

---

## 8. 完成确认清单

```
- [ ] sqc/reconstruction/basis.py 完整实现
- [ ] sqc/reconstruction/{wiener,hammerstein,numerical_inverse,cryoscope}.py 完整实现
- [ ] sqc/calibration/{qubit_frequency,flux_response}.py 完整实现
- [ ] sqc/calibration/base.py:CalibrationTable.evaluate/inverse 实现
- [ ] sqc/experiments/cryoscope.py 完整实现
- [ ] src/analysis.py:Analysis 改为 facade
- [ ] src/protocal.py:Calibration 改为 facade
- [ ] src/protocal.py case 5 改为转发 CryoscopeExperiment
- [ ] tests/baselines/ 新增 3 个 pkl
- [ ] tests/regression/generate_baselines.py 包含 cryoscope/lm/transient_calib
- [ ] tests/integration/ 至少新增 4 个文件
- [ ] pytest tests/unit -v 全部 pass
- [ ] pytest tests/integration -v 全部 pass
- [ ] pytest tests/regression -m regression -v 全部 pass(关键!)
- [ ] python web_demo.py 仍然跑通(若有 Cryoscope tab 也跑通)
- [ ] Simulation.ipynb 前 25 个 cell 可运行,包括 Cryoscope/LM 反演相关 cell
```

---

**下一步**:[phase_4_handbook.md](phase_4_handbook.md) — ControlLine + DistortionModel + 预失真。
