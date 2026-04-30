# Phase 1 Handbook — sqc/ 骨架 + ABC + 数据结构 + src/ 镜像

> 前置阅读:[_refactor_plan.md](_refactor_plan.md) §4–§8、[phase_0_handbook.md](phase_0_handbook.md)  
> 估计工时:3–5 天  
> 触发条件:Phase 0 已完成,baseline 测试全部通过  
> 完成标志:`sqc/` 完整目录建立,所有 ABC 定义,所有数据结构定义,`src/` 镜像不破坏任何 baseline

---

## 1. 目标

完成全栈重构的"骨架阶段"。本 phase 不实现任何业务逻辑,只:

1. 创建 `sqc/` 完整目录树。
2. 在所有应有 ABC 的位置写好抽象基类,实现处 `raise NotImplementedError`。
3. 定义所有数据结构(`@dataclass`)。
4. 把当前 `src/qubit.py` 中的 `TransmonQubit` 等已实现的物理对象**搬到** `sqc/devices/transmon.py`,并修复 D1 债务(qubit_in_mag 副作用)。
5. 实现 `src/` 永久镜像层,使 `from src.qubit import TransmonQubit` 完全等价于 `from sqc.devices.transmon import TransmonQubit`。
6. 物理回归测试 100% 通过。

**本 phase 完成后**,重构有了承载未来工作的"地基",但实验逻辑仍主要在 `src/protocal.py` 中。Phase 2 会迁移已实现的实验逻辑;Phase 3+ 会迁移 Track B 完成的功能。

---

## 2. 前置条件

| # | 条件 | 验证方式 |
|---|---|---|
| 2.1 | Phase 0 完成 | `pytest tests/regression -m regression` 全部 pass |
| 2.2 | git 工作区干净 | `git status` |
| 2.3 | 当前 baseline 已锁定在 git | `git ls-files tests/baselines/` 列出所有 pkl |

---

## 3. 任务清单

任务编号即建议执行顺序。粗体任务是关键路径,不可跳过。

### 3.1 创建 sqc/ 完整目录骨架

**任务**:在项目根创建以下目录与空 `__init__.py` 文件。

```
sqc/__init__.py
sqc/devices/__init__.py
sqc/hardware/__init__.py
sqc/control/__init__.py
sqc/simulation/__init__.py
sqc/experiments/__init__.py
sqc/calibration/__init__.py
sqc/reconstruction/__init__.py
sqc/workflows/__init__.py
```

`sqc/__init__.py` 内容:

```python
"""sqc — Superconducting Quantum Control simulation framework.

See idea/refactor/_refactor_plan.md for full architecture.
"""
__version__ = "0.1.0"
```

每个子包的 `__init__.py` 暂时为空(P5 时酌情暴露顶层快捷 import,见主方案 §13.2 Q1)。

**验收**:`python -c "import sqc; print(sqc.__version__)"` 输出 `0.1.0`。

### 3.2 定义所有数据结构

按主方案 §5 的定义实现。**每个 dataclass 必须有 docstring 和类型注解**。

#### 3.2.1 `sqc/devices/transmon.py` 中的 QubitSpec

参考主方案 §5.1。**关键约束**:`QubitSpec` 是 `@dataclass(frozen=True)`,`frequency()`、`anharmonicity()`、`sensitivity()` 是方法而非属性,且不修改对象。

```python
# sqc/devices/transmon.py
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from qutip import basis, destroy, num, qeye, Qobj


@dataclass(frozen=True)
class QubitSpec:
    """Pure parameter description of a flux-tunable Transmon qubit.
    
    Holds NO experimental state. Methods are pure: they take parameters
    and return values without mutating self.
    """
    name: str
    EC: float                 # rad·GHz
    EJ: float                 # rad·GHz at zero flux
    T1: float                 # ns
    T2: float                 # ns
    flux_bias: float = 0.0    # Φ₀
    n_levels: int = 3
    
    def EJ_at(self, flux: float | None = None) -> float:
        f = self.flux_bias if flux is None else flux
        return self.EJ * abs(math.cos(math.pi * f))
    
    def frequency(self, flux: float | None = None) -> float:
        """f₀₁(Φ) = √(8 EJ(Φ) EC) - EC, in rad·GHz."""
        return float(np.sqrt(8 * self.EJ_at(flux) * self.EC) - self.EC)
    
    def anharmonicity(self) -> float:
        return -self.EC
    
    def sensitivity(self, flux: float | None = None,
                    delta: float = 1e-6) -> float:
        f = self.flux_bias if flux is None else flux
        f_plus = self.frequency(f + delta)
        f_minus = self.frequency(f - delta)
        return (f_plus - f_minus) / (2 * delta)
    
    def with_flux(self, new_flux: float) -> "QubitSpec":
        """Return a new QubitSpec with flux_bias=new_flux. Does not mutate self."""
        return QubitSpec(
            name=self.name, EC=self.EC, EJ=self.EJ,
            T1=self.T1, T2=self.T2,
            flux_bias=new_flux, n_levels=self.n_levels,
        )
    
    @staticmethod
    def optimal_work_point() -> float:
        """Φ at which |df/dΦ| is maximum, in units of Φ₀."""
        return float(np.arctan(np.sqrt(2)) / np.pi)
```

#### 3.2.2 `sqc/control/waveform.py` 中的 Waveform 与 CompositeWaveform

参考主方案 §5.2。`Waveform` 是 mutable dataclass(因为 truncate 等会修改 samples)。

```python
# sqc/control/waveform.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any
import numpy as np


@dataclass
class Waveform:
    """Generic time-domain waveform, semantically neutral.
    
    For physical interpretation (e.g., flux signal), use FluxSignal.
    """
    t_list: np.ndarray
    samples: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        self.t_list = np.asarray(self.t_list, dtype=float)
        self.samples = np.asarray(self.samples, dtype=float)
        if self.t_list.shape != self.samples.shape:
            raise ValueError(
                f"shape mismatch: t_list {self.t_list.shape} vs "
                f"samples {self.samples.shape}"
            )
    
    @property
    def duration(self) -> float:
        return float(self.t_list[-1] - self.t_list[0])
    
    @property
    def n_points(self) -> int:
        return len(self.t_list)
    
    def value_at(self, t: float) -> float:
        """Sample-and-hold at time t. Returns 0 if t out of range."""
        if t < self.t_list[0] or t > self.t_list[-1]:
            return 0.0
        idx = int(np.argmin(np.abs(self.t_list - t)))
        return float(self.samples[idx])
    
    def truncate(self, t_start: float, t_end: float) -> "Waveform":
        """Return a NEW waveform with samples zeroed outside [t_start, t_end]."""
        mask = (self.t_list >= t_start) & (self.t_list <= t_end)
        new_samples = np.where(mask, self.samples, 0.0)
        return type(self)(
            t_list=self.t_list.copy(),
            samples=new_samples,
            metadata=dict(self.metadata),
        )
    
    def copy(self) -> "Waveform":
        return type(self)(
            t_list=self.t_list.copy(),
            samples=self.samples.copy(),
            metadata=dict(self.metadata),
        )
    
    def plot(self, ax=None, **kwargs):
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 4))
        ax.plot(self.t_list, self.samples, **kwargs)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Amplitude")
        ax.grid(True)
        return ax


@dataclass
class CompositeWaveform(Waveform):
    """Concatenation of multiple Waveforms in time."""
    components: list[Waveform] = field(default_factory=list)
    
    @classmethod
    def from_components(cls, components: list[Waveform]) -> "CompositeWaveform":
        t_list = []
        samples = []
        offset = 0.0
        for w in components:
            t_list.extend(t + offset for t in w.t_list)
            samples.extend(w.samples)
            if len(w.t_list) > 0:
                offset = t_list[-1] + 1e-9   # tiny gap to avoid duplicates
        return cls(
            t_list=np.array(t_list),
            samples=np.array(samples),
            components=list(components),
        )
```

#### 3.2.3 `sqc/control/flux_signal.py` 中的 FluxSignal

```python
# sqc/control/flux_signal.py
from __future__ import annotations
from .waveform import Waveform


class FluxSignal(Waveform):
    """Flux signal, samples in units of Φ₀."""
    pass


# Backward-compat alias for src/signal.py:Signal
Signal = FluxSignal
```

注:`Signal` 别名仅供 `src/signal.py` 镜像使用。新代码应直接用 `FluxSignal`。

#### 3.2.4 其余数据结构

按主方案 §5.3–§5.7 实现:

| 数据结构 | 文件 |
|---|---|
| `ControlLine` | `sqc/hardware/control_line.py` (P4 完整实现,P1 仅 dataclass 骨架) |
| `PulseSequence` | `sqc/control/sequence.py`(P1 仅 dataclass) |
| `ExperimentResult` | `sqc/simulation/result.py` |
| `MeasurementTrace` | `sqc/simulation/result.py` |
| `TransferMatrix` | `sqc/hardware/transfer_matrix.py` (P5 完整实现,P1 仅 dataclass 骨架) |
| `CalibrationTable` | `sqc/calibration/base.py` |

P1 阶段 ControlLine、TransferMatrix 只需要把 dataclass 字段定义出来,方法体 `raise NotImplementedError("implemented in PhaseN")`。

### 3.3 定义所有 ABC

每个 ABC **必须**:
1. 有 docstring 说明用途和子类必须实现的不变量。
2. 标注 `@abstractmethod`。
3. 不导入下层不应该依赖的模块(防止循环导入)。

#### 3.3.1 `sqc/devices/base.py` Device ABC

```python
# sqc/devices/base.py
from __future__ import annotations
from abc import ABC, abstractmethod
from qutip import Qobj


class Device(ABC):
    """Base class for all physical devices (qubits, resonators, couplers).
    
    Subclasses MUST be parameter-only objects: they describe physics but
    do not store experimental state. To represent a device under perturbation
    (e.g., qubit at non-zero flux), use a method that returns a NEW Device.
    """
    
    @property
    @abstractmethod
    def name(self) -> str: ...
    
    @abstractmethod
    def hilbert_dim(self) -> int: ...
    
    @abstractmethod
    def hamiltonian_static(self) -> Qobj:
        """Time-independent Hamiltonian without driving or external fields."""
    
    @abstractmethod
    def collapse_operators(self) -> list[Qobj]: ...
```

#### 3.3.2 其余 ABC

按主方案 §6 写出。每一处都遵循同一模板:

```python
# sqc/<layer>/base.py
from __future__ import annotations
from abc import ABC, abstractmethod

class XxxBase(ABC):
    """One-line summary.
    
    Detailed contract:
    - Invariant 1: ...
    - Invariant 2: ...
    """
    
    @abstractmethod
    def method_one(self, ...) -> ReturnType: ...
```

具体清单:

| 文件 | ABC | 关键方法 |
|---|---|---|
| `sqc/hardware/distortion.py` | `DistortionModel` | apply, step_response, impulse_response, frequency_response |
| `sqc/hardware/readout.py` | `ReadoutModel` | measure |
| `sqc/control/pulse.py` | `PulseBase` | hamiltonian, t_list, frame |
| `sqc/simulation/runner.py` | `RunnerBase`, `MesolveRunner`, `SlidingMeasurementRunner` | run |
| `sqc/experiments/base.py` | `Experiment` | build_sequence, run |
| `sqc/calibration/base.py` | `Calibration` | calibrate |
| `sqc/reconstruction/base.py` | `Reconstruction` | reconstruct |
| `sqc/workflows/base.py` | `Workflow` | run |

P1 阶段所有 abstract 方法在子类中暂时 `raise NotImplementedError`,具体实现留给 P2–P5。

### 3.4 迁移 TransmonQubit 到 sqc/devices/transmon.py

这是 P1 唯一的实质迁移工作。必须**保持物理行为完全等价**(回归测试通过),同时**消除 D1 债务**(qubit_in_mag 副作用)。

#### 3.4.1 新版 TransmonQubit 设计

**关键变化**:
- 新 `TransmonQubit` 仍然是一个有 mutable state 的类(为兼容 src/ 的 import),但**所有 `qubit_in_mag` 等会修改自身的方法**,内部调用 `HamiltonianBuilder.build()` 并显式存储结果到自身,**保留旧行为**。
- 同时暴露 `QubitSpec`,新代码应优先使用 `QubitSpec`。
- `TransmonQubit` 与 `QubitSpec` 之间的桥接通过 `TransmonQubit.spec()` 方法。

```python
# sqc/devices/transmon.py (续 §3.2.1)

class TransmonQubit(Device):
    """Backward-compatible Transmon qubit class.
    
    Wraps QubitSpec and maintains mutable state for legacy callers
    (src/protocal.py, Simulation.ipynb). New code should use QubitSpec
    + HamiltonianBuilder directly.
    """
    
    def __init__(self, EC, EJ, T1, T2, flux=0.0, state=0, n_levels=3,
                 name="Q"):
        self._spec = QubitSpec(
            name=name, EC=EC, EJ=EJ, T1=T1, T2=T2,
            flux_bias=flux, n_levels=n_levels,
        )
        # Legacy attributes (preserved verbatim for backward compatibility)
        self.EC = EC
        self.EJ = EJ * abs(math.cos(math.pi * flux))   # at current flux
        self.EJ_0 = EJ
        self.flux = flux
        self.n_levels = n_levels
        self.T1 = T1
        self.T2 = T2
        self.frequency = self._spec.frequency()
        self.anharmonicity = self._spec.anharmonicity()
        self.a = destroy(n_levels)
        self.a_dag = self.a.dag()
        self.n = num(n_levels)
        self.I = qeye(n_levels)
        self.hamiltonian = self.get_hamiltonian()
        self.c_ops = self.get_collapse_operators()
        # Init state
        if isinstance(state, int) and 0 <= state < n_levels:
            self.state = basis(n_levels, state)
        elif isinstance(state, Qobj) and state.dims == [[n_levels], [1]]:
            self.state = state.unit()
        else:
            self.state = basis(n_levels, 0)
    
    # --- new accessor ---
    def spec(self) -> QubitSpec:
        return self._spec
    
    # --- Device ABC implementation ---
    @property
    def name(self) -> str:
        return self._spec.name
    
    def hilbert_dim(self) -> int:
        return self.n_levels
    
    def hamiltonian_static(self) -> Qobj:
        return self.hamiltonian
    
    def collapse_operators(self) -> list[Qobj]:
        return self.c_ops
    
    # --- Legacy methods (verbatim port from src/qubit.py) ---
    
    def calculate_frequency(self) -> float:
        return self._spec.frequency()
    
    def calculate_anharmonicity(self) -> float:
        return self._spec.anharmonicity()
    
    def frequency_sensitivity(self, flux, delta_flux=1e-6) -> float:
        return self._spec.sensitivity(flux=flux, delta=delta_flux)
    
    def get_hamiltonian(self) -> Qobj:
        n_levels = self.n_levels
        n = self.n
        H_0 = (-self.EJ + 0.25 * self.EC) * qeye(n_levels)
        H_1 = self.frequency * (n + 0.5 * qeye(n_levels))
        H_2 = (self.anharmonicity / 2) * (n * n - n)
        return H_0 + H_1 + H_2
    
    def get_hamiltonian_rwa(self, omega_d) -> Qobj:
        n = self.n
        n_levels = self.n_levels
        Delta = self.frequency - omega_d
        return Delta * n + (self.anharmonicity / 2) * (n * n - n)
    
    def get_collapse_operators(self) -> list[Qobj]:
        gamma_1 = 1.0 / self.T1
        gamma_phi = 1.0 / self.T2 - 0.5 * gamma_1
        return [
            np.sqrt(gamma_1) * self.a,
            np.sqrt(gamma_phi) * self.n,
        ]
    
    def generate_1f_noise(self, t_lists, amplitude, f_min, f_max):
        # Verbatim port from src/qubit.py:128
        from sqc.simulation.noise import generate_1f_noise
        return generate_1f_noise(t_lists, amplitude, f_min, f_max)
    
    def calculate_state_projection(self, target_state):
        # verbatim from src/qubit.py:151
        ...
    
    def qubit_under_mag(self, Phi_signal, is_noise=False):
        # verbatim from src/qubit.py:164
        ...
    
    def qubit_under_mag_hamiltonian(self, qubit_t, t_list,
                                     frame=0, omega_d=None):
        # verbatim from src/qubit.py:185
        ...
    
    def qubit_in_mag(self, Phi_signal, frame=0, omega_d=None):
        """Legacy interface; wraps HamiltonianBuilder.build and stores
        results as instance attributes.
        
        Sets self.isinmag, self.mag_signal, self.freq_coeffs, self.H_list.
        """
        from sqc.simulation.hamiltonian import HamiltonianBuilder
        H_list, _ = HamiltonianBuilder.build(
            qubit=self._spec, flux_signal=Phi_signal,
            pulse=None, frame="lab" if frame == 0 else "rotating",
            omega_d=omega_d,
        )
        self.isinmag = True
        self.mag_signal = Phi_signal
        # Extract freq_coeffs from H_list[1][1] for backward compat
        self.freq_coeffs = np.asarray(H_list[1][1])
        self.H_list = H_list
    
    def change_flux(self, flux):
        """Legacy mutator. New code should use QubitSpec.with_flux()."""
        self.flux = flux
        self.EJ = self.EJ_0 * abs(math.cos(math.pi * flux))
        self._spec = self._spec.with_flux(flux)
        self.frequency = self._spec.frequency()
        self.anharmonicity = self._spec.anharmonicity()
        self.hamiltonian = self.get_hamiltonian()
    
    def optimal_work_point(self):
        return QubitSpec.optimal_work_point() * np.pi   # legacy returned in radians
    
    def ideal_gate(self, theta, phi):
        # verbatim from src/qubit.py:260
        ...
    
    def simulate_gate(self, theta, phi, T, sigma):
        # verbatim from src/qubit.py:269
        ...
```

**关键细节**:
- `qubit_in_mag` 不再自己写 freq_coeffs 计算逻辑,而是委托给 `HamiltonianBuilder.build`。
- `optimal_work_point()` 旧版返回 `np.arctan(np.sqrt(2))`(单位:弧度,**不是 Φ₀**),新版静态方法 `QubitSpec.optimal_work_point()` 返回单位 Φ₀,因此 legacy method 需要乘以 π 还原。**此处务必小心**,通过 baseline 测试验证。

#### 3.4.2 实现 HamiltonianBuilder

```python
# sqc/simulation/hamiltonian.py
from __future__ import annotations
from typing import Literal
import numpy as np
import math
from qutip import qeye, num, destroy, Qobj


class HamiltonianBuilder:
    """Builds time-dependent Hamiltonian list-format from device + flux + pulse.
    
    Pure functional: never mutates input objects.
    """
    
    @staticmethod
    def build(qubit, flux_signal, pulse,
              frame: Literal["lab", "rotating"] = "rotating",
              omega_d: float | None = None) -> tuple[list, np.ndarray]:
        """Returns (H_list, t_global) for QobjEvo consumption.
        
        H_list = [[op0, coeff_array_or_scalar], ...].
        
        If flux_signal is None and pulse is None, returns just [[H_static, 1]].
        If flux_signal is given, time-dependent freq_coeffs are computed.
        If pulse is given, pulse.hamiltonian is appended.
        """
        # Resolve QubitSpec: accept TransmonQubit or QubitSpec
        if hasattr(qubit, "spec"):
            spec = qubit.spec()
        else:
            spec = qubit
        
        n_levels = spec.n_levels
        n_op = num(n_levels)
        
        if flux_signal is None:
            t_list = np.linspace(0, 100, 100)
            freq_coeffs = np.full(len(t_list), spec.frequency())
        else:
            t_list = np.asarray(flux_signal.t_list)
            samples = np.asarray(flux_signal.samples
                                 if hasattr(flux_signal, "samples")
                                 else flux_signal.signal)
            freq_coeffs = np.array([
                spec.frequency(spec.flux_bias + s) for s in samples
            ])
        
        if frame == "rotating":
            if omega_d is None:
                omega_d = spec.frequency()
            freq_coeffs = freq_coeffs - omega_d
        
        H_list = [
            spec.anharmonicity() * 0.5 * (n_op * n_op - n_op),
        ]
        if frame == "lab":
            H_list.append([n_op + 0.5 * qeye(n_levels), freq_coeffs])
        else:
            H_list.append([n_op, freq_coeffs])
        
        if pulse is not None:
            for op, coeff in pulse.hamiltonian:
                H_list.append([op, coeff])
        
        return H_list, t_list
```

**关键约束**:
- 不修改 `spec`、`flux_signal`、`pulse` 中的任何字段。
- 返回的 H_list 与 src/qubit.py:qubit_in_mag 产生的 H_list **数值完全相同**(回归测试 transient_default、ramsey_default 必须通过)。

### 3.5 迁移 Cavity → Resonator,Coupled_System → CoupledSystem

逐字符 port 到 `sqc/devices/resonator.py` 和 `sqc/devices/chip.py`,只改类名和 docstring。所有方法保持原样(包括内部循环结构)。

### 3.6 迁移单/双比特门函数到 sqc/control/gates.py

把 `src/qubit.py` 中模块级的 `ideal_iSWAP`、`simulate_iSWAP`、`ideal_CZ`、`simulate_CZ` 移到 `sqc/control/gates.py`。在 `src/qubit.py` 镜像层重新导出。

注:`Coupled_System.simulate_iSWAP` 是**实例方法**,迁移到 `sqc/devices/chip.py: CoupledSystem.simulate_iSWAP`(保持为方法)。

### 3.7 迁移 Signal/CompositeSignal 行为到 sqc/control/

P1 阶段**不**做"拆分 Signal 为 Waveform + FluxSignal"的全面迁移(那会破坏太多 import)。改为以下方式:

#### 3.7.1 方案

`sqc/control/flux_signal.py` 提供 `FluxSignal` 类,**直接 inherit 自 Waveform**,但额外接受一个 type-based 工厂构造参数,以兼容旧 `Signal(type=N, t_list=..., **kwargs)` 调用方式:

```python
# sqc/control/flux_signal.py (P1 完整版)
from __future__ import annotations
from typing import Any
import numpy as np

from .waveform import Waveform, CompositeWaveform


class FluxSignal(Waveform):
    """Flux signal in units of Φ₀.
    
    Backward-compatible with src.signal.Signal: supports type-based
    construction.
    """
    
    def __init__(self, type=0, t_list=None, **kwargs):
        # Build samples first using legacy logic
        self._type = type
        self._params = self._fill_default_params(kwargs)
        if type == 6:
            self._build_basis_functions(t_list)
        samples = self._generate(type, t_list, self._params)
        # Initialize parent
        super().__init__(
            t_list=np.asarray(t_list),
            samples=samples,
            metadata={"type": type, "params": dict(kwargs)},
        )
    
    @property
    def signal(self) -> np.ndarray:
        """Backward-compat alias for samples."""
        return self.samples
    
    @signal.setter
    def signal(self, value):
        self.samples = np.asarray(value, dtype=float)
    
    @property
    def type(self) -> int:
        return self._type
    
    @property
    def params(self) -> dict:
        return self._params
    
    @property
    def basis_functions(self) -> list:
        return getattr(self, "_basis_functions", [])
    
    # Legacy methods
    def update_signal(self, **kwargs):
        # verbatim port from src/signal.py:192
        for k, v in kwargs.items():
            self._params[k] = v
        self.samples = self._generate(self._type, self.t_list, self._params)
    
    def truncate(self, t_start, t_end):
        # legacy: in-place. Override Waveform.truncate (which returns new).
        mask = (self.t_list < t_start) | (self.t_list > t_end)
        self.samples[mask] = 0.0
    
    def value_at(self, t):
        # verbatim port from src/signal.py:211
        ...
    
    def copy(self) -> "FluxSignal":
        # verbatim port
        ...
    
    # Internal helpers (verbatim port from src/signal.py)
    @staticmethod
    def _fill_default_params(params: dict) -> dict: ...
    @staticmethod
    def _generate(type, t_list, params) -> np.ndarray: ...
    def _build_basis_functions(self, t_list): ...


# Aliases
Signal = FluxSignal


class CompositeSignal(CompositeWaveform):
    """Backward-compat for src.signal.CompositeSignal."""
    
    def __init__(self, signals: list[FluxSignal]):
        super().__init__(
            t_list=np.array([]),
            samples=np.array([]),
            components=list(signals),
        )
        self.signals = list(signals)
        self.t_list = self._compute_t_list()
        self.samples = self._compute_samples()
    
    @property
    def signal(self) -> np.ndarray:
        return self.samples
    
    def _compute_t_list(self) -> np.ndarray: ...
    def _compute_samples(self) -> np.ndarray: ...
```

**保留语义**:旧 `Signal(type=2, ...)`、`Signal(type=8, signal=arr)`、`signal.signal`、`signal.t_list`、`signal.value_at(t)`、`signal.truncate(...)`、`signal.update_signal(...)`、`signal.params` 全部行为不变。

### 3.8 迁移 Pulse 到 sqc/control/pulse.py

逐字符 port `src/pulse.py:Pulse` 和 `CompositePulse`。两件事:

1. **去掉 `get_kernel` 方法**(D2 债务)。kernel 计算的逻辑挪到 `sqc/reconstruction/kernel.py: KernelEstimator`(P2 实现;P1 留 ABC 与 stub `raise NotImplementedError`)。
2. 在镜像层 `src/pulse.py` 中,如果有代码调用了旧 `pulse.get_kernel(qubit)`,通过 facade 兼容:`def get_kernel(self, qubit): from sqc.reconstruction.kernel import KernelEstimator; t_samples, kernel = KernelEstimator().estimate(self, qubit); self.t_samples = t_samples; self.kernel = kernel`

注:**P1 完成时 `KernelEstimator` 还是 stub**,所以 `get_kernel` 在 P1 后会 `raise NotImplementedError`。这暂时会破坏 case 4。

**解决**:P1 阶段的 `Pulse.get_kernel`(以及 `CompositePulse.get_kernel`)**逐字符保留** legacy 实现在 `sqc/control/pulse.py` 内,**作为 deprecated 方法**。P2 把它内化到 `KernelEstimator` 后,把 `Pulse.get_kernel` 改为转发调用。

### 3.9 迁移 sequence factory 到 sqc/control/sequence.py

`create_pulse`、`create_ramsey_pulse`、`create_diff_echo_pulse`、`create_echo_pulse`、`create_cpmg_pulse`、`create_cryoscope_pulse` 全部 port 到 `sqc/control/sequence.py`,逐字符。

### 3.10 实现 src/ 永久镜像层

按主方案 §8.2 模板实现。**所有 src/ 文件改写为镜像 + facade**。

#### 3.10.1 src/qubit.py

```python
"""src.qubit — sqc.devices.transmon 的兼容镜像。

新代码应直接使用 sqc.devices.* 。本文件为兼容历史代码而存在,
不应包含任何业务逻辑。

See idea/refactor/_refactor_plan.md §8.
"""
from sqc.devices.transmon import TransmonQubit, QubitSpec
from sqc.devices.resonator import Resonator as Cavity
from sqc.devices.chip import CoupledSystem as Coupled_System
from sqc.control.gates import (
    ideal_iSWAP,
    simulate_iSWAP_simple as simulate_iSWAP,
    ideal_CZ,
    simulate_CZ,
)

__all__ = [
    "TransmonQubit", "QubitSpec",
    "Cavity", "Coupled_System",
    "ideal_iSWAP", "simulate_iSWAP",
    "ideal_CZ", "simulate_CZ",
]
```

#### 3.10.2 src/signal.py

```python
"""src.signal — sqc.control.flux_signal 的兼容镜像。"""
from sqc.control.flux_signal import FluxSignal as Signal, CompositeSignal

__all__ = ["Signal", "CompositeSignal"]
```

#### 3.10.3 src/pulse.py

```python
"""src.pulse — sqc.control.pulse + sqc.control.sequence 的兼容镜像。"""
from sqc.control.pulse import Pulse, CompositePulse
from sqc.control.sequence import (
    create_pulse,
    create_ramsey_pulse,
    create_diff_echo_pulse,
    create_echo_pulse,
    create_cpmg_pulse,
    create_cryoscope_pulse,
)

__all__ = [
    "Pulse", "CompositePulse",
    "create_pulse", "create_ramsey_pulse", "create_diff_echo_pulse",
    "create_echo_pulse", "create_cpmg_pulse", "create_cryoscope_pulse",
]
```

#### 3.10.4 src/protocal.py 和 src/analysis.py

P1 阶段**不动**(它们本身是 Track A 后续 phase 的工作内容)。保留原文件。

### 3.11 简单的 noise.py 模块

```python
# sqc/simulation/noise.py
from __future__ import annotations
import numpy as np


def generate_1f_noise(t_list, amplitude, f_min, f_max,
                      seed: int | None = None) -> np.ndarray:
    """1/f noise generator. Verbatim from TransmonQubit.generate_1f_noise."""
    if seed is not None:
        np.random.seed(seed)
    dt = t_list[1] - t_list[0]
    n = len(t_list)
    freqs = np.fft.fftfreq(n, dt)
    spectrum = np.zeros(n, dtype=complex)
    for i in range(1, n // 2):
        f = abs(freqs[i])
        if f_min <= f <= f_max:
            spectrum[i] = (amplitude / np.sqrt(f) *
                           (np.random.normal() + 1j * np.random.normal()))
    spectrum[n // 2 + 1:] = np.conj(spectrum[1:n // 2][::-1])
    return np.fft.ifft(spectrum).real
```

### 3.12 ExperimentResult 与 result 工具

```python
# sqc/simulation/result.py
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import pickle

import numpy as np


@dataclass
class ExperimentResult:
    data: dict[str, np.ndarray] = field(default_factory=dict)
    axes: dict[str, np.ndarray] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)
    config: dict = field(default_factory=dict)
    
    def save(self, path: str | Path) -> None:
        with open(path, "wb") as f:
            pickle.dump(self, f, protocol=5)
    
    @classmethod
    def load(cls, path: str | Path) -> "ExperimentResult":
        with open(path, "rb") as f:
            obj = pickle.load(f)
        if not isinstance(obj, cls):
            raise TypeError(f"Loaded object is {type(obj).__name__}, "
                            f"not ExperimentResult")
        return obj


@dataclass(frozen=True)
class MeasurementTrace:
    axis: np.ndarray
    p_e: np.ndarray
    p_e_iq: tuple[np.ndarray, np.ndarray] | None = None
    metadata: dict = field(default_factory=dict)


def extract_expectation(result, e_ops_index: int = 0) -> np.ndarray:
    """Extract expectation values from a QuTiP result object."""
    if hasattr(result, "expect"):
        if isinstance(result.expect, list) and len(result.expect) > e_ops_index:
            return np.asarray(result.expect[e_ops_index])
        return np.asarray(result.expect)
    raise ValueError("Result has no .expect attribute")


def extract_population(result, level: int) -> np.ndarray:
    """Extract population of a Fock level from result.states."""
    from qutip import basis, expect
    if not hasattr(result, "states"):
        raise ValueError("Result has no .states attribute (set store_states=True)")
    populations = []
    for state in result.states:
        n_levels = state.dims[0][0]
        proj = basis(n_levels, level) * basis(n_levels, level).dag()
        populations.append(expect(proj, state))
    return np.asarray(populations)
```

---

## 4. 验收标准

### 4.1 结构完整性

```bash
# 4.1.1 所有目录与 __init__.py 存在
find sqc -type d -exec ls {} \;

# 4.1.2 所有 ABC 定义无语法错误
python -c "
from sqc.devices.base import Device
from sqc.hardware.distortion import DistortionModel
from sqc.hardware.readout import ReadoutModel
from sqc.control.pulse import PulseBase
from sqc.simulation.runner import RunnerBase, MesolveRunner, SlidingMeasurementRunner
from sqc.experiments.base import Experiment
from sqc.calibration.base import Calibration
from sqc.reconstruction.base import Reconstruction
from sqc.workflows.base import Workflow
print('all ABCs ok')
"

# 4.1.3 所有数据结构可实例化
python -c "
from sqc.devices.transmon import QubitSpec
from sqc.control.waveform import Waveform, CompositeWaveform
from sqc.control.flux_signal import FluxSignal
from sqc.simulation.result import ExperimentResult, MeasurementTrace
import numpy as np
spec = QubitSpec(name='Q0', EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000)
print('QubitSpec OK, frequency =', spec.frequency())
"
```

### 4.2 镜像等价性(关键)

```bash
# 4.2.1 旧 import 仍然工作
python -c "
from src.qubit import TransmonQubit, Cavity, Coupled_System
from src.signal import Signal, CompositeSignal
from src.pulse import (Pulse, CompositePulse, create_pulse,
                       create_ramsey_pulse, create_diff_echo_pulse,
                       create_echo_pulse, create_cpmg_pulse,
                       create_cryoscope_pulse)
print('src/* imports OK')
"

# 4.2.2 镜像与 sqc 是同一对象
python -c "
import src.qubit as a
import sqc.devices.transmon as b
assert a.TransmonQubit is b.TransmonQubit, 'TransmonQubit not aliased'
print('alias check OK')
"
```

### 4.3 物理回归

```bash
# 4.3.1 所有 Phase 0 的 baseline 测试仍然通过
pytest tests/regression -m regression -v
```

期望:全部 pass,**不需要重新生成 baseline**。

### 4.4 端到端

```bash
# 4.4.1 web_demo 启动并跑 Ramsey
python web_demo.py
# 浏览器打开,选 Ramsey,运行,确认有结果

# 4.4.2 Notebook 至少前 5 个 cell 可运行
# (手动验证)
```

### 4.5 单元测试新增

P1 引入的新模块必须有 unit test。最低集合(`tests/unit/`):

| 文件 | 测试 |
|---|---|
| `test_qubit_spec.py` | QubitSpec 不可变性、frequency/anharmonicity/sensitivity 数值 |
| `test_waveform.py` | Waveform 构造、value_at、truncate 返回新对象、CompositeWaveform 拼接 |
| `test_flux_signal.py` | type=0..8 的 FluxSignal 与旧 Signal 输出相同 samples |
| `test_hamiltonian_builder.py` | HamiltonianBuilder.build 输出与旧 qubit_in_mag 数值相等 |
| `test_imports.py` | 所有 ABC、数据结构、镜像 import 不报错 |

每个 test 文件至少 3 个测试用例。

---

## 5. 测试要求

### 5.1 必须新增的测试

| 测试文件 | 关键 case |
|---|---|
| `tests/unit/test_qubit_spec.py` | QubitSpec frozen 检查、frequency 数值、optimal_work_point 返回 Φ₀ |
| `tests/unit/test_waveform.py` | t_list/samples shape 检查、truncate 不修改原对象、value_at 越界返回 0 |
| `tests/unit/test_flux_signal.py` | type=2 (sinusoidal) 与 type=8 (custom) 与 src/signal.py 输出对齐 |
| `tests/unit/test_hamiltonian_builder.py` | flux=None 情况;rotating + omega_d=qubit.frequency 情况;flux_signal 给定时与旧 qubit_in_mag 输出一致 |
| `tests/unit/test_pulse_mirror.py` | `from src.pulse import Pulse` 与 `from sqc.control.pulse import Pulse` 是同一对象 |

### 5.2 必须保持的测试

`tests/regression/test_physics_baseline.py` 中**所有用例继续 pass**。

---

## 6. 风险与回滚

| # | 风险 | 缓解 |
|---|---|---|
| 6.1 | TransmonQubit 迁移漏掉某个属性,导致 src/protocal.py 在某个 case 中 AttributeError | tests/unit/test_imports.py + tests/regression 端到端覆盖;若漏属性,baseline 测试会失败 |
| 6.2 | `optimal_work_point()` 单位混淆(rad vs Φ₀) | test_qubit_spec.py 显式断言两个版本各自的单位 |
| 6.3 | HamiltonianBuilder 在 frame=lab 时与旧实现差异 | test_hamiltonian_builder.py 同时覆盖 frame=lab 和 frame=rotating |
| 6.4 | Signal type=6 的 basis_functions 字段在镜像里丢失 | 显式在 FluxSignal 中保留 self._basis_functions 属性 + 测试 |
| 6.5 | sqc/__init__.py 引入循环依赖 | 子包 __init__.py 保持空;只在叶子模块 import |
| 6.6 | pickle 反序列化时找不到 sqc.* 类(旧 baseline pickle 含 src.* 类引用) | baseline 用 numpy 数组而非 Signal 对象,无引用问题 |

**回滚策略**:本 phase 改动较大但 src/ 保留 facade,可分两级回滚:

- **轻度回滚**:某个 sqc 模块出问题,把 src/*.py 恢复为原版(git checkout),sqc/ 保留作为 dead code。
- **完全回滚**:`git rm -r sqc/` + `git checkout src/` 即可回到 P0 状态。

---

## 7. 输出物

### 7.1 文件清单(新增)

```
sqc/
├── __init__.py
├── devices/
│   ├── __init__.py
│   ├── base.py
│   ├── transmon.py
│   ├── resonator.py
│   ├── coupler.py                  (空 ABC,P5 实现)
│   └── chip.py
├── hardware/
│   ├── __init__.py
│   ├── control_line.py             (dataclass + ABC 占位)
│   ├── distortion.py               (ABC 占位)
│   ├── transfer_matrix.py          (dataclass + ABC 占位)
│   ├── electronics.py              (空,P5)
│   └── readout.py                  (ABC + IQ_readout_legacy facade)
├── control/
│   ├── __init__.py
│   ├── waveform.py
│   ├── flux_signal.py
│   ├── pulse.py
│   ├── sequence.py
│   ├── schedules.py                (空,P5)
│   └── gates.py
├── simulation/
│   ├── __init__.py
│   ├── hamiltonian.py
│   ├── runner.py                   (ABC + stub,P2 实现)
│   ├── noise.py
│   └── result.py
├── experiments/
│   ├── __init__.py
│   └── base.py                     (Experiment ABC)
├── calibration/
│   ├── __init__.py
│   └── base.py                     (Calibration ABC + CalibrationTable)
├── reconstruction/
│   ├── __init__.py
│   ├── base.py                     (Reconstruction ABC)
│   └── kernel.py                   (KernelEstimator,P1 verbatim port,P2 整合)
└── workflows/
    ├── __init__.py
    └── base.py                     (Workflow ABC)
```

### 7.2 文件清单(修改)

```
src/qubit.py                        (改为镜像)
src/signal.py                       (改为镜像)
src/pulse.py                        (改为镜像)
src/__init__.py                     (无变化,保留)
```

### 7.3 文件清单(不改)

```
src/protocal.py                     (P2 改)
src/analysis.py                     (P2/P3 改)
web_demo.py                         (永久不改)
Simulation.ipynb                    (永久不改)
```

### 7.4 接口快照

| 符号 | 来源 | 调用范例 |
|---|---|---|
| `QubitSpec(name, EC, EJ, T1, T2, flux_bias=0, n_levels=3)` | sqc.devices.transmon | `QubitSpec("Q0", 2*pi*0.2, 2*pi*15, 1e4, 8e3)` |
| `TransmonQubit(EC, EJ, T1, T2, flux=0, state=0, n_levels=3, name="Q")` | sqc.devices.transmon | 与 src/qubit.py 完全等价 |
| `Waveform(t_list, samples, metadata={})` | sqc.control.waveform | |
| `FluxSignal(type=N, t_list=..., **params)` | sqc.control.flux_signal | 与 src/signal.py:Signal 等价 |
| `HamiltonianBuilder.build(qubit, flux_signal, pulse, frame, omega_d) -> (H_list, t_list)` | sqc.simulation.hamiltonian | 静态方法,无副作用 |
| `ExperimentResult(data, axes, metadata, config)` | sqc.simulation.result | `.save(path)` / `.load(path)` |
| `extract_expectation(result, idx)` | sqc.simulation.result | |
| `extract_population(result, level)` | sqc.simulation.result | |

---

## 8. 完成确认清单

```
- [ ] sqc/ 目录骨架完整(9 个子包 + 至少 25 个 .py)
- [ ] 所有 ABC 文件无语法错误,可被 import
- [ ] 所有数据结构可实例化(QubitSpec/Waveform/FluxSignal/ExperimentResult/...)
- [ ] HamiltonianBuilder.build 实现并通过测试
- [ ] sqc/devices/transmon.py 完整实现 TransmonQubit + QubitSpec
- [ ] sqc/control/{waveform,flux_signal,pulse,sequence,gates}.py 全部完成迁移
- [ ] src/qubit.py / signal.py / pulse.py 改为镜像
- [ ] src/protocal.py / analysis.py 保持原状(不改)
- [ ] tests/unit/ 至少 5 个新文件,每个至少 3 个用例
- [ ] pytest tests/unit -v 全部 pass
- [ ] pytest tests/regression -m regression -v 全部 pass(关键!)
- [ ] python web_demo.py Ramsey/Rabi/Transient 三个功能均跑通
- [ ] Simulation.ipynb 前 10 个 cell 可运行(手动验证)
```

**触发 P2** 的条件:上述清单全部打勾,且 commit message 引用本清单。

---

**下一步**:[phase_2_handbook.md](phase_2_handbook.md) — 已实现协议的实验对象化与反向 wrapper。
