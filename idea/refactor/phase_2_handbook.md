# Phase 2 Handbook — 已实现协议的实验对象化 + KernelEstimator 去重

> 前置阅读:[_refactor_plan.md](_refactor_plan.md) §6.5、§7.4、[phase_1_handbook.md](phase_1_handbook.md)  
> 估计工时:4–6 天  
> 触发条件:Phase 1 完成 + Track B 已修复 LM 收敛(_TODO_master.md 0.1–0.3 完成,**仅 0.3 是硬依赖**)  
> 完成标志:`RamseyExperiment` / `DiffEchoExperiment` / `TransientSensingExperiment` / `RabiExperiment` 实现并通过 baseline 回归;`KernelEstimator` 在三处 kernel 实现去重

---

## 1. 目标

把当前 `src/protocal.py` 中 case 0/1/2/4 的实现**对象化**到 `sqc/experiments/`,消除 D3、D4 债务。同时:

1. 实现 `KernelEstimator`,把 `src/pulse.py` 和 `src/analysis.py` 中三份重复的 kernel 计算合并为一个(消除 D2 债务)。
2. 实现 `MesolveRunner` 和 `SlidingMeasurementRunner`,把 `single_measurement` 和 `sliding_measrement` 从 `Protocal` 中抽离(并修正拼写)。
3. 实现 `IQReadoutModel`,把 `IQ_readout()` 类化(消除 D4 债务)。
4. 在 `src/protocal.py` 中,case 0/1/2/4 改为转发调用新 Experiment 类(facade 模式)。
5. 物理回归测试 100% 通过。

**本 phase 完成后**,case 0/1/2/4 在 sqc/ 中有"原生实现",src/protocal.py 是 facade。case 3/5/6/7/8 仍由 src/protocal.py 处理(P3 内化)。

---

## 2. 前置条件

| # | 条件 | 验证方式 |
|---|---|---|
| 2.1 | Phase 1 完成 | sqc/ 完整骨架,镜像测试通过 |
| 2.2 | (强烈推荐)Track B `_TODO_master.md 0.1` 完成 | case 1 死代码清理;不影响 baseline |
| 2.3 | Track B `_TODO_master.md 0.2` 完成 | kernel 自动校准;**完成后必须重生成 transient_default.pkl** |
| 2.4 | (Optional but recommended)Track B `_TODO_master.md 0.3` 完成 | LM 收敛修复 |
| 2.5 | git 工作区干净,P1 commit 已合并 | `git status` |

---

## 3. 任务清单

### 3.1 实现 KernelEstimator (D2 债务)

#### 3.1.1 接口

```python
# sqc/reconstruction/kernel.py
from __future__ import annotations
from dataclasses import dataclass

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.devices.transmon import TransmonQubit, QubitSpec
from sqc.control.pulse import CompositePulse, PulseBase
from sqc.control.flux_signal import FluxSignal


@dataclass
class KernelEstimator:
    """Estimate the control kernel of a pulse (or composite pulse)
    by perturbing with a narrow Gaussian stimulus at each time point.
    
    Replaces the three duplicate implementations in:
      - src/pulse.py:Pulse.get_kernel
      - src/pulse.py:CompositePulse.get_kernel
      - src/analysis.py:Analysis.get_kernel
    """
    stim_amplitude: float = 0.0215      # default matches legacy
    stim_width: float = 3.0             # ns
    auto_calibrate: bool = False        # if True, adjust amplitude per qubit
    
    def estimate(self, pulse: PulseBase | CompositePulse,
                 qubit: TransmonQubit) -> tuple[np.ndarray, np.ndarray]:
        """Return (t_samples, kernel) for the given pulse + qubit.
        
        Numerically equivalent to legacy CompositePulse.get_kernel.
        """
        t_list = np.asarray(pulse.t_list)
        n_levels = qubit.n_levels
        psi_e = basis(n_levels, 1)
        
        # Static H_0
        if pulse.frame == 0:
            H_0 = QobjEvo(qubit.get_hamiltonian())
        else:
            H_0 = QobjEvo(qubit.get_hamiltonian_rwa(qubit.frequency))
        
        H_pulse = QobjEvo(pulse.hamiltonian, tlist=t_list, order=1)
        H_base = H_0 + H_pulse
        result_base = mesolve(H_base, qubit.state, t_list, [],
                              e_ops=[psi_e * psi_e.dag()])
        p_e_base = result_base.expect[0][-1]
        
        kernel = np.zeros(len(t_list))
        for i, t_i in enumerate(t_list):
            stim = FluxSignal(
                type=3, t_list=t_list,
                amplitude=self._amplitude_for(qubit),
                center=t_i, width=self.stim_width,
            )
            stim_area = np.trapezoid(stim.signal, stim.t_list)
            qubit_t = qubit.qubit_under_mag(stim)
            H_stim = QobjEvo(
                qubit.qubit_under_mag_hamiltonian(
                    qubit_t, stim.t_list, pulse.frame, pulse.omega_d
                ),
                tlist=stim.t_list, order=1,
            )
            H_total = H_0 + H_pulse + H_stim
            res = mesolve(H_total, qubit.state, t_list, [],
                          e_ops=[psi_e * psi_e.dag()])
            kernel[i] = (res.expect[0][-1] - p_e_base) / stim_area
        
        return t_list.copy(), kernel
    
    def _amplitude_for(self, qubit: TransmonQubit) -> float:
        if not self.auto_calibrate:
            return self.stim_amplitude
        # Auto-calibration: pick amplitude such that the stimulus produces
        # a frequency shift ~ 1% of qubit anharmonicity. (Algorithm pending,
        # see _TODO_master.md 0.2.)
        kappa = qubit.frequency_sensitivity(qubit.flux)
        if kappa == 0:
            return self.stim_amplitude
        target_shift = 0.01 * abs(qubit.anharmonicity)
        return target_shift / abs(kappa)
```

#### 3.1.2 替换 src/pulse.py 中的 get_kernel

`sqc/control/pulse.py: Pulse.get_kernel` 和 `CompositePulse.get_kernel` 改为转发:

```python
def get_kernel(self, qubit):
    """Backward-compat. Internally delegates to KernelEstimator.
    
    Sets self.t_samples and self.kernel as side effects (legacy).
    """
    from sqc.reconstruction.kernel import KernelEstimator
    estimator = KernelEstimator()
    self.t_samples, self.kernel = estimator.estimate(self, qubit)
    return self.t_samples, self.kernel
```

src/pulse.py:Analysis.get_kernel 同样转发(注:`Analysis` 是 P3 才类化的对象,P2 阶段先在 src/analysis.py 中改 `Analysis.get_kernel` 转发即可)。

### 3.2 实现 MesolveRunner

```python
# sqc/simulation/runner.py
from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from qutip import QobjEvo, mesolve, basis

from sqc.simulation.result import ExperimentResult


class RunnerBase(ABC):
    @abstractmethod
    def run(self, *args, **kwargs) -> ExperimentResult: ...


@dataclass
class MesolveRunner(RunnerBase):
    options: dict | None = None
    
    def run(self, H_list, psi0, t_list, c_ops, e_ops,
            store_states: bool = False) -> ExperimentResult:
        opts = dict(self.options or {})
        opts.setdefault("store_states", store_states)
        H = QobjEvo(H_list, tlist=t_list, order=1)
        result = mesolve(H, psi0, t_list, c_ops, e_ops, options=opts)
        return ExperimentResult(
            data={"expect": np.array(result.expect)},
            axes={"t": np.asarray(t_list)},
            metadata={"runner": "MesolveRunner"},
        )
```

### 3.3 实现 SlidingMeasurementRunner

直接 port 自 `src/protocal.py:Protocal.sliding_measrement`(244–321 行)的逻辑,但:
- 拆为独立类,不依赖 `Protocal`。
- 输入:qubit、flux_signal、composite_pulse,可选 scan_list。
- 输出:`ExperimentResult`,含 `axes["scan"]` 和 `data["p_e"]`。

```python
# sqc/simulation/runner.py (续)
from sqc.devices.transmon import TransmonQubit
from sqc.control.flux_signal import FluxSignal
from sqc.control.pulse import CompositePulse


@dataclass
class SlidingMeasurementRunner(RunnerBase):
    """Slides a control pulse across a flux signal, measures p_e at each delay.
    
    Replaces Protocal.sliding_measrement (note legacy spelling).
    """
    options: dict | None = None
    
    def run(self, qubit: TransmonQubit, flux_signal: FluxSignal,
            control_pulse: CompositePulse,
            scan_list: np.ndarray | None = None) -> ExperimentResult:
        # ... port from src/protocal.py:252-321 verbatim,
        # but renamed to "scan_list" instead of legacy "scan_list".
        ...
        return ExperimentResult(
            data={"p_e": p_e},
            axes={"scan": scan_list},
            metadata={"runner": "SlidingMeasurementRunner"},
        )
```

### 3.4 实现 IQReadoutModel

```python
# sqc/hardware/readout.py
from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.devices.transmon import TransmonQubit
from sqc.control.sequence import create_ramsey_pulse
from sqc.simulation.result import ExperimentResult


class ReadoutModel(ABC):
    @abstractmethod
    def measure(self, *args, **kwargs) -> dict[str, float]: ...


@dataclass
class IdealProjectiveReadout(ReadoutModel):
    """Projective measurement onto |1⟩."""
    
    def measure(self, state, qubit: TransmonQubit) -> dict[str, float]:
        psi_e = basis(qubit.n_levels, 1)
        proj = psi_e * psi_e.dag()
        from qutip import expect
        return {"p_e": float(expect(proj, state).real)}


@dataclass
class IQReadoutModel(ReadoutModel):
    """IQ readout via two Ramsey sequences with π/2 phase offset.
    
    Replaces src/protocal.py:IQ_readout.
    """
    tau: float = 20.0          # ns, free precession
    t_rabi: np.ndarray | None = None
    
    def measure(self, qubit: TransmonQubit, **extra) -> dict[str, float]:
        """Run two Ramsey sequences (φ₂=0 and φ₂=π/2) and return I, Q.
        
        Requires qubit.H_list and qubit.mag_signal to be populated
        (call qubit.qubit_in_mag(...) first).
        """
        if self.t_rabi is None:
            self.t_rabi = np.linspace(0, 10, 20)
        omega_d = qubit.frequency
        ctrl_I = create_ramsey_pulse(self.t_rabi, self.tau,
                                      omega_d=omega_d,
                                      phase1=np.pi/2, phase2=0.0)
        ctrl_Q = create_ramsey_pulse(self.t_rabi, self.tau,
                                      omega_d=omega_d,
                                      phase1=np.pi/2, phase2=np.pi/2)
        H_I = (QobjEvo(ctrl_I.hamiltonian, tlist=ctrl_I.t_list, order=1) +
               QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1))
        H_Q = (QobjEvo(ctrl_Q.hamiltonian, tlist=ctrl_Q.t_list, order=1) +
               QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1))
        psi_e = basis(qubit.n_levels, 1)
        result_I = mesolve(H_I, qubit.state, ctrl_I.t_list, [],
                           e_ops=[psi_e * psi_e.dag()])
        result_Q = mesolve(H_Q, qubit.state, ctrl_Q.t_list, [],
                           e_ops=[psi_e * psi_e.dag()])
        return {
            "p_e_I": float(result_I.expect[0][-1]),
            "p_e_Q": float(result_Q.expect[0][-1]),
        }


# Backward-compat function for src/protocal.py
def IQ_readout_legacy(qubit, type, **kwargs):
    """Drop-in replacement for src/protocal.py:IQ_readout."""
    if type in (2, 3):  # cryoscope calib / measurement
        readout = IQReadoutModel(
            tau=kwargs.get("tau", 20),
            t_rabi=kwargs.get("t_rabi", np.linspace(0, 10, 20)),
        )
        result = readout.measure(qubit)
        return result["p_e_I"], result["p_e_Q"]
    else:
        raise NotImplementedError(f"IQ_readout type={type} not supported")
```

### 3.5 实现 RabiExperiment

```python
# sqc/experiments/rabi.py
from __future__ import annotations
from dataclasses import dataclass

import numpy as np
from qutip import QobjEvo, basis

from sqc.devices.transmon import TransmonQubit
from sqc.control.sequence import create_pulse
from sqc.simulation.runner import MesolveRunner
from sqc.simulation.result import ExperimentResult
from sqc.experiments.base import Experiment


@dataclass
class RabiExperiment(Experiment):
    qubit: TransmonQubit
    t_rabi: np.ndarray | None = None
    omega_d: float | None = None
    
    def __post_init__(self):
        if self.t_rabi is None:
            self.t_rabi = np.linspace(0, 40, 1000)
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency
    
    def build_sequence(self):
        return create_pulse(
            self.qubit, frame=1, type=1, t_list=self.t_rabi,
            omega_d=self.omega_d, phase=0.0,
        )
    
    def run(self) -> ExperimentResult:
        psi_e = basis(self.qubit.n_levels, 1)
        H_0 = self.qubit.get_hamiltonian_rwa(self.qubit.frequency)
        H_pulse = self.build_sequence()
        H_rabi = QobjEvo(H_0) + QobjEvo(H_pulse, tlist=self.t_rabi)
        runner = MesolveRunner()
        result = runner.run(
            H_list=[H_0] + [[op, c] for op, c in H_pulse.ops_lst()],
            psi0=self.qubit.state, t_list=self.t_rabi,
            c_ops=[], e_ops=[psi_e * psi_e.dag()],
        )
        result.metadata["experiment"] = "RabiExperiment"
        result.metadata["qubit_spec"] = self.qubit.spec()
        return result
```

(注:上面 H_pulse 的处理需要根据 create_pulse 实际返回类型调整;此处仅示意结构。)

### 3.6 实现 RamseyExperiment

最重要的一个。**关键**:必须复现 `src/protocal.py:46-99` 的行为(返回 Phi、tau_list、p_e_list)。

```python
# sqc/experiments/ramsey.py
from __future__ import annotations
from dataclasses import dataclass, field

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.devices.transmon import TransmonQubit
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.simulation.result import ExperimentResult
from sqc.experiments.base import Experiment


@dataclass
class RamseyExperiment(Experiment):
    """Ramsey protocol with optional flux signal.
    
    Replaces Protocal.evolve case 1.
    Returns ExperimentResult with:
      - data["p_e"]: shape (len(tau_list),)
      - data["flux_samples"]: snapshot of flux signal
      - axes["tau"]: tau_list
      - axes["t_flux"]: flux signal time axis
    """
    qubit: TransmonQubit
    flux_signal: FluxSignal | None = None
    omega_d: float | None = None
    t_rabi: np.ndarray = field(default_factory=lambda: np.linspace(0, 20, 40))
    tau_list: np.ndarray = field(default_factory=lambda: np.linspace(0, 250, 500))
    t_global: np.ndarray = field(default_factory=lambda: np.linspace(-50, 300, 700))
    phase1: float = 0.0
    phase2: float = 0.0
    
    def __post_init__(self):
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency
        if self.flux_signal is None:
            # Default test signal (matches Protocal case 1)
            self.flux_signal = FluxSignal(
                type=2, t_list=np.linspace(0, 250, 500),
                amplitude=0.001, frequency=0.01,
                rise=10, fall=10, center=100, noise_level=0.0,
            )
    
    def build_sequence(self):
        # Sequence is built per-tau in run(), not statically
        return None
    
    def run(self) -> ExperimentResult:
        # 1. Couple flux to qubit (legacy side effect; HamiltonianBuilder
        # produces the same H_list)
        self.qubit.qubit_in_mag(self.flux_signal, frame=1, omega_d=self.omega_d)
        
        psi_e = basis(self.qubit.n_levels, 1)
        p_e_list = np.zeros(len(self.tau_list))
        for i, tau in enumerate(self.tau_list):
            ctrl = create_ramsey_pulse(
                self.t_rabi, tau, omega_d=self.omega_d,
                phase1=self.phase1, phase2=self.phase2,
            )
            ctrl.t_list = ctrl.t_list - self.t_rabi[-1]
            H = (QobjEvo(self.qubit.H_list,
                          tlist=self.qubit.mag_signal.t_list, order=1) +
                 QobjEvo(ctrl.hamiltonian, tlist=ctrl.t_list, order=1))
            result = mesolve(H, self.qubit.state, self.t_global, [],
                             e_ops=[psi_e * psi_e.dag()])
            p_e_list[i] = result.expect[0][-1]
        
        return ExperimentResult(
            data={
                "p_e": p_e_list,
                "flux_samples": self.flux_signal.signal.copy(),
            },
            axes={
                "tau": self.tau_list.copy(),
                "t_flux": self.flux_signal.t_list.copy(),
            },
            metadata={
                "experiment": "RamseyExperiment",
                "qubit_spec": self.qubit.spec(),
                "omega_d": self.omega_d,
            },
            config={
                "phase1": self.phase1, "phase2": self.phase2,
                "t_rabi": self.t_rabi.copy(),
                "t_global": self.t_global.copy(),
            },
        )
```

### 3.7 实现 DiffEchoExperiment 和 TransientSensingExperiment

完全类比 RamseyExperiment 的写法,逐字 port `src/protocal.py:117-141` (case 2) 和 `src/protocal.py:144-163` (case 4)。

**关键**:
- DiffEchoExperiment 的默认参数与 case 2 完全一致(k=5, t_list=linspace(0,100,200), Phi type=3 amplitude=0.01...)。
- TransientSensingExperiment 的默认参数与 case 4 完全一致,且**调用 KernelEstimator** 计算 kernel。

### 3.8 改写 src/protocal.py 的 Protocal 类为 facade

```python
# src/protocal.py (P2 之后)
"""src.protocal — facade over sqc.experiments and sqc.calibration.

Legacy API preserved verbatim:
  Protocal(type=N).evolve(qubit) -> tuple-of-things-matching-old-behavior

Implementation delegates to sqc.experiments.{rabi,ramsey,...}.
"""
from __future__ import annotations

import numpy as np
from qutip import basis

from sqc.devices.transmon import TransmonQubit
from sqc.control.flux_signal import FluxSignal as Signal, CompositeSignal
from sqc.control.pulse import Pulse, CompositePulse
from sqc.control.sequence import (
    create_pulse, create_ramsey_pulse, create_diff_echo_pulse,
    create_cpmg_pulse, create_cryoscope_pulse,
)
from sqc.hardware.readout import IQ_readout_legacy as IQ_readout

# P2 internalized
from sqc.experiments.rabi import RabiExperiment
from sqc.experiments.ramsey import RamseyExperiment
from sqc.experiments.echo import DiffEchoExperiment
from sqc.experiments.transient import TransientSensingExperiment


class Protocal:
    """Legacy protocol class, now a thin facade."""
    
    def __init__(self, type=0, **kwargs):
        self.type = type
        self.params = kwargs
    
    def initialize(self, qubit, state=0):
        # Verbatim from old src/protocal.py:23-44
        ...
    
    def evolve(self, qubit):
        match self.type:
            case 0:
                exp = RabiExperiment(qubit=qubit)
                result = exp.run()
                # Legacy returned just `result`; mimic that
                return result   # may need raw mesolve result; see below
            case 1:
                exp = RamseyExperiment(qubit=qubit)
                result = exp.run()
                # Legacy returned (Phi, tau_list, p_e_list)
                return (
                    exp.flux_signal,
                    result.axes["tau"],
                    result.data["p_e"].tolist(),
                )
            case 2:
                exp = DiffEchoExperiment(qubit=qubit)
                result = exp.run()
                # Legacy returned (Phi, tau_list, p_e_list, k, t_int)
                return (
                    exp.flux_signal,
                    result.axes["tau"],
                    result.data["p_e"].tolist(),
                    exp.k,
                    exp.t_int,
                )
            case 3:
                pass   # CPMG stub, P3+
            case 4:
                exp = TransientSensingExperiment(qubit=qubit)
                result = exp.run()
                # Legacy returned (t_samples, kernel, scan_list, delta_p, p_e, Phi, control_pulse)
                return (
                    result.axes["t_samples"], result.data["kernel"],
                    result.axes["scan"], result.data["delta_p"],
                    result.data["p_e"], exp.flux_signal,
                    exp.control_pulse,
                )
            case 5:
                # Cryoscope: P3 internalized, leave delegate to CryoscopeExperiment
                from sqc.experiments.cryoscope import CryoscopeExperiment
                exp = CryoscopeExperiment(qubit=qubit)
                result = exp.run()
                return (
                    exp.trunc_list, result.data["varphi"],
                    exp.flux_signal, [result.data["p_e_I"], result.data["p_e_Q"]],
                )
            case _:
                raise ValueError(f"Unknown protocol type {self.type}")
    
    # Helper methods kept as shims (delegate to runners)
    def single_measurement(self, qubit, Phi, ctrl, t_delay,
                           qubit_t=None, H=None, t_evole=None, index=None):
        """Backward-compat. Delegates to MesolveRunner."""
        ...
    
    def sliding_measrement(self, qubit, Phi, ctrl):
        """Backward-compat (legacy spelling preserved)."""
        from sqc.simulation.runner import SlidingMeasurementRunner
        runner = SlidingMeasurementRunner()
        result = runner.run(qubit, Phi, ctrl)
        return result.axes["scan"], result.data["p_e"].tolist()


# Calibration is P3-internalized; for P2, leave the original src/protocal.py
# Calibration class and IQ_readout function intact below this line.
class Calibration:
    # ... (verbatim original until P3) ...
    pass
```

**必读注释**:`Protocal.evolve` 必须返回与旧版**完全相同**的元组结构(顺序、类型),否则 web_demo.py 和 Notebook 会崩。Phase 2 的 baseline 测试就是验证这一点。

### 3.9 更新 src/analysis.py:Analysis.get_kernel 转发

```python
# src/analysis.py
class Analysis:
    def get_kernel(self, control_pulse, qubit):
        """Now delegates to KernelEstimator."""
        from sqc.reconstruction.kernel import KernelEstimator
        return KernelEstimator().estimate(control_pulse, qubit)
    
    # ... rest verbatim until P3 ...
```

---

## 4. 验收标准

### 4.1 Experiment 类工作正常

```bash
python -c "
import numpy as np
from sqc.devices.transmon import TransmonQubit
from sqc.experiments.ramsey import RamseyExperiment

q = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000)
exp = RamseyExperiment(qubit=q)
result = exp.run()
assert 'p_e' in result.data
assert len(result.data['p_e']) == len(result.axes['tau'])
print('Ramsey OK, p_e shape:', result.data['p_e'].shape)
"
```

### 4.2 KernelEstimator 等价性

```bash
python -c "
import numpy as np
from sqc.devices.transmon import TransmonQubit
from sqc.control.sequence import create_ramsey_pulse
from sqc.reconstruction.kernel import KernelEstimator

q = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000)
ctrl = create_ramsey_pulse(np.linspace(0,10,20), tau=0.0, omega_d=q.frequency)
ke = KernelEstimator()
t, k = ke.estimate(ctrl, q)
print('kernel shape:', k.shape, 'sum:', k.sum())
"
```

### 4.3 物理回归(关键)

```bash
pytest tests/regression -m regression -v
```

**全部 pass**。如果 Track B 已重生成 transient_default.pkl,则用新版本回归。

### 4.4 web_demo 与 Notebook

- `python web_demo.py` 启动并跑 Rabi/Ramsey/Transient 三个协议。
- `Simulation.ipynb` 至少前 15 个 cell 可运行。

### 4.5 单元测试新增

| 文件 | 用例 |
|---|---|
| `tests/unit/test_kernel_estimator.py` | KernelEstimator 与 legacy CompositePulse.get_kernel 输出一致 |
| `tests/unit/test_runners.py` | MesolveRunner 输出 ExperimentResult 形状正确;SlidingMeasurementRunner 与旧 sliding_measrement 输出 p_e 等价 |
| `tests/unit/test_iq_readout.py` | IQReadoutModel.measure 与 IQ_readout(qubit, type=3) 输出一致 |
| `tests/integration/test_ramsey_experiment.py` | RamseyExperiment.run 与 Protocal(type=1).evolve 输出 p_e 等价 |
| `tests/integration/test_diff_echo_experiment.py` | 同上 |
| `tests/integration/test_transient_experiment.py` | 同上 |

---

## 5. 测试要求

### 5.1 必须新增的等价性测试

每个 Experiment 类都要有"与 baseline 比对"的测试。例:

```python
# tests/integration/test_ramsey_experiment.py
import numpy as np
import pytest
from sqc.devices.transmon import TransmonQubit
from sqc.experiments.ramsey import RamseyExperiment
from tests.conftest import load_baseline, assert_array_close


pytestmark = pytest.mark.integration


def test_ramsey_experiment_matches_baseline(qubit_default):
    exp = RamseyExperiment(qubit=qubit_default)
    result = exp.run()
    bl = load_baseline("ramsey_default")
    assert_array_close(result.data["p_e"], bl["p_e_list"], name="p_e")
    assert_array_close(result.axes["tau"], bl["tau_list"], name="tau")
```

### 5.2 必须保持的回归

`tests/regression/test_physics_baseline.py` 全部 pass(因为 src/protocal.py 现在是 facade,旧 API 行为通过 Experiment 类间接验证)。

---

## 6. 风险与回滚

| # | 风险 | 缓解 |
|---|---|---|
| 6.1 | RamseyExperiment 的默认参数与 case 1 不完全一致,导致 p_e 偏差 | 严格按 src/protocal.py:62-66 复制 amplitude=0.001、frequency=0.01、rise=fall=10、center=100、noise_level=0 |
| 6.2 | TransientSensingExperiment 中 KernelEstimator 与旧 control_pulse.get_kernel 数值不等(浮点积分顺序) | 测试中允许 rtol=1e-5(略宽于 1e-6),并记录在 commit message |
| 6.3 | facade Protocal.evolve 返回元组形状错误,web_demo.py 解包失败 | 在 P2 PR 必须显式跑 web_demo.py;tests/integration/test_legacy_facade.py 验证返回元组 |
| 6.4 | sliding_measrement 的拼写改动后 src/protocal.py 中调用失效 | 保留方法名 `sliding_measrement`(拼写错误)在 Protocal 类中;新代码用 `SlidingMeasurementRunner` |
| 6.5 | KernelEstimator 把 stim_amplitude 默认值改成 auto_calibrate 后 transient_default 数值变化 | P2 阶段 `auto_calibrate=False`,保持 0.0215 默认值;auto_calibrate 是 P3 或更晚的事 |

**回滚**:
- 单 commit 回滚:任意 baseline 失败,revert 该 commit。
- Phase 级回滚:`git checkout master -- src/protocal.py` + 删除 `sqc/experiments/{rabi,ramsey,echo,transient}.py`,即可回到 P1 状态。

---

## 7. 输出物

### 7.1 文件清单(新增)

```
sqc/experiments/rabi.py
sqc/experiments/ramsey.py
sqc/experiments/echo.py                          (DiffEchoExperiment + EchoExperiment)
sqc/experiments/transient.py
sqc/experiments/legacy.py                        (可选:Protocal facade,若 src/protocal.py 太复杂)
sqc/simulation/runner.py                         (P1 是 stub,P2 完整实现)
sqc/hardware/readout.py                          (P1 是 ABC,P2 完整实现 + IQ_readout_legacy)
sqc/reconstruction/kernel.py                     (P1 仅有,P2 完整实现 KernelEstimator)
tests/unit/test_kernel_estimator.py
tests/unit/test_runners.py
tests/unit/test_iq_readout.py
tests/integration/test_ramsey_experiment.py
tests/integration/test_diff_echo_experiment.py
tests/integration/test_transient_experiment.py
tests/integration/test_legacy_facade.py          (验证 Protocal.evolve 返回元组形状)
```

### 7.2 文件清单(修改)

```
src/protocal.py                                  (Protocal 改为 facade,case 0/1/2/4 内化;case 3/5 留旧实现)
src/analysis.py                                  (Analysis.get_kernel 转发到 KernelEstimator)
sqc/control/pulse.py                             (Pulse.get_kernel / CompositePulse.get_kernel 转发)
```

### 7.3 接口快照

| 符号 | 来源 |
|---|---|
| `KernelEstimator(stim_amplitude=0.0215, stim_width=3.0, auto_calibrate=False)` | sqc.reconstruction.kernel |
| `KernelEstimator.estimate(pulse, qubit) -> (t_samples, kernel)` | 同上 |
| `MesolveRunner(options=None)` | sqc.simulation.runner |
| `SlidingMeasurementRunner(options=None)` | 同上 |
| `IQReadoutModel(tau=20, t_rabi=None)` | sqc.hardware.readout |
| `IdealProjectiveReadout()` | 同上 |
| `RabiExperiment(qubit, t_rabi=None, omega_d=None)` | sqc.experiments.rabi |
| `RamseyExperiment(qubit, flux_signal=None, omega_d=None, t_rabi, tau_list, t_global, phase1=0, phase2=0)` | sqc.experiments.ramsey |
| `DiffEchoExperiment(qubit, flux_signal=None, ...)` | sqc.experiments.echo |
| `TransientSensingExperiment(qubit, flux_signal=None, ...)` | sqc.experiments.transient |

---

## 8. 与 Track B 的协调点

| Track B 任务 | P2 影响 |
|---|---|
| _TODO 0.1 case 1 死代码清理 | 完成后,RamseyExperiment 不需要 mirror 死代码 |
| _TODO 0.2 kernel 自动校准 | 完成后,KernelEstimator 加 auto_calibrate=True 路径;**重生成 transient_default.pkl** |
| _TODO 0.3 LM 收敛 | 与 P2 无直接耦合,但建议先完成,以免 P3 阻塞 |
| _TODO 1.1 Cryoscope (case 6/7) | P2 留 case 5/6/7 在旧 protocal.py 中,P3 内化 |
| _TODO 1.2 瞬态标定 (case 8) | 同上 |

---

## 9. 完成确认清单

```
- [ ] sqc/experiments/{rabi,ramsey,echo,transient}.py 全部实现
- [ ] sqc/simulation/runner.py 完整实现 (MesolveRunner + SlidingMeasurementRunner)
- [ ] sqc/hardware/readout.py 完整实现 (IQReadoutModel + IQ_readout_legacy)
- [ ] sqc/reconstruction/kernel.py 完整实现 KernelEstimator
- [ ] sqc/control/pulse.py 中 get_kernel 转发到 KernelEstimator
- [ ] src/protocal.py:Protocal.evolve case 0/1/2/4 改为 facade
- [ ] src/protocal.py:Protocal.sliding_measrement 转发到 SlidingMeasurementRunner
- [ ] src/protocal.py:Protocal.evolve case 3/5 保持原状
- [ ] src/protocal.py:Calibration 类保持原状(P3 内化)
- [ ] src/analysis.py:Analysis.get_kernel 转发到 KernelEstimator
- [ ] tests/unit/ 至少新增 3 个文件
- [ ] tests/integration/ 至少新增 4 个文件
- [ ] pytest tests/unit -v 全部 pass
- [ ] pytest tests/integration -v 全部 pass
- [ ] pytest tests/regression -m regression -v 全部 pass(关键!)
- [ ] python web_demo.py Rabi/Ramsey/Transient 三个功能均跑通
- [ ] Simulation.ipynb 前 15 个 cell 可运行
```

---

**下一步**:[phase_3_handbook.md](phase_3_handbook.md) — Track B 成果(LM/Cryoscope/瞬态标定)内化。
