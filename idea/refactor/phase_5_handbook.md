# Phase 5 Handbook — 多 qubit Z-crosstalk + Cavity 表征 论文级 demo

> 前置阅读:[_refactor_plan.md](_refactor_plan.md) §1.3、§4–§6、§15.4、§15.5,[phase_4_handbook.md](phase_4_handbook.md)  
> 物理参考(本地 PDF):[`./Gao 等 - 2021 - Practical Guide for Building Superconducting Quantum Devices.pdf`](./Gao%20%E7%AD%89%20-%202021%20-%20Practical%20Guide%20for%20Building%20Superconducting%20Quantum%20Devices.pdf) §V.D(双比特门三类),§V.E + Fig. 14 + Fig. 16(残余 ZZ Eq. 75–85),§V.F + Fig. 17(cavity 表征:number splitting/Ramsey revival/Wigner Eq. 87–89)  
> 估计工时:5–7 天(基础)+ 3–5 天(可选 cavity 表征扩展)  
> 触发条件:Phase 4 完成,DistortionModel + ControlLine + PredistortionDesigner 已稳定  
> 完成标志:`TransferMatrix` 完整实现 + 双 qubit `ChipTopology` + Z-crosstalk demo workflow + 双 qubit baseline;**(可选)** Cavity number splitting / Ramsey revival / Wigner tomography 三件套;**`src/` 仍未被 Track A 修改**

---

## 1. 目标

把项目推进到展示**核心应用价值**的阶段:多 qubit 串扰建模与补偿。本 phase 主要新增,但需要把现有 `Coupled_System` 升级为**通用 ChipTopology** + 双 qubit Z-crosstalk demo。

**物理设定**:在 qubit A 的 Z 线上施加 flux pulse,通过电感/容性耦合,qubit B 也感受到一个寄生磁通 Φ_B(t) = H_BA(ω) · V_A(ω)。本 phase 演示如何:
1. 在 qubit B 上**重建** Φ_B(t)。
2. 提取 H_BA(ω)。
3. 设计补偿,使 qubit B 的寄生相位下降。

主要交付:
1. `TransferMatrix` 完整实现(P1 只是 dataclass 占位)。
2. `ChipTopology` — 多 qubit + coupler + cavity + control lines 的拓扑容器。
3. `ZCrosstalkWorkflow` — 端到端串扰提取 + 补偿 demo。
4. 双 qubit `Coupled_System` 在新框架下的兼容封装(可选)。
5. Notebook 级 demo:展示 H_BA(t) 重建、补偿前后寄生相位对比。

---

## 2. 前置条件

| # | 条件 | 验证 |
|---|---|---|
| 2.1 | Phase 4 完成 | DistortionModel + ControlLine + Predistortion 工作 |
| 2.2 | (推荐) Track B `_TODO_master.md 4.x` 中相关多 qubit 模型已设想好物理参数 | 双 qubit 间耦合 g_AB 选合理值 |

---

## 3. 任务清单

### 3.1 ChipTopology 完整实现

**文件**:`sqc/devices/chip.py`

```python
# sqc/devices/chip.py (P5 增量)
from __future__ import annotations
from dataclasses import dataclass, field

import numpy as np
from qutip import tensor, qeye, destroy, basis, Qobj

from sqc.devices.transmon import TransmonQubit
from sqc.devices.resonator import Resonator
from sqc.devices.base import Device


@dataclass
class ChipTopology:
    """Multi-qubit chip topology: qubits + couplers + resonators + control lines.
    
    For backward compat, CoupledSystem (P1's port of Coupled_System) is a
    special case: 2 qubits + 1 multimode resonator.
    """
    qubits: list[TransmonQubit]
    resonators: list[Resonator] = field(default_factory=list)
    couplings: dict[tuple[str, str], float] = field(default_factory=dict)
    # ^ key: (qubit_name, resonator_name) → coupling strength g (rad/ns)
    control_lines: dict[str, "ControlLine"] = field(default_factory=dict)
    transfer_matrix: "TransferMatrix | None" = None
    
    def hilbert_dim(self) -> int:
        d = 1
        for q in self.qubits:
            d *= q.n_levels
        for r in self.resonators:
            for n in r.n_levels:
                d *= n
        return d
    
    def lift_qubit_op(self, op: Qobj, qubit_name: str) -> Qobj:
        """Embed a single-qubit operator into the chip Hilbert space."""
        ops = []
        for q in self.qubits:
            if q.name == qubit_name:
                ops.append(op)
            else:
                ops.append(qeye(q.n_levels))
        for r in self.resonators:
            for n in r.n_levels:
                ops.append(qeye(n))
        return tensor(*ops)
    
    def hamiltonian_static(self) -> Qobj:
        """Sum of static qubit + resonator + coupling Hamiltonians."""
        ...
    
    def collapse_operators(self) -> list[Qobj]:
        ...
    
    @classmethod
    def from_legacy_coupled_system(cls, coupled_system) -> "ChipTopology":
        """Adapter: build ChipTopology from legacy Coupled_System."""
        ...
```

### 3.2 TransferMatrix 完整实现

**文件**:`sqc/hardware/transfer_matrix.py`

```python
# sqc/hardware/transfer_matrix.py (P5 完整版)
from __future__ import annotations
from dataclasses import dataclass, field

import numpy as np

from sqc.control.waveform import Waveform
from sqc.control.flux_signal import FluxSignal


@dataclass
class TransferMatrix:
    """Multi-qubit Z-line transfer matrix H_ji(ω).
    
    Phi_j(ω) = Σ_i H_ji(ω) V_i(ω)
    
    elements[(target, source)] is a complex array on frequency_axis.
    """
    elements: dict[tuple[str, str], np.ndarray] = field(default_factory=dict)
    frequency_axis: np.ndarray = field(default_factory=lambda: np.zeros(0))
    time_axis: np.ndarray | None = None
    
    @property
    def sources(self) -> list[str]:
        return sorted({s for (_, s) in self.elements.keys()})
    
    @property
    def targets(self) -> list[str]:
        return sorted({t for (t, _) in self.elements.keys()})
    
    def H_ji(self, target: str, source: str) -> np.ndarray:
        return self.elements[(target, source)]
    
    def diagonal(self) -> dict[str, np.ndarray]:
        """Self-response H_ii."""
        return {t: self.elements[(t, t)] for t in self.targets if (t, t) in self.elements}
    
    def off_diagonal(self) -> dict[tuple[str, str], np.ndarray]:
        """Cross-talk H_ji (j ≠ i)."""
        return {(t, s): H for (t, s), H in self.elements.items() if t != s}
    
    def apply(self, source_voltages: dict[str, Waveform]) -> dict[str, FluxSignal]:
        """Compute on-chip flux at each target, given AWG voltages on each source."""
        # FFT each source, multiply by appropriate H_ji, sum, IFFT
        targets = self.targets
        result = {}
        for j in targets:
            t_list = None
            phi_j_omega = None
            for i, V_i in source_voltages.items():
                if (j, i) not in self.elements:
                    continue
                if t_list is None:
                    t_list = V_i.t_list
                V_omega = np.fft.fft(V_i.samples)
                # Interpolate H_ji onto V_omega's frequency axis
                omega_grid = 2 * np.pi * np.fft.fftfreq(len(V_i.samples),
                                                         t_list[1] - t_list[0])
                H_ji = self._interpolate(self.elements[(j, i)], omega_grid)
                if phi_j_omega is None:
                    phi_j_omega = H_ji * V_omega
                else:
                    phi_j_omega += H_ji * V_omega
            if phi_j_omega is None:
                continue
            phi_j = np.real(np.fft.ifft(phi_j_omega))
            result[j] = FluxSignal(
                type=8, t_list=t_list, signal=phi_j,
            )
        return result
    
    def _interpolate(self, H_arr: np.ndarray, omega_target: np.ndarray) -> np.ndarray:
        re = np.interp(omega_target, self.frequency_axis, H_arr.real)
        im = np.interp(omega_target, self.frequency_axis, H_arr.imag)
        return re + 1j * im
    
    @classmethod
    def from_dc_matrix(cls, dc_matrix: np.ndarray, source_names: list[str],
                        target_names: list[str]) -> "TransferMatrix":
        """Build a frequency-flat TransferMatrix from a DC crosstalk matrix.
        
        dc_matrix[j, i] = H_ji(ω=0) for all ω.
        """
        omega = np.linspace(-np.pi, np.pi, 1024)
        elements = {}
        for j, t in enumerate(target_names):
            for i, s in enumerate(source_names):
                elements[(t, s)] = dc_matrix[j, i] * np.ones_like(omega)
        return cls(elements=elements, frequency_axis=omega)
```

### 3.3 ZCrosstalkWorkflow

**文件**:`sqc/workflows/z_crosstalk.py`

```python
# sqc/workflows/z_crosstalk.py
from __future__ import annotations
from dataclasses import dataclass, field

import numpy as np

from sqc.devices.transmon import TransmonQubit
from sqc.devices.chip import ChipTopology
from sqc.control.waveform import Waveform
from sqc.control.flux_signal import FluxSignal
from sqc.hardware.control_line import ControlLine
from sqc.hardware.distortion import DistortionModel
from sqc.hardware.transfer_matrix import TransferMatrix
from sqc.experiments.transient import TransientSensingExperiment
from sqc.experiments.cryoscope import CryoscopeExperiment
from sqc.reconstruction.wiener import WienerReconstruction
from sqc.calibration.predistortion import PredistortionDesigner
from sqc.workflows.base import Workflow


@dataclass
class ZCrosstalkWorkflow(Workflow):
    """End-to-end Z-crosstalk demo:
    
    Setup: 2-qubit chip. Qubit A has its own Z line. Apply flux pulse to A.
    Qubit B (parked at max-sensitivity point) feels parasitic flux.
    
    Steps:
    1. Apply flux pulse on A's Z line through ground-truth crosstalk.
    2. Reconstruct Φ_B(t) from B via TransientSensingExperiment + Wiener.
    3. Extract H_BA(ω) by deconvolving Φ_B(t) with the input pulse.
    4. Design compensation pulse on A's Z line so that B's parasitic flux is cancelled.
    5. Verify: with compensation, B's parasitic phase φ_B is reduced.
    """
    chip: ChipTopology
    flux_pulse_on_A: Waveform
    true_transfer_matrix: TransferMatrix       # ground truth
    qubit_A_name: str = "QA"
    qubit_B_name: str = "QB"
    
    def run(self) -> dict:
        # 1. Apply ground-truth distortion + crosstalk
        source_voltages = {self.qubit_A_name: self.flux_pulse_on_A}
        on_chip_fluxes = self.true_transfer_matrix.apply(source_voltages)
        phi_A = on_chip_fluxes[self.qubit_A_name]
        phi_B_true = on_chip_fluxes[self.qubit_B_name]
        
        # 2. Reconstruct Φ_B(t) via transient sensing on B
        qubit_B = next(q for q in self.chip.qubits if q.name == self.qubit_B_name)
        # Park B at max-sensitivity point
        qubit_B = TransmonQubit(
            EC=qubit_B.EC, EJ=qubit_B.EJ_0, T1=qubit_B.T1, T2=qubit_B.T2,
            flux=qubit_B.spec().optimal_work_point(), n_levels=qubit_B.n_levels,
            name=qubit_B.name,
        )
        exp = TransientSensingExperiment(
            qubit=qubit_B, flux_signal=phi_B_true,
        )
        result = exp.run()
        
        # Use Wiener to reconstruct Φ_B
        wiener = WienerReconstruction(lambda_reg=1e-3)
        phi_B_reconstructed = wiener.reconstruct(
            result, kernel=result.data["kernel"],
            dt=phi_B_true.t_list[1] - phi_B_true.t_list[0],
        )
        
        # 3. Extract H_BA(ω) by frequency-domain deconvolution
        omega = 2 * np.pi * np.fft.fftfreq(
            len(self.flux_pulse_on_A.samples),
            self.flux_pulse_on_A.t_list[1] - self.flux_pulse_on_A.t_list[0],
        )
        Phi_B_omega = np.fft.fft(phi_B_reconstructed.samples)
        V_A_omega = np.fft.fft(self.flux_pulse_on_A.samples)
        H_BA_estimated = Phi_B_omega / (V_A_omega + 1e-12)
        
        # 4. Design compensation: target Φ_B should be 0
        # Compensation pulse on A's line: V_A_comp such that
        # Φ_B(comp) = -Φ_B(uncompensated)
        # i.e., V_comp_A = -Φ_B_estimated / H_BA_estimated
        # But this requires also accounting for H_AA, see below.
        # Simplified: assume H_AA ≈ 1, so compensation is straightforward.
        designer = PredistortionDesigner(method="frequency_inverse",
                                          regularization=1e-3)
        # Build the inverse matrix (only for the BA off-diagonal element)
        ...
        
        # 5. Verify by applying with compensation
        ...
        
        return {
            "phi_A": phi_A,
            "phi_B_true": phi_B_true,
            "phi_B_reconstructed": phi_B_reconstructed,
            "H_BA_estimated": H_BA_estimated,
            "H_BA_true": self.true_transfer_matrix.H_ji(
                self.qubit_B_name, self.qubit_A_name
            ),
            "fit_error_dB": self._fit_error(H_BA_estimated, ...),
            # Compensation results
            "compensation_pulse": ...,
            "phi_B_after_compensation": ...,
            "parasitic_phase_uncompensated": ...,
            "parasitic_phase_compensated": ...,
            "compensation_factor": ...,
        }
    
    @staticmethod
    def _fit_error(H_est, H_true) -> float:
        """RMS dB error of estimated vs true transfer function."""
        return 20 * np.log10(np.abs(H_est - H_true).mean() / np.abs(H_true).mean())
```

### 3.4 双 qubit Notebook demo (可选,不影响验收)

**新建**:`Simulation_z_crosstalk.ipynb` 或在 `Simulation.ipynb` 增加新 section。

包含 cell:
1. 构造 2-qubit ChipTopology + ground-truth TransferMatrix。
2. 跑 ZCrosstalkWorkflow。
3. 画图:phi_A、phi_B(真实)、phi_B(重建)、H_BA(频域真 vs 估)、补偿前后的寄生相位。

### 3.5 (可选扩展) Cavity 表征三件套 — 对应 Gao 2021 §V.F

如果未来要把项目扩展到 cavity-based bosonic qubit(论文 §VI.B),P5 可顺手实现 cavity 表征三件套。**这部分与磁通传感主线正交,优先级低,但实现门槛低,且物理上漂亮。**

#### 3.5.1 NumberSplittingExperiment

参考论文 §V.F + Fig. 17(b)。在 cavity 中位移到 |β⟩,然后用 transmon 的 spectroscopy 看到一系列等间距 χ 峰(每个对应一个 Fock 态)。物理:由色散耦合 $H_{int} = -(\chi/2) a^\dagger a \sigma_z$,transmon 频率被偏移 $-n \chi$,$n$ 是 cavity photon 数。

```python
# sqc/experiments/cavity_spectroscopy.py
@dataclass
class NumberSplittingExperiment(Experiment):
    """Probe transmon spectrum when cavity is displaced; reveals
    photon-number-resolved peaks separated by chi.
    
    Replicates Gao 2021 §V.F Fig. 17(b).
    """
    qubit: TransmonQubit
    cavity: Resonator
    g: float                         # cavity-qubit coupling (rad/ns)
    beta: complex = 1.0              # displacement amplitude
    omega_scan: np.ndarray | None = None
    
    def run(self) -> ExperimentResult:
        # 1. Displace cavity to |beta>
        # 2. Spectroscopy on transmon while measuring P_e
        # 3. Output: spectrum showing peaks at omega_q - n*chi for n=0,1,2,...
        ...
```

**验收**:输出谱中相邻峰间距 ≈ $\chi = 2 E_C (g/\Delta)^2$(Gao 2021 Eq. 35),误差 < 5%。

#### 3.5.2 RamseyRevivalExperiment

参考论文 §V.F Eq. (87) $P_e = \frac{1}{2}\{1 + e^{-2|\beta|^2 \sin^2(\chi t/2)} \cos(|\beta|^2 \sin\chi t)\}$。在 transmon 上做 Ramsey,但 cavity 被显式放在大相干态。$P_e$ 随时间出现 revival(在 $t = 2n\pi/\chi$ 处)。

```python
@dataclass
class RamseyRevivalExperiment(Experiment):
    """Transmon Ramsey with cavity in coherent state |beta>.
    Revival times reveal chi precisely.
    
    Replicates Gao 2021 §V.F Fig. 17(d).
    """
    qubit: TransmonQubit
    cavity: Resonator
    beta: float = 2.0                # |beta|, "large" coherent state
    t_max: float = 1000.0            # ns, scan range
    n_t: int = 200
    
    def run(self) -> ExperimentResult: ...
```

**验收**:revival 峰间距 = $2\pi/\chi$,提取的 $\chi$ 与 number splitting 一致(误差 < 5%)。

#### 3.5.3 WignerTomographyWorkflow

参考论文 §V.F Eq. (89) $W(\alpha) = (2/\pi) \mathrm{Tr}[D_\alpha^\dagger \rho D_\alpha P]$。通过对 cavity 做扫描位移 $D_\alpha$ + 奇偶测量(用 transmon 作 ancilla),得到 Wigner 函数。

```python
# sqc/workflows/wigner_tomography.py
@dataclass
class WignerTomographyWorkflow(Workflow):
    """Wigner-function tomography via parity measurement.
    
    Implements Gao 2021 §V.F Fig. 17(e):
    1. Prepare cavity state (e.g., Fock |1> via SNAP/OCT)
    2. Displace cavity by alpha (scan grid)
    3. Map photon-parity onto transmon Y/2 - C_phi(t=π/χ) - Y/2
    4. Single-shot transmon readout → P(even) - P(odd) = parity
    5. W(alpha) = (2/π) parity
    """
    qubit: TransmonQubit
    cavity: Resonator
    initial_cavity_state: Qobj         # e.g., Fock |1> or coherent
    alpha_grid: np.ndarray             # 2D complex grid
    
    def run(self) -> dict:
        # ... per-alpha parity measurement ...
        # output: 2D Wigner image
        ...
```

**验收**:对 |0⟩ 输入,Wigner 函数等于二维高斯;对 |1⟩,中心点为负值(展示 photon Fock 态的非经典性)。

**注**:这三个实验在论文 §V.F 中是 cavity-based logical qubit 的标准表征工具。本项目当前不做 logical qubit,但实现它们的成本低(每个 < 200 行),**是论文级 demo 的天然候选**。如果时间允许,推荐在 P5 完成后立刻做。

### 3.6 (可选) Coupled_System 兼容包装

如果 P5 实质上替代了 `src/qubit.py:Coupled_System`,可在 `sqc/devices/chip.py` 中加:

```python
# sqc/devices/chip.py
class CoupledSystem(ChipTopology):
    """Backward-compat: 2 qubits + 1 multimode resonator.
    
    Direct port of src/qubit.py:Coupled_System.
    """
    
    def __init__(self, qubit1, qubit2, cavity):
        super().__init__(
            qubits=[qubit1, qubit2],
            resonators=[cavity],
        )
        self.qubit1 = qubit1
        self.qubit2 = qubit2
        self.cavity = cavity
        # ... port of all methods from src/qubit.py:Coupled_System ...
```

P1 已经把 Coupled_System 迁移过去了,P5 主要是给它加 ChipTopology 父类。

---

## 4. 验收标准

### 4.1 TransferMatrix 基本功能

```bash
python -c "
import numpy as np
from sqc.hardware.transfer_matrix import TransferMatrix
from sqc.control.waveform import Waveform

tm = TransferMatrix.from_dc_matrix(
    dc_matrix=np.array([[1.0, 0.05], [0.03, 1.0]]),
    source_names=['QA', 'QB'],
    target_names=['QA', 'QB'],
)
t = np.linspace(0, 100, 1000)
V_A = Waveform(t_list=t, samples=np.where(t > 10, 1.0, 0.0))
V_B = Waveform(t_list=t, samples=np.zeros_like(t))
fluxes = tm.apply({'QA': V_A, 'QB': V_B})
phi_B = fluxes['QB'].samples[-1]
print(f'phi_B at end = {phi_B:.4f} (should be ~0.03)')
"
```

期望:`phi_B at end ≈ 0.03`(±0.001)。

### 4.2 端到端 Z-crosstalk demo

```bash
python -c "
... 跑 ZCrosstalkWorkflow,断言:
- |H_BA_estimated - H_BA_true| / |H_BA_true| < 0.1 (10%)
- compensation_factor > 5 (寄生相位降至原来的 1/5 以下)
"
```

### 4.3 物理回归

```bash
pytest tests/regression -m regression -v
```

P5 引入新 baseline:
- `z_crosstalk_default.pkl`(ZCrosstalkWorkflow 输出 + 关键 metric)

### 4.4 web_demo

(可选)增加 "Z-Crosstalk" tab,展示串扰提取 + 补偿 demo。

---

## 5. 测试要求

### 5.1 必须新增的测试

| 文件 | 用例 |
|---|---|
| `tests/unit/test_chip_topology.py` | hilbert_dim 计算、lift_qubit_op tensor 顺序、hamiltonian_static 维度 |
| `tests/unit/test_transfer_matrix.py` | from_dc_matrix → apply 在 DC 极限下等价于矩阵乘法;H_ji / diagonal / off_diagonal 正确;FFT/IFFT 卷积一致性 |
| `tests/integration/test_z_crosstalk_workflow.py` | 完整 ZCrosstalkWorkflow,验证 H_BA 提取误差 < 10%,补偿因子 > 5 |
| `tests/regression/test_z_crosstalk_baseline.py` | 与 baseline 比对 |

### 5.2 必须保持的回归

所有已有 baseline 仍 pass。

---

## 6. 风险与回滚

| # | 风险 | 缓解 |
|---|---|---|
| 6.1 | 双 qubit Hilbert 空间维度爆炸(3³ × 多模 cavity = 数百维) | 限制 cavity n_levels=2、qubit n_levels=3,Hilbert dim ≤ 36 |
| 6.2 | TransferMatrix.apply 中 FFT 边界效应导致 phi_B 失真 | 在 ZCrosstalkWorkflow 加 padding 选项 |
| 6.3 | qubit B 在最优灵敏度点引入额外退相干 | T2 noise 在 ZCrosstalkWorkflow 关闭(设大数),让物理结果纯净;后续可加噪 |
| 6.4 | H_BA 在 |V_A(ω)| ≈ 0 处估计错 | regularization 加在分母:H_BA = (Phi_B · conj(V_A)) / (|V_A|² + λ²) |
| 6.5 | 补偿后 H_AA 不为 1,导致 V_A 自身变化 | demo 假设 H_AA = 1(已校准);若不为 1,先用 P4 PredistortionDesigner 处理 H_AA |

**回滚**:本 phase 全部新增,删除 sqc/devices/chip.py 中 ChipTopology 类 + sqc/hardware/transfer_matrix.py 完整版 + sqc/workflows/z_crosstalk.py 即可。CoupledSystem 兼容包装回退到 P1 版本。

---

## 7. 输出物

### 7.1 文件清单(新增)

```
sqc/devices/chip.py                     (扩展 ChipTopology + CoupledSystem 父类化)
sqc/hardware/transfer_matrix.py         (P1 占位 → P5 完整)
sqc/workflows/z_crosstalk.py
tests/unit/test_chip_topology.py
tests/unit/test_transfer_matrix.py
tests/integration/test_z_crosstalk_workflow.py
tests/regression/test_z_crosstalk_baseline.py
tests/baselines/z_crosstalk_default.pkl
Simulation_z_crosstalk.ipynb            (可选)
```

### 7.2 文件清单(修改)

```
sqc/devices/chip.py                     (CoupledSystem 加 ChipTopology 父类)
tests/regression/generate_baselines.py  (添加 z_crosstalk_default 生成)
docs/architecture.md                    (可选:加 Z-crosstalk 章节)
```

### 7.3 接口快照

| 符号 | 来源 |
|---|---|
| `ChipTopology(qubits, resonators, couplings, control_lines, transfer_matrix)` | sqc.devices.chip |
| `ChipTopology.lift_qubit_op(op, qubit_name) -> Qobj` | 同上 |
| `ChipTopology.hamiltonian_static() -> Qobj` | 同上 |
| `TransferMatrix(elements, frequency_axis, time_axis)` | sqc.hardware.transfer_matrix |
| `TransferMatrix.apply(source_voltages) -> dict[str, FluxSignal]` | 同上 |
| `TransferMatrix.H_ji(target, source) -> np.ndarray` | 同上 |
| `TransferMatrix.diagonal() / off_diagonal()` | 同上 |
| `TransferMatrix.from_dc_matrix(dc_matrix, source_names, target_names)` | 同上 |
| `ZCrosstalkWorkflow(chip, flux_pulse_on_A, true_transfer_matrix, qubit_A_name, qubit_B_name)` | sqc.workflows.z_crosstalk |
| `ZCrosstalkWorkflow.run() -> dict` | 同上 |

---

## 8. 完成确认清单

```
- [ ] sqc/devices/chip.py 扩展 ChipTopology 与 CoupledSystem 父类化
- [ ] sqc/hardware/transfer_matrix.py 完整实现
- [ ] sqc/workflows/z_crosstalk.py 完整实现
- [ ] tests/baselines/z_crosstalk_default.pkl 生成
- [ ] tests/unit 新增 2 个文件,各 ≥ 3 个用例
- [ ] tests/integration 新增 1 个文件
- [ ] tests/regression 新增 1 个文件
- [ ] pytest tests/unit -v 全部 pass
- [ ] pytest tests/integration -v 全部 pass
- [ ] pytest tests/regression -m regression -v 全部 pass
- [ ] ZCrosstalkWorkflow.run() 输出 H_BA 估计误差 < 10%
- [ ] ZCrosstalkWorkflow.run() 补偿因子 > 5
- [ ] (可选) Simulation_z_crosstalk.ipynb 端到端跑通
- [ ] (可选) web_demo.py 加 Z-Crosstalk tab 并跑通
```

---

## 9. 后续(P5 之后,本方案不覆盖)

P5 完成后,本方案描述的全栈重构基本结束。后续可能的扩展(对应 _TODO_master.md 阶段四):

| # | 任务 | 备注 |
|---|---|---|
| ext.1 | 噪声集成 | 已存设计文档 idea/TODO_noise.md,可作为下一份 handbook |
| ext.2 | 并行优化 sliding_measrement | 已存设计文档 idea/parallel_optimization_plan.md |
| ext.3 | CPMG 频谱重建 | 已存理论 idea/_sensing theory.md §3 |
| ext.4 | Volterra 二阶非线性核 | 已存 idea/nonlinear_response_extension.md |
| ext.5 | 压缩感知采样 | 已存 idea/innovation_proposals.md |
| ext.6 | 全局 Protocal 拼写修正 | 影响 web_demo + notebook,需配合 phase 6+ |

每一项扩展若需要,可按本方案的 phase handbook 模板新写一份独立文档。

---

## 10. 全方案完结确认

完成 P5 时,执行者应在主项目 README 或 idea/refactor/_refactor_plan.md 末尾追加一段:

```
## 重构完结记录

- Phase 0 完成:commit <SHA>,date <YYYY-MM-DD>
- Phase 1 完成:commit <SHA>,date <YYYY-MM-DD>
- Phase 2 完成:commit <SHA>,date <YYYY-MM-DD>
- Phase 3 完成:commit <SHA>,date <YYYY-MM-DD>
- Phase 4 完成:commit <SHA>,date <YYYY-MM-DD>
- Phase 5 完成:commit <SHA>,date <YYYY-MM-DD>

总新增代码:~7500 行(估)
总测试用例:~80 个
物理回归 baseline:8 个
src/ 镜像策略:永久维护
```

至此,Sensing-Project 已具备:
- 真实 cQED 全栈架构
- 三大主线(波形重建/qubit 标定/波形预失真)的端到端实现
- 物理结果不变的回归测试保护
- 可扩展到多 qubit、多 control line、多失真源的 demo 框架

后续科研工作可以直接在 sqc/ 中扩展,旧 src/ 永远可工作。
