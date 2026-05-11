# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## Project overview

Quantum sensing simulation platform for superconducting Transmon qubits. The platform simulates time-dependent magnetic field sensing via qubit-frequency transduction, supporting multiple sensing protocols. All simulations use natural units (ħ=1) with Fock basis expansion via QuTiP.

## Workflow — before writing any code

1. **Restate your understanding.** What problem are you solving? What is the deliverable? Flag any assumptions you made. If you see a better technical approach, say so directly — the user will decide.
2. **Ask clarifying questions** (at most 3 at a time) until you are 100% confident about:
   - The real goal the user wants to achieve (not just what they said literally).
   - Unspoken constraints or preferences — tech stack, performance requirements, existing code that must remain compatible, parts of the system that must not be touched.
   - Your planned implementation — the core idea and why you chose it.
3. **Do not write code or modify any files** until the user explicitly signals approval (e.g. "可以开始", "go ahead", "proceed").

## Development commands

```bash
# Install dependencies (conda environment, per .vscode/settings.json)
pip install -r requirements.txt

# Core dependency check
python -c "import qutip, numpy, matplotlib; print('ok')"

# Verify project imports
python -c "from src.qubit import TransmonQubit; from src.signal import Signal; print('ok')"

# Run the Gradio web demo
python web_demo.py

# There is no test suite or linting configuration in this project.
```

## Architecture: data flow

```
Signal (signal.py)  →  qubit.qubit_in_mag(Signal)  →  Protocol.evolve(qubit)  →  Analysis methods
                         ↓                                    ↓                        ↓
                  Qubit frequency shifts            mesolve() with              Deconvolution /
                  under external flux               control pulses              numerical inversion
```

## Core modules (`src/`)

### [qubit.py](src/qubit.py) — Qubit physical models
- **`TransmonQubit`** — Central class. Holds EC, EJ, T1, T2, flux, state. Computes frequency and anharmonicity analytically. Key methods:
  - `qubit_in_mag(Phi_signal)` — accepts a `Signal` and pre-computes `self.freq_coeffs` and `self.H_list` (time-dependent Hamiltonian in QuTiP list format). This is the main interface for sensing simulations.
  - `qubit_under_mag(Phi_signal)` — older, slower per-timestep variant; creates a new qubit object at each time point. Prefer `qubit_in_mag()`.
  - `frequency_sensitivity(flux)` — numerical dω/dΦ via central difference.
  - Gate simulation: `simulate_gate()` with DRAG pulses, `ideal_gate()`.
- **`Cavity`** — Multi-mode cavity model with tensor-product operators.
- **`Coupled_System`** — Two qubits coupled via a tunable coupler to a multi-mode cavity. Used for two-qubit gate simulation (iSWAP, CZ).

### [signal.py](src/signal.py) — Magnetic flux signals
- **`Signal`** — Time-domain signal with 8 types (0=zero, 1=constant, 2=sinusoidal, 3=Gaussian, 4=asymmetric impulse, 5=double-peak, 6=basis-expanded, 7=complex wave-packet, 8=custom). Key property: `self.signal` (numpy array), `self.t_list`. Methods: `value_at(t)`, `truncate()`, `update_signal()`, `plot()`.
- **`CompositeSignal`** — Concatenates multiple Signals sequentially.

### [pulse.py](src/pulse.py) — Control pulses
- **`Pulse`** — Single control pulse. Accepts a `Signal` as the Rabi envelope (`Omega`). Builds Hamiltonian in either lab frame (`frame=0`) or rotating frame (`frame=1`), with or without RWA. Returns QuTiP list-format `[op, coeff_array]`.
- **`CompositePulse`** — Sequences multiple Pulses with proper time alignment.
- **Factory functions**: `create_pulse()`, `create_ramsey_pulse()`, `create_echo_pulse()`, `create_cpmg_pulse()`, `create_diff_echo_pulse()`, `create_cryoscope_pulse()`.
- **`get_kernel()`** on `Pulse`/`CompositePulse` computes the control kernel by perturbing with a narrow Gaussian stimulus at each time point.

### [protocal.py](src/protocal.py) — Sensing protocols
- **`Protocal`** (intentional spelling, used project-wide). `type` selects the protocol:
  - `0` Rabi, `1` Ramsey, `2` Differential echo, `3` CPMG (stub), `4` Transient field with sliding measurement, `5` Cryoscope.
- `evolve(qubit)` runs the protocol end-to-end: creates test signals, constructs Hamiltonians, calls `mesolve()`.
- `sliding_measrement()` — slides a control pulse across a signal, used by protocol type 4. Pre-computes qubit time series and Hamiltonians.
- `single_measurement()` — evaluates p_e at one delay time.
- **`Calibration`** — For frequency-flux and phase-flux calibration curves.
- **`IQ_readout()`** — I/Q demodulation by running two Ramsey sequences with π/2 phase offset.

### [analysis.py](src/analysis.py) — Data analysis and field reconstruction
- **`Analysis`** class with methods for extracting physical quantities from mesolve results:
  - `get_signal_from_ramsey_by_iq()` / `by_unwrap()` — Ramsey phase unwrapping and B-field derivation.
  - `get_signal_from_diff_echo()` — Direct formula B = -φ/(2k·κ·t_int).
  - `wiener_deconvolution()` — Linear Wiener deconvolution (FFT-based) of kernel from Δp measurements.
  - `hammerstein_wiener_deconvolution()` — Nonlinear block model: B → ω(B) → φ → Δp.
  - `numerical_inverse()` — Levenberg-Marquardt optimization using full density matrix simulation with basis function parameterization (B-spline, Fourier, Legendre). Includes adjoint-method Jacobian.
  - `get_signal_from_cryoscope()` — Phase-differential method with calibration-curve inversion.
- **Helper functions**: `generate_basis_functions()`, `basis_function_decomposition()`, `forward_simulation()`, `compute_jacobian()` (adjoint), `compute_jacobian_finite_difference()`, `levenberg_marquardt()`.

## Key conventions and gotchas

- **Class name**: `Protocal` is intentionally misspelled everywhere — do not "fix" it.
- **Units**: Frequencies/energies in GHz, time in ns, flux in Φ₀. ħ=1 throughout.
- **Hamiltonian format**: QuTiP list format `[[H0, coeffs], [H1, coeffs_array], ...]` is used pervasively for time-dependent Hamiltonians. `QobjEvo(H, tlist=t_list)` wraps these.
- **`qubit_in_mag()` must be called before simulation** whenever the signal changes — it pre-computes `freq_coeffs` and `H_list` used by `mesolve()`.
- **Kernel computation** uses a hardcoded stimulus amplitude (0.0215) and width (3 ns) in both `Pulse.get_kernel()` and `Analysis.get_kernel()`. These need manual adjustment when qubit parameters change.
- **LM inversion** (numerical_inverse) relies on `qubit_in_mag()` being called before each forward pass. The adjoint Jacobian (`compute_jacobian`) has known convergence issues vs. finite difference — see project TODO 0.3.

## Notebook and results

- [`Simulation.ipynb`](Simulation.ipynb) is the main experimental notebook for protocol testing and visualization.
- [`result/`](result/) contains output data and plotting scripts organized by experiment type (amplitude_scan, lambda_scan, convergence, basis_comparison, etc.).
- [`idea/`](idea/) contains design documents and implementation plans (markdown). The master TODO is [`idea/_TODO_master.md`](idea/_TODO_master.md).
- [`web_demo.py`](web_demo.py) is a Gradio web interface exposing protocols 0, 1, and 4.

---

## Refactoring activity (Track A) — hard rules

> **何时适用**:当你的任务涉及 [`idea/refactor/`](idea/refactor/) 中的任何文档,或者需要在 [`sqc/`](sqc/) / [`src_mirror/`](src_mirror/) 目录下增删改文件时,本节规则**全部生效**。其他常规开发(在 `src/` 内做功能扩展、写 notebook 等)不受本节约束。

### R1. 不可违反的硬约束:`src/` 永远不变

- Track A(重构)**永远不修改、不删除、不重命名** `src/` 下的任何文件。
- 任何重构相关 PR 中,`git diff master -- src/` **必须为空**。如果你发现自己即将修改 `src/`,**立即停手并报告**。
- 所有 facade、wrapper、镜像导出代码,**一律放在新建的 `src_mirror/` 目录下**(与 `src/` 平级)。
- 详细原因见 [`idea/refactor/_refactor_plan.md`](idea/refactor/_refactor_plan.md) §8.8。

### R2. 启动顺序(每个 phase 子代理或新 session 开始时必读)

按顺序读以下文件,**不准跳过**:
1. [`idea/refactor/_refactor_plan.md`](idea/refactor/_refactor_plan.md) §0、§8、§13、§15(全局认知 + 兼容层契约 + 物理对应)
2. [`idea/refactor/_handoff_state.md`](idea/refactor/_handoff_state.md)(上一阶段交接状态)
3. 当前要执行的 [`idea/refactor/phase_N_handbook.md`](idea/refactor/)(具体任务清单)
4. 物理参考(按需查阅):[`idea/refactor/Gao 等 - 2021 - Practical Guide for Building Superconducting Quantum Devices.pdf`](idea/refactor/) — 主方案 §15 给出按章节定位

### R3. 测试纪律

- 任何代码改动后:`pytest tests/regression -m regression` 必须 100% 通过。
- 任何 commit 前:`git diff --quiet master -- src/` 必须返回 0(R1 的自动化检查)。
- 物理回归 baseline 容限:`rtol=1e-6, atol=1e-9`。**不准为了让测试通过而放宽这个容限**;如果数值偏移,先调查根因。
- 内化新功能(如 LM、Cryoscope)时,新增 baseline pickle 必须在 `tests/baselines/` 下,生成命令进入 `tests/regression/generate_baselines.py`。

### R4. 命名与拼写

- `Protocal`(故意错拼)在 `src/` 中**永远不动**;`src_mirror/protocal.py:Protocal` 同样保留这个拼写以保持 API 一致。
- `sqc/` 中所有新代码用 `Protocol` / `Experiment` / `Calibration` 等**正确拼写**的类名。
- `sliding_measrement`(故意错拼,`src/protocal.py:252`)在 `src_mirror/` facade 中保留;`sqc/simulation/runner.py:SlidingMeasurementRunner` 用正确拼写。

### R5. 子代理执行模式

本项目用 **subagent-per-phase** 方式执行(详见主方案 §13 与 [`.Codex/agents/refactor-phase-executor.md`](.Codex/agents/refactor-phase-executor.md))。如果你正以子代理身份运行:
- 你**不能中途询问用户**。所有决策点要么自行决定(选最保守),要么 abort 并在 handoff 中标 `DECISION_NEEDED:`。
- 拿不准时**绝不放宽 R1/R3**;宁可 abort 也不污染主分支。
- 完成后**必须更新** `idea/refactor/_handoff_state.md`,这是与下一个子代理的唯一交接载体。

### R6. 物理参数 sanity check

任何 PR 引入新的 qubit/cavity 参数,先核对论文 §III.B 推荐范围(详见主方案 §15.7):
- `EJ/h ∈ [10, 25] GHz`,`EC/h ∈ [160, 400] MHz`,`EJ/EC ∈ [40, 80]`
- `f_01 ∈ [4, 8] GHz`,`α/h ∈ [200, 300] MHz`
- 偏离推荐范围**必须在 commit message 中说明物理动机**。

### R7. Track B(功能开发)与 Track A(重构)的区分

- **Track B**(在 `src/` 内做功能开发,如实现 Cryoscope case 6/7、瞬态标定 case 8 等)**不受本节 R1 约束** — Track B 可以自由修改 `src/`,这是 `src/` 的正当演化。
- **Track A**(本次重构)受 R1–R6 全部约束。
- 如何分辨:看 PR 是否动到 `sqc/` / `src_mirror/` / `idea/refactor/`。动了就是 Track A,不动就是 Track B。
- Track B 的进度跟踪在 [`idea/_TODO_master.md`](idea/_TODO_master.md);Track A 的进度跟踪在 [`idea/refactor/_handoff_state.md`](idea/refactor/_handoff_state.md)。
