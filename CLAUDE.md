# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

# Run the v2.0 visualization demo
python web_demo_v2.py

# Run tests (pytest)
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m pytest tests/ -v
```

### Conda environment path

The project runs in conda env `qutip-env` at the path below. Use this **full Python path** in all bash commands (`python` may not be in PATH):

```
C:\Users\21034\anaconda3\envs\qutip-env\python.exe
```

Example usage:

```bash
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m pytest tests/ -v
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m tests.regression.generate_baselines
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -c "from sqc.config import CONFIG; print(CONFIG.awg.dt)"
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

本项目用 **subagent-per-phase** 方式执行(详见主方案 §13 与 [`.claude/agents/refactor-phase-executor.md`](.claude/agents/refactor-phase-executor.md))。如果你正以子代理身份运行:
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

### R8. 重大改动前的安全协议

涉及"重大改动"时,按以下顺序执行,**确保任何时刻都能回溯到改前状态**:

1. **预先 commit 安全点**:执行 `git add -A && git commit -m "safety checkpoint: before XXX"`。如有未 commit 改动,**先 commit** 才能开始重大改动;有 untracked 文件就先 `git add` 或 `.gitignore`。
2. **告知用户回溯路径**:明确告诉用户"如果出问题可 `git reset --hard <SHA>` 回退到此点",并给出具体 SHA。
3. **改动过程中持续测试**:每改完一个模块就跑相应单元测试,**不要**累积 10 个改动再一起测——出错难以定位。出错时立刻定位,不前进到下一个改动。
4. **改完整体测试**:跑 `pytest tests/ -v` 完整套件,确认全绿才声明完成。任何测试失败必须先解决,不能"先 commit 再说"。
5. **出文档**:如改动影响 API/约定/默认值,**必须**更新 `docs/architecture.md` 对应章节(详见 R11)。

**什么算"重大改动"**——满足任意一条:
- 修改 5 个以上 `sqc/` 文件
- 改变任意 baseline pickle 的数值
- 修改 `sqc.config.CONFIG` 默认值
- 引入新的全局约定(如时间轴推导规则)
- 改变实验类/重建算法的默认行为
- 任何会让现有测试失败的改动

**反例**:仅改 docstring、仅加新文件不动旧文件、仅修 typo 等不算重大改动,可直接 commit。

### R9. 全局配置规范(v2.0+)

所有时间轴**必须**从 `sqc.config.CONFIG.awg.dt` 派生:

- 使用 `np.arange(start, stop, dt)`,**禁止** `np.linspace(start, stop, N)` 用于时间轴
- 模块内部辅助常量用 `_` 前缀(如 `_GT = CONFIG.awg.dt`、`_gap(duration)`),不在 `__all__` 中
- 外部代码用 `CONFIG.pulse.t_rabi.copy()` / `CONFIG.pulse.t_global.copy()` / `CONFIG.pulse.make_time(start, end)`
- 非时间轴的离散数组(如 `h_list = linspace(-0.03, 0.03, 21)` 的 flux 扫描点)可用 `linspace`——只有"时间"参数受此约束
- 修改 `CONFIG` 默认值前先查 `docs/architecture.md` §10 看影响范围
- 修改 `CONFIG` 默认值后**必须**重新生成 baseline:`"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m tests.regression.generate_baselines`
- 重新生成 baseline 必须在 commit message 显式说明:`"intentional physics change: regenerated baselines due to ..."`

### R10. 依赖与配置文件政策

- 只用以下依赖:**numpy / scipy / qutip / matplotlib / dataclasses(标准库)/ pytest(测试)/ gradio(demo)**
- **不引入** pydantic / attrs / yaml / toml / click / pyyaml 等新依赖,除非用户明确批准
- **不引入** yaml / toml / .ini 配置文件——所有参数通过 Python 对象传递(`sqc.config.CONFIG`)
- 引入新依赖时必须:
  1. 在 PR 描述里说明物理动机/技术必要性
  2. 更新 `requirements.txt` 并钉死版本范围(如 `pydantic>=2.0,<3.0`)
  3. 更新 `docs/architecture.md` §1.4 技术栈表
- 任何 `pip install <new_package>` 之前**先询问用户**

### R11. 文档同步要求

满足以下情况之一时,**必须**同步更新 `docs/architecture.md`:

- 新增 `sqc/` 文件或模块
- 改变任何公共类/方法的签名
- 新增/修改 `sqc.config.CONFIG` 字段
- 引入新的约定(命名、单位、回归基线、时间轴规范等)
- 改变扩展规范(影响 §7 扩展指南)
- 重构改变了 R1–R10 中任意一条

**更新原则**:"只增不减"——保留所有原文档内容,新增内容追加到合适章节。绝不删除既有章节。

**版本号递增**:
- `v1.x` → `v1.(x+1)`:增加章节/扩展示例(minor)
- `v1.x` → `v2.0`:有 breaking change(如改变默认值、改变 API 签名,需大改用户使用方式)

**变更记录**:在 `docs/architecture.md` 末尾"文档维护"表追加一行,记录日期+主要变更摘要。
