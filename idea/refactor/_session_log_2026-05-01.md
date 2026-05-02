# Refactoring Session Log — 2026-05-01

> 完整记录了从 P0 到 P5 的全部重构过程，包括每个 sub-agent 的输入、输出、关键决策和验收结果。

---

## 会话上下文

- **分支**: 项目重建-v2 (基于 master fork)
- **Python**: qutip-env conda 环境 (qutip 5.2.0, numpy 2.3.1)
- **核心约束**: R1 (永不修改 src/), R3 (rtol=1e-6, atol=1e-9)
- **执行模式**: subagent-per-phase (`refactor-phase-executor`)

---

## Phase 0 — 测试基线 + 工具基础设施

### Agent 输入

```
Execute Phase 0: 建立测试基础设施。创建 tests/ 目录、pytest.ini、
conftest.py、baseline 生成脚本、物理回归测试、smoke 测试。
环境: conda python, n_levels 注意与生产代码一致。
```

### Agent 关键操作 (89 tool calls, ~99K tokens)

1. 读取 CLAUDE.md、_refactor_plan.md §0/§8/§13/§15、_handoff_state.md、phase_0_handbook.md
2. 验证 git 状态与环境 (qutip 5.2.0, numpy 2.3.1)
3. 创建 `pytest.ini` (4 个 marker: regression/integration/unit/slow)
4. 创建 `tests/` 完整目录骨架 (8 个 .py + README)
5. 创建 `tests/conftest.py` (qubit_default fixture, assert_array_close, RTOL=1e-6, ATOL=1e-9)
6. 创建 `tests/regression/generate_baselines.py` (4 个 baseline 生成函数)
7. 运行 baseline 生成: `python -m tests.regression.generate_baselines`
8. 创建 `tests/regression/test_physics_baseline.py` (4 个回归测试)
9. 创建 `tests/unit/test_smoke.py` (6 个冒烟测试)
10. 钉死 `requirements.txt` 版本: qutip>=5.0,<6.0; numpy>=1.24,<3.0; scipy>=1.11,<2.0; matplotlib>=3.7,<4.0; gradio>=4.0,<5.0; pytest>=7.4,<9.0
11. **关键决策**: n_levels=2 (因为 create_*_pulse() 默认使用 2-level 算符，用 3 会维度不匹配)

### Agent 输出

```
Status: DONE
Commits: ce3b7a9, ab444f6, 9c5c5f6, 1cd021f
A2 unit: 6 passed
A3 regression: 4 passed (24.59s)
DECISION_NEEDED: n_levels=2 vs handbook n_levels=3
DECISION_NEEDED: Qt GUI crash in Phi.plot() during headless tests
DECISION_NEEDED: gradio upgraded 3.24.1 → 4.44.1
```

### 用户验收
✅ P0 验收通过

---

## Phase 1 — sqc/ 骨架 + ABC + 数据结构 + src_mirror/ 镜像

### Agent 输入

```
Execute Phase 1: 创建 sqc/ 完整目录、所有 ABC、数据结构、
TransmonQubit 迁移、src_mirror/ 镜像层。CRITICAL: 永不修改 src/，
使用 src_mirror/ 作为 facade。P0 baselines 为 n_levels=2，
生产代码默认 n_levels=3。
```

### Agent 关键操作 (98 tool calls, ~142K tokens)

1. 启动序列: 读 CLAUDE.md → _refactor_plan.md §0/§8/§13/§15 → _handoff_state.md → phase_1_handbook.md
2. 创建 sqc/ 目录骨架 (9 个子包, 34 个 .py 文件)
3. 定义所有 ABC: Device, DistortionModel, ReadoutModel, PulseBase, RunnerBase, Experiment, Calibration, Reconstruction, Workflow
4. 定义数据结构: QubitSpec(@dataclass frozen=True), Waveform, CompositeWaveform, FluxSignal, ControlLine, TransferMatrix, PulseSequence, ExperimentResult, MeasurementTrace, CalibrationTable
5. **核心迁移**: TransmonQubit 从 src/qubit.py 逐字 port 到 sqc/devices/transmon.py
   - 内部包装 QubitSpec，保留所有 legacy 方法签名
   - `qubit_in_mag()` 委托给 HamiltonianBuilder.build()
   - `optimal_work_point()` 旧版返回弧度，新版返回 Φ₀（legacy 方法乘 π 还原）
6. HamiltonianBuilder: 纯函数，静态方法，不修改任何输入
   - 验证: H_list[0], H_list[1][0], H_list[1][1], freq_coeffs 与旧 qubit_in_mag 全等 (np.allclose)
7. 迁移 Cavity→Resonator, Coupled_System→CoupledSystem
8. 迁移 Signal→FluxSignal (type 0-8 支持), Pulse, 6 个序列工厂函数, 4 个门函数
9. **src_mirror/ 镜像层** (替代 handbook §3.10):
   - src_mirror/qubit.py → re-export from sqc.devices.transmon
   - src_mirror/signal.py → re-export from sqc.control.flux_signal
   - src_mirror/pulse.py → re-export from sqc.control.pulse + sqc.control.sequence
   - **src/ 一行未改** (R1 严格遵守)
10. 5 个新单元测试文件: test_qubit_spec(9), test_waveform(8), test_flux_signal(11), test_hamiltonian_builder(5), test_imports(10)

### Agent 输出

```
Status: DONE
Commits: 0efc938, 57fb362
A2 unit: 53 passed
A3 regression: 4 passed (36s)
镜像等价性: src_mirror.TransmonQubit is sqc.devices.transmon.TransmonQubit ✓
A7 physics sanity: f01=4.699 GHz, EJ/EC=75.0 ✓
```

### 用户验收
✅ P1 验收通过

---

## Phase 2 — 实验对象化 + KernelEstimator + Facade

### Agent 输入

```
Execute Phase 2: 实现 4 个 Experiment 类 (Rabi/Ramsey/DiffEcho/Transient)、
KernelEstimator (三合一去重)、MesolveRunner/SlidingMeasurementRunner、
IQReadoutModel、src_mirror/protocal.py facade。
CRITICAL: src.protocal.Protocal(type=N).evolve() 必须与
src_mirror.protocal.Protocal(type=N).evolve() 在 case 0/1/2/4 上数值等价。
```

### Agent 关键操作 (122 tool calls, ~193K tokens)

1. **实现 sqc/experiments/**:
   - RabiExperiment: port case 0 (Rabi 振荡)
   - RamseyExperiment: port case 1 (Ramsey 干涉) — 最关键，默认参数完全匹配旧版
   - DiffEchoExperiment: port case 2 (差分回波)
   - TransientSensingExperiment: port case 4 (瞬态传感滑动测量)
2. **实现 KernelEstimator**: 合并 Pulse.get_kernel / CompositePulse.get_kernel / Analysis.get_kernel 三处重复实现。stim_amplitude=0.0215, stim_width=3.0 保持兼容
3. **实现 MesolveRunner + SlidingMeasurementRunner**: 完整 port single_measurement 和 sliding_measrement
4. **实现 IQReadoutModel**: port IQ_readout() 模块级函数为类方法
5. **创建 src_mirror/protocal.py facade**:
   - Protocal.evolve() 返回与旧版完全相同的元组
   - Case 0: Result 对象, Case 1: (Phi, tau_list, p_e_list), Case 2: (Phi, tau_list, p_e_list, k, t_int), Case 4: (t_samples, kernel, scan_list, delta_p, p_e, Phi, control_pulse)
6. **Bug 修复**:
   - FluxSignal type-4 公式: `exp(-t/fall) - exp(-t/rise)` → 正确双指数
   - PulseBase ABC `@property @abstractmethod` 阻止子类实例化 → 改为普通 @abstractmethod
   - create_ramsey_pulse / create_diff_echo_pulse 缺少 qubit 参数 → 添加
7. **等价性测试**: test_protocal_mirror.py 验证 case 0/1/2/4 数值等价

### Agent 输出

```
Status: DONE
Commits: cdd32a3 (用户手动 commit)
A2 unit: 67 passed
A3 regression: 4 passed
A4 equivalence: 4 passed (case 0/1/2/4)
A6 integration: 6 passed
Total: 81 passed
```

### 用户验收
✅ P2 验收通过

---

## Phase 3a — basis + Wiener/RamseyIQ/Hammerstein 内化

### Agent 输入

```
Execute Phase 3a: basis 模块 + Wiener/RamseyIQ/RamseyUnwrap/DiffEcho
reconstruction 类 + HammersteinWiener + src_mirror/analysis.py facade。
不依赖 Track B 的部分。LM/Cryoscope 留 NotImplementedError。
```

### Agent 关键操作 (85 tool calls, ~151K tokens)

1. **sqc/reconstruction/basis.py**: 逐字 port generate_basis_functions / basis_function_decomposition / regularization_matrix (Fourier + B-spline + Legendre)
2. **sqc/reconstruction/wiener.py**: 4 个 Reconstruction 子类
   - WienerReconstruction(lambda_reg) — 线性 Wiener 反卷积
   - RamseyIQReconstruction(qubit) — IQ 解调 → arcsin → B(τ)
   - RamseyUnwrapReconstruction(qubit, k_span=3) — 相位解缠绕
   - DiffEchoReconstruction(qubit, t_int, k) — 差分回波直接公式
   - 双重 sensitivity 检测: QubitSpec (新) 和 legacy TransmonQubit
3. **sqc/reconstruction/hammerstein.py**: HammersteinWienerReconstruction (块结构非线性反演，内联 Wiener + 反色散)
4. **CalibrationTable 扩展**: 添加 inputs/outputs/kind/fit_params 字段 + evaluate()/inverse() 插值方法 (cubic → quadratic → linear 自动回退)
5. **src_mirror/analysis.py facade** (新文件，R1 合规):
   - 已迁移: get_expectation_values, get_population, get_signal_from_ramsey_by_iq, get_signal_from_ramsey_by_unwrap, get_signal_from_diff_echo, get_kernel, wiener_deconvolution, hammerstein_wiener_deconvolution
   - 延迟 (NotImplementedError): numerical_inverse, get_signal_from_cryoscope, get_h_from_phi
6. 31 个新测试 (14 basis + 9 calibration_table + 8 analysis mirror equivalence)

### Agent 输出

```
Status: DONE
Commits: 5488e20, 844e875, 2eafda9, 6bed1af
A2 unit: 90 passed (+23)
A4 equivalence: 12 passed (+8 analysis mirror)
Total: 112 passed
```

### 用户验收
✅ P3a 验收通过

---

## Phase 3b — LMReconstruction 完整内化

### Agent 输入

```
Execute Phase 3b: Track B 0.3 (LM 收敛修复) 已完成。
Port src/analysis.py:272-825 → sqc/reconstruction/numerical_inverse.py。
LMReconstruction + 伴随 Jacobian + 有限差分 Jacobian + LM 优化循环。
更新 src_mirror/analysis.py numerical_inverse → LMReconstruction。
生成 lm_default.pkl baseline。
```

### Agent 关键操作 (85 tool calls, ~134K tokens)

1. **sqc/reconstruction/numerical_inverse.py** (~550 行):
   - LMReconstruction dataclass with reconstruct() public API
   - _forward_simulation() — 每个测量延迟运行 mesolve
   - _compute_jacobian_adjoint() — 伴随状态法 (solve_ivp)
   - _compute_jacobian_fd() — 有限差分 fallback
   - _levenberg_marquardt() — 主 LM 优化循环
   - _build_h_for_signal() — 每个延迟的 Hamiltonian list 构造器
2. **src_mirror/analysis.py**: numerical_inverse() 不再是 NotImplementedError → 委托 LMReconstruction
3. 模块级 facade: forward_simulation, compute_jacobian, compute_jacobian_fd, levenberg_marquardt
4. **lm_default.pkl baseline**: 最小 LM baseline (n_basis=5, max_iter=2, 11 time points)
5. 5 个 LM 单元测试 + 1 个 LM 回归测试

### Agent 输出

```
Status: DONE
Commits: 4cd6bc7, 3f37ede, a8ec1f2
A2 unit: 95 passed (+5 LM tests)
A3 regression: 5 passed (+1 LM baseline)
Total: 118 passed
Known caveat: LM baseline signal weak (amplitude=0.01, p_e ~0.5)
```

### 用户验收
✅ P3b 验收通过

---

## Phase 3c (PARTIAL) — Cryoscope/Calibration 内化

### Agent 输入

```
Execute P3c PARTIAL: Track B 1.1/1.2 未完成。
完成无依赖部分: CryoscopeExperiment (case 5), CryoscopeReconstruction 骨架,
FluxResponseCalibration (ramsey only), QubitFrequencyCalibration,
Calibration facade。Track B 依赖部分留 NotImplementedError stub。
务必在文档和交接文档中注明特殊情况。
```

### Agent 关键操作 (92 tool calls, ~137K tokens)

1. **sqc/experiments/cryoscope.py**: 完整 port case 5 (truncation scan + IQ readout → varphi)
2. **sqc/reconstruction/cryoscope.py**: 类骨架 (reconstruct 方法 raise NotImplementedError("Track B 1.1"))
3. **sqc/calibration/flux_response.py**: method="ramsey" 实现 (扫描 DC flux steps → fit f(Φ)), method="cryoscope"/"transient" stub
4. **sqc/calibration/qubit_frequency.py**: QubitFrequencyCalibration (case 0, FFT-based frequency fitting) + TransientFrequencyCalibration stub
5. **src_mirror/protocal.py**: Calibration facade (type 0 → QubitFrequencyCalibration, type 1 → FluxResponseCalibration(ramsey), type 2/3 → NotImplementedError)
6. **src_mirror/analysis.py**: get_h_from_phi 完整实现 (纯插值), get_signal_from_cryoscope deferred
7. **Bug 修复**: IQReadoutModel 现在将 qubit 传给 create_ramsey_pulse (修复 n_levels 不匹配)
8. 所有 NotImplementedError 消息包含 "Track B X.Y" 引用

### Agent 输出

```
Status: PARTIAL
Commit: aab4045
A2 unit: 103 passed
A4 equivalence: 14 passed (+2: case 5 + get_h_from_phi)
Total: 128 passed
Stubs: 6 处 NotImplementedError (均标注 Track B 引用)
```

### 用户验收
✅ P3c 验收通过（接受 PARTIAL 状态）

---

## Track B 1.3 in src_mirror/ — DistortionModel 实现

### 主 session 直接实现

在 `src_mirror/distortion.py` 中创建完整的 DistortionModel 实现（R1 合规，不修改 src/）：

1. **DistortionModel ABC**: apply(waveform, dt), step_response(t), impulse_response(t), frequency_response(omega), apply_to_waveform(wf)
2. **SingleExponentialDistortion**: h(t) = (1-amp)δ(t) + (amp/τ)exp(-t/τ)Θ(t)
   - IIR 离散化 (bilinear transform)
   - 解析频率响应: H(ω) = (1-amp) + amp/(1+jωτ)
3. **MultiExponentialDistortion**: K 个指数尾巴的并行和
   - 每路独立 bilinear + 直接路径求和
4. **FIRDistortion**: 通用 FIR 滤波器 (lfilter)
5. **IIRDistortion**: 通用 IIR 滤波器 (b/a coefficients)
6. **CustomTransferDistortion**: 自定义复频域 H(ω) (FFT 乘法)

### Smoke Test 结果

```
Test 1 - Step settles to: 0.999996 ✓
Test 2 - Step response max diff: 4.99e-02 (bilinear approximation)
Test 3 - MultiExp settles to: 0.997 ✓ (after bug fix)
Test 4 - FIR impulse sum*dt: 1.0000 ✓
Test 5 - IIR step settles to: 1.000000 ✓
Test 6 - Custom(identity) max err: 4.44e-16 ✓
Test 7 - Waveform apply: OK ✓
Test 8 - Freq response max err: 0.00e+00 ✓
```

### Bug 修复
- MultiExponentialDistortion.apply(): 原始 SOS cascade 实现错误，改为正确的并行 bilinear 求和

### Commit: 07fdc9b

---

## Phase 4 — ControlLine + DistortionModel + Predistortion

### Agent 输入

```
Execute Phase 4: Track B 1.3 完成 (src_mirror/distortion.py)。
Internalize DistortionModel 到 sqc/hardware/、实现 ControlLine、
TransferFunctionCalibration、PredistortionDesigner、
PredistortionValidationWorkflow。
```

### Agent 关键操作 (102 tool calls, ~130K tokens)

1. **sqc/hardware/distortion.py**: 完整实现 (替换 P1 stub) — 5 个子类 + apply_to_waveform()
2. **sqc/hardware/control_line.py**: 完整实现 — apply() + predistort() + 向后兼容物理字段
3. **sqc/calibration/transfer_function.py**: TransferFunctionCalibration (阶跃响应测量 + scipy curve_fit 多模型拟合 + to_distortion_model())
4. **sqc/calibration/predistortion.py**: PredistortionDesigner
   - 自动检测: 指数模型 → 解析 IIR 逆 (bilinear)
   - 通用模型 → 频域逆: H_inv = conj(H) / (|H|² + λ²)
   - analytical_inverse 方法: SingleExp → 完美抵消
5. **sqc/workflows/predistortion_validation.py**: 端到端验证工作流
6. **src_mirror/distortion.py**: 改为 re-export from sqc.hardware.distortion (facade)
7. 45 个新测试 (24 distortion + 9 control_line + 12 designer), 9 个集成测试, 1 个回归测试

### 预失真性能

| 模型 | 改善倍数 |
|---|---|
| SingleExponentialDistortion | ~8.5×10¹² (解析逆完美抵消) |
| MultiExponentialDistortion | ~1× (频域逆，bilinear warping 残留) |

### Agent 输出

```
Status: DONE
Commits: fdc55c8, d21194b, 235831f
A2 unit: 148 passed (+45)
A3 regression: 6 passed (+1 predistortion)
Total: 183 passed
```

### 用户验收
✅ P4 验收通过

---

## Phase 5 — TransferMatrix + ChipTopology + Z-Crosstalk

### Agent 输入

```
Execute Phase 5 (FINAL): TransferMatrix 完整实现、ChipTopology 扩展、
ZCrosstalkWorkflow、双 qubit Z-crosstalk demo。
```

### Agent 关键操作 (87 tool calls, ~163K tokens)

1. **sqc/hardware/transfer_matrix.py**: 完整 FFT 实现
   - apply(): source voltages → FFT → multiply H_ji(ω) → sum → IFFT → on-chip fluxes
   - H_ji(), diagonal(), off_diagonal() 查询
   - from_dc_matrix() 工厂 (频率平坦)
   - 自动 resampling 处理不匹配的频率网格
   - 零填充处理不同长度波形
2. **sqc/devices/chip.py**: ChipTopology 类
   - lift_qubit_op() — 将单 qubit 算符嵌入完整张量积 Hilbert 空间
   - hamiltonian_static() — 所有器件哈密顿量之和
   - collapse_operators() — 所有 Lindblad 算符
   - from_legacy_coupled_system() — 适配现有 CoupledSystem
3. **sqc/workflows/z_crosstalk.py**: ZCrosstalkWorkflow
   - Step 1: 应用真实 TransferMatrix → on-chip fluxes
   - Step 2: 在 qubit B 上运行 TransientSensingExperiment 测量寄生磁通
   - Step 3: Wiener 反卷积重建 Φ_B(t)
   - Step 4: 正则化频域解卷积提取 H_BA(ω): H = (Φ_B · conj(V_A)) / (|V_A|² + λ²)
   - Step 5: 设计与验证补偿
4. 37 个新单元测试 + 7 个集成测试 + 1 个回归测试

### Agent 输出

```
Status: DONE
Commits: 1877732, dac77bc
A2 unit: 185 passed (+37)
A3 regression: 7 passed (+1 z_crosstalk)
A4 equivalence: 14 passed
A6 integration: 24 passed (+9)
Total: 230 passed
```

### 用户验收
✅ P5 验收通过 — 全方案重构完成

---

## 最终交付 — 文档 + 前端 Demo

### 主 session 直接实现

1. **docs/architecture.md** (~500 行, 10 章节):
   - 平台概览、物理背景、全栈架构、模块详解 (每个子包)
   - 数据流、快速开始、扩展指南 (新实验/算法/模型)
   - API 参考 (核心类索引 + 关键函数签名)
   - 测试与回归机制、设计原则与约定

2. **web_demo_v2.py** (6 个 Tab 的交互式 Gradio 应用):
   - Tab 0: 🏗️ Architecture — cQED 六层可视化架构图
   - Tab 1: ⚛️ Qubit Config — 参数滑块 + 实时 sanity check
   - Tab 2: 📈 Waveform Designer — 7 种信号类型交互设计
   - Tab 3: 🧪 Protocol Runner — 5 种协议一键执行
   - Tab 4: 🔄 Reconstruction — Wiener + Kernel
   - Tab 5: 🔧 Distortion — 预失真闭环验证 (4 子图)
   - Tab 6: 📖 About — 版本信息

### Commit: d46c25e

---

## 重构完整统计

### Commits (共 20 个)

```
d46c25e Final: technical documentation + visual frontend demo
dac77bc P5: update handoff state
1877732 P5: TransferMatrix + ChipTopology + ZCrosstalkWorkflow
235831f P4: update handoff state
d21194b P4: calibration + predistortion + validation workflow + tests
fdc55c8 P4: internalize DistortionModel + ControlLine
07fdc9b Track B 1.3: DistortionModel implementations in src_mirror/
aab4045 P3c (PARTIAL): CryoscopeExperiment + Calibration facade + stubs
a8ec1f2 P3b: update handoff state
3f37ede P3b: add LM baseline + regression test
4cd6bc7 P3b: LMReconstruction internalization
6bed1af P3a: update handoff state
2eafda9 P3a: src_mirror/analysis.py facade
844e875 P3a: Wiener + Hammerstein reconstruction classes
5488e20 P3a: basis module
cdd32a3 P2: Experiment classes + KernelEstimator + runners + facade
57fb362 P1: update handoff state
0efc938 P1: sqc/ skeleton + ABCs + data structures + src_mirror/
1cd021f P0: update handoff state
9c5c5f6 Restore src/__pycache__ tracking (R1 compliance)
ab444f6 Add .gitignore, remove __pycache__ from tracking
ce3b7a9 P0: test infrastructure + baselines + smoke tests
```

### 文件统计

| 类别 | 数量 |
|---|---|
| sqc/ 新增文件 | ~50 .py |
| src_mirror/ 新增文件 | 6 |
| tests/ 新增文件 | ~25 |
| docs/ 新增文件 | 1 |
| web_demo_v2.py | 1 |
| 修改文件 | ~10 |
| src/ 修改 | **0** |

### 测试统计

| 套件 | 最终数量 |
|---|---|
| 单元测试 (unit) | 185 |
| 回归测试 (regression) | 7 |
| 等价性测试 (equivalence) | 14 |
| 集成测试 (integration) | 24 |
| **总计** | **230** |

### Sub-agent 消耗

| Phase | Agent | Tool Calls | Tokens | Duration |
|---|---|---|---|---|
| P0 | a2cde3aaf5b873688 | 89 | ~99K | 33 min |
| P1 | aaf408813cbe44b33 | 98 | ~142K | 56 min |
| P2 | acddd8109621c321e | 122 | ~193K | 105 min |
| P3a | a2a4a80513a4d1377 | 85 | ~151K | 53 min |
| P3b | a956069ec5079aa82 | 85 | ~134K | 46 min |
| P3c | ae2e5a2d550face6c | 92 | ~137K | 159 min |
| P4 | a23d281f6f2afe1be | 102 | ~130K | 152 min |
| P5 | a5f26a15306c24ead | 87 | ~163K | 321 min |
| **Total** | — | **760** | **~1.15M** | **~15.5 hr** |

---

## 已知问题与待办

### DECISION_NEEDED (来自各 phase)

1. **n_levels=2 vs 3**: P0 baselines 锚定在 n_levels=2, 生产代码默认 n_levels=3
2. **Qt GUI crash**: Phi.plot() 在 headless 环境触发 Windows 异常
3. **LM baseline 信号弱**: 需要更强的测试信号
4. **MultiExp predistortion 性能**: 频域逆与 bilinear filter 不完美匹配 (~1x improvement)

### Track B 待办 (阻塞 P3c 完整)

| 任务 | 状态 |
|---|---|
| 0.1 case 1 死代码清理 | ○ |
| 0.2 kernel 自动校准 | ○ |
| 0.3 LM 收敛修复 | ✓ (P3b 前完成) |
| 1.1 Cryoscope case 6/7 | ○ — 阻塞 P3c-full |
| 1.2 瞬态频率标定 case 8 | ○ — 阻塞 P3c-full |
| 1.3 DistortionModel | ✓ (P4 前完成) |

### P3c Stub 清单

| 位置 | 阻塞条件 |
|---|---|
| CryoscopeReconstruction.reconstruct() | Track B 1.1 |
| FluxResponseCalibration._calibrate_cryoscope() | Track B 1.1 |
| FluxResponseCalibration._calibrate_transient() | Track B 1.2 |
| TransientFrequencyCalibration.calibrate() | Track B 1.2 |
| Calibration type 2, 3 | Track B 1.1/1.2 |
| get_signal_from_cryoscope() | Track B 1.1 |

---

## 会话文件变更总览

```
NEW:     sqc/                        (~50 .py)
NEW:     src_mirror/                 (6 .py)
NEW:     tests/                      (~25 .py + 7 .pkl)
NEW:     docs/architecture.md
NEW:     web_demo_v2.py
NEW:     pytest.ini
MODIFIED: requirements.txt
MODIFIED: idea/refactor/_handoff_state.md
MODIFIED: .gitignore
UNCHANGED: src/                      (ALL files — R1 compliance)
UNCHANGED: web_demo.py
UNCHANGED: Simulation.ipynb
```

---

*日志结束。下一次 session 可从 `idea/refactor/_handoff_state.md` 读取当前进度。*
