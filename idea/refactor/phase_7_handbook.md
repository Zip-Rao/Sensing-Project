# Phase 7 — 统一 mesolve 时间轴到 t_global

> **状态**: 规划阶段，待实施。解决 sqc/ 内时间轴"全局画布 + 局部脉冲"混用导致的反演伪影与 API 不一致。
> **前置**: P0–P6 全部完成（推荐 P6a 先做，因 CONFIG 参数接线影响本方案的迁移模板）。
> **配套文件**:
> - 本方案详细 plan: [unified-kindling-metcalfe.md](../../../.claude/plans/unified-kindling-metcalfe.md)
> - 全局配置规范: [CLAUDE.md §R9](../../CLAUDE.md)
> - 重构总方案: [_refactor_plan.md](_refactor_plan.md) §8 (兼容层契约)、§13 (子代理执行模式)
> - 当前交接状态: [_handoff_state.md](_handoff_state.md)

---

## §1 问题诊断: 时间轴混用的三个致命问题

### 1.1 不同实验 mesolve 跑在不同轴上

| 实验 | mesolve 的 `t_list` | 来源 |
|---|---|---|
| RamseyExperiment / EchoExperiment / Calibration | `CONFIG.pulse.t_global` (-50→400ns, 900pt) | [sqc/experiments/ramsey.py:126](sqc/experiments/ramsey.py) |
| CryoscopeExperiment (via IQReadoutModel) | `ctrl_I.t_list` (0→2·t_rabi+tau, ~140pt) | [sqc/hardware/readout.py:134,148-149](sqc/hardware/readout.py) |
| TransientSensingExperiment (via SlidingMeasurementRunner) | `np.linspace(t_start, t_end, N)` | [sqc/simulation/runner.py:317](sqc/simulation/runner.py) |

**后果**: 改 `CONFIG.pulse.t_global` 不影响 Cryoscope；改 IQReadoutModel 的 tau 不影响 Ramsey。用户在 notebook 里调参数时行为不一致。

### 1.2 没有统一的"把对象放到 global"机制

- 每个实验调用方手动 `ctrl.t_list -= self.t_rabi[-1]` ([ramsey.py:111](sqc/experiments/ramsey.py)) 做偏移
- `CompositePulse.get_t_list()` 用 `curr += duration + 1e-9` 拼时间轴 ([pulse.py:354-368](sqc/control/pulse.py))—— `1e-9` 与 `dt=0.5ns` 的网格哲学冲突
- SlidingMeasurementRunner 用 `np.searchsorted` + `np.linspace` 重采样——违反 R9
- `HamiltonianBuilder.build` 在 `flux_signal=None` 时 fallback 到 `np.linspace(0, 100, 100)` ([hamiltonian.py:78](sqc/simulation/hamiltonian.py))

**后果**: 每个新实验都要手写对齐代码，没有可复用的 helper，容易出错且无单元测试覆盖。

### 1.3 Cryoscope 反演的静默截断 bug（用户报告）

用户用 Cryoscope 重建 wave-packet 正弦信号：把 `trunc_list` 延长到 100ns 时，t > 60ns 的重建曲线 flatline 到 0，完全丢失 t≈70ns 的深谷。

**根因链**：
1. Cryoscope 构造 `FluxSignal` 用 `t_list=make_time(0, 80)` ([cryoscope.py:66](sqc/experiments/cryoscope.py))
2. 用户只延长了 `trunc_list` 但没延长 `flux_signal.t_list`
3. `truncate()` 是零掩码 (mask `t > t_d`)，当 `t_d > t_list[-1]=80` 时 mask 全 False → truncate 是 no-op → 所有超长 trunc 点测得同一个总相位常数 → dφ/dt = 0 → flatline
4. 框架没有 `trunc_list.max() > flux_signal.t_list[-1]` 的边界检查，**silent failure**

---

## §2 目标状态: 统一的全局时间轴模式

### 2.0 核心语义: `trigger` + 自带时间轴

每个 `Pulse` / `FluxSignal` 持有两个东西：

- **`trigger: float`** — 该对象在 global 时间轴上的**触发时刻**（即 "local t=0 对应的 global t"）。物理直觉对标 AWG 触发点。
- **`t_list: np.ndarray`** — 该对象**自己的局部时间轴**，永远从 0 起（典型 `arange(0, duration, dt)`）。

**拼接规则**：把对象映射到 global 时，对每个 global 时刻 `t_g`，先算 local 时间 `t_loc = t_g - trigger`，若 `t_loc ∈ [t_list[0], t_list[-1]]` 则取该对象 samples/coeff 的插值，否则置 0。

举例：一个 10 ns 的 π/2 脉冲，`trigger=30, t_list=[0, 0.5, ..., 9.5]`，在 global 轴上**只在 t∈[30, 39.5] 内有非零贡献**，其余处为 0。Ramsey 的两个 π/2 各自带 trigger（如 `0` 和 `t_rabi[-1]+tau`）。

### 2.1 统一模板：所有实验用同一段代码

```python
t_global = self.t_global  # 来自 CONFIG.pulse.t_global

# 1) flux 投影到 t_global（FluxSignal 自带 trigger）
flux_samples_global = self.flux_signal.samples_on(t_global)
flux_global = FluxSignal(type=8, t_list=t_global, signal=flux_samples_global, trigger=0.0)
self.qubit.qubit_in_mag(flux_global, frame=1, omega_d=self.omega_d)

# 2) 控制脉冲投影到 t_global（各子脉冲自带 trigger）
ctrl = create_xxx_pulse(..., trigger=...)
H_ctrl = ctrl.hamiltonian_on(t_global)

# 3) 合并 H，统一在 t_global 上 mesolve
H = (QobjEvo(self.qubit.H_list, tlist=t_global, order=1)
     + QobjEvo(H_ctrl, tlist=t_global, order=1))
result = mesolve(H, self.qubit.state, t_global, [], e_ops=[...])
```

**关键点**：`qubit.H_list` 的 `tlist` 也是 `t_global`（因为传入 `flux_global` 的 t_list 就是 t_global），两个 QobjEvo 项的 tlist 一致，QuTiP 不再跨界外推。

### 2.2 干掉的东西

| 旧代码 | 替换 |
|---|---|
| `ctrl.t_list -= self.t_rabi[-1]` 手动偏移 | `create_ramsey_pulse(..., trigger=t_rabi[-1])` — 第一个 π/2 的 trigger 声明第二个 π/2 的等待 |
| `CompositePulse.get_t_list()` 累加 + `1e-9` | `CompositePulse.hamiltonian_on(t_global)` 收集各子脉冲投影 |
| `IQReadoutModel` 用 `ctrl_I.t_list` | 投影到 `t_global` 后用 `t_global` 跑 mesolve |
| `np.linspace(...)` 生成时间轴 | `CONFIG.pulse.make_time(start, end)` 或 `np.arange(start, end, dt)` |
| `HamiltonianBuilder` fallback `np.linspace(0,100,100)` | `CONFIG.pulse.make_time(0, 100)` |

---

## §3 新增 API

### 3.1 `FluxSignal` 扩展

```python
class FluxSignal(Waveform):
    def __init__(self, ..., trigger: float = 0.0):
        self.trigger = trigger
        # 其余不变；samples 仍然定义在 local t_list 上（永远从 0 起）

    def samples_on(self, t_global: np.ndarray) -> np.ndarray:
        """投影到 global：每个 global 时刻 t_g 取 local 时间 t_loc=t_g-trigger 的插值。
        [trigger, trigger + t_list[-1]] 之外置 0。
        """
        out = np.zeros(len(t_global), dtype=float)
        t_loc = t_global - self.trigger
        mask = (t_loc >= self.t_list[0]) & (t_loc <= self.t_list[-1])
        out[mask] = np.interp(t_loc[mask], self.t_list, self.samples)
        return out
```

### 3.2 `Pulse` 扩展

```python
class Pulse:
    def __init__(self, ..., trigger: float = 0.0):
        self.trigger = trigger
        # 其余不变；hamiltonian/coeffs 仍然定义在 local t_list 上

    def hamiltonian_on(self, t_global: np.ndarray) -> list:
        """返回 H_list，coeffs 零填充到 len(t_global)。"""
        H_local = self.hamiltonian
        result = []
        t_loc = t_global - self.trigger
        mask = (t_loc >= self.t_list[0]) & (t_loc <= self.t_list[-1])
        for op, coeff_local in H_local:
            if np.iscomplexobj(coeff_local):
                coeff_global = np.zeros(len(t_global), dtype=complex)
                coeff_global[mask] = (
                    np.interp(t_loc[mask], self.t_list, coeff_local.real)
                    + 1j * np.interp(t_loc[mask], self.t_list, coeff_local.imag)
                )
            else:
                coeff_global = np.zeros(len(t_global), dtype=float)
                coeff_global[mask] = np.interp(t_loc[mask], self.t_list, coeff_local)
            result.append([op, coeff_global])
        return result
```

### 3.3 `CompositePulse` 简化

设计哲学：**CompositePulse 不再"拼接时间轴"，只是子脉冲的集合**。每个子脉冲自带 trigger。

```python
class CompositePulse:
    def hamiltonian_on(self, t_global: np.ndarray) -> list:
        """收集所有子脉冲的全局贡献。各子脉冲 trigger 是绝对值（相对 global t=0）。"""
        merged = []
        for pulse in self.pulses:
            merged.extend(pulse.hamiltonian_on(t_global))
        return merged
```

- 删除 `get_t_list()` 中的 `curr += duration + 1e-9` 脏 trick
- 旧的 `get_t_list()` / `get_hamiltonian()` 保留为 deprecated

### 3.4 Factory 函数签名变化

所有 `create_*_pulse` 函数加 `trigger=0.0` 参数：

```python
def create_ramsey_pulse(t_rabi, tau, omega_d, phase1, phase2, qubit, trigger=0.0):
    pulse1 = Pulse(..., phase=phase1, trigger=trigger)
    pulse2 = Pulse(..., phase=phase2, trigger=trigger + t_rabi[-1] + tau)
    return CompositePulse([pulse1, pulse2])
```

调用方再也不需要 `ctrl.t_list -= self.t_rabi[-1]`。

---

## §4 受影响文件（基于完整审计）

### Tier B: 底层 API 扩展（3 文件）

| 文件 | 改动 |
|---|---|
| [sqc/control/pulse.py](sqc/control/pulse.py) | 加 `trigger` + `hamiltonian_on`；`CompositePulse.hamiltonian_on` 调用子脉冲投影；保留旧方法 deprecated |
| [sqc/control/flux_signal.py](sqc/control/flux_signal.py) | 加 `trigger` + `samples_on` |
| [sqc/control/sequence.py](sqc/control/sequence.py) | 所有 `create_*_pulse` factory 加 `trigger=0.0` 参数 |

### Tier A1: 实验层迁移（7 文件）

| 文件 | 改动 |
|---|---|
| [sqc/experiments/rabi.py](sqc/experiments/rabi.py) | mesolve 改用 t_global；rabi pulse 设 trigger=0 |
| [sqc/experiments/ramsey.py](sqc/experiments/ramsey.py) | 删 `ctrl.t_list -= self.t_rabi[-1]` hack；用新模板 |
| [sqc/experiments/echo.py](sqc/experiments/echo.py) | 用新模板，确认无 offset hack |
| [sqc/experiments/cryoscope.py](sqc/experiments/cryoscope.py) | flux_signal 默认 t_list 延长到 `make_time(0, 100)`；用新模板；**加 trunc_list 边界 sanity check** |
| [sqc/experiments/delay_ramsey.py](sqc/experiments/delay_ramsey.py) | 删手搓 `t_sig` + windowed signal 逻辑；flux 用 `trigger=t_d` 直接定位 |
| [sqc/experiments/pi_pulse_comp.py](sqc/experiments/pi_pulse_comp.py) | mesolve 改用 t_global |
| [sqc/experiments/transient.py](sqc/experiments/transient.py) | 通过 `SlidingMeasurementRunner` 间接迁移 |

### Tier A2–A6: 硬件/标定/重建/仿真/Workflow 层（9 文件）

| 文件 | 改动 |
|---|---|
| [sqc/hardware/readout.py](sqc/hardware/readout.py) | `IQReadoutModel.measure` 改用 t_global 模板；删 `t_evolve = ctrl_I.t_list` |
| [sqc/calibration/frequency.py](sqc/calibration/frequency.py) | 4 处 mesolve 全部迁移到 t_global 模板 |
| [sqc/reconstruction/transient.py](sqc/reconstruction/transient.py) | mesolve 改用 t_global；同步 forward sim 时间轴逻辑 |
| [sqc/reconstruction/ramsey.py](sqc/reconstruction/ramsey.py) | mesolve 改用 t_global |
| [sqc/reconstruction/kernel.py](sqc/reconstruction/kernel.py) | kernel 估计的 mesolve 改用 t_global；stimulus 用 trigger 定位 |
| [sqc/reconstruction/delay_ramsey.py](sqc/reconstruction/delay_ramsey.py) | 与 `DelayRamseyExperiment` 配套迁移 |
| [sqc/simulation/runner.py](sqc/simulation/runner.py) | `SlidingMeasurementRunner` 去 `np.linspace` (line 158,239,317)；改用 t_global |
| [sqc/simulation/hamiltonian.py](sqc/simulation/hamiltonian.py) | line 78 fallback `np.linspace(0,100,100)` → `CONFIG.pulse.make_time(0, 100)` |
| [sqc/workflows/z_crosstalk.py](sqc/workflows/z_crosstalk.py) | 内部组合实验配合新 API |
| [sqc/workflows/predistortion_validation.py](sqc/workflows/predistortion_validation.py) | 内部组合实验配合新 API |

### Tier C: 审计但不迁移（6 文件，非时间轴用途）

[sqc/devices/{transmon,chip,base,coupler,resonator}.py](sqc/devices/) — 不直接调 mesolve；[sqc/hardware/{distortion,electronics,transfer_matrix}.py](sqc/hardware/) — `np.linspace` 用于参数扫描非时间轴，不违反 R9。

### 测试影响

**直接需修复的单元测试**（API 形态变化或 t_list 长度断言）：
- [tests/unit/test_cryoscope_experiment.py](tests/unit/test_cryoscope_experiment.py) — `phi.t_list[60:10:-1]` 切片范围
- [tests/unit/test_delay_ramsey_experiment.py](tests/unit/test_delay_ramsey_experiment.py) — 长度断言
- [tests/unit/test_flux_signal.py](tests/unit/test_flux_signal.py) — CompositeSignal 合并测试
- [tests/unit/test_waveform.py](tests/unit/test_waveform.py) — CompositeWaveform 长度/起止断言
- [tests/unit/test_hamiltonian_builder.py](tests/unit/test_hamiltonian_builder.py) — 数值变化
- [tests/unit/test_lm_reconstruction.py](tests/unit/test_lm_reconstruction.py) — `meas_start` 计算
- [tests/unit/test_control_line.py](tests/unit/test_control_line.py) — 形状断言
- [tests/unit/test_transfer_matrix.py](tests/unit/test_transfer_matrix.py) — 固定长度断言
- [tests/unit/test_runners.py](tests/unit/test_runners.py) — 等价性测试

**重新生成的 baseline**（mesolve 数值变化）：全部 7 个 pkl。
**对应回归测试**：`test_physics_baseline.py` (5 子)、`test_predistortion_baseline.py`、`test_z_crosstalk_baseline.py`、`test_analysis_mirror.py`。

### Notebook 影响

[Simulation_sqc.ipynb](../../Simulation_sqc.ipynb) 14+ cells 触碰时间轴。用户调用代码**不需要改**（API 向后兼容），但需重跑 cells 视觉确认。重点验证 cell 15 (Cryoscope 主流程) 和 cell 29-31 (阶跃响应)。

---

## §5 执行计划（10 个 Phase，每个末尾 commit checkpoint）

### P7.0: 安全检查点 + 改前基线

```bash
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m pytest tests/ -v > /tmp/pre_refactor_results.txt
git add -A
git commit -m "safety checkpoint: before unified global time axis refactor (P7)"
```
记录此 SHA 作为回滚锚点。

### P7.1: Tier B — 底层 API 扩展（不破坏旧调用）

文件: `pulse.py`, `flux_signal.py`, `sequence.py`

- 给 `FluxSignal` 加 `trigger` + `samples_on`
- 给 `Pulse` 加 `trigger` + `hamiltonian_on`
- `CompositePulse.hamiltonian_on` = extend 子脉冲投影
- **不动** `get_t_list()` 的 `1e-9` trick（留待 P7.5 清）
- 工厂函数加 `trigger=0.0` 形参

新写测试: `tests/unit/test_global_embedding.py`。
验收: `pytest tests/unit/test_global_embedding.py tests/unit/test_flux_signal.py tests/unit/test_waveform.py -v` 全绿。
Commit: `"feat(P7): add trigger + global embedding API to Pulse/FluxSignal"`

### P7.2: Tier A2 + A5 — 底层 mesolve 消费者

文件: `readout.py`, `runner.py`, `hamiltonian.py`

- `IQReadoutModel.measure` 改用 t_global 模板
- `HamiltonianBuilder.build` line 78 fallback 修复

验收: `pytest tests/unit/test_runners.py tests/unit/test_hamiltonian_builder.py -v`。
Commit: `"refactor(P7): IQReadoutModel uses t_global; fix HamiltonianBuilder linspace"`

### P7.3: Tier A1 — 实验层逐个迁移

按依赖顺序：Rabi → Ramsey → Echo → Cryoscope → DelayRamsey → PiPulseComp → Transient (+ SlidingMeasurementRunner)

每个子步骤跑对应单测，全绿再进下一个。Cryoscope 顺带修复默认 `flux_signal.t_list` 延长到 100ns + 加 trunc 边界 sanity check。

Commit: 每个实验一次：`"refactor(P7): migrate {Name} to t_global"`

### P7.4: Tier A4 — 重建层迁移

文件: `reconstruction/{transient,ramsey,kernel,delay_ramsey}.py`

验收: `pytest tests/unit/test_lm_reconstruction.py -v`。
Commit: `"refactor(P7): migrate reconstruction layer to t_global"`

### P7.5: R9 清理

- `SlidingMeasurementRunner` 的 `np.linspace` → `arange`（如 P7.3 未清完）
- `CompositePulse.get_t_list()` 干掉 `1e-9` 分隔
- `hamiltonian.py` line 78 最终确认

Commit: `"refactor(P7): clean R9 violations — remove np.linspace from time axes"`

### P7.6: Tier A6 — 标定 + Workflow

文件: `calibration/frequency.py`, `workflows/{z_crosstalk,predistortion_validation}.py`

验收: `pytest tests/unit -v` 全套。
Commit: `"refactor(P7): migrate calibration + workflow layer to t_global"`

### P7.7: 测试修复

集中修 P7.1–P7.6 期间发现的单元测试失败（长度断言、切片范围等）。借此机会把 `len(t_list) == N` 类脆断言改为 `>=` 或 `% dt == 0` 等更鲁棒的形态断言。

Commit: `"test(P7): fix unit tests affected by t_global unification"`

### P7.8: 重新生成 baseline + 回归测试

```bash
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m tests.regression.generate_baselines
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m pytest tests/regression tests/equivalence -v
```

7 个 baseline pkl 全部刷新。Equivalence 测试必要时把容差从 1e-6 放宽到 1e-4。

Commit: `"test(P7): regenerate baselines for unified t_global

intentional physics change per CLAUDE.md R9: mesolve time axis unified to
CONFIG.pulse.t_global for all experiments. All 7 regression baselines
regenerated; values shift due to integration density change."`

### P7.9: 人工验证 Notebook + 文档 + 最终检查

1. `Simulation_sqc.ipynb` 重启 kernel 从头跑；重点验证 Cell 15 (Cryoscope 60ns 后无 flatline)、Cell 29-31 (阶跃响应形态正常)
2. 更新 `docs/architecture.md`：
   - §10 新增「v2.0：统一 global 时间轴策略」章节
   - §7 扩展指南补 §2.1 的标准模板
   - 末尾文档维护表追加：`2026-05-15 v1.x → v2.0：统一 mesolve 时间轴到 t_global；Pulse/FluxSignal 加 trigger 参数`
3. 最终验证：
```bash
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m pytest tests/ -v
git diff --stat master..HEAD
git diff master -- src/   # 必须为空（R1）
grep -r "np.linspace" sqc/ --include="*.py" | grep -v "# OK"   # 只剩非时间轴用途
```

---

## §6 验证标准

| # | 标准 |
|---|---|
| 1 | `pytest tests/ -v` 全绿（含重新生成的 7 个 regression baseline） |
| 2 | `git diff master -- src/` 输出为空（R1） |
| 3 | `grep -r "np.linspace" sqc/ --include="*.py"` 只剩非时间轴用途 |
| 4 | 用户在 notebook 用正弦/wave-packet 跑 Cryoscope，60ns 后不再翘尾或 flatline |
| 5 | 用户构造新实验无需 `ctrl.t_list -= offset`，只需 `Pulse(..., trigger=30)` + `pulse.hamiltonian_on(t_global)` |
| 6 | `docs/architecture.md` §10 有新章节说明统一策略，版本号 → v2.0 |

---

## §7 风险与缓解

| 风险 | 缓解 |
|---|---|
| t_global 默认 900pt vs 当前 ~140pt，cryoscope 跑慢 5-6× | 用户已知接受。若不可忍受可缩短默认 t_global——单独 PR |
| `qubit.qubit_in_mag` 在 `src/` 中（R1 不可改），但需要它在 global 轴上工作 | 调用前把 flux 投影到 t_global，构造 `FluxSignal(type=8, t_list=t_global, signal=...)` 后再传 |
| 改动期间无 baseline 回归保护 | 用单元测试逐 Phase 卡，最后一步集中重生成 |
| `Pulse.hamiltonian` 中复数 coeffs（rotating frame）的 interp | real/imag 分别 `np.interp` 后合成 |
| 改 `CompositePulse.get_t_list` 去 `1e-9` 可能影响 `get_kernel`（deprecated） | 该路径已 deprecated，P7.5 整段 deprecate；若回归测试触发再单独修 |
| `IQReadoutModel.measure` 接口变化可能影响 `IQ_readout_legacy` | 通过现有签名保持不变，内部转发；测试覆盖 |
| 范围比预期大（~17 文件迁移），可能引入意外耦合 | 严格按 P7.1→P7.9 顺序，每 Phase 跑测试 + commit checkpoint |
| Workflow 层可能有未注意到的 mesolve 依赖路径 | P7.6 先纯审计后再动 |
| 单元测试中 `len(t_list) == specific_number` 类硬断言较多 | P7.7 集中处理，改为鲁棒形态断言 |
| Notebook 数值微偏可能让用户误以为重构破坏物理 | P7.9 重启 kernel 重跑；约定 "形态正常 + 数值偏差 < 1e-4 相对" 为通过线 |
| equivalence 测试容差不够 | 必要时放宽到 1e-4，commit message 解释 |

---

## §8 不做的事

- **不改 `src/`** (R1 硬约束)
- **不改 `Pulse.get_kernel` / `CompositePulse.get_kernel`**（已 deprecated）
- **不改 `web_demo.py` / `web_demo_v2.py`**（demo 用上层 Experiment API，除非 import 失败否则不动）
- **不修改 `CONFIG.pulse.t_global` 的默认范围**（900pt 性能问题作为独立 PR）
- **不修改 Tier C 文件的 `np.linspace`**（参数扫描用途，非时间轴，不违反 R9）

---

> **更新记录**: 2026-05-15 — 初次编写。基于用户 Cryoscope 反演翻车诊断 + sqc/ 全代码库时间轴审计。配套详细 plan 在 `.claude/plans/unified-kindling-metcalfe.md`。
