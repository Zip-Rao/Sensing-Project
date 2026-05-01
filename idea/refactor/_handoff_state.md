# 重构进度交接状态 (Track A handoff state)

> **作用**:这是各 phase 子代理之间唯一的"接力棒"。每个 phase 完成时**必须**更新本文件。  
> **规则**:本文件由 [.claude/agents/refactor-phase-executor.md](../../.claude/agents/refactor-phase-executor.md) 维护;主 session(用户)只读,不直接编辑(除非要修正子代理的笔误)。  
> **配套阅读**:[`_refactor_plan.md`](_refactor_plan.md) §3(双 Track 编织)、§13(假设)。

---

## 上次更新

| 字段 | 值 |
|---|---|
| 完成 phase | P0 |
| 完成日期 | 2026-05-01 |
| commit SHA | 9c5c5f6 (final), ce3b7a9 (P0 main), ab444f6 (.gitignore) |
| 执行者 (人/agent id) | refactor-phase-executor (P0 run) |
| 本次 token 实际消耗 | ~80K |

---

## 已完成 phase 列表

- [x] **P0** — 测试基线 + 工具基础设施 (2026-05-01, commits: ce3b7a9, ab444f6, 9c5c5f6)
- [ ] **P1** — sqc/ 骨架 + ABC + 数据结构 + src_mirror/ 镜像
- [ ] **P2** — 已实现协议 (case 1/2/4) 实验对象化 + KernelEstimator 去重 + src_mirror/protocal.py 创建
- [ ] **P3a** — basis 模块 + Wiener / RamseyIQ / DiffEcho / HammersteinWiener 内化 + src_mirror/analysis.py 创建
- [ ] **P3b** — LMReconstruction 完整内化 (含伴随 Jacobian)
- [ ] **P3c** — Cryoscope/Calibration 内化 + src_mirror/protocal.py 中 Calibration 类补全
- [ ] **P4** — ControlLine + DistortionModel + PredistortionDesigner + workflow
- [ ] **P5** — TransferMatrix + 双 qubit Z-crosstalk + (可选) Cavity 表征三件套

---

## 当前 git 状态

| 字段 | 值 |
|---|---|
| 当前分支 | `项目重建-v2` |
| 最近 commit | `9c5c5f6` (Restore src/__pycache__ tracking) |
| `git rev-parse HEAD:src` | `ecc3758842f7c85e00758c67e607a9e37cd3be2b` |
| `git diff --quiet master -- src/` 是否返回 0 | ✓ (src/*.py files unchanged; __pycache__ diffs fixed in commit 9c5c5f6) |
| 未合并到 master 的 refactor 分支 | `项目重建-v2` |

---

## 测试状态

| 测试套件 | 上次结果 | 用时 |
|---|---|---|
| `pytest tests/unit -v` | 6 passed, 0 failed | 0.04s |
| `pytest tests/regression -m regression` | 4 passed, 0 failed | 24.59s |
| `pytest tests/equivalence` | (P1 后填写) | — |
| `pytest tests/integration` | (P2 后填写) | — |

baseline pickle 清单(`tests/baselines/` 内):
- [x] `qubit_static.pkl` (P0)
- [x] `ramsey_default.pkl` (P0)
- [x] `diff_echo_default.pkl` (P0)
- [x] `transient_default.pkl` (P0)
- [ ] `cryoscope_default.pkl` (P3c)
- [ ] `lm_default.pkl` (P3b)
- [ ] `transient_calib_default.pkl` (P3c)
- [ ] `predistortion_default.pkl` (P4)
- [ ] `z_crosstalk_default.pkl` (P5)

---

## DECISION_NEEDED / 已知问题

1. **n_levels discrepancy**: The handbook specifies `n_levels=3` (conftest.py fixture, generate_baselines.py, test_physics_baseline.py), but the actual working code (Simulation.ipynb, web_demo.py) uses `n_levels=2`. Using `n_levels=3` causes dimension mismatch errors because `create_ramsey_pulse()` and similar functions in `src/pulse.py` default to `n_levels=2` and don't accept a qubit parameter. P0 was completed with `n_levels=2` to match the working production code. The upgrade to `n_levels>=3` should happen in P1/P2 when sqc/ pulse constructors are refactored to accept qubit-level info.

   **Impact**: All baseline pickle files are anchored at `n_levels=2`. When P1 upgrades to `n_levels=3`, baselines will need regeneration with explicit justification in commit message.

2. **Qt GUI crash in regression tests**: `test_ramsey_default_baseline` triggers a Windows fatal exception (code 0xc0000139) from matplotlib's Qt backend because `src/protocal.py:65` calls `Phi.plot()` which creates a GUI figure. The test still PASSES (pytest catches the exception). This is a pre-existing production code issue — `Protocal.evolve()` has side-effect plot calls in a headless context. Not fixed in P0 (would require modifying src/). Recommend setting `matplotlib.use('Agg')` at module top of `test_physics_baseline.py` in P2+ or adding a `--no-plot` flag to `Protocal` in the sqc/ refactor.

3. **gradio upgrade**: Upgraded from gradio 3.24.1 to 4.44.1 per handbook's `gradio>=4.0,<5.0` pin. `web_demo.py` imports verified OK with gradio 4.x. Full end-to-end GUI test deferred to user.

---

## Track B 当前进度(供 Track A 决定何时启动 P3)

参考 [`_TODO_master.md`](../_TODO_master.md)。Track B 的进度状态由 Track B 维护,Track A 只读取以判断 P3 是否可以启动。

| Track B 任务 | 状态 (○/△/✓) | 影响 Phase |
|---|---|---|
| 0.1 case 1 死代码清理 | ○ | P2(无强依赖) |
| 0.2 kernel 自动校准 | ○ | P2(若完成需重生 transient_default.pkl) |
| 0.3 LM 收敛修复 | ○ | **P3b(硬依赖)** |
| 1.1 Cryoscope (case 6/7) | ○ | **P3c(硬依赖)** |
| 1.2 瞬态频率标定 (case 8) | ○ | **P3c(硬依赖)** |
| 1.3 DistortionModel (src/ 内基础实现) | ○ | **P4(硬依赖)** |

P3a 不依赖 Track B,可与 P2 之后立即启动。
P3b/c 必须等对应 Track B 任务完成才能启动。

---

## 下一 phase 启动前的自动检查清单

子代理在每个 phase 启动时必须能勾完以下条目,否则 abort。

### 启动 P0 之前
- [x] git 工作区干净 (尚未开始,默认满足)
- [x] `python -c "from src.qubit import TransmonQubit; from src.protocal import Protocal; print('ok')"` 能跑通
- [x] `python -c "import qutip, numpy, scipy, matplotlib; print('ok')"` 能跑通
- [x] `requirements.txt` 已审查,需补的依赖列出

### 启动 P1 之前
- [x] P0 已完成,本文件"已完成 phase"中 P0 已勾选
- [x] `pytest tests/regression -m regression` 全部通过 (4/4)
- [x] `tests/baselines/` 下至少有 4 个 pkl
- [ ] git 工作区干净 (P0 结束后有未提交的 handoff_state.md 更新;P1 启动前需提交)

### 启动 P2 之前
- [ ] P1 已完成
- [ ] `from src_mirror.qubit import TransmonQubit` 能 import 成功
- [ ] `git diff --quiet master -- src/` 返回 0
- [ ] (推荐) Track B 0.3 (LM 收敛修复) 完成 — 若未完成,P2 内不内化 LM,留到 P3b

### 启动 P3a 之前
- [ ] P2 已完成
- [ ] `tests/equivalence/test_protocal_mirror.py` 全部通过

### 启动 P3b 之前
- [ ] P3a 已完成
- [ ] **Track B 0.3 LM 收敛修复已完成**(硬依赖)
- [ ] src/analysis.py 中 LM 在 baseline 参数下能稳定收敛(由 Track B 验证)

### 启动 P3c 之前
- [ ] P3a 已完成
- [ ] **Track B 1.1 Cryoscope case 6/7 已完成**(硬依赖)
- [ ] **Track B 1.2 瞬态标定 case 8 已完成**(硬依赖)
- [ ] src/protocal.py 中 case 5/6/7/8 在 baseline 参数下能稳定输出

### 启动 P4 之前
- [ ] P3a/b/c 全部完成(或至少 P3a + 必要的部分)
- [ ] **Track B 1.3 DistortionModel 已在 src/ 内实现**(硬依赖)

### 启动 P5 之前
- [ ] P4 已完成
- [ ] `PredistortionValidationWorkflow.run()` 改善 factor > 10
- [ ] (可选)cavity 三件套需求已确认

---

## 历次执行记录(append-only,最新在最上)

每个子代理执行完后,在此追加一条记录,**不修改前面的记录**。

| 日期 | Phase | 执行者 | commit SHA | 状态 | token 消耗 (估) | 备注 |
|---|---|---|---|---|---|---|
| 2026-05-01 | P0 | refactor-phase-executor | 9c5c5f6 | ✅ DONE | ~80K | n_levels=2 (see DECISION_NEEDED #1); 4 baselines generated; src/ anchor: ecc37588 |

---

## 文档版本

| 字段 | 值 |
|---|---|
| 本文件初始化版本 | v1.0 |
| 初始化日期 | 2026-05-01 |
| 配套主方案版本 | `_refactor_plan.md` v1.1 |
