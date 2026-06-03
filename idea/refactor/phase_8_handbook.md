# Phase 8 — Filter function 代码适配：零遮盖 → 核函数卷积

> **状态**: 规划阶段，待实施。将 delay Ramsey 和 Ramsey sensing 的"脉冲期 flux 零遮盖"近似替换为完整的 filter function 卷积模型。
> **前置**: P0–P7 全部完成（P7 的时间轴统一是理想前提，但 P8 可与 P7 并行实施——P8 的 filter function 计算和卷积重建不依赖 t_global）。
> **与 P10 的关系（执行顺序约定，方案 A）**: **P10 先实施，P8 后实施**。具体 gate 与契约见 §1.0。
> **理论配套**: [`_sensing theory.md`](../_sensing%20theory.md) "核函数与卷积测量：统一理论框架"节、"瞬态磁场协议 → 从哈密顿量到卷积结构"节。
> **配套文件**:
> - 重构总方案: [_refactor_plan.md](_refactor_plan.md) §8、§13
> - 当前交接状态: [_handoff_state.md](_handoff_state.md)
> - 上游 phase: [phase_10_kernel_extension_handbook.md](phase_10_kernel_extension_handbook.md)（核函数三维扩展，P8 启动的前置）

---

## §1.0 P8 ↔ P10 协调（方案 A：P10 先 → P8 后）

### 决策

**P8 必须在 P10 的关键里程碑完成后才启动**。原因：

1. **P8 是 KernelEstimator 的消费者** — P8 §3 Tier A4 已明确"`kernel.py` 无改动"；P8 只调用 estimator。
2. **避免 κ workaround 扩散** — 若 P8 早于 P10 实施，`delay_ramsey.py` / `ramsey.py` 的 Wiener 反卷积只能复用 [frequency.py:294-302](../../sqc/calibration/frequency.py#L294-L302) 的 flux kernel + κ workaround anti-pattern。P10 先实施后，P8 可以直接用 `mode='omega'`，新代码从一开始就干净。
3. **P10 的默认行为 byte-equivalent** — P10 设计上保证 `KernelEstimator()`（默认参数）在任何时刻都与当前实现等价；这让 P8 可以放心调用，不会因 P10 中间状态破坏 baseline。

### P8 启动 gate

启动 P8 之前，**必须确认以下 P10 里程碑已合入 master**：

| P10 里程碑 | 给 P8 的能力 | 是否阻塞 P8 启动 |
|-----------|------------|----------------|
| **Phase 10.1**（mode 维度 + Virtual Z） | `KernelEstimator(mode='omega', method='exp', order=1)` 可用 | **阻塞**（P8.3/P8.4 用此路径） |
| **Phase 10.5**（frequency.py 迁移，消除 κ workaround） | `pulse.get_kernel()` shim 稳定，调用模式可参考 | **强烈建议先合入**（P8 沿用一致的调用模式） |
| Phase 10.2（sim 方法） | sim mode | 不阻塞（P8 默认走 exp） |
| Phase 10.3（高阶 + 序列化） | order≥2、`KernelResult.save/load` | 不阻塞（P8 默认 order=1） |
| Phase 10.4（Hammerstein-Volterra） | 高阶反卷积 | 不阻塞（P8 用 order=1 Wiener） |

**最小 gate**：**Phase 10.1 已合入** → P8 可启动；**Phase 10.5 已合入** → P8 才进入 P8.3/P8.4 的标定 + 反卷积层。

### P8 实施时对 P10 的硬约定

1. **优先使用 `mode='omega'`** — P8.3 / P8.4 的 Wiener 反卷积**必须**走 omega kernel + `_omega_to_flux` 反演路径，而不是 flux kernel + 手动 κ。与 P10.5 frequency.py 的做法保持一致。
2. **不要在 P8 代码里写新的 κ workaround** — 如果发现协议层需要 κ 换算才能完成 wiener 反卷积，**停手并报告**（说明 P10 设计有 gap，需要走 Phase 11 立项扩展 estimator，而不是在 P8 里 patch）。
3. **kernel 缓存 key 必须包含 `KernelEstimator` 的非默认字段** — §6.1 风险 1 提到的 transfer_function hash 之外，还要加上 `mode`、`method`、`order`、`virtual_z_impl` 等字段。**实操**：直接用 `KernelResult.metadata` 字典序列化后做 hash 当 cache key。
4. **不动 `sqc/reconstruction/kernel.py`** — P8 实施者若发现需要在 estimator 里加新功能，**停手并升级到 Phase 11 单独立项**。P10 之后的 estimator 是 P8 的依赖，不是 P8 的修改对象。
5. **不动 `KernelResult` 数据结构** — 同上；P8 只读不写。

### P10 给 P8 的契约（P10 实施者必须保证）

1. **默认参数 byte-equivalent**：`KernelEstimator()`（无参数）在 P10 的任何子 phase 完成后，都等价于当前 [sqc/reconstruction/kernel.py](../../sqc/reconstruction/kernel.py) 的行为。
2. **`pulse.get_kernel()` shim 始终返回 `(t_samples, kernel_ndarray)`**。
3. **`KernelResult.k1` 始终是 1D ndarray**（P8 不需要区分 order=1 vs order≥2）。
4. **`reconstruction/__init__.py` 只增不减**，P10 新导出的符号不能改名或删除。

### 并发场景

短期内不期望 P8 和 P10 并行开发。若出现：

- P8 分支已切出（基于 P10 启动前的 master）→ P10 合入后 P8 必须 rebase 到包含 P10.1 + P10.5 的 master，再继续
- 严禁 P8 抢先合入 master（会污染 P10.5 frequency.py 迁移的 review）

具体协议层 vs P9（控制路由）的 merge 协调见下文 §6.1。

---

## §1 问题诊断

### 1.1 哪些协议受影响

只有 flux 信号**连续存在、无法与脉冲在时间上分离**的协议需要 filter function：

| 协议 | 当前做法 | 问题 |
|------|---------|------|
| **瞬态磁场** (`transient.py`) | ✅ 已使用 kernel + 卷积 | — |
| **delay Ramsey** (`delay_ramsey.py`) | ❌ 脉冲期 flux 强制 zero-out | 不符合物理实际；τ_R 短时脉冲边沿贡献显著 |
| **Ramsey sensing** (`ramsey.py`) | ❌ 用调制函数 $y(t)$，忽略脉冲边沿 | 当 τ_R 与 T_π/2 可比时，灵敏度过渡区不可忽略 |

以下协议**不受影响**（flux 可通过时序设计限制在自由演化期内）：

| 协议 | 原因 |
|------|------|
| Cryoscope | flux 在第一 π/2 之后触发、第二 π/2 之前截断 |
| Echo / CPMG | τ ≫ T_π，脉冲边沿贡献可忽略 |
| π 脉冲补偿 | 关心共振条件（二值判断），不依赖脉冲期相位细节 |
| Ramsey 频率标定 | τ_R ∼ μs ≫ T_π/2，且信号为准静态 |

### 1.2 delay Ramsey 的具体症状

当前 `DelayRamseyExperiment.run()`（[sqc/experiments/delay_ramsey.py:134-143](sqc/experiments/delay_ramsey.py)）：

```python
signal = np.zeros(len(t_sig), dtype=float)
free_mask = (t_sig >= self.t_rabi[-1]) & (t_sig <= self.t_rabi[-1] + self.tau_R)
signal[free_mask] = np.array([...])  # 只在自由演化窗放 flux
```

`DelayRamseyCalibration.calibrate()`（[sqc/reconstruction/delay_ramsey.py:146-149](sqc/reconstruction/delay_ramsey.py)）同样 zero-out。

**后果**：
1. 标定曲线 φ_cal(z) = κ·τ_R·z 假设了灵敏度全在自由演化期——实际 filter function 面积 ≠ τ_R
2. 测量步重建 Φ_tail = φ / (κ·τ_R) 忽略脉冲边沿对相位的贡献
3. 时间映射偏移 bug（flux signal time 未减 t_rabi[-1]）叠加

### 1.3 Ramsey sensing 的具体症状

`RamseyExperiment` 目前只用简单调制函数 y(t) = 1 在自由演化期。当 τ_R 短、信号快时，脉冲边沿区域的灵敏度过渡不可忽略。

---

## §2 目标状态

### 2.1 核心变更

| 组件 | 旧行为 | 新行为 |
|------|--------|--------|
| 实验层 | 脉冲期 flux zero-out；只输出 φ | 连续施加 flux；输出 φ + 可选 `kernel` |
| 标定层 | φ_cal(z) = κ·τ_R·z（仅自由演化） | 标定也连续加 flux；标定曲线的斜率 = κ·∫k(t)dt |
| 重建层 | Φ = φ / k（除法） | 默认保持除法（向后兼容）；新增 `method="wiener"` 反卷积路径 |
| filter function | 不存在于非 transient 协议中 | `KernelEstimator` 复用到 delay Ramsey / Ramsey |

### 2.2 用户 API（预期）

```python
# delay Ramsey — 默认行为不变（向后兼容）
exp = DelayRamseyExperiment(qubit=qubit, flux_signal=tail, tau_R=10.0)
result = exp.run()  # zero-out，与现有 baseline 一致

# delay Ramsey — 启用 filter function
exp = DelayRamseyExperiment(qubit=qubit, flux_signal=tail, tau_R=10.0,
                            use_filter_function=True)
result = exp.run()  # 连续 flux，result.data["kernel"] 包含 filter function

# 重建 — 反卷积路径
recon = DelayRamseyReconstruction(method="wiener", lambda_reg=5.0)
flux = recon.reconstruct(result)  # 使用 result.data["kernel"] 做 Wiener 反卷积
```

### 2.3 filter function 的存储与复用

- `KernelEstimator`（[sqc/reconstruction/kernel.py](sqc/reconstruction/kernel.py)）已可计算任意 pulse 的 kernel
- 对于给定的 `(t_rabi, tau_R, omega_d, qubit_spec)`，kernel 是确定的——应缓存，不必每次 run 重算
- kernel 存入 `ExperimentResult.data["kernel"]`，重建时直接取用

---

## §3 受影响文件

### Tier A1: 实验层

| 文件 | 改动内容 |
|------|---------|
| [sqc/experiments/delay_ramsey.py](sqc/experiments/delay_ramsey.py) | 新增 `use_filter_function: bool = False`；True 时去掉 zero-out，连续施加 flux_window；修复时间映射偏移 `t - t_rabi[-1]`；计算并存储 kernel |
| [sqc/experiments/ramsey.py](sqc/experiments/ramsey.py) | 新增 `use_filter_function: bool = False`；True 时输出 kernel |

### Tier A4: 重建层

| 文件 | 改动内容 |
|------|---------|
| [sqc/reconstruction/delay_ramsey.py](sqc/reconstruction/delay_ramsey.py) | `DelayRamseyReconstruction` 新增 `method: Literal["slope", "wiener"] = "slope"`；`method="wiener"` 时从 result 取 kernel 做 Wiener 反卷积；新增 `lambda_reg` 参数。**必须用 `KernelEstimator(mode='omega', method='exp', order=1)`**（P10 协调约定，§1.0）；不允许复制 frequency.py 的 κ workaround |
| [sqc/reconstruction/ramsey.py](sqc/reconstruction/ramsey.py) | 同 delay_ramsey，新增 wiener 反卷积路径。同样**走 omega kernel + `_omega_to_flux` 反演**，不写 κ workaround |
| [sqc/reconstruction/kernel.py](sqc/reconstruction/kernel.py) | **无改动**（P10 已实施完毕，estimator 已扩展）；P8 仅作为消费者调用 |

### Tier A5: 标定层

| 文件 | 改动内容 |
|------|---------|
| [sqc/reconstruction/delay_ramsey.py](sqc/reconstruction/delay_ramsey.py) `DelayRamseyCalibration` | `calibrate()` 新增 `use_filter_function: bool = False`；True 时标定信号连续跨整个序列（不 zero-out），标定曲线斜率为 κ·∫k dt 而非 κ·τ_R |

### 文档

| 文件 | 改动内容 |
|------|---------|
| [docs/architecture.md](docs/architecture.md) | 新增 P8 协议扩展说明；更新 delay Ramsey + Ramsey 的 API 文档 |
| [idea/_sensing theory.md](../_sensing%20theory.md) | 已完成（"核函数与卷积测量：统一理论框架"） |

---

## §4 执行计划

### P8.0: 安全检查点

- [ ] `git add -A && git commit -m "safety checkpoint: before P8 filter function"`
- [ ] `pytest tests/regression -m regression` 确认全绿

### P8.1: delay Ramsey — 实验层

- [ ] `DelayRamseyExperiment` 加 `use_filter_function: bool = False`
- [ ] 修复时间映射偏移：`value_at(t_fall + t_d + t - t_rabi[-1])`
- [ ] `use_filter_function=True` 时：去掉 free_mask zero-out，将整个 flux 窗口连续施加；flux 信号窗口改为 `signal = flux_signal.samples` 对 `t_sig + t_start` 插值
- [ ] 计算 kernel（调用 `KernelEstimator`，若缓存命中则跳过）并存入 result
- [ ] `use_filter_function=False` 保持现有行为（向后兼容）

### P8.2: delay Ramsey — 标定层

- [ ] `DelayRamseyCalibration` 加 `use_filter_function: bool = False`
- [ ] True 时：标定信号 `z` 连续跨整个 t_sig（不 zero-out 脉冲期）
- [ ] 标定曲线输出 φ(z) 而非 φ_shifted(z)，线性拟合的斜率对应 κ·∫k dt

### P8.3: delay Ramsey — 重建层

- [ ] `DelayRamseyReconstruction` 加 `method: Literal["slope", "wiener"] = "slope"` 和 `lambda_reg: float`
- [ ] `method="slope"`：保持当前标定斜率除法
- [ ] `method="wiener"`：从 result 取 kernel，调用 Wiener 反卷积（可复用 `sqc/reconstruction/transient.py` 的 Wiener 逻辑，或抽公共函数）
- [ ] 新增 `_via_wiener()` 方法

### P8.4: Ramsey sensing — 实验 + 重建

- [ ] `RamseyExperiment` 加 `use_filter_function: bool = False`
- [ ] True 时：连续施加 flux 信号，输出 kernel
- [ ] `RamseyReconstruction` 加 `method="wiener"` 路径
- [ ] 修复时间映射（同理减 t_rabi[-1]）

### P8.5: 公共 Wiener 反卷积函数提取（可选优化）

- [ ] 检查 `transient.py`、`delay_ramsey.py`、`ramsey.py` 的 Wiener 逻辑是否可抽公共函数
- [ ] 若三处逻辑相同，抽取 `_wiener_deconvolve(phi, kernel, dt, lambda_reg) -> FluxSignal` 到 `sqc/reconstruction/base.py` 或独立模块

### P8.6: 测试

- [ ] 新增 `tests/unit/test_delay_ramsey_filter_function.py`
  - `use_filter_function=True` 与 `False` 在 τ_R ≫ T_π/2（慢变 tail）下结果接近（rtol=1e-2）
  - kernel 面积验证：∫k dt ≈ τ_R + 2/Ω
  - Wiener 反卷积 vs 标定斜率法对比
- [ ] 新增 `tests/unit/test_ramsey_filter_function.py`（类似）
- [ ] 确保现有 `use_filter_function=False` 测试全绿
- [ ] 生成新的 baseline（仅当 `use_filter_function=True` 且确定为正确行为后）

### P8.7: 文档 + Notebook

- [ ] 更新 `docs/architecture.md` 对应章节
- [ ] 在 `Simulation_sqc.ipynb` 的 delay Ramsey section 加 `use_filter_function=True` vs `False` 对比 cell

---

## §5 验证标准

- [ ] `use_filter_function=False` 所有现有测试全绿（向后兼容）
- [ ] `use_filter_function=True` + 慢变 tail：与 `False` 结果一致（rtol=1e-2）
- [ ] filter function 面积 ∫k dt 与理论值一致
- [ ] Wiener 反卷积重建的 tail 与真实 tail 的 RMSE 显著优于简单除法
- [ ] `pytest tests/regression -m regression` 全绿

---

## §6 依赖与风险

### 依赖

- **P10（核函数三维扩展）**：**硬依赖**。P10.1（mode 维度上线）是 P8 启动的最小 gate；P10.5（frequency.py 迁移）是 P8.3/P8.4 标定+反卷积的强烈建议前置。详见 §1.0。
- **P7（时间轴统一）**：非硬依赖。P8 可独立实施；若 P7 先完成，P8 的触发对齐代码可简化。
- **P9（级联预失真 + 协议驱动测量）**：非硬依赖，详见下文 §6.1。
- **KernelEstimator**：依赖 P10 之后的版本（带 `mode='omega'` 等扩展），不再使用 P10 前的单点 API。
- **Wiener 反卷积**：`transient.py` 中已有实现，可直接复用或提取。

### 6.1 与 P9 的协调

P9（[phase_9_handbook.md](phase_9_handbook.md)）在 `delay_ramsey.py` / `ramsey.py` 有文件级交集，但**行级隔离**：

- P8 改 `run()` **中段**：flux modulation 逻辑（zero-out vs 连续施加 + kernel 计算）
- P9.B 改 `run()` **开头**：`flux = self._route_flux(flux)`（注入 `control_line` 失真）

两者改的位置错开，merge 不冲突。

**情形 1：P8 先实施，P9 后实施** — P8 无需为 P9 做调整；P9.B 在 `run()` 开头插 `_route_flux()`，位置早于 P8 的 zero-out / 连续施加分歧点。

**情形 2：P9 先实施，P8 后实施** — P8 实施时 `run()` 开头已有 `_route_flux()`。**P8 需注意**：
- `use_filter_function=True` 分支里计算 kernel 的 flux signal 应当是 `_route_flux()` **之后** 的版本（即经 `control_line` 畸变后到达 qubit 的真实信号）——这是物理正确做法
- 实操：把 kernel 计算位置放在 `_route_flux()` 之后；若 P8 当前实现是基于"理想 flux signal"求 kernel，与 P9 合并后 kernel 会自动变成基于"畸变 flux signal"——这是符合物理的，无需额外代码调整
- 仅在 P8 的 baseline test 中需要确认：当调用方**不传** `control_line`（即旧行为）时，kernel 与 P9 实施前一致

**情形 3：同 phase 并行** — 建议 P8.1 / P8.4 完成后再做 P9.B.3，减少 review 复杂度。

### 风险

1. **kernel 缓存一致性**：kernel 依赖 `(t_rabi, tau_R, omega_d, qubit_spec)`，若 qubit 参数中途改变需刷新缓存。**P10 后约束**：cache key 还必须包含 `KernelEstimator` 的非默认字段（`mode`、`method`、`order`、`virtual_z_impl`）和 `control_line.transfer_function` hash（若 P9 已实施）。**实操**：用 `KernelResult.metadata` 字典序列化后做 hash 当 cache key
2. **性能**：kernel 计算需对每个时间点做两次 mesolve，耗时与脉冲长度成正比——但对 20 点 t_rabi 的 Ramsey 序列（~3×20=60 点），开销可接受
3. **数值稳定性**：Wiener 反卷积在小信号（深 tail 区）可能放大噪声；λ 正则化参数需合理默认值

---

## §7 不做的事

- **不修改** Echo、CPMG、Cryoscope、π 脉冲补偿的 flux 处理逻辑（它们不需要）
- **不引入**新的依赖库
- **不删除**现有的 zero-out 路径（保持向后兼容）
- **不修改** `src/` 下任何文件（R1）
- **不统一**所有协议到一个"万能反卷积框架"——保持每个协议的重建类独立
