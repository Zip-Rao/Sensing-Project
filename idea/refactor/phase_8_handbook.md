# Phase 8 — Filter function 代码适配：零遮盖 → 核函数卷积

> **状态**: 规划阶段，待实施。将 delay Ramsey 和 Ramsey sensing 的"脉冲期 flux 零遮盖"近似替换为完整的 filter function 卷积模型。
> **前置**: P0–P7 全部完成（P7 的时间轴统一是理想前提，但 P8 可与 P7 并行实施——P8 的 filter function 计算和卷积重建不依赖 t_global）。
> **理论配套**: [`_sensing theory.md`](../_sensing%20theory.md) "核函数与卷积测量：统一理论框架"节、"瞬态磁场协议 → 从哈密顿量到卷积结构"节。
> **配套文件**:
> - 重构总方案: [_refactor_plan.md](_refactor_plan.md) §8、§13
> - 当前交接状态: [_handoff_state.md](_handoff_state.md)

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
| [sqc/reconstruction/delay_ramsey.py](sqc/reconstruction/delay_ramsey.py) | `DelayRamseyReconstruction` 新增 `method: Literal["slope", "wiener"] = "slope"`；`method="wiener"` 时从 result 取 kernel 做 Wiener 反卷积；新增 `lambda_reg` 参数 |
| [sqc/reconstruction/ramsey.py](sqc/reconstruction/ramsey.py) | 同 delay_ramsey，新增 wiener 反卷积路径 |
| [sqc/reconstruction/kernel.py](sqc/reconstruction/kernel.py) | 无改动（已可复用）；确认 `KernelEstimator` 的默认参数从 CONFIG 读取 |

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

- **P7（时间轴统一）**：非硬依赖。P8 可独立实施；若 P7 先完成，P8 的触发对齐代码可简化。
- **KernelEstimator**：已存在且可复用。
- **Wiener 反卷积**：`transient.py` 中已有实现，可直接复用或提取。

### 风险

1. **kernel 缓存一致性**：kernel 依赖 `(t_rabi, tau_R, omega_d, qubit_spec)`，若 qubit 参数中途改变需刷新缓存
2. **性能**：kernel 计算需对每个时间点做两次 mesolve，耗时与脉冲长度成正比——但对 20 点 t_rabi 的 Ramsey 序列（~3×20=60 点），开销可接受
3. **数值稳定性**：Wiener 反卷积在小信号（深 tail 区）可能放大噪声；λ 正则化参数需合理默认值

---

## §7 不做的事

- **不修改** Echo、CPMG、Cryoscope、π 脉冲补偿的 flux 处理逻辑（它们不需要）
- **不引入**新的依赖库
- **不删除**现有的 zero-out 路径（保持向后兼容）
- **不修改** `src/` 下任何文件（R1）
- **不统一**所有协议到一个"万能反卷积框架"——保持每个协议的重建类独立
