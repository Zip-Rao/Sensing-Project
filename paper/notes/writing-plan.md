# 短脉冲序列参与频率标定：当前写作计划

> 本文件是写作决策和证据边界，不是论文正文。正式文字仍受
> [`../../WRITING_BOUNDARY.md`](../../WRITING_BOUNDARY.md) 约束。2026-08-11
> 修改前的完整讨论稿已归档至
> `paper/archive/pre-v2-state-machine-2026-08-11/notes/writing-plan.md`。

## 1. 核心问题

全文围绕一个限定问题展开：短脉冲序列能否在明确的局部有效范围内
参与超导量子比特频率标定。

响应核、局部三阶反演、阻尼割线更新和 drive-reference tracking 是回答
这一问题的机制，不应扩展成“短脉冲普遍替代 Ramsey 或 spectroscopy”的
结论。

## 2. 当前协议

正式状态机为：

```text
Acquire -> Track -> Verify -> Lock
             ^        |        |
             |        |        |
             +--------+        |
                 Reacquire <----+

Any state -- budget/interlock/cancel --> SafeStop / SafeHold
```

状态职责：

- **Acquire**：宽范围测量，建立绝对频率参考和不确定度。
- **Track**：唯一允许修改物理磁通偏置的状态；使用局部短脉冲估计器、
  阻尼割线更新和 drive tracking。
- **Verify**：冻结候选偏置，使用独立频率测量确认；Track 的连续命中不能
  计入 Verify。
- **Lock**：保留已验证偏置并进行长期低成本监测；不直接调节偏置。
- **Reacquire**：局部范围、参考或置信度丢失后重新建立绝对参考。
- **SafeStop**：预算耗尽、interlock、取消或不可恢复故障后的安全保持。

Lock 中的低成本 monitor 和周期 Ramsey audit 都先进入 Verify。只有大幅跳变
或参考丢失直接进入 Reacquire。

## 3. 判据与单位

理论公式统一使用角频率，普通频率只用于图表和数值报告：

```text
f = omega / (2 pi)
delta_f = Delta / (2 pi)
epsilon_f = epsilon_omega / (2 pi)
```

不得把 rad GHz 数值直接标为 MHz。

主要统计量：

```text
U_tra = |r_tra| + z_(1-beta) sigma_tra
U_ver = |r_ver| + z_(1-beta) sigma_ver
U_loc = |r_loc| + z_(1-beta) sigma_loc
```

阈值顺序保持：

```text
epsilon_hold < epsilon_final < epsilon_enter < Delta_val
```

转换规则：

- Acquire 或 Track 的候选满足 `U <= epsilon_enter` 后进入 Verify。
- Verify 在同一冻结偏置上连续 `N_verify` 次满足
  `U_ver <= epsilon_final` 后进入 Lock。
- Verify 未通过但局部条件仍可靠时返回 Track；歧义或越界时进入 Reacquire。
- Lock 的 clear/suspect/large-jump 阈值与 `epsilon_hold` 的物理目标分开定义。
- Verify attempts、Verify shots、总命令、总 shots、solver calls、wall time 和
  Reacquire 次数都受预算约束。

当前确定性后端返回 `sigma=0`，并标记
`uncertainty_source="deterministic_zero"`。这不是实验置信区间。真机执行器必须
明确 `z_(1-beta)`、sigma 估计方法、覆盖率和多重检验策略。

## 4. 四类数据

全文严格区分：

1. **response-calibration data**：选择探针设置，拟合 `G1`、`G3` 和 nuisance
   corrections。
2. **response-validation data**：独立确定 `Delta_val`，检查分支连续性和局部
   反演误差。
3. **closed-loop independent-verification data**：在冻结候选偏置上判断 Verify；
   不得复用 Track convergence streak。
4. **long-term monitoring data**：进入 Lock 后的 monitor 与 Ramsey audit 时间序列。

解析 transmon 色散用于确定性仿真的 ground-truth scoring，不属于第三或第四类
实验数据。

## 5. 章节分工

### Abstract / Introduction

- 保留“探究可能性”的限定叙事。
- 数值残差称为 analytically checked residual。
- 不称为 independent Verify，也不称 experimental speedup。

### Section II: System and calibration task

- II.A：器件、目标和角频率记号。
- II.B：scanned Ramsey、fixed-delay Ramsey 和 bipartite short-pulse probe。
- II.C：完整 V2 状态机、误差量、阈值层级、数据独立性和评价接口。
- Fig. 1(c) 承担状态机总览，不另加第四面板。

### Section III: Short-pulse frequency inference and calibration

- III.A：短脉冲差分响应和 `G1/G3`。
- III.B：局部一阶与三阶估计器、fold 和 validated interval。
- III.C：候选进入、blind step、secant sensitivity 和 flux-bias update。
- III.D：一般下一步失谐、fixed-drive 对照、drive tracking 和范围守卫。
- 本章只实例化 Track；Verify、Lock 和 Reacquire 的总职责不重复展开。

旧 `04_feedback_protocol.tex` 与 III.C--III.D 重复，不再由 `main.tex` 输入；文件
保留，修改前版本也在 archive 中。

### Numerical results

- F2/F4/F5 只作为确定性数值证据。
- F4 是获得初始 seed 后的 local Track 案例，不是完整 V2 状态机运行结果。
- 984/4000 是 instrumented solver-call comparison，不是实验 acquisition cost。
- 0.00519 MHz 是解析色散复算残差，不是独立 Verify 测量。

### Discussion / Conclusion

- 说明局部有效范围、fold、drive tracking 成立条件和成本定义。
- 明确尚无有限 shots、readout error、drift、pulse distortion 和 leakage 验证。
- 不宣称长期稳定、闭环带宽、PSD 抑制或 Allan stability。

### Appendix

- 给出数值 benchmark 参数和 solver-call accounting。
- 给出状态机阈值、`N_verify`、最大 Verify attempts/shots、monitor/audit 周期、
  数据角色、预算和恢复语义。
- 区分 benchmark 的 0.016-MHz common tolerance 与 V2 Verify threshold。

## 6. 图表与篇幅

- 主文理论部分目标 4--4.5 页，硬上限约 5 页。
- 删除重复 Section IV，为 Verify/Lock 定义释放版面。
- Results、Discussion 和 Conclusion 不因状态机扩展而压缩。
- 优先维持全文 11 页；若增加到 12 页，新增内容应主要来自必要 Appendix，
  不能让结果图进一步脱离对应文字。
- Fig. 1 组合图必须使用最新 `fig1d-v8-editable.svg`，最终交付稳定命名的
  editable SVG、PDF 和 PNG。

## 7. 当前证据边界

已经支持：

- 短脉冲局部频率信息和三阶校正的确定性数值案例。
- seeded Track、drive tracking、范围/灵敏度/置信守卫。
- Acquire/Track/Verify/Lock/Reacquire/SafeStop 软件状态机。
- 命令/事件、ReasonCode、预算预检、取消、checkpoint 和 at-least-once 恢复。
- 旧 workflow 接口兼容。

尚未支持为论文实证结论：

- 真机有限-shot Verify 统计。
- 正的 hardware monitor/audit 周期及长期运行数据。
- 实验 acquisition-time speedup。
- 长期稳定性、闭环带宽、PSD 或 Allan stability。
- 对模型误差、pulse distortion、readout error 和 leakage 的系统鲁棒性。

## 8. 交付检查

- [ ] `eq:target-entry-condition` 唯一定义且全部引用解析。
- [ ] `eq:drive-tracking` 只有一个 active label。
- [ ] `mckay2017efficient` 引文解析。
- [ ] 无旧 `Acquire--Track--Lock--Reacquire` 正文。
- [ ] 无 `N_lock` 作为 Verify 次数。
- [ ] 组合 Fig. 1 使用最新 panel (d)。
- [ ] `main.log` 和 `build/main.log` 无 undefined/multiply-defined/citation warning。
- [ ] 最终 PDF 逐页检查，理论不超过约 5 页，结果图文关系可接受。
- [ ] 相关测试、完整 pytest、Ruff、双语 Sphinx 和 `git diff --check` 有明确记录。
