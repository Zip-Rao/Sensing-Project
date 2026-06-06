# Phase 11 Handbook — 瞬态测频 order≥3 反演的 G₃ 修复(双路线可选)

> 日期: 2026-06-06 | 状态: **DRAFT / 待批准**(尚未开 PR,未动源码)
> 起草: 主 session | 关联:
> [phase_10_kernel_theory_and_fix.md](phase_10_kernel_theory_and_fix.md) §3,
> [transient_frequency_theory.md](transient_frequency_theory.md) §6,
> [`result/transient_error/高阶核函数差异根因.md`](../../result/transient_error/高阶核函数差异根因.md)
> 必读启动顺序见 [CLAUDE.md](../../CLAUDE.md) R2。

---

## 0. 一句话目标

修复 `sqc/calibration/frequency.py:_measure_frequency_transient` 中 order≥3
cubic-Newton 用错 G₃ 的概念 bug:把"对角核积分 $G_3^{\rm diag}$"换成"真·立方
Taylor 系数 $G_3^{\rm Taylor}$"。提供 A/B 两条路线,仿真时可自由选择。

---

## 1. 问题陈述(根因)

P10 把核标定做对了,但 `frequency.py` 反演恒定 detuning 时用错了物理对象:

| 量 | 定义 | 数值(EC=0.2·2π, EJ=15·2π, 方波 π/2=10ns) | 反演恒定 Δ 是否该用 |
|---|---|---|---|
| $G_3^{\rm diag}=\int k_3^{\rm diag}(t)\,dt$ | 单时刻探针三阶差分 | $\approx5.4$,$\lvert G_3^{\rm diag}/G_1\rvert\approx0.86$–$1.0$ | ❌ 否 |
| $G_3^{\rm Taylor}=\dfrac{d^3p_{\rm diff}}{d\Delta^3}\big\rvert_0$ | $p_{\rm diff}=G_1\Delta+\frac16 G_3\Delta^3$ 的立方系数 | $\approx1057$,$G_3/G_1\approx-167\,\mathrm{ns}^2$ | ✅ 是 |

恒定 Δ 的三阶响应是**三重时间积分** $\iiint k_3(t_1,t_2,t_3)$,对角项只是零测度切片;
两者差约 $(2T_{\pi/2})^2$。$G_1$ 两路径一致(−6.31),故仅三阶受影响。数值证明见
`result/transient_error/diag_vs_taylor_G3.py`(纯 numpy)与 `verify_G3_qutip.py`(真核)。

---

## 2. 两条修复路线(本 handbook 的核心:都做成可选项)

### 路线 A — 直接拟合 $p_{\rm diff}(\Delta)$(推荐默认)

标定阶段施加一组**已知** $\Delta$(仿真 trivial),测 $p_{\rm diff}$,对奇多项式
$\{\Delta,\Delta^3,\Delta^5\}$ 拟合,立方系数 ×6 得 $G_3^{\rm Taylor}$。

- 标定一次性成本: ~2×15 次 mesolve;结果作为脉冲序列固有常数缓存复用。
- 测量时成本: **零额外**(仍只测一个 $p_{\rm diff}$,Newton 复用缓存的 $G_1,G_3$)。
- 优点: 定义上即正确,不依赖对易子推导,无三重积分对称因子歧义。
- 缺点: 标定需能施加已知失谐(真机=已知 Δ 的 Ramsey 标定,本就要做)。

### 路线 B — 多时间核三重积分

把 `KernelEstimator` 三阶部分从单时刻探针升级为三时刻探针,显式算
$k_3(t_1,t_2,t_3)$ 再做时序楔形区三重积分(P10 §2.4 的 $k_n^{\rm ord}$ 已给定义)。

- 标定一次性成本: $O(N^3)$ 次 mesolve(N≈40 → 数千次)。
- 优点: 与现有核框架统一,纯 VZ 探针,真机无需额外施加已知 Δ。
- 缺点: 贵;三重积分对称因子($1/3!$、$t_1{\le}t_2{\le}t_3$ 排序)易差常数,
  **必须用路线 A 的值标定校验**。

### 选择接口建议

`FrequencyMeasurement` / `_measure_frequency_transient` 增加参数
`g3_source: Literal["fit", "kernel_triple", "diag_legacy"] = "fit"`:
- `"fit"` → 路线 A(默认)
- `"kernel_triple"` → 路线 B
- `"diag_legacy"` → 保留旧 $G_3^{\rm diag}$ 行为(仅向后兼容/诊断,文档标注其物理不正确)

$k_3^{\rm diag}$ **保留**——它回答的是另一个合法问题(单点 detuning 脉冲三阶灵敏度),
只是不再默认用于反演恒定失谐。

---

## 3. 任务清单

- [ ] T1. 在 `sqc/calibration/frequency.py` 加 `g3_source` 参数;`"fit"` 路线实现
      奇多项式拟合(扫已知 Δ),产出 $G_3^{\rm Taylor}$,喂入现有 Newton(Newton 本身不改)。
- [ ] T2. 缓存机制:$G_1,G_3^{\rm Taylor}$ 作脉冲序列固有量缓存,测量时复用(对齐核"测一次复用")。
- [ ] T3.(可选)路线 B:`KernelEstimator` 三时刻探针 + 三重积分;用 T1 的值断言一致。
- [ ] T4. 更正 `transient_frequency_theory.md` §6:区分 $G_3^{\rm diag}$ 与 $G_3^{\rm Taylor}$,
      安全区改为实测 $\Delta^\*\approx17$ MHz(§R11 文档同步)。
- [ ] T5. 回归: `pytest tests/regression -m regression` 100% 通过;`git diff --quiet master -- src/`。
- [ ] T6. 新增单测: order=3 在安全区 RMS 显著优于 order=1(用 `result/transient_error` 的判据)。

---

## 4. 硬约束提醒(R1/R3/R8/R11)

- **不动 `src/`**(R1);所有改动在 `sqc/`。
- 若改 order=3 默认行为 → 属"重大改动"(R8),先 commit 安全点,逐模块测试。
- 若数值 baseline 变化 → R3 不准放宽 `rtol=1e-6`;在 commit message 标
  `intentional physics change`。
- 改 API 签名 / 新约定 → 必须同步 `docs/architecture.md`(R11)。

---

## 5. 验收

1. `g3_source="fit"` 时,order=3 在 $\lvert\Delta\rvert\lesssim0.8\Delta^\*$ 内
   RMS 误差较 order=1 降 ≥5×(对照 `瞬态测频误差来源分析.md` §4.2)。
2. 若实现路线 B,其 $G_3$ 与路线 A 的 $-167\,\mathrm{ns}^2$ 相对偏差 <5%。
3. 全回归绿;`src/` diff 为空;文档已同步。
