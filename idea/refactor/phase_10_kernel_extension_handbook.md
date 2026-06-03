# Phase 10 — 核函数体系扩展：flux/omega 双单位 × sim/exp 双方法 × 高阶 Volterra

> **状态**: 规划阶段，待实施。本 phase 把当前的"单一 flux-based, 单一 exp 方法, 单一 1 阶"核函数升级为三维正交的统一框架。
>
> **前置**: P0–P7 + P9 完成（P9 非硬依赖，但建议在 P9 cascade 预失真稳定后实施以避免 baseline 二次重做）。**注意 P8 是 P10 的下游**，不是前置。
>
> **与 P8 的关系（执行顺序约定，方案 A）**: **P10 先实施，P8 后实施**。理由与契约见 §0.1。
>
> **理论配套**: [`_sensing theory.md`](../_sensing%20theory.md) "核函数与卷积测量：统一理论框架"节、"瞬态磁场协议"节、"频率标定"节。
>
> **配套文件**:
> - 重构总方案: [_refactor_plan.md](_refactor_plan.md) §7.3 (KernelEstimator)、§14.1 (debt D2)
> - 当前交接状态: [_handoff_state.md](_handoff_state.md)
> - 上游 phase: [phase_8_handbook.md](phase_8_handbook.md)（filter function 适配，P10 完成后才启动）

---

## §0.1 P10 ↔ P8 协调（方案 A：P10 先 → P8 后）

### 决策

**P10 必须先于 P8 实施**。原因：

1. **P8 是 KernelEstimator 的消费者** — P8 handbook §3 Tier A4 已明确"`kernel.py` 无改动（已可复用）"；P8 只调用 estimator，不改其内部。
2. **避免 κ workaround 扩散** — 若 P8 先实施，`delay_ramsey.py` / `ramsey.py` 的 Wiener 反卷积只能用 flux kernel + κ workaround（复制 [frequency.py:294-302](../../sqc/calibration/frequency.py#L294-L302) 的 anti-pattern）。P10 先实施后，P8.3/P8.4 可以直接用 `mode='omega'` 路径，新代码从一开始就干净。
3. **不阻塞 P8 的语义任务** — P8 的核心改动是协议层（zero-out → 连续 flux），与 KernelEstimator 内部解耦无关；等 P10 只需要等 §9 Phase 10.1 完成（mode 维度上线）即可启动。

### 启动 P8 的最小前置（gate）

P8 可以在 P10 的以下里程碑完成后立刻启动：

| P10 里程碑 | 给 P8 的能力 | P8 是否需要 |
|-----------|------------|------------|
| **Phase 10.1**（mode 维度 + Virtual Z） | `KernelEstimator(mode='omega', method='exp', order=1)` 可用 | **必需**（P8.3/P8.4 用此路径） |
| **Phase 10.2**（method 维度） | sim 方法 | 可选（P8 默认走 exp，不依赖 sim） |
| **Phase 10.3**（高阶 + 序列化） | order≥2、`KernelResult.save/load` | 可选（P8 默认 order=1，缓存机制不依赖 save/load） |
| **Phase 10.4**（Hammerstein-Volterra） | 高阶反卷积 | 可选（P8 用 order=1 Wiener 即可） |
| **Phase 10.5**（frequency.py 迁移） | 消除 κ workaround、`pulse.get_kernel()` shim 稳定 | **强烈建议**（P8 调用习惯应与 frequency.py 一致） |

**最小 gate**：**Phase 10.1 + 10.5 完成** → 启动 P8。

### P10 在实施期间对 P8 的契约

P10 的 6 个子 phase 在演进过程中，**必须保证**：

1. **默认参数 byte-equivalent**：`KernelEstimator()`（无参数）在任何时刻都等价于当前 [sqc/reconstruction/kernel.py](../../sqc/reconstruction/kernel.py) 的行为。这保证 P8 启动后可以无需修改先跑通"老 zero-out + 旧 estimator 默认"路径。
2. **`pulse.get_kernel()` shim 不破坏返回签名**：始终返回 `(t_samples, kernel_ndarray)`，保 P8 可以选择走 shim 或新 API。
3. **`KernelResult.k1` 永远是 1D ndarray**：P8 取 `result.k1` 时不需要区分 order=1 vs order≥2。
4. **`reconstruction/__init__.py` 只增不减**：P10 新导出的符号（`KernelResult`、`make_kernel_estimator_for`）不能改名或删除，否则 P8 引用会断。

### P8 实施时对 P10 的约定

P8 启动时（即 P10.1 + P10.5 完成后），P8 实施者必须：

1. **优先使用 `mode='omega'`** — P8.3 / P8.4 的 Wiener 反卷积走 omega kernel + `_omega_to_flux` 反演，而不是 flux kernel + 手动 κ。这与 P10.5 frequency.py 的做法一致。
2. **不要在 P8 代码里出现新的 κ workaround** — 若发现需要 κ 换算，停手并报告（说明 P10 设计有 gap）。
3. **kernel 缓存 key 必须包含 `KernelEstimator` 的所有非默认字段** — P8 handbook §6.1 风险 1 提到 transfer_function hash；P10 上线后还需加 `mode`、`method`、`order`、`virtual_z_impl` 等字段进 cache key（建议直接用 `KernelResult.metadata` 序列化后做 hash）。
4. **不动 `kernel.py`** — P8 实施者若发现需要在 estimator 里加新功能，停手并升级到 Phase 11 单独立项；不允许通过修改 estimator 来满足 P8 的协议层需求。

### 并发场景说明（如果未来同时跑两个分支）

短期内不期望 P8 和 P10 并行开发，但若出现：

- P8 分支已切出（基于 P10 启动前的 master）→ P10 完成后 P8 必须 rebase 到包含 P10.1 + P10.5 的 master，再继续
- 反向不允许（P10 不能基于"含 P8 改动"的 master，因为 P8 改的 `experiments/`、`reconstruction/{delay_ramsey,ramsey}.py` 在 P10 不动它们的前提下不影响 P10 推进——但 review 复杂度会上升）

---

## §0 阅读路径

按以下顺序阅读，**不准跳过**：

1. **§1 问题诊断**——为什么现有 KernelEstimator 不够用
2. **§2 三维设计空间**——flux/omega、sim/exp、阶数 1/N
3. **§3 理论推导**（最长一节，含 Volterra 展开、单位换算、sim/exp 等价性证明）
4. **§4 API 设计**——统一 KernelEstimator + 配套 reconstruction
5. **§5 实现细节**——数值差分、刺激形状、对角化提取
6. **§6 反卷积策略**——一阶 Wiener vs 高阶 Hammerstein-Wiener
7. **§7 受影响文件清单**
8. **§8 测试与 baseline 策略**
9. **§9 分阶段实施 plan**（PR 拆分）
10. **§10 Open questions / 后续 TODO**

---

## §1 问题诊断

### 1.1 当前 KernelEstimator 的三个隐含假设

当前实现 [`sqc/reconstruction/kernel.py`](../../sqc/reconstruction/kernel.py:53-149) 把**三个独立的设计选择硬编码在了同一段代码里**：

| 假设 | 代码体现 | 物理意义 |
|------|---------|---------|
| **A1：单位是 flux** | [L114-120](../../sqc/reconstruction/kernel.py#L114) 用 `FluxSignal(type=3, amplitude=...Φ₀)` | 刺激物理量是磁通 $\delta\Phi$ |
| **A2：方法是"实验模拟"** | [L124-136](../../sqc/reconstruction/kernel.py#L124) 调 `qubit.qubit_under_mag()` + `qubit_under_mag_hamiltonian()` | 经过完整 $\omega(\Phi)$ 非线性 + 多能级 |
| **A3：阶数是 1**（线性近似） | [L144](../../sqc/reconstruction/kernel.py#L144) `(p_e_stim - p_e_base) / stim_area` 是单边差分，只提取一阶 | $\Delta p_e$ 对刺激幅度严格线性 |

这三个假设**互相独立**，但代码把它们粘在一起，导致：

- 想做**纯理论核函数验证**（不依赖 qubit 色散）→ 必须改代码
- 想做**频率标定**（直接得 $\delta\omega$ 不绕道）→ 必须乘 $\kappa$（[sqc/calibration/frequency.py:294-302](../../sqc/calibration/frequency.py#L294-L302) 的 workaround）
- 想做**大信号高阶建模** → 框架内无此能力

### 1.2 现有代码已经在打补丁

[sqc/calibration/frequency.py:294-302](../../sqc/calibration/frequency.py#L294-L302) 是活化石：

```python
# Unit conversion: the code's kernel is built with a *flux* stimulus
# (FluxSignal with stim_area in Φ₀·ns), so G_diff has units 1/Φ₀ ...
# Convert via κ = dω/dΦ:
#     G_freq = G_diff / κ   ⇒   Δω = p_diff · κ / G_diff
kappa = qubit.frequency_sensitivity(qubit.flux + flux)
delta_omega = p_diff * kappa / G_diff
```

这段代码用 flux kernel 算频率：1) 必须乘 $\kappa$；2) $\kappa$ 是局部线性化值，大信号下不准。**正确的做法是直接用 omega-based kernel**。

### 1.3 sim/exp 边界混淆

[`_sensing theory.md` "瞬态磁场协议" §3.3](../_sensing%20theory.md) 明确把核函数计算分两种方法：

- **数值仿真**：用 $H_{\text{stim}}^{\text{sim}} = \hbar\delta\omega(t-t_j)a^\dagger a$（频率刺激，$a^\dagger a$ 形式，2 能级时与 $-\sigma_z/2$ 差常数项）
- **实验测量**：用 Virtual Z（相位平移）或真实磁通脉冲

当前代码**实际是"实验测量法"的磁通子版本**，但被命名为 `KernelEstimator` 没有体现这个选择。用户想做"纯理论仿真"必须手写一个新类。

### 1.4 高阶近似的物理动机

[`_sensing theory.md` "瞬态磁场协议" §6.1](../_sensing%20theory.md) 给出了一阶近似失效的量化判据 $\phi_{\rm acc}=|\delta\omega|_{\max}T_{\rm eff}\ll 1$ 及非线性响应的三种来源。这暗示：

- 当前 1 阶核函数 $k_1(t)$ 在 $\delta\omega \to 0$ 极限下严格成立
- $\delta\omega$ 不小时（大信号探测、强偏置工作点远离 sweet spot），$\Delta p_e$ 对 $\delta\omega$ 含 $O(\delta\omega^2), O(\delta\omega^3), \ldots$ 偏差
- 想用更大的刺激幅度提高信噪比 → 必须用高阶核函数补偿

### 1.5 总结：三个维度必须解耦

```
当前实现:  (flux, exp, order=1)  一个固定点
目标实现:  (flux|omega, sim|exp, order=1..N)  三维连续选择
```

---

## §2 三维设计空间

### 2.1 三个正交维度

```
维度 1: mode    ∈ {'flux', 'omega'}        刺激物理量
维度 2: method  ∈ {'sim', 'exp'}            模拟方式
维度 3: order   ∈ {1, 2, 3, ...}           Volterra 阶数
```

### 2.2 合法组合矩阵

| (mode, method) | qubit | order=1 | order≥2 |
|----------------|-------|---------|---------|
| (omega, sim) | 可选 | ✓ 笔记"仿真"法标准实现 | ✓ 对角 Volterra 提取 |
| (omega, exp) | 必需 | ✓ Virtual Z 实测对应 | ✓ Virtual Z 幅度扫描 |
| (flux, sim) | ❌ 非法 | — | — |
| (flux, exp) | 必需 | ✓ **当前默认行为** | ✓ flux 幅度扫描 |

非法组合 `(flux, sim)` 的物理含义：磁通刺激必须经过 $\omega(\Phi)$ dispersion → 必须有 qubit。如果想要"纯理论 flux kernel"，等价于 `(omega, sim)` 的结果乘以一个标量 $\kappa$——用户自己换算更清楚，不要在 API 里造一个伪法律组合。

### 2.3 默认值选择

```python
KernelEstimator(
    mode='flux',        # 默认 'flux'（向后兼容当前行为）
    method='exp',       # 默认 'exp'（向后兼容）
    order=1,            # 默认 1 阶（向后兼容）
)
```

**关键设计点**：所有默认值匹配当前实现，**已有用户代码完全不需改**。

---

## §3 理论推导

> 本节是 handbook 的核心。所有后续实现细节都从这里推出。

### 3.1 基本设定

旋转坐标系下，二能级 Transmon 哈密顿量（RWA 近似）：

$$H(t) = \frac{1}{2}\delta\omega(t)\,\sigma_z + H_{\text{ctrl}}(t) \tag{3.1}$$

其中 $H_{\text{ctrl}}(t)$ 是控制脉冲（已知），$\delta\omega(t)$ 是**待测频率信号**（小量微扰）。

频率信号和磁通信号的关系（Transmon 色散）：

$$\delta\omega(t) = \omega_q(\Phi_w + \delta\Phi(t)) - \omega_q(\Phi_w) = \kappa\,\delta\Phi(t) + \frac{1}{2}\kappa'\,\delta\Phi(t)^2 + \cdots \tag{3.2}$$

其中 $\kappa = d\omega/d\Phi|_{\Phi_w}$，$\kappa' = d^2\omega/d\Phi^2|_{\Phi_w}$。

### 3.2 Volterra 级数展开（核函数的严格定义）

测量结果 $\Delta p_e(t_d)$ 是输入 $\delta\omega(\cdot)$ 的**泛函**。对其做 Volterra 级数展开（关于 $\delta\omega = 0$ 点），采用与 [`_sensing theory.md` §6.2](../_sensing%20theory.md) 一致的约定——$1/n!$ 因子写在求和号前、$k_n$ 不含 $1/n!$：

$$\boxed{\Delta p_e(t_d) = \sum_{n=1}^{\infty} \frac{1}{n!} \int\!\!\cdots\!\!\int k_n^{(\omega)}(t_1-t_d, \ldots, t_n-t_d)\,\delta\omega(t_1)\cdots\delta\omega(t_n)\,dt_1\cdots dt_n} \tag{3.3}$$

其中 **n 阶核函数**定义为泛函导数（无额外系数）：

$$k_n^{(\omega)}(\tau_1, \ldots, \tau_n) = \left.\frac{\delta^n p_e}{\delta\omega(\tau_1)\cdots\delta\omega(\tau_n)}\right|_{\delta\omega=0} \tag{3.4}$$

#### 3.2.1 性质

1. **对称性**：$k_n$ 对参数置换完全对称（$\delta^n p_e / \delta\omega(\tau_1)\delta\omega(\tau_2) = \delta^n p_e / \delta\omega(\tau_2)\delta\omega(\tau_1)$，对所有置换成立）。
2. **因果性**：脉冲序列在 $[0, T]$ 期间作用，$k_n(\tau_1, \ldots, \tau_n) = 0$ 若任一 $\tau_i \notin [0, T]$。
3. **量纲**：$[k_n^{(\omega)}] = \text{rad}^{-n}$（无量纲）。

#### 3.2.2 一阶项（当前实现）

对应当前 KernelEstimator 提取的核函数：

$$\Delta p_e^{(1)}(t_d) = \int k_1^{(\omega)}(t-t_d)\,\delta\omega(t)\,dt = (k_1^{(\omega)} * \delta\omega)(t_d) \tag{3.5}$$

$n=1$ 时 $1/1! = 1$，一阶核与笔记 §1.2 的 $k_1 = i\langle 0|[W,Q]|0\rangle$ 完全一致。

#### 3.2.3 高阶项

二阶项含**两个时间变量**的非对角元，写入 $\tfrac12$ 因子：

$$\Delta p_e^{(2)}(t_d) = \frac{1}{2}\iint k_2^{(\omega)}(t_1-t_d, t_2-t_d)\,\delta\omega(t_1)\delta\omega(t_2)\,dt_1 dt_2 \tag{3.6}$$

完整 $k_2$ 是二维函数，存储成本 $O(N^2)$。$k_2$ 的算符表达见笔记 §6.3：$k_2(t_>,t_<) = -\langle 0|[W(t_<),[W(t_>),Q]]|0\rangle$。

#### 3.2.4 对角化近似（本 phase 采用）

当 $\delta\omega(t)$ 在**记忆时间** $\tau_{\text{mem}}$（脉冲序列时长）内基本恒定（准静态近似），可做对角化。此时 $[\delta\omega(t)]^n$ 在非对角元上的积分退化为同一点取值：

$$\Delta p_e^{(n)}(t_d) \approx \frac{1}{n!}\int k_n^{(\omega)}(t-t_d, t-t_d, \ldots, t-t_d)\,[\delta\omega(t)]^n\,dt \tag{3.7}$$

记**对角核** $k_n^{\text{diag}}(\tau) \equiv k_n^{(\omega)}(\tau, \tau, \ldots, \tau)$，则：

$$\boxed{\Delta p_e(t_d) \approx \sum_{n=1}^{N} \frac{1}{n!}\int k_n^{\text{diag}}(t-t_d)\,[\delta\omega(t)]^n\,dt} \tag{3.8}$$

**这是本 phase 实际实现的近似**。物理图像：$k_1^{\text{diag}}, k_2^{\text{diag}}, k_3^{\text{diag}}, \ldots$ 都是 1D 函数，分别捕获"线性响应"、"二阶非线性"、"三阶非线性"。

> **记号约定**：下文用 $k_n^{\text{diag}}$（或简写 $\tilde{k}_n$，等价于笔记的 $k_n(t,\ldots,t)$ 对角元）表示对角核，用 $k_n$（无上标）表示完整 $n$ 变量核。笔记 §6.6 的 $k_n^{(\omega)}$、$k_n^{(\Phi)}$ 指完整核；压缩到对角时需显式注明。

#### 3.2.5 对角化的适用条件

(3.8) 严格成立的条件：

$$\tau_{\text{signal}} \ll \tau_{\text{mem}} \quad \text{或} \quad \tau_{\text{signal}} \gg \tau_{\text{mem}} \tag{3.9}$$

即信号变化时间尺度和核函数记忆时间尺度**充分分离**。在两者可比时（典型瞬态测量），(3.8) 是 leading-order 近似，残差是 $k_2$ 的非对角元——这些项的贡献和 $|d\delta\omega/dt| \cdot \tau_{\text{mem}}$ 同阶。

对于 transmon sensing 的典型场景：脉冲时长 $\tau_p \sim 10$ ns，信号变化时间尺度若 $\geq 100$ ns，对角化是良好近似。

### 3.3 一阶核函数的解析形式（验证基准）

#### 3.3.1 闭式推导

从 (3.1) 出发，把 $\delta\omega(t)\sigma_z/2$ 当微扰。互作用表象下的演化算符到一阶：

$$U_I(T, 0) = I - \frac{i}{2}\int_0^T U_{\text{ctrl}}^\dagger(t', 0)\sigma_z U_{\text{ctrl}}(t', 0)\,\delta\omega(t')\,dt' + O(\delta\omega^2) \tag{3.10}$$

代入 $p_e = |\langle 1| U_{\text{ctrl}}(T,0) U_I(T,0) |0\rangle|^2$，对 $\delta\omega$ 取一阶变分：

$$\boxed{k_1^{(\omega)}(t) = \Im\Bigl[\langle 0|U_{\text{ctrl}}^\dagger(T,0)\,\sigma_z\,U_{\text{ctrl}}(T,t)\,\sigma_z\,U_{\text{ctrl}}(t,0)|0\rangle\Bigr]} \tag{3.11}$$

#### 3.3.2 理想 $\pi/2$-$\pi/2$ 方波脉冲

笔记 "瞬态磁场协议" §2(b) 给出闭式解（方波包络，$k_1(t) = -\tfrac12\sin\theta_y(t)$ 在第一个脉冲内等）：

$$k_1^{(\omega)}(t) = \begin{cases} \sin[\Omega(\tau_p/2 - |t|)] & |t| < \tau_p/2 \\ 0 & \text{otherwise} \end{cases} \tag{3.12}$$

**单测必须验证**：sim 模式 + 方波包络下数值结果与 (3.12) 在 $10^{-6}$ 量级吻合。

### 3.4 高阶核函数的推导

#### 3.4.1 二阶项

继续展开 (3.10) 到二阶：

$$U_I^{(2)}(T,0) = -\frac{1}{4}\int_0^T dt_1 \int_0^{t_1} dt_2\,U_{\text{ctrl}}^\dagger(t_1,0)\sigma_z U_{\text{ctrl}}(t_1, t_2)\sigma_z U_{\text{ctrl}}(t_2,0)\,\delta\omega(t_1)\delta\omega(t_2) \tag{3.13}$$

(3.13) 中的 $1/4$ 来自 $W=Z_c/2$（$W(t_1)W(t_2)=\tfrac14 Z_c(t_1)Z_c(t_2)$，再利用 $U_{\rm ctrl}(t_1,0)U_{\rm ctrl}^\dagger(t_2,0)=U_{\rm ctrl}(t_1,t_2)$ 化简）。完整的 $k_2^{(\omega)}$ 算符表达式见笔记 [`_sensing theory.md` §6.3](../_sensing%20theory.md)：$k_2(t_>,t_<)=-\langle 0|[W(t_<),[W(t_>),Q]]|0\rangle$。**实现层不直接用闭式**，而是用 §3.4.2 的数值幅度扫描法提取对角元。

#### 3.4.2 对角核 $k_n^{\text{diag}}$ 与刺激响应的关系

对一个**位于 $t_j$ 的窄高斯刺激**（宽度 $\sigma$，幅度 $\varepsilon$，面积 $\phi_\varepsilon = \varepsilon\sqrt{2\pi}\sigma$）：

$$\delta\omega(t) = \varepsilon\,e^{-(t-t_j)^2/(2\sigma^2)} \tag{3.14}$$

代入对角化 Volterra 展开 (3.8)：

$$\Delta p_e(t_j) = k_1^{\text{diag}}(t_j)\phi_\varepsilon + \frac{1}{2!}k_2^{\text{diag}}(t_j)\,\phi_\varepsilon^2\,c_2' + \frac{1}{3!}k_3^{\text{diag}}(t_j)\,\phi_\varepsilon^3\,c_3' + \cdots \tag{3.15}$$

其中 $c_n' = \int [e^{-x^2/(2\sigma^2)}]^n dx \,/\, (\sqrt{2\pi}\sigma)^n \cdot n!$ 是把 $[\delta\omega(t)]^n$ 的面积归一化到 $\phi_\varepsilon^n$ 时的积分修正因子（精确表达式见 §5.3）。

**关键观察**：$\Delta p_e$ 是 $\phi_\varepsilon$（或等价地 $\varepsilon$）的多项式。**在每个 $t_j$ 处扫描 $\varepsilon$，多项式拟合，提取系数 → 得对角核 $k_n^{\text{diag}}(t_j)$。**

### 3.5 单位与换算

#### 3.5.1 omega-based vs flux-based

omega 单位（频率刺激），对角近似下：

$$\Delta p_e = \sum_n \frac{1}{n!}\int k_n^{(\omega),\text{diag}}(t)\,[\delta\omega(t)]^n\,dt, \quad [k_n^{(\omega)}] = \text{rad}^{-n} = 1 \tag{3.16}$$

flux 单位（磁通刺激），对角近似下：

$$\Delta p_e = \sum_n \frac{1}{n!}\int k_n^{(\Phi),\text{diag}}(t)\,[\delta\Phi(t)]^n\,dt, \quad [k_n^{(\Phi)}] = \Phi_0^{-n} \cdot \text{ns}^{-(n-1)} \tag{3.17}$$

> 完整（非对角）形式带有 $1/n!$ 因子与多变量积分，见 (3.3)。本节后续的 (3.18)–(3.19) 对完整核与对角核均成立。

#### 3.5.2 线性色散下的换算

若 $\delta\omega = \kappa\,\delta\Phi$（色散是线性的）：

$$k_n^{(\Phi)}(t_1,\ldots,t_n) = \kappa^n\,k_n^{(\omega)}(t_1,\ldots,t_n) \tag{3.18}$$

#### 3.5.3 非线性色散下的混合

若 $\delta\omega = \kappa\delta\Phi + \frac{1}{2}\kappa'\delta\Phi^2 + \frac{1}{6}\kappa''\delta\Phi^3 + \cdots$，则 flux 核是色散非线性 + qubit 动力学非线性的**混合产物**。把色散展开代入 $\omega$ 核的 Volterra 级数，再按 $\delta\Phi$ 的幂次重新收集即得。系数已用符号微分验证，与笔记 [`_sensing theory.md` §6.6](../_sensing%20theory.md) 一致：

$$k_1^{(\Phi)}(t) = \kappa\,k_1^{(\omega)}(t) \tag{3.19a}$$

$$k_2^{(\Phi)}(t_1,t_2) = \kappa^2\,k_2^{(\omega)}(t_1,t_2) + \kappa'\,k_1^{(\omega)}(t_1)\,\delta(t_1-t_2) \tag{3.19b}$$

$$k_3^{(\Phi)}(t_1,t_2,t_3) = \kappa^3\,k_3^{(\omega)} + \kappa\kappa'\!\!\sum_{\text{3 对称项}}\!\!k_2^{(\omega)}\!\circ\!\delta + \kappa''\,k_1^{(\omega)}\!\circ\!\delta\!\circ\!\delta \tag{3.19c}$$

其中 (3.19c) 的显式展开为 $\kappa\kappa'[k_2^{(\omega)}(t_1,t_2)\delta(t_2-t_3) + k_2^{(\omega)}(t_2,t_3)\delta(t_1-t_3) + k_2^{(\omega)}(t_1,t_3)\delta(t_1-t_2)] + \kappa''k_1^{(\omega)}(t_1)\delta(t_1-t_2)\delta(t_1-t_3)$。

对角核（令所有 $t_i$ 相等）从 (3.19b) 直接读出：$k_2^{(\Phi),\text{diag}} = \kappa^2 k_2^{(\omega),\text{diag}} + \kappa' k_1^{(\omega)}\delta(0)$，其中 $\delta(0)$ 在离散实现中为 $1/dt$。

> **版本说明**：本文早期版本 (v2.0) 的换算系数为 $\tilde{k}_2^{(\Phi)} = \kappa^2\tilde{k}_2^{(\omega)} + \tfrac12\kappa'\tilde{k}_1^{(\omega)}$、$\tilde{k}_3^{(\Phi)} = \kappa^3\tilde{k}_3^{(\omega)} + \tfrac32\kappa\kappa'\tilde{k}_2^{(\omega)} + \tfrac16\kappa''\tilde{k}_1^{(\omega)}$。这些系数源自对角核记号下的换算错误——$\tfrac12\kappa'$ 应为 $\kappa'$（缺失的因子来自 $\tfrac12\kappa'\delta\Phi^2$ 中 $1/2$ 与 Volterra 对称化 $1/2!$ 的抵消不对等），$\tfrac32$ 与 $\tfrac16$ 对应地也需修正。本节 (3.19a–c) 采用经过符号微分验证的正确系数，与笔记 [`_sensing theory.md` §6.6](../_sensing%20theory.md) 的修正说明完全一致。代码中若出现上述旧系数，应立即替换。

#### 3.5.3.1 实现路径：方案 A（默认，本 phase 采用）

**用户决策**：本 phase 不通过 (3.19) 做后处理换算，而是**直接在 flux 域多项式扫描提取** $k_n^{(\Phi),\text{diag}}$。

理由：
1. **实验侧不可逆**——实验只能在 flux 域施加刺激（AWG 输出磁通脉冲）。要让 sim 和 exp 模式 API 一致，flux 高阶核必须能在 flux 域直接扫描。
2. **不需要 qubit 解析色散导数**——`TransmonQubit` 暴露的是 $\omega(\Phi)$ 数值方法，没有 $\kappa', \kappa''$ 的解析 API。引入 $\kappa', \kappa''$ 计算会触发 R10/R11（新公共 API、必须同步文档）。
3. **多项式扫描自然包含两类非线性**——在 $t_j$ 处扫描 $\varepsilon \in [-A, +A]$ 时，$\Delta p_e$ 拟合到 $\varepsilon^2$ 的系数**同时**包含 (3.19b) 右边两项；无需分离。

**约束代价**：flux 高阶核**不能拆解**成 qubit 动力学非线性 vs dispersion 非线性的独立贡献。若用户需要这种解耦，应改走 `mode='omega'` 路径（得到纯 $k_n^{(\omega)}$，与 dispersion 无关）。

方案 A 的实现伪代码：

```python
# 在 _extract_kn (flux 模式) 内部：
for t_j in t_samples:
    delta_p = []
    for eps in eps_grid:                # eps 是 flux 幅度 (Φ₀)
        stim = FluxSignal(type=3, amplitude=eps, center=t_j, width=sigma)
        H_total = build_total_hamiltonian(qubit, pulse, stim)
        p_e = mesolve(H_total, ...).expect[0][-1]
        delta_p.append(p_e - p_base)
    # phi_eps_grid = flux stimulus areas (Φ₀·ns)
    coeffs = polyfit_no_constant(phi_eps_grid, delta_p, N)
    # 直接得到 flux 对角核，不经过 (3.19) 换算
    for n in range(1, N+1):
        kernels[n-1][t_j] = coeffs[n-1] / c_n_flux(n, sigma)
```

**与方案 B/C 的区别**：方案 B/C 在 handbook 早期讨论中曾被列为备选——B 用 omega 核 + 后处理，C 显式传入 dispersion 导数。这两种路径**不在本 phase 实现**；用户需要这种能力时通过 `mode='omega'` 直接得到 $k_n^{(\omega),\text{diag}}$，可在 notebook 里自行按 (3.19) 后处理。

#### 3.5.4 物理解释

| 阶数 | $k_n^{(\omega)}$ 含义 | $k_n^{(\Phi)}$ 含义 |
|------|------------------------------|------------------------------|
| n=1 | qubit 对 $\delta\omega$ 的线性响应（**纯脉冲动力学**） | dispersion × qubit 线性响应 |
| n=2 | qubit 动力学的纯二阶非线性 | dispersion 二阶（$\kappa'$）+ dispersion×qubit 交叉项 |
| n=3 | qubit 动力学的纯三阶非线性 | dispersion 三阶（$\kappa''$）+ 多重交叉项 |

**结论**：omega 单位下的 $k_n^{(\omega)}$ **物理意义最清晰**——纯粹反映 qubit 动力学。flux 单位下的 $k_n^{(\Phi)}$ 是"工程方便量"，把所有非线性打包在一起，方便直接做 flux 信号反演。

### 3.6 sim 方法与 exp 方法的等价性

#### 3.6.1 sim 方法（笔记 "瞬态磁场协议" §3.3）

$H_{\text{stim}}^{\text{sim}} = \hbar\delta\omega(t-t_j) a^\dagger a$，**直接用频率单位**，不需要 qubit 的 $\omega(\Phi)$ 关系（但需要 $\alpha, n_{\text{levels}}$ 以含 leakage）。$a^\dagger a$ 在 2 能级时与 $-\sigma_z/2$ 差一个常数项，不影响差分结果。

#### 3.6.2 exp 方法 (a)：Virtual Z

$H_{\text{stim}}^{\text{vz}}$ 等价于在 $t_j$ 处**瞬时**对 qubit 施加 $\sigma_z\phi_z/2$ 的旋转。数学上等价于：

$$\delta\omega(t) = \phi_z \cdot \delta(t-t_j) \tag{3.20}$$

代入 (3.5) 直接得：

$$\Delta p_e(t_j) = \phi_z\,\tilde{k}_1^{(\omega)}(t_j) + O(\phi_z^2) \tag{3.21}$$

→ **Virtual Z 法和 sim 法对**一阶**核函数给出相同结果**（同 omega 单位）。

#### 3.6.3 exp 方法 (b)：磁通刺激（当前实现）

$H_{\text{stim}}^{\text{flux}}$ = 经过 `qubit_under_mag()` 注入的完整时变 Hamiltonian。在 Φ → ω 局部线性时等价于 sim 法乘 $\kappa$；非线性时按 (3.19) 修正。

#### 3.6.4 三种方法对 $k_1$ 的等价表

| 方法 | 单位 | 与 sim 的等价关系（线性色散区） |
|------|------|-------------------------------|
| sim（$\sigma_z$ 频率刺激） | $k_1^{(\omega)}$ | reference |
| exp Virtual Z | $k_1^{(\omega)}$ | 完全相同（含多能级 leakage） |
| exp 磁通刺激 | $k_1^{(\Phi)}$ | $k_1^{(\Phi)} = \kappa\,k_1^{(\omega)}$ |

**单测必须验证**：sim 模式 + 2 能级 + 方波包络下 $k_1^{(\omega)}$ 和理论解 (3.12) 一致；exp 磁通模式下 $k_1^{(\Phi)} / \kappa$ 与 sim $k_1^{(\omega)}$ 在 $10^{-3}$ 相对误差内一致（差异来自多能级修正）。

### 3.7 Volterra 反卷积

#### 3.7.1 单阶 Wiener（当前实现）

只用 $\tilde{k}_1$，FFT 域除法 + Tikhonov 正则：

$$\hat{X}(\omega) = \frac{\hat{K}_1^*(\omega)}{|\hat{K}_1(\omega)|^2 + \lambda^2}\hat{Y}(\omega) \tag{3.22}$$

→ 重建 $\delta\omega(t)$（omega 单位）或 $\Phi(t)$（flux 单位）。

#### 3.7.2 高阶 Hammerstein-Wiener 迭代

(3.8) 是 Hammerstein 结构（线性核函数 + 输入端多项式非线性）。用 fixed-point 迭代：

$$\delta\omega^{(0)} = \text{Wiener}_1(\Delta p_e) \tag{3.23a}$$

$$\delta\omega^{(k+1)} = \text{Wiener}_1\Bigl(\Delta p_e - \sum_{n=2}^{N}\,\tilde{k}_n * [\delta\omega^{(k)}]^n\Bigr) \tag{3.23b}$$

每次迭代用一阶 wiener 反卷积"残差"。收敛准则：$\|\delta\omega^{(k+1)} - \delta\omega^{(k)}\|_2 < \text{tol}$。

#### 3.7.3 适用条件

(3.23) 收敛要求 $|\sum_{n\geq 2} \tilde{k}_n * (\delta\omega)^n| \ll |\tilde{k}_1 * \delta\omega|$，即**高阶贡献是 leading-order 一阶项的 perturbation**。这正好是高阶核函数的有效区间——若高阶项已经 dominant，对角化近似 (3.8) 本身就失效，需要全 Volterra。

---

## §4 API 设计

### 4.1 KernelEstimator 类

```python
# sqc/reconstruction/kernel.py

@dataclass
class KernelEstimator:
    """Estimate the control kernel (Volterra series) of a pulse sequence.

    Three orthogonal design dimensions:
    - `mode`: physical quantity of the stimulus
    - `method`: simulation strategy
    - `order`: highest Volterra order extracted

    Parameters
    ----------
    mode : {'flux', 'omega'}
        'flux' → kernel acts on δΦ(t), units 1/(Φ₀ⁿ·ns).
        'omega' → kernel acts on δω(t), units rad⁻ⁿ (dimensionless).
    method : {'sim', 'exp'}
        'sim' → σ_z perturbation, no qubit dispersion (requires `kappa`
        only for flux mode).
        'exp' → full qubit Hamiltonian; uses `qubit.qubit_under_mag()`
        (flux mode) or Virtual Z (omega mode).
    order : int
        Highest Volterra order. Default 1 (linear, backward-compatible).
        order ≥ 2 enables polynomial-fit amplitude scan.

    Stimulus parameters
    -------------------
    stim_amplitude : float
        Reference amplitude. Units depend on mode (Φ₀ or rad·GHz).
    stim_width : float
        Gaussian σ (ns) for omega mode, or width parameter for flux mode.
    n_amp_samples : int
        Number of amplitudes scanned per t_j when order ≥ 2. Default 5.
    amp_scan_factor : float
        Amplitude scan range = [-factor, +factor] × stim_amplitude.
        Default 1.0.
    extract_off_diagonal : bool
        If True (default False), compute full off-diagonal k_n(t_i, t_j, ...)
        instead of diagonal-only tilde{k}_n(t). Stimulus count grows as
        O(M^n) where M = len(t_samples). Only enable when (a) signal
        timescale ≈ pulse timescale (breaking quasi-static approximation),
        OR (b) high-precision cross-time-correlation analysis required.
        See §3.4 for theoretical context.

    Sim-mode-only parameters
    ------------------------
    n_levels : int
        Hilbert space dimension when method='sim'. Default 2.
    anharmonicity : float
        Used when n_levels ≥ 3 and method='sim'. Default 0.
    kappa : float | None
        Linear dispersion d ω/dΦ (GHz/Φ₀). REQUIRED for sim mode + flux
        unit. Ignored otherwise.

    Exp-mode-only parameters
    ------------------------
    auto_calibrate : bool
        Adjust stim_amplitude to target Δp_e ~ 1% (avoids saturation).
    virtual_z_impl : {'math', 'hardware'}
        Virtual-Z stimulus implementation for omega-mode + exp-method:
        'math' (default) — insert narrow σ_z impulse via Hamiltonian list.
        'hardware' — rebuild CompositePulse with sub-pulse phase shift.
        See §5.2.2 for trade-offs.
    deprecation_warn_legacy : bool
        If True (default True), `pulse.get_kernel()` emits DeprecationWarning
        and forwards to KernelEstimator. Set False to silence in tests.

    Returns (from `estimate`)
    -------------------------
    KernelResult dataclass containing:
      - t_samples : ndarray of shape (M,)
      - kernels : list[ndarray], length = order
                  kernels[n-1] = diagonal n-th order kernel \tilde{k}_n(t)
      - metadata : dict with mode, method, units, etc.
    """
    mode: Literal['flux', 'omega'] = 'flux'
    method: Literal['sim', 'exp'] = 'exp'
    order: int = 1
    stim_amplitude: float = field(default_factory=lambda: CONFIG.reconstruction.stim_amplitude)
    stim_width: float = field(default_factory=lambda: CONFIG.reconstruction.stim_width)
    n_amp_samples: int = 5
    amp_scan_factor: float = 1.0
    extract_off_diagonal: bool = False
    n_levels: int = 2
    anharmonicity: float = 0.0
    kappa: float | None = None
    auto_calibrate: bool = False
    virtual_z_impl: Literal['math', 'hardware'] = 'math'
    deprecation_warn_legacy: bool = True

    def estimate(
        self,
        pulse,
        qubit=None,
        t_samples: np.ndarray | None = None,
    ) -> KernelResult:
        self._validate_inputs(qubit)
        if self.method == 'sim':
            return self._estimate_sim(pulse, qubit, t_samples)
        else:
            return self._estimate_exp(pulse, qubit, t_samples)
```

### 4.2 KernelResult dataclass

```python
@dataclass
class KernelResult:
    t_samples: np.ndarray
    kernels: list[np.ndarray]    # length == order; kernels[0] is k_1, etc.
    mode: Literal['flux', 'omega']
    method: Literal['sim', 'exp']
    order: int
    stim_amplitude: float
    units: str                    # e.g. '1/(Φ₀·ns)' or 'dimensionless'

    @property
    def k1(self) -> np.ndarray:
        """First-order kernel (backward-compat shortcut)."""
        return self.kernels[0]
```

### 4.3 合法性校验

`_validate_inputs` 实现下表：

| (mode, method, qubit, kappa) | 行为 |
|-----------------------------|------|
| (flux, sim, *, *) | **RAISE** `ValueError("flux+sim illegal; use omega+sim or supply kappa")` |
| (flux, exp, None, *) | **RAISE** `ValueError("flux+exp requires qubit")` |
| (omega, sim, None, _) | OK（纯理论） |
| (omega, sim, qubit, _) | OK（取 n_levels, anharmonicity 自动） |
| (omega, exp, None, *) | **RAISE** `ValueError("omega+exp requires qubit for Virtual Z")` |
| (omega, exp, qubit, *) | OK |

### 4.4 reconstruction 端配对

```python
# sqc/reconstruction/transient.py 扩展

@dataclass
class TransientReconstructor:
    method: Literal['wiener', 'hammerstein_volterra', 'lm']
    kernel_mode: Literal['flux', 'omega'] = 'flux'  # 与 KernelResult.mode 必须匹配
    max_volterra_iter: int = 5                       # only for 'hammerstein_volterra'
    volterra_tol: float = 1e-4

    def reconstruct(self, measurement, kernel: KernelResult, **kwargs):
        # method='wiener' → 只用 kernel.k1
        # method='hammerstein_volterra' → 用 kernel.kernels[:order]
        ...
```

### 4.5 工厂函数（便利接口）

```python
def make_kernel_estimator_for(use_case: str, **overrides) -> KernelEstimator:
    """Convenience factory matching common use cases.

    use_case ∈ {
      'waveform_reconstruction',   # flux + exp + order=1（默认）
      'frequency_calibration',     # omega + exp + order=1
      'theory_validation',         # omega + sim + order=1
      'large_signal_sensing',      # flux + exp + order=2
    }
    """
```

---

## §5 实现细节

### 5.1 sim 方法实现

#### 5.1.1 Hamiltonian 构造

```python
def _estimate_sim(self, pulse, qubit, t_samples):
    n = self.n_levels if qubit is None else qubit.n_levels
    a = destroy(n)
    
    # H_0 = anharmonicity only (no flux dependence)
    alpha = self.anharmonicity if qubit is None else qubit.anharmonicity
    if n == 2:
        H_0 = 0 * qeye(2)
    else:
        H_0 = (alpha / 2.0) * (a.dag() * a.dag() * a * a)
    
    # Pulse Hamiltonian (already in rotating frame)
    H_pulse = QobjEvo(pulse.hamiltonian, tlist=pulse.t_list, order=1)
    
    # Baseline
    psi_e = basis(n, 1)
    state_init = basis(n, 0)
    p_e_base = mesolve(
        H_0 + H_pulse, state_init, pulse.t_list, [],
        e_ops=[psi_e * psi_e.dag()],
    ).expect[0][-1]
    
    # Kernel extraction
    if self.order == 1:
        return self._extract_k1_sim(pulse, H_0, H_pulse, p_e_base, t_samples)
    else:
        return self._extract_kn_sim(pulse, H_0, H_pulse, p_e_base, t_samples)
```

#### 5.1.2 一阶提取（双边对称差分，消除偶次误差）

```python
def _extract_k1_sim(self, pulse, H_0, H_pulse, p_e_base, t_samples):
    # 单位转换：sim 模式刺激单位是 rad·GHz（频率）
    # 在 omega 单位下直接返回；在 flux 单位下用 kappa 换算
    
    kernel = np.zeros(len(t_samples))
    sigma = self.stim_width
    
    for i, t_i in enumerate(t_samples):
        # Symmetric +/- stimulus (eliminates O(ε²) error)
        stim_arr_plus = self.stim_amplitude * np.exp(
            -((pulse.t_list - t_i) ** 2) / (2 * sigma ** 2)
        )
        stim_arr_minus = -stim_arr_plus
        
        H_stim_plus = QobjEvo([sigmaz_n(n) / 2, stim_arr_plus],
                              tlist=pulse.t_list, order=1)
        H_stim_minus = QobjEvo([sigmaz_n(n) / 2, stim_arr_minus],
                               tlist=pulse.t_list, order=1)
        
        p_plus = mesolve(H_0 + H_pulse + H_stim_plus, ...).expect[0][-1]
        p_minus = mesolve(H_0 + H_pulse + H_stim_minus, ...).expect[0][-1]
        
        stim_area = self.stim_amplitude * sigma * np.sqrt(2 * np.pi)
        kernel[i] = (p_plus - p_minus) / (2 * stim_area)
    
    if self.mode == 'flux':
        kernel = kernel * self.kappa  # ω → Φ 换算
    
    return KernelResult(t_samples, [kernel], self.mode, 'sim', 1, ...)
```

注意 `sigmaz_n(n)` = $\text{diag}(0, 1, 2, ..., n-1) \cdot 2 - (n-1) = $ 等效的多能级 $\sigma_z$。在 n=2 时等于 `sigmaz()`，n≥3 时反映"频率刺激加在 $a^\dagger a$"上。具体形式：

$$\sigma_z^{(n)} = \sum_k (2k - (n-1)) |k\rangle\langle k| \quad \text{或者更标准地：用} \quad H_{\text{stim}} = \delta\omega \cdot a^\dagger a \tag{5.1}$$

**推荐用 $a^\dagger a$ 形式**（笔记 "瞬态磁场协议" §3.3）：在 n=2 时 $a^\dagger a = (I - \sigma_z)/2$，差一个常数项，对差分无影响。高阶提取 (§5.1.3) 统一使用 `a.dag() * a`。

#### 5.1.3 高阶提取（多项式拟合）

```python
def _extract_kn_sim(self, pulse, H_0, H_pulse, p_e_base, t_samples):
    """Polynomial-fit amplitude scan for diagonal k_1, ..., k_N."""
    N = self.order
    M = self.n_amp_samples
    
    # Symmetric amplitude grid
    eps_grid = np.linspace(
        -self.amp_scan_factor, self.amp_scan_factor, M
    ) * self.stim_amplitude
    
    kernels = [np.zeros(len(t_samples)) for _ in range(N)]
    sigma = self.stim_width
    
    for i, t_i in enumerate(t_samples):
        delta_p = np.zeros(M)
        
        for j, eps in enumerate(eps_grid):
            stim_arr = eps * np.exp(
                -((pulse.t_list - t_i) ** 2) / (2 * sigma ** 2)
            )
            H_stim = QobjEvo([a.dag() * a, stim_arr],
                             tlist=pulse.t_list, order=1)
            p_e = mesolve(H_0 + H_pulse + H_stim, ...).expect[0][-1]
            delta_p[j] = p_e - p_e_base
        
        # Fit polynomial: Δp_e = a_1·ε + a_2·ε² + ... + a_N·ε^N
        # (no constant term: a_0 = 0 by definition of Δp_e)
        coeffs = _fit_polynomial_no_constant(eps_grid, delta_p, N)
        
        # Convert coefficients to diagonal kernel values via §3.4.2 formula
        for n in range(1, N + 1):
            c_n = _gauss_integral_factor(n, sigma)  # see §5.3
            kernels[n-1][i] = coeffs[n-1] / c_n
    
    if self.mode == 'flux':
        # Apply §3.5.3 conversion (assumes linear dispersion for now)
        kernels = [k * (self.kappa ** (n+1)) for n, k in enumerate(kernels)]
    
    return KernelResult(t_samples, kernels, self.mode, 'sim', N, ...)
```

### 5.2 exp 方法实现

#### 5.2.1 flux 模式（保留当前行为）

基本沿用当前 [kernel.py:111-149](../../sqc/reconstruction/kernel.py#L111) 的逻辑：

- 用 `FluxSignal(type=3)` 构造刺激
- 用 `qubit.qubit_under_mag()` + `qubit_under_mag_hamiltonian()` 注入
- order=1 时单边差分（保 baseline 不变），order≥2 切换到对称扫描

**唯一改动**：order=1 时改成**双边差分**作为可选项（默认仍单边以保 baseline）。

#### 5.2.2 omega 模式 via Virtual Z

Virtual Z 的核心是在 $t_j$ 处对 qubit 注入相位刺激 $\phi_z$。两种实现给出**数学等价**的 $\Delta p_e$ vs $\phi_z$，因此核函数一致——但执行路径和数值代价不同。**两种实现都作为可选项**，由 `virtual_z_impl` 字段切换；默认走 `'math'`。

##### (a) `virtual_z_impl='hardware'`：硬件式相位重建

模拟真实硬件做 VZ 的方式：**在 $t_j$ 之后所有微波 sub-pulse 上叠加全局相位 $\phi_z$**。

```python
def _apply_virtual_z_hardware(composite_pulse, t_j, phi_z):
    """Rebuild CompositePulse with phase shift on sub-pulses starting after t_j."""
    new_subs = []
    for sub in composite_pulse.pulses:
        if sub.trigger >= t_j:
            new_subs.append(sub.with_phase_shift(phi_z))  # requires Tier A2 API
        elif sub.trigger + sub.duration > t_j:
            # t_j 落在 sub-pulse 内部 → 拒绝并提示用户切到 'math'
            raise ValueError(
                f"t_j={t_j} falls inside sub-pulse at trigger={sub.trigger}; "
                f"hardware VZ requires t_j on sub-pulse boundary. "
                f"Use virtual_z_impl='math' instead."
            )
        else:
            new_subs.append(sub)
    return CompositePulse(new_subs)
```

**优点**
- 与实验完全一致——同一段 control sequence 在仿真和实验里 indistinguishable
- 自动正确处理 sub-pulse 边界与多段串联

**缺点**
- $t_j$ 只能落在 sub-pulse 边界（否则抛错让用户切回 `'math'`）
- 每个 $t_j$ 要重建一次 CompositePulse → $O(M)$ 次额外对象构造
- **依赖 `Pulse.with_phase_shift()` API**（Tier A2 新增，仅 hardware VZ 需要）

##### (b) `virtual_z_impl='math'`（默认）：瞬时 σ_z 冲激

数学等价做法：在 Hamiltonian 里直接插入瞬时 $\sigma_z$ 冲激：

$$H_{\text{stim}}^{\text{vz}}(t) = \frac{\phi_z}{2}\,\delta(t - t_j)\,\sigma_z \tag{5.5}$$

数值实现用窄高斯近似 $\delta$（宽度 $\sigma_t \sim 2$ dt，面积归一化为 1）：

```python
def _apply_virtual_z_math(pulse, t_list, t_j, phi_z, n_levels, frame, omega_d):
    """Build σ_z impulse Hamiltonian term for Virtual Z."""
    sigma_t = 2.0 * CONFIG.awg.dt                  # 窄高斯宽度
    coeff = phi_z / (sigma_t * np.sqrt(2 * np.pi)) * np.exp(
        -0.5 * ((t_list - t_j) / sigma_t) ** 2
    )                                              # ∫ coeff dt ≈ φ_z
    # σ_z in n-level Fock basis: diag(1, -1, 0, ..., 0)
    sigma_z_op = qutip.Qobj(np.diag([1.0, -1.0] + [0.0] * (n_levels - 2)))
    return [0.5 * sigma_z_op, coeff]               # QuTiP list-format term
```

**优点**
- 与 pulse 结构解耦——单段、composite、空闲段都按同一路径处理
- $t_j$ 可取**任意时刻**，时间分辨率只受 `dt` 限制
- 实现简单，仅向 `H_list` 追加一项；不依赖 Pulse class 新 API

**缺点**
- 窄高斯近似 $\delta$ 引入 $O(\sigma_t)$ 量级近似误差（典型 < 0.1%，$\sigma_t = 2 \cdot dt$ 时）
- 不直接对应硬件——实验里 VZ 是 frame redefinition，物理上没有真正的 $\sigma_z$ 瞬时脉冲；解释 $t_j$ 落在 sub-pulse 内部的结果时要小心
- 高泄漏机制下 $\sigma_z$ 在 Fock 空间的截断需手动校验（上面用 $\text{diag}(1, -1, 0, ..., 0)$ 是 2 能级嵌入式定义）

##### 切换策略

| 场景 | 推荐 `virtual_z_impl` |
|------|----------------------|
| 高时间分辨率扫描，纯仿真验证 | `'math'`（默认） |
| 端到端与实验比对，控制序列已含硬件 VZ | `'hardware'` |
| `t_j` 落在 sub-pulse 内部 | 必须 `'math'`（`'hardware'` 抛错） |
| 高能级 leakage 严重的强驱动 | `'math'` + 校验 $\sigma_z$ 矩阵定义 |

##### 主提取流程

```python
def _extract_k1_exp_vz(self, pulse, qubit, t_samples):
    """Virtual Z phase modulation for omega-based kernel."""
    kernel = np.zeros(len(t_samples))
    phi_z = self.stim_amplitude  # 此处单位是 rad（相位增量）
    
    for i, t_i in enumerate(t_samples):
        if self.virtual_z_impl == 'hardware':
            pulse_p = _apply_virtual_z_hardware(pulse, t_i, +phi_z)
            pulse_m = _apply_virtual_z_hardware(pulse, t_i, -phi_z)
            p_plus = _simulate_pulse(pulse_p, qubit)
            p_minus = _simulate_pulse(pulse_m, qubit)
        else:  # 'math'
            H_vz_p = _apply_virtual_z_math(pulse, pulse.t_list, t_i, +phi_z,
                                            qubit.n_levels, pulse.frame, pulse.omega_d)
            H_vz_m = _apply_virtual_z_math(pulse, pulse.t_list, t_i, -phi_z,
                                            qubit.n_levels, pulse.frame, pulse.omega_d)
            p_plus = _simulate_pulse_with_extra(pulse, qubit, H_vz_p)
            p_minus = _simulate_pulse_with_extra(pulse, qubit, H_vz_m)
        
        kernel[i] = (p_plus - p_minus) / (2 * phi_z)
    
    return KernelResult(t_samples, [kernel], 'omega', 'exp', 1, ...)
```

#### 5.2.3 偶次误差消除（笔记 "瞬态磁场协议" §3.3 末尾）

可选 `symmetric_axis=True` flag：每次刺激跑两个旋转轴 $R_x$ 和 $R_{-x}$，差分：

$$p_e^{\text{sym}} = \frac{p_e|_{R_x} - p_e|_{R_{-x}}}{2} \tag{5.2}$$

减少"偶次脉冲不对称"导致的系统偏差。

### 5.3 数值积分因子

对角化 Volterra 展开 (3.8) 中，$[\delta\omega(t)]^n$ 对窄高斯刺激 $\delta\omega(t)=\varepsilon\,e^{-(t-t_j)^2/(2\sigma^2)}$ 的积分为

$$I_n \equiv \int_{-\infty}^{\infty} [\delta\omega(t)]^n\,dt = \varepsilon^n \int_{-\infty}^{\infty} e^{-n x^2/(2\sigma^2)}\,dx = \varepsilon^n\,\sigma\sqrt{\frac{2\pi}{n}} \tag{5.3}$$

若以刺激面积 $\phi_\varepsilon = \varepsilon\sqrt{2\pi}\,\sigma$ 为自变量，则 $\varepsilon = \phi_\varepsilon/(\sqrt{2\pi}\,\sigma)$，代入得

$$I_n = \phi_\varepsilon^n \cdot (2\pi)^{(1-n)/2}\,\sigma^{1-n}\,n^{-1/2} \tag{5.4}$$

代入 (3.8)，对角核 $k_n^{\text{diag}}$ 与测量响应 $\Delta p_e$ 的关系为

$$\Delta p_e(t_j) = \sum_{n=1}^{N} \frac{1}{n!}\,k_n^{\text{diag}}(t_j)\,I_n \tag{5.5}$$

**实现建议**：直接在代码中以 $\varepsilon$（而非 $\phi_\varepsilon$）为扫描变量，拟合 $\Delta p_e = \sum_n a_n \varepsilon^n$ 得系数 $a_n$，再由 $a_n = k_n^{\text{diag}} \cdot \sigma\sqrt{2\pi/n}\,/\,n!$ 反解 $k_n^{\text{diag}}$。这样避免了 $\phi_\varepsilon$ 表达式中 $(2\pi)^{(1-n)/2}\sigma^{1-n}$ 的繁琐幂次。

### 5.4 多项式拟合的数值稳定性

`_fit_polynomial_no_constant`：

```python
def _fit_polynomial_no_constant(x, y, N):
    """Fit y = a_1 x + a_2 x^2 + ... + a_N x^N (no constant term)."""
    # Use Vandermonde matrix V[i, n] = x[i]^(n+1) for n = 0..N-1
    V = np.vander(x, N=N+1)[:, :-1][:, ::-1]   # drop constant column
    # NOTE: condition number degrades for high N. Use orthogonal basis.
    if N >= 4:
        # Use Chebyshev or Hermite basis for stability
        coeffs = _fit_orthogonal_basis(x, y, N)
    else:
        coeffs, *_ = np.linalg.lstsq(V, y, rcond=None)
    return coeffs
```

**关键**：当 `order >= 4` 时 Vandermonde 矩阵病态，必须用正交基（Chebyshev / Hermite）拟合后再换回幂级数系数。本 phase 推荐 `order` 上限 3。

### 5.5 Legacy `pulse.get_kernel()` 兼容桥（方案 A）

为最小化对 [sqc/calibration/frequency.py](../../sqc/calibration/frequency.py)、[sqc/experiments/transient.py](../../sqc/experiments/transient.py) 等多处调用点的冲击，保留 `Pulse.get_kernel()` / `CompositePulse.get_kernel()` 作为 backward-compat shim：

```python
# sqc/control/pulse.py
import warnings

class Pulse:
    def get_kernel(
        self,
        qubit,
        t_samples=None,
        *,
        _suppress_deprecation: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """[DEPRECATED] Legacy kernel estimation interface.

        Internally forwards to KernelEstimator(mode='flux', method='exp',
        order=1). For new code, use KernelEstimator directly to access
        omega kernel, higher-order Volterra, sim mode, etc.
        """
        if not _suppress_deprecation:
            warnings.warn(
                "Pulse.get_kernel() is deprecated; use "
                "sqc.reconstruction.KernelEstimator(mode='flux', method='exp', "
                "order=1).estimate(pulse, qubit) instead. "
                "This shim will be removed in Phase 11.",
                DeprecationWarning,
                stacklevel=2,
            )
        from sqc.reconstruction.kernel import KernelEstimator
        est = KernelEstimator(mode='flux', method='exp', order=1,
                              deprecation_warn_legacy=False)
        result = est.estimate(self, qubit, t_samples=t_samples)
        return result.t_samples, result.k1  # 旧返回签名 (ndarray, ndarray)
```

**关键约定**：
- 旧返回签名 `(t_samples, kernel_array)` **必须保持**——`frequency.py:281` 等下游代码按位置解包
- 新返回的 `KernelResult` 通过 estimator 直接调用获取；shim 不暴露
- DeprecationWarning 用 `stacklevel=2` 让用户看到自己的调用点
- 单测里需要静默时显式传 `_suppress_deprecation=True`，避免污染测试输出

**Phase 10.5 迁移示例**（frequency.py 内部代码改动）：

```python
# 旧（保留作为兼容入口供其他下游调用，但 frequency.py 自己迁移）：
# t_samples, kernel_flux = ctrl_x.get_kernel(qubit)
# kappa = qubit.frequency_sensitivity(qubit.flux + flux)
# delta_omega = p_diff * kappa / G_diff      # ← κ workaround

# 新：
from sqc.reconstruction.kernel import KernelEstimator
result = KernelEstimator(mode='omega', method='exp', order=1).estimate(ctrl_x, qubit)
G_freq = np.trapezoid(result.k1, result.t_samples)
delta_omega = p_diff / G_freq                # 直接频率单位，无 κ
```

### 5.6 性能估计

| 配置 | 每个 $t_j$ 的 mesolve 次数 | 总 mesolve 次数（M=200 个 t_j） |
|------|----------------------------|--------------------------------|
| order=1（单边差分） | 1 | 200 |
| order=1（双边差分） | 2 | 400 |
| order=2（5 点扫描） | 5 | 1000 |
| order=3（7 点扫描） | 7 | 1400 |
| order=2 + `extract_off_diagonal=True` | M=200 | 200 × 200 = 4×10⁴ |
| `virtual_z_impl='hardware'` 额外开销 | +0 mesolve | +M 次 CompositePulse 重建（$O(M)$ Python，$\ll$ mesolve 耗时） |

对应当前一次 KernelEstimator.estimate() 在典型 transmon + 10 ns 脉冲 + 200 点时间轴下耗时约 30 s。order=3 约 7×30 = 210 s。可接受范围内。

`extract_off_diagonal=True` 时 order=2 已经达 ~6000 s（100 min），仅适合**高精度诊断**，不适合频繁调用。**不并行化**——Q7 已决策。

---

## §6 反卷积策略

### 6.1 一阶 Wiener（保持现有 [transient.py:_reconstruct_wiener](../../sqc/reconstruction/transient.py#L124)）

输入 `KernelResult` 时只用 `kernel.k1`，**忽略高阶项**。这是默认行为，保持现有 baseline 不变。

### 6.2 高阶 Hammerstein-Volterra

新增 `method='hammerstein_volterra'`：

```python
def _reconstruct_hammerstein_volterra(
    self, measurement, kernel: KernelResult, dt: float,
) -> FluxSignal:
    """Iterative Volterra inversion (3.23)."""
    Y = np.asarray(measurement.data['delta_p'], dtype=float)
    K = kernel.kernels  # list of k_1, k_2, ..., k_N
    N = len(K)
    
    # Initial guess: linear wiener
    X = _wiener_step(Y, K[0], dt, self.lambda_reg)
    
    for it in range(self.max_volterra_iter):
        # Build nonlinear correction
        correction = np.zeros_like(Y)
        for n in range(2, N + 1):
            correction += np.convolve(K[n-1], X ** n, mode='same') * dt
        
        # Wiener-deconvolve residual
        X_new = _wiener_step(Y - correction, K[0], dt, self.lambda_reg)
        
        if np.linalg.norm(X_new - X) / np.linalg.norm(X) < self.volterra_tol:
            X = X_new
            break
        X = X_new
    
    return FluxSignal(type=8, t_list=..., signal=X)
```

### 6.3 omega → Φ 反演（仅 omega kernel 路径）

当 `kernel.mode == 'omega'` 且用户最终需要 flux：

```python
def _omega_to_flux(self, delta_omega, qubit, branch='positive'):
    """Invert ω(Φ) → Φ. Handles sign ambiguity via continuity."""
    EC, EJ, freq = _get_qubit_params(qubit)
    omega_total = delta_omega + freq
    ratio = np.clip((omega_total + EC) ** 2 / (8 * EC * EJ), 0.0, 1.0)
    phi_abs = (1.0 / np.pi) * np.arccos(ratio)
    
    # Sign disambiguation:
    # - Reference: previous DC offset trap fix (commits 30a26d1, 3bc2a9c)
    # - Strategy: enforce continuity from initial condition + monotonicity
    return _resolve_sign(phi_abs, branch)
```

**这是 P10 必须解决的 DC offset trap 衍生问题**——与 P8 的 cryoscope 解决方案保持一致。

---

## §7 受影响文件清单

### Tier A1: 核心改动（必须）

| 文件 | 改动内容 | 风险 |
|------|---------|------|
| [sqc/reconstruction/kernel.py](../../sqc/reconstruction/kernel.py) | 重写 KernelEstimator 支持三维选择；新增 KernelResult；保持 `estimate(pulse, qubit)` 默认行为不变 | 中（broadly 修改单元，但向后兼容） |
| [sqc/config.py](../../sqc/config.py) ReconstructionConfig | 新增 `n_amp_samples`、`amp_scan_factor`、`max_volterra_iter`、`volterra_tol` 字段 | 低（仅添加字段） |
| [sqc/reconstruction/transient.py](../../sqc/reconstruction/transient.py) | 新增 `method='hammerstein_volterra'`；`reconstruct` 接受 `KernelResult` 或旧 ndarray | 中 |
| [sqc/reconstruction/__init__.py](../../sqc/reconstruction/__init__.py) | 导出 `KernelResult`、新工厂 `make_kernel_estimator_for` | 低 |

### Tier A2: 下游调用者适配

| 文件 | 改动 |
|------|------|
| [sqc/experiments/transient.py:118-123](../../sqc/experiments/transient.py#L118) | 调用方式不变（`pulse.get_kernel()` 走兼容路径），存入 `result.data["kernel"]` 仍是 ndarray；新增可选 `result.data["kernel_full"] = KernelResult` |
| [sqc/calibration/frequency.py:274-302](../../sqc/calibration/frequency.py#L274) | **重大优化**：用 `KernelEstimator(mode='omega', method='exp', order=1)` 替代当前 flux kernel + κ workaround；消除 L294-302 的 unit conversion 注释 |
| [sqc/calibration/waveform.py](../../sqc/calibration/waveform.py) | 间接受益，无直接改动 |
| [sqc/control/pulse.py](../../sqc/control/pulse.py) | **仅 `virtual_z_impl='hardware'` 需要**：新增 `Pulse.with_phase_shift(phi_z) -> Pulse` 方法（返回相位偏移后的新 Pulse 副本，不改原对象）；默认 `'math'` 路径不依赖此 API |
| [sqc/control/pulse.py](../../sqc/control/pulse.py) `pulse.get_kernel()` | 改为 backward-compat shim：内部转发到 `KernelEstimator(mode='flux', method='exp', order=1)`，emit `DeprecationWarning`（可通过 estimator 的 `deprecation_warn_legacy=False` 抑制） |
| Notebooks (`Simulation_sqc.ipynb`, `closed_loop_calibration_test.ipynb`) | 老调用方式继续工作；推荐迁移到新 API |

### Tier A3: 测试

| 文件 | 改动 |
|------|------|
| [tests/unit/test_kernel_estimator.py](../../tests/unit/test_kernel_estimator.py) | 扩展为 ~12 个 case：每种 (mode, method, order) 组合 + 合法性校验 |
| `tests/unit/test_kernel_volterra.py`（新增） | 高阶 Volterra 提取 + 反卷积的单测 |
| [tests/regression/test_physics_baseline.py](../../tests/regression/test_physics_baseline.py) | 添加 sim 模式 vs 解析公式 (3.12) 的对照 |
| `tests/equivalence/test_kernel_sim_vs_exp.py`（新增） | 小信号区 sim 和 exp 一致性测试 |

### Tier A4: 文档

| 文件 | 改动 |
|------|------|
| [docs/architecture.md](../../docs/architecture.md) | 新增 §X 核函数体系；按 R11 版本号 v2.x → v2.(x+1) |
| [idea/_sensing theory.md](../_sensing%20theory.md) | 在"核函数与卷积测量"节后加"高阶 Volterra 展开"小节；在"核函数的数值计算"小节明确 sim/exp 单位区别 |
| [.claude/agents/sqc-module-guide-knowledge.md](../../.claude/agents/sqc-module-guide-knowledge.md) | 更新 KernelEstimator 描述 |

---

## §8 测试与 baseline 策略

### 8.1 等价性测试（关键）

**T1: sim + omega + order=1 vs 解析公式 (3.12)**

```python
def test_sim_omega_order1_matches_analytic_pi2():
    """sim 模式应该精确再现 sin 形理论核函数。"""
    pulse = _ideal_pi2_pi2_square_pulse(Omega=0.5, tau_p=10.0)
    est = KernelEstimator(mode='omega', method='sim', order=1, n_levels=2)
    result = est.estimate(pulse)
    k_analytic = _analytic_sine_kernel(result.t_samples, Omega=0.5, tau_p=10.0)
    np.testing.assert_allclose(result.k1, k_analytic, rtol=1e-3, atol=1e-5)
```

**T2: sim vs exp 在线性区一致**

```python
def test_sim_vs_exp_linear_regime():
    """小刺激下 sim 和 exp 给出相同的 k_1（差 κ 因子）。"""
    qubit = _standard_qubit()
    pulse = _standard_pi2_pulse()
    
    k_sim = KernelEstimator(mode='omega', method='sim',
                             n_levels=qubit.n_levels,
                             anharmonicity=qubit.anharmonicity).estimate(pulse)
    k_exp = KernelEstimator(mode='omega', method='exp').estimate(pulse, qubit)
    
    np.testing.assert_allclose(k_sim.k1, k_exp.k1, rtol=1e-2)
```

**T3: flux vs omega（同 mode='exp'）的换算**

```python
def test_flux_omega_kappa_conversion():
    qubit = _standard_qubit()
    kappa = qubit.frequency_sensitivity(qubit.flux)
    
    k_flux = KernelEstimator(mode='flux', method='exp').estimate(pulse, qubit)
    k_omega = KernelEstimator(mode='omega', method='exp').estimate(pulse, qubit)
    
    # In linear dispersion regime: k_flux = κ * k_omega
    np.testing.assert_allclose(k_flux.k1, kappa * k_omega.k1, rtol=5e-2)
```

**T4: 高阶项的提取精度**

```python
def test_order2_extraction_polynomial_response():
    """在已知二次响应系统下，order=2 应该恢复 k_2。"""
    # Construct a synthetic system where Δp_e = a·ε + b·ε² is known
    ...
```

### 8.2 回归测试

- 现有 `test_transient_default_baseline` **保持不变**（默认 `mode='flux', method='exp', order=1` 与当前实现等价）
- 新增 `test_kernel_extension_baselines` 锁定新 sim/omega/高阶路径的数值

### 8.3 baseline 是否需要重做

| 测试 | 是否需要重做 baseline |
|------|----------------------|
| `test_qubit_static_properties` | 否 |
| `test_ramsey_default_baseline` | 否 |
| `test_diff_echo_default_baseline` | 否 |
| `test_transient_default_baseline` | **否**（默认行为不变） |
| `test_lm_default_baseline` | 否 |
| `test_predistortion_default_baseline` | 否 |

**结论**：本 phase 严格向后兼容，**不需要重做任何现有 baseline**。新功能的 baseline 全部以新文件形式添加。

### 8.4 性能基准

| 用例 | 当前 | 实施后（order=1 默认） | 实施后（order=2） |
|------|------|----------------------|-----------------|
| `test_transient_default_baseline` | 30 s | 30 s（无变化） | N/A |
| `test_kernel_omega_sim_order1` | — | ~15 s | — |
| `test_kernel_omega_sim_order2` | — | — | ~75 s |

**性能预算**：不允许默认路径 (order=1) 退化超过 5%。

---

## §9 分阶段实施 plan

### Phase 10.1 — 解耦 mode 维度（小到中，~250 行）

**目标**：把 `mode='flux'|'omega'` 拆出来，order=1, method='exp' 固定；同时实现两种 Virtual Z（`'math'` 默认 + `'hardware'` 可选）。

- 在 KernelEstimator 加 `mode` 字段，默认 `'flux'`
- 加 `virtual_z_impl` 字段，默认 `'math'`
- `mode='omega'` + `virtual_z_impl='math'`：实现 `_apply_virtual_z_math`（QuTiP list-format σ_z 冲激）
- `mode='omega'` + `virtual_z_impl='hardware'`：实现 `_apply_virtual_z_hardware` + `Pulse.with_phase_shift()` 新 API（Tier A2）
- `mode='flux'` 完全复用当前代码
- 测试 T1（分立的 omega kernel 正确性）、T3（κ 换算）、新增 T_VZ（两种 VZ 在边界点 $t_j$ 给出 ≤ 0.1% 偏差）

**PR 大小**：~250 行代码 + ~5 个单测。**baseline 无变化**。

### Phase 10.2 — 解耦 method 维度（中，~250 行）

**目标**：把 `method='sim'|'exp'` 拆出来。

- 在 KernelEstimator 加 `method` 字段，默认 `'exp'`
- 实现 `_estimate_sim`（直接 σ_z 刺激，不调 `qubit_under_mag()`）
- 实现 sim 模式的合法性校验
- 测试 T1（sim + omega 与解析公式对照）、T2（sim vs exp 一致性）

**PR 大小**：~250 行代码 + ~4 个单测。**baseline 无变化**。

### Phase 10.3 — 高阶提取（中，~250 行）

**目标**：实现 `order >= 2` 的对角 Volterra 提取 + 完整非对角可选 + KernelResult 序列化。

- 多项式拟合 `_fit_polynomial_no_constant`
- 对角元提取主循环 `_extract_kn_sim`、`_extract_kn_exp`（方案 A，直接 flux 域扫描）
- `extract_off_diagonal=True` 路径：$n$ 维网格 + n 元多项式拟合；下游只能配 `method='lm'`
- KernelResult 数据结构最终化 + `save(path)` / `load(path)`（`numpy.savez`）
- 测试 T4（已知二次响应系统）、T_off_diag（非对角恢复双 $\delta\omega$ 脉冲的 $k_2(t_1, t_2)$）

**PR 大小**：~250 行代码 + ~4 个单测。

### Phase 10.4 — 反卷积升级（小到中，~150 行）

**目标**：实现 `method='hammerstein_volterra'` 反卷积。

- `_reconstruct_hammerstein_volterra` 函数
- `_omega_to_flux` 含符号歧义处理
- 端到端测试：先施加已知 $\delta\omega$，正向计算 $\Delta p_e$，反卷积恢复 $\delta\omega$，比对

**PR 大小**：~150 行代码 + ~2 个端到端测试。

### Phase 10.5 — 下游迁移（小，~50 行）

**目标**：在 [sqc/calibration/frequency.py](../../sqc/calibration/frequency.py) 切换到 omega kernel，消除 L294-302 的 κ workaround。

**PR 大小**：~50 行 + 频率标定专项测试。

### Phase 10.6 — 文档（小）

- 更新 docs/architecture.md、_sensing theory.md、agent knowledge file
- 写入 R11 变更记录

---

## §10 Open questions / 后续 TODO

### Q1: 完整非对角 Volterra 是否实现 — **RESOLVED：作为可选项**

当前 phase 默认只做对角元（`extract_off_diagonal=False`）。完整 $k_2(t_1, t_2), k_3(t_1, t_2, t_3), \ldots$ 在 §4.1 中通过 `extract_off_diagonal=True` 启用。

启用条件（用户文档须明确列出）：
- 信号变化时间尺度 ~ 脉冲时长（破坏 (3.9) 的对角化近似）
- 需要研究 qubit 动力学的 cross-time correlation
- 高精度/高时间分辨率诊断场景

代价：
- 刺激次数 $O(M^n)$，$n$ 阶非对角 + $M=200$ 时间点 → $n=2$ 已达 $4\times10^4$ 次 mesolve；$n=3$ 不实际
- 存储 $O(M^n)$ 浮点数

**实施层约定**：`extract_off_diagonal=True` 时 KernelResult 的 `kernels[n-1]` 变为 n-维 ndarray（shape = `(M,) * n`）。下游 `TransientReconstructor` 在 `extract_off_diagonal=True` 时只能用 `method='lm'`（LM 优化器原生接受 n-D 核函数）；Wiener 和 Hammerstein-Volterra 路径不支持非对角，加 `ValueError` 提示。

### Q2: 非线性 dispersion 在 flux 高阶里的处理 — **RESOLVED：方案 A**

§3.5.3.1 已确定：本 phase 用**方案 A**——直接在 flux 域多项式扫描提取 $\tilde{k}_n^{(\Phi)}$，**不经 (3.19) 换算**。

要点回顾：
- 不需要 qubit 解析 $\kappa', \kappa''$；只用数值 mesolve 输出
- $\tilde{k}_n^{(\Phi)}$ 同时打包"qubit 动力学非线性"和"色散非线性"两类贡献，无法解耦
- 需要解耦的用户走 `mode='omega'` 路径，自行按 (3.19) 后处理（notebook level，不进 sqc/）

**实施层约定**：`KernelEstimator(mode='flux', order=N)` 直接对 `eps_grid`（flux 单位）做 polyfit；`metadata` 字段写入 `"flux_extraction_method": "direct_polyfit_diag_only"`。

### Q3: Virtual Z 实现方式 (a) 还是 (b) — **RESOLVED：两者均作为可选项**

§5.2.2 已展开：
- `virtual_z_impl='math'`（默认）—— 瞬时 $\sigma_z$ 冲激，与 pulse 结构解耦
- `virtual_z_impl='hardware'`（可选）—— 重建 CompositePulse 的 sub-pulse 相位，与硬件 1:1 对应

切换策略与限制见 §5.2.2 表格。`'hardware'` 依赖 `Pulse.with_phase_shift()` 新 API（Tier A2 实现）；`'math'` 不需要任何 Pulse class 新方法。

**实施层约定**：Phase 10.1 同时实现两种路径；Phase 10.1 的单测必须验证 `'math'` 和 `'hardware'`（在边界点 $t_j$）给出 $\leq 0.1\%$ 相对偏差。

### Q4: KernelResult 的序列化 + legacy `pulse.get_kernel()` 迁移 — **RESOLVED：方案 A，不并行**

**序列化**：Phase 10.3 引入 `KernelResult.save(path)` / `load(path)`，用 `numpy.savez`（不需要 pickle，避免类版本兼容性问题）。metadata dict 序列化为 `.npz` 内 0-d object array。

**legacy 接口迁移（Plan A）**：
- 保留 `pulse.get_kernel(qubit, ...)` 作为 backward-compat shim（被 [sqc/calibration/frequency.py:281](../../sqc/calibration/frequency.py#L281)、[sqc/experiments/transient.py:118-123](../../sqc/experiments/transient.py#L118) 等多处调用）
- 内部转发到 `KernelEstimator(mode='flux', method='exp', order=1).estimate(pulse, qubit)`
- 加 `DeprecationWarning`，指向新 API；用户可通过 `KernelEstimator(deprecation_warn_legacy=False)` 抑制（仅供测试用）
- **不引入并行化**——`n_jobs` 字段从设计中删除；性能优化留待 future phase 若 profiling 显示是瓶颈再考虑

**实施层约定**：Phase 10.5 把 [sqc/calibration/frequency.py](../../sqc/calibration/frequency.py) 切换到新 API（消除 κ workaround），同时让 [sqc/experiments/transient.py](../../sqc/experiments/transient.py) 仍走兼容路径直到 Phase 11 统一迁移。

### Q5: 与 `sqc/calibration/frequency.py` 的具体迁移路径 — **RESOLVED：与 Q4 合并**

已并入 Q4 决策：Plan A（保留兼容入口，deprecation warning），Phase 10.5 让 [sqc/calibration/frequency.py:274-302](../../sqc/calibration/frequency.py#L274) 切换到 `KernelEstimator(mode='omega', method='exp', order=1)`，消除 κ workaround。

### Q6: ExperimentResult.data["kernel"] 的兼容

当前 `transient.py` 把 kernel 存为 ndarray。新 API 需要 `KernelResult`。**建议**：

- `result.data["kernel"]` 保持是 ndarray（k_1 数组），向后兼容
- 新增 `result.data["kernel_full"]` 存 KernelResult dataclass

### Q7: 性能优化预留 — **RESOLVED：暂不并行化**

**用户决策**：本 phase 不引入并行化（`multiprocessing.Pool` 或 GPU）。

理由：
- order=3 + 7 点扫描 + 200 时间点 ≈ 210 s（§5.5），可接受
- multiprocessing 在 Windows 下 pickle qubit/QobjEvo 有兼容性风险
- 真正瓶颈应先用 profiling 定位（可能是 mesolve 内部而非 Python 循环开销）

后续若需要加速，优先考虑：
1. **采样优化**：自适应 `n_amp_samples`（动态精度）
2. **批量 mesolve**：QuTiP 5.x 的 `parallel_map` 后端（若稳定）
3. **粗粒度并行**：notebook 层面分 t_samples 块独立跑（用户手动）

GPU 加速复杂度过大，**不在本 phase 或近期 phase 考虑**。

---

## 附录 A: 三种方法的实施伪代码对比

### A.1 sim + omega + order=1

```
For each t_j:
  +stim = ε·gaussian(t-t_j, σ)
  +H_stim = QobjEvo([a†a, +stim])
  +p_+ = mesolve(H_0 + H_pulse + H_stim, ...)
  -stim = -ε·gaussian(t-t_j, σ)
  -H_stim = QobjEvo([a†a, -stim])
  +p_- = mesolve(H_0 + H_pulse + H_stim_neg, ...)
  k_1(t_j) = (p_+ - p_-) / (2 ε √(2π) σ)
```

### A.2 exp + omega + order=1（Virtual Z）

```
For each t_j:
  pulse_+ = apply_virtual_z(pulse, t_j, +φ_z)
  pulse_- = apply_virtual_z(pulse, t_j, -φ_z)
  p_+ = mesolve(H_0 + pulse_+.H, ...)
  p_- = mesolve(H_0 + pulse_-.H, ...)
  k_1(t_j) = (p_+ - p_-) / (2 φ_z)
```

### A.3 exp + flux + order=1（当前实现）

```
For each t_j:
  stim = FluxSignal(type=3, amp=A, center=t_j, width=σ)
  qubit_t = qubit.qubit_under_mag(stim)
  H_stim = QobjEvo(qubit.qubit_under_mag_hamiltonian(qubit_t, ...))
  p_+ = mesolve(H_0 + H_pulse + H_stim, ...)
  k_1(t_j) = (p_+ - p_base) / (A · σ · √(2π))   # 单边
```

### A.4 sim + omega + order=2

```
For each t_j:
  For ε in [-A, -A/2, 0, A/2, A]:
    stim = ε·gaussian(t-t_j, σ)
    H_stim = QobjEvo([a†a, stim])
    p_e[ε] = mesolve(H_0 + H_pulse + H_stim, ...)
    Δp_e[ε] = p_e[ε] - p_base
  
  # Fit Δp_e = a_1·ε + a_2·ε²  (fit in ε, not φ_ε — see §5.3)
  a_1, a_2 = polyfit_no_constant(ε_grid, Δp_e, order=2)
  
  # Convert: a_n = (1/n!)·k_n^diag · σ·√(2π/n)  (from §5.3 eq. 5.5)
  k_1(t_j) = a_1 / (σ·√(2π))          # n=1: I₁ = ε·σ·√(2π)
  k_2(t_j) = 2·a_2 / (σ·√(π))         # n=2: I₂ = ε²·σ·√(π), factor 1/2!
```

---

## 附录 B: 词汇表

| 术语 | 含义 |
|------|------|
| Volterra 级数 | 非线性输入-输出系统的泛函展开 |
| 对角核 $\tilde{k}_n$ | $k_n(\tau, \tau, \ldots, \tau)$，多变量核函数的对角元 |
| 准静态近似 | 信号变化时间尺度远大于核函数记忆时间 |
| Hammerstein 模型 | 输入端非线性 + 线性核函数 |
| Wiener 模型 | 线性核函数 + 输出端非线性 |
| Virtual Z | AWG 相位平移实现的瞬时 σ_z 旋转 |
| sim 模式 | 直接 σ_z 频率刺激，不依赖 qubit 色散 |
| exp 模式 | 通过 qubit 完整 Hamiltonian 注入刺激，含 leakage 等 |

---

## §11 实施承诺

实施本 phase 时严格遵守 [CLAUDE.md](../../CLAUDE.md) 的 R1–R11：

- **R1**: 不动 `src/`
- **R3**: 所有 baseline 保持数值不变（默认路径与当前实现 byte-equivalent）
- **R8**: 实施前 commit safety checkpoint；分 Phase 10.1～10.6 逐步推进
- **R9**: 时间轴仍从 `CONFIG.awg.dt` 派生
- **R10**: 不引入新依赖
- **R11**: 完成后更新 docs/architecture.md，版本号 v2.x → v2.(x+1)

实施前必须读：
1. [_refactor_plan.md](_refactor_plan.md) §0、§8、§13、§15
2. [_handoff_state.md](_handoff_state.md)
3. 本 handbook 完整 §3 理论部分
