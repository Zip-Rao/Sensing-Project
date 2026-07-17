# 核函数与测频精度 · 统一理论笔记

> 日期: 2026-06-07 | 这是**核函数体系 + 瞬态测频精度**的唯一权威理论笔记。
> 代码行级对应: [transient_frequency_theory.md](refactor/transient_frequency_theory.md)(frequency.py 逐行) ·
> [kernel_and_frequency_architecture.md](refactor/kernel_and_frequency_architecture.md)(API/架构)
> 分析记录(证据): [高阶核函数差异根因.md](../result/transient_error/高阶核函数差异根因.md) ·
> [瞬态测频误差来源分析.md](../result/transient_error/瞬态测频误差来源分析.md)
> 代码: `sqc/reconstruction/kernel.py`, `sqc/calibration/frequency.py`
> 验证/图: [frequency_method_comparison.ipynb](../kernel/frequency_method_comparison.ipynb) ·
> [verify_full_kernel.py](../kernel/verify_full_kernel.py) · 流程图 `result/flowchart/kernel_multidim.tex`, `QSL_highorder.tex`
> 基准参数: 2 能级, EC=0.2·2π, EJ=15·2π, π/2 脉冲 10 ns, dt=0.5 ns

---

## 0. 全局图景

```
                         ┌─ sim (Heisenberg 对易子, 精确, 需 H(t))
核函数 k_n(t₁..tₙ) ──────┤
   (Volterra 核)         └─ exp (Virtual-Z + 有限差分, 可测量)
        │
        │ 积分
        ▼
   标量 G_α = ∫…∫ k_α        ──► 瞬态测频:  p_diff = G₁Δ + (1/6)G₃Δ³ + …
                                          反演 Δ → f = ω_d − Δ
                                          翻折点 Δ_fold 封顶 → 大失谐改用 Ramsey
```

一句话:**核函数描述 qubit 对失谐的 n 阶响应;把它积分成标量 G_α 喂给瞬态测频;测频精度由 G_α 的算法精度 + 泰勒翻折点共同决定。**

---

## 1. 核函数理论

### 1.1 Volterra 定义

末态 |e⟩ 布居对失谐微扰 δω(t) 的泛函展开:

$$p_e = p_0 + \sum_n \int\!\cdots\!\int k_n(t_1,\dots,t_n)\,\delta\omega(t_1)\cdots\delta\omega(t_n)\,d^n t$$

$k_n$ = n 阶 Volterra 核 = $\delta^n p_e/\delta\omega^n$。一阶 $k_1(t)$ 是线性响应核。

### 1.2 阶数 · 对角 vs 非对角

- **对角元** $k_n^{\rm diag}(t)=k_n(t,t,\dots,t)$:同一时刻的 n 阶响应,1 维数组。
- **非对角(完整)** $k_n(t_1,\dots,t_n)$:不同时刻,n 维对称张量(shape `(M,)*n`)。
- 用途分流:**逐时核**(k(t))→ Wiener/LM 重建时变场;**积分标量** $G_α=\int…\int k_α$ → 测频(对恒定失谐)。

### 1.3 两条计算路线

| | **sim**(Heisenberg) | **exp**(测量式) |
|---|---|---|
| 原理 | 一次 `sesolve`→U(t),嵌套对易子 | Virtual-Z 探针 + 有限差分 |
| 公式 | $k_n=i^n\langle0\|[W(t_n),[\cdots,[W(t_1),Q]]]\|0\rangle$ \| $\partial^n p_e/\partial\phi_z^n$,5 点 FD stencil |
|  | $W=U^\dagger\sigma_z U/2,\;Q=U^\dagger(T)MU(T)$ | |
| 精度 | 机器精度 | $O(h^4)/O(h^2)$ + 有限-σ_t 偏差 |
| 非对角 | ✅ 一次 sesolve 出整个张量(`_heisenberg_kernels_offdiag`, order≤3) | ✅ Phase 14: 混合 FD (`_extract_kn_offdiag_exp`, order≤3, omega + flux; 代价 ~O(M³) mesolve, 建议 M≤8) |
| 适用 | 仅仿真(需知 H(t)) | 任何可测 p_e(含真实硬件) |

代码:`KernelEstimator(mode∈{flux,omega}, method∈{sim,exp}, order, extract_off_diagonal)`。

### 1.4 关键性质

- **零失谐 Y-X Ramsey:$k_2^{\rm diag}\equiv0$**(脉冲 Y-X 正交对称严格保证,非近似);$k_3^{\rm diag}=-k_1$(2 能级对易子周期-2)。Gaussian 包络有限上升时间产生 ~4% 残余 k₂。
- 历史 bug:多项式拟合提取高阶核病态(Vandermonde),已被 5 点 FD stencil 取代(见 [phase_10_kernel_theory_and_fix.md](refactor/phase_10_kernel_theory_and_fix.md))。

---

## 2. 核函数估计的精度

### 2.1 exp 高阶对角偏差(两个根因)

1. **σ_t 涂抹**:VZ 高斯探针宽度 σ_t 有限,把完整非对角核在 σ_t 球内卷积平均;非对角"体积"巨大(实测 $G_3^{\rm full}\approx190\times G_3^{\rm diag}$)且沿时序楔形陡峭 → 对角估计偏低。默认 σ_t=2·dt 时 **k₃ 偏低 ~15%**(k₁ 仅 0.7%)。
2. **FD 容差**:stencil 除以 $h^n$,默认 mesolve 容差 ~1e-8,使 noise$/h^3\sim$1e-2/点主导高阶。

### 2.2 σ_t 旋钮 + Richardson 外推(Phase 12 修复)

- FD mesolve 收紧到 `atol=1e-12, rtol=1e-10`(消除噪声主导)。
- 新增 `probe_sigma_t`(默认 None→2·dt,零回归);`richardson=True` 对多个 σ_t 外推 σ_t→0。
- 效果(Y-X):$G_3/G_3^{\rm sim}$ 从 **0.844(2·dt)→0.922(dt)→0.961(Richardson)**。σ_t<dt 网格欠采样会崩,故 σ_t≥dt。

### 2.3 非对角(sim + exp)提取

- **sim**:`extract_off_diagonal=True` 经 `_heisenberg_kernels_offdiag` 一次 sesolve + numpy 对易子产出 $k_2(M,M)$、$k_3(M,M,M)$(order≤3)。验证:对角切片 $k_3(t,t,t)=-k_1$(RMSE~1e-6),$k_2$ 对称,$G_2^{\rm full}\approx0$。
- **exp**(Phase 14):`_extract_kn_offdiag_exp` 混合偏导数 FD → 同 shape 的 n 维张量,从可测 $p_e$ 恢复,无需知 $H(t)$。支持 omega + flux 双模式;含 memoized 求值缓存、对称楔形填充。代价 ~O(M³) mesolve(M=时间点数),建议 M≤8。
- `KernelResult.off_diagonal` 标记 + n 维 save/load;Wiener/Hammerstein 收到 n 维核抛 ValueError(仅 LM 可消费)。

### 2.4 G_freq 路线:为何一阶测频比旧结果更准

一阶测频精度全看 $G_{\rm freq}=\int k_1$ 算多准。P10.5 把它从 flux 路线换成 omega 路线:

| 路线 | 误差来源 | 实测相对误差(bias=0.9) |
|---|---|---|
| **flux 高斯核 / κ**(旧) | 高斯近似δ + κ 局部线性化 + 单边差分(三重叠加) | **186%**(符号都翻) |
| **omega 虚拟Z核**(现) | 双边差分、免 κ、模拟理想 δ | **2.6%** |

→ 一阶比旧 output 准的根因是核算法升级,不是反演公式变化。

---

## 3. 瞬态测频管线

### 3.1 线性反演

τ=0 正交 Ramsey(R_y–R_{±x}),差分布居:

$$p_{\rm diff}=\frac{p_x-p_{-x}}{2}=G_{\rm freq}\,\Delta+O(\Delta^3),\qquad
\boxed{f_{\rm meas}=\omega_d-\frac{p_{\rm diff}}{G_{\rm freq}}}$$

$p_{\rm diff}$ 对 Δ 是**奇函数**(偶次项差分抵消)→ 一阶精确到 $O(\Delta^3)$。$G_{\rm freq}=\int (k_1^x-k_1^{-x})/2\,dt$。

### 3.2 符号约定(关键)

omega 虚拟Z 核以 $+\phi_z=-\Delta$ 约定 → $G_{\rm freq}>0$,$\delta\omega=p_{\rm diff}/G_{\rm freq}=-\Delta$,故 $f=\omega_d-\delta\omega=\omega_d+\Delta$ 正确。**所有高阶修正必须在同一 $\delta\omega=-\Delta$ 约定下**(见 §4.3)。

---

## 4. 高阶测频

### 4.1 ⚠️ G₃ 是 Taylor 系数(三重积分),不是对角核积分

对**全程恒定** Δ 的反演:

$$G_3^{\rm Taylor}=\frac{d^3p_{\rm diff}}{d\Delta^3}\Big|_0=\iiint k_3(t_1,t_2,t_3)\,d^3t\;\neq\;\int k_3^{\rm diag}\,dt$$

二者差 $\sim(2T_{\pi/2})^2$(实测 $G_3^{\rm diag}$ 小约 **170×**、符号也不同)。**用对角核做三次修正几乎不生效**——历史 bug,根因见 [高阶核函数差异根因.md](../result/transient_error/高阶核函数差异根因.md)。

### 4.2 两条产 G₃ 的正确路线

| `g3_source` | 算法 | 代码 |
|---|---|---|
| `"fit"`(默认) | 奇多项式拟合 $p_{\rm diff}(\Delta)$,**自适应 delta_max**(收缩扫描至 G₁ 收敛) | `_calibrate_g3_taylor` |
| `"kernel_full"` | $\iiint k_3$(Heisenberg sim,免扫描) | `_calibrate_g3_kernel_full` |

二者互为交叉验证(实测 |G₃| 吻合 ~10%)。**对角核 `diag_legacy` 已移除**。
> **delta_max 坑**:旧默认 0.08 GHz 超出线性区 → 拟合被条纹饱和带偏,$G_1$ 偏低 0.6×、$G_3$ 全错。`None`=自适应修复;`FrequencyMeasurement.g3_delta_max` 可手动覆盖。

### 4.3 符号 bug 修复(Phase 12.1)

旧版三次 Newton 用 $G_1^{\rm fit}$(Δ 约定,负)与线性 $G_{\rm freq}$(δω 约定,正)**符号相反** → 返回 $\omega_d-\Delta$(误差 ≈ $-2\Delta$,比线性更差)。现统一 δω=-Δ:Route A 取 $(-G_1^{\rm fit},-G_3^{\rm Taylor})$,Route B 用原始核积分。Newton 解 $p_{\rm diff}=G_{\rm lin}\delta\omega+\tfrac16 G_3\delta\omega^3$,$f=\omega_d-\delta\omega$。

### 4.4 为何高阶扩展线性区

- 线性法系统误差 = 被忽略的三次项 $\propto\Delta^3$;三次法吃掉它,残余 $\propto\Delta^5$。
- $\Delta\ll1$ 时 $\Delta^5\ll\Delta^3$ → 同容差下可容忍更大 Δ → 安全区变宽。
- 安全区边界 $\Delta_n\sim{\rm tol}^{1/(n+2)}$:order-n 误差 $\propto\Delta^{n+2}$,每加一奇次项再外推一截。
- 实测(bias=0.9,有效区):order3-fit 比 order1 **~10×**(0.15 vs 1.39 MHz)。

---

## 5. 翻折点(fold)理论

### 5.1 起源

$p_{\rm diff}\simeq A\sin(\varphi)$,$\varphi=\Delta\cdot T_{\rm eff}$($T_{\rm eff}$=核加权脉冲时长,τ=0 时≈两 π/2 脉冲)。在 **φ=π/2 翻折**:极值后回卷、多值不可逆。实测极值 ≈±0.5,单调区 |Δ|≲0.02 GHz。

### 5.2 位置 · 调控 · 守恒

$$\boxed{\Delta_{\rm fold}\approx\frac{\pi}{2T_{\rm eff}}}\quad(\text{10 ns π/2}\approx0.02\text{ GHz})$$

缩短 π/2 脉冲推远翻折,$\Delta_{\rm fold}\propto1/T_{\pi/2}$,但**量程–灵敏度守恒**:

$$\boxed{G_1\cdot\Delta_{\rm fold}\approx\text{const}\;(\approx0.124)}$$

| π/2 时长 | Δ_fold | G₁ | G₁·Δ_fold |
|---|---|---|---|
| 10 ns | 0.020 GHz | −6.22 | 0.124 |
| 5 ns | 0.040 GHz | −3.07 | 0.123 |
| 2.5 ns | 0.088 GHz | −1.42 | 0.125 |

### 5.3 翻折点 = 泰勒收敛半径

$p_{\rm diff}=G_1\Delta+\tfrac16 G_3\Delta^3+\tfrac1{120}G_5\Delta^5+\cdots$ 的收敛半径就是 Δ_fold;**越过它任何有限阶失效**。
- 实用安全区:**|Δ|≲0.5·Δ_fold**(此处约 ≲10 MHz)。
- 三次截断模型自身翻折于 $|\delta\omega|=\sqrt{|2G_1/G_3|}\approx\Delta_{\rm fold}$;逼近时 Newton 的 $f'\to0$ 跳伪根(见过 −55 MHz 尖刺)。

### 5.4 缩脉冲的物理下限

泄漏(Ω→|α|=EC,2.5 ns 时 Ω/|α|≈0.5 需 DRAG)、RWA(Ω≪ω₀₁)、AWG 网格(T_{π/2}≫dt)。

---

## 6. 方法决策 + 与 Ramsey 对比

| 失谐范围 | 方法 | 理由 |
|---|---|---|
| \|Δ\| ≲ 0.5·Δ_fold(小失谐精修) | transient `order=3`(`fit`/`kernel_full`) | ~10× 精度,2 次 mesolve/点,极快 |
| \|Δ\| ~ Δ_fold | 缩短 π/2 脉冲推远翻折 | 量程↑(灵敏度同比↓) |
| \|Δ\| 大 / 未知 | **Ramsey 双扫 FFT** | 读振荡频率,无歧义、无界 |
| 兼顾量程+分辨 | 短/长脉冲组合解卷绕(CRT 式) | 短脉冲粗定位 + 长脉冲精测 |
| 任何情况 | **勿用** `diag_legacy`(已移除) | 错误物理对象 |

**Ramsey 无翻折**因为它测的是条纹**频率**(周期),不是有界的 $p_{\rm diff}$。
> **Ramsey 实现坑**:τ-扫的 `t_global` 必须覆盖完整序列(脉冲 + 最大 τ),否则 `hamiltonian_on` 截断丢弃 τ 处第二个 π/2 脉冲,扫描失效(误差卡 ~1 MHz 且加长 τ 不改善;覆盖后达 ~kHz)。

---

## 7. 实测数据汇总 + 复现 + 代码索引

**基准(10 ns π/2, 2 能级):** Δ_fold≈0.02 GHz, p_diff 极值±0.5, G₁≈−6.2, G₁·Δ_fold≈0.124, G₃ᵀᵃʸˡᵒʳ≈+1035(G₃/G₁≈−165 ns²), G₃ᵈⁱᵃᵍ 小 ~170×。

**测频精度(frequency_method_comparison.ipynb,RMS MHz):**

| bias | Ramsey | order1 | fit | kernel_full |
|---|---|---|---|---|
| 0.0(Δ<1MHz) | 0.002 | 0.015 | 0.0003 | 0.012 |
| 0.9(Δ<14MHz) | 0.0065 | 1.00 | 0.19 | 0.40 |

**核心公式速查:**
- 核(sim): $k_n=i^n\langle0|\,\mathrm{ad}_W^n Q\,|0\rangle$
- 测频: $f=\omega_d-p_{\rm diff}/G_{\rm freq}$(线性);$p_{\rm diff}=G_1\delta\omega+\tfrac16 G_3\delta\omega^3$(三次,δω=-Δ)
- G₃ 对象: $G_3^{\rm Taylor}=\iiint k_3$(≠ ∫k₃_diag)
- 翻折: $\Delta_{\rm fold}\approx\pi/(2T_{\rm eff})$;守恒 $G_1\Delta_{\rm fold}=$const;安全区 $|\Delta|\lesssim0.5\Delta_{\rm fold}$;误差 $\propto\Delta^{n+2}$

**代码索引:**
- `sqc/reconstruction/kernel.py` — `KernelEstimator`(mode/method/order/extract_off_diagonal, σ_t, richardson), `_heisenberg_kernels[_offdiag]`, `_extract_kn_*`
- `sqc/calibration/frequency.py` — `FrequencyMeasurement`(method, order, g3_source, g3_delta_max), `_measure_frequency_transient`, `_calibrate_g3_taylor`, `_calibrate_g3_kernel_full`
- `docs/architecture.md §4.6.2` — kernel 架构(Phase 12)
- 复现脚本: `kernel/verify_full_kernel.py`, `kernel/verify_transient_highorder.py`
- 流程图: `result/flowchart/kernel_multidim.tex`(多维核形成), `QSL_highorder.tex`(高阶闭环测频)
