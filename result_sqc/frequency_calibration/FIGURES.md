# 块 2 图注说明(频率标定 `frequency_calibration/`)

块 2 是 sqc **新增能力**,`src/` 无对应基准,故只画 sqc 自身物理正确性图
(与解析真值对比,不叠加 src)。每张图按四段式:**展示 / 如何读 / 物理含义 / 关键数值**。

> **目录约定(v2)**:每张图独占一个子文件夹 `F<n>_<slug>/`,内含该图的
> **脚本 + 数据(.npz)+ 图(.png/.pdf)**。子夹清单:
> `F1_flux_curve/`、`F2_transient_accuracy/`、`F3_step_method/`、
> `F4_hybrid_convergence/`、`F5_ramsey_vs_transient/`。

通用设置:Transmon `EC=0.2·2π, EJ=10·2π GHz, T1=100μs, T2=50μs`,基准工作点
`Φ_opt=0.9553`;`FrequencyMeasurement.measure(flux=offset)` 中 flux 为相对
`Φ_opt` 的偏置。⭐ = 论文重点。

---

## F1 — 频率-磁通标定曲线 f(Φ):Ramsey vs transient

图文件:`F1_flux_curve/f_phi_curve_sqc.png/.pdf`,数据:`F1_flux_curve/f_phi_curve_sqc.npz`

- **展示**:上图 f01 vs 磁通 Φ —— 解析曲线(黑实线)+ Ramsey 测点(蓝圆)+
  transient 测点(红空心方);下图残差(测量−解析,MHz)。扫描 Φ_opt±0.03,21 点。
- **如何读**:测点落在黑线上即准;残差子图看两方法偏离解析多少 MHz。
- **物理含义**:transmon 磁通色散 f(Φ)=√(8·E_J(Φ)·E_C)−E_C 是器件标定基元。
  两种测频法各自恢复它:
  - **Ramsey**:τ 扫描 + FFT,鲁棒但成本高(标准方法)。
  - **transient**:τ=0 正交读出 + 核灵敏度,便宜,但依赖弱信号线性(+三次)近似。
  图展示便宜的 transient 探针能在其有效窗口内复现 Ramsey 标定,窗口外近似失效。
- **关键数值**:全程 MAE Ramsey=0.019 MHz,transient=4.37 MHz。**但 transient
  误差集中在扫描边缘**:中心窗口(Φ≈0.94–0.975)transient 与 Ramsey/解析几乎重合
  (残差 <~1 MHz);左端 Φ=0.925 处残差达 28 MHz。
- **物理解读(可写入报告)**:驱动频率固定在 Φ_opt,偏置越大失谐 Δ 越大;
  单扫 transient 假设 |Δ|<0.1 GHz,大偏置(尤其左侧色散更陡)超出近似有效域 →
  误差骤增。这**定量刻画了 transient 测频的精度-成本-有效域权衡**。

---

## F2 ⭐ — transient 测频精度 vs 失谐:线性 vs 三次修正

图文件:`F2_transient_accuracy/transient_accuracy_sqc.png/.pdf`,
数据:`F2_transient_accuracy/transient_accuracy_sqc.npz`

- **展示**:横轴为**真实失谐** Δ=f(Φ)−f(Φ_opt)(MHz,由 DC 磁通偏置注入,
  驱动固定在 Φ_opt),扫 Φ_opt±0.028(13 点)。上图 measured−true f01(有符号,
  MHz);下图 |error| 对数 y。三条曲线:`order=1`(线性,蓝圆)、
  `order=3 g3_source="fit"`(奇多项式拟合 G₃,红空心方)、
  `order=3 g3_source="kernel_full"`(全离对角 ∭k₃ dt³,绿空心三角)。绿色阴影为
  三次可信窗口 |Δ|≲17.6 MHz(折叠 turnover),虚线为 1 MHz 参考。
- **如何读**:曲线越贴 0(上)/越低(下)越准。看窗口内 order=3 相对 order=1
  下压多少;看两条 order=3(fit vs kernel_full)是否重合(G₃ 两条独立路径的交叉验证);
  看窗口边缘处三条曲线是否重新合拢(Newton 折叠回退)。
- **物理含义**:transient 测频用正交 Ramsey 差分读出 p_diff 反演失谐。线性估计
  `δω=p_diff/G_freq` 在大 |Δ| 处因条纹折叠而**系统性偏低**;三次 Newton
  `p_diff=G₁·δω+(G₃/6)·δω³` 展宽可信域。G₃ 有两个**等价来源**:p_diff(Δ) 的奇多项式
  拟合(fit)与全离对角核三重积分(kernel_full)——两者一致,即
  [[transient-g3-sim-fullkernel]] 记录的 G1_fit≈0.6× 偏差的解决。折叠点
  |δω_fold|=√(−2G_lin/G₃) 之外三次不可逆,`_solve_cubic_detuning` 优雅回退到线性。
- **关键数值**:G₃ 标定 G1=−6.28、G3=1024.8(G3/G1=−163.3 ns²,fit 与 kernel_full
  一致);Δ 范围 −32.5–+17.0 MHz(**不对称**——Φ_opt 左侧色散更陡,同样 ±0.028
  偏置在低磁通侧产生更大 |Δ|,13 点中 10 点落在折叠窗口内)。窗口内 MAE:
  order=1=1.27 MHz,order=3(fit)=0.94 MHz,order=3(kernel_full)=0.97 MHz;
  Δ≈0 附近 order=3(fit)|error| 低至 ~10⁻³ MHz。窗口外(|Δ|≳18 MHz)三者收敛
  (order=3 回退线性)。
- **注意(诚实标注)**:正 Δ 折叠边缘(Δ≈+17 MHz)处三次略微**过冲**、|error|
  短暂超过线性——这是折叠边界的真实行为,非 bug,保留未做特殊处理。

---

## F3 — 闭环三步法收敛对比:secant / bisection / gradient

图文件:`F3_step_method/step_method_comparison_sqc.png/.pdf`,
数据:`F3_step_method/step_method_comparison_sqc.npz`

- **展示**:`SinglePointFrequencyCalibration` 用同一目标 f_target=f_opt−20 MHz
  跑三种根搜索步法。左图 |残差|=|f_q−f_target|(MHz,对数 y)vs 迭代,三条曲线
  secant(蓝圆)/ bisection(红方)/ gradient(绿三角),灰虚线为收敛容限
  ε=1e-4 rad·GHz≈0.016 MHz。右图 bisection 的括号宽度 V_hi−V_lo(对数 y)vs
  迭代,叠加理想 2⁻ⁿ 折半黑点线。
- **如何读**:左图看谁最快压到虚线下方(收敛速度)、曲线斜率(收敛阶);右图看
  实测括号宽度是否落在理想折半线上(bisection 的确定性 O(log₂) 行为)。
- **物理含义**:器件标定中把 f_q(V) 调到目标频率是把 r(V)=f_q(V)−f_target 求根。
  三种步法各有取舍:**secant** 超线性(用前两点割线,1–3 步),最快但需括号;
  **bisection** 每步折半,确定性 O(log₂(range/ε)),慢但最稳、趋势最干净(适合可视化);
  **gradient** 阻尼割线 + best-point 跟踪,**无需预括号**(只需 V_seed),鲁棒性介于两者。
  逐点测频由内部 `FrequencyMeasurement`(此处 ramsey 双扫,`f_artificial=None`,
  |Δ| 无界)完成。
- **关键数值**:target −20 MHz(根在偏置 offset≈−0.016);括号(偏置坐标)
  (−0.05, +0.01) 跨根。收敛步数 **secant=2、gradient=5、bisection=9**(教科书序);
  终残差全部 <ε:secant=0.009 MHz、bisection=0.012 MHz、gradient=0.015 MHz。
  bisection 括号宽度实测与理想 2⁻ⁿ 折半线**逐点重合**。
- **关键坑位(已解决,务必记住)**:`V_a/V_b/V_seed` 是相对 OPTIMAL_FLUX 的
  **偏置(offset)**,不是绝对磁通——`SinglePointFrequencyCalibration` 把 V 直接
  转给 `measure(flux=V)`,后者叠加到构造工作点上(见 [[sqc-frequency-measurement-api]]
  "flux is an OFFSET")。曾误用绝对磁通括号 (0.90, 0.98) → 实际总磁通 1.855 越界 →
  FFT 找不到峰 → 测频回退返回 f_opt → 残差冻结在 +20 MHz 不收敛。改回偏置坐标后
  三法全部快速收敛。
- **轴对齐说明**:各步法 history 起始 iter 不同(gradient 记 iter 0 种子点,
  括号法从 iter 1 起),已各自重标到 0 基迭代轴,故左图比较的是**收敛速率**而非
  绝对测量次数(secant 在进入主循环前已做 2 次括号端点测量)。

---

## F4 ⭐ — 闭环测频成本解剖:驱动跟随 + kernel_full 才是提速来源(4.1×)

图文件:`F4_hybrid_convergence/freq_calibration_convergence_sqc.png/.pdf`(收敛)
与 `freq_calibration_cost_sqc.png/.pdf`(成本),
数据:`F4_hybrid_convergence/freq_calibration_sqc.npz`

> ⚠️ **结论两次修正,以本节为准**:
> 1. 原计划标题"混合闭环**提速**"**不成立**:预设的 `sweet+fit` 混合方案在
>    任何有用容限下都**不快于**纯 Ramsey(实测 0.8×,即更贵)。
> 2. 但换配置后**提速真实存在**:`drive_policy="track"` + `g3_source="kernel_full"`
>    的**单阶段 transient**(不需要 Ramsey 精修)达到 **4.1× 更便宜且 3.5× 更准**。
> 3. ~~曾一度记为 333×~~ —— 那是**计数漏了 sesolve** 的错误。补上 `KernelEstimator`
>    内部的 sesolve 后真实开销是 984 次(不是 12 次)。**以 4.1× 为准。**

- **展示**:目标 f_target=f_opt−20 MHz,三种策略同起点(V_seed=0 偏置)。
  **收敛图**:|残差| 对数 y vs 迭代 —— 混合两相底色(coarse 绿 / fine 蓝)+
  switch 黑虚线,叠加 track+kernel_full(紫三角)与纯 Ramsey(红空心方)。
  **成本图**(三面板,x 轴一律为**实测**解算调用数):
  (a) |残差| vs 累计真实成本;(b) 各策略相对纯 Ramsey 的**加速比 vs 残差阈值**
  (向右容限更严),break-even=1×;(c) 收敛所需真实解算总数柱状图 + 各自的
  **独立核验残差**。
- **如何读**:(a) 同一残差高度上越靠左越省——紫线全程最左;(b) 读加速比:紫线
  全程在 1× 上方(3.7–7.3×),蓝线在容限收紧后**跌到 1× 以下**(变亏);
  (c) 直接读总账与真实精度。
- **物理含义**:混合思路是"便宜探针跑长途、贵探针收尾"。实测发现**真正的瓶颈不在
  分工,而在测量是否落在有效域内**:
  - **预设 sweet+fit 混合为何不赚**:gradient 步法有**固定 bootstrap 开销**
    (种子测量 + 盲探步 + 割线预热 ≈5 次测量),且迭代数**几乎与初始距离无关**——
    实测纯 Ramsey 在 −20 MHz 需 5 步、在 −60 MHz 仅需 6 步。所以 fine 阶段无论
    种子多好都要 ~5 次昂贵测量,coarse 的"领先"换不成节省。
  - **track+kernel_full 为何赚**:它把驱动 ω_d 每轮预测到下一测量点的 f_q,使
    **实测失谐始终 <1 MHz**(见收敛表 `delta` 列:−10.3 → −0.5 → −0.14 → +0.003),
    深居 transient 线性窗内 → 测量可信 → 梯度可信 → **6 步单调收敛,无需 Ramsey**。
    即"提速"来自**让便宜探针一直准**,而不是"便宜探针+贵探针接力"。
- **关键数值(全部经解析 f(Φ) 独立核验残差)**:

  | 策略 | 迭代 | 真实解算调用 | 真实残差 | vs 纯 Ramsey | 墙钟 |
  |---|---|---|---|---|---|
  | hybrid sweet+fit | 7+5 | 5214 | 0.00214 MHz | **0.8×**(更贵) | 49 s |
  | **transient track+kernel_full** | **6** | **984** | 0.00519 MHz | **4.1×** | **8 s** |
  | ramsey-only(基线) | 5 | 4000 | 0.01797 MHz | 1.0× | 46 s |

  按阈值的真实加速比(hybrid / track+kf):10 MHz → 4.2× / **7.3×**;
  5 MHz → 2.0× / **4.9×**;1 MHz → 0.9× / **3.7×**;
  ε=0.016 MHz → **0.8×** / **4.1×**。
  真实单次成本:transient(track+kf)≈164、Ramsey 双扫 ≈800 → 比值 ≈**4.9×**
  (**不是** workflow 记账暗示的 400×)。
- **⚠️ 成本记账陷阱(务必知道)**:`FrequencyCalibrationWorkflow` 的 `_stage_cost`
  把一次 transient 测量记作**固定 2**,严重低估:
  - `g3_source="fit"` + track:ω_d 每轮变化使 G₃ 缓存键 `(t_rabi, ω_d, tag)` 失效,
    每轮重标定 ~80 次 → 实测 1432 次 vs 记账 32 次(**低估 45×**)。
  - `g3_source="kernel_full"`:核估计走 **sesolve**(不在 mesolve 计数里),
    实测 984 次 vs 记账 12 次(**低估 82×**)。
  故本图 x 轴一律用**实测**计数(同时 patch `sqc.calibration.frequency.mesolve`、
  `sqc.reconstruction.kernel.mesolve`、`qutip.sesolve`)。sesolve 与 mesolve
  **按 1:1 计数**是保守做法(sesolve 是更便宜的调用),所以上表加速比是**下界**;
  墙钟 8 s vs 46 s(**5.8×**)是更贴近实际的收益。
- **coarse 相为何震荡(与 F2 的关联)**:收敛图里 coarse 在 20→4.5 MHz 间来回
  游走 7 步(非单调)。因为它从 **20 MHz 失谐**起步,而 F2 已量化 transient 的
  折叠有效边界是 **17.6 MHz** —— coarse 一开始就在有效域之外,测量有偏,梯度搜索
  自然游走。这是 F2 结论在闭环里的直接体现。
- **度量方法学(重要)**:**总成本不是公平指标**——两次运行停在不同终精度
  (混合 0.0050 MHz vs 纯 Ramsey 0.0151 MHz,混合多跑一步换来 3× 精度)。
  公平指标是"**首次达到给定残差的成本**",即成本图 (b)。另注:成本 trace 只对
  history 行计费,而 secant 步法在主循环前的 2 次括号端点测量**不入账**,故
  混合与基线均用 gradient 步法以保证同类比较。
- **配置说明(非调参凑结果)**:fine 阶段的盲探步长取 `switch_residual/灵敏度`
  =5 MHz/1202 MHz·Φ₀⁻¹≈0.0042 Φ₀,而非预设默认的 0.01 Φ₀。默认值在此处等于
  **12 MHz 盲踢**,会把搜索踢得比进入时更远,反而让混合总成本**高于**纯 Ramsey
  (实测 0.83×)。按"进入尺度"定标探步是 fine 阶段应有之义,故用显式 stages 覆盖
  (预设的扁平参数不会把步长转发给 fine 阶段)。
- **🔴 track 必须配 kernel_full,不能配 fit(机制性不兼容,已实测)**:
  `_fit_g3_at_delta_max` 的做法是把量子比特放在**甜点**、外加人工失谐扫 ±Δ,
  再拟合**奇多项式** `p_diff=c₁Δ+c₃Δ³+c₅Δ⁵` —— 这**假设 p_diff 关于 Δ=0 对称**。
  track 把 ω_d 移开后,扫描窗口中心被整体偏置 (f_sweet−ω_d),奇函数假设破裂:
  实测 G₁ 在迭代间取值 −6.28 / −0.99 / **+0.51 / +4.01** / …**反复穿零变号**,
  测量随之崩坏(残差跳到 63 / 165 / 124 MHz),V 早已锁定却仍乱跳,16 步不收敛,
  终点真实残差 4.07 MHz,白花 1432 次解算。
  `kernel_full` 直接用**实际 ω_d 的脉冲**算 `∭k₃ dt³`,无对称性假设,故与 track 相容。
- **结论(可写入报告)**:闭环测频的提速不来自"粗精两段分工",而来自
  **让廉价探针始终工作在其线性有效域内**——即 `drive_policy="track"` 的驱动跟随,
  外加一个不依赖对称性假设的 G₃ 来源(`kernel_full`)。这条链把 F2(折叠边界
  17.6 MHz)与 F5(order=3 的成本优势)在闭环场景里合并成一个可用结论。
- **改进方向(未实施)**:混合若仍要用,fine 阶段应改用**能利用好种子的步法**
  (如用已标定灵敏度做单次 Newton 步),而非盲探 bootstrap;另可给 `_stage_cost`
  加上 sesolve/G₃ 重标定的真实计费,避免记账误导。

---

## F5 — Ramsey vs transient:精度-成本帕累托前沿

图文件:`F5_ramsey_vs_transient/ramsey_vs_transient_sqc.png/.pdf`,
数据:`F5_ramsey_vs_transient/ramsey_vs_transient_sqc.npz`

- **与 F1 的分工**:F1 画"精度 vs **磁通**"(便宜探针能否复现色散曲线);
  F5 画"精度 vs **成本**"(各策略花多少次 mesolve 换到多少精度)。F5 的磁通窗口
  **刻意限制在 transient 有效域内**(|Δ|≲12 MHz,offset ±0.010,5 点),
  以隔离成本变量——折叠失效已由 F2 刻画,不在此重复。
- **展示**:6 种测频**策略**各为一个点。(a) MAE vs 每次测量的 mesolve 调用数;
  (b) MAE vs 每次测量墙钟秒数;均双对数,黑虚线为帕累托前沿。
  策略含 Ramsey 双扫/单扫 × τ 网格密度(Δτ=0.5/1/2 ns),以及 transient order=1/3。
- **如何读**:越靠**左下**越好(更便宜且更准)。落在别人右上方的点即被**支配**
  (花更多成本却更不准)。黑虚线连出的三点是不可被支配的最优选择。
- **物理含义**:两族测频法的成本结构完全不同。Ramsey 靠 τ 扫描+FFT,成本
  ∝ 扫描点数(双扫再 ×2),精度受 FFT 分辨率限制;transient 靠 τ=0 正交读出+
  核灵敏度,单次仅 2 次 mesolve,精度受弱信号近似阶数限制——**加三次修正
  (order=3)几乎不增边际成本却大幅提精度**,这正是 F2 的结论在成本维度的回报。
- **关键数值(实测 mesolve 计数,非估算)**:

  | 策略 | MAE (MHz) | mesolve/次 | 边际 | 秒/次 |
  |---|---|---|---|---|
  | Ramsey 双扫 Δτ=0.5 ns | **0.006** | 800 | 800 | 8.68 |
  | Ramsey 单扫 Δτ=0.5 ns | 0.023 | 400 | 400 | 4.41 |
  | Ramsey 单扫 Δτ=1 ns | 0.028 | 200 | 200 | 2.20 |
  | Ramsey 单扫 Δτ=2 ns | 0.030 | 100 | 100 | 1.06 |
  | Transient order=1 | 0.192 | 2 | 2 | 0.18 |
  | **Transient order=3 (fit)** | **0.014** | 18 | **2** | 0.38 |

- **⭐ 头条结论**:**transient order=3 帕累托支配全部三种 Ramsey 单扫配置**——
  MAE 0.014 MHz 优于 400 次求解的单扫(0.023 MHz),成本却低 **22×**(摊销)
  至 **200×**(边际)。仅 Ramsey 双扫更准(0.006 MHz),代价是 **44×**(摊销)/
  **400×**(边际)的求解量。即:要极致精度用双扫,其余场合 order=3 是更优选择。
- **摊销说明(重要)**:order=3 的 18 次是**摊销值**——G₃ 标定一次性花 80 次
  mesolve(模块级缓存,按 `(t_rabi, ω_d)` 键),边际成本仅 2 次。5 个磁通点摊销
  得 (80+2×5)/5=18;点数越多越接近 2(如 F1 的 21 点扫描 ≈5.8 次/点)。
  真实实验中这笔 setup 每种脉冲/驱动配置只付一次。
- **实现说明**:成本为**实测**——用上下文管理器把 `sqc.calibration.frequency.mesolve`
  换成计数包装(该模块 `from qutip import mesolve`,故必须替换**模块属性**,
  patch `qutip.mesolve` 太晚无效)。τ 网格按 R9 从 `CONFIG.awg.dt` 派生
  (`np.arange`,步长取 dt 整数倍)。粗 τ 网格会压低 FFT Nyquist=1/(2Δτ),
  故最粗只到 Δτ=2 ns(Nyquist 0.25 GHz),仍安全高于单扫的 f_artificial=0.1 GHz 峰。
- **坑位**:`measure(flux=...)` 收的是**相对工作点的偏置**,stub 原写
  `flux=OPTIMAL_FLUX+offset`(重复叠加,越界)——与 F3/F4 同一个坑,已修
  (见 [[sqc-frequency-measurement-api]])。
