# 短脉冲序列参与标定：论文写作计划

> **文件性质**：本文件是用户正式写作时使用的章节指南与讨论记录，不是论文正文，也不是已经确认的论文提纲。正式文字由用户亲自撰写；AI 仅对用户提供的文字进行润色、逻辑检查、证据核对和格式整理。开始使用前先阅读 [`../../WRITING_BOUNDARY.md`](../../WRITING_BOUNDARY.md)。

## A. 使用说明

本计划按照论文各章节组织。后续讨论形成的新决定，应补充到对应章节，不另建相互冲突的提纲。

状态标记：

- **[基本原则]**：贯穿全文，后续写作不能偏离。
- **[已确认]**：已经由用户明确决定。
- **[暂定]**：当前采用，但保留在后续证据或排版讨论中调整的可能。
- **[建议]**：当前可行方案，仍可在讨论后调整。
- **[待讨论]**：尚未决定，正式写作前需要确认。
- **[需代码调整]**：论文比较或实验依赖的能力尚未完成实现与验证。
- **[待实施]**：方案已确定或暂定，但对应代码、流程或数据产出尚未完成。
- **[待实验]**：定义已建立，具体数值或证据必须由真实实验冻结。
- **[已有素材]**：项目中已有数据、图片或参考资料。

每个章节统一使用以下栏目：

1. 本章任务；
2. 要回答的问题；
3. 建议的正文组织；
4. 图表安排；
5. 已有素材；
6. 完成检查；
7. 待讨论事项。

## B. 全过程基本叙事

**[基本原则]**

> 探究短脉冲序列参与标定的可能。

这一基本叙事不仅约束论文结构，还应持续决定：研究问题如何提出、材料如何取舍、证据如何排序、术语如何使用，以及结论延伸到哪里。响应核、drive tracking、闭环算法和具体数值结果只能作为可能的机制或证据，不能自行取代这一基本叙事。

## B.1 2026-08-06 第二、三章问答决策

以下决定覆盖本文件中更早的冲突性候选表述：

- **[已确认]** 原 Section III/IV 合并为一个理论与协议主体章节；第三章正式名称仍为 **[暂定]**。
- **[已确认]** Section II 采用三个功能小节：II.A 定义器件与标定任务，II.B 定义测量探针与基线，II.C 给出标定 workflow 和评价框架。目标与成功条件先于总体流程，避免先讲方法、后讲问题。
- **[已确认]** Section II 与 Section III 形成“完整架构--核心实例化”关系：Section II 说明闭环包含什么、各模块承担什么职责；Section III 只对短脉冲 local TRACK 链路给出响应模型、反演、磁通更新和 drive tracking 的技术实现。
- **[已确认]** 第三章保留 III.A--III.D 四个小节，drive tracking 单列 III.D。
- **[已确认]** `$G_3$` 留在主文；主文直接使用准静态标量响应的奇次展开 `$p_{\mathrm d}(\Delta)=G_1\Delta+G_3\Delta^3/3!+\mathcal O(\Delta^5)$`，不为未单独界定或使用的余项引入 `$R_5$`。
- **[已确认]** 主文同时保留一阶基线估计器和三阶主估计器；`$G_1$` 与 `$G_3$` 是零失谐处的一阶和三阶响应导数。若由响应核计算 `$G_3$`，正确对象是完整三时间积分。
- **[已确认]** 主文只简述 model route 与 Virtual-Z measurement route；一般 Volterra 泛函展开、嵌套对易子、时间排序、有限差分模板和三维采样细节进入 Supplemental Material。
- **[暂定]** 真机 `$G_3$` 采用离线 Virtual-Z 采样，在闭环中复用或插值；采样网格、插值变量、复用范围和重标定条件待实验确定。
- **[已确认]** 理想主模型保持奇次结构；实验零点偏置和分支不平衡作为 nuisance parameters 在 calibration set 标定，冻结后在独立 validation set 上检验。
- **[已确认]** 区分 `response fold` 与更严格的 `validated local operating interval`；后者的数值边界由独立验证数据确定。
- **[已确认]** 主文只完整介绍无预括号的阻尼割线/Newton 型偏置更新；不得称为 gradient descent，也不得写成对称双点 bootstrap。
- **[已确认]** Fig. 1(d) 和 Section II.C 总览完整 `ACQUIRE--TRACK--LOCK--REACQUIRE` 架构；Section III.C--III.D 不重新介绍流程，而是形式化 TRACK 的进入/退出接口、磁通更新和 drive-reference tracking。实际阈值进入 Appendix/实验方案。
- **[暂定]** ACQUIRE/REACQUIRE 使用双扫 Ramsey；TRACK 使用割线预测驱动参考；LOCK 将驱动参考置于目标频率。
- **[已确认]** LOCK 只承担达到目标后的确认与短时保持；长期频率稳定、闭环带宽、PSD 和 Allan deviation 放 Discussion/Outlook。
- **[已确认]** coarse-to-fine `transient -> Ramsey` 负面结果主要进入 Supplemental Material，正文只概述限制性结论。
- **[暂定]** fixed-drive 作为正文机制对照；最终位置随图数与实验结果确认。
- **[暂定]** Fig. 2 使用“差分响应 + 估计误差”两面板，也可在压缩时将响应改为 inset。
- **[已确认]** 不单设标定机制图；Fig. 1(d) 承担完整状态机，并区分偏置更新与驱动参考更新。
- **[已确认]** 资源结果同时报告 cold-start、amortized maintenance 和 local-tracking 三种成本，不能用局部 tracking 优势替代完整协议优势。
- **[已确认]** 固定延迟 Ramsey 进入正文强基线；其双正交 I/Q 单点实现属于 **[需代码调整]**，详见 II.B。

## B.2 2026-08-06 控制变量层级决策

以下决定统一第二、三章中主文物理坐标、实验执行量和现有数值变量的记号；它们覆盖此前引入通用控制命令 `$u$` 的候选方案：

- **[已确认]** Section II--III 主文统一使用 `$\Phi$` 作为控制器更新的 calibrated flux-bias coordinate，不再引入 `$u$`。为避免过度声称，`$\Phi$` 表示经电压--磁通换算得到的 flux-equivalent bias，而不是仪器直接测得的片上真实磁通。
- **[已确认]** 主文首次定义 `$\Phi$` 时只用一句话说明：实验通过偏置电压 `$V_{\mathrm b}$` 和单独标定的电压--磁通换算实现该坐标；换算关系、零点、周期、不确定度、漂移、DAC 分辨率和限幅进入 Appendix。
- **[已确认]** 标定目标写成寻找 `$\Phi_\star$` 使 `$|\omega_{01}(\Phi_\star)-\omega_{\mathrm{tar}}|\leq\epsilon_\omega$`。控制器可使用已校准的执行坐标 `$\Phi$`，但不获得定量色散 `$\omega_{01}(\Phi)$` 或其导数；已知 `$V_{\mathrm b}\mapsto\Phi$` 不等于已知 `$\Phi\mapsto\omega_{01}$`。
- **[已确认]** 旋转系失谐写成 `$\Delta(\Phi)=\omega_{01}(\Phi)-\omega_d$`。Section III 的割线灵敏度、偏置更新和 drive tracking 分别使用 `$S_\Phi=\partial\omega_{01}/\partial\Phi$`、`$\Phi_{n+1}-\Phi_n$` 和相邻实测点估计。
- **[已确认]** Fig. 1(a) 与主文闭环轨迹使用 `$\Phi/\Phi_0$` 或相对磁通偏置 `$\delta\Phi/\Phi_0$` 作为横坐标。实验点由原始 `$V_{\mathrm b}$` 数据换算得到；原始电压坐标及换算验证放入 Supplemental Material，主文图不强制设置电压副轴。
- **[平台查证记录]** 代码中的 `V_seed`、`V_a`、`V_b` 及历史字段 `V` 虽以 `V` 命名，实际由 `FrequencyMeasurement.measure(flux=...)` 作为相对磁通偏置使用，单位为 `$\Phi_0$`，不是伏特。论文正文统一改写为 `$\Phi_{\mathrm{seed}}$`、磁通边界和磁通步长；代码字段名只在复现说明中交代。
- **[成本边界]** 若实验闭环依赖预先获得的 `$V_{\mathrm b}\mapsto\Phi$` 换算，其标定成本必须计入 cold-start setup，或在复用时明确计入 amortized setup，不得把该执行坐标写成无成本先验。
- **[范围边界]** 本文只需处理闭环相关的准静态电压--磁通换算。控制线动态传递函数、脉冲失真和电压波形到片上瞬态磁通的完整模型放入 Appendix/Supplemental Material；除非实验表明它们决定主要结论，否则不扩展为第二章主线。

## B.3 2026-08-06 频率记号决策

- **[已确认]** Section II--III 的 Hamiltonian、响应核、估计器、标定目标、残差、灵敏度、停止阈值和 drive tracking 统一使用角频率 `$\omega$`、失谐 `$\Delta$` 与阈值 `$\epsilon_\omega$`，不在理论公式中交替使用 `$f$` 和 `$\omega$`。
- **[已确认]** 普通频率只作为实验与数值报告量，统一定义 `$f_{01}=\omega_{01}/(2\pi)$`、`$f_d=\omega_d/(2\pi)$`、`$\delta f=\Delta/(2\pi)$` 和 `$\epsilon_f=\epsilon_\omega/(2\pi)$`。图、表和正文数值可以使用 GHz/MHz，但必须明确已经除以 `$2\pi$`。
- **[已确认]** 理论资源指标写成 `$T_{\mathrm{hit}}(\epsilon_\omega)$`；当图表用 MHz 给出共同阈值时，图注可等价写 `$T_{\mathrm{hit}}(\epsilon_f)$`，但不得在同一推导中混用两种阈值。

## B.4 2026-08-07 短序列转角与方法材料分工

- **[已确认]** Bipartite short-pulse probe 使用一般转角 `$\alpha$`，序列写成 `$R_y(\alpha)$` 与 `$R_{\pm x}(\alpha)$`；`$\pi/2$` 只是可采用的参考设置，不是该序列的定义。
- **[已确认]** 转角 `$\alpha$` 是探针设计参数。在固定最大 Rabi 频率和脉冲形状下，减小 `$\alpha$` 可缩短序列并收窄时间响应，但会减弱布居响应；因此它调节时间分辨率、灵敏度和采集成本之间的权衡，并会影响 `$G_1$`、`$G_3$` 与 validated local operating interval。
- **[待实验]** 主比较采用的 `$\alpha$` 应在 calibration set 上按预先规定的精度--时间准则选择，随后冻结并在 validation set 上评价；Results 应至少给出转角选择依据或 `$\alpha$` 消融，避免把单一 `$\pi/2$` 设置写成方法本身。
- **[已确认]** `$\alpha$` 保留给短序列转角；若正文需要 Transmon 非谐性，改用 `$\alpha_{\mathrm{anh}}$`，置信分位数改写为 `$z_{1-\beta}$`，避免符号重用。
- **[已确认]** Section III 是本文提出的 short-pulse local TRACK 方法。`Methods` 不再作为含混的泛称：可复现但不构成创新主线的实验和数值步骤放入主文后的 `Appendix: Experimental and numerical methods`；长篇响应核推导、额外对照和扩展数据放入独立 Supplemental Material。候选稿中的交叉引用统一指向 Appendix 或 Supplemental Material，不能只写未定位的 `the Methods`。

## B.5 2026-08-07 Section III.A 响应展开与证据边界

以下决定覆盖本文件中关于主文一般 Volterra 展开、`$R_5$` 和时间对角线对照的更早候选表述：

- **[已确认]** 标量频率估计采用准静态条件：单个短序列内的失谐变化相对序列响应可忽略，写作 `$\Delta(t)\simeq\Delta$`；这不排除重复测量或反馈迭代之间的漂移。快速序列内波动若不可忽略，将破坏标量三次模型，作为成立条件或限制讨论。
- **[已确认]** 奇次结构来自理想 phase-cycled 两分支的交换对称性 `$p_{+X}(-\Delta)=p_{-X}(\Delta)$`，从而 `$p_{\mathrm d}(-\Delta)=-p_{\mathrm d}(\Delta)$`；中心有限差分只是响应导数的一种数值/实验估计手段，不是奇对称性的物理来源。
- **[已确认]** 主文不必先展示无限阶 Volterra 泛函展开。正文从准静态条件与分支交换对称性直接得到标量奇次展开，并保留 `$G_1$`、`$G_3$`，因为它们分别决定局部灵敏度、三次非线性修正和 response fold。
- **[已确认]** 主文可给出 `$G_1=\left.\partial_\Delta p_{\mathrm d}\right|_0$` 与 `$G_3=\left.\partial_\Delta^3p_{\mathrm d}\right|_0$` 的操作性定义，并在介绍 kernel route 时给出 `$G_1=\int k_1$` 与 `$G_3=\iiint k_3$`。核 `$k_m$` 已由完整短脉冲探针及其差分输出定义，无需使用上标 `$(\mathrm d)$`。
- **[证据边界]** 当前主结果没有 full three-time `$G_3$` 与 time-diagonal `$k_3(t,t,t)$` 近似的直接对照，因此“时间对角线不充分”不作为主结果或突出贡献。若主文保留核路径，只陈述并使用数学上正确的完整三重积分；对角线区别、推导与采样细节进入 Supplemental Material。
- **[证据边界]** 当前可由结果支持的表述是：kernel-derived `$G_3$` 与独立 odd-polynomial-fit 路径在局部响应/估计结果上交叉验证。Abstract、Introduction 和 Discussion 不得突出 full-vs-diagonal 优势，除非后续加入直接对照数据。
- **[待讨论] 全文统筹** 当前 Section III.A 草稿暂时保留 `calibration set` 与 `validation set` 的写法，并分别叙述 `$G_1/G_3$` 的获得和 `$b_0$` 的零点修正；本轮不提前合并或重命名。主要未决问题是这两个 set 在正文中首次出现时尚缺少数据来源、采集时机、划分原则和职责定义，且容易与 Section II 的 closed-loop independent success verification 混淆。
- **[待讨论] 后续专项审查** 完成 Section III 后文和 Results/Appendix 相关内容后，汇总全文中 `calibration set`、`validation set`、`calibration data`、`validation data`、`independent validation` 与 `independent verification` 的全部用法，再统一决定：（1）是否在 Section II.C 或 III.A 首次定义；（2）是否改称 `response-calibration data` 与 `response-validation data`；（3）参数拟合、偶分量诊断、`$\Delta_{\mathrm{val}}$` 确定和闭环成功核验分别由哪组数据承担；（4）`frozen` 是否统一改为更正式的 `held fixed`。在该专项审查完成前，不把当前 set 术语视为最终定稿。

## B.6 2026-08-08 Section III 第二版同步与初步审查

- **[已确认]** 当前 `paper/sections/03_kernel_correction.tex` 是第三章第二版正式审查对象，已写至 III.D；`~~~` 中的旧候选稿不再代表当前章节进度。
- **[已确认]** 第三章四个小节采用 `Short-pulse differential response`、`Local frequency estimator`、`Flux-bias update` 和 `Drive-reference tracking`。Section II.C 定义完整架构与评价原则；III.C 只形式化 TRACK 的进入/退出接口和物理磁通更新；III.D 只解释固定 drive 的局部范围问题及 drive-reference tracking。
- **[已确认]** TRACK 到 LOCK 的转换与 LOCK 内确认分开：环内残差满足 target-entry condition 后进入 LOCK；连续 `$N_{\mathrm{lock}}$` 次独立核验用于在 LOCK 内确认成功，而不是作为“进入 LOCK”的条件。确认失败且局部条件仍成立时返回 TRACK，否则进入 REACQUIRE。
- **[已确认]** 偏置更新使用一个实测 seed 加一个有界单侧 blind step 初始化，随后由最近两个实测点估计 secant sensitivity。`$\chi_\Phi\in\{-1,+1\}$` 只表示已知局部单调方向，blind step 只使用 `$\operatorname{sgn}(\widehat r_{\mathrm{seed}})$`，不使用定量斜率。seed 统一写 `$\Phi_{\mathrm{seed}}$`，避免与磁通量子 `$\Phi_0$` 混淆。
- **[已确认]** 主文把控制律写成 damped secant--Newton proposal；最大单步和安全磁通范围用正文文字说明，精确 clipping 约定与参数进入 Appendix。不得沿用代码内部历史名称 `gradient update`，也不必在正文强调“not gradient descent”。
- **[已确认]** III.D 必须明确机制：drive tracking 将下一轮失谐从累计频率位移转化为一步频率预测误差；它降低越界风险，但在割线预测、置信条件或局部反演失效时不保证仍在范围内。
- **[已确认]** 第三章作为方法章，应在更新职责、状态依赖和失败边界处结束。fixed-drive 消融安排、Section IV 的检验顺序、当前数值实现只覆盖 seeded TRACK 的证据声明，统一移到 Results/Discussion，不在方法章提前评价结果。
- **[待实施]** 第二版 III.D 末尾现有两组重复的“flux-bias 与 drive-reference 职责”和 workflow 收束段。最后加入的 workflow 段位置正确，可作为第三章结尾；下一轮正文整理时只保留一组，并把 LOCK 的独立确认与 REACQUIRE 失败路径合并在同一收束段中。
- **[待实施]** `paper/main.tex` 仍在第三章后载入旧的 `sections/04_feedback_protocol`。该文件重复 III.C--III.D，使用旧 `$V$`、`$\omega_q$` 和 `gradient` 术语，并造成 `eq:drive-tracking` 重复标签及 Results 章节顺延。正式整理时必须移除该输入或明确改作其他用途，不得继续作为独立协议章。
- **[待实施]** III.B 当前同时写“branch/root failure 回退线性估计”和“越界触发 REACQUIRE”，失败政策互相冲突。线性回退仅可在其自身独立验证区间内使用；否则 root/branch failure 应拒绝局部估计并进入 REACQUIRE。
- **[待实施]** III.D 当前先以 fixed-drive 特例 `$\omega_{d,n}$` 定义 `$\Delta_{n+1}$`，随后又称将 tracking 公式代入“该定义”。下一轮应先给一般定义 `$\Delta_{n+1}=\omega_{01}(\Phi_{n+1})-\omega_{d,n+1}$`，再令 fixed drive 满足 `$\omega_{d,n+1}=\omega_{d,n}$`，避免代入对象不一致。
- **[待实施]** Section II 仍需补回控制器不可使用定量 `$\omega_{01}(\Phi)$` 的明确限制，修正 `$\epsilon_\omega$`、失谐首次定义、Fig. 1 标签拼写及遗留占位字符，保证第二章向第三章提供完整接口。

## B.7 2026-08-08 Fig. 1(d) workflow 图示决策

- **[已确认]** Fig. 1(d) 采用紧凑状态图，状态框只保留 `ACQUIRE`、`TRACK`、`LOCK` 和 `REACQUIRE`；各状态使用的 probe 不在框内重复标注，而由 Fig. 1 caption 和 Section II.C 说明：`ACQUIRE/REACQUIRE` 使用 scanned Ramsey，`TRACK` 使用 short-pulse local probe，`LOCK` 使用 independent frequency verification。
- **[已确认]** 图中保留 `ACQUIRE -> TRACK` 与 `REACQUIRE -> TRACK` 的 local-entry transition、`TRACK -> LOCK` 的 target-entry transition、`TRACK -> REACQUIRE` 的 local-guard failure、`LOCK -> TRACK` 的 failed-verification/local-valid path、`LOCK -> REACQUIRE` 的 failed-verification/local-invalid path，以及 `LOCK` 内 `$N_{\mathrm{lock}}$` 次连续独立核验后的成功确认。
- **[已确认]** `TRACK` 在目标尚未达到且局部条件仍成立时继续迭代；为控制图面密度，不单独画 self-loop，由 caption 或正文说明。
- **[已确认]** `TRACK` 框内只以最小视觉编码区分 `$\Phi_{n+1}$` 的 physical flux-bias update 和 `$\omega_{d,n+1}$` 的 drive-reference update；具体职责仍由 caption 和 Section III.C--III.D 定义。
- **[待讨论]** 若 `ACQUIRE` 或 `REACQUIRE` 的宽范围独立测量已经满足 target-entry condition，是否允许直接进入 `LOCK`，尚未确认，因此当前 Fig. 1(d) 不增加该直达路径。

## B.8 2026-08-08 Section II--III 综合审查与待办

本节记录第三章第二版完成后对 Section II--III 进行的源码、图表、引用和状态机综合审查。以下条目不重设 B.1--B.7 已确认的章节结构与方法决定；`[建议]` 和 `[待讨论]` 仍需用户逐项确认，确认后再更新对应正式正文。

### B.8.1 编译结构与 Section II 接口

- [ ] **[待实施]** 从 `paper/main.tex` 移除旧 `sections/04_feedback_protocol` 的正式输入。该文件重复 III.C--III.D，使用旧 `$V$`、`$\omega_q$` 和 `gradient` 术语，造成 `eq:drive-tracking` 重复标签，并使 Results 顺延为 Section V；旧文件是否归档或删除可在停止编译后另行决定。
- [ ] **[待实施]** 将标定成功条件统一为 `$|\omega_{01}(\Phi_\star)-\omega_{\mathrm{tar}}|\leq\epsilon_\omega$`，修复当前裸文本 `epsilon`，并在普通频率图表中只使用已定义的 `$\epsilon_f=\epsilon_\omega/(2\pi)$`。
- [x] **[已实施]** 在 II.A 明确区分物理上存在、可用于模拟 ground truth 或独立核验的 `$\omega_{01}(\Phi)$`，与控制器运行时不可调用的定量色散及其导数；已知 `$V_{\mathrm b}\mapsto\Phi$` 不等于已知 `$\Phi\mapsto\omega_{01}$`。
- [x] **[已实施]** 在旋转系 Hamiltonian 附近就地定义 `$\Delta(\Phi)=\omega_{01}(\Phi)-\omega_d$`，保持 Section II--III 全部理论公式使用角频率。
- [ ] **[部分完成]** 已修复正式正文中的 Fig. 1(b)/(c) 标签、`room-temperatue`、句号后缺空格、主谓一致和单复数，并拆分 REVTeX 双栏中溢出的 `$G_1/G_3$` 公式；`********Figure~\ref{fig:system-protocal}(a)` 作为用户保留的占位行暂不修改。
- [ ] **[写作提醒]** 正文已经通过 `as detailed/specified/given in Appendix/Supplemental Material` 作出实质复现承诺，后续撰写时必须逐项兑现。Appendix 至少覆盖电压--磁通换算、Ramsey scan/I--Q reconstruction、`$\tau_0$`、response calibration 与 coefficient reuse/recalibration、blind step、精确 clipping、置信参数和 `$N_{\mathrm{lock}}$`；Supplemental 至少覆盖一般 response-function 推导、时间排序与嵌套对易子、Virtual-Z 有限差分和完整三时间采样。写作时须同步清理旧 `damped-gradient`、`$A_1/A_3$` 术语，并严格区分 protocol specification 与实际执行证据：固定延迟 I/Q Ramsey 和自动 REACQUIRE/独立 LOCK 核验在实现前不得写成已完成方法或结果。
- [ ] **[待实施]** 在 II.C 补齐评价框架接口：`$T_{\mathrm{hit}}(\epsilon_\omega)$`、独立核验原则以及 cold-start、amortized maintenance、local-tracking 三类成本；固定延迟双正交 I/Q Ramsey 在代码完成前只能写成计划基线，不能写成已有结果。

### B.8.2 Section III 响应、反演与更新

- [ ] **[待实施]** III.A 按 B.1/B.5 的材料分工，从准静态条件和 `$+X/-X$` 分支交换对称性直接建立标量奇次展开；一般 Volterra 泛函展开、时间排序和嵌套对易子移入 Supplemental Material。
- [ ] **[待实施]** 在 III.A 开头说明 `$k_m(\ldots;\alpha)$`、`$G_m(\alpha)$` 和 `$\Delta_{\mathrm{val}}(\alpha)$` 的依赖；`$\alpha$` 在 response-calibration data 上按预定准则选择并 held fixed 后，正文才省略其显式依赖。
- [x] **[已实施]** 已合并 III.A 两次 coefficient-route 与 Supplemental 介绍。direct-response/odd-polynomial-fit route 与 kernel route 只说明一次；kernel route 下区分 calibrated-control-model calculation 与 offline Virtual-Z measurement sampling。McKay 2017 只支撑 Virtual-Z 操作，full-three-time 与 time-diagonal 的区别未提升为当前主结果。
- [ ] **[待实施]** 明确 response-calibration data、response-validation data 与 closed-loop independent verification 的来源、采集时机和职责；参数拟合、偶分量诊断、`$\Delta_{\mathrm{val}}$` 确定和最终成功核验不得由同一数据概念含混承担。
- [ ] **[方法边界；需代码调整]** III.A 作为推荐方法保留 `$b_0$` 常数零点修正、偶分量诊断和 branch-contrast imbalance 保护，不因当前代码尚未实现而降低方法层级。需补充 model-mismatch 后的拒绝/重标定政策，并避免把偶分量唯一归因于 branch imbalance；代码与实验流程应向上实现 `$b_0$` 标定、对称失谐点的偶分量检验、必要的 branch-resolved contrast check 及其阈值测试。Results/Discussion 必须另行说明当前证据覆盖到哪一层，不得把推荐方法全部写成已经执行。
- [x] **[已实施]** III.B 已取消 cubic root/branch failure 后的运行时线性回退。线性估计器保留为 first-order baseline；只有当离线 response calibration 判定 `$G_3$` 在规定精度下可忽略、且线性估计器具有自己的独立验证区间时，才可预先选用。cubic convergence、branch continuity 或 range validation 失败时拒绝局部估计并进入 REACQUIRE。
- [x] **[已实施]** 已将 `$\Delta_{\mathrm{val}}<|\Delta_{\mathrm{fold}}|$` 写成只在 `$G_1G_3<0$` 且实 fold 存在时成立的附加条件；没有实 fold 时仍保留由 validation data 确定的 validated interval，并明确 cubic response 失去一一对应性而非正向响应变成多值函数。
- [ ] **[待实施]** 明确 `$z_{1-\beta}$`、`$\sigma_{\Delta,n}$` 的统计语义和覆盖率；区分响应斜率 `$G_1$` 太小与磁通灵敏度 `$S_\Phi$` 异常两类失败。
- [ ] **[待实施]** 在 III.C 区分 `proposal` 与 clipping 后实际施加的 `$\Phi_{n+1}$`，并确认 III.D 的 drive prediction 使用实际施加值。
- [ ] **[待实施]** III.D 先定义一般 `$\Delta_{n+1}=\omega_{01}(\Phi_{n+1})-\omega_{d,n+1}$`，再令 fixed drive 满足 `$\omega_{d,n+1}=\omega_{d,n}$`；合并 III.C--III.D 重复的更新职责和 workflow 收束段。

### B.8.3 Fig. 1 与后续方法图

- [ ] **[建议]** Fig. 1 的四面板框架已经足够，不增加第五面板，也不把核函数图塞入 Fig. 1(c) 角落。Fig. 1 只承担 physical task、probe timing、local differential response 和完整 workflow 四项总览职责。
- [ ] **[待实施]** Fig. 1(a) 展示 `$f_{01}(\Phi)$`、seed/target/迭代点及 `withheld from controller`；不得把 `$\Delta_{\mathrm{val}}$` 无条件画成固定的磁通区间，除非明确它相对于某个 `$\omega_d$` 的映射。
- [ ] **[待实施]** Fig. 1(b) 承担 scanned Ramsey、fixed-delay I/Q Ramsey 和 `$R_y(\alpha)\to R_{\pm x}(\alpha)$` 两分支定义；`$\pi/2$` 只能作为参考设置。
- [ ] **[建议]** Fig. 1(c) 优先展示原始 `$p_{\mathrm d}(\Delta)$`、独立 validation points、零点斜率和 validated interval，不叠加 linear/cubic/fold/estimator-error 全部细节。若横轴使用 MHz，写 `$\delta f=\Delta/(2\pi)$`，不能把角频率 `$\Delta$` 直接标成 MHz。
- [ ] **[建议]** Fig. 2(a) 补充较完整的原始差分响应、线性/三次模型、fold 与 validated interval，Fig. 2(b) 展示一阶/三阶估计误差；现有 `transient_accuracy.pdf` 的两个面板均为误差，尚不能单独承担原始响应证据。
- [ ] **[建议]** 若 kernel route 保留为可复现方法，在 Supplemental Fig. S1 展示 `$k_1(t)$`、明确定义的 `$k_3$` 切片或边缘积分，以及积分收敛到 `$G_1/G_3$` 的过程；避免用任意三时间切片装饰主图或暗示已有 full-vs-diagonal 实证比较。

### B.8.4 文献引用补充

- [ ] **[待实施]** II.A 的 transmon/控制背景补 `koch2007charge,krantz2019quantum`；慢频率漂移和旧 map 失准补 `paladino2014noise,vepsalainen2022improving`。
- [ ] **[待实施]** II.B 的旋转系驱动和 cQED 读出补 `krantz2019quantum,blais2021circuit`；Ramsey 原理补 `ramsey1950molecular`；fixed-delay I/Q、相位缠绕与反馈背景使用 `vepsalainen2022improving`。
- [ ] **[待实施]** short-pulse 两相移序列使用 `herb2024quantum`；有限控制时间下的 time-resolution/control-design 背景可使用 `herb2024quantum,lin2025application`，但不得让这些文献背书本文的三阶反演或三时间 `$G_3$`。
- [ ] **[待实施]** II.C 的校准编排和闭环背景补 `kelly2018physical,wittler_integrated_2021,vepsalainen2022improving`；校准时间和资源背景补 `rol2017restless,werninghaus_high-speed_2021`。
- [ ] **[待实施]** III.A 保留 `kubo1957statistical,mukamel1995principles` 作为一般 response-function 背景；准静态慢漂移可补 `paladino2014noise,vepsalainen2022improving`；Virtual-Z 操作补 `mckay2017efficient`。
- [ ] **[引用边界]** 不为 `$\Delta_n/r_n$` 双误差定义、cubic inversion、response fold、blind step、secant--Newton proposal、drive-reference predictor 或本文四状态转换规则强行附加既有文献。这些属于本文定义、初等推导或方法构造。McKay 2017 只支撑 Virtual-Z 操作，不自动支撑本文的三时间有限差分核采样。

### B.8.5 四状态机与实现证据边界

- [ ] **[已确认待实施]** 四状态机围绕两个误差量组织：`$\Delta_n$` 决定局部测频是否可信，`$r_n$` 决定相对固定目标的标定进度；drive-reference update 控制前者，flux-bias update 减小后者，二者不得混用。
- [ ] **[已确认待实施]** ACQUIRE/REACQUIRE 的退出条件是宽范围频率估计及其不确定度足以初始化 drive 并满足 local-entry condition，不是“已经接近目标”。
- [ ] **[待实施]** TRACK 到 LOCK 的 target-entry condition 应保证将 drive 切换到 `$\omega_{\mathrm{tar}}$` 后仍处于局部可信范围，例如 `$|\widehat r_n|+z_{1-\beta}\sigma_{r,n}\leq\Delta_{\mathrm{lock}}$` 且 `$\Delta_{\mathrm{lock}}\leq\Delta_{\mathrm{val}}$`；`$\Delta_{\mathrm{lock}}$` 不等于最终成功容差 `$\epsilon_\omega$`。
- [ ] **[已确认待实施]** LOCK 不继续承担常规磁通精调：它把 drive 置于目标频率并进行独立确认和短时保持。核验失败但局部条件仍有效时返回 TRACK；只有局部反演、范围或置信条件失效时才进入 REACQUIRE。
- [ ] **[已确认待实施]** 连续 `$N_{\mathrm{lock}}$` 次命中必须来自 independent verification，而不是同一环内短序列 estimator 的 `converge_streak`；角频率残差 `$r_n$` 与 `$\epsilon_\omega$` 比较，使用 MHz 时先定义 `$r_{f,n}=r_n/(2\pi)$` 再与 `$\epsilon_f$` 比较。
- [ ] **[待实施]** 运行时范围保护优先使用 independently validated `$\Delta_{\mathrm{val}}$` 或冻结的 `$\Delta_{\mathrm{guard}}\leq\Delta_{\mathrm{val}}$`。代码中的 `rho=0.6` 与 `linear_range` 可作为 Appendix 的实现参数，不能未经验证直接替代主文的 validated interval。
- [ ] **[证据边界]** 当前代码支持顺序 stage、`sweet/track/target` drive policy、`out_of_range` 标志和环内 `converge_streak`，但 `out_of_range` 只暴露给调用方，尚无自动 REACQUIRE 分支；`converge_streak` 也不是独立核验。Results/Discussion 应称其为 seeded TRACK 的现有证据和推荐四状态框架，不得声称端到端全自动状态机已经实现或验证。

### B.8.6 建议处理顺序

1. [ ] 先修复 `main.tex`、Section II 接口、标签和公式等编译/定义阻断项。
2. [ ] 再统一状态转换、失败政策、LOCK 独立核验和 `$\Delta/r$` 单位。
3. [ ] 整理 III.A--III.D 的材料分工、validated interval、proposal/applied update 与章末收束。
4. [ ] 冻结 Fig. 1(c)、Fig. 2 和 Supplemental kernel figure 的证据分工。
5. [ ] 按引用边界补文献，并逐条核对 Appendix/Supplemental 是否兑现正文交叉引用。

## C. 暂定章节导航

> **注意**：以下章节名称和拆分方式均为 **[建议]**，用于组织写作计划，不代表正式论文结构已经确定。

| 章节 | 暂定作用 | 主要图表 | 当前状态 |
|---|---|---|---|
| Section I. Introduction | 提出为什么要探究短脉冲序列参与标定 | 通常不单独放结果图 | 待后续讨论 |
| Section II. System and calibration task | 定义器件与任务、测量探针、完整 workflow 和评价框架 | Fig. 1 问题、探针与协议图 | 三小节结构已确认 |
| Section III. Short-pulse frequency inference and calibration | 实例化 local TRACK：短脉冲响应、局部反演、磁通更新和局部范围维持 | 响应与反演图、机制对照 | 结构已确认；名称暂定 |
| Section IV. Results and comparison | 展示精度、范围、收敛和资源成本 | 性能与基线比较图 | 待后续讨论 |
| Section V. Discussion | 解释可能性成立的条件、限制和可推广边界 | 原则上不新增核心结果 | 待后续讨论 |
| Section VI. Conclusion | 回答研究问题，不引入新结果 | 无 | 待后续讨论 |
| Appendix: Experimental and numerical methods | 放置复现实验与数值结果所必需的步骤和参数 | 参数表、流程定义 | 内容待写；位置已确认 |
| Supplemental Material | 放置长推导、额外对照和扩展数据 | 补充推导、补充图 | 内容待写 |
| Title / Abstract | 在正文稳定后概括全文 | 无 | 最后撰写 |

---

## Section I. Introduction

### I.1 本章任务

**[待讨论]** 建立研究背景并逐步收束到核心问题：短脉冲序列除了执行控制操作，是否还能在标定过程中提供有用信息。

### I.2 要回答的问题

- 为什么量子比特需要重复标定？
- 传统频率测量和标定流程解决了什么问题？
- 为什么值得考察短脉冲期间产生的参数响应？
- 现有工作已经做到哪里，本文实际要检验的缺口是什么？
- 本文的研究问题、范围和证据类型是什么？

### I.3 建议的正文组织

**[待讨论]** 后续逐段确定。暂不沿用任何既有 AI 草稿的段落结构。

建议至少区分：研究背景、已有测频与反馈方法、短脉冲响应带来的问题、本文研究问题与文章安排。

### I.4 图表安排

- Introduction 原则上不单独放性能结果图。
- Fig. 1 可以在 Introduction 末尾被首次引用，但图的主要解释放在 Section II。

### I.5 已有素材

- Vepsäläinen et al. (2022)：闭环频率稳定与 Ramsey 测频背景。
- Sung et al. (2021)：器件和物理模型的图文组织参考，不代表本文需要采用相同器件介绍。

### I.6 完成检查

- [ ] Introduction 提出的核心问题能够在 Conclusion 中被直接回答。
- [ ] 没有预先宣称短脉冲方案必然优于传统方法。
- [ ] 没有让 response kernel、反馈算法或某个数值指标替代全过程基本叙事。

### I.7 待讨论事项

- **[待讨论]** 文献脉络如何分组。
- **[待讨论]** Introduction 中对“标定”的定义需要多具体。
- **[待讨论]** 本文贡献在引言末尾采用问题式还是结论式表达。

---

## Section II. System and calibration task

> **[结构已确认；名称暂定]** 章节名暂用 `System and calibration task`，采用三个小节：`II.A Flux-tunable transmon and calibration task`、`II.B Ramsey and bipartite short-pulse probes` 与 `II.C Calibration workflow and evaluation framework`。三节依次回答“标定什么”“频率信息从哪里来”“这些信息如何进入完整流程并被评价”。

### II.1 本章任务

在进入短脉冲响应推断与 TRACK 技术实现前，用紧凑篇幅建立完整研究框架。本章不是教材式器件介绍，也不提前推导控制律；应让读者理解：

1. 标定对象是磁通可调 Transmon 的角跃迁频率 `$\omega_{01}(\Phi)$`；控制器更新校准后的磁通偏置坐标 `$\Phi$`，目标量是 `$\omega_{\mathrm{tar}}$`；
2. 本文研究闭环频率标定，而非锁定后抑制漂移的频率稳定；
3. 控制器在不使用定量 `$\omega_{01}(\Phi)$` 关系或其导数的条件下，通过测量和磁通偏置更新寻找目标工作点；
4. 短脉冲序列的频率信息来自脉冲作用期间的失谐响应，而不是自由演化时间扫描；
5. 完整 workflow 由宽范围捕获、局部跟踪、目标确认和失败恢复组成；第三章只实例化其中的短脉冲 local TRACK 链路；
6. “参与标定”是指短脉冲差分响应提供局部频率信息，并进入后续偏置更新，不预先宣称其替代全部传统测频方法。

本文包含真实器件实验，并以一致的两能级有效模型组织主文中的理论与实验。多能级效应、泄漏、`$T_1/T_2$`、SPAM、有限采样和开放系统模型原则上放入 Appendix 或 Supplemental Material；若后续实验表明其中某项决定主结论，再提升到主文。

### II.2 II.A Flux-tunable transmon and calibration task

#### II.2.1 器件、磁通偏置坐标与静态色散

主文保留一般非对称 SQUID 的色散关系：

$$
E_J(\Phi)=E_{J,\Sigma}\sqrt{\cos^2(\pi\Phi/\Phi_0)+d^2\sin^2(\pi\Phi/\Phi_0)},
$$

$$
\hbar\omega_{01}\simeq\sqrt{8E_CE_J(\Phi)}-E_C.
$$

需要明确区分：当前数值模拟采用 `$d=0$`；实验曲线使用光谱数据拟合参数；完整实验 `$f_{01}(\Phi)$` 曲线只作为 withheld ground truth 和独立核验，不作为控制器输入。Transmon 的主要器件引用采用 Koch et al. (2007)。主文将 `$\Phi$` 定义为 calibrated flux-bias coordinate；实验中该坐标由偏置电压 `$V_{\mathrm b}$` 通过单独标定的换算实现。正文只保留一句边界说明，完整电压换算和不确定度进入 Appendix。

主文 Hamiltonian 移到 II.B 开头，作为 Ramsey 与短脉冲探针共享的测量动力学，而不是从标定算法直接引入。Hamiltonian 使用角频率，并统一记号：

$$
\omega_{01}=2\pi f_{01},\qquad \omega_d=2\pi f_d,\qquad
\Delta(\Phi)=\omega_{01}(\Phi)-\omega_d.
$$

公式使用角频率；图中的频率、失谐和误差使用 GHz 或 MHz，并在图注或坐标轴上标清是否除以 `$2\pi$`。

#### II.2.2 标定任务、先验和物理成功条件

本节把 frequency calibration 定义为 calibrated flux coordinate 中的闭环搜索：从初始偏置出发，通过后文定义的测量、频率估计和有界磁通更新寻找目标工作点。这里只定义任务，不展开具体残差、初始化和更新公式。

允许控制器使用的弱先验仅包括安全磁通区间、局部单调方向以及安全步长和驱动限制；不得使用定量 `$\omega_{01}(\Phi)$` 或已知导数。calibration 是到达目标工作点的过程；stabilization 是到达后抑制漂移的任务，不属于本文主要范围。

标定成功定义为

$$
\lvert \omega_{01}(\Phi_\star)-\omega_{\mathrm{tar}}\rvert\leq\epsilon_\omega.
$$

该式定义物理目标，而不是控制器的自我停止条件。内部估计残差与独立成功核验的关系放 II.C；`$\widehat r_n$`、盲探测步和割线更新只在 III.C 定义。

### II.3 II.B Ramsey and bipartite short-pulse probes

II.B 开头先定义 measurement probe 为“固定磁通偏置下的控制序列加布居读出”，再引入共同两能级 Hamiltonian 和 `$\Delta(\Phi)$`。随后依次介绍 Ramsey 与 bipartite short-pulse probes。测量链固定写为“驱动演化把 `$\Delta$` 编码为末态布居；measurement 获得布居；estimator 返回 `$\widehat\Delta$` 和 `$\widehat\omega_{01}$`”。

#### II.3.1 Ramsey 基线

Section II 只保留通用 Ramsey 响应：

$$
p_R(\Delta;\tau,\varphi)=\frac{1}{2}[1+C\cos(\Delta\tau+\varphi)].
$$

设置两个职责不同的 Ramsey 基线：扫描 Ramsey 用于宽范围捕获和独立核验；固定延迟 Ramsey 是主要局部实验基线。Vepsäläinen et al. (2022) 只支撑固定延迟 Ramsey 测量与频率稳定背景，不作为本项目阻尼割线优化器的来源。FFT、相位解缠、采样数优化和分支选择细节不放在本章。

**[需代码调整]** 当前 `FrequencyMeasurement(method="ramsey")` 只实现了
`$\tau$` 扫描 FFT 路径；现有 `DelayRamseyExperiment` 用于控制线尾波重建，
不能直接作为本文的固定延迟单点测频基线。正文主比较前需要新增固定延迟
Ramsey 测频模式，采用同一 `\tau_0` 下的双正交 I/Q 分支，由两路布居通过
`atan2` 恢复局部相位并反演失谐。`\tau_0` 应在 calibration set 上按单位时间
Fisher 信息或固定精度资源成本选择，随后冻结并在独立 validation set 上评价。
实现任务至少包括：`FrequencyMeasurement` 接口、shots/读出噪声支持、闭环调用
接口、符号与相位缠绕测试，以及与双分支短脉冲方法在相同 shots 和真实时间模型
下的公平比较。该路径完成并验证前，不得把固定延迟 Ramsey 写成已有数值结果。

#### II.3.2 Bipartite short-pulse probe

正式名称使用 `bipartite zero-delay Ramsey sequence`，行文简称 `bipartite short-pulse sequence`；两个分支称为 `phase-cycled +X and -X branches`。两段脉冲采用一般转角 `$\alpha$`，即 `$R_y(\alpha)$` 与 `$R_{\pm x}(\alpha)$`，而不是把序列固定定义为 `$\pi/2$` 脉冲。文献关系表述为该序列“following”或“as analyzed by Herb and Degen”，不声称 Herb and Degen (2024) 首次提出该序列。

在固定最大 Rabi 频率和脉冲形状下，减小 `$\alpha$` 会缩短序列、收窄时间响应并减弱布居信号。正文应明确把 `$\alpha$` 作为决定时间分辨率、灵敏度和采集成本的探针设计参数；主比较所用转角的选择和冻结规则放入 Appendix，转角消融或精度--时间权衡进入 Results。`$\pi/2$` 仅作为 Fig. 1(b) 或实验参数表中的一个参考设置。

两个分支的末态布居记为 `$p_{+X}$` 和 `$p_{-X}$`，实际差分观测量定义为

$$
p_{\mathrm d}=\frac{p_{+X}-p_{-X}}{2}.
$$

本章说明失谐在有限脉冲持续时间内改变末态布居，并可给出局部可辨识条件

$$
\left.\partial_{\Delta}p_{\mathrm d}\right|_{\Delta=0}\neq0.
$$

`$G_1/G_3$`、响应核、三次反演、翻折与详细有效性推导移到 Section III。需要持续区分脉冲序列、原始分支布居、差分观测量和由其得到的频率估计。

II.B 结尾给出“参与标定”的操作性定义：在已验证的局部工作区间内，短脉冲差分响应产生可用于估计当前频率的观测量，该估计随后参与偏置更新；这不等同于全局无歧义测频，也不等同于频率稳定。

### II.4 II.C Calibration workflow and evaluation framework

#### II.4.1 完整架构与模块接口

本节只给 workflow overview，不给控制器推导。一次 calibration cycle 在功能上由 measurement、estimation 和 tuning 组成：measurement 返回布居数据，estimation 返回 `$\widehat\Delta_n$` 与 `$\widehat\omega_{01,n}$`，tuning 改变 `$\Phi_n$`。Fig. 1(d) 展示完整 `ACQUIRE--TRACK--LOCK--REACQUIRE` 架构，并分开标出改变物理工作点的 flux-bias update 与改变测量参考的 drive-reference update。

逐状态只说明职责：scanned Ramsey 用于 ACQUIRE/REACQUIRE 和独立核验；短脉冲局部估计器用于 TRACK；fixed-delay Ramsey 是 TRACK 的局部强基线；LOCK 负责目标确认和短时保持。状态转换不等式、残差、割线更新和 drive tracking 公式留在 Section III。

#### II.4.2 独立核验与评价维度

闭环使用的局部估计器不能单独证明自身成功。模拟以 withheld transmon dispersion 核验，实验以独立频率测量核验。Section II 只建立这一非循环验证原则，并说明后续评价采用共同的精度阈值、局部有效范围和真实实验时间；不提前给出性能结论。

主要资源指标暂定为 `$T_{\mathrm{hit}}(\epsilon_\omega)$`，即独立核验首次以规定置信度确认误差进入阈值所需的墙钟时间。图表使用 MHz 时以 `$\epsilon_f=\epsilon_\omega/(2\pi)$` 表示同一阈值。核验频率、计时边界、cold-start/amortized/local-tracking 成本和置信规则进入 Section IV/Appendix；求解器调用次数只作数值诊断。

#### II.4.3 向 Section III 的交接

第二章结尾必须指出：完整 workflow 将 local TRACK 隔离为短脉冲可能参与标定的关键阶段。要实现这一点，必须同时解决两个问题：短脉冲非线性局部响应如何可靠反演，以及磁通偏置移动时如何让后续测量继续留在局部有效范围。Section III 分别以三阶局部估计和 drive-reference tracking 解决这两个耦合问题。

### II.5 Fig. 1 安排

Fig. 1 定义物理问题和标定协议，不展示最终精度、收敛速度或资源优势；主文 Fig. 1 不放芯片显微照片或完整低温测控链路。

| 面板 | 已确定或暂定内容 | 证据与边界 |
|---|---|---|
| Fig. 1(a) | 实验 `$f_{01}(\Phi)$` 数据与拟合，标出参考点、目标点和闭环实际访问的迭代点 | 完整曲线明确标注为 withheld from optimizer，只作物理背景与独立核验 |
| Fig. 1(b) | 扫描/固定延迟 Ramsey 与 bipartite short-pulse 时序，重点突出可调转角 `$\alpha$` 和 phase-cycled `$+X/-X$` 分支 | 不展开 Ramsey 数据处理算法；`$\pi/2$` 只可作为参考设置 |
| Fig. 1(c) | 只画 `$p_{\mathrm d}(\Delta)$`：理论/数值响应曲线叠加独立实验验证点及误差棒 | 标出零失谐、局部斜率和一个 `validated local operating interval`；`$p_{+X}$`、`$p_{-X}$` 移至补充材料，不在本图分别画线 |
| Fig. 1(d) | `ACQUIRE--TRACK--LOCK--REACQUIRE` 完整状态机 | **[已确认]** 同一面板区分磁通偏置更新与视觉次要的 drive-reference update 虚线路径；不另设标定机制图。图中不暗示自动重捕获已完成验证 |

Fig. 1(c) 不在第二章分别标注线性区、三次修正区和翻折区，这些边界留给 Section III。实验数据按职责拆分：calibration set 用于确定零点和响应系数；validation set 提供 Fig. 1(c) 的主实验点；在预测验证集前冻结 nuisance parameters，避免同一数据同时承担拟合与验证。

#### II.5.1 设备、能级与补充材料分工

- 主文只在公式或 Fig. 1 附近保留服务于失谐定义的 `$\lvert0\rangle,\lvert1\rangle$` 两能级表示。
- 完整多能级结构、泄漏、非谐性、开放系统参数和 SPAM/有限采样模型放 Appendix 或 Supplemental Material。
- 芯片照片、等效电路、低温测控链路和详细参数表放 Appendix 或 Supplemental Material；只有后续证据表明某项直接决定主结论时才提升到主文。

### II.6 参考与原创性边界

- Koch et al. (2007)：磁通可调 Transmon 的主要器件与色散参考。
- Krantz et al. (2019) 与 Blais et al. (2021)：超导量子比特驱动、旋转系有效模型和 cQED 读出的一般背景；不替代 Koch 2007 对 transmon 色散的直接支撑。
- Paladino et al. (2014) 与 Vepsäläinen et al. (2022)：低频噪声、慢频率漂移以及闭环频率反馈的依据；用于说明既有频率--磁通映射可能不足以满足当前目标容差，不声称 `$V_{\mathrm b}\mapsto\Phi$` 执行坐标同时失效。
- Ramsey (1950) 与 Vepsäläinen et al. (2022)：分别支撑 separated-field Ramsey 原理和 superconducting-transmon 中固定延迟 Ramsey 的局部频率估计与反馈使用。
- Herb and Degen (2024)：短脉冲序列及其线性感知分析的直接参考；使用“following/analyzed by”，不作“first proposed by”的优先权判断。
- Vepsäläinen et al. (2022)：固定延迟 Ramsey 测量及稳定背景。其闭环在已建立频率-磁通关系附近工作，与本文“不向优化器提供定量色散关系、通过迭代到达目标频率”的标定任务不同。
- Kubo (1957) 与 Mukamel (1995)：一般线性/非线性 response-function 展开的理论背景；不作为本文 `$G_3$` 三时间核的现成推导来源，详细三阶推导仍由 Supplemental Material 承担。
- McKay et al. (2017)：Virtual-`$Z$` 操作的直接来源；只支撑相位更新实现，不自动支撑本文的有限差分核采样方案。

### II.7 已有素材

- `report/img/fig_freq_calib.png`：只作内容清单参考；信息过密，不直接用于 Fig. 1。
- `result_sqc/frequency_calibration/report_figures/RF1_flux_curve_report.*` 与 `F1_flux_curve/f_phi_curve_sqc.*`：可作为 Fig. 1(a) 数据与版式基础，需加入实验拟合、目标点、迭代点和 withheld 标注。
- `result_sqc/frequency_calibration/F2_transient_accuracy/transient_accuracy_sqc.*`：适合后续精度与范围图，不能代替 Fig. 1(c) 的原始 `$p_{\mathrm d}(\Delta)$` 响应。
- `paper/figures/flux_curve.pdf` 和 `paper/figures/transient_accuracy.pdf`：既有 AI 草稿使用的冻结图，仅作素材，不自动决定正式图序。

### II.8 完成检查

- [ ] II.A--II.C 三小节形成“系统与任务--测量信息--workflow 与评价”的连续节奏。
- [ ] 不看后文，读者能够区分 measurement、calibration、tuning 和 stabilization，并说清标定量、受控量与目标量。
- [ ] 定量 `$\omega_{01}(\Phi)$` 关系未进入控制器；图中以 `$f_{01}=\omega_{01}/(2\pi)$` 报告的完整色散曲线只承担 withheld ground truth 和独立核验。
- [ ] 角频率/普通频率记号与 GHz/MHz 单位一致。
- [ ] 扫描 Ramsey 与固定延迟 Ramsey 的不同基线职责清楚。
- [ ] 短脉冲序列的 `$+X/-X$` 分支、相位和 `$p_{\mathrm d}$` 定义完整。
- [ ] 原始分支布居、差分响应和频率估计没有混为一谈。
- [ ] Fig. 1(c) 的 calibration set、validation set 和冻结参数规则明确。
- [ ] 没有在 Results 前提前宣称精度、成本或闭环优势。
- [ ] 主文两能级模型与实验器件、数值模型及补充材料中的非理想因素边界清楚。
- [ ] Section II 只定义完整架构和接口，没有提前重复 Section III 的残差、控制律和状态转换公式。
- [ ] 章末明确提出局部反演与局部范围维持两个问题，并自然交给 Section III。

### II.9 后续需要最终确认

- 章节最终英文名称；
- Fig. 1(d) 是否保留视觉次要的 drive-reference update 虚线路径；
- `$T_{\mathrm{hit}}(\epsilon_\omega)$` 的置信度、独立核验频率及墙钟时间计入规则；
- 实验 calibration/validation 数据的具体采集、划分和 nuisance-parameter 冻结方案。

---

## Section III. Short-pulse frequency inference and calibration

> **[已确认合并；名称暂定]** 原 Section III 和 Section IV 合并为一个顶层理论与协议主体章节，用小节串联“响应包含信息—局部反演—进入标定—维持有效范围”的连续证据链。它不是对 Section II workflow 的重复，而是对 Fig. 1(d) 中 local TRACK 核心链路的技术实例化。后续不再按旧 `.tex` 的 `Response-kernel calibration` 与独立反馈章节组织；第三章英文名称在全文主张稳定后最终确认。

### III.1 本章任务

本章开头先从 Section II 的完整架构中隔离两个耦合瓶颈：第一，短脉冲响应只在局部可逆且存在非线性，必须建立带边界的可靠反演；第二，磁通更新会移动量子比特频率，固定驱动会使下一次测量离开该局部范围。第三章围绕这两个问题依次回答：

1. 短脉冲差分响应中是否存在可用于频率标定的信息；
2. 如何在已知局部灵敏度时把该响应转化为频率估计；
3. 局部估计如何进入偏置更新；
4. 如何使多轮测量保持在短脉冲探针的有效范围内；
5. 哪些结论来自短脉冲响应，哪些依赖外部更新策略。

本章在论证地位上是核心方法章，但不把每个组成件都声称为原创：短脉冲序列和一阶感知关系有直接文献基础，普通割线更新是标准数值工具；本文需要论证的是三阶局部响应扩展、带验证边界的频率反演，以及它们与 model-free drive-reference tracking 的标定集成。

**[已确认]** 小节结构：

- `III.A Short-pulse differential response`
- `III.B Local frequency estimator`
- `III.C Calibration protocol`
- `III.D Maintaining the local operating range`

合并的目的不是增加方法内容，而是避免把核函数和控制算法误写成两个独立贡献。代码模块的边界不自动决定论文的章节边界。

### III.2 III.A Short-pulse differential response

#### III.2.1 本节任务

从 Section II 定义的两条短脉冲序列和差分读出出发，说明有限脉冲持续时间内的失谐如何改变末态布居，并使差分信号在参考工作点附近具有可辨识的频率响应。

#### III.2.2 建议的正文组织

第一段：承接 Section II 已定义的脉冲时序、两个分支布居和差分观测量，说明标量频率估计限制在准静态条件：单个短序列内 `$\Delta(t)\simeq\Delta$`，但允许重复测量和反馈迭代之间存在漂移。

第二段：利用理想 phase-cycled 分支的交换关系 `$p_{+X}(-\Delta)=p_{-X}(\Delta)$` 得到 `$p_{\mathrm d}(-\Delta)=-p_{\mathrm d}(\Delta)$`，从物理对称性而不是中心有限差分推出零阶和偶数阶消失，并写成

$$
p_{\mathrm d}(\Delta)=G_1\Delta+\frac{G_3}{3!}\Delta^3+\mathcal O(\Delta^5).
$$

第三段：先以 `$G_1=\left.\partial_\Delta p_{\mathrm d}\right|_0$` 与 `$G_3=\left.\partial_\Delta^3p_{\mathrm d}\right|_0$` 定义实际反演所需系数，再简述由已知 Hamiltonian 计算的 model route 与通过 Virtual-Z 探针采样末态布居的 measurement route。若正文展示核接口，使用 `$G_1=\int k_1dt$` 与 `$G_3=\iiint k_3dt_1dt_2dt_3$`，不再使用核上标 `$(\mathrm d)$`，也不在此处展开一般 `$G_m$`。

第四段：解释 `response fold` 与 `validated local operating interval` 的区别，明确后者由独立 validation data 确定，且局部响应不等同于全局测频能力。

实验中保持主理论的奇次模型；零失谐偏置 `$b_0$`、分支对比度差异和其他 nuisance parameters 在 calibration set 上估计并冻结，以 `$\widetilde p_{\mathrm d}=p_{\mathrm d}^{\mathrm{meas}}-b_0$` 进入反演。显著偶次残差作为模型失配诊断，不自动并入主估计器。

#### III.2.3 核函数的写作位置

**[已确认]** 不设置独立的 `Kernel theory` 大节。主文从准静态标量响应和分支交换对称性引出 `$G_1/G_3$`；一般 Volterra 展开、嵌套对易子、时间排序、有限差分模板和三维采样实现进入 Supplemental Material。

Herb 2024 已在弱信号、线性响应条件下提出感知核关系。引用仅支撑该范围；完整三阶三时间核需要独立高阶响应依据或本文 Supplemental Material 的自洽推导。

$$
G_1=\int k_1(t)\,dt .
$$

实验路径 **[暂定]** 为离线 Virtual-Z 采样：在代表性驱动频率或工作点获得 `$G_1,G_3$`，闭环中复用或插值。采样网格、插值变量、复用范围、重标定条件及摊销成本待实验确定；不得把离线采样写成逐轮实时三维核测量。

引用与原创性边界：

- **[引用边界]** Herb 2024 可支撑一阶、弱信号、线性响应下的感知核表述。
- **[引用边界]** Herb 2024 不能直接作为完整三阶三时间核 `$k_3(t_1,t_2,t_3)$` 的来源。
- **[已确认]** `$G_3=\iiint k_3\,dt_1dt_2dt_3$` 留在主文；其详细推导与实现进入 Supplemental Material。
- **[证据边界]** 主结果没有 full-vs-diagonal 的直接对照；正文不把排除 `$k_3(t,t,t)$` 写成已验证优势或主要贡献。

#### III.2.4 图表安排

- 若 Fig. 1(c) 已完整展示差分响应和局部可逆区，本节直接引用，不重复绘图。
- 若 Fig. 1(c) 仅作示意，则 Fig. 2(a) 展示 `$p_{\mathrm d}(\Delta)$`、零点斜率、线性/三次近似、翻折点和安全工作区。
- 本节不提前展示闭环成本或性能优势。

### III.3 III.B Local frequency estimator

#### III.3.1 本节任务

说明如何由短脉冲差分响应得到局部频率估计，并明确反演方法、可信范围、翻折和失败判据。

#### III.3.2 建议的正文组织

第一段：给出一阶基线估计器 `$\widehat\Delta_1=\widetilde p_{\mathrm d}/G_1$`，主文统一使用物理失谐 `$\Delta=\omega_{01}-\omega_d$`。代码内部 Virtual-Z 的 `$\delta\omega=-\Delta$` 约定只在 Appendix/Supplemental 说明。

第二段：说明一阶估计误差如何随失谐增加，以及奇次非线性如何导致系统偏差和最终翻折。

第三段：给出作为主方法的三次反演：

$$
p_{\mathrm d}=G_1\Delta+\frac{G_3}{6}\Delta^3 .
$$

正文说明从线性估计初始化 Newton 求根，并选择与零失谐连续的局部分支。当前实现的 fold guard 在越界根或退化三次项下回退到一阶估计；必须同时说明该保护不恢复全局可辨识性。嵌套对易子、时间排序和三时间数值积分放 Supplemental Material。

第四段：只使用 `response fold` 和 `validated local operating interval` 两层术语。fold 是三次模型的数学边界；validated interval 是满足截断误差、bias/RMSE、置信度和独立验证要求的实际工作区间，通常更严格。

#### III.3.3 图表与证据

**[暂定]** Fig. 2 以“局部估计的准确度和边界”为主题：

- Fig. 2(a)：原始差分响应，若 Fig. 1 未完整承担；
- Fig. 2(b)：一阶与三次反演的估计误差；
- 若图数需压缩，将原始响应改为估计误差面板的 inset；`fit` 与 `kernel_full` 的完整交叉验证优先进入 Supplemental Material。

**[证据现状]** 当前 F2 数据表明三次修正在约 `$|\Delta|\lesssim17.6\,\mathrm{MHz}$` 的折叠窗口内降低误差，但边界附近存在过冲，窗口外回退到线性。这支撑“局部改进及其边界”，不支撑全局优势。

### III.4 III.C Calibration protocol

#### III.4.1 本节任务

承接 Section II.C 已给出的完整 `ACQUIRE--TRACK--LOCK--REACQUIRE` 架构，形式化 local TRACK 的进入/退出接口与磁通更新；不重新介绍整套 workflow。严格区分测量、偏置更新和测量参考更新。

必须定义：

- 标定目标：目标频率或目标工作点；
- 测量输出：当前频率估计；
- 残差：当前频率与目标频率之差；
- 磁通偏置更新：改变 `$\Phi$`，使真实量子比特频率靠近目标；
- 驱动参考更新：改变下一轮测量的 `$\omega_d$`，使短脉冲探针保持在局部有效区。

偏置命令更新改变被标定对象的工作点；驱动参考更新改变测量坐标系或参考。二者不能合并为同一个更新机制。

**[已确认]** 主文给出四状态的完整协议设计和符号化转换条件：

- `ACQUIRE -> TRACK`：宽范围方法给出的估计及其不确定度落入 validated local operating interval；
- `TRACK -> LOCK`：连续 `$N_{\mathrm{lock}}$` 次满足目标残差阈值；
- `TRACK/LOCK -> REACQUIRE`：测量越界、根不连续、置信度不足、局部斜率异常或连续失败；
- `REACQUIRE -> TRACK`：重新获得满足局部区间条件的绝对频率估计。

`$\Delta_{\mathrm{val}}$`、置信系数、连续确认次数和实际阈值放 Appendix/实验方案，由 calibration/validation data 冻结。当前代码支持顺序多阶段编排和 `out_of_range`/predicate 接口，但尚未形成自动循环状态机；正文必须区分“完整协议设计”与“已验证的局部 tracking 结果”。

#### III.4.2 偏置更新策略

**[已确认]** 主文只完整介绍实际主协议使用的一种无预括号的阻尼割线/Newton 型更新。代码名称虽为 `gradient`，论文中不称为 gradient descent，因为其更新依赖局部割线斜率和最佳点跟踪。

主文需要写清：

- 初始种子处先完成一次频率测量，再按残差符号执行有界的单侧盲探测步；后续相邻实测点才提供割线斜率。不得误写成对称双点 bootstrap；
- 步长如何阻尼或限制；
- 收敛条件和失败条件；
- 每轮短脉冲估计如何产生下一次偏置更新。

主文所有控制器公式使用 `$\Phi$`：相邻实测点给出 `$\widehat S_{\Phi,n}=(\widehat\omega_n-\widehat\omega_{n-1})/(\Phi_n-\Phi_{n-1})$`，磁通偏置按阻尼、限步长的割线/Newton 型规则更新。实验中的电压命令由 Supplemental Material 给出的 `$V_{\mathrm b}\mapsto\Phi$` 换算实现；代码中的 `V` 字段仍按磁通偏置解释。

普通 secant、bisection 以及三算法的完整比较建议放 Supplemental Material。它们用于说明控制器选择的取舍，不占据主文叙事中心。

#### III.4.3 与 Ramsey 和混合策略的关系

- **[暂定]** 双扫 Ramsey 承担 ACQUIRE/REACQUIRE 和独立核验；固定延迟双正交 I/Q Ramsey 是 TRACK 阶段的正文强基线。后者尚需代码实现，见 II.B。
- **[证据现状]** “短脉冲粗调 + Ramsey 精调”的预设混合策略没有表现出成本优势。它可作为负面或限制性对照，不定义为本文主方法。
- **[已确认]** coarse-to-fine `transient -> Ramsey` 负面结果主要进入 Supplemental Material，正文只概述其限制性结论。

### III.5 III.D Maintaining the local operating range

#### III.5.1 本节任务

解释局部短脉冲估计器如何参与跨越更大总调节范围的多轮标定。

建议的论证顺序：

1. 固定驱动时，器件频率移动会使失谐增大并最终离开局部可逆区；
2. 每轮更新驱动参考，使下一次短脉冲探针重新靠近零失谐；
3. 参考更新只维持测量有效性，不直接完成器件偏置标定；
4. 用 fixed-drive 对照检验该机制是否必要。

#### III.5.2 驱动更新策略

**[已确认]** 主文采用 model-free secant sensitivity 型 drive tracking，利用已经获得的频率估计和磁通变化 `$\Phi_{n+1}-\Phi_n$` 更新下一轮驱动参考。预测式使用 `$\widehat S_{\Phi,n}$`；依赖器件色散模型的 model-based tracking 放 Supplemental Material。

逐状态驱动策略 **[暂定]** 为：ACQUIRE/REACQUIRE 使用双扫 Ramsey 建立绝对频率参考；TRACK 用割线灵敏度预测下一工作点的驱动参考；LOCK 将驱动参考置于目标频率。LOCK 只承担达到目标后的连续确认和短时保持，长期频率稳定、动态带宽、PSD 和 Allan deviation 放 Discussion/Outlook。

**[证据现状]** 当前 F4 支持 `drive_policy="track"` 与 `g3_source="kernel_full"` 的组合。fixed drive 或与 tracking 不相容的 `fit` 路径不能自动得到同样结论。正式表述应把该组合写成当前协议的成立条件，而不是短脉冲方法的普遍性质。

**[平台查证记录]** F4 的冻结设置字段为 `V_seed=0`，边界 `(-0.05, 0.01) Phi_0`、首个盲探测步 `0.01 Phi_0`、阻尼 `0.8`、最大步长 `0.02 Phi_0`、目标频移 `-20 MHz`。其中 `V_seed` 和历史 `V` 实际表示相对磁通偏置 `$\delta\Phi/\Phi_0$`，不是伏特。该结果验证的是获得可用初始频率种子后的局部 tracking 阶段，不验证完整的 capture/track/lock/reacquire 状态机；第三章讨论时不得扩大证据范围。

**[暂定]** fixed-drive 作为正文机制对照，用于检验 drive tracking 的必要性；最终位置随图数和实验结果确认。当前 F4 未设置 `linear_range/rho`，不得写成已由 range guard 自动拦截越界。

### III.6 可检验预测与 Section IV 接口

第三章结尾不提前报告数值优势，而是把方法推导收束为三个由 Section IV 检验的有限预测：

1. 在 independently validated local operating interval 内，保留 `$G_3$` 的三阶估计器相对于一阶估计器应降低由响应非线性造成的模型偏差；该预测不延伸到 fold 外或全局测频。
2. 在相邻点割线预测有效时，drive-reference tracking 相对于 fixed drive 应减小下一轮测量失谐，并延长短脉冲估计器停留在局部有效区的能力；该预测不等同于长期频率稳定。
3. 当分支连续性、置信度或局部斜率条件失效时，局部反演不应外推，协议必须返回宽范围 reacquisition。

Section IV 分别用 estimator-error、tracking/fixed-drive 和完整资源核算检验上述预测。第三章只提供机制、公式、成立条件和失败边界；定量性能与统计显著性留给 Results。

### III.7 主文与 Supplemental Material 分工

主文建议保留：

- 差分观测量及局部奇函数结构；
- 一阶灵敏度 `$G_1$` 的操作性定义；
- 三阶系数 `$G_3$`、实际采用的三次反演及其 fold；
- 实际采用的反演式与可信区；
- 主偏置更新和驱动参考更新的职责；
- 标定流程成立所需的条件和失败边界。

Supplemental Material 建议放置：

- 偏置电压 `$V_{\mathrm b}$` 到 calibrated flux-bias coordinate `$\Phi$` 的换算，包括电压零点、磁通周期、拟合或标定方法、换算不确定度与漂移；
- 原始电压坐标与 `$\Phi/\Phi_0$` 坐标的对应验证，以及 DAC 分辨率、可用范围和限幅规则；
- 电压--磁通标定在 cold-start 与复用场景下的成本计入和摊销口径；
- 一般 Volterra/核函数理论的详细回顾和一阶核完整推导；
- 高阶核的嵌套对易子、时间排序和完整三时间积分；
- full three-time 与 time-diagonal 收缩的数学区别；在没有直接结果对照时不把该区别提升为主文性能主张；
- `fit` 与 `kernel_full` 两条 `$G_3$` 路径的实现和交叉验证；
- Virtual-Z、脉冲离散化和多能级设置；
- secant、bisection、阻尼割线/Newton 型更新的完整比较；
- model-based drive tracking、参数敏感性和额外失败案例；
- 成本计数、缓存和预标定成本的摊销规则。
- 离线 Virtual-Z 三时间采样的网格、有限差分、插值、复用和重标定规则。
- 完整状态机实现、状态转换阈值敏感性以及额外失效/重捕获案例。

### III.8 完成检查

- [ ] 本章先建立短脉冲响应中的信息，再讨论反演和闭环使用。
- [ ] Herb 2024 的引用没有越过一阶线性感知核边界。
- [ ] `$G_3$` 的来源、必要性和主文位置有明确依据。
- [ ] 专项核对全文 calibration/validation 数据术语：首次定义、采集时机、数据独立性和各自职责明确，并与闭环 independent verification 严格区分。
- [ ] 局部有效性没有被写成全局测频能力。
- [ ] 偏置更新和驱动参考更新定义清楚且没有混为一谈。
- [ ] “参与标定”的具体位置能在流程中被指出。
- [ ] Ramsey、fixed drive 和 hybrid 的角色与当前证据一致。
- [ ] 算法选择没有取代短脉冲序列参与标定的基本叙事。
- [ ] 开头明确说明为何 local TRACK 是完整 workflow 中需要细化的关键环节。
- [ ] 结尾输出可由 Section IV 检验的有限预测，而不是提前宣称性能优势。

### III.9 写作时需要定案的事项

- 第三章最终英文名称；
- fixed-drive 对照最终留正文还是 Supplemental Material；
- Fig. 2 保留双面板，还是将响应压缩为估计误差图的 inset；
- 离线 Virtual-Z 的采样网格、插值变量、复用范围和重标定条件；
- ACQUIRE/REACQUIRE 的宽范围方法及逐状态驱动策略最终版本；
- validated local operating interval 和状态转换阈值的实验数值；
- 自动循环状态机在投稿前需要完成到何种实现与验证等级。

---

## Section IV. Results and comparison

### IV.1 本章任务

**[待讨论]** 用数据依次回答短脉冲方法是否准确、在哪个范围有效、能否参与完整标定，以及与基线相比需要多少资源。

### IV.2 建议的结果顺序

1. 局部信息是否存在：展示 `$p_{\mathrm d}(\Delta)$`、局部斜率和可逆区；
2. 反演是否准确：比较一阶与三次估计误差并显示折叠边界；
3. 闭环能否收敛：展示主偏置更新下的残差轨迹；
4. 驱动跟随是否必要：比较 tracking 与 fixed drive；
5. 不同策略承担什么角色：双扫 Ramsey 负责捕获，固定延迟双正交 I/Q Ramsey 是局部强基线，一阶短脉冲是内部基线，三阶短脉冲 + drive tracking 是主方法；
6. 资源结论是否公平：在共同置信度和残差阈值下分别报告 cold-start、amortized maintenance 和 local-tracking 成本，并说明求解调用数、墙钟时间和实验资源不是同一指标；
7. 负面结果或失败案例。

**[待讨论]** 该顺序需要根据最终主张和实验数据调整。

**[证据现状]** 当前结果支持把 `track + kernel_full` 作为数值主候选协议，把 fixed drive 作为暂定机制对照，把双扫 Ramsey-only 作为已有宽范围基线，把 coarse-to-fine hybrid 作为补充材料中的负面结果。固定延迟双正交 I/Q Ramsey 尚待实现，完成前不得写成已有比较结果。

### IV.3 图表与已有素材

- `F1_flux_curve/f_phi_curve_sqc.*`：频率-磁通估计范围。
- `F2_transient_accuracy/transient_accuracy_sqc.*`：局部估计精度与翻折边界。
- `F3_step_method/step_method_comparison_sqc.*`：更新方法比较。
- `F4_hybrid_convergence/`：收敛与成本。
- `F5_ramsey_vs_transient/ramsey_vs_transient_sqc.*`：精度-资源比较。

现有图只作为数据和设计素材；正式图序、指标和结论需重新讨论。

候选图表分工：

| 图 | 目的 | 候选内容 | 建议位置 |
|---|---|---|---|
| Fig. 1 | 定义物理问题和测量协议 | 频率-磁通关系、Ramsey/短脉冲时序、差分响应示意、最简标定流程 | Section II |
| Fig. 2 | 建立局部频率推断 | 原始差分响应、一阶/三次误差、可信窗口与翻折 | Section III 或 Results 起始 |
| Fig. 3 | 验证标定机制 | tracking/fixed-drive 对照（暂定）；不重复 Fig. 1(d) 的流程示意 | Section IV |
| Fig. 4 | 验证闭环表现 | 主协议、固定延迟 Ramsey 强基线、宽范围 Ramsey 与共同阈值成本 | Section IV |
| Fig. 5（可选） | 给出精度-成本边界 | Ramsey 与短脉冲策略的 Pareto 比较 | Section IV 或 Supplemental |

若需压缩 PRA 主文，不新增独立流程图；算法全比较、coarse-to-fine 负结果和扩展 Pareto 图优先移至 Supplemental Material。

### IV.4 完成检查

- [ ] 每个结论都能指向具体图、数据或统计量。
- [ ] 所有策略在共同误差定义和资源口径下比较。
- [ ] 数值求解成本没有被写成实验时间优势。
- [ ] 正面结果和失效范围同时呈现。
- [ ] tracking 的效果通过 fixed-drive 消融而不是仅靠单条成功轨迹说明。
- [ ] 当前协议的成立条件没有被写成短脉冲序列的普遍性质。

### IV.5 待讨论事项

- **[暂定]** 实验主资源指标采用 `$T_{\mathrm{hit}}(\epsilon_\omega)$`，即固定置信度下首次经独立核验达到阈值的真实墙钟时间；图表用 `$\epsilon_f=\epsilon_\omega/(2\pi)$` 报告 MHz 阈值。仍需确定置信度、核验频率和计时边界。
- **[已确认]** 同时报告三种成本：`$T_{\mathrm{cold}}=T_{\mathrm{setup}}+T_{\mathrm{run}}$`、`$T_{\mathrm{amortized}}=T_{\mathrm{run}}+T_{\mathrm{setup}}/N_{\mathrm{reuse}}$`，以及只隔离 TRACK 阶段的 local-tracking cost。三者不得互相替代。
- setup 至少包括脉冲/读出准备和离线 Virtual-Z 核标定；run 包括 acquire、track、lock/verify 和实际发生的 reacquire。
- 理论/数值与真实器件实验共同构成证据；仍需确定各结果面板的具体数据来源及 calibration/validation 实施方案。
- 求解器调用次数只作为数值诊断，不作为实验时间或主资源指标。
- coarse-to-fine 的负面结果主要进入 Supplemental Material。

---

## Section V. Discussion

### V.1 本章任务

**[待讨论]** 回到全过程基本叙事，解释短脉冲序列参与标定的成立条件、实际意义、局限和可推广边界。

### V.2 要回答的问题

- 当前证据实际证明了哪一种“可能”？
- 哪些条件来自短脉冲响应本身，哪些来自外部标定流程？
- 与传统 Ramsey 是替代、互补还是分层协作关系？
- 哪些结论只适用于当前频率标定案例？
- 实验噪声、泄漏、漂移和硬件带宽会改变什么？

### V.3 完成检查

- [ ] Discussion 解释结果而不重复 Results。
- [ ] 推广结论不超出当前证据。
- [ ] 明确区分数值可行性和实验优势。

### V.4 待讨论事项

- 是否讨论其他标定任务。
- response kernel 和控制策略在 Discussion 中各占多少权重。
- 长期频率稳定、闭环带宽、频率噪声 PSD 和 Allan deviation 只作为讨论或展望，不作为当前 LOCK 状态已经验证的主结果。

---

## Section VI. Conclusion

### VI.1 本章任务

**[待讨论]** 直接回答 Introduction 提出的研究问题，不加入新数据或新机制。

### VI.2 完成检查

- [ ] 首句回到短脉冲序列参与标定的可能。
- [ ] 总结证据、成立条件和限制。
- [ ] 不把单一案例写成普遍结论。
- [ ] 不把数值资源降低写成已验证的实验加速。

---

## Appendix and Supplemental Material

### M.1 本部分任务

主文后的 Appendix 保存复现实验和数值结果所必需、但会中断主体论证的技术步骤；独立 Supplemental Material 保存长推导、额外对照和扩展数据。Appendix 不是 Section III 所指的 proposed method。

### M.2 候选内容

- Transmon 参数和 Hamiltonian 细节；
- 偏置电压 `$V_{\mathrm b}$` 到 calibrated flux-bias coordinate `$\Phi$` 的换算：电压零点、磁通周期、拟合/标定方法、换算不确定度和时间漂移；
- Fig. 1(a) 实验点从原始 `$V_{\mathrm b}$` 到 `$\Phi/\Phi_0$` 的转换验证，以及 DAC 分辨率、执行范围和限幅；
- 电压--磁通换算的 cold-start 标定成本与跨轮次复用时的摊销规则；
- 脉冲包络、一般转角 `$\alpha$` 的定义与选择准则、采样步长和读出模型；
- response kernel 推导或数值标定过程；
- 一阶核的完整推导，高阶核的嵌套对易子、时间排序和三时间积分；
- `fit` 与 `kernel_full` 两条 `$G_3$` 路径的实现与交叉验证；
- 完整多能级、泄漏和噪声分析；
- 芯片照片、等效电路和低温测控链路；
- 控制器参数、额外基线和稳健性测试；
- secant、bisection、阻尼割线/Newton 型更新的完整比较；
- model-based drive tracking、参数敏感性和额外失败案例；
- 成本计数、缓存和预标定成本的摊销规则；
- 数据与代码可复现说明。

### M.3 完成检查

- [ ] 主文中的每个关键结果都能在方法或补充材料中找到复现依据。
- [ ] 补充材料承载细节，但不隐藏主张成立所必需的关键证据。

---

## Title and Abstract

### T.1 写作顺序

**[建议]** 在正文、图序和 Discussion 稳定后最后撰写。

### T.2 完成检查

- [ ] 标题准确反映核心问题，不把某个机制误写成全文中心。
- [ ] 摘要包含背景、问题、方法、主要证据、边界和意义。
- [ ] 摘要中的每个定量结果都能在正文中找到对应证据。

---

## D. 建议的实际写作顺序

章节在论文中的展示顺序与实际写作顺序可以不同。当前建议：

1. Section II：先完成系统、标定任务和 Fig. 1 定义；
2. Section III.A-III.B：写短脉冲响应、核函数的必要接口和局部反演；
3. Section IV：整理局部精度、闭环、消融与资源比较结果；
4. Section III.C-III.D：根据结果确定偏置更新与 drive tracking 的正文权重；
5. Section V：讨论成立条件与边界；
6. Section I：在主线稳定后写 Introduction；
7. Section VI：写 Conclusion；
8. Appendix / Supplemental；
9. 最后写 Title 和 Abstract。

## E. 全局待讨论事项

- 论文最终只讨论频率标定案例，还是在 Discussion 中谨慎延伸到其他标定任务？
- **[已确认]** 正式论文包含真实器件实验，并采用理论/数值与实验共同支撑的证据结构。
- **[待实施]** 各结果面板的理论/数值与实验分工，以及 calibration/validation 数据的采集、划分和参数冻结方案。
- **[暂定]** 主要实验资源指标采用 `$T_{\mathrm{hit}}(\epsilon_\omega)$`，图表可等价报告 `$T_{\mathrm{hit}}(\epsilon_f)$`；具体统计、独立核验与计时口径由 Section IV/Appendix 和实验方案冻结。
- **[暂定]** fixed-drive 对照最终留正文还是 Supplemental Material。
- **[暂定]** Fig. 2 双面板或响应 inset 的最终版式。
- **[暂定]** 离线 Virtual-Z 核采样、插值、复用与重标定方案。
- **[暂定]** ACQUIRE/REACQUIRE 方法和逐状态驱动策略。
- **[待实施]** 固定延迟双正交 I/Q Ramsey 强基线与自动循环状态机。
- **[待实验]** validated local operating interval、状态转换阈值和三种成本口径的真实数据。

---

## F. Section II--III 连续候选初稿（2026-08-06）

> **文件性质**：以下内容是依据当前已确认写作决策形成的英文参考初稿，附在计划中供用户逐段审核；它不是正式论文正文，不自动替换 paper/sections 下的既有文件。任何措辞、公式取舍和章节标题仍由用户最终决定。
>
> **当前引用边界**：Koch、Krantz、Blais、Herb 和 Vepsäläinen 已在本地工作区核验；Paladino、Ramsey、Kubo 与 McKay 已通过 Crossref/Semantic Scholar 核验 DOI 和内容范围；Mukamel 专著已核验书目信息但尚未本地入库。Herb 2024 只支撑一阶、弱信号感知核及相关短序列分析，Kubo/Mukamel 只提供一般响应理论背景，均不作为完整三阶三时间核的现成来源。第三阶推导由本文 Supplemental Material 自洽承担。

> **同步提示（2026-08-07）**：以下候选稿的 Section III.A 尚未按 B.5 收缩，仍保留一般 Volterra 展开、`$k_m^{(\mathrm d)}$`、`$R_5$` 和 time-diagonal 强调。下一轮修改以 B.5 与 III.2 的最新安排为准，不应把这些旧表述视为已确认正文。

~~~latex
\section{System and calibration task}
\label{sec:system}

\subsection{Flux-tunable transmon and calibration task}
\label{sec:system-objective}

We consider a flux-tunable transmon controlled through a calibrated flux-bias
coordinate $\Phi$. Experimentally, $\Phi$ is controlled through a bias voltage
$V_{\mathrm b}$ generated by room-temperature electronics and expressed in
flux units using a separately calibrated voltage-to-flux relation. Details of
this conversion are provided in Appendix~\ref{app:methods}. Throughout the main
text, $\Phi$ denotes this calibrated flux-equivalent bias rather than a direct
measurement of the on-chip flux.

In the transmon regime, the flux dispersion of the qubit angular transition
frequency is approximately given by
%
\begin{equation}
\hbar\omega_{01}(\Phi)
\simeq
\sqrt{8E_CE_J(\Phi)}-E_C,
\qquad
E_J(\Phi)
=
E_{J,\Sigma}
\sqrt{
\cos^2\left(\pi\Phi/\Phi_0\right)
+
d^2\sin^2\left(\pi\Phi/\Phi_0\right)
},
\label{eq:transmon-dispersion}
\end{equation}
%
where $E_C$ is the charging energy, $E_{J,\Sigma}$ is the total Josephson
energy, $d$ characterizes the junction asymmetry, and $\Phi_0$ is the flux
quantum~\cite{koch2007charge,krantz2019quantum}. The numerical calculations
use the symmetric limit $d=0$.

The calibration objective is to locate a flux bias $\Phi_\star$ at which the
qubit angular transition frequency matches a prescribed target
$\omega_{\mathrm{tar}}$ within a tolerance $\epsilon_\omega$,
%
\begin{equation}
\left|
\omega_{01}(\Phi_\star)-\omega_{\mathrm{tar}}
\right|
\leq
\epsilon_\omega,
\label{eq:calibration-success}
\end{equation}
%
We consider this task in situations where an up-to-date quantitative relation
$\omega_{01}(\Phi)$ is unavailable. Low-frequency noise can produce slow,
temporally correlated shifts of the qubit frequency, so the frequency
predicted by a previously characterized map may no longer be accurate at the
target tolerance~\cite{paladino2014noise,vepsalainen2022improving}. The prior
information available to the closed-loop calibration is restricted to a safe
flux interval, the local monotonic direction, and bounds on the flux-bias and
drive-frequency updates. Here calibration denotes locating the target working
point, whereas stabilization denotes suppressing subsequent drift after that
point has been established and is not the primary task considered here.

Figure~\ref{fig:system-protocol}(a) represents this task on the flux
dispersion by marking the initial bias, the target frequency, and the local
operating region. Its horizontal coordinate is the calibrated flux bias
$\Phi/\Phi_0$; experimental points are converted from the applied
$V_{\mathrm b}$ using the conversion defined in Appendix~\ref{app:methods},
with the corresponding raw-voltage representation provided in the
Supplemental Material.

\begin{figure*}[t]
    \centering
    % PLACEHOLDER: replace with the final four-panel system/protocol figure.
    \fbox{%
        \parbox[c][42mm][c]{0.96\textwidth}{%
            \centering
            \textbf{Placeholder for Fig.~1}\\[1ex]
            (a) flux dispersion and calibration target;\quad
            (b) probe timing;\\[0.5ex]
            (c) bipartite differential observable;\quad
            (d) calibration protocol
        }%
    }
    \caption{Schematic overview of the physical system, measurement probes,
    and calibration protocol. (a) Qubit flux dispersion in the calibrated
    coordinate $\Phi/\Phi_0$, with the initial bias, target frequency, and
    local operating region indicated. (b) Timing diagrams for scanned Ramsey,
    fixed-delay Ramsey, and the bipartite short-pulse probe, with the tunable
    short-pulse rotation angle $\alpha$ indicated. (c) Construction
    of the differential observable $p_{\mathrm d}$ from the phase-cycled
    $+X$ and $-X$ branches, together with its locally monotonic response about
    resonance. (d) The \textsc{Acquire--Track--Lock--Reacquire} protocol,
    distinguishing flux-bias updates from drive-reference updates.}
    \label{fig:system-protocol}
\end{figure*}

\subsection{Ramsey and bipartite short-pulse probes}
\label{sec:system-probes}

At each calibration iteration, the qubit is interrogated at a fixed flux bias
by a control sequence followed by population readout. We refer to this
pulse-and-readout operation as a measurement probe. The probe maps the
detuning between the qubit at the current flux bias and the applied microwave
drive onto measured population data. We consider two probe families: Ramsey
probes and bipartite short-pulse probes. The driven evolution in both families
is described by an effective two-level Hamiltonian. In a frame rotating at
the microwave drive frequency $\omega_d$ and under the rotating-wave
approximation~\cite{krantz2019quantum,blais2021circuit},
%
\begin{equation}
\frac{H(t;\Phi)}{\hbar}
=
\frac{\Delta(\Phi)}{2}\sigma_z
+
\frac{\Omega(t)}{2}
\left[
\cos\varphi(t)\sigma_x
+
\sin\varphi(t)\sigma_y
\right],
\qquad
\Delta(\Phi)
=
\omega_{01}(\Phi)-\omega_d,
\label{eq:rotating-frame-hamiltonian}
\end{equation}
%
where $\Omega(t)$ and $\varphi(t)$ are the drive envelope and phase. Population
readout does not measure detuning directly: the control sequence encodes
$\Delta$ in the final-state populations, from which the estimator obtains
$\widehat\Delta$ and hence
$\widehat\omega_{01}=\omega_d+\widehat\Delta$. Equations use angular
frequencies; values plotted in hertz are divided by $2\pi$.

The timing structures of the three probes are compared in
Fig.~\ref{fig:system-protocol}(b). Ramsey sensing includes an explicit
free-evolution interval, whereas the bipartite short-pulse probe encodes
detuning during the finite control pulses themselves.

For an ideal Ramsey sequence with free-evolution time $\tau$ and analysis
phase $\varphi$, the final population can be written
as~\cite{ramsey1950molecular,vepsalainen2022improving}
%
\begin{equation}
p_R(\Delta;\tau,\varphi)
=
\frac{1}{2}
\left[
1+C\cos(\Delta\tau+\varphi)
\right],
\label{eq:ramsey-response}
\end{equation}
%
where $C$ is the Ramsey contrast. Scanning $\tau$ provides a comparatively
wide capture range and is used for acquisition, reacquisition, and independent
verification. Once a local frequency reference is available, a fixed-delay
Ramsey measurement can instead infer the phase at a selected delay $\tau_0$.
Two orthogonal analysis phases provide the quadratures needed to recover the
local phase, but the inferred detuning remains restricted by phase wrapping.
We therefore treat scanned Ramsey and fixed-delay Ramsey as distinct
baselines: the former is a wide-capture method, whereas the latter is a local
frequency-estimation method~\cite{vepsalainen2022improving}. Scan processing,
phase reconstruction, and the selection of $\tau_0$ are specified in
Appendix~\ref{app:methods}.

The short-pulse probe removes the free-evolution interval altogether.
Following the bipartite zero-delay Ramsey sequence analyzed by Herb and Degen
\cite{herb2024quantum}, we use a common preparation pulse of rotation angle
$\alpha$ followed immediately by one of two phase-cycled analysis pulses of
the same angle,
%
\begin{equation}
R_y(\alpha)
\longrightarrow
R_{+x}(\alpha)
\quad\text{or}\quad
R_y(\alpha)
\longrightarrow
R_{-x}(\alpha).
\label{eq:bipartite-sequences}
\end{equation}
%
Here $\alpha=\int\Omega(t)\,\dd t$ for the corresponding calibrated pulse
envelope; $\alpha=\pi/2$ is one reference setting rather than a defining
constraint. At fixed maximum Rabi rate and pulse shape, reducing $\alpha$
shortens the sequence and narrows its temporal response, but also reduces the
population contrast. The angle therefore tunes the tradeoff among temporal
resolution, detuning sensitivity, and acquisition cost. Its selection is
specified in Appendix~\ref{app:methods} and tested in Sec.~IV.

We refer to these as the $+X$ and $-X$ branches. If $p_{+X}$ and $p_{-X}$
are their final excited-state populations, the measured differential
observable is
%
\begin{equation}
p_{\mathrm d}
=
\frac{p_{+X}-p_{-X}}{2}.
\label{eq:differential-population}
\end{equation}
%
The individual branch populations, their difference, and the frequency
estimate obtained by inverting that difference are distinct quantities.
Figure~\ref{fig:system-protocol}(c) makes this measurement chain explicit and
indicates schematically the locally monotonic response used for inversion;
the quantitative response and its validated range are established separately.

Although Eq.~\eqref{eq:bipartite-sequences} contains no independently varied
free-evolution interval, detuning acts during the finite control pulses in
Eq.~\eqref{eq:rotating-frame-hamiltonian}. It therefore changes the final
branch populations and produces a differential response near resonance. A
necessary condition for local frequency identification is
%
\begin{equation}
\left.
\frac{\partial p_{\mathrm d}}{\partial\Delta}
\right|_{\Delta=0}
\neq 0.
\label{eq:local-identifiability}
\end{equation}
%
The experimental response is calibrated on one data set and evaluated on an
independent validation set. Zero-point offsets, branch imbalance, and other
nuisance parameters are frozen before the validation data are analyzed.

Equation~\eqref{eq:local-identifiability} establishes only local
identifiability. It does not imply that $p_{\mathrm d}$ provides an
unambiguous frequency over an unrestricted detuning range. In this work, a
short-pulse sequence is said to participate in calibration when its measured
differential response yields a validated local frequency estimate and that
estimate subsequently enters the flux-bias update. This definition does
not assume that the short-pulse probe replaces wide-range Ramsey acquisition
or that local calibration is equivalent to long-term stabilization.

\subsection{Calibration workflow and evaluation framework}
\label{sec:calibration-workflow}

At the functional level, a calibration cycle combines measurement,
estimation, and tuning. A measurement returns population data from a probe at
fixed $\Phi_n$; estimation maps those data to $\widehat\Delta_n$ and
$\widehat\omega_{01,n}$; and tuning changes the flux bias. Repeated cycles
locate the target working point defined by
Eq.~\eqref{eq:calibration-success}.

This modular separation is consistent with established closed-loop frequency
feedback and integrated calibration frameworks for superconducting
qubits~\cite{vepsalainen2022improving,wittler_integrated_2021}.

We use the term closed-loop calibration controller for the classical feedback
logic that, given the current frequency estimate and protocol state, selects
the next flux bias $\Phi_{n+1}$ and drive reference $\omega_{d,n+1}$. The
controller does not include the quantum probe, the frequency estimator, or the
independent verification measurement.

Figure~\ref{fig:system-protocol}(d) summarizes the complete
\textsc{Acquire--Track--Lock--Reacquire} workflow. Scanned Ramsey supplies a
wide-range frequency reference during \textsc{Acquire} and
\textsc{Reacquire} and provides independent verification. During
\textsc{Track}, the bipartite short-pulse probe supplies the population data
for the local estimator developed in Sec.~III; fixed-delay Ramsey serves as
the corresponding local baseline. \textsc{Lock} confirms arrival at the
target and provides only short-term holding. The diagram distinguishes the
flux-bias update, which moves the physical working point, from the
drive-reference update, which changes the reference of the next measurement.
The state-transition conditions and both update rules are developed in
Sec.~III rather than specified at this architectural level.

Successful calibration is not declared solely from the estimator used in the
feedback loop. We verify Eq.~\eqref{eq:calibration-success} independently,
using the withheld transmon dispersion in simulation and an independent
frequency measurement in experiment. Experimental calibration time is itself
a relevant resource in superconducting
processors~\cite{rol2017restless,werninghaus_high-speed_2021}. We therefore denote by
$T_{\mathrm{hit}}(\epsilon_\omega)$ the elapsed experimental time at which
this verification first confirms, at a prescribed confidence, that the error
is within $\epsilon_\omega$. The verification schedule, timing convention,
and cold-start and amortized cost accounting are specified in
Appendix~\ref{app:methods} and evaluated in Sec.~IV.

This workflow isolates local tracking as the stage in which short-pulse
measurements may participate in calibration between wide-range acquisitions.
Doing so requires both a reliable inverse of their nonlinear local response
and a way to keep successive measurements within the validated operating
range as the flux bias moves. Section~III develops these two coupled elements.


\section{Short-pulse frequency inference and calibration}
\label{sec:inference-calibration}

We now instantiate the local \textsc{Track} branch of
Fig.~\ref{fig:system-protocol}(d). Two coupled problems determine whether the
short-pulse probe can serve this role. First, its differential response is
only locally invertible and becomes nonlinear away from resonance, so a
frequency estimator must include both a nonlinear correction and an explicit
validity boundary. Second, changing the flux bias moves the qubit frequency;
with a fixed drive, successive measurements can therefore leave that local
range even when the device moves toward the target.

The following sections address these problems in sequence. We first derive
the local differential response and its inverse, then use the resulting
frequency estimate in a bounded flux-bias update, and finally update the drive
reference to maintain the measurement range. This construction instantiates
the measurement, estimation, and tuning blocks of Sec.~II.C without claiming
that the short-pulse probe replaces wide-range acquisition or independent
verification.

\subsection{Short-pulse differential response}
\label{sec:differential-response}

The response depends on the selected rotation angle $\alpha$. In general the
kernels, integrated coefficients, and validated interval may be written as
$k_m^{(\mathrm d)}(\ldots;\alpha)$, $G_m(\alpha)$, and
$\Delta_{\mathrm{val}}(\alpha)$. Once $\alpha$ has been selected on the
calibration set and frozen for validation, we suppress this parameter in the
notation below.

To quantify the frequency information contained in
Eq.~\eqref{eq:differential-population}, we regard each branch population as a
functional of a time-dependent detuning $\Delta(t)$. Expanding around the
resonant control evolution gives the response-function
expansion~\cite{kubo1957statistical,mukamel1995principles}
%
\begin{equation}
p_{\mathrm d}[\Delta]
=
p_{\mathrm d}^{(0)}
+
\sum_{m=1}^{\infty}
\frac{1}{m!}
\int_0^T\!\dd t_1\cdots\int_0^T\!\dd t_m\,
k_m^{(\mathrm d)}(t_1,\ldots,t_m)
\prod_{j=1}^{m}\Delta(t_j),
\label{eq:differential-volterra}
\end{equation}
%
where $k_m^{(\mathrm d)}$ is the symmetrized differential response kernel and
$T$ is the total sequence duration. The first-order term is the linear
sensing response considered in Ref.~\cite{herb2024quantum}. The derivation of
the higher-order kernels, including time ordering and nested commutators, is
given in the Supplemental Material.

For frequency estimation at a fixed flux bias, the detuning is constant
during a single short sequence, $\Delta(t)=\Delta$. Defining
%
\begin{equation}
G_m
=
\int_0^T\!\dd t_1\cdots\int_0^T\!\dd t_m\,
k_m^{(\mathrm d)}(t_1,\ldots,t_m)
\label{eq:integrated-response-coefficients}
\end{equation}
%
reduces Eq.~\eqref{eq:differential-volterra} to an ordinary local expansion
in $\Delta$. For the ideal phase-cycled sequence,
$p_{\mathrm d}(-\Delta)=-p_{\mathrm d}(\Delta)$, so the zero-order and even
terms cancel. Truncating after third order gives
%
\begin{equation}
p_{\mathrm d}(\Delta)
=
G_1\Delta
+
\frac{G_3}{3!}\Delta^3
+
R_5(\Delta),
\qquad
R_5(\Delta)=O(\Delta^5).
\label{eq:short-pulse-cubic-response}
\end{equation}
%
The required coefficients are
%
\begin{equation}
G_1
=
\int_0^T k_1^{(\mathrm d)}(t)\,\dd t,
\qquad
G_3
=
\iiint_{[0,T]^3}
k_3^{(\mathrm d)}(t_1,t_2,t_3)
\,\dd t_1\dd t_2\dd t_3.
\label{eq:g1-g3}
\end{equation}
%
The cubic coefficient is a full three-time integral. Restricting it to
coincident times, for example by integrating only
$k_3^{(\mathrm d)}(t,t,t)$, discards unequal-time contributions and does not
give the cubic response to a constant detuning.

We consider two routes for obtaining these coefficients. In the model route,
the kernels are calculated from the control Hamiltonian and the resonant
propagator. In the measurement route, calibrated Virtual-$Z$
perturbations~\cite{mckay2017efficient} are inserted at selected pulse times
and finite differences of the measured branch populations sample the required
response derivatives. The latter
permits $G_1$ and $G_3$ to be calibrated without supplying the controller
with a quantitative $\omega_{01}(\Phi)$ model. The offline calibration and its
cost accounting are specified in Appendix~\ref{app:methods}; the finite-
difference construction and three-time sampling scheme are given in the
Supplemental Material. The calibrated coefficients may be reused or
interpolated only over a range established on the calibration set, with the
reuse range and recalibration condition frozen before validation. A complete
three-dimensional measurement is not repeated at every feedback iteration.

Experimental imperfections can break the ideal odd symmetry. We therefore
define the corrected observable
%
\begin{equation}
\widetilde p_{\mathrm d}
=
\mathcal C\!\left(
p_{+X}^{\mathrm{meas}},p_{-X}^{\mathrm{meas}};
\bm\theta_{\mathrm{nuis}}
\right),
\label{eq:corrected-differential-response}
\end{equation}
%
where $\mathcal C$ reduces to $(p_{+X}-p_{-X})/2$ in the ideal case. The
nuisance parameters $\bm\theta_{\mathrm{nuis}}$, including the zero-point
offset, relative branch contrast, and readout corrections, are determined on
the calibration set and frozen before validation. A statistically significant
even residual is treated as evidence of model mismatch rather than being
absorbed automatically into the ideal odd response model.

\subsection{Local frequency estimator}
\label{sec:local-estimator}

The first-order frequency estimate follows from the local slope:
%
\begin{equation}
\widehat\Delta_1
=
\frac{\widetilde p_{\mathrm d}}{G_1},
\qquad
\widehat\omega_{01,1}
=
\omega_d+\widehat\Delta_1.
\label{eq:linear-frequency-estimator}
\end{equation}
%
This estimator provides a transparent baseline but neglects the cubic and
higher-order terms in Eq.~\eqref{eq:short-pulse-cubic-response}. Its
systematic error therefore grows as the probe moves away from resonance.

The main local estimator retains the third-order term and solves
%
\begin{equation}
G_1\widehat\Delta_3
+
\frac{G_3}{6}\widehat\Delta_3^3
=
\widetilde p_{\mathrm d}.
\label{eq:cubic-frequency-estimator}
\end{equation}
%
Newton iteration is initialized with $\widehat\Delta_1$, and the accepted
solution is the branch connected continuously to $\Delta=0$. When
$G_1G_3<0$, the cubic response has stationary points at
%
\begin{equation}
\Delta_{\mathrm{fold}}
=
\sqrt{-\frac{2G_1}{G_3}}.
\label{eq:response-fold}
\end{equation}
%
Beyond this fold, multiple detunings can produce the same differential
population. The implementation rejects a root that leaves the connected local
branch. A degenerate cubic coefficient may justify the linear estimator only
when that estimator has its own validated interval; otherwise a failed root or
branch test rejects the local estimate and triggers reacquisition. This guard
prevents selection of a disconnected root; it does not restore global
identifiability.

We distinguish the response fold from the validated local operating interval,
%
\begin{equation}
|\Delta|
\leq
\Delta_{\mathrm{val}}
,\qquad
\Delta_{\mathrm{val}}<\Delta_{\mathrm{fold}}
\quad (G_1G_3<0).
\label{eq:validated-interval}
\end{equation}
%
The fold is a mathematical boundary of the truncated cubic response, whereas
$\Delta_{\mathrm{val}}$ is determined from independent validation data. When
the cubic model has no real stationary point, the second condition in
Eq.~\eqref{eq:validated-interval} is absent, but the validation boundary still
applies. It must satisfy prescribed limits on truncation error, estimator
bias, root-mean-square error, and confidence under the experimental noise
model. Measurements outside this interval trigger wide-range reacquisition
rather than extrapolation of the local inverse.

\subsection{Calibration protocol}
\label{sec:calibration-protocol}

The state roles and data flow of the full workflow were defined in Sec.~II.C.
Here we formalize the interfaces required to enter, execute, and leave its
local \textsc{Track} branch. In \textsc{Acquire}, a scanned Ramsey measurement
supplies an absolute frequency estimate and its uncertainty. The protocol
enters \textsc{Track} only when the inferred
measurement detuning is compatible with the validated interval, for example
when
%
\begin{equation}
|\widehat\Delta_n|
+
z_{1-\beta}\sigma_{\Delta,n}
\leq
\Delta_{\mathrm{val}},
\label{eq:track-entry-condition}
\end{equation}
%
where $\sigma_{\Delta,n}$ is the uncertainty of the detuning estimate and
$z_{1-\beta}$ is a confidence coefficient fixed by the experimental protocol.
In \textsc{Track}, the short-pulse estimator supplies the frequency estimate
used by the bias controller. The protocol enters \textsc{Lock} after the
target condition is satisfied for $N_{\mathrm{lock}}$ consecutive
verification measurements. A range violation, a discontinuous root,
insufficient confidence, an anomalous local slope, or repeated measurement
failure sends \textsc{Track} or \textsc{Lock} to \textsc{Reacquire}. A new
wide-range estimate then permits a return to \textsc{Track}.

The estimated angular-frequency residual is
%
\begin{equation}
\widehat r_n
=
\widehat\omega_{01,n}-\omega_{\mathrm{tar}}.
\label{eq:estimated-feedback-residual}
\end{equation}
%
At the seed bias $\Phi_0^{(\mathrm{seed})}$, the controller first performs one frequency
measurement. Because no quantitative local slope is yet available, the
second point is generated by a bounded one-sided exploratory step. If
$\chi_\Phi$ denotes the known local monotonic direction,
%
\begin{equation}
\Phi_1
=
\operatorname{clip}_{[\Phi_{\min},\Phi_{\max}]}
\left[
\Phi_0^{(\mathrm{seed})}
-
\chi_\Phi\operatorname{sgn}(\widehat r_0)
\delta\Phi_{\mathrm{blind}}
\right].
\label{eq:blind-step}
\end{equation}
%
This initialization consists of one measured seed followed by a blind step,
not a symmetric two-point bootstrap.

After two neighboring measurements, the local composite sensitivity is
estimated by
%
\begin{equation}
\widehat S_{\Phi,n}
=
\frac{
\widehat\omega_{01,n}-\widehat\omega_{01,n-1}
}{
\Phi_n-\Phi_{n-1}
}.
\label{eq:secant-sensitivity}
\end{equation}
%
The controller estimates this slope directly from neighboring measurements
and is not supplied with the quantitative dispersion or its derivative. The
next flux bias follows a damped, step-limited secant--Newton update,
%
\begin{equation}
\Phi_{n+1}
=
\operatorname{clip}_{[\Phi_{\min},\Phi_{\max}]}
\left[
\Phi_n
+
\operatorname{clip}_{[-\delta\Phi_{\max},\,\delta\Phi_{\max}]}
\left(
-\eta\frac{\widehat r_n}{\widehat S_{\Phi,n}}
\right)
\right],
\label{eq:flux-bias-update}
\end{equation}
%
where $\eta$ is the damping factor. Because its local sensitivity is inferred
from two measured iterates, Eq.~\eqref{eq:flux-bias-update} is a damped
secant--Newton step rather than gradient descent.

The flux-bias update changes the physical working point and reduces the
residual relative to the fixed target. It must be
distinguished from the drive-reference update, which changes only the
coordinate system of the next measurement.

\subsection{Maintaining the local operating range}
\label{sec:drive-tracking}

If the microwave drive is held fixed while the flux bias moves the qubit,
the next measurement detuning is
%
\begin{equation}
\Delta_{n+1}
=
\omega_{01}(\Phi_{n+1})-\omega_{d,n}.
\label{eq:fixed-drive-detuning}
\end{equation}
%
A calibration trajectory can then leave the validated interval even when
each bias update moves the device toward the target. The failure occurs
because a local measurement is being forced to span the full frequency
displacement.

During \textsc{Track}, we instead predict the transition frequency at the next
flux bias using the measured secant sensitivity and set
%
\begin{equation}
\omega_{d,n+1}
=
\widehat\omega_{01,n}
+
\widehat S_{\Phi,n}
\left(
\Phi_{n+1}-\Phi_n
\right).
\label{eq:drive-tracking}
\end{equation}
%
This reference update uses only estimates already produced by the loop and
does not evaluate Eq.~\eqref{eq:transmon-dispersion}. If the
local prediction is accurate, the next short-pulse measurement is performed
near zero detuning even though the absolute qubit frequency continues to move
toward $\omega_{\mathrm{tar}}$.

Equations~\eqref{eq:flux-bias-update} and \eqref{eq:drive-tracking} serve
different objectives. The former changes the physical working point; the
latter changes the measurement reference. Updating $\omega_d$ does not
redefine the target or conceal the calibration residual, because
$\widehat\omega_{01,n}=\omega_{d,n}+\widehat\Delta_n$ is always compared with
the fixed $\omega_{\mathrm{tar}}$.

The drive-reference policy depends on the protocol state. \textsc{Acquire}
and \textsc{Reacquire} use scanned Ramsey to establish an absolute frequency
reference. \textsc{Track} uses Eq.~\eqref{eq:drive-tracking} to predict the
reference at the next command. In \textsc{Lock}, the drive is placed at the
target frequency and the residual is checked repeatedly before calibration
is declared complete. Here \textsc{Lock} denotes confirmation and short-term
holding around the reached target; it does not claim demonstrated long-term
stabilization, feedback bandwidth, noise suppression, or improved Allan
deviation.

A fixed-drive implementation provides the direct mechanism control for
Eq.~\eqref{eq:drive-tracking}: both protocols use the same short-pulse
estimator and flux-bias update, but only one updates the measurement
reference. A model-based reference computed from a calibrated
$\omega_{01}(\Phi)$ relation is a supplementary comparison rather than the main
protocol, because it introduces the quantitative prior excluded from the
calibration task.

The construction above yields three testable predictions. Over the portion of
the independently validated interval in which cubic truncation remains
accurate, retaining $G_3$ should reduce systematic nonlinear bias relative to
the first-order inverse. When the measured secant slope predicts the next
working point reliably, drive-reference tracking should keep the next
measurement closer to resonance than a fixed drive. When branch continuity,
confidence, or local-slope conditions fail, the local inverse should not be
extrapolated and the workflow should return to wide-range reacquisition.
Section~IV tests these predictions through the estimator-error,
tracked-versus-fixed-drive, and resource comparisons.

The complete four-state workflow remains the practical architecture for
combining wide-range acquisition, local short-pulse tracking, target
confirmation, and recovery. The technical development in this section and
the current numerical implementation concern the seeded local
\textsc{Track} stage. Automatic state transitions and successful
reacquisition under all failure modes require separate thresholds and
validation before the complete state machine can be described as
experimentally demonstrated.
~~~

### F.1 候选稿下一轮审核重点

- Section II.A 的联合 Transmon 色散公式是否保留 $E_J(\Phi)$ 的显式形式，或移入 Appendix；
- Virtual-Z measurement route 的正文表述是否与最终实验实施完全一致；
- `$\alpha$` 的选择准则及其精度--时间权衡是否获得充分的实验或数值证据；
- $T_{\mathrm{hit}}$、$\Delta_{\mathrm{val}}$、$z_{1-\beta}$ 和 $N_{\mathrm{lock}}$ 的实验定义冻结后，替换当前符号化叙述；
- 第三章最终标题以及 fixed-drive 对照的主文位置。
