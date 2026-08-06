# `result_sqc/` — sqc API 成果图产出

本目录用**新 sqc API** 产出成果图,按能力分为**三大块**:

```
result_sqc/
├── reconstruction/          # 【块 1】波形重建 — 复现 result/ 冻结基准,叠加 src 逐点对比
│   ├── different_signals/
│   ├── D1_waveform_reconstruction/  (图 D1:双峰/复杂波形 + Wiener 正则化)
│   ├── amplitude_scan/
│   ├── lambda_scan/
│   ├── basis_comparison/
│   ├── convergence/
│   ├── D2_resolution_sensitivity/   (图 D2:时间分辨率 + 灵敏度)
│   ├── D3_nonlinearity_boundary/    (图 D3:局部线性边界 + LM 对照)
│   └── workpoint_comparison/        (result/ 无对应,新增对照)
├── frequency_calibration/   # 【块 2】频率标定 — sqc 新能力,无 src 基准
│   ├── F1_flux_curve/           # 每张图独占一夹:脚本 + 数据 + 图
│   ├── F2_transient_accuracy/
│   ├── F3_step_method/
│   ├── F4_hybrid_convergence/
│   └── F5_ramsey_vs_transient/
└── predistortion/           # 【块 3】波形预失真 — sqc 新能力,无 src 基准
    ├── P1_distortion_gallery/
    ├── P2_correction/
    ├── P3_transfer_fit/
    ├── P4_improvement_vs_strength/
    └── P5_protocol_driven/
```

> **块 2/3 目录约定(v2)**:与块 1「按主题」分夹不同,块 2/3 每张图独占一个
> `F<n>_<slug>/` 或 `P<n>_<slug>/` 子夹,内含该图的**脚本 + 数据(.npz)+
> 图(.png/.pdf)**。脚本因此深两级(`sys.path` 上跳 `../..` 找 `_common`),
> 其 `SUBDIR` 指向对应子夹,产出直接落在图旁。

**块 1(波形重建)** 有冻结 `result/` 基准:
- `result/` 保持不动(R1 冻结),作为**基准 (ground truth)**。
- 每个脚本用 sqc API 重跑同一物理实验,并在产出中**叠加 src 基准曲线**做逐点对比,
  目标是**数值尽可能贴合原图**。
- `workpoint_comparison/` 是 `result/` 没有的额外对照(甜点 / 线性区 / 非线性区)。

**块 2/3(频率标定、波形预失真)** 是 sqc **新增能力,`src/` 里没有对应基准**:
- 只画 sqc 自身的物理正确性图(收敛曲线、校正前后对比),**不做 src 叠加**。
- 数据全部来自对应 workflow 的 `run()` 返回 dict(见各节表格“sqc 源”)。

> 状态标记: ✅ 骨架已搭 / 🔨 待实现数据+绘图 / ⏳ 待澄清源

## 通用配方 (src → sqc)

`result/` 的数据类脚本几乎都是同一套 src 流程，对应到 sqc：

| 步骤 | src (基准) | sqc (复现) |
|---|---|---|
| 量子比特 | `src.qubit.TransmonQubit(EC=0.2·2π, EJ=10·2π, T1=100e3, T2=50e3)` | `sqc.devices.transmon.TransmonQubit(同参)` |
| 工作点 | `qubit.optimal_work_point()` = 0.9553 | 同值(`_common.OPTIMAL_FLUX`) |
| 前向测量 | `src.protocal.Protocal(4).evolve()` | `sqc.experiments.transient.TransientSensingExperiment(...).run()` |
| 信号 | `src.signal.Signal(type=2/4/5/7)` | `sqc.control.flux_signal.FluxSignal(type=2/4/5/7)`(类型定义一致) |
| Wiener | `src.analysis.Analysis.wiener_deconvolution(λ=5)` | `TransientReconstruction(method="wiener", lambda_reg=5)` |
| LM | `Analysis.numerical_inverse(basis, n_basis, λ, max_iter)` | `TransientReconstruction(method="lm", control_pulse=exp.control_pulse, qubit=q, basis_type=..., n_basis=..., lambda_reg=..., max_iter=...)` → `(signal, history)` |

Wiener 通道已验证与 src 逐位一致(见 `tests/integration/test_transient_experiment.py::test_sqc_wiener_matches_src_wiener_exactly`)。

## 目录与预计产出

### 块 1 — 波形重建 `reconstruction/`(有 src 基准,叠加对比)

> 📄 **每张图的报告级图注**(展示内容/如何读/物理含义/关键数值)见
> [`reconstruction/FIGURES.md`](reconstruction/FIGURES.md)。

| 子目录 | 对应 `result/` | 预计产出文件 | src 源脚本 | 状态 |
|---|---|---|---|---|
| `reconstruction/different_signals/` | 同名 | `comparison_sqc.png/.pdf`(叠加src)、`reconstruction_only_sqc.png/.pdf`(仅sqc,原波形+重建)、`comparison_metrics_sqc.txt`、`check_{sine,step,double_peak,complex}_sqc.png`、`signal_{name}_sqc.npz` | `result/different_signals/generate_single_signal_data.py` | ✅ |
| `reconstruction/D1_waveform_reconstruction/` | **新增(图 D1)** | `d1_waveform_regularization_sqc.png/.pdf/.npz`、metadata、metrics；高斯 `π/2` 控制下的双峰/复杂波形两法对比 + 类阶跃 Wiener λ 扫描 | 独立重算量子动力学、响应核及 LM | ✅ |
| `reconstruction/amplitude_scan/` | `signal_amp/amplitude_scan/` | `amplitude_scan_results_sqc.png/.pdf`(仅sqc,报告候选)、`amplitude_scan_vs_src_sqc.png/.pdf`(叠加src)、`amplitude_{val}_sqc.npz` | `result/signal_amp/amplitude_scan/generate_amplitude_scan_data.py` | ✅ |
| `reconstruction/lambda_scan/` | 同名 | `{wiener,lm}_lambda_scan_sqc.png/.pdf`(仅sqc)、`{wiener,lm}_lambda_scan_vs_src_sqc.png/.pdf`(叠加src)、`lambda_scan_sqc.npz` | `result/lambda_scan/plot_lambda_scan.py` | ✅ |
| `reconstruction/basis_comparison/` | 同名 | `basis_comparison_sqc.png/.pdf`(仅sqc,src无npz)、`basis_comparison_sqc.npz` | `result/basis_comparison/generate_basis_comparison_data.py` | ✅ |
| `reconstruction/convergence/` | 同名 | `lm_convergence_sqc.png`、`lm_residual_heatmap_sqc.png`、`lm_coefficient_norm_sqc.png`、`lm_history_sqc.npz`(仅sqc) | `result/convergence/plot_lm_convergence.py` | ✅ |
| `reconstruction/workpoint_comparison/` | **新增(无 src 对应)** | `workpoint_comparison_sqc.png/.pdf`(仅sqc)、`workpoint_metrics_sqc.npz` | — | ✅ |
| `reconstruction/D2_resolution_sensitivity/` | **新增(无 src 对应,图 D2)** | `resolution_sensitivity_sqc.png/.pdf`(Wiener 主图)、`time_resolution_sqc.npz`、`sensitivity_sqc.npz`、metadata、metrics | — | ✅ |
| `reconstruction/D3_nonlinearity_boundary/` | **替代旧硬编码非线性图(图 D3)** | `nonlinearity_boundary_sqc.png/.pdf`(仅sqc,报告候选)、高斯控制幅值/工作点缓存、`nonlinearity_boundary_sqc.npz`、4 个 `lm_amplitude_*_sqc.npz`、metadata、metrics | 统一采用 10 ns 截断高斯 `pi/2` 控制独立重算；旧脚本仅作构图参考 | ✅ |

### 块 2 — 频率标定 `frequency_calibration/`(sqc 新能力,无 src 基准)

每张图一个脚本,共享 `SUBDIR="frequency_calibration"`。⭐ = 论文重点。

> 📄 图注说明见 [`frequency_calibration/FIGURES.md`](frequency_calibration/FIGURES.md)。
> F1 已把 Ramsey vs transient 对比并入(原 F5 动机);F5 脚本仍可独立跑权衡图。

| 图 | 子夹 / 脚本 | sqc 源 | 预计产出(均在子夹内) | 成本 | 状态 |
|---|---|---|---|---|---|
| F1 f(Φ) 标定曲线(Ramsey vs transient) | `F1_flux_curve/plot_flux_response_curve_sqc.py` | `FrequencyMeasurement` 两 method 逐点 `measure(flux=offset)` | `f_phi_curve_sqc.png/.pdf/.npz`(解析 f(Φ)+两法测点+残差) | 中 | ✅ |
| **F2 ⭐ transient 测频精度** | `F2_transient_accuracy/plot_transient_accuracy_sqc.py` | `FrequencyMeasurement(method="transient")` order=1/3,g3_source fit/kernel_full | `transient_accuracy_sqc.png/.pdf/.npz`(测频误差 vs 失谐,3 曲线:线性+两条三次;窗口内 MAE 1.27→0.94 MHz) | 高 | ✅ |
| F3 三步法收敛对比 | `F3_step_method/plot_step_method_comparison_sqc.py` | `SinglePointFrequencyCalibration` × secant/bisection/gradient | `step_method_comparison_sqc.png/.pdf/.npz`(\|残差\| log-y vs 迭代 + bisection 括号折半;secant=2/gradient=5/bisection=9 步) | 中高 | ✅ |
| **F4 ⭐ 闭环测频成本解剖**(提速来自 track+kernel_full,**不是**粗精分工) | `F4_hybrid_convergence/plot_hybrid_convergence_sqc.py` | `FrequencyCalibrationWorkflow` × 3 策略:hybrid sweet+fit / **transient track+kernel_full** / ramsey-only | `freq_calibration_convergence_sqc.png/.pdf`、`freq_calibration_cost_sqc.png/.pdf`(3 面板)、`.npz`(实测解算计数;track+kf **4.1× 更省且 3.5× 更准**,预设 hybrid 反而 0.8×) | 高 | ✅ |
| F5 ramsey vs transient 权衡(精度-成本帕累托) | `F5_ramsey_vs_transient/plot_ramsey_vs_transient_sqc.py` | `FrequencyMeasurement` 两 method × τ 网格密度,6 策略 | `ramsey_vs_transient_sqc.png/.pdf/.npz`(MAE vs 实测 mesolve 数/墙钟;**order=3 支配全部 Ramsey 单扫**) | 中 | ✅ |

### 块 3 — 波形预失真 `predistortion/`(sqc 新能力,无 src 基准)

每张图一个脚本,共享 `SUBDIR="predistortion"`。⭐ = 论文重点。

> 📄 图注说明见 [`predistortion/FIGURES.md`](predistortion/FIGURES.md)。
> **P1 坑位**:指数族 `frequency_response(ω)` 是物理 rad/ns,而 `FIRDistortion`/
> `IIRDistortion` 是归一化 dt=1(`distortion.py:597`/`:745`)——同轴混画必须传 `ω·dt`。
> **P2/P4 坑位(重要)**:当拟合模型类 == 真值模型类(如 single_exp 真值 + single_exp
> 拟合)时,`curve_fit` 精确回收参数 → 逆滤波器是**代数精确抵消** → improvement ≈ 3e12,
> 是同义反复。真实数字要靠**模型失配**(P2 的 mismatch 场景 ×108)。**P4 骨架现在扫的是
> matched single_exp,会得到一条 ~1e12 的水平线,必须改成失配才有趋势** —— 详见 P2 节
> 实现注记 4。**(已在 P4 中处理:改为 matched/mismatch 双场景扫描。)**
> **P4 结论(反直觉,注意别写错)**:improvement **随失真强度下降**,∝A^−1.06 ——
> uncorrected 误差线性于 A(A^1.0000),失配残差是二阶(A^2.06)。**预失真对弱失真
> 最有效**。另有两个硬约束:settling 10⁻³ 指标在 **A≈0.06–0.08 之间失效**;AWG 摆幅
> 在 A=0.5 时需 **2.82×** 理想值。
> **P5 坑位(重要)**:① 4 个协议里**只有 cryoscope** 经公开路径可用
> (delay_ramsey 衰减到 0 且 `t_max` 无效;transient 崩 `TypeError`;pi_pulse 值域
> 0..0.2)。② cryoscope 的 rmse 3.7e-3 **不是精度指标** —— 理想无失真线的对照实验给出
> **5.4e-3(更大)**,残差由协议伪影主导。**表征精度前必须先跑 null 对照**。
> ③ 可用区间双边受限:step_amplitude ≲0.05、窗口固定 20.5..90 ns
> (`trunc_list` 硬编码切片,与 `t_max` 无关且不被转发)。
> **P3 坑位**:①「标定窗长」在解析路径上是**死轴** —— 无噪声阶跃被任意有限弧唯一确定,
> `t_max/τ` 从 0.2 扫到 20 相对误差都 <1e-9,所以 P4 引入失配只剩「模型失配 / 测量噪声」
> 两条路。② cascade 的参考幅度是**部分分式留数** (0.0480, 0.0305),**不是**逐级幅度
> (0.05, 0.03),搞错会误判拟合偏了 4%。③ rms vs 拟合阶数 K 在 K>K_true 处**非单调**
> (局部极小),**不可**用 argmin 做阶数选择。

| 图 | 子夹 / 脚本 | sqc 源 | 预计产出(均在子夹内) | 成本 | 状态 |
|---|---|---|---|---|---|
| P1 失真模型图谱 | `P1_distortion_gallery/plot_distortion_gallery_sqc.py` | `sqc.hardware.distortion.*` 的 `step_response`/`frequency_response`(6 模型) | `distortion_gallery_sqc.png/.pdf/.npz`(2×2:阶跃全窗 + 尾偏差 log-y + Bode 幅度/相位;含直通平台 inset) | 低 | ✅ |
| P2 校正前后波形 | `P2_correction/plot_correction_before_after_sqc.py` | `PredistortionValidationWorkflow` × 2 场景(matched / mismatch) | `predistortion_correction_sqc.png/.pdf`(2×2:两场景波形+残差 log-y+AWG 实播)、`predistortion_metrics_sqc.npz` | 低中 | ✅ |
| P3 传递函数拟合质量 | `P3_transfer_fit/plot_transfer_fit_quality_sqc.py` | `WaveformCalibration(method="transfer_function")` × 3 真值 × K=1..5 | `transfer_fit_quality_sqc.png/.pdf/.npz`(2×2:测阶跃 vs 拟合 + 残差 log-y + rms vs 阶数 K + 参数回收误差) | 低 | ✅ |
| **P4 ⭐ 改善 vs 失真强度** | `P4_improvement_vs_strength/plot_improvement_vs_strength_sqc.py` | `PredistortionValidationWorkflow` 扫 A(12 点)× 2 场景 | `improvement_vs_strength_sqc.png/.pdf/.npz`(2×2:rmse 幂律 + **improvement ∝1/A** + settling 阈值 + AWG 动态范围) | 中 | ✅ |
| **P5 ⭐ 协议驱动测量** | `P5_protocol_driven/plot_protocol_driven_measurement_sqc.py` | `WaveformCalibration(measurement_protocol="cryoscope")` + 理想线对照 | `protocol_driven_measurement_sqc.png/.pdf/.npz`(2×2:量子 vs 解析阶跃 + **伪影对照实验** + 工作区间 + τ 扫描) | 高 | ✅ |

### 暂不纳入 / 待澄清

| `result/` 图 | 原因 |
|---|---|
| `ramsey/field_reconstruction*.png` | ⏳ 未定位到明确生成脚本(候选 `appendix/protocol_demo.py`),定位后再补 `ramsey/` 子目录 |
| `signal_amp/transmon_nonlinearity.png` | ✅ 已由 `reconstruction/D3_nonlinearity_boundary/nonlinearity_boundary_sqc.png` 替代；Wiener、工作点和 LM 数据均在统一高斯控制下独立重算，不采用原脚本硬编码 RMSE |
| `transient_error/*` | 🔶 原脚本为**独立解析 2 能级模型**(手写 2×2 expm)，不走 src/sqc API，属物理验证，不在本次 API 复现范围 |
| `flowchart/*.pdf` | ❌ 手绘/TikZ 示意图，非数据图 |

## 运行方式(实现后)

```bash
# 块 1 波形重建(全部 ✅;含 LM 的节较慢,均带 npz 缓存,重画秒级)
python result_sqc/reconstruction/different_signals/generate_signal_data_sqc.py all       # 数据+check图(含LM,~25min)
python result_sqc/reconstruction/different_signals/generate_signal_data_sqc.py sine --no-lm  # 只 wiener(快)
python result_sqc/reconstruction/different_signals/plot_signal_comparison_panel_sqc.py    # 2x4 vs-src 总览+metrics
python result_sqc/reconstruction/different_signals/plot_reconstruction_only_sqc.py        # sqc-only 报告候选
python result_sqc/reconstruction/D1_waveform_reconstruction/generate_d1_waveform_regularization_sqc.py # D1 高斯包络正式图(有缓存时秒级；--force 全量重算)
python result_sqc/reconstruction/workpoint_comparison/generate_workpoint_comparison_sqc.py   # 工作点对比(纯wiener,快)
python result_sqc/reconstruction/amplitude_scan/generate_amplitude_scan_sqc.py            # 幅度扫描(纯wiener,快)
python result_sqc/reconstruction/lambda_scan/plot_lambda_scan_sqc.py                      # λ 扫描(含LM~24min;缓存后秒出)
python result_sqc/reconstruction/basis_comparison/generate_basis_comparison_sqc.py        # 基函数对比(LMx3)
python result_sqc/reconstruction/convergence/plot_lm_convergence_sqc.py                   # LM 收敛(含LM~6min)
# 块 1 图 D2 时间分辨率 + 灵敏度(数据/绘图分离;LM 实测 ~300s/次,是唯一成本)
python result_sqc/reconstruction/D2_resolution_sensitivity/generate_time_resolution_sqc.py --quick --no-lm  # 冒烟(~30s)
python result_sqc/reconstruction/D2_resolution_sensitivity/generate_sensitivity_sqc.py --quick --no-lm      # 冒烟(~40s)
python result_sqc/reconstruction/D2_resolution_sensitivity/generate_time_resolution_sqc.py    # D2-A 正式(16 间距 + 7 LM,~40min)
python result_sqc/reconstruction/D2_resolution_sensitivity/generate_sensitivity_sqc.py        # D2-B 正式(12 幅度×3 shot×32 种子 + 18 LM,~95min)
python result_sqc/reconstruction/D2_resolution_sensitivity/plot_resolution_sensitivity_sqc.py # D2 2×2 报告图(只读缓存,秒级)
# 块 1 图 D3 非线性边界(8 个 Wiener 幅值点 + 29 个工作点 + 4 个 LM 点；缓存后秒级)
python result_sqc/reconstruction/D3_nonlinearity_boundary/generate_nonlinearity_boundary_sqc.py
# 块 2 频率标定(每图独占一子夹,⭐重点)
python result_sqc/frequency_calibration/F1_flux_curve/plot_flux_response_curve_sqc.py        # F1
python result_sqc/frequency_calibration/F2_transient_accuracy/plot_transient_accuracy_sqc.py # F2 ⭐
python result_sqc/frequency_calibration/F4_hybrid_convergence/plot_hybrid_convergence_sqc.py # F4 ⭐
# 块 3 波形预失真(每图独占一子夹,⭐重点)
python result_sqc/predistortion/P4_improvement_vs_strength/plot_improvement_vs_strength_sqc.py # P4 ⭐
python result_sqc/predistortion/P2_correction/plot_correction_before_after_sqc.py             # P2
```

所有脚本从仓库根目录运行(`sys.path` 已在 `_common.py` 处理)。
`reconstruction/` 下脚本深两级(上跳 `../..` 找 `_common`);块 2/3 现按「每图一夹」
也深两级(上跳 `../..`),各脚本 `SUBDIR` 指向自己的 `F<n>_*/` 或 `P<n>_*/` 子夹。

## 关键坑位 (fresh context 必读)

0. **type=4 阶跃信号的 sqc/src 幅度不一致(已核实)**:`FluxSignal(type=4,
   amplitude=0.01, center=100, rise=10, fall=10)` 在 sqc 侧峰值 ≈0.0037,
   在 src 侧 `Signal(type=4,…)` 峰值 ≈0.0111(约 3×)。两侧重建都忠实跟随
   **各自的** original(sqc wiener 1.14×、src 0.97×),故各自 pipeline 无误
   —— 差异在 type=4「非对称脉冲」的信号定义本身。看 `comparison_sqc.png`
   面板 (b)/(f) 时,src 虚线比 sqc 实线高约 3× 是**信号定义差异,不是重建
   误差**。其余三种(sine/double_peak 差 ~10%、complex 近乎一致)对得上。
   **不要**为了让 step 对齐去改 `sqc/control/flux_signal.py`(核心模块)或
   `src/`(R1 冻结)。

1. **结果字典的键位置**：`TransientSensingExperiment.run()` 返回的
   `ExperimentResult` 里 —
   - `res.data`: `p_e`, `delta_p`, `kernel`, `flux_samples`
   - `res.axes`: `scan`, `t_samples`, `t_flux`
   (真值波形是 `res.data["flux_samples"]`，时间轴是 `res.axes["t_flux"]`。)

2. **LM 路径需要适配器（重要）**：`_reconstruct_lm` 读的是
   `data["p_meas"]` / `axes["t_signal"]`，与 `TransientSensingExperiment`
   产出的 `p_e` / `t_flux` **键名不一致**（sqc 既有问题，commit 059bbe68，
   与本次 wiener /dt 工作无关）。直接对 transient 结果用 `method="lm"` 会
   `KeyError: 'p_meas'`（`SensingWorkflow` 的 lm 分支同样中招）。
   **必须先过 `_common.adapt_for_lm(res)`**。已验证:适配后 LM ratio≈1.0，
   返回 `(FluxSignal, history)` 元组，单信号约 1–2 分钟(全密度矩阵反演,慢)。

3. **Wiener 已验证与 src 逐位一致**，且都含 `/dt`。不要再"修" `/dt`
   (见 architecture.md v2.17；曾误判并回退)。λ 默认已是 5.0。

4. **控制脉冲**:LM 需要 `exp.control_pulse`,该属性在 `exp.run()` **之后**
   才被填充,故必须先 `run()` 再构造 LM reconstructor。

5. **D2(`D2_resolution_sensitivity`)坑位**:
   - **LM 实测 ~300 s/次**,不是坑位 2 里写的 1–2 min。规划 LM 点数时按 5 min/点
     算,否则预算会差 3–5 倍(D2 全套 LM 25 次 ≈ 2 小时)。
   - **LM 的 Δt_min 远高于 Wiener**:fourier `n_basis=100` + λ=100 过度平滑,
     Δt=10/12/14 ns 处 LM 只给出**单峰**(Cv 未定义),而 Wiener 已分辨两峰。
     选 LM 验证点必须往**大间距**方向取才能把 LM 自己的穿越点括住。
   - **真值自身有几何天花板**:两个 σ=5 ns 高斯在间距 ≤ 2σ=10 ns 时**输入本身
     就是单峰**,任何传感器都不可能分辨。报告 Δt_min 必须同时给出这条天花板,
     否则会把信号几何限制误记成传感器分辨率。
   - **`FluxSignal(type=5)` 不能做该扫描**:它把峰间距和单峰宽度都绑到 `width`。
     必须用 `type=8` 显式自定义数组。
   - **η_Φ 不可算**:`sqc.config` 有 `pulse.t_rabi_duration` 但**没有读出/复位
     时间**,占空比未知 → `eta_phi` 一律存 NaN,只报 A_min。**不得**拿求解器
     墙钟时间冒充采集时间。
   - **A_min 依赖估计器增益**:Wiener 幅度估计在扫描区间上端 `Ahat/A ≈ 1.2`
     (过恢复)。A_min 是**该估计器的**检测阈,不是增益校正后的物理幅度。

## src 基准函数签名 (做叠加对比时调用)

```python
from src.qubit import TransmonQubit      # (EC, EJ, T1, T2, flux, state, n_levels); .change_flux(); .optimal_work_point()
from src.signal import Signal            # (type, t_list, **params)  type 同 FluxSignal
from src.protocal import Protocal        # Protocal(4).initialize(q,state=0); .evolve(q) ->
                                         #   (t_samples, kernel, scan_list, delta_p, p_e, Phi_signal, control_pulse)
from src.analysis import Analysis        # .wiener_deconvolution(delta_p, kernel, dt, lambdas=5.0) -> (t, recon)
                                         # .numerical_inverse(q, control_pulse, p_e, t_list, Phi,
                                         #     basis_type, n_basis, lambdas, max_iter, tol) -> (B_opt, history)
```
`_common.make_src_qubit(flux)` 已封装 src 侧量子比特构造。
