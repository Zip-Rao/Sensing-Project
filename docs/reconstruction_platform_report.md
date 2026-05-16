# 重建平台技术报告 — `sqc/` 框架

> 内部技术文档 · 面向后续维护者与协作者
> 对应代码版本：`sqc/` v0.1.0（refactor 已完成 Phase 1–5）
> Notebook 主导引：[`Simulation_sqc.ipynb`](../Simulation_sqc.ipynb)
> 配套架构文档：[`docs/architecture.md`](architecture.md)

---

## 0. 导言

### 0.1 平台定位

本平台 (`sqc/`) 是一套面向超导 Transmon qubit 的"端到端时变磁通感知 + 控制线标定 + AWG 预失真"仿真框架。所有时间使用 ns、所有频率使用 GHz·(2π)，自然单位 ħ=1。底层物理求解器为 [QuTiP](https://qutip.org/) `mesolve`。

整套平台围绕一个**统一时间量子** `dt = 1 / awg.sample_rate` 组织（[`sqc/config.py:51-56`](../sqc/config.py#L51-L56)），所有时间轴一律用 `np.arange(start, stop, dt)` 派生（项目硬约束 R9）。

### 0.2 三条主线

| 主线 | 目标 | 输入 | 输出 | 主要模块 |
|------|------|------|------|----------|
| **A. 波形重建** | 反演未知磁通波形 Φ(t) | qubit 测量数据 (p_e, p_e_I, p_e_Q) | `FluxSignal` (Φ vs t) | `sqc/experiments/`, `sqc/reconstruction/` |
| **B. 频率标定** | 建立 f₀₁(Φ) 曲线 / 单点 f₀₁ | qubit + flux 扫描 | `CalibrationTable` (Φ→f 或 V→f) | `sqc/calibration/` |
| **C. 预失真** | 反求 AWG 输入使片上波形等于目标 | 测得的 H(ω) 或阶跃响应 | 逆 `DistortionModel` + 预补偿 AWG 波形 | `sqc/hardware/`, `sqc/calibration/waveform.py`, `sqc/workflows/` |

三条主线在数据上互相喂养：

```
                  ┌────────────────────────────────────┐
                  │    主线 B（频率标定）              │
                  │  Ramsey / closed-loop / transient  │
                  └──────────┬─────────────────────────┘
                             │  CalibrationTable(f_phi, f01)
                             ▼
  ┌──────────────────────────────────────────────────┐
  │  主线 A（波形重建）                              │
  │  Cryoscope / DelayRamsey / πPulseComp /Transient │
  └──────────┬───────────────────────────────────────┘
             │  Φ(t) 重建结果 / kernel(t) / transfer fn
             ▼
  ┌──────────────────────────────────────────────────┐
  │  主线 C（预失真）                                │
  │  WaveformCalibration → PredistortionDesigner     │
  │  → PredistortionValidationWorkflow / ZCrosstalk  │
  └──────────────────────────────────────────────────┘
```

### 0.3 共用基础设施

所有协议、所有重建、所有标定共享 4 个抽象基类与一个全局配置单例：

| 类 / 模块 | 文件 | 作用 |
|-----------|------|------|
| `CONFIG` (singleton) | [`sqc/config.py:233`](../sqc/config.py#L233) | 全局配置：AWG、qubit 默认、pulse 时间轴、reconstruction 默认超参 |
| `Experiment` (ABC) | [`sqc/experiments/base.py:10`](../sqc/experiments/base.py#L10) | 每个协议的统一接口：`build_sequence()` + `run() → ExperimentResult` |
| `Reconstruction` (ABC) | [`sqc/reconstruction/base.py:11`](../sqc/reconstruction/base.py#L11) | 重建算法接口：`reconstruct(measurement, ...) → FluxSignal` |
| `Calibration` (ABC) | [`sqc/calibration/base.py:167`](../sqc/calibration/base.py#L167) | 标定接口：`calibrate() → CalibrationTable` |
| `CalibrationTable` (dataclass) | [`sqc/calibration/base.py:19`](../sqc/calibration/base.py#L19) | 标定结果容器，含 `evaluate(x)` 正向 + `inverse(y)` 反向插值 |
| `Workflow` (ABC) | [`sqc/workflows/base.py:11`](../sqc/workflows/base.py#L11) | 高层工作流：编排多个 experiment + calibration + reconstruction |
| `IQReadoutModel` | [`sqc/hardware/readout.py:74`](../sqc/hardware/readout.py#L74) | I/Q 双通道读出（两次 Ramsey 序列，π/2 phase offset） |

`CalibrationTable.evaluate / inverse` 实现在 [`sqc/calibration/base.py:63-164`](../sqc/calibration/base.py#L63-L164)，统一使用 scipy `interp1d` cubic/quadratic/linear 自动降级。

### 0.4 文档导航

- §1 = 主线 A（波形重建，4 个协议）
- §2 = 主线 B（频率标定，2 个标定器）
- §3 = 主线 C（预失真，5 个失真模型 + 设计器 + 端到端 workflow）
- 附录 A = 公式速查
- 附录 B = 参考文献
- 附录 C = 测试与回归索引

每个协议的子章节统一按以下顺序写：

> 物理原理与目标 → 数学公式 → 实验层代码（experiments/）→ 重建层代码（reconstruction/）→ 配套标定（如有）→ Notebook 调用示例 → 关键超参与约束 → 局限与适用范围

---

# 主线 A · 波形重建

> 共四个子协议：**Cryoscope**、**Delay Ramsey**、**π 脉冲补偿**、**瞬态磁场**。它们都是把未知磁通 Φ(t) 反演出来，但探测窗口与精度权衡不同。

## 1.1 共用基石：相位 ⇄ 频率 ⇄ 磁通

Transmon 的基态-激发态频率随磁通的依赖关系（Gao 2021 §III.B Eq. 13–14）：

$$
\omega_Q(\Phi) \;=\; \sqrt{8\,E_J\,|\cos(\pi\Phi/\Phi_0)|\,E_C}\;-\;E_C
$$

代码实现：[`sqc/reconstruction/dispersion.py:11`](../sqc/reconstruction/dispersion.py#L11) `qubit_inverse_frequency()`，给定 Δω 反求 Φ：

$$
\cos(\pi\Phi)\;=\;\frac{(f_{\text{target}}+E_C)^2}{8\,E_J\,E_C},\qquad
\Phi\;=\;\frac{1}{\pi}\arccos(\mathrm{clip}(\text{ratio},0,1))
$$

```python
# sqc/reconstruction/dispersion.py:11-41
def qubit_inverse_frequency(dphi_dt: np.ndarray, qubit) -> np.ndarray:
    f_q = qubit.frequency
    f_target = f_q + dphi_dt
    EC = qubit.EC
    EJ0 = getattr(qubit, "EJ_0", qubit.EJ)
    ratio = (f_target + EC) ** 2 / (8.0 * EC * EJ0)
    ratio = np.clip(ratio, 0.0, 1.0)
    total_flux = np.arccos(ratio) / np.pi
    bias = getattr(qubit, "flux_bias", getattr(qubit, "flux", 0.0))
    return total_flux - bias
```

Ramsey-类协议（Cryoscope、Delay Ramsey）的相位累积公式（Gao 2021 §V.B）：

$$
\varphi(t) \;=\; \int_0^{t}\!\!\Delta\omega(\tau)\,\mathrm d\tau
\;\approx\; \int_0^t\!\big[\omega_Q(\Phi_{\mathrm{bias}}+h(\tau))-\omega_d\big]\,\mathrm d\tau
$$

由此 `dφ/dt = Δω(t)`，再通过上面的 dispersion 反演就得到 h(t)。

I/Q 双通道相位提取在 [`sqc/hardware/readout.py:97-160`](../sqc/hardware/readout.py#L97-L160)，跑两次 Ramsey 序列（phase2=0 → I 通道；phase2=π/2 → Q 通道）：

$$
p_e^I \;=\; \tfrac12\big[1+\cos\varphi\big],\qquad p_e^Q \;=\; \tfrac12\big[1+\sin\varphi\big]
$$

$$
\varphi \;=\; \mathrm{atan2}(0.5-p_e^I,\; p_e^Q-0.5)
$$

```python
# sqc/hardware/readout.py:122-160 (节选)
ctrl_I = create_ramsey_pulse(t_rabi, tau, omega_d=omega_d,
                             phase1=np.pi / 2, phase2=0.0, qubit=qubit)
ctrl_Q = create_ramsey_pulse(t_rabi, tau, omega_d=omega_d,
                             phase1=np.pi / 2, phase2=np.pi / 2, qubit=qubit)
# ... mesolve I 通道、Q 通道
return {"p_e_I": float(result_I.expect[0][-1]),
        "p_e_Q": float(result_Q.expect[0][-1])}
```

---

## 1.2 Cryoscope — 截断扫描重建波形

> **Notebook 对应**: §5（基础 demo）、§10.2.1（阶跃响应测量）
> **理论**: Rol 2020 Cryoscope (Appl. Phys. Lett. 116, 054001) · Gao 2021 §V.E

### 1.2.1 物理原理

施加待测的磁通波形 Φ(t)，在不同截断时刻 t_d 处把信号切掉，使 qubit 在 [0, t_d] 区间内自由演化，累积相位：

$$
\varphi(t_d) \;=\; \int_0^{t_d}\Delta\omega(\Phi(\tau))\,\mathrm d\tau
$$

数值微分得到瞬时频率失谐：

$$
\frac{\mathrm d\varphi}{\mathrm d t_d} \;=\; \Delta\omega\big(\Phi(t_d)\big) \;=\; 2\pi\,\Delta f\big(h(t_d)\big)
$$

最后用**两种反演策略**之一映射 Δω → h：

- **`"calibration"`** —— 用 `CalibrationTable.inverse()` 查 φ↔h 标定表（本质是 φ(h) = κ·h·τ 的非线性版本）。
- **`"response"`** —— 用解析 Transmon 色散 (§1.1) 反演 Δω → Φ。

### 1.2.2 实验层

[`sqc/experiments/cryoscope.py:26`](../sqc/experiments/cryoscope.py#L26) `CryoscopeExperiment`：

```python
@dataclass
class CryoscopeExperiment(Experiment):
    qubit: object
    flux_signal: FluxSignal | None = None
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())
    tau: float = field(default_factory=lambda: CONFIG.reconstruction.cryoscope_tau)  # 100 ns
    trunc_list: np.ndarray | None = None
    omega_d: float | None = None
```

`run()` 流程 ([`sqc/experiments/cryoscope.py:76-137`](../sqc/experiments/cryoscope.py#L76-L137))：

1. 构造 [`IQReadoutModel(tau, t_rabi, omega_d)`](../sqc/hardware/readout.py#L73)
2. 对每个截断时刻 `trunc` ∈ `trunc_list`：
   - `phi_truncated = flux_signal.copy(); phi_truncated.truncate(0, trunc)` —— 把信号截到 t=trunc
   - `qubit.qubit_in_mag(phi_truncated, frame=1, omega_d=...)` —— 让 qubit 看到截断后的磁通
   - `readout.measure(qubit)` → `p_e_I`, `p_e_Q`
3. 相位提取：

```python
# sqc/experiments/cryoscope.py:117-118
varphi = np.arctan2(0.5 - p_e_I, p_e_Q - 0.5)[::-1]
varphi = np.unwrap(varphi)
```

输出 `ExperimentResult(data={"varphi", "p_e_I", "p_e_Q"}, axes={"trunc"})`。

### 1.2.3 重建层

[`sqc/reconstruction/cryoscope.py:28`](../sqc/reconstruction/cryoscope.py#L28) `CryoscopeReconstruction`：

```python
@dataclass
class CryoscopeReconstruction(Reconstruction):
    tau: float = field(default_factory=lambda: CONFIG.reconstruction.cryoscope_tau)
    inversion: Literal["calibration", "response"] = "calibration"
    calibration: CalibrationTable | None = None
    qubit: object | None = None
    use_sg_filter: bool = False
    sg_window: int = 7
    sg_poly: int = 2
```

核心算法 ([`sqc/reconstruction/cryoscope.py:65-105`](../sqc/reconstruction/cryoscope.py#L65-L105))：

```python
def reconstruct(self, measurement, kernel=None, calibration=None, dt=None) -> FluxSignal:
    varphi = np.asarray(measurement.data["varphi"], dtype=float)
    trunc = np.asarray(measurement.axes["trunc"], dtype=float)
    if len(trunc) > 1 and trunc[0] > trunc[-1]:
        trunc = trunc[::-1]
    if dt is None:
        dt = float(trunc[1] - trunc[0])

    if self.use_sg_filter:
        from scipy.signal import savgol_filter
        dphi_dt = savgol_filter(varphi, window_length=self.sg_window,
                                polyorder=self.sg_poly, deriv=1, delta=dt)
    else:
        dphi_dt = np.gradient(varphi, trunc)

    if self.inversion == "calibration":
        phi_equiv = dphi_dt * self.tau              # 模拟 τ-长 Ramsey 相位
        h_recon = cal.inverse(phi_equiv)            # CalibrationTable.inverse()
    elif self.inversion == "response":
        h_recon = qubit_inverse_frequency(dphi_dt, self.qubit)

    return FluxSignal(type=8, t_list=trunc, signal=h_recon)
```

**关键点**：`phi_equiv = dphi_dt * tau` 把瞬时频率失谐 Δω 换算成一个等价的"τ-长 Ramsey 相位"，再走 [`CalibrationTable.inverse()`](../sqc/calibration/base.py#L105) 查表。这里 τ 必须与标定时使用的 τ 严格一致（标定时同样跑 τ-长 square pulse + IQ Ramsey）。

### 1.2.4 配套标定 `CryoscopeCalibration`

[`sqc/reconstruction/cryoscope.py:113`](../sqc/reconstruction/cryoscope.py#L113) `CryoscopeCalibration`：

```python
@dataclass
class CryoscopeCalibration(Calibration):
    qubit: object
    h_list: np.ndarray | None = None        # 扫的 flux 列表，默认 linspace(-0.03,0.03,51)
    tau: float = 100.0                       # 与重建器必须一致
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())
```

`calibrate()` 流程 ([`sqc/reconstruction/cryoscope.py:141-180`](../sqc/reconstruction/cryoscope.py#L141-L180))：

1. 对每个 h ∈ h_list 构造 τ-长 square pulse：
   ```python
   t_sig = CONFIG.pulse.make_time(0, t_total)   # t_total = 2*t_pi2 + tau
   signal[mask] = float(h)                       # mask 是 [t_pi2, t_pi2+tau]
   Phi = FluxSignal(type=8, t_list=t_sig, signal=signal)
   ```
2. IQ Ramsey 读取 → `p_e_I, p_e_Q`
3. 计算原始相位 `varphi_raw = arctan2(0.5 - p_e_I, p_e_Q - 0.5)`
4. **关键：基于 Transmon 解析模型反推 2π 跳变次数** ([`sqc/reconstruction/cryoscope.py:182-194`](../sqc/reconstruction/cryoscope.py#L182-L194))：
   ```python
   omega_q = np.sqrt(8.0 * EJ0 * np.abs(np.cos(np.pi * total_flux)) * EC) - EC
   varphi_theory = (omega_q - omega_d) * self.tau
   n_wraps = np.round((varphi_theory - varphi_raw) / (2.0 * np.pi))
   return varphi_raw + 2.0 * np.pi * n_wraps
   ```
   这是把 `arctan2` 输出的 ±π 包裹放回正确的分支，**用物理模型驱动 unwrap** 而非纯数值 unwrap——后者在大相位（>π）时容易跳错。
5. 返回 `CalibrationTable(kind="phi_h", inputs=h_list, outputs=unwrapped_phi)`。

### 1.2.5 Notebook 调用示例

完整 pipeline（[`Simulation_sqc.ipynb`](../Simulation_sqc.ipynb) §5，cell 12）：

```python
# === Cryoscope: 标定 → 测量 → 重建 ===
qubit = TransmonQubit(
    EC=0.2 * 2 * np.pi, EJ=10.0 * 2 * np.pi,
    T1=100.0e3, T2=50.0e3,
    flux=np.arctan(np.sqrt(2)) / np.pi,   # 这是 legacy "甜点"工作点
    n_levels=3,
)

# Step 1: 标定
cal = CryoscopeCalibration(qubit=qubit,
                            h_list=np.linspace(-0.03, 0.03, 21),
                            tau=50.0)
cal_table = cal.calibrate()

# Step 2: 测量
test_signal = FluxSignal(type=3, t_list=CONFIG.pulse.make_time(0, 80),
                          amplitude=0.01, center=40, width=5)
cryo_exp = CryoscopeExperiment(qubit=qubit, flux_signal=test_signal, tau=50.0)
meas = cryo_exp.run()

# Step 3: 重建（两种 inversion 对比）
recon_cal  = CryoscopeReconstruction(calibration=cal_table, tau=50.0, inversion="calibration")
recon_resp = CryoscopeReconstruction(qubit=qubit,            tau=50.0, inversion="response")
h_cal  = recon_cal.reconstruct(meas)
h_resp = recon_resp.reconstruct(meas)
```

**§10.2.1 阶跃响应特殊用法**（cell 29）：把 flux 工作点从甜点改到敏感点（`flux=0.25`），因为甜点处 dω/dΦ = 0（二次响应），φ(h) 是对称抛物线，1D inverse 区分不出 ±h。敏感点处 φ(h) 单调，可以反演。

```python
qubit_cryo = TransmonQubit(... flux=0.25, n_levels=3)   # ←敏感点
step_amp = 0.001   # 选小幅度保证 |φ_max| < π，避免 unwrap 风险
```

### 1.2.6 约束与陷阱

| 约束 | 出处 | 说明 |
|------|------|------|
| τ 必须与标定一致 | [`cryoscope.py:55-56` vs `:132`](../sqc/reconstruction/cryoscope.py#L55) | 否则 `phi_equiv = dphi_dt * tau` 量纲对不上 |
| flux 工作点必须敏感（dω/dΦ ≠ 0） | notebook cell 29 | 甜点 φ(h) 是抛物线，无法 1D 反演 |
| `|φ_max| < π` | notebook §10.2.1 | 否则需要更激进的 unwrap（基于物理模型那种） |
| `arctan2` 顺序：`(0.5−p_e^I, p_e^Q−0.5)` | [`cryoscope.py:117`](../sqc/experiments/cryoscope.py#L117) | 与 src 历史保持一致；调换顺序会得到 −φ |
| `trunc_list` 默认倒序 `[140:20:-1]` | [`cryoscope.py:70`](../sqc/experiments/cryoscope.py#L70) | 与 `src/protocal.py case 5` legacy 一致；输出 `varphi` 再 `[::-1]` 翻回 |

---

## 1.3 Delay Ramsey — 衰减拖尾测量

> **Notebook 对应**: §7
> **物理**: 用短 Ramsey 序列在阶跃下降沿后**滑动**采样，重建 flux pulse 的拖尾。

### 1.3.1 物理原理

施加方波 flux 脉冲，其下降沿后存在指数拖尾（IIR-type 失真）。在下降沿后延迟 t_d 处插入一段长度为 τ_R 的 Ramsey 序列：

$$
\varphi(t_d) \;=\; \int_{t_d}^{t_d+\tau_R}\!\!\Delta\omega(\tau)\,\mathrm d\tau
\;\approx\; \tau_R\cdot\Delta\omega\big(\bar\Phi_{\mathrm{tail}}(t_d)\big)
$$

近似式假设 τ_R 足够短，拖尾在窗内近似常数。此式直接给出 **"τ_R-窗口平均"** 而非瞬时值（窗口越长，时间分辨率越差）。

线性化（小 h 极限）：

$$
\varphi(t_d) \;\approx\; \tau_R \cdot \kappa \cdot \bar h(t_d)
$$

其中 κ = dω/dΦ 是 frequency sensitivity（[`src/qubit.py`](../src/qubit.py) `frequency_sensitivity`）。

### 1.3.2 实验层

[`sqc/experiments/delay_ramsey.py:40`](../sqc/experiments/delay_ramsey.py#L40) `DelayRamseyExperiment`：

```python
@dataclass
class DelayRamseyExperiment(Experiment):
    qubit: object
    flux_signal: FluxSignal | None = None
    t_d_list: np.ndarray | None = None
    tau_R: float = field(default_factory=lambda: CONFIG.reconstruction.delay_ramsey_tau)  # 20 ns
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())
    omega_d: float | None = None
    run_baseline: bool = False    # 默认 off：零信号时 arctan2(0,0) 不稳
    t_fall: float = 0.0           # 下降沿位置，t_d 相对此参考
```

`run()` 流程 ([`sqc/experiments/delay_ramsey.py:93-203`](../sqc/experiments/delay_ramsey.py#L93-L203)) 关键步骤：

1. 时间网格 `t_sig = CONFIG.pulse.make_time(0, 2*t_pi2 + tau_R)`
2. **可选 baseline** —— 零 flux 跑一次记录 `p_e_I_base, p_e_Q_base`（用 I/Q 复数空间减，避开 `arctan2(0,0)`）
3. 对每个 t_d ∈ t_d_list：
   - 把 flux 窗口化到自由演化区间 `[t_pi2, t_pi2+τ_R]`（**π/2 脉冲期间必须置零**，否则 flux 会失谐脉冲本身）：
     ```python
     # sqc/experiments/delay_ramsey.py:134-143
     signal = np.zeros(len(t_sig), dtype=float)
     free_mask = (t_sig >= t_pi2_end) & (t_sig <= t_pi2_end + tau_R)
     signal[free_mask] = np.array(
         [self.flux_signal.value_at(self.t_fall + t_d + float(t))
          for t in t_sig[free_mask]],
         dtype=float,
     )
     ```
   - `qubit.qubit_in_mag(phi_windowed, frame=1, omega_d=omega_d)`
   - `readout.measure(qubit)` → `p_e_I, p_e_Q`
4. 相位提取 + **复数空间 baseline 相减**（[`delay_ramsey.py:157-178`](../sqc/experiments/delay_ramsey.py#L157-L178)）：
   ```python
   varphi_base = float(np.arctan2(0.5 - p_e_I_base, p_e_Q_base - 0.5))
   varphi_shifted = np.arctan2(0.5 - p_e_I, p_e_Q - 0.5) - varphi_base
   varphi_shifted = (varphi_shifted + np.pi) % (2 * np.pi) - np.pi   # wrap 回 [-π,π]
   varphi = np.unwrap(varphi_shifted)
   ```
   这一步比直接 unwrap 原始相位更稳——大 t_d 处信号衰减，原始 `arctan2(0,0)` 在 ~3π/4 噪声底飘移，会污染 unwrap 的 2π 分支判断。

### 1.3.3 重建层

[`sqc/reconstruction/delay_ramsey.py:26`](../sqc/reconstruction/delay_ramsey.py#L26) `DelayRamseyReconstruction`：

```python
@dataclass
class DelayRamseyReconstruction(Reconstruction):
    inversion: Literal["response", "calibration"] = "response"
    qubit: object | None = None
    calibration: CalibrationTable | None = None
    tau_R: float | None = field(default_factory=lambda: CONFIG.reconstruction.delay_ramsey_tau)
```

两种反演路径：

**A. `inversion="response"`** —— 用 Transmon 解析色散反演（[`delay_ramsey.py:65-72`](../sqc/reconstruction/delay_ramsey.py#L65-L72)）：

```python
def _via_response(self, varphi, t_axis, tau) -> FluxSignal:
    phi = np.asarray(varphi, dtype=float).copy()
    if np.median(phi) > 0:
        phi = phi - 2.0 * np.pi             # 防止整体平移到错的分支
    h = qubit_inverse_frequency(phi / tau, self.qubit)
    return FluxSignal(type=8, t_list=t_axis, signal=h)
```

**B. `inversion="calibration"`** —— 用预先标定的 φ(z) 表反查（[`delay_ramsey.py:74-83`](../sqc/reconstruction/delay_ramsey.py#L74-L83)）：

```python
def _via_calibration(self, varphi, t_axis, cal) -> FluxSignal:
    phi = np.asarray(varphi, dtype=float).copy()
    cal_center = 0.5 * (cal.outputs.min() + cal.outputs.max())
    shift = cal_center - np.median(phi)
    n2pi = np.round(shift / (2.0 * np.pi))
    phi = phi + n2pi * (2.0 * np.pi)        # 把测量相位平移到标定表的覆盖区
    flux = cal.inverse(phi)
    return FluxSignal(type=8, t_list=t_axis, signal=flux)
```

### 1.3.4 配套标定 `DelayRamseyCalibration`

[`sqc/reconstruction/delay_ramsey.py:91`](../sqc/reconstruction/delay_ramsey.py#L91) `DelayRamseyCalibration` —— 与重建器拍配套，扫描已知 flux 高度 z，拟合斜率 k：

$$
\varphi_{\mathrm{cal}}(z) \;=\; k\cdot z \;+\; \text{intercept},\qquad k = \tau_R\cdot\kappa
$$

`calibrate()` ([`delay_ramsey.py:126-175`](../sqc/reconstruction/delay_ramsey.py#L126-L175)) 关键步骤：

1. 先在零 flux 跑一次 baseline，得 `varphi_base`
2. 对每个 z ∈ z_list（默认 `linspace(-0.02, 0.02, 41)`）施加方波 flux，跑 IQ Ramsey 测 φ
3. 同样做复数空间 baseline 相减 + wrap + unwrap
4. 一阶多项式拟合 `k, intercept = np.polyfit(z_list, varphi, 1)`
5. 返回 `CalibrationTable(kind="phi_z", inputs=z_list, outputs=varphi_unwrapped, fit_params={"k": k, ...})`

### 1.3.5 Notebook 调用示例

完整 pipeline（[`Simulation_sqc.ipynb`](../Simulation_sqc.ipynb) §7，cell 16）：

```python
qubit = TransmonQubit(EC=0.2*2*np.pi, EJ=10.0*2*np.pi, T1=100e3, T2=50e3,
                      flux=0.25, n_levels=2)
tau_R = 20.0

# Step 1: 构造一个带指数拖尾的方波（待测信号）
t_list = CONFIG.pulse.make_time(0, 100)
sig = np.zeros(len(t_list))
sig[(t_list >= 10) & (t_list <= 40)] = 0.008
sig[t_list > 40] = 0.008 * np.exp(-(t_list[t_list > 40] - 40) / 20)
test_flux = FluxSignal(type=8, t_list=t_list, signal=sig)

# Step 2: 标定 φ(z) 斜率
cal = DelayRamseyCalibration(qubit=qubit, z_list=np.linspace(-0.010, 0.010, 41),
                              tau_R=tau_R, t_rabi=np.linspace(0, 10, 20))
cal_table = cal.calibrate()
k = cal_table.fit_params["k"]    # rad/Φ₀

# Step 3: 测量（启用 baseline 相减）
exp = DelayRamseyExperiment(qubit=qubit, flux_signal=test_flux,
                             t_d_list=CONFIG.pulse.make_time(0, 60)[::2],
                             tau_R=tau_R, t_fall=40.0, run_baseline=True)
meas = exp.run()

# Step 4: 重建（用标定查表反演）
recon = DelayRamseyReconstruction(inversion="calibration", calibration=cal_table)
tail_flux = recon.reconstruct(meas)
```

### 1.3.6 约束与陷阱

| 约束 | 出处 | 说明 |
|------|------|------|
| 必须从下降沿后开始扫 `t_d` | `t_fall=40.0` notebook | 否则 flux 还在 hold 阶段，测的不是拖尾 |
| `τ_R` 必须远短于拖尾时间常数 | 物理 | 否则窗口平均把拖尾抹平；典型 τ_R = 10–30 ns |
| 自由演化窗内才有 flux，π/2 期间必须置零 | [`delay_ramsey.py:134-143`](../sqc/experiments/delay_ramsey.py#L134) | 否则 flux 失谐 π/2 脉冲 |
| baseline 相减必须在复数空间做 | [`delay_ramsey.py:171-176`](../sqc/experiments/delay_ramsey.py#L171) | 否则大 t_d 处的 `arctan2(0,0)` 噪声会污染 unwrap |
| 重建出的是**窗口平均** ≠ 瞬时 Φ | 推导式 | notebook §7 提供了 τ_R-window moving average 作为真值参考 |

---

## 1.4 π 脉冲补偿 — 二维扫描的"自归零"测量

> **Notebook 对应**: §8
> **物理**: 拖尾使 qubit 失谐；施加一个反向补偿 flux z 使其归零，π 脉冲翻转概率 P_e 达极大。

### 1.4.1 物理原理

设拖尾 flux 在 [τ, τ+T_π] 区间内平均为 ⟨Φ_tail⟩(τ)，叠加一个常数补偿 z 后总失谐为 (κ·(⟨Φ_tail⟩+z))。π 脉冲只有在失谐 ≈ 0 时才能完美翻转 |0⟩→|1⟩。扫描 z 找极大：

$$
z^*(\tau) \;=\; \arg\max_z P_e(\tau,z)
\quad\Longrightarrow\quad
\Phi_{\mathrm{tail}}(\tau) \;\approx\; -z^*(\tau)
$$

**注意**：由于 π 脉冲宽度 T_π > 0，`z*` 实际上是 [τ, τ+T_π] 窗口的反向平均，**不是瞬时拖尾**。指数拖尾的衰减时间常数仍能被保留，但幅度会被窗口平均稀释。

### 1.4.2 实验层

[`sqc/experiments/pi_pulse_comp.py:42`](../sqc/experiments/pi_pulse_comp.py#L42) `PiPulseCompensationExperiment`：

```python
@dataclass
class PiPulseCompensationExperiment(Experiment):
    qubit: object
    flux_signal: FluxSignal | None = None
    tau_list: np.ndarray | None = None       # delay 扫描
    z_list: np.ndarray | None = None         # 补偿幅度扫描
    T_pi: float = field(default_factory=lambda: CONFIG.reconstruction.pi_pulse_T_pi)  # 10 ns
    omega_bias: float | None = None
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())
    t_fall: float = 0.0
```

`run()` 流程（[`pi_pulse_comp.py:97-203`](../sqc/experiments/pi_pulse_comp.py#L97-L203)）核心：

```python
for i, tau in enumerate(self.tau_list):
    for j, z in enumerate(self.z_list):
        # 拖尾 + 补偿 z（整个 π 脉冲期间）
        tail = np.array(
            [self.flux_signal.value_at(self.t_fall + tau + float(t))
             for t in t_sig], dtype=float)
        signal = tail + float(z)

        phi_composite = FluxSignal(type=8, t_list=t_sig, signal=signal)
        self.qubit.qubit_in_mag(phi_composite, frame=1, omega_d=self.omega_bias)

        # 同时加 π 脉冲（在 t_rabi 时间窗内）
        H_pi = create_pulse(self.qubit, frame=1, type=1, t_list=t_sig,
                            omega_d=self.omega_bias, phase=0.0, angle=np.pi)
        H_total = (
            QobjEvo(self.qubit.H_list, tlist=self.qubit.mag_signal.t_list, order=1)
            + H_pi
        )
        result = mesolve(H_total, self.qubit.state, t_sig, [],
                          e_ops=[psi_e * psi_e.dag()])
        p_e_2d[i, j] = float(result.expect[0][-1])
```

提取 `z*` 时用**抛物线顶点**做子分辨率插值（[`pi_pulse_comp.py:167-182`](../sqc/experiments/pi_pulse_comp.py#L167-L182)）：

```python
for i in range(n_tau):
    j = int(np.argmax(p_e_2d[i]))
    if 0 < j < n_z - 1:
        pl, pc, pr = p_e_2d[i, j-1], p_e_2d[i, j], p_e_2d[i, j+1]
        denom = pr - 2.0 * pc + pl
        if abs(denom) > 1e-15:
            z_star[i] = self.z_list[j] - 0.5 * dz * (pr - pl) / denom
        else:
            z_star[i] = self.z_list[j]
    else:
        z_star[i] = self.z_list[j]
```

这避免了 z_list 分辨率太粗造成的"台阶状"重建结果。

### 1.4.3 重建层

[`sqc/reconstruction/pi_pulse_comp.py:21`](../sqc/reconstruction/pi_pulse_comp.py#L21) `PiPulseCompReconstruction`：

```python
@dataclass
class PiPulseCompReconstruction(Reconstruction):
    def reconstruct(self, measurement, **kwargs) -> FluxSignal:
        z_star = np.asarray(measurement.data["z_star"], dtype=float)
        t_axis = np.asarray(measurement.axes["tau"], dtype=float)
        return FluxSignal(type=8, t_list=t_axis, signal=-z_star)
```

—— 实际上**没有反演**：拖尾 = -z\*，已经在物理意义上是 1:1 的。

### 1.4.4 Notebook 调用示例

完整 pipeline（[`Simulation_sqc.ipynb`](../Simulation_sqc.ipynb) §8，cell 18）：

```python
qubit = TransmonQubit(EC=0.2*2*np.pi, EJ=10.0*2*np.pi, T1=100e3, T2=50e3,
                      flux=0.25, n_levels=2)
omega_bias = qubit.frequency

# 与 §7 共用同一个待测信号（一段方波 + 拖尾）
t_list = CONFIG.pulse.make_time(0, 100)
test_signal_arr = np.zeros(len(t_list))
test_signal_arr[(t_list >= 10) & (t_list <= 40)] = 0.008
test_signal_arr[t_list > 40] = (
    0.008 * np.exp(-(t_list[t_list > 40] - 40) / 20)
    + 0.008 * np.sin(2*np.pi * (t_list[t_list > 40] - 40) / 100)
)
test_flux = FluxSignal(type=8, t_list=t_list, signal=test_signal_arr)

# 2D 扫描
ppc_exp = PiPulseCompensationExperiment(
    qubit=qubit, flux_signal=test_flux,
    tau_list=CONFIG.pulse.make_time(0, 60)[::4],
    z_list=np.linspace(-0.012, 0.008, 17),
    T_pi=10.0, omega_bias=omega_bias,
    t_rabi=np.linspace(0, 10, 20),
    t_fall=40.0,
)
result = ppc_exp.run()

# 重建：Φ_tail = -z*
recon = PiPulseCompReconstruction()
tail_flux = recon.reconstruct(result)
```

### 1.4.5 约束与陷阱

| 约束 | 出处 | 说明 |
|------|------|------|
| z 扫描范围必须覆盖真实 -⟨Φ_tail⟩ | 物理 | 否则 argmax 在边界，子分辨率插值无效 |
| z 分辨率 dz 决定原始分辨率 | [`pi_pulse_comp.py:170`](../sqc/experiments/pi_pulse_comp.py#L170) | 抛物线插值能改善 ~5×，但仍受 `dz` 限制 |
| 测的是窗口平均 ≠ 瞬时 | 1.4.1 推导 | 与 Delay Ramsey 同样的窗口效应，但窗长是 T_π 而非 τ_R |
| 复杂度 O(n_tau · n_z) | 算法 | 单次 mesolve × 21 × 17 ≈ 357 次；最贵 |
| 不需要相位 unwrap | 无 | 直接看 P_e 峰，相位整圈跳变不影响 argmax |

---

## 1.5 瞬态磁场感知 — 滑动 Ramsey + 反卷积

> **Notebook 对应**: §4, §10.2.2
> **物理**: 把 Ramsey 序列当一个"探针"滑过未知信号，测得到的是信号与探针 kernel 的**卷积**。反卷积恢复信号。

### 1.5.1 物理原理

控制脉冲的"等效核函数" k(t) 描述了它对每个时刻磁通的灵敏度。在小信号近似下：

$$
\Delta p_e(t_d) \;=\; (k \ast B)(t_d) \;=\; \int k(t_d - \tau)\,B(\tau)\,\mathrm d\tau
$$

—— 即测量信号是真实信号 B(t) 与核 k(t) 的卷积。

反演策略三选一：

| 方法 | 公式 | 何时用 |
|------|------|--------|
| **Wiener** | $\hat B(\omega) = \dfrac{H^*(\omega)\,\Delta P(\omega)}{|H(\omega)|^2 + \lambda^2}$ | 线性区间，最快 |
| **Hammerstein-Wiener** | Wiener 出 Δω → 解析 Transmon 反推 Φ | 大信号、非线性 |
| **Levenberg-Marquardt** | 基函数参数化 B + 全密度矩阵正向仿真 + 牛顿-高斯 | 数据贵、精度要求高 |

### 1.5.2 核函数估计 `KernelEstimator`

[`sqc/reconstruction/kernel.py:27`](../sqc/reconstruction/kernel.py#L27) 统一替换了 src 的三处重复实现（`Pulse.get_kernel`、`CompositePulse.get_kernel`、`Analysis.get_kernel`）：

```python
@dataclass
class KernelEstimator:
    stim_amplitude: float = field(default_factory=lambda: CONFIG.reconstruction.stim_amplitude)  # 0.0215
    stim_width: float = field(default_factory=lambda: CONFIG.reconstruction.stim_width)          # 3.0 ns
    auto_calibrate: bool = False
```

`estimate(pulse, qubit)` 流程 ([`kernel.py:53-149`](../sqc/reconstruction/kernel.py#L53-L149))：

1. 跑一次基线（无刺激）→ `p_e_base`
2. 对每个时刻 t_i 注入窄 Gaussian flux 刺激（amplitude=0.0215, width=3 ns, center=t_i），跑一次 → `p_e_stim`
3. `kernel[i] = (p_e_stim - p_e_base) / stim_area`，其中 `stim_area = ∫ stim(t) dt`

这是 **Mehmet Canturk 式** 的"用窄 Gaussian 当 δ-函数"做线性响应函数测量。`auto_calibrate=True` 时按 `kappa = dω/dΦ` 自动调 amplitude 让频率漂移 ~ 1% × |anharmonicity|（[`kernel.py:151-176`](../sqc/reconstruction/kernel.py#L151-L176)）。

### 1.5.3 实验层

[`sqc/experiments/transient.py:24`](../sqc/experiments/transient.py#L24) `TransientSensingExperiment`：

```python
@dataclass
class TransientSensingExperiment(Experiment):
    qubit: object
    flux_signal: FluxSignal | None = None
    flux_signal_zero: FluxSignal | None = None    # 0-flux 参考
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())
    omega_d: float | None = None
    scan_list: np.ndarray | None = None
```

`run()` 流程 ([`transient.py:89-144`](../sqc/experiments/transient.py#L89-L144))：

1. 构造 Ramsey 控制脉冲 `create_ramsey_pulse(t_rabi, tau=0.0, omega_d, qubit)`（**tau=0** 关键 —— 让两段 π/2 紧贴）
2. 用 `SlidingMeasurementRunner` 把 Ramsey 在 flux_signal 上滑动 → `p_e(scan)`
3. 同样在零 flux 上滑动 → `p_e_base(scan)`
4. `delta_p = p_e - p_e_base`
5. 同时调用 `control_pulse.get_kernel(qubit)` 计算 `kernel` 与 `t_samples`

输出 `ExperimentResult(data={"p_e", "delta_p", "kernel", "flux_samples"}, axes={"scan", "t_samples", "t_flux"})`。

### 1.5.4 重建层

[`sqc/reconstruction/transient.py:45`](../sqc/reconstruction/transient.py#L45) `TransientReconstruction` 支持三种 method：

#### A. Wiener 反卷积

[`transient.py:124-157`](../sqc/reconstruction/transient.py#L124-L157)：

```python
def _reconstruct_wiener(self, measurement, kernel, dt=None, **__) -> FluxSignal:
    delta_p = np.asarray(measurement.data["delta_p"], dtype=float)
    if dt is None:
        scan = np.asarray(measurement.axes["scan"])
        dt = scan[1] - scan[0]

    N_del = len(delta_p); N_ker = len(kernel)
    N = N_del - N_ker + 1
    N_fft = N_del

    p_pad = np.zeros(N_fft); p_pad[:N_del] = delta_p
    k_pad = np.zeros(N_fft); k_pad[:N_ker] = kernel

    Y = np.fft.fft(p_pad)
    H = np.fft.fft(k_pad)
    G = np.conj(H) / (np.abs(H) ** 2 + self.lambda_reg ** 2)   # ← Wiener 滤波器
    X_w = Y * G / dt
    x_rec = np.real(np.fft.ifft(X_w))[:N]

    return FluxSignal(type=8,
                       t_list=np.linspace(0, (N - 1) * dt, N),
                       signal=x_rec)
```

正则化 `λ` 由 `CONFIG.reconstruction.lambda_reg = 10.0` 默认控制；**小信号、阶跃响应时常需要调到 ~50** 抑制 ringing（notebook §10.2.2 cell 31）。

#### B. Hammerstein-Wiener

[`transient.py:163-182`](../sqc/reconstruction/transient.py#L163-L182)：先用 Wiener 得到等效频率失谐 ω(t)，再用解析 Transmon 色散反推磁通：

```python
def _reconstruct_hammerstein(self, measurement, kernel, dt=None, **__) -> FluxSignal:
    omega_signal = self._reconstruct_wiener(measurement, kernel, dt=dt)
    omega = np.asarray(omega_signal.signal, dtype=float)

    EC, EJ, freq = _get_qubit_params(self.qubit)
    ratio = np.clip((omega + freq + EC) ** 2 / (8 * EC * EJ), 0.0, 1.0)
    B = (1.0 / np.pi) * np.arccos(ratio)

    bias = getattr(self.qubit, "flux_bias", getattr(self.qubit, "flux", 0.0))
    B = B - bias    # 减掉工作点 bias，得 h = Φ − Φ_bias
    return FluxSignal(type=8, t_list=omega_signal.t_list.copy(), signal=B)
```

#### C. Levenberg-Marquardt

[`transient.py:188-490`](../sqc/reconstruction/transient.py#L188-L490) —— 把 B(t) 用基函数展开 `B(t) = Σ_k b_k φ_k(t)`，最小化残差 `||p_meas - p_sim(b)||²` 求 b：

$$
b^{(n+1)} \;=\; b^{(n)} + \Big(J^{\mathsf T}J + \mu I + \lambda R\Big)^{-1}\,J^{\mathsf T}\,\big(p_{\text{meas}}-p_{\text{sim}}(b^{(n)})\big)
$$

其中 R 是 smoothness 正则化矩阵（`regularization_matrix`，[`sqc/reconstruction/basis.py:138-170`](../sqc/reconstruction/basis.py#L138-L170)）。

基函数生成在 [`basis.py:20-90`](../sqc/reconstruction/basis.py#L20-L90)，三种选项 `bspline / fourier / legendre`。`fourier` 对应物理上常见的"频谱稀疏"先验，`bspline` 对应"光滑"先验，`legendre` 对应"低多项式"先验。

Jacobian 有两种实现：

- **伴随法（adjoint, 默认 `use_adjoint=True`）** —— [`transient.py:301-390`](../sqc/reconstruction/transient.py#L301-L390)。物理上是把测量算符 μ_0 = |1⟩⟨1| 反向演化到 t_evolve[0]，再与正向演化的密度矩阵 ρ(t) 做迹运算：
  $$
  J_{ik} \;=\; \int_0^{T}\!\!\mathrm{Tr}\Big[\lambda(t)\,\big[G,\rho(t)\big]\Big]\,\frac{\partial\Delta\omega}{\partial\Phi}\,\phi_k(t)\,\mathrm dt
  $$
  其中 G = `qubit.n`（数算符）。
- **有限差分** —— [`transient.py:396-420`](../sqc/reconstruction/transient.py#L396-L420)。慢但稳。

正向仿真在 [`_forward_simulation`](../sqc/reconstruction/transient.py#L272-L295)：对每个 t_delay 单独跑一次 `mesolve(H_total, qubit.state, t_evolve, ...)`。

> **已知问题**：adjoint Jacobian 与有限差分的偏差比预期大（master TODO 0.3），但 LM 仍能收敛——只是迭代次数偏多。生产环境推荐先用 Wiener / Hammerstein 做 warm start，再用 LM 精修。

### 1.5.5 Notebook 调用示例

**基础 demo**（§4，cell 10）：

```python
qubit = TransmonQubit(EC=0.2*2*np.pi, EJ=10.0*2*np.pi, T1=100e3, T2=50e3,
                      flux=0.9553, n_levels=2)   # 几乎到 sweet spot 附近
trans_exp = TransientSensingExperiment(qubit=qubit)
trans_result = trans_exp.run()
kernel = trans_result.data["kernel"]

rec = TransientReconstruction(method="wiener",
                               lambda_reg=CONFIG.reconstruction.lambda_reg)
B_rec_signal = rec.reconstruct(trans_result, kernel=kernel)
```

**阶跃响应测量**（§10.2.2，cell 31）：

```python
qubit_sens = qubit_cryo   # 复用 §10.2.1 的敏感点 qubit (flux=0.25)
t_scan = CONFIG.pulse.make_time(0, 70)
trans_exp = TransientSensingExperiment(
    qubit=qubit_sens,
    flux_signal=on_chip_flux,    # 来自 §10.2.1 的失真后片上磁通
    scan_list=t_scan,
)
trans_result = trans_exp.run()
kernel = trans_result.data["kernel"]

lambda_reg = 50.0    # 阶跃信号弱→强正则化抑 ringing
rec_wiener = TransientReconstruction(method="wiener", lambda_reg=lambda_reg)
rec_hammer = TransientReconstruction(method="hammerstein",
                                       lambda_reg=lambda_reg, qubit=qubit_sens)

B_wiener = rec_wiener.reconstruct(trans_result, kernel=kernel)
B_hammer = rec_hammer.reconstruct(trans_result, kernel=kernel)
```

### 1.5.6 约束与陷阱

| 约束 | 出处 | 说明 |
|------|------|------|
| Wiener 假设线性 | 推导 1.5.1 | 大信号下相位被压（cosθ ≠ 1−θ²/2），需要 Hammerstein |
| `λ_reg` 决定噪声/分辨率权衡 | `CONFIG.reconstruction.lambda_reg` | 噪声大或信号弱时 λ ↑（如 §10.2.2 用 50.0） |
| kernel 必须重新算 | `TransientSensingExperiment.run()` 末尾 | qubit 参数变了 kernel 也变 |
| stim_amplitude 0.0215 是 legacy 默认 | [`kernel.py:46`](../sqc/reconstruction/kernel.py#L46) | qubit 灵敏度变化时考虑 `auto_calibrate=True` 或显式覆盖 |
| LM 启动需要 `qubit_in_mag` 已调用 | [`transient.py:411,419,437`](../sqc/reconstruction/transient.py#L411) | adjoint 计算依赖 `qubit.freq_coeffs` |
| LM Jacobian 偏差已知 | TODO 0.3 | 不影响收敛，影响一阶 step 大小 |

---

## 1.6 四种波形重建方法对比

> 来源：notebook §10.2.3（cell 33）总结 + 本文档对各方法的代码层分析。

| 维度 | Cryoscope | Delay Ramsey | π-pulse Comp. | Transient (Wiener) |
|------|-----------|--------------|---------------|--------------------|
| **测量量** | φ(t_trunc) | φ(t_d) | P_e(τ, z) | Δp_e(scan) |
| **解相位** | atan2 + 物理-driven unwrap | atan2 + baseline 减 + unwrap | argmax + 抛物线插值 | 不需要 |
| **时间分辨率** | 一个截断采样 | τ_R 移动窗 | T_π 移动窗 | AWG dt |
| **幅度分辨率** | 标定查表 / 解析反演 | k·z + 标定查表 / 解析反演 | dz · 0.5×插值 | λ_reg / SNR |
| **测一次的代价** | 1× mesolve / trunc | 2× mesolve / t_d (含 baseline) | n_z × mesolve / τ | 2 × n_scan mesolve (含 baseline) + kernel |
| **典型应用** | 阶跃响应（>10 ns IIR 尾） | 方波拖尾 | 方波拖尾、复杂波形 | 高频成分、宽带瞬态 |
| **工作点** | 敏感点（dω/dΦ≠0） | 敏感点 | 敏感点 | 任意（但越敏感越好） |
| **代码** | `experiments/cryoscope.py` + `reconstruction/cryoscope.py` | `experiments/delay_ramsey.py` + `reconstruction/delay_ramsey.py` | `experiments/pi_pulse_comp.py` + `reconstruction/pi_pulse_comp.py` | `experiments/transient.py` + `reconstruction/transient.py` |

**Notebook §10.2.3 给出的实证比较**（阶跃响应同一组 qubit/失真/step_amp 上）：

| 方法 | RMSE (Φ₀) | 上升沿 (ns) | 采样点 | 备注 |
|------|-----------|-------------|--------|------|
| 真实值 | — | ~3 ns | — | `dist.step_response()` |
| Cryoscope | ~1×10⁻⁵ | 与真实接近 | 60+ trunc | 相位精度高、需要 |φ_max|<π |
| Transient (Wiener) | ~5×10⁻⁵ | 略宽 | 70+ scan | 时间分辨率好、噪声敏感 |

**最佳实践**：Cryoscope 处理 >10 ns 的 IIR-type 尾，Transient 处理 <10 ns 的 FIR-type 快瞬态——二者互补。

---

# 主线 B · Qubit 频率标定

> 目标：建立 f₀₁ 在 (Φ, V) 参数空间上的"刻度"。两个标定器类、四种 method、一个调度器。

## 2.1 标定层概览

```
                       sqc.calibration
                              │
        ┌─────────────────────┼─────────────────────────────────┐
        │                     │                                 │
   CalibrationTable      Calibration (ABC)               CalibrationScheduler
   (base.py:19)          (base.py:167)                   (scheduler.py:19)
        │                     │                                 │
   evaluate(x)        ┌───────┴───────┬─────────────┐    register / run / run_next
   inverse(y)         │               │             │
                      │               │             │
              FluxResponseCal   SinglePointFreqCal  WaveformCal
              (frequency.py)   (frequency.py)      (waveform.py)
              ├ "ramsey"        ├ "ramsey"          ├ "transfer_function"
              └ "transient"*    ├ "closed_loop"     └ "predistortion"
                                └ "transient"
                                  ├ measure="ramsey"
                                  └ measure="transient"
              * Track B 1.2 待实现           Δω = p_diff / G_α
```

注：`FluxResponseCalibration / SinglePointFrequencyCalibration` 是频率主线的核心；`WaveformCalibration / PredistortionDesigner` 形式上也属于 calibration 模块，但属于主线 C（§3.3 详述）。

`CalibrationScheduler` 是一个轻量的"控制室"：

- **注册** ([`scheduler.py:61-102`](../sqc/calibration/scheduler.py#L61-L102))：维护任务名 → Calibration 类的映射，及依赖关系。
- **运行** ([`scheduler.py:108-172`](../sqc/calibration/scheduler.py#L108-L172))：`run(name)` 检查依赖 → 自动注入依赖结果（如 `flux_response_ramsey` 的 inputs 自动作为 `closed_loop` 的 V_a/V_b bracket）→ 执行。
- **DAG 自动化** ([`scheduler.py:261-316`](../sqc/calibration/scheduler.py#L261-L316))：`check_state / maintain / auto_calibrate` 是 Kelly 2018 风格的 stub，留给未来。

预注册任务列表（[`scheduler.py:82-102`](../sqc/calibration/scheduler.py#L82-L102)）：

| 任务名 | 类 | 依赖 |
|--------|----|------|
| `flux_response_ramsey` | `FluxResponseCalibration(method="ramsey")` | — |
| `frequency_ramsey` | `SinglePointFrequencyCalibration(method="ramsey")` | — |
| `frequency_closed_loop` | `SinglePointFrequencyCalibration(method="closed_loop")` | `flux_response_ramsey` |
| `waveform_transfer_function` | `WaveformCalibration(method="transfer_function")` | — |
| `waveform_predistortion` | `WaveformCalibration(method="predistortion")` | `waveform_transfer_function` |

## 2.2 共用基石：Ramsey FFT 频率拟合

[`sqc/calibration/frequency.py:28`](../sqc/calibration/frequency.py#L28) `_fit_ramsey_frequency()` —— 跑一组 Ramsey 扫 τ，FFT 找 detuning 峰：

```python
def _fit_ramsey_frequency(qubit, omega_d, tau_list, t_rabi, t_global, flux=0.0):
    # 1. 把 qubit 钉在指定 flux
    Phi = FluxSignal(type=1 if flux != 0.0 else 0, t_list=t_sig,
                      amplitude=float(flux))
    qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

    # 2. 扫 tau 跑 Ramsey
    for i, tau in enumerate(tau_list):
        ctrl = create_ramsey_pulse(t_rabi, tau, omega_d=omega_d,
                                    phase1=0.0, phase2=0.0, qubit=qubit)
        ctrl.t_list = ctrl.t_list - t_rabi[-1]
        H = (QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1)
             + QobjEvo(ctrl.hamiltonian, tlist=ctrl.t_list, order=1))
        result = mesolve(H, qubit.state, t_global, [], e_ops=[psi_e * psi_e.dag()])
        p_e_vals[i] = float(result.expect[0][-1])

    # 3. FFT + 二次插值精确找峰
    p_centered = p_e_vals - np.mean(p_e_vals)
    fft_vals = np.abs(rfft(p_centered, n=2048))
    freqs = rfftfreq(2048, d=dt_val)
    peak_idx = int(np.argmax(fft_vals))

    # 子-bin 二次插值（关键：把分辨率从 1/n_fft 提升到连续）
    y1, y2, y3 = fft_vals[peak_idx-1], fft_vals[peak_idx], fft_vals[peak_idx+1]
    denom = 2 * (y1 + y3 - 2 * y2)
    delta_bin = (y1 - y3) / denom
    detuning_hz = (peak_idx + delta_bin) * df_bin
    return float(omega_d) + 2 * np.pi * detuning_hz
```

**子-bin 二次插值**（[`frequency.py:103-117`](../sqc/calibration/frequency.py#L103-L117)）是这个标定器的核心——把"找峰的分辨率"从 `1/T_total` 提到连续值，FFT bin 数 2048 时 ~kHz 级。

## 2.3 `FluxResponseCalibration` — 扫 Φ 拟 f(Φ)

[`sqc/calibration/frequency.py:125`](../sqc/calibration/frequency.py#L125) `FluxResponseCalibration`：

```python
@dataclass
class FluxResponseCalibration(Calibration):
    qubit: object
    method: Literal["ramsey", "transient"] = "ramsey"
    h_list: np.ndarray | None = None          # default linspace(-0.03, 0.03, 51)
    tau: float = 100.0
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())
```

`_calibrate_ramsey()` ([`frequency.py:170-192`](../sqc/calibration/frequency.py#L170-L192))：

1. 对每个 h ∈ h_list 调用 `_fit_ramsey_frequency(qubit, omega_d, tau_list, t_rabi, t_global, flux=h)`
2. 输出 `CalibrationTable(kind="f_phi", inputs=h_list, outputs=f01_list)`

`_calibrate_transient()` 是 stub —— **Track B 1.2 未实现**（在 sqc 框架中保留接口位置，实际算法待补）。

**用途**：

- 给 closed-loop 标定提供 `V_a, V_b` bracket（scheduler 自动注入）。
- 作为 Cryoscope/Delay Ramsey 标定的"宽视野" sanity check。

## 2.4 `SinglePointFrequencyCalibration` — 单点 f₀₁

[`sqc/calibration/frequency.py:207`](../sqc/calibration/frequency.py#L207) 三种 method：

| method | 输出 | 何时用 |
|--------|------|--------|
| `"ramsey"` | f₀₁ at flux=0 | 单点测量、bare freq 校验 |
| `"closed_loop"` | f₀₁ at 目标值附近的电压 V_opt | 把 qubit 调到 target frequency |
| `"transient"` | (与 `"ramsey"` 同接口) | stub，Track B 1.2 |

### 2.4.1 `"ramsey"` 单次

[`frequency.py:303-321`](../sqc/calibration/frequency.py#L303-L321)：直接调一次 `_fit_ramsey_frequency(..., flux=0.0)` 包成 `CalibrationTable(kind="f01")`。

### 2.4.2 `"closed_loop"` 反馈控制

[`frequency.py:327-352`](../sqc/calibration/frequency.py#L327-L352) `_calibrate_closed_loop()`：把"找 V 使 f(V) = f_target"当成根求解问题 r(V) = f_Q(V) − f_target = 0，两种算法：

**A. Secant 法（默认 `step_method="secant"`）** ([`frequency.py:355-400`](../sqc/calibration/frequency.py#L355-L400))：

$$
V_{n+1} \;=\; V_n \;-\; r_n\cdot\frac{V_n - V_{n-1}}{r_n - r_{n-1}}
$$

附加 **bracket tightening (regula falsi)**（[`frequency.py:393-398`](../sqc/calibration/frequency.py#L393-L398)）：每步把 [V_lo, V_hi] 裹紧到当前估计点，防止 secant 跳出有效区间。

**B. Bisection 法（`step_method="bisection"`）** ([`frequency.py:403-459`](../sqc/calibration/frequency.py#L403-L459)）：

$$
V_n \;=\; \tfrac12 (V_{\mathrm{lo}}+V_{\mathrm{hi}}),
\qquad
\text{bracket width} \;\sim\; \frac{V_b - V_a}{2^n}
$$

特别处理 **even f(Φ)**（甜点处 f 对称，r(V_lo) 和 r(V_hi) 同号）：自动在中点劈半 bracket（[`frequency.py:418-432`](../sqc/calibration/frequency.py#L418-L432)）。

### 2.4.3 `measure_method="transient"`：核函数法测频率

[`frequency.py:504-582`](../sqc/calibration/frequency.py#L504-L582) `_measure_frequency_transient()` —— closed-loop 的每次迭代不一定都跑 Ramsey FFT（贵），可以用瞬态核函数法只跑 2 次测量得到 detuning：

构造两个正交 Ramsey 序列：`ctrl_x = (R_y, R_x)` 与 `ctrl_mx = (R_y, R_{-x})`，τ=0，跑 mesolve 得 p_x, p_mx。差分：

$$
p_{\text{diff}} \;=\; \tfrac12 (p_x - p_{-x})
$$

理论上（[_sensing theory.md](../note/_sensing%20theory.md) §瞬态磁场协议）：

$$
p_{\text{diff}} \;\approx\; G_\alpha \cdot \Delta\omega,\qquad
G_\alpha \;=\; \int k(t)\,\mathrm dt
$$

所以 `delta_omega = p_diff / G_alpha`：

```python
# sqc/calibration/frequency.py:543-582 (节选)
ctrl_x = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=omega_d,
                              phase1=np.pi/2, phase2=0.0, qubit=qubit)
ctrl_mx = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=omega_d,
                               phase1=np.pi/2, phase2=np.pi, qubit=qubit)
# ... mesolve 两次得 p_x, p_mx
p_diff = (p_x - p_mx) / 2.0

# kernel sensitivity G_α = ∫ k(t) dt
ctrl_x.get_kernel(qubit)
G_alpha = float(np.trapezoid(ctrl_x.kernel, ctrl_x.t_samples))

delta_omega = p_diff / G_alpha
return float(omega_d + delta_omega)
```

**优势**：每次迭代仅 2 次 mesolve + 1 次 kernel 计算（kernel 可以 caching），而 Ramsey FFT 需要扫 `len(tau_list)` 次。

### 2.4.4 Notebook 调用示例 —— 闭环单点频率标定

完整 pipeline（[`Simulation_sqc.ipynb`](../Simulation_sqc.ipynb) §9，cells 20–22）：

**Step 1: 粗标定 f(Φ) 给 bracket**（cell 20）：

```python
qubit = TransmonQubit(EC=0.2*2*np.pi, EJ=10.0*2*np.pi,
                      T1=100e3, T2=50e3, flux=0.0, n_levels=2)
flux_cal = FluxResponseCalibration(qubit=qubit, method="ramsey",
                                    h_list=np.linspace(-0.03, 0.03, 15), tau=100.0)
flux_table = flux_cal.calibrate()
V_a, V_b = float(np.min(flux_table.inputs)), float(np.max(flux_table.inputs))

# 选个目标：Φ=0.025 处的频率
f_interp = interp1d(flux_table.inputs, flux_table.outputs, kind="cubic")
f_target = float(f_interp(0.025))
```

**Step 2: 两种闭环（Ramsey FFT vs Transient kernel） Bisection**（cell 21）：

```python
# 方法 A: Ramsey FFT + Bisection
cal_ramsey = SinglePointFrequencyCalibration(
    qubit=qubit, method="closed_loop", measure_method="ramsey",
    step_method="bisection",
    tau_list=CONFIG.pulse.make_time(0, 500),    # 长 τ 提升 FFT 分辨率
    f_target=f_target,
    epsilon_f=1e-4 * 2 * np.pi,
    V_a=V_a, V_b=V_b, max_iter=20,
)
res_r = cal_ramsey.calibrate()

# 方法 B: Transient kernel + Bisection
cal_transient = SinglePointFrequencyCalibration(
    qubit=qubit2, method="closed_loop", measure_method="transient",
    step_method="bisection",
    f_target=f_target,
    epsilon_f=1e-4 * 2 * np.pi,
    V_a=V_a, V_b=V_b, max_iter=20,
)
res_t = cal_transient.calibrate()
```

**Step 3: 收敛对比可视化**（cell 22）展示 6 张子图：

1. (a) f(Φ) 曲线 + target 标注
2. (b) V vs iter 收敛轨迹
3. (c) f vs iter
4. (d) |residual| log 衰减
5. (e) bracket width log 衰减（验证 1/2ⁿ 理论线）
6. (f) Summary 表（迭代次数 + 单次测量代价 + 总测量数）

**关键观察**（notebook 最后一行打印）：transient kernel 每次迭代仅需 2 次单点测量，比 Ramsey FFT 快 `len(tau_list)/2` 倍。精度相当。

### 2.4.5 约束与陷阱

| 约束 | 出处 | 说明 |
|------|------|------|
| Bisection 需要 r(V_a)·r(V_b) < 0 | [`frequency.py:418`](../sqc/calibration/frequency.py#L418) | 甜点处 f(Φ) 偶函数，自动 split bracket 算法已内置 |
| Secant 可能跳出 bracket | [`frequency.py:377-380`](../sqc/calibration/frequency.py#L377-L380) | 设了 fallback 到 midpoint |
| Ramsey 分辨率 ≈ 1/(2·tau_max) | FFT | tau_long=500 → ~1 MHz 极限 |
| Transient kernel 假设线性 | 推导 2.4.3 | 大 detuning 时 |Δω·tau_R| > 1，公式失真 |
| Kernel 必须在 qubit 当前 flux 下重新算 | [`frequency.py:570-576`](../sqc/calibration/frequency.py#L570) | 否则 G_α 不对 |

---

# 主线 C · 预失真

> 目标：测出控制线的传递函数 H(ω)，设计逆滤波器，让 AWG → 片上波形误差最小。

## 3.1 失真模型层

[`sqc/hardware/distortion.py`](../sqc/hardware/distortion.py) 定义了 5 个 LTI 失真模型，全部实现 `DistortionModel` ABC（[`distortion.py:31`](../sqc/hardware/distortion.py#L31)）的 4 个接口：`apply()`, `step_response()`, `impulse_response()`, `frequency_response()`。

### 3.1.1 `SingleExponentialDistortion` — 单指数尾

[`distortion.py:105`](../sqc/hardware/distortion.py#L105)。这是 Cryoscope 标定中最常见的 IIR-type 失真：

$$
h(t) \;=\; (1-\alpha)\,\delta(t) \;+\; \frac{\alpha}{\tau}\,e^{-t/\tau}\,\Theta(t)
$$

$$
H(s) \;=\; \frac{1 + s\tau(1-\alpha)}{1 + s\tau},\qquad
H(\omega) \;=\; (1-\alpha) \;+\; \frac{\alpha}{1 + j\omega\tau}
$$

$$
s(t) \;=\; 1 \;-\; \alpha\,e^{-t/\tau} \quad (\text{step response})
$$

代码 ([`distortion.py:152-179`](../sqc/hardware/distortion.py#L152-L179))：用 bilinear transform 把连续时间 H(s) 离散化成一阶 IIR `(b, a)`，再用 scipy `lfilter` 应用。

### 3.1.2 `MultiExponentialDistortion` — K 个指数并联

[`distortion.py:187`](../sqc/hardware/distortion.py#L187)：

$$
H(\omega) \;=\; (1-\sum_k\alpha_k) + \sum_k\frac{\alpha_k}{1+j\omega\tau_k},\qquad
s(t) \;=\; 1 - \sum_k\alpha_k\,e^{-t/\tau_k}
$$

`apply()` ([`distortion.py:218-237`](../sqc/hardware/distortion.py#L218-L237))：DC gain 直传 + K 个并联 IIR 累加。

### 3.1.3 FIR / IIR / Custom

| 类 | 描述 | 应用 |
|----|------|------|
| `FIRDistortion` | y[n] = Σ b[k]·x[n-k] | 短脉冲响应 |
| `IIRDistortion` | y[n] = Σ b[k]·x[n-k] − Σ a[k]·y[n-k] | 反馈系统 |
| `CustomTransferDistortion` | 用户提供 H(ω) on grid | 从 Cryoscope/Transient 测得的非参数 H |

`apply_to_waveform()` ([`distortion.py:84-97`](../sqc/hardware/distortion.py#L84-L97)) 是统一适配器：接收 `Waveform`，调底层 `apply(samples, dt)`，返回新 `Waveform`。

## 3.2 控制线封装 `ControlLine`

[`sqc/hardware/control_line.py:19`](../sqc/hardware/control_line.py#L19) 把 `DistortionModel` 装到一条物理线上：

```python
@dataclass
class ControlLine:
    name: str                                # "Z0"
    kind: Literal["xy", "z", "readout"]
    source: str                              # "AWG0:CH0"
    target: str                              # "Q0"
    transfer_function: Optional["DistortionModel"] = None
    # 物理参数（可选）
    impedance: float = 50.0
    attenuation_db: float = 20.0
    delay: float = 0.0
```

两个核心方法：

- **`apply(awg_waveform)`** ([`control_line.py:65-99`](../sqc/hardware/control_line.py#L65-L99)) —— 正向模型：AWG → 片上。如果有 `transfer_function`，调 `apply_to_waveform()` 加失真；如果有 `delay > 0`，做样本级时移。
- **`predistort(target, designer)`** ([`control_line.py:101-120`](../sqc/hardware/control_line.py#L101-L120)) —— 反向：目标片上 → 所需 AWG，委托 `PredistortionDesigner`。

## 3.3 标定层 `WaveformCalibration` + `PredistortionDesigner`

### 3.3.1 `WaveformCalibration` —— 统一入口

[`sqc/calibration/waveform.py:27`](../sqc/calibration/waveform.py#L27) `WaveformCalibration` 有两种 method：

**A. `method="transfer_function"`** ([`waveform.py:103-119`](../sqc/calibration/waveform.py#L103-L119))：
1. 用 `distortion.step_response(t)` 直接得到阶跃响应
2. 用 `_fit_step_response()` 拟合成参数化模型
3. 返回 `CalibrationTable(kind="transfer_function", inputs=t, outputs=step, fit_params=...)`

`fit_type` 四种：

| fit_type | 拟合函数 | 物理模型 |
|----------|----------|----------|
| `"single_exp"` | `s(t) = 1 - amp·exp(-t/tau)` | `SingleExponentialDistortion` |
| `"multi_exp"` | `s(t) = 1 - Σ amp_k·exp(-t/tau_k)` | `MultiExponentialDistortion` |
| `"fir"` | impulse h ≈ ds/dt 截 N taps | `FIRDistortion` |
| `"iir"` | 1 阶 IIR 双线性 | `IIRDistortion` |

`_fit_multi_exp` ([`waveform.py:163-214`](../sqc/calibration/waveform.py#L163-L214))：**先串行剥离 K 次单指数**（每次拟合后从残差里减掉），再一次性 `curve_fit` 联合优化所有 2K 参数。这种 "greedy 初始化 + 联合精修" 比直接 K-参数初始化更稳。

`to_distortion_model()` ([`waveform.py:239-269`](../sqc/calibration/waveform.py#L239-L269))：把标定结果包成对应的 `DistortionModel` 子类。

**B. `method="predistortion"`** ([`waveform.py:275-296`](../sqc/calibration/waveform.py#L275-L296))：把 transfer_model 喂给 `PredistortionDesigner`，得到逆滤波器，包成 `CalibrationTable(kind="predistortion", fit_params={"inverse_model": ...})`。

### 3.3.2 `PredistortionDesigner` —— 三种逆滤波器设计

[`sqc/calibration/waveform.py:304`](../sqc/calibration/waveform.py#L304) `PredistortionDesigner`：

```python
@dataclass
class PredistortionDesigner:
    method: Literal["auto", "fir_inverse", "iir_inverse", "frequency_inverse"] = "auto"
    n_taps: int = 64
    regularization: float = 1e-4
```

`design(transfer_model, dt)` ([`waveform.py:326-351`](../sqc/calibration/waveform.py#L326-L351)) 根据 `method` 派遣到三条路径：

**A. `iir_inverse`（指数失真的解析逆）** ([`waveform.py:360-397`](../sqc/calibration/waveform.py#L360-L397))：

对 `SingleExponentialDistortion(amp, tau)`：原始 $H(s)=\dfrac{1+s\tau(1-\alpha)}{1+s\tau}$，**逆**为：

$$
H^{-1}(s) \;=\; \frac{1+s\tau}{1+s\tau(1-\alpha)}
$$

```python
# sqc/calibration/waveform.py:373-383
def _single_exp_to_iir_inverse(model, dt):
    amp, tau = model.amplitude, model.tau
    b_cont = [tau, 1.0]
    a_cont = [tau * (1.0 - amp), 1.0]
    b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
    return IIRDistortion(b_coeffs=b, a_coeffs=a)
```

对 `MultiExponentialDistortion`：直接用频域 Wiener 逆（[`waveform.py:385-397`](../sqc/calibration/waveform.py#L385-L397)），因为多个并联指数的 IIR 解析逆是高阶有理函数，不稳定。

**B. `fir_inverse`（FIR 频域 Wiener 逆）** ([`waveform.py:401-418`](../sqc/calibration/waveform.py#L401-L418))：

构造 δ 输入 → 跑正向 → 得经验 impulse response → FFT → Wiener 逆 → IFFT 截 N taps：

```python
n_fft = 4096
imp = np.zeros(n_fft); imp[0] = 1.0 / dt
h_fwd = transfer_model.apply(imp, dt)
H_emp = np.fft.fft(h_fwd)
H_inv = np.conj(H_emp) / (np.abs(H_emp) ** 2 + self.regularization ** 2)
h_inv_full = np.fft.ifft(H_inv).real
h_trunc = h_inv_full[:self.n_taps]
return FIRDistortion(taps=h_trunc)
```

**C. `frequency_inverse`（直接频域）** ([`waveform.py:422-435`](../sqc/calibration/waveform.py#L422-L435))：

调 `transfer_model.frequency_response(omega_grid)` 拿 H(ω)，Wiener 逆后包成 `CustomTransferDistortion`。

**`"auto"` 选择规则** ([`waveform.py:353-356`](../sqc/calibration/waveform.py#L353-L356))：类名含 "Exponential" 走 `iir_inverse`，否则走 `frequency_inverse`。

### 3.3.3 稳定性检查

[`waveform.py:459-465`](../sqc/calibration/waveform.py#L459-L465) `check_pole_stability(b, a)`：

$$
\text{所有} \;|p_i| < 1\;\text{（其中 }p_i\text{ 是 a(z) 的根）} \quad\Longleftrightarrow\quad \text{IIR 稳定}
$$

```python
@staticmethod
def check_pole_stability(b_coeffs, a_coeffs) -> bool:
    roots = np.roots(a_coeffs)
    return bool(np.all(np.abs(roots) < 1.0 - 1e-10))
```

—— 用于 IIR 解析逆设计后做 sanity check。

## 3.4 端到端工作流 `PredistortionValidationWorkflow`

[`sqc/workflows/predistortion_validation.py:28`](../sqc/workflows/predistortion_validation.py#L28) 把 §3.1–§3.3 串成 6 步：

```python
@dataclass
class PredistortionValidationWorkflow(Workflow):
    target_waveform: Waveform
    true_distortion: object       # DistortionModel
    designer: Optional[object] = None    # PredistortionDesigner, 默认 auto
    control_line_params: dict = field(default_factory=dict)
```

`run()` 流程（[`predistortion_validation.py:58-118`](../sqc/workflows/predistortion_validation.py#L58-L118)）：

1. **建控制线** —— `ControlLine(name, kind, ..., transfer_function=true_distortion)`
2. **无补偿测量** —— `on_chip_uncorrected = line.apply(target_waveform)`
3. **标定 H(ω)** —— `WaveformCalibration(distortion=true_distortion, method="simulation", fit_type=auto_inferred).calibrate()`
4. **设计逆** —— `inverse_model = designer.design(measured_model, dt)`
5. **预失真** —— `awg_predistorted = inverse_model.apply_to_waveform(target)`
6. **再次正向** —— `on_chip_corrected = line.apply(awg_predistorted)`
7. **指标对比** —— RMSE before/after + settling time before/after

`_infer_fit_type` ([`predistortion_validation.py:139-150`](../sqc/workflows/predistortion_validation.py#L139-L150))：根据真实失真类名自动选 fit_type。

指标计算 ([`predistortion_validation.py:156-195`](../sqc/workflows/predistortion_validation.py#L156-L195))：

```python
@staticmethod
def _rmse(measured, target):
    return float(np.sqrt(np.mean((measured.samples - target.samples) ** 2)))

@staticmethod
def _settling_time(w, target, tolerance=0.001):
    error = np.abs(w.samples - target.samples)
    bad = np.where(error > tolerance)[0]
    if len(bad) == 0: return 0.0
    return float(w.t_list[bad[-1]])
```

返回 dict 含：`target, on_chip_uncorrected, on_chip_corrected, awg_predistorted, inverse_model, measured_model, metrics`。

## 3.5 Z-线交叉串扰 `ZCrosstalkWorkflow`

[`sqc/workflows/z_crosstalk.py:25`](../sqc/workflows/z_crosstalk.py#L25) `ZCrosstalkWorkflow` —— 多 qubit 场景下，Z 控制线之间的电磁串扰使得"驱动 QA 时 QB 也感受到寄生磁通"。

物理：定义交叉传递矩阵 H_BA(ω) = Φ_B(ω) / V_A(ω)。**用 transient sensing 在 QB 上反演 φ_B(t)**，FFT 出 H_BA。

`run()` 流程（[`z_crosstalk.py:81-142`](../sqc/workflows/z_crosstalk.py#L81-L142)）：

1. **应用真实传递矩阵** —— `on_chip_fluxes = true_transfer_matrix.apply({QA: V_A})`，得到 `phi_A, phi_B_true`
2. **QB 上 transient sensing** ([`z_crosstalk.py:148-180`](../sqc/workflows/z_crosstalk.py#L148-L180))：
   - 把 QB 移到 sensitivity 最大点（`_make_qubit_at_optimal`, [`z_crosstalk.py:297-337`](../sqc/workflows/z_crosstalk.py#L297-L337)）
   - 跑 `TransientSensingExperiment(qubit=qubit_B_at_opt, flux_signal=phi_B_true)`
   - 用 `TransientReconstruction(method="wiener")` 反演 → `phi_B_reconstructed`
3. **频域反卷积 H_BA** ([`z_crosstalk.py:186-225`](../sqc/workflows/z_crosstalk.py#L186-L225))：

$$
\hat H_{BA}(\omega) \;=\; \frac{\hat\Phi_B(\omega)\cdot V_A^*(\omega)}{|V_A(\omega)|^2 + \lambda^2}
$$

```python
Phi_B_omega = np.fft.fft(phi_B_resampled)
V_A_omega = np.fft.fft(V_A_arr)
H_BA = Phi_B_omega * np.conj(V_A_omega) / (np.abs(V_A_omega) ** 2 + lam**2)
```

4. **真值对比** —— `H_BA_true = true_transfer_matrix.H_ji(QB, QA)`，计算 fit_error_dB
5. **补偿验证** ([`z_crosstalk.py:249-291`](../sqc/workflows/z_crosstalk.py#L249-L291))：
   - 寄生相位 ≈ ∫|φ_B_true| dt（未补偿）
   - 补偿后残差 φ_B_true − φ_B_reconstructed
   - compensation_factor = 未补偿 / 已补偿

返回 dict 含 `phi_A, phi_B_true, phi_B_reconstructed, H_BA_estimated, H_BA_true, fit_error_dB, parasitic_phase_uncompensated/compensated, compensation_factor, phi_B_after_compensation`。

## 3.6 Notebook 端到端示例（§10）

[`Simulation_sqc.ipynb`](../Simulation_sqc.ipynb) §10 是一套完整的预失真验证流程，5 个子章节：

### 3.6.1 §10.0 — 设置（cell 24）

```python
qubit = TransmonQubit(EC=0.2*2*np.pi, EJ=10.0*2*np.pi,
                     T1=100e3, T2=50e3, flux=0.0, n_levels=2)

# 故意夸大失真以可视化清晰（真实系统 amp ~ 1-10%）
dist = SingleExponentialDistortion(amplitude=0.3, tau=30.0)

line = ControlLine(name="Z0", kind="z", source="AWG0:CH0", target="Q0",
                    transfer_function=dist)
```

### 3.6.2 §10.1 — Rabi Chevron 失真诊断（cell 26）

跑两个 Chevron（detuning × pulse duration 的 P_e 二维图）：

- **理想**: 解析公式 `p_e = (Ω/Ω_eff)²·sin²(½·Ω_eff·t)`
- **失真**: 用 `dist.apply(Omega·ones, dt)` 得失真后的 Rabi 包络，再 mesolve

差分图 (b−a) 揭示失真在 |Δ|≤Ω 附近最显著。**附带证据**：累积 ∫Ω(t)dt vs t 的曲线在初期更平、晚期赶上——这正是阶跃响应的镜像。

### 3.6.3 §10.2 — 两种阶跃响应测量方法

- **§10.2.1 Cryoscope**（cell 29）：用 §1.2 的 pipeline 测 step。flux=0.25 敏感点 + step_amp=0.001 + tau=50ns。
- **§10.2.2 Transient**（cell 31）：用 §1.5 的 pipeline 测同一个 step。`lambda_reg=50.0`（增加正则避免 ringing）。
- **§10.2.3 对比**（cell 33）：可视化 RMSE、上升时间、采样代价。

### 3.6.4 §10.3 — 预失真滤波器设计（cell 35）

三种方法并跑：

```python
designer_iir  = PredistortionDesigner(method="iir_inverse")
designer_fir  = PredistortionDesigner(method="fir_inverse", n_taps=32, regularization=1e-3)
designer_freq = PredistortionDesigner(method="frequency_inverse", regularization=1e-4)

inv_iir  = designer_iir.design(dist, dt=dt)
inv_fir  = designer_fir.design(dist, dt=dt)
inv_freq = designer_freq.design(dist, dt=dt)
```

6 张图诊断：
- (a) |H(ω)| 原始
- (b) |H_inv(ω)| 三方法
- (c) 级联 |H_inv·H| ≈ 1 验证（IIR 在所有频段几乎完美 = 1，FIR 高频偏离）
- (d) 级联相位 ≈ 0 验证
- (e) FIR taps 图（看尾部能量分布）
- (f) 预失真后阶跃响应：理想 / 未校正 / IIR 校正 / FIR 校正

定量：未校正 RMSE → IIR 校正 RMSE 提升 ~50–100× (notebook 实测)。

### 3.6.5 §10.4 — 三层验证（cells 38, 40, 42）

**§10.4.1 阶跃响应再次 Cryoscope 测量**（cell 38）：

预失真 AWG → 失真线 → 片上 → Cryoscope 测出 h_cryo2，对比未补偿的 h_cryo。RMSE before/after 改善因子直接量化。

**§10.4.2 Chevron 复检**（cell 40）：

预失真后的驱动 Chevron 图 (b) vs 未补偿 (a) vs 理想 (c)：肉眼可见 (b) 几乎与 (c) 一致。

**§10.4.3 端到端工作流**（cell 42）：

```python
target = Waveform(t_list=t_wf,
                   samples=np.where((t_wf > 30) & (t_wf < 100), 1.0, 0.0))
wf = PredistortionValidationWorkflow(
    target_waveform=target,
    true_distortion=dist,
    designer=PredistortionDesigner(method="auto"),
)
result = wf.run()
m = result["metrics"]
print(f"RMSE: {m['rmse_uncorrected']:.6f} → {m['rmse_corrected']:.6f}")
print(f"Settling: {m['settling_uncorrected_ns']:.1f} → {m['settling_corrected_ns']:.1f} ns")
```

> notebook 实测的 RMSE 改善因子 ~30–100×，settling time 从 ~200 ns → 0 ns（在容差 1e-3 内）。

## 3.7 约束与陷阱

| 约束 | 出处 | 说明 |
|------|------|------|
| 单指数 IIR 逆解析最优 | §3.3.2.A | 但多指数时不稳定，必须走频域 |
| FIR 逆需 N_taps ≥ τ/dt | 物理 | 否则截断误差大，notebook 用 64 |
| `regularization` 决定噪声/精度权衡 | [`waveform.py:323`](../sqc/calibration/waveform.py#L323) | 默认 1e-4 偏激进；噪声大时调到 1e-3 |
| `check_pole_stability` 必须做 | §3.3.3 | 不稳 IIR 在某些 (b, a) 下会发散 |
| Predistortion 必须在 AWG dt 同一时基 | [`waveform.py:281`](../sqc/calibration/waveform.py#L281) | 否则采样不齐，加重 aliasing |
| `simulation` method == `transfer_function` | [`waveform.py:90`](../sqc/calibration/waveform.py#L90) | alias 保留，未来可能合并 |
| Z-crosstalk 要 QB 在 sensitivity 最大点 | [`z_crosstalk.py:299-321`](../sqc/workflows/z_crosstalk.py#L299) | 否则 Δω 太小测不出 |

---

# 附录

## 附录 A · 公式速查

| 量 | 公式 | 代码位置 |
|----|------|----------|
| Transmon f(Φ) | $\omega_Q = \sqrt{8E_J|\cos(\pi\Phi)|E_C} - E_C$ | [`src/qubit.py`](../src/qubit.py) `TransmonQubit.frequency` |
| Transmon Φ(ω) 反演 | $\Phi = \pi^{-1}\arccos\!\big[(f+E_C)^2/(8E_JE_C)\big]$ | [`reconstruction/dispersion.py:11`](../sqc/reconstruction/dispersion.py#L11) |
| Ramsey 相位 | $\varphi(\tau)=\int_0^\tau\Delta\omega(t)\,dt$ | [`hardware/readout.py:97`](../sqc/hardware/readout.py#L97) |
| IQ 相位提取 | $\varphi=\mathrm{atan2}(\tfrac12-p_e^I,p_e^Q-\tfrac12)$ | [`experiments/cryoscope.py:117`](../sqc/experiments/cryoscope.py#L117) |
| Cryoscope dφ/dt | $d\varphi/dt = \Delta\omega(h(t))$ | [`reconstruction/cryoscope.py:86`](../sqc/reconstruction/cryoscope.py#L86) |
| Delay Ramsey 斜率 | $\varphi_{\mathrm{cal}}(z)=\tau_R\kappa z$ | [`reconstruction/delay_ramsey.py:162`](../sqc/reconstruction/delay_ramsey.py#L162) |
| π pulse comp | $\Phi_{\mathrm{tail}}(\tau)=-z^*(\tau)$ | [`reconstruction/pi_pulse_comp.py:43`](../sqc/reconstruction/pi_pulse_comp.py#L43) |
| Wiener 反卷积 | $\hat B(\omega)=\frac{H^*\Delta P}{|H|^2+\lambda^2}$ | [`reconstruction/transient.py:149`](../sqc/reconstruction/transient.py#L149) |
| Hammerstein 反演 | Wiener → ω → 解析 Φ | [`reconstruction/transient.py:163`](../sqc/reconstruction/transient.py#L163) |
| LM 更新 | $\Delta b=(J^TJ+\mu I+\lambda R)^{-1}J^Tr$ | [`reconstruction/transient.py:461`](../sqc/reconstruction/transient.py#L461) |
| Adjoint Jacobian | $J_{ik}=\int\mathrm{Tr}[\lambda(t)[G,\rho]]\partial\Delta\omega/\partial\Phi\cdot\phi_k\,dt$ | [`reconstruction/transient.py:301`](../sqc/reconstruction/transient.py#L301) |
| Secant 法 | $V_{n+1}=V_n - r_n\frac{V_n-V_{n-1}}{r_n-r_{n-1}}$ | [`calibration/frequency.py:374`](../sqc/calibration/frequency.py#L374) |
| Bisection 误差 | $|r_n| \sim (V_b-V_a)/2^n$ | [`calibration/frequency.py:437`](../sqc/calibration/frequency.py#L437) |
| Transient frequency | $\Delta\omega = p_{\mathrm{diff}}/G_\alpha,\;G_\alpha=\int k(t)dt$ | [`calibration/frequency.py:567-581`](../sqc/calibration/frequency.py#L567) |
| 单指数失真 | $H(\omega)=(1-\alpha)+\alpha/(1+j\omega\tau)$ | [`hardware/distortion.py:177`](../sqc/hardware/distortion.py#L177) |
| 单指数解析逆 | $H^{-1}(s)=(1+s\tau)/(1+s\tau(1-\alpha))$ | [`calibration/waveform.py:373`](../sqc/calibration/waveform.py#L373) |
| 频域 Wiener 逆 | $H^{-1}=\frac{H^*}{|H|^2+\lambda^2}$ | [`calibration/waveform.py:432`](../sqc/calibration/waveform.py#L432) |
| Z-crosstalk H_BA | $H_{BA}(\omega)=\Phi_B(\omega)V_A^*/(|V_A|^2+\lambda^2)$ | [`workflows/z_crosstalk.py:223`](../sqc/workflows/z_crosstalk.py#L223) |

## 附录 B · 参考文献

1. **Gao et al. 2021** — *A Practical Guide for Building Superconducting Quantum Devices*. PRX Quantum 2, 040202. 文件: [`idea/refactor/Gao 等 - 2021 - Practical Guide for Building Superconducting Quantum Devices.pdf`](../idea/refactor/Gao%20等%20-%202021%20-%20Practical%20Guide%20for%20Building%20Superconducting%20Quantum%20Devices.pdf)
   - §III.B — Transmon 参数推荐范围
   - §III.D — 控制线传递函数
   - §V.B — Ramsey 相位累积
   - §V.E — Cryoscope distortion 标定
2. **Rol et al. 2020** — *Time-domain characterization and correction of on-chip distortion of control pulses in a quantum processor*. Appl. Phys. Lett. 116, 054001. (Cryoscope 原始论文)
3. **Vepsäläinen et al. 2022** — *Improving qubit coherence using closed-loop feedback*. Nat. Commun. 13, 1932. (Closed-loop frequency calibration 参考)
4. **Kelly et al. 2018** — *Physical qubit calibration on a directed acyclic graph*. arXiv:1803.03226. (Calibration DAG 设计参考；scheduler 的 `check_state/maintain/auto_calibrate` 是这套设计的 stub)

辅助资料：

- [`note/_sensing theory.md`](../note/_sensing%20theory.md) — 项目内部物理理论笔记，含瞬态磁场协议、核函数策略推导
- [`idea/refactor/_refactor_plan.md`](../idea/refactor/_refactor_plan.md) — sqc/ 框架重构主方案
- [`docs/architecture.md`](architecture.md) — 框架架构总览

## 附录 C · 测试与回归索引

### 单元测试（`tests/unit/`）

| 模块 | 测试文件 |
|------|----------|
| Cryoscope Experiment | [`test_cryoscope_experiment.py`](../tests/unit/test_cryoscope_experiment.py) |
| Delay Ramsey Experiment | [`test_delay_ramsey_experiment.py`](../tests/unit/test_delay_ramsey_experiment.py) |
| π-pulse Comp Experiment | [`test_pi_pulse_comp_experiment.py`](../tests/unit/test_pi_pulse_comp_experiment.py) |
| LM Reconstruction | [`test_lm_reconstruction.py`](../tests/unit/test_lm_reconstruction.py) |
| Basis Functions | [`test_basis_module.py`](../tests/unit/test_basis_module.py) |
| Kernel Estimator | [`test_kernel_estimator.py`](../tests/unit/test_kernel_estimator.py) |
| IQ Readout | [`test_iq_readout.py`](../tests/unit/test_iq_readout.py) |
| Calibration Table | [`test_calibration_table.py`](../tests/unit/test_calibration_table.py) |
| Distortion Models | [`test_distortion_models.py`](../tests/unit/test_distortion_models.py) |
| Predistortion Designer | [`test_predistortion_designer.py`](../tests/unit/test_predistortion_designer.py) |
| Control Line | [`test_control_line.py`](../tests/unit/test_control_line.py) |
| Transfer Matrix | [`test_transfer_matrix.py`](../tests/unit/test_transfer_matrix.py) |
| Chip Topology | [`test_chip_topology.py`](../tests/unit/test_chip_topology.py) |

### 回归测试（`tests/regression/`）

| 测试 | 内容 |
|------|------|
| [`test_physics_baseline.py`](../tests/regression/test_physics_baseline.py) | 物理回归 baseline (rtol=1e-6, atol=1e-9，**不可放宽**) |
| [`test_predistortion_baseline.py`](../tests/regression/test_predistortion_baseline.py) | 预失真端到端 baseline |
| [`test_z_crosstalk_baseline.py`](../tests/regression/test_z_crosstalk_baseline.py) | Z-crosstalk workflow baseline |
| [`generate_baselines.py`](../tests/regression/generate_baselines.py) | baseline 生成入口（修改 CONFIG 后必须重跑） |

### 运行测试

```bash
# 完整测试套件
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m pytest tests/ -v

# 仅回归
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m pytest tests/regression -m regression

# 单个模块
"C:\Users\21034\anaconda3\envs\qutip-env\python.exe" -m pytest tests/unit/test_cryoscope_experiment.py -v
```

## 附录 D · CONFIG 默认值（关键超参）

| 字段 | 默认值 | 影响 |
|------|--------|------|
| `CONFIG.awg.sample_rate` | 2.0 GSa/s | `dt = 0.5 ns` 派生 |
| `CONFIG.pulse.t_rabi_duration` | 10.0 ns | π 脉冲长度 |
| `CONFIG.pulse.t_global_end` | 400.0 ns | 仿真窗口 |
| `CONFIG.transmon.EC` | 0.2 GHz·(2π) | E_C/h = 200 MHz |
| `CONFIG.transmon.EJ` | 15.0 GHz·(2π) | E_J/h = 15 GHz |
| `CONFIG.reconstruction.lambda_reg` | 10.0 | Wiener 正则化 |
| `CONFIG.reconstruction.stim_amplitude` | 0.0215 | kernel 探针 |
| `CONFIG.reconstruction.stim_width` | 3.0 ns | kernel 探针宽度 |
| `CONFIG.reconstruction.lm_n_basis` | 100 | LM 基函数数 |
| `CONFIG.reconstruction.lm_basis_type` | `"fourier"` | LM 基类型 |
| `CONFIG.reconstruction.lm_lambda` | 100.0 | LM 平滑正则化 |
| `CONFIG.reconstruction.cryoscope_tau` | 100.0 ns | Cryoscope 标定 τ |
| `CONFIG.reconstruction.delay_ramsey_tau` | 20.0 ns | Delay Ramsey τ_R |
| `CONFIG.reconstruction.pi_pulse_T_pi` | 10.0 ns | π pulse 宽度 |

**修改 CONFIG 后**：必须重新生成 baseline（`tests/regression/generate_baselines.py`），并在 commit message 明确标注（CLAUDE.md R9）。

---

*文档版本 v1.0*
*生成日期 2026-05-14*
*维护：sqc/ 框架重构 Track A*
