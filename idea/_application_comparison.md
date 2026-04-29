# 量子传感协议的三大应用场景：方案设计与方法对比

## 0. 总览

本文档围绕仿真平台现有能力，系统性地将量子传感协议映射到三个核心应用场景，并与主流实验方法进行对比分析。

| 应用场景 | 核心问题 | 主流方法 | 本平台方法 | 关键区别 |
|----------|----------|----------|------------|----------|
| **频率标定** | 确定 qubit 工作频率 $f_Q(\Phi)$ | 标准 Ramsey | Ramsey (case 1) | 基本等价，可扩展到非线性标定 |
| **波形重建** | 还原片上实际磁通波形 $\Phi(t)$ | Cryoscope（截断 Ramsey） | 滑动 Ramsey + 反卷积 (case 4) | 信息提取方式不同 |
| **波形预失真** | 设计补偿滤波器 $h_{\text{filt}}$ | Cryoscope + IIR/FIR 滤波 | 重建 + 逆滤波设计 | 标定-补偿-验证闭环 |

---

## 1. 应用一：Qubit 频率标定

### 1.1 问题定义

确定 Transmon qubit 的跃迁频率 $f_{01}$ 及其随磁通的依赖关系 $f_Q(\Phi)$，包括：

- **点频标定**：在给定工作点 $\Phi_w$ 处精确确定 $f_{01}$
- **磁通-频率映射**：扫描 $\Phi$，建立完整的 $f_Q(\Phi)$ 曲线
- **灵敏度标定**：确定 $\kappa = \partial f_Q / \partial \Phi|_{\Phi_w}$

### 1.2 主流方法：标准 Ramsey 干涉

#### 原理

1. 设置驱动频率 $\omega_d$ 略偏离 $f_{01}$
2. 执行 Ramsey 序列：$(\pi/2)_X \to \text{free}(\tau) \to (\pi/2)_X$
3. 扫描 $\tau$，观察振荡：$p_e(\tau) = \frac{1}{2}(1 - \cos(\Delta \cdot \tau))$
4. 拟合振荡频率得 $\Delta = f_{01} - \omega_d/(2\pi)$

#### 实验流程

```
粗标定: 频谱扫描 → 找到共振峰 → 确定 f_01 ± 10 MHz
↓
精标定: Ramsey 干涉 → 拟合 Δ → f_01 精度 ~ kHz
↓
磁通映射: 扫描 Φ，每个点做 Ramsey → f_Q(Φ) 曲线
```

### 1.3 本平台实现

#### 已有能力

平台的 **case 1**（Ramsey 测量）已经实现了标准 Ramsey 干涉：

```python
# protocal.py case 1
# 创建 Ramsey 脉冲序列，扫描 tau，测量 p_e(tau)
control_pulse = create_ramsey_pulse(t_rabi, tau, omega_d, phase1, phase2)
# → 返回 tau_list, p_e_list
```

`analysis.py` 中的相位提取方法可以从 $p_e(\tau)$ 还原出 $\varphi(\tau)$：
- `get_signal_from_ramsey_by_iq()`：IQ 双通道提取（精确）
- `get_signal_from_ramsey_by_unwrap()`：单通道 Viterbi 解缠绕

#### 频率标定的具体实现

**点频标定**：恒定偏置下的 Ramsey

```python
# 设置 qubit 在工作点 Φ_w
qubit = TransmonQubit(EC, EJ, T1, T2, flux=Phi_w, n_levels=2)
omega_d = qubit.frequency - detuning  # 设置一个已知的小失谐

# 执行标准 Ramsey（case 1，但信号为零）
Phi_zero = Signal(type=0, t_list=tau_list)  # 零信号
# Ramsey 序列下扫描 tau
# p_e(tau) = 0.5 * (1 - cos(Δ·tau))
# 拟合得 Δ → f_01 = omega_d + Δ
```

**磁通-频率映射**：在标定步中扫描 Φ

```python
# 扫描外加磁通
Phi_list = np.linspace(-0.3, 0.3, 61)
f_list = []
for Phi in Phi_list:
    qubit.change_flux(Phi)
    # 做 Ramsey → 拟合得 f_01(Phi)
    f_list.append(f_01_fitted)
# f_list 即为 f_Q(Φ) 曲线
```

**灵敏度标定**：对 $f_Q(\Phi)$ 曲线做数值微分

$$\kappa(\Phi_w) = \frac{f_Q(\Phi_w + \delta\Phi) - f_Q(\Phi_w - \delta\Phi)}{2\delta\Phi}$$

### 1.4 对比分析

| 对比维度 | 标准 Ramsey | 本平台 (case 1) |
|----------|------------|----------------|
| 原理 | 完全一致 | 完全一致 |
| 相位提取 | 通常拟合 $\cos$ 振荡 | 支持 IQ 双通道和 Viterbi 解缠绕 |
| 精度 | 受 $T_2$ 限制 | 同上，额外受数值精度影响 |
| 磁通映射 | 逐点 Ramsey | 逐点 Ramsey |
| 非线性标定 | 需要额外分析 | `frequency_sensitivity()` 已实现 |

**结论**：对于频率标定，本平台的 case 1 与实验完全等价。平台的额外优势在于可以方便地研究非线性标定（即 $f_Q(\Phi)$ 的高阶项对后续波形重建的影响），以及 IQ 双通道精确相位提取。

### 1.5 可扩展方向

1. **Cryoscope 标定曲线**：固定 $T_{\text{bias}}$，扫描 $h$，建立 $\varphi(h)$ 曲线——这实际上就是 $f_Q(\Phi)$ 的一种采样方式（见 `cryoscope_implementation.md`）
2. **非线性映射校正**：当工作在甜点附近时，$f_Q(\Phi) \propto \Phi^2$，标准线性反演会产生误差。可利用 `hammerstein_wiener_deconvolution()` 进行非线性校正

---

## 2. 应用二：磁通波形重建

### 2.1 问题定义

已知 AWG 输出了一个理想波形 $V_{\text{in}}(t)$，由于控制线路中各种元件（偏置器、低通滤波器、阻抗失配、趋肤效应等）引入的线性动力学失真，到达芯片的实际磁通波形为：

$$\Phi_Q(t) = h * V_{\text{in}}(t)$$

其中 $h$ 是系统冲激响应。目标：从量子测量中还原 $\Phi_Q(t)$。

### 2.2 主流方法：Cryoscope

#### 2.2.1 核心思想

Cryoscope 利用 **截断 Ramsey**（truncated Ramsey）实验，以 AWG 的时间分辨率逐点采样片上磁通波形。

#### 2.2.2 脉冲序列

```
(π/2)_Y → 待测波形（截断到 τ） → free(T_sep - τ) → (π/2)_{X,Y} → 读出
```

对于两个相邻截断点 $\tau$ 和 $\tau + \Delta\tau$：

$$\overline{\Delta f}_R(\tau) = \frac{\varphi_{\tau+\Delta\tau} - \varphi_\tau}{2\pi \cdot \Delta\tau}$$

这直接给出了 $\tau$ 时刻的瞬时频率偏移，进而通过 $f_Q(\Phi)$ 的反函数还原磁通。

#### 2.2.3 关键特性

- **时间分辨率**：等于 AWG 采样间隔 $\Delta\tau = 1/f_{\text{AWG}}$
- **精度**：阶跃响应可达 $\pm 0.1\%$
- **无需反卷积**：差分操作直接给出瞬时值
- **误差来源**：关断瞬态的差异，但被甜点处的二次非线性强烈抑制（$|\varepsilon|/\overline{\Delta f}_R \lesssim 10^{-2} \sim 10^{-3}$）

#### 2.2.4 信号处理流程

```
IQ 数据 → arctan2 → φ(τ) → demod → unwrap → SG 微分 → Δf(τ) → f_Q⁻¹ → Φ(τ)
```

### 2.3 本平台方法：滑动 Ramsey + 反卷积 (case 4)

#### 2.3.1 核心思想

使用固定结构的 Ramsey 脉冲序列，通过改变脉冲序列与磁通信号之间的 **时间延迟** $t_{\text{delay}}$，进行等效时间采样。测量结果是磁场波形与脉冲核函数的卷积，需要通过反卷积还原。

#### 2.3.2 脉冲序列

```
Ramsey: (π/2)_X → free(τ) → (π/2)_X
       ↕
信号:  ← t_delay → Φ(t) 沿时间轴平移
```

#### 2.3.3 数学模型

单次测量的激发态概率变化为：

$$\delta p_e(t_{\text{delay}}) = \int K(t - t_{\text{delay}}) \cdot \kappa \cdot \Phi(t) \, dt$$

其中 $K(t)$ 是 Ramsey 脉冲序列的核函数（由脉冲序列的调制函数决定）。这是一个卷积关系：

$$\delta p_e = K * (\kappa \Phi)$$

#### 2.3.4 波形还原流程

```
方法1 (线性): 滑动测量 → δp_e(t_delay) → Wiener 反卷积 → Φ(t)
方法2 (非线性): 滑动测量 → Hammerstein-Wiener 反卷积 → Φ(t)
方法3 (全局优化): 滑动测量 → p_e(t_delay) → LM 数值反演 → Φ(t)
```

### 2.4 核心区别与对比

#### 2.4.1 信息提取方式

| 维度 | Cryoscope | 本平台 case 4 |
|------|-----------|--------------|
| **采样方式** | 截断波形，扫描截断点 $\tau$ | 固定脉冲，扫描延迟 $t_{\text{delay}}$ |
| **原始数据** | $\varphi(\tau) = \int_0^\tau \Delta f(t) dt$ （相位积分） | $\delta p_e(t_d) = \int K(t-t_d) \kappa\Phi(t) dt$ （卷积） |
| **信息维度** | 每点携带 $[0, \tau]$ 的完整相位历史 | 每点仅反映脉冲窗口内的加权平均 |
| **波形还原** | 差分 → 瞬时值，**无需反卷积** | **需要反卷积或优化** |
| **核函数** | 无（微分即可） | 需要标定核函数 $K(t)$ |

**物理本质的差异**：

- **Cryoscope** 测量的是波形的 **积分**（从 0 到 $\tau$ 的相位积累），通过微分得到瞬时值。每个数据点包含从 $t=0$ 到截断点的全部历史信息。

- **Case 4** 测量的是波形与核函数的 **卷积**（滑动窗口加权平均）。每个数据点仅包含脉冲窗口附近的局部信息。

#### 2.4.2 性能对比

| 性能指标 | Cryoscope | 本平台 case 4 |
|----------|-----------|--------------|
| **时间分辨率** | $\Delta\tau$（AWG 步长） | $\sim$ 脉冲宽度（$\tau_{\text{Ramsey}}$ 决定） |
| **精度** | 0.1%（实验验证） | 依赖反卷积正则化 |
| **噪声放大** | SG 微分有一定放大，但可控 | 反卷积放大高频噪声 |
| **非线性处理** | 通过 $f_Q^{-1}$ 完整处理 | Hammerstein-Wiener 或 LM |
| **计算复杂度** | 低（微分 + 查表） | 高（反卷积/优化迭代） |
| **对脉冲形状的要求** | 只需精确的 $\pi/2$ 脉冲 | 需要已知的精确核函数 |
| **适用工作点** | 甜点最优（二次抑制瞬态误差） | 最大灵敏度点最优 |

#### 2.4.3 优劣势总结

**Cryoscope 的优势**：
1. 无需反卷积，避免了噪声放大和正则化参数选择
2. 时间分辨率由 AWG 硬件决定，不受脉冲宽度限制
3. 在甜点工作时误差被自然抑制
4. 信号处理流程标准化，易于自动化

**Case 4（滑动 Ramsey）的优势**：
1. 脉冲序列结构固定，实验实现更简单（不需要逐点改变波形截断）
2. 核函数框架具有一般性——更换脉冲序列（Echo, CPMG, 差分回波）等价于更换核函数，整个反卷积框架不变
3. LM 数值反演可以自然处理非线性效应
4. 可以利用不同脉冲序列的频率选择性（如 CPMG 的带通特性）优化特定频段的重建

**Case 4 的劣势**：
1. 反卷积是病态问题，需要精心选择正则化参数
2. 核函数标定本身需要额外实验
3. 分辨率受限于核函数宽度

### 2.5 整合方案：将两种方法统一到同一框架

两种方法可以统一理解为 **不同的调制函数** 下的相位积累：

$$\varphi = \int y(t') \cdot \Delta\omega(t') \, dt'$$

| 方法 | 调制函数 $y(t')$ | 测量量 |
|------|-----------------|--------|
| Cryoscope | $y(t') = \mathbb{1}_{[0, \tau]}(t')$（截断窗口） | $\varphi(\tau) = \int_0^\tau \Delta\omega(t')dt'$ |
| 滑动 Ramsey | $y(t') = y_{\text{Ramsey}}(t' - t_d)$（平移核） | $\delta p_e(t_d) \propto \int y(t'-t_d) \Delta\omega(t') dt'$ |
| 差分回波 | $y(t') = M(t, t')$（差分调制） | $p_e(t) \propto B(t)$ |
| CPMG | $y(t') = y_{\text{CPMG}}(t')$（带通滤波） | 频域采样 |

**本平台的独特价值**：能够在同一仿真框架下比较所有这些方法，这是单一实验方法无法实现的。

### 2.6 本平台实现 Cryoscope 的路径

在现有框架基础上实现 Cryoscope 非常自然——它本质上是 case 1（Ramsey）的变体，只是：

1. 信号从"外部未知信号"变为"已知但截断的待测波形"
2. 扫描变量从 $\tau$（等待时间）变为截断点位置
3. 相位提取方式相同（IQ 或 unwrap）

具体实现见 `cryoscope_implementation.md`。

---

## 3. 应用三：波形预失真

### 3.1 问题定义

控制线路的传递函数 $H(\omega)$（或等价地，冲激响应 $h(t)$）导致实际到达芯片的波形偏离理想。目标是设计预失真滤波器 $h_{\text{filt}}$，使得：

$$h_{\text{filt}} * h * V_{\text{in}}(t) \approx V_{\text{in}}(t)$$

即 $h_{\text{filt}} \approx h^{-1}$。

### 3.2 主流方法：Cryoscope + IIR/FIR 迭代校正

#### 3.2.1 工作流程（Rol et al. 2020）

```
第1步: Cryoscope 测量阶跃响应 s(t)
         ↓
第2步: 设计 IIR 滤波器 → 校正慢动态 (>30 ns)
       s(t) = g·(1 + A·exp(-t/τ_IIR)) · u(t)
       优化 A, τ_IIR, g → 3~5 个 IIR 段
         ↓
第3步: Cryoscope 验证 IIR 效果 → 慢动态校正到 ±0.1%
         ↓
第4步: 设计 FIR 滤波器 → 校正快动态 (<30 ns)
       72 个系数 (30 ns), CMA-ES 优化
         ↓
第5步: Cryoscope 最终验证 → 全时间尺度 ±0.1%
         ↓
第6步: Chevron 实验独立验证
```

#### 3.2.2 IIR 滤波器设计

每个 IIR 段校正一个指数衰减分量：

$$s_{\text{corrected}}(t) = s(t) \cdot \frac{1}{1 + A \cdot e^{-t/\tau}}$$

对应的 z 变换 IIR 滤波器：

$$H_{\text{IIR}}(z) = \frac{1 - \alpha z^{-1}}{1 - \beta z^{-1}}$$

其中 $\alpha, \beta$ 由 $A, \tau$ 确定。通常需要 3–5 个这样的 IIR 段级联。

#### 3.2.3 FIR 滤波器设计

FIR 滤波器由 72 个系数描述（对应 30 ns），参数化为 40 个变量，通过 CMA-ES 全局优化算法求解：

$$\min \| h_{\text{FIR}} * h_{\text{IIR}} * s(t) - u(t) \|$$

### 3.3 本平台实现方案

#### 3.3.1 方案 A：基于 Cryoscope 的标准流程

直接在仿真中复现 Rol et al. 的流程。这是最直接的方案，实现路径已在 `cryoscope_implementation.md` 中描述。

```
1. 在仿真中引入失真模型（指数衰减、带限等）
2. 用 Cryoscope (case 7) 测量阶跃响应
3. 拟合 IIR 参数，设计补偿滤波器
4. 施加补偿后再次用 Cryoscope 验证
5. 迭代直到满足精度要求
```

**优势**：与实验流程完全一致，仿真结果可直接指导实验参数选择。

#### 3.3.2 方案 B：基于 Case 4 滑动测量的间接方法

利用现有的 case 4 框架，通过以下路径实现预失真：

```
1. 施加理想阶跃波形（包含失真模型）
2. 用 case 4 滑动测量 → Wiener 反卷积重建实际波形
3. 计算传递函数：H(ω) = FFT(重建波形) / FFT(理想波形)
4. 设计逆滤波器：H_filt(ω) = 1/H(ω)（带正则化）
5. 施加预失真后的波形再次测量验证
```

**优势**：利用现有的反卷积基础设施，无需额外实现 Cryoscope。

**劣势**：反卷积本身引入误差，预失真的精度受限于重建精度。

#### 3.3.3 方案 C：全局数值优化方法

利用 `analysis.py` 中的 LM 数值反演框架，直接优化预失真波形：

```
1. 定义目标：施加 V_predist(t) 后，测量结果应等于理想波形的预期测量结果
2. 参数化 V_predist(t) 为基函数展开
3. 正向模拟：V_predist → 失真 → 片上波形 → 量子测量 → p_e
4. 目标函数：min ||p_e_simulated - p_e_target||²
5. LM 优化迭代
```

**优势**：
- 自动处理非线性效应
- 不需要分离 IIR/FIR 两步
- 可以直接优化门保真度而非波形形状

**劣势**：
- 计算量大
- 可能陷入局部最优
- 对实验不太实用（太慢）

### 3.4 失真模型设计

为了在仿真中测试预失真方案，需要引入可控的失真模型。建议在 `signal.py` 中新增失真功能：

```python
class DistortionModel:
    """控制线路失真模型"""
    
    def __init__(self, distortion_type='exponential', **params):
        """
        distortion_type:
            'exponential': s(t) = 1 + A*exp(-t/τ)  （指数过冲/欠冲）
            'lowpass':     H(ω) = 1/(1 + jωτ)      （一阶低通）
            'bandlimit':   H(ω) = sinc(ω/ω_bw)     （带限）
            'multi_exp':   s(t) = 1 + Σ A_i*exp(-t/τ_i)  （多指数）
            'custom':      用户自定义传递函数
        """
        self.type = distortion_type
        self.params = params
    
    def apply(self, signal):
        """对信号施加失真"""
        if self.type == 'exponential':
            A = self.params.get('A', 0.1)      # 过冲幅度
            tau = self.params.get('tau', 40)    # 时间常数 (ns)
            t = signal.t_list
            kernel = np.zeros_like(t)
            kernel[0] = 1
            for i in range(1, len(t)):
                dt = t[i] - t[i-1]
                kernel[i] = A/tau * np.exp(-t[i]/tau)
            distorted = np.convolve(signal.signal, kernel * (t[1]-t[0]), 
                                     mode='full')[:len(t)]
            return distorted
        
        elif self.type == 'multi_exp':
            A_list = self.params.get('A_list', [0.05, 0.03, 0.02])
            tau_list = self.params.get('tau_list', [30, 100, 300])
            t = signal.t_list
            dt = t[1] - t[0]
            distorted = signal.signal.copy()
            for A, tau in zip(A_list, tau_list):
                decay = A * np.exp(-t / tau)
                correction = np.convolve(signal.signal, decay * dt, 
                                          mode='full')[:len(t)]
                distorted += correction
            return distorted
    
    def get_step_response(self, t_list):
        """返回理论阶跃响应"""
        if self.type == 'exponential':
            A = self.params.get('A', 0.1)
            tau = self.params.get('tau', 40)
            return 1 + A * np.exp(-t_list / tau)
        elif self.type == 'multi_exp':
            s = np.ones_like(t_list)
            for A, tau in zip(self.params['A_list'], self.params['tau_list']):
                s += A * np.exp(-t_list / tau)
            return s
    
    def design_iir_correction(self, s_measured, t_list, n_stages=3):
        """
        从测量的阶跃响应设计 IIR 校正滤波器
        拟合 s(t) = g * (1 + Σ A_i * exp(-t/τ_i)) 
        """
        from scipy.optimize import curve_fit
        
        def model(t, *params):
            g = params[0]
            result = np.ones_like(t) * g
            for i in range((len(params)-1)//2):
                A_i = params[1 + 2*i]
                tau_i = params[2 + 2*i]
                result += g * A_i * np.exp(-t / tau_i)
            return result
        
        # 初始猜测
        p0 = [1.0]
        for i in range(n_stages):
            p0.extend([0.05, 30 * (i+1)])  # A, tau 初始值
        
        popt, _ = curve_fit(model, t_list, s_measured, p0=p0, maxfev=10000)
        
        return popt  # [g, A1, τ1, A2, τ2, ...]
```

### 3.5 对比分析

| 维度 | Cryoscope + IIR/FIR | Case 4 + 逆滤波 | LM 全局优化 |
|------|---------------------|-----------------|-------------|
| **精度** | 0.1%（实验验证） | 取决于反卷积质量 | 理论上最优 |
| **计算量** | 低 | 中 | 高 |
| **非线性处理** | 通过 $f_Q^{-1}$ 查表 | Hammerstein-Wiener | 全密度矩阵自动处理 |
| **迭代次数** | 2-3 次 | 需要更多迭代 | 10-50 次 |
| **实验可行性** | ✓ 已验证 | ✓ 可行 | ✗ 太慢 |
| **适用范围** | 线性失真 | 线性失真 | 线性+非线性 |

---

## 4. 统一仿真框架设计

### 4.1 三个应用的代码映射

```
应用一：频率标定
├── protocal.py case 1 (Ramsey)     → 点频标定
├── protocal.py case 6 (Cryoscope标定) → φ(h) 曲线 = 频率-磁通映射
└── analysis.py frequency_sensitivity() → κ 标定

应用二：波形重建
├── protocal.py case 4 (滑动Ramsey)  → 卷积测量 + Wiener/LM 反演
├── protocal.py case 7 (Cryoscope)   → 截断Ramsey + 微分反演
├── protocal.py case 2 (差分回波)    → 差分相位 + 直接反演
└── analysis.py 各种反演方法

应用三：波形预失真
├── signal.py DistortionModel       → 引入失真 (新增)
├── 方案A: case 7 测量 → IIR/FIR 设计 → case 7 验证
├── 方案B: case 4 测量 → 逆滤波设计 → case 4 验证
└── 方案C: numerical_inverse() 直接优化预失真波形
```

### 4.2 推荐实施路线

```
Phase 1: 基础建设（已完成 ✓）
├── [✓] Ramsey 协议 (case 1)
├── [✓] 滑动测量 (case 4)
├── [✓] 差分回波 (case 2)
├── [✓] Wiener 反卷积
├── [✓] LM 数值反演
└── [✓] 核函数标定

Phase 2: Cryoscope 实现
├── [○] case 6: 标定协议
├── [○] case 7: 截断波形测量
├── [○] 分析方法: SG微分 + 标定反查
└── [○] 失真模型类

Phase 3: 应用对比
├── [○] 频率标定: Ramsey vs Cryoscope 标定
├── [○] 波形重建: Cryoscope vs case 4 vs 差分回波
│   ├── 精度对比 (不同信号类型)
│   ├── 噪声鲁棒性对比
│   ├── 时间分辨率对比
│   └── 计算效率对比
└── [○] 预失真: IIR/FIR vs 逆滤波 vs LM
    ├── 校正精度
    ├── 迭代收敛速度
    └── Chevron 验证
```

### 4.3 对比实验设计

#### 实验 1: 频率标定精度

| 测试条件 | 方法 A: Ramsey 拟合 | 方法 B: Cryoscope 标定曲线 |
|----------|--------------------|-----------------------|
| 恒定偏置 h=0.01 Φ₀ | 标准 cos 拟合 → Δf | φ(h) 查表 → Δf |
| 扫描 h, 建立 f(Φ) 曲线 | 逐点 Ramsey 拟合 | 单次标定扫描 |
| **指标**: 与理论值的偏差 | | |

#### 实验 2: 波形重建精度

| 测试信号 | Cryoscope | Case 4 + Wiener | Case 4 + LM | 差分回波 |
|----------|-----------|----------------|-------------|---------|
| 高斯脉冲 (type=3) | ✓ | ✓ | ✓ | ✓ |
| 双指数脉冲 (type=4) | ✓ | ✓ | ✓ | ✓ |
| 方波 (type=1) | ✓ | ✓ | ✓ | ✓ |
| 复杂波包 (type=7) | ✓ | ✓ | ✓ | ✓ |
| **指标**: RMSE, 峰值误差, 时间分辨率 | | | | |

#### 实验 3: 预失真效果

| 失真类型 | Cryoscope + IIR/FIR | Case 4 + 逆滤波 | LM 直接优化 |
|----------|--------------------|--------------|-----------| 
| 单指数 (A=0.1, τ=40 ns) | ✓ | ✓ | ✓ |
| 多指数 (3 段) | ✓ | ✓ | ✓ |
| 带限 (BW=500 MHz) | ✓ | ✓ | ✓ |
| **指标**: 校正后残差, 迭代次数, 计算时间 | | | |

---

## 5. 本平台的独特价值

### 5.1 相比纯实验

1. **方法可控性**：可以精确控制失真类型和幅度，系统性地研究各方法的适用边界
2. **真值可知**：仿真中原始波形已知，可以计算精确误差
3. **无实验限制**：可以测试极端参数（如非常大的失真、非常短的脉冲等）
4. **方法对比**：在完全相同的条件下对比不同方法，消除实验系统差异

### 5.2 相比纯理论

1. **全密度矩阵模拟**：考虑了 $T_1$、$T_2$ 退相干，不是理想化的 Bloch 方程
2. **非线性效应**：Transmon 的 $f_Q(\Phi)$ 非线性和多能级泄露效应被自然包含
3. **有限脉冲宽度**：$\pi/2$ 脉冲不是理想的瞬间操作，有限宽度效应被正确模拟
4. **数值反演验证**：可以验证各种反演算法（Wiener、LM、直接公式）的实际表现

### 5.3 创新点

1. **统一框架**：首次在同一仿真平台上实现并对比 Cryoscope、滑动 Ramsey、差分回波、CPMG 四种波形重建方法
2. **核函数视角**：将所有方法统一为不同调制函数下的相位积累，建立方法间的数学联系
3. **非线性反演**：LM + 伴随法雅可比矩阵的全密度矩阵反演方法可以处理大信号情况，这是实验中常见但理论中常忽略的
4. **闭环验证**：标定 → 预失真 → 验证的完整闭环，可以量化每一步的误差传递

---

## 6. 下一步工作建议

### 优先级排序

1. **[高]** 实现 Cryoscope 协议 (case 6, 7)：这是对比研究的核心基准
2. **[高]** 实现失真模型类 `DistortionModel`：三个应用都需要
3. **[中]** 对比实验 2 (波形重建)：最有学术价值的对比
4. **[中]** 预失真闭环验证 (实验 3)
5. **[低]** 频率标定对比 (实验 1)：差异较小，优先级低
6. **[低]** CPMG 频域重建 (case 3 实现)：补充频域视角

### 预期成果

- 一组系统性的方法对比数据
- 各方法在不同信号类型、不同噪声水平、不同工作点下的适用条件总结
- 预失真设计的最优工作流程推荐
- （可选）发表一篇仿真对比综述文章
