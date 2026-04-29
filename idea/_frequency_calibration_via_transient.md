# 瞬态磁场协议用于 Qubit 频率标定

## 1. 核心思想

Case 4 的滑动 Ramsey 测量模型：

$$\delta p_e(t_d) = \int K(t - t_d) \cdot \delta\omega(\Phi(t)) \, dt$$

标准应用中 $\Phi(t)$ 未知、$\delta\omega(\Phi) = \kappa \Phi$（线性近似），目标是还原 $\Phi(t)$。

**频率标定的逆问题**：施加 **已知信号** $\Phi(t)$，目标变为提取 **未知的频率响应** $\delta\omega(\Phi)$。

这是一个从"传感"到"标定"的视角翻转：同一个实验框架，仅仅交换了已知量和未知量。

---

## 2. 三种方案

### 2.1 方案一：恒定信号幅度扫描（朴素方案）

#### 原理

施加一系列恒定信号 $\Phi(t) = h_i$，对每个 $h_i$ 做滑动测量。

当 Ramsey 核函数完全落在信号区域内时（$t_d$ 远离边界）：

$$\delta p_e(t_d) \approx \delta\omega(h_i) \cdot \underbrace{\int K(t) dt}_{K_0} = \delta\omega(h_i) \cdot K_0$$

扫描 $h_i \in [h_{\min}, h_{\max}]$，逐点提取：

$$\delta\omega(h_i) = \frac{\delta p_e^{\text{plateau}}(h_i)}{K_0}$$

#### 评价

- **等价于**逐点 Ramsey + 拟合，没有本质优势
- 每个 $h_i$ 需要一次完整滑动扫描（冗余，因为只用了 plateau 区域的值）
- 唯一优势：可以顺便看到信号上升/下降沿处的瞬态响应

### 2.2 方案二：斜坡信号——单次扫描非线性标定（推荐）

#### 原理

施加线性斜坡信号：

$$\Phi(t) = \Phi_{\min} + \alpha \cdot t, \quad \alpha = \frac{\Phi_{\max} - \Phi_{\min}}{T}$$

斜坡将 **时间轴映射到磁通轴**：$t \leftrightarrow \Phi$。因此滑动测量扫描 $t_d$ 时，Ramsey 窗口在不同的磁通值处采样：

$$\delta p_e(t_d) = \int K(t - t_d) \cdot \delta\omega(\Phi_{\min} + \alpha \cdot t) \, dt$$

做变量替换 $\Phi = \Phi_{\min} + \alpha t$：

$$\delta p_e(t_d) = \frac{1}{\alpha} \int K\left(\frac{\Phi - \Phi_d}{\alpha}\right) \cdot \delta\omega(\Phi) \, d\Phi$$

其中 $\Phi_d = \Phi_{\min} + \alpha \cdot t_d$。

这依然是一个卷积关系，但卷积核变成了 $\tilde{K}(\Phi) = K(\Phi/\alpha)/\alpha$（在磁通域上的核），待还原函数变成了 $\delta\omega(\Phi)$。

#### 反卷积

直接复用现有的 Wiener 反卷积框架：

$$\delta\omega(\Phi) = \mathcal{F}^{-1}\left[\frac{\tilde{K}^*(\omega)}{\left|\tilde{K}(\omega)\right|^2 + \lambda^2} \cdot \delta \tilde{P}_e(\omega)\right]$$

其中所有变换在磁通域上进行。

#### 关键优势

| 维度 | 标准 Ramsey 标定 | 斜坡 + 滑动 Ramsey |
|------|-----------------|-------------------|
| 测量次数 | $M \times N$（M 个磁通点 × N 个 $\tau$ 点） | $N$（单次滑动扫描） |
| 得到什么 | 逐点 $f_Q(\Phi_i)$ | 连续曲线 $\delta\omega(\Phi)$ |
| 非线性 | 每点独立，自动包含 | 反卷积可处理 |
| 磁通分辨率 | 由扫描步长决定 | 由核函数宽度和斜坡斜率共同决定 |

#### 磁通分辨率

在磁通域上的有效分辨率：

$$\delta\Phi \approx \alpha \cdot \Delta t_K$$

其中 $\Delta t_K$ 为 Ramsey 核函数的等效宽度（约等于 Ramsey 自由演化时间 $\tau_{\text{Ramsey}}$）。

选择 $\alpha$ 时需要平衡：
- $\alpha$ 大 → 磁通范围大，但分辨率低
- $\alpha$ 小 → 分辨率高，但范围小

#### 实际约束

1. **斜坡斜率上限**：频率变化速率不能太快，否则 Ramsey 在一个核函数宽度内经历的频率变化过大，模糊了信息
   
   $$\alpha \cdot \Delta t_K \cdot \kappa_{\max} \ll \frac{1}{\Delta t_K} \quad \Rightarrow \quad \alpha \ll \frac{1}{\kappa_{\max} \cdot \Delta t_K^2}$$

2. **斜坡斜率下限**：需要在可用时间 $T \lesssim T_2$ 内覆盖目标磁通范围

3. **线性近似失效**：当信号幅度大时，$\delta\omega(\Phi)$ 的非线性意味着核函数实际上不再是常数——这恰好就是我们要标定的非线性效应

### 2.3 方案三：Hammerstein-Wiener 自洽标定

#### 原理

利用已有的 Hammerstein-Wiener 框架（`analysis.py`），但反转标定和测量的角色：

**标准用法**（波形重建）：
```
已知: K(t), f_Q(Φ) 关系
未知: Φ(t)
模型: Φ → δω(Φ) [非线性] → K*δω [线性卷积] → δp_e
```

**标定用法**（频率标定）：
```
已知: K(t), Φ(t)  (施加已知信号)
未知: δω(Φ) 关系
模型: Φ(t) → δω(?) [待标定] → K*δω [线性卷积] → δp_e
```

将 $\delta\omega(\Phi)$ 参数化为多项式展开：

$$\delta\omega(\Phi) = \sum_{n=1}^{N} c_n \cdot (\Phi - \Phi_w)^n$$

其中 $c_1 = \kappa$（线性灵敏度），$c_2 = \kappa'/2$（二阶项），等等。

#### 优化目标

$$\min_{\{c_n\}} \sum_{t_d} \left| \delta p_e^{\text{meas}}(t_d) - \int K(t-t_d) \cdot \sum_n c_n [\Phi(t) - \Phi_w]^n \, dt \right|^2$$

这是一个关于 $\{c_n\}$ 的线性最小二乘问题（因为 $\Phi(t)$ 已知，非线性在 $\Phi$ 上而非 $c_n$ 上）！

#### 具体实现

定义 $\Phi_n(t) = [\Phi(t) - \Phi_w]^n$，则：

$$\delta p_e(t_d) = \sum_n c_n \underbrace{\int K(t-t_d) \cdot \Phi_n(t) \, dt}_{R_n(t_d)}$$

$R_n(t_d)$ 是 $K$ 与 $\Phi_n$ 的卷积，可以预计算。于是：

$$\delta \mathbf{p}_e = \mathbf{R} \cdot \mathbf{c}$$

其中 $\mathbf{R} \in \mathbb{R}^{M \times N}$ 是响应矩阵，$\mathbf{c} = [c_1, c_2, \ldots, c_N]^T$。

求解（带 Tikhonov 正则化）：

$$\hat{\mathbf{c}} = (\mathbf{R}^T \mathbf{R} + \lambda \mathbf{I})^{-1} \mathbf{R}^T \delta\mathbf{p}_e$$

#### 优势

1. **一次测量提取多阶标定系数**：从 $c_1$（线性灵敏度）到 $c_N$（高阶非线性）
2. **计算简单**：线性最小二乘，无需迭代
3. **信号设计灵活**：理论上任何已知信号都可以，但信号应覆盖足够的幅度范围
4. **自然包含核函数效应**：不需要假设核函数为 delta 函数

#### 最优信号设计

为了最好地标定非线性系数，信号的选择很重要：

| 信号类型 | 特点 | 适合标定的阶数 |
|----------|------|--------------|
| 恒定 $h$ × 多个幅度 | 每次只采样一个 $\Phi$ 值 | 任意，但效率低 |
| 线性斜坡 | 连续扫描 $\Phi$ | 1-3 阶（分辨率受限） |
| 正弦 $A\sin(\omega t)$ | 重复扫描 $[-A, A]$ | 奇数阶 |
| 多频正弦 | 丰富的幅度分布 | 高阶 |
| 高斯脉冲（不同幅度） | 覆盖完整幅度范围 | 任意 |

**推荐**：使用 **多幅度高斯脉冲序列** 或 **慢斜坡信号**：
- 高斯脉冲提供平滑的幅度覆盖
- 慢斜坡简单直接，与方案二互补

---

## 3. 与 Cryoscope 标定的对比

Cryoscope 的标定步（见 `cryoscope_implementation.md` 第 2.3 节）也是施加已知方波偏置 $h$，扫描 $h$ 建立 $\varphi(h)$ 曲线。

| 维度 | Cryoscope 标定 | 瞬态磁场斜坡标定 (方案二) | Hammerstein 标定 (方案三) |
|------|---------------|-------------------------|------------------------|
| 信号 | 恒定方波，扫描幅度 | 线性斜坡，单次扫描 | 任意已知信号 |
| 测量次数 | $M$（每个 $h$ 一次） | 1 次滑动扫描 | 1 次滑动扫描 |
| 输出 | $\varphi(h)$ 查找表 | $\delta\omega(\Phi)$ 连续曲线 | $\{c_n\}$ 多项式系数 |
| 处理方法 | 插值反查 | Wiener 反卷积 | 线性最小二乘 |
| 非线性处理 | 自然包含（直接测量） | 通过反卷积+非线性反演 | 显式参数化 |
| 精度 | 受 $h$ 采样密度限制 | 受核函数分辨率限制 | 受模型阶数限制 |
| 计算量 | 低（插值） | 中（FFT） | 低（矩阵求逆） |

---

## 4. 方案三的完整代码实现

### 4.1 `analysis.py` 新增

```python
def calibrate_frequency_response(self, qubit, control_pulse, Phi_known, 
                                   delta_p, scan_list, n_order=4, lambdas=1e-3):
    """
    利用瞬态磁场协议标定 qubit 的频率-磁通非线性响应
    
    模型: δω(Φ) = Σ c_n · (Φ - Φ_w)^n,  n = 1, ..., n_order
    测量: δp_e(t_d) = ∫ K(t-t_d) · δω(Φ(t)) dt = Σ c_n · R_n(t_d)
    
    参数:
        qubit:          TransmonQubit 对象
        control_pulse:  Ramsey 脉冲序列（需已计算 kernel）
        Phi_known:      已知的施加信号 (Signal 对象)
        delta_p:        测量的概率变化 (array)
        scan_list:      滑动延迟时间列表 (array)
        n_order:        多项式展开阶数
        lambdas:        正则化参数
    
    返回:
        coeffs:         标定系数 [c_1, c_2, ..., c_n_order]
        freq_response:  函数 δω(Φ)
        R_matrix:       响应矩阵（可用于诊断）
    """
    # 获取核函数
    if not hasattr(control_pulse, 'kernel'):
        control_pulse.get_kernel(qubit)
    kernel = np.array(control_pulse.kernel)
    t_kernel = control_pulse.t_samples
    
    # 获取已知信号在各时刻的值
    Phi_w = qubit.flux
    t_signal = Phi_known.t_list
    phi_values = Phi_known.signal  # Φ(t) - 已知
    
    M = len(scan_list)       # 测量点数
    N = n_order              # 未知系数数
    
    # 构建响应矩阵 R[m, n] = ∫ K(t - t_d[m]) · (Φ(t) - Φ_w)^(n+1) dt
    R = np.zeros((M, N))
    
    dt_signal = t_signal[1] - t_signal[0]
    
    for n in range(N):
        # Φ_n(t) = (Φ(t) - Φ_w)^(n+1)
        Phi_n = phi_values ** (n + 1)
        
        # R_n(t_d) = K * Φ_n (卷积)
        conv = np.convolve(kernel, Phi_n, mode='full') * dt_signal
        # 对齐到 scan_list
        R[:, n] = conv[:M] if len(conv) >= M else np.pad(conv, (0, M - len(conv)))
    
    # 求解正则化最小二乘
    delta_p_vec = np.array(delta_p)
    I = np.eye(N)
    coeffs = np.linalg.solve(R.T @ R + lambdas * I, R.T @ delta_p_vec)
    
    # 构建频率响应函数
    def freq_response(Phi):
        """δω(Φ) = Σ c_n · Φ^n"""
        result = np.zeros_like(np.asarray(Phi, dtype=float))
        for n in range(N):
            result += coeffs[n] * np.asarray(Phi) ** (n + 1)
        return result
    
    return coeffs, freq_response, R
```

### 4.2 `protocal.py` 新增 case

```python
case 8:  # 瞬态磁场协议用于频率标定
    # 参数
    Phi_max = self.params.get('Phi_max', 0.05)      # 斜坡最大磁通 (Φ0)
    Phi_min = self.params.get('Phi_min', -0.05)      # 斜坡最小磁通 (Φ0)
    T_ramp = self.params.get('T_ramp', 200.0)        # 斜坡持续时间 (ns)
    n_points = self.params.get('n_points', 400)       # 信号采样点数
    t_rabi = self.params.get('t_rabi', np.linspace(0, 10, 20))
    signal_type = self.params.get('signal_type', 'ramp')  # 'ramp' 或 'gaussian_set'
    
    omega_d = qubit.frequency
    t_list = np.linspace(0, T_ramp, n_points)
    
    # 创建已知标定信号
    if signal_type == 'ramp':
        Phi_cal = Signal(type=0, t_list=t_list)
        Phi_cal.signal = np.linspace(Phi_min, Phi_max, n_points)
    elif signal_type == 'sine':
        amp = (Phi_max - Phi_min) / 2
        Phi_cal = Signal(type=2, t_list=t_list, amplitude=amp, 
                         frequency=1.0/T_ramp)
    
    Phi_cal.plot()
    Phi_0 = Signal(type=1, t_list=t_list, amplitude=0.0)
    
    # Ramsey 脉冲
    control_pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=omega_d)
    
    # 滑动测量
    print("Calibration: sliding measurement with known signal...")
    scan_list, p_e = self.sliding_measrement(qubit, Phi_cal, control_pulse)
    scan_list_base, p_e_base = self.sliding_measrement(qubit, Phi_0, control_pulse)
    
    # 核函数
    control_pulse.get_kernel(qubit)
    t_samples, kernel = control_pulse.t_samples, control_pulse.kernel
    delta_p = np.array(p_e) - np.array(p_e_base)
    
    return t_samples, kernel, scan_list, delta_p, p_e, Phi_cal, control_pulse
```

### 4.3 测试 Notebook

```python
import numpy as np
import matplotlib.pyplot as plt
from src.qubit import TransmonQubit
from src.signal import Signal
from src.protocal import Protocal
from src.analysis import Analysis

# ============ 初始化 ============
qubit = TransmonQubit(
    EC=0.2 * 2*np.pi, EJ=10.0 * 2*np.pi,
    T1=100e3, T2=50e3,
    flux=np.arctan(np.sqrt(2))/np.pi,  # 最优灵敏度点
    state=0, n_levels=2
)
kappa_true = qubit.frequency_sensitivity(qubit.flux)
print(f"True κ = {kappa_true:.4f} rad·GHz/Φ0")

# ============ 方案二：斜坡标定 ============
Phi_max = 0.03
T_ramp = 200.0
n_pts = 400

# 使用 case 8 (或直接调用 case 4 的框架)
cal = Protocal(type=8, Phi_max=Phi_max, Phi_min=-Phi_max, T_ramp=T_ramp,
               n_points=n_pts, signal_type='ramp',
               t_rabi=np.linspace(0, 10, 20))
cal.initialize(qubit)
t_samples, kernel, scan_list, delta_p, p_e, Phi_cal, control_pulse = cal.evolve(qubit)

# Wiener 反卷积：从 delta_p 还原 δω(t)
analysis = Analysis()
dt = scan_list[1] - scan_list[0]
t_rec, omega_rec = analysis.wiener_deconvolution(delta_p, kernel, dt, lambdas=0.01)

# 将时间轴映射到磁通轴
alpha = 2 * Phi_max / T_ramp
Phi_axis = -Phi_max + alpha * t_rec

# 理论曲线
import math
def delta_omega_theory(Phi):
    """精确的 δω(Φ) = f_Q(Φ_w + Φ) - f_Q(Φ_w)"""
    Phi_w = qubit.flux
    EJ_shifted = qubit.EJ_0 * abs(math.cos(math.pi * (Phi_w + Phi)))
    f_shifted = np.sqrt(8 * EJ_shifted * qubit.EC) - qubit.EC
    return f_shifted - qubit.frequency

Phi_theory = np.linspace(-Phi_max, Phi_max, 200)
omega_theory = np.array([delta_omega_theory(p) for p in Phi_theory])

plt.figure(figsize=(10, 5))
plt.plot(Phi_theory, omega_theory / (2*np.pi), 'k-', linewidth=2, label='Theory $\\delta\\omega(\\Phi)$')
plt.plot(Phi_axis, omega_rec / (2*np.pi), 'r--', label='Ramp + Wiener')
plt.plot(Phi_theory, kappa_true * Phi_theory / (2*np.pi), 'b:', label='Linear approx $\\kappa\\Phi$')
plt.xlabel('Flux offset $\\delta\\Phi$ ($\\Phi_0$)')
plt.ylabel('Frequency shift $\\delta f$ (GHz)')
plt.legend()
plt.grid()
plt.title('Frequency calibration via transient field protocol')
plt.show()

# ============ 方案三：Hammerstein 多项式标定 ============
coeffs, freq_func, R = analysis.calibrate_frequency_response(
    qubit, control_pulse, Phi_cal, delta_p, scan_list,
    n_order=4, lambdas=1e-3
)

print("Calibrated coefficients:")
print(f"  c1 (κ, linear):     {coeffs[0]:.4f}  (true: {kappa_true:.4f})")
if len(coeffs) > 1:
    print(f"  c2 (quadratic):     {coeffs[1]:.6f}")
if len(coeffs) > 2:
    print(f"  c3 (cubic):         {coeffs[2]:.6f}")

omega_hw = freq_func(Phi_theory)
plt.figure(figsize=(10, 5))
plt.plot(Phi_theory, omega_theory / (2*np.pi), 'k-', linewidth=2, label='Theory')
plt.plot(Phi_theory, omega_hw / (2*np.pi), 'g--', linewidth=2, label=f'Polynomial fit (order {len(coeffs)})')
plt.plot(Phi_theory, kappa_true * Phi_theory / (2*np.pi), 'b:', label='Linear $\\kappa\\Phi$')
plt.xlabel('Flux offset $\\delta\\Phi$ ($\\Phi_0$)')
plt.ylabel('Frequency shift $\\delta f$ (GHz)')
plt.legend()
plt.grid()
plt.title('Hammerstein polynomial calibration')
plt.show()
```

---

## 5. 总结：为什么用瞬态磁场协议做频率标定

### 5.1 标准 Ramsey 标定的局限

标准方法需要**逐点扫描**：在每个磁通点 $\Phi_i$ 处独立做一组完整的 Ramsey 实验（扫描 $\tau$、拟合振荡频率）。如果要标定 $M$ 个磁通点，需要 $M \times N_\tau$ 次测量。

### 5.2 瞬态协议的核心优势

1. **单次扫描获取连续曲线**（方案二）：用斜坡信号把时间轴映射到磁通轴，一次滑动扫描同时采样所有磁通值。测量次数从 $M \times N_\tau$ 降为 $N_{\text{delay}}$。

2. **显式提取非线性项**（方案三）：将 $\delta\omega(\Phi)$ 参数化为多项式，对已知信号做卷积后，问题退化为线性最小二乘。可以一步提取 $\kappa$、$\kappa'$、$\kappa''$ 等所有阶的标定系数。

3. **与波形重建/预失真共享基础设施**：标定和重建使用完全相同的实验框架（滑动测量 + 核函数），只需切换已知/未知量。这意味着频率标定可以作为波形重建的 **前置步骤**，两者共享硬件和数据处理流程。

### 5.3 适用场景

| 场景 | 推荐方案 |
|------|---------|
| 快速确认 $f_{01}$ | 标准 Ramsey（最简单） |
| 精确标定 $\kappa$ | 方案二（斜坡）或方案三 |
| 标定非线性项 $\kappa', \kappa''$ | 方案三（多项式拟合） |
| 标定 + 波形重建一体化 | 方案三 + case 4（共享核函数） |
| 动态频率跟踪（频率漂移） | 方案二（重复斜坡扫描） |
