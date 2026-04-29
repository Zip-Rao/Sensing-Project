# Cryoscope 协议实现方案

## 1. 协议概述

基于 Rol et al., Appl. Phys. Lett. 116, 054001 (2020) 的 Cryoscope 方法，结合本仿真平台的 Transmon qubit 模型，设计两步测量协议，用于**在片上重建任意波形**或**表征控制线路失真**。

### 1.1 与原始文献的区别

原始 Cryoscope 使用差分 Ramsey 相位直接反演瞬时频率，再通过已知的 $f_Q(\Phi_Q)$ 关系还原磁通。本方案将其拆分为两步：

1. **标定步（Calibration）**：建立方波偏置高度 $h$ 与积累相位 $\varphi$ 的映射关系 $\varphi(h)$
2. **测量步（Measurement）**：施加待测波形并扫描截断延迟 $t_{\text{delay}}$，得到 $\varphi(t_{\text{delay}})$

结合两者即可还原 $h(t_{\text{delay}})$，即波形的时域波形。

### 1.2 应用场景

- **波形探测**：还原 AWG 输出的实际波形形状（如方波、高斯脉冲等）
- **失真表征**：测量控制线路的阶跃响应，量化 IIR/FIR 失真

---

## 2. 理论推导

### 2.1 Transmon 频率与磁通的关系

Transmon qubit 的频率依赖于 SQUID 环路中的磁通 $\Phi_Q$：

$$f_Q(\Phi_Q) \approx \frac{1}{\hbar}\left(\sqrt{8E_J E_C \left|\cos\left(\pi\frac{\Phi_Q}{\Phi_0}\right)\right|} - E_C\right)$$

定义失谐量：

$$\Delta f_Q(t) = f_{\max} - f_Q(\Phi_Q(t))$$

其中 $f_{\max}$ 为磁通甜点处的最大频率。

### 2.2 磁通灵敏度

在工作点 $\Phi_0$ 处定义磁通灵敏度：

$$\kappa(\Phi_0) = \frac{\partial f_Q}{\partial \Phi}\bigg|_{\Phi=\Phi_0}$$

对于小偏置 $h$（以 $\Phi_0$ 为单位），频率变化近似为：

$$\Delta f_Q \approx \kappa \cdot h + \frac{1}{2}\kappa' \cdot h^2 + \cdots$$

在甜点（$\Phi_0 = 0$）处，$\kappa = 0$，一阶灵敏度为零，频率变化为二次的：

$$\Delta f_Q \approx \alpha \cdot h^2, \quad \alpha = \frac{4\pi^2 E_J E_C}{f_{\max} + E_C}$$

在非甜点工作点处（如最优灵敏度点 $\Phi = \arctan(\sqrt{2})/\pi$），$\kappa \neq 0$，频率变化以线性项为主。

### 2.3 第一步：标定 — 建立 $\varphi(h)$ 关系

#### 脉冲序列

```
(π/2)_Y → z-bias 方波（高度 h, 持续时间 T_bias） → (π/2)_{X 或 Y} → 读出
```

详细时序：

1. $(\pi/2)_Y$ 脉冲：将 qubit 从 $|0\rangle$ 旋转到 $\frac{1}{\sqrt{2}}(|0\rangle + |1\rangle)$
2. **z 方向偏置方波**：在自由演化期间施加恒定磁通偏置 $h$，持续时间 $T_{\text{bias}}$，qubit 积累相位
3. $(\pi/2)_{X}$ 和 $(\pi/2)_{Y}$ 脉冲：分别提取 $\langle X \rangle$ 和 $\langle Y \rangle$ 分量

#### 相位积累

在自由演化期间，qubit 积累的相位为：

$$\varphi(h) = 2\pi \int_0^{T_{\text{bias}}} \Delta f_Q(h) \, dt = 2\pi \cdot \Delta f_Q(h) \cdot T_{\text{bias}}$$

由于 $h$ 为常数，因此：

$$\varphi(h) = 2\pi \cdot \left[f_{\max} - f_Q(h)\right] \cdot T_{\text{bias}}$$

#### IQ 读出提取相位

通过两次测量（最后 $\pi/2$ 脉冲分别沿 $X$ 和 $Y$ 轴），得到：

$$\langle X \rangle = \cos(\varphi), \quad \langle Y \rangle = \sin(\varphi)$$

即：

$$p_e^{(X)}(h) = \frac{1}{2}\left(1 - \cos\varphi(h)\right), \quad p_e^{(Y)}(h) = \frac{1}{2}\left(1 + \sin\varphi(h)\right)$$

从而：

$$\varphi(h) = \text{atan2}\left(2p_e^{(Y)} - 1, \; 1 - 2p_e^{(X)}\right)$$

#### 扫描 $h$ 建立标定曲线

扫描 $h \in [h_{\min}, h_{\max}]$，对每个 $h$ 执行上述 Ramsey 实验，得到 $\varphi(h)$ 的标定曲线。

**关键**：标定步的 $T_{\text{bias}}$ 应与测量步中 AWG 的时间步长 $\Delta\tau$ 相匹配，使得后续测量的每个时间窗口积累的相位可以直接用标定曲线反演。

### 2.4 第二步：测量 — 扫描 $t_{\text{delay}}$ 还原波形

#### 脉冲序列

```
(π/2)_Y → 待测波形（截断到 t_delay + Δτ） → (π/2)_{X 或 Y} → 读出
```

详细时序：

1. $(\pi/2)_Y$ 脉冲：初始化叠加态
2. **待测波形**：施加截断到时刻 $t_{\text{delay}} + \Delta\tau$ 的待测磁通波形 $\Phi_{\text{test}}(t)$，其中波形在 $t > t_{\text{delay}} + \Delta\tau$ 后截断为零
3. $(\pi/2)_{X/Y}$ 脉冲 + 读出

#### 差分相位提取

对两个相邻截断时刻 $t_{\text{delay}} = \tau$ 和 $\tau + \Delta\tau$，分别测量积累相位 $\varphi_\tau$ 和 $\varphi_{\tau+\Delta\tau}$，则：

$$\Delta\varphi = \varphi_{\tau+\Delta\tau} - \varphi_\tau = 2\pi \int_\tau^{\tau+\Delta\tau} \Delta f_Q(\Phi_{\text{test}}(t)) \, dt$$

当 $\Delta\tau$ 足够小时：

$$\Delta\varphi \approx 2\pi \cdot \Delta f_Q(\Phi_{\text{test}}(\tau)) \cdot \Delta\tau$$

这给出了在时刻 $\tau$ 处的瞬时频率偏移 $\overline{\Delta f}_R$：

$$\overline{\Delta f}_R(\tau) = \frac{\Delta\varphi}{2\pi \cdot \Delta\tau} = \frac{\varphi_{\tau+\Delta\tau} - \varphi_\tau}{2\pi \cdot \Delta\tau}$$

#### 误差抑制（甜点处的二次非线性）

文献中指出，差分操作的关键优势在于：截断后的关断瞬态（turn-off transient）对相位的贡献几乎完全抵消。这是因为在甜点处，$\Delta f_Q(\Phi_Q) \propto \Phi_Q^2$，当磁通快速回到甜点时，二次非线性强烈抑制了关断瞬态对相位的贡献。误差量级：

$$|\varepsilon| / \overline{\Delta f}_R \lesssim 10^{-2} \sim 10^{-3}$$

### 2.5 波形还原

对于每个扫描时刻 $\tau$，得到瞬时频率偏移 $\overline{\Delta f}_R(\tau)$。将其与标定步建立的 $\varphi(h)$ 关系对应：

$$\Delta\varphi(\tau) = \varphi(h_{\text{eff}}(\tau))$$

其中 $h_{\text{eff}}(\tau)$ 为该时刻的等效偏置高度。通过反查标定曲线：

$$h_{\text{eff}}(\tau) = \varphi^{-1}\left(\Delta\varphi(\tau)\right)$$

这里 $\varphi^{-1}$ 为标定曲线的反函数，可通过插值实现。

**对于线性响应区（非甜点工作点）**：

$$h(\tau) \approx \frac{\overline{\Delta f}_R(\tau)}{\kappa}$$

**对于二次响应区（甜点工作点）**：

需要通过完整的 $f_Q(\Phi_Q)$ 关系进行数值反演：

$$h(\tau) = f_Q^{-1}\left(f_{\max} - \overline{\Delta f}_R(\tau)\right)$$

### 2.6 信号处理细节

#### Savitzky-Golay 微分滤波器

为减少数值微分放大噪声，使用二阶 Savitzky-Golay 滤波器对 $\varphi(\tau)$ 进行局部多项式拟合后求导：

$$\overline{\Delta f}_R(\tau) = \frac{1}{2\pi} \cdot \text{SG-derivative}\left[\varphi(\tau)\right]$$

#### 相位解缠绕

当波形导致大失谐（$\Delta f_Q \gg 1/\Delta\tau$）时，相邻时间步之间的相位差可能超过 $\pi$，导致混叠。此时需要：

1. 先用 FFT 确定主频 $f_{\text{demod}}$，去调制
2. 对残余相位进行 `unwrap`
3. 最终频率 = SG 导数 + $f_{\text{demod}}$ + $n \cdot f_{\text{Nyquist}}$（如需要）

在本仿真中，由于信号幅度可控，一般不需要处理 Nyquist 跳跃。

---

## 3. 实现方案

### 3.1 与现有代码模块的对应关系

| 功能 | 模块 | 现有/新增 |
|------|------|-----------|
| Ramsey 脉冲序列 + z-bias | `pulse.py` → `create_cryoscope_calibration_pulse()` | **新增** |
| Ramsey 脉冲序列 + 截断波形 | `pulse.py` → `create_cryoscope_measurement_pulse()` | **新增** |
| 标定协议执行 | `protocal.py` → `case 6` | **新增** |
| 测量协议执行 | `protocal.py` → `case 7` | **新增** |
| 相位提取（IQ） | `analysis.py` → `get_signal_from_ramsey_by_iq()` | **可复用** |
| 标定曲线拟合/插值 | `analysis.py` → `build_calibration_curve()` | **新增** |
| 波形还原 | `analysis.py` → `reconstruct_waveform_cryoscope()` | **新增** |
| z-bias 方波信号 | `signal.py` → `Signal(type=1)` | **可复用** |

### 3.2 Step 1: 标定脉冲构建

新增 `create_cryoscope_calibration_pulse()` 函数：

```python
def create_cryoscope_calibration_pulse(t_rabi, T_bias, omega_d, phase_analysis='X'):
    """
    创建 Cryoscope 标定用 Ramsey 脉冲序列
    
    序列结构：(π/2)_Y → free(T_bias) → (π/2)_{X or Y}
    
    参数:
        t_rabi:   π/2 脉冲时间数组 (ns)
        T_bias:   z-bias 方波持续时间 (ns)，即自由演化时间
        omega_d:  驱动频率 (GHz)
        phase_analysis: 分析脉冲相位, 'X' 对应 phase=0, 'Y' 对应 phase=π/2
    
    返回:
        CompositePulse 对象
    
    说明:
        z-bias 方波不在脉冲序列中，而是作为 Signal 在 protocal 中与 qubit_in_mag 结合
        此处脉冲序列只包含微波驱动部分
    """
    phase2 = 0.0 if phase_analysis == 'X' else np.pi / 2
    return create_ramsey_pulse(t_rabi, tau=T_bias, omega_d=omega_d, 
                                phase1=np.pi/2, phase2=phase2)
```

### 3.3 Step 2: 标定协议（`protocal.py` case 6）

```python
case 6:  # Cryoscope 标定
    # 参数
    T_bias = self.params.get('T_bias', 20.0)          # z-bias 持续时间 (ns)
    h_list = self.params.get('h_list', np.linspace(-0.1, 0.1, 101))  # 偏置高度扫描
    t_rabi = self.params.get('t_rabi', np.linspace(0, 10, 20))
    omega_d = qubit.frequency
    
    p_e_X_list = []
    p_e_Y_list = []
    psi_e = basis(qubit.n_levels, 1)
    
    for h in h_list:
        for phase_label in ['X', 'Y']:
            # 1. 创建 z-bias 方波信号
            Phi_bias = Signal(type=1, t_list=np.linspace(0, T_bias, max(int(T_bias*2), 40)),
                              amplitude=h)
            
            # 2. 创建 Ramsey 脉冲
            phase2 = 0.0 if phase_label == 'X' else np.pi / 2
            control_pulse = create_ramsey_pulse(t_rabi, tau=T_bias, omega_d=omega_d,
                                                phase1=np.pi/2, phase2=phase2)
            control_pulse.t_list = control_pulse.t_list - t_rabi[-1]
            
            # 3. 将 qubit 置于 z-bias 磁场下
            qubit.qubit_in_mag(Phi_bias, frame=1, omega_d=omega_d)
            
            # 4. 构建总哈密顿量并演化
            t_global = np.linspace(-t_rabi[-1]-5, T_bias + t_rabi[-1]+5, 
                                    max(500, int((T_bias + 2*t_rabi[-1]+10)*5)))
            H = QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1) \
              + QobjEvo(control_pulse.hamiltonian, tlist=control_pulse.t_list, order=1)
            result = mesolve(H, qubit.state, t_global, [], e_ops=[psi_e * psi_e.dag()])
            p_e = result.expect[0][-1]
            
            if phase_label == 'X':
                p_e_X_list.append(p_e)
            else:
                p_e_Y_list.append(p_e)
    
    # 5. 提取相位
    p_e_X = np.array(p_e_X_list)
    p_e_Y = np.array(p_e_Y_list)
    cos_phi = 1 - 2 * p_e_X      # <X> = cos(φ)
    sin_phi = 2 * p_e_Y - 1      # <Y> = sin(φ)
    phi_cal = np.arctan2(sin_phi, cos_phi)
    phi_cal = np.unwrap(phi_cal)
    
    return h_list, phi_cal, p_e_X, p_e_Y
```

### 3.4 Step 3: 测量协议（`protocal.py` case 7）

```python
case 7:  # Cryoscope 波形测量
    # 参数
    Phi_test = self.params.get('Phi_test')   # 待测波形 Signal 对象
    t_rabi = self.params.get('t_rabi', np.linspace(0, 10, 20))
    T_sep = self.params.get('T_sep', None)   # 两个 π/2 脉冲的间隔
    omega_d = qubit.frequency
    
    # 扫描截断延迟
    tau_list = Phi_test.t_list.copy()
    dt = tau_list[1] - tau_list[0]  # 时间步长 = AWG 分辨率
    
    if T_sep is None:
        T_sep = tau_list[-1] + 100  # 比最大截断时间多 100 ns
    
    p_e_X_list = []
    p_e_Y_list = []
    psi_e = basis(qubit.n_levels, 1)
    
    for tau in tau_list:
        for phase_label in ['X', 'Y']:
            # 1. 截断波形到 tau
            t_trunc = Phi_test.t_list[Phi_test.t_list <= tau]
            if len(t_trunc) < 2:
                t_trunc = Phi_test.t_list[:2]
            sig_trunc_values = Phi_test.signal[:len(t_trunc)]
            
            # 创建截断后的复合信号：截断波形 + 零信号（到 T_sep）
            Phi_trunc = Signal(type=1, t_list=t_trunc, amplitude=0.0)
            Phi_trunc.signal = sig_trunc_values  # 直接覆盖信号值
            
            T_remaining = T_sep - tau
            Phi_zero = Signal(type=0, t_list=np.linspace(0, max(T_remaining, 1), 
                              max(int(T_remaining*2), 20)))
            Phi_full = CompositeSignal([Phi_trunc, Phi_zero])
            
            # 2. 创建 Ramsey 脉冲
            phase2 = 0.0 if phase_label == 'X' else np.pi / 2
            control_pulse = create_ramsey_pulse(t_rabi, tau=T_sep, omega_d=omega_d,
                                                phase1=np.pi/2, phase2=phase2)
            control_pulse.t_list = control_pulse.t_list - t_rabi[-1]
            
            # 3. qubit 在截断波形下演化
            qubit.qubit_in_mag(Phi_full, frame=1, omega_d=omega_d)
            
            # 4. 演化并测量
            t_global = np.linspace(-t_rabi[-1]-5, T_sep + t_rabi[-1]+5, 
                                    max(1000, int((T_sep + 2*t_rabi[-1]+10)*5)))
            H = QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1) \
              + QobjEvo(control_pulse.hamiltonian, tlist=control_pulse.t_list, order=1)
            result = mesolve(H, qubit.state, t_global, [], e_ops=[psi_e * psi_e.dag()])
            p_e = result.expect[0][-1]
            
            if phase_label == 'X':
                p_e_X_list.append(p_e)
            else:
                p_e_Y_list.append(p_e)
    
    # 5. 提取相位
    p_e_X = np.array(p_e_X_list)
    p_e_Y = np.array(p_e_Y_list)
    cos_phi = 1 - 2 * p_e_X
    sin_phi = 2 * p_e_Y - 1
    phi_meas = np.arctan2(sin_phi, cos_phi)
    phi_meas = np.unwrap(phi_meas)
    
    return Phi_test, tau_list, phi_meas, p_e_X, p_e_Y
```

### 3.5 Step 4: 分析与波形还原（`analysis.py`）

新增两个函数：

```python
def build_calibration_curve(self, h_list, phi_cal):
    """
    根据标定数据构建 φ(h) 和 h(φ) 的插值函数
    
    参数:
        h_list:  偏置高度数组 (Φ0)
        phi_cal: 对应的积累相位数组 (rad)
    
    返回:
        phi_of_h: 插值函数 φ(h)
        h_of_phi: 反函数 h(φ)
    """
    from scipy.interpolate import interp1d
    
    phi_of_h = interp1d(h_list, phi_cal, kind='cubic', fill_value='extrapolate')
    
    # 构建反函数：需要 φ(h) 单调
    # 如果在甜点工作，φ(h) ∝ h² 不单调 → 只取正半支
    # 如果在非甜点工作，φ(h) 近似线性 → 直接反查
    sorted_idx = np.argsort(phi_cal)
    phi_sorted = phi_cal[sorted_idx]
    h_sorted = h_list[sorted_idx]
    
    # 去除重复值
    unique_mask = np.diff(phi_sorted, prepend=-np.inf) > 1e-12
    h_of_phi = interp1d(phi_sorted[unique_mask], h_sorted[unique_mask], 
                         kind='cubic', fill_value='extrapolate')
    
    return phi_of_h, h_of_phi


def reconstruct_waveform_cryoscope(self, qubit, tau_list, phi_meas, 
                                     h_of_phi=None, T_bias=None,
                                     sg_window=5, sg_poly=2):
    """
    从 Cryoscope 测量数据还原波形
    
    参数:
        qubit:     TransmonQubit 对象
        tau_list:  截断延迟时间数组 (ns)
        phi_meas:  各截断时刻的积累相位 (rad), 已 unwrap
        h_of_phi:  标定反函数 h(φ) (可选, 若不提供则使用解析公式)
        T_bias:    标定步的 bias 持续时间 (ns), 用于计算 Δφ 对应的 h
        sg_window: Savitzky-Golay 窗口大小 (奇数)
        sg_poly:   Savitzky-Golay 多项式阶数
    
    返回:
        tau_list:  时间数组 (ns)
        h_recon:   重建的偏置高度数组 (Φ0)
    """
    from scipy.signal import savgol_filter
    
    dt = tau_list[1] - tau_list[0]
    
    # 方法一：使用 SG 滤波器做微分，得到瞬时频率偏移
    if sg_window >= len(tau_list):
        sg_window = max(3, len(tau_list) // 2 * 2 - 1)
    
    dphi_dt = savgol_filter(phi_meas, sg_window, sg_poly, deriv=1, delta=dt)
    delta_f = dphi_dt / (2 * np.pi)  # 瞬时频率偏移 (GHz)
    
    if h_of_phi is not None and T_bias is not None:
        # 方法二：使用标定曲线直接反查
        # 每个 Δτ 窗口的相位差
        delta_phi = np.diff(phi_meas, prepend=0)
        # delta_phi 对应于 T_bias = dt 时间内积累的相位
        # 如果标定步的 T_bias != dt，需要缩放
        delta_phi_scaled = delta_phi * T_bias / dt
        h_recon = np.array([h_of_phi(dp) for dp in delta_phi_scaled])
    else:
        # 方法三：使用解析公式反演（适用于非甜点工作点）
        kappa = qubit.frequency_sensitivity(qubit.flux)
        if abs(kappa) > 1e-6:
            # 线性区域：Δf ≈ κ·h
            h_recon = delta_f / kappa
        else:
            # 甜点处：需要完整的 f_Q(Φ) 反函数
            # f_Q = sqrt(8*EJ*EC*|cos(π*Φ)|) - EC
            # Δf = f_max - f_Q → f_Q = f_max - Δf
            f_Q = qubit.frequency - delta_f  # f_max - Δf
            ratio = (f_Q + qubit.EC)**2 / (8 * qubit.EC * qubit.EJ_0)
            ratio = np.clip(ratio, -1, 1)
            h_recon = np.arccos(ratio) / np.pi
    
    return tau_list, h_recon
```

---

## 4. 完整代码实现

### 4.1 `pulse.py` 新增函数

```python
def create_cryoscope_ramsey_pulse(t_rabi, T_free, omega_d, phase_analysis='X'):
    """
    Cryoscope 用 Ramsey 脉冲：(π/2)_Y → free(T_free) → (π/2)_{X or Y}
    
    第一个 π/2 脉冲沿 Y 轴 (phase = π/2)，
    分析脉冲根据 phase_analysis 选择 X 或 Y。
    """
    phase2 = 0.0 if phase_analysis == 'X' else np.pi / 2
    return create_ramsey_pulse(t_rabi, tau=T_free, omega_d=omega_d,
                                phase1=np.pi/2, phase2=phase2)
```

### 4.2 `protocal.py` 新增 case

```python
case 6:  # Cryoscope 标定：扫描 h 得到 φ(h)
    T_bias = self.params.get('T_bias', 20.0)
    h_list = self.params.get('h_list', np.linspace(-0.05, 0.05, 51))
    t_rabi = self.params.get('t_rabi', np.linspace(0, 10, 20))
    omega_d = qubit.frequency
    
    phi_cal_list = []
    p_e_X_list = []
    p_e_Y_list = []
    psi_e = basis(qubit.n_levels, 1)
    
    for h in h_list:
        p_e_pair = {}
        for phase_label in ['X', 'Y']:
            # z-bias 方波信号
            Phi_bias = Signal(type=1, 
                              t_list=np.linspace(0, T_bias, max(int(T_bias*2), 40)),
                              amplitude=h)
            
            # Ramsey 脉冲
            phase2 = 0.0 if phase_label == 'X' else np.pi / 2
            control_pulse = create_ramsey_pulse(
                t_rabi, tau=T_bias, omega_d=omega_d,
                phase1=np.pi/2, phase2=phase2
            )
            control_pulse.t_list = control_pulse.t_list - t_rabi[-1]
            
            # qubit 在 z-bias 下
            qubit.qubit_in_mag(Phi_bias, frame=1, omega_d=omega_d)
            
            # 构建总 H 并演化
            t_start = control_pulse.t_list[0] - 5
            t_end = control_pulse.t_list[-1] + 5
            n_pts = max(500, int((t_end - t_start) * 5))
            t_global = np.linspace(t_start, t_end, n_pts)
            
            H = (QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1) 
               + QobjEvo(control_pulse.hamiltonian, tlist=control_pulse.t_list, order=1))
            result = mesolve(H, qubit.state, t_global, [], 
                             e_ops=[psi_e * psi_e.dag()])
            p_e_pair[phase_label] = result.expect[0][-1]
        
        p_e_X_list.append(p_e_pair['X'])
        p_e_Y_list.append(p_e_pair['Y'])
    
    p_e_X = np.array(p_e_X_list)
    p_e_Y = np.array(p_e_Y_list)
    cos_phi = 1 - 2 * p_e_X
    sin_phi = 2 * p_e_Y - 1
    phi_cal = np.unwrap(np.arctan2(sin_phi, cos_phi))
    
    return h_list, phi_cal, p_e_X, p_e_Y

case 7:  # Cryoscope 波形测量：扫描 t_delay 得到 φ(t_delay)
    Phi_test = self.params.get('Phi_test')
    t_rabi = self.params.get('t_rabi', np.linspace(0, 10, 20))
    T_sep_extra = self.params.get('T_sep_extra', 100)
    omega_d = qubit.frequency
    
    tau_list = Phi_test.t_list.copy()
    T_sep = tau_list[-1] + T_sep_extra
    
    p_e_X_list = []
    p_e_Y_list = []
    psi_e = basis(qubit.n_levels, 1)
    
    for i, tau in enumerate(tau_list):
        p_e_pair = {}
        for phase_label in ['X', 'Y']:
            # 截断待测波形到 tau
            mask = Phi_test.t_list <= tau
            if np.sum(mask) < 2:
                mask[:2] = True
            t_trunc = Phi_test.t_list[mask]
            sig_trunc = Phi_test.signal[mask]
            
            # 截断波形 + 后续零信号 → CompositeSignal
            Phi_trunc = Signal(type=0, t_list=t_trunc)
            Phi_trunc.signal = sig_trunc.copy()
            
            T_remain = T_sep - t_trunc[-1]
            Phi_zero = Signal(type=0, t_list=np.linspace(0, max(T_remain, 1),
                              max(int(T_remain * 2), 20)))
            Phi_full = CompositeSignal([Phi_trunc, Phi_zero])
            
            # Ramsey 脉冲
            phase2 = 0.0 if phase_label == 'X' else np.pi / 2
            control_pulse = create_ramsey_pulse(
                t_rabi, tau=T_sep, omega_d=omega_d,
                phase1=np.pi/2, phase2=phase2
            )
            control_pulse.t_list = control_pulse.t_list - t_rabi[-1]
            
            # qubit 在截断波形下演化
            qubit.qubit_in_mag(Phi_full, frame=1, omega_d=omega_d)
            
            # 演化
            t_start = control_pulse.t_list[0] - 5
            t_end = control_pulse.t_list[-1] + 5
            n_pts = max(1000, int((t_end - t_start) * 5))
            t_global = np.linspace(t_start, t_end, n_pts)
            
            H = (QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1)
               + QobjEvo(control_pulse.hamiltonian, tlist=control_pulse.t_list, order=1))
            result = mesolve(H, qubit.state, t_global, [],
                             e_ops=[psi_e * psi_e.dag()])
            p_e_pair[phase_label] = result.expect[0][-1]
        
        p_e_X_list.append(p_e_pair['X'])
        p_e_Y_list.append(p_e_pair['Y'])
        
        if (i + 1) % 20 == 0:
            print(f"Cryoscope measurement: {i+1}/{len(tau_list)} points done")
    
    p_e_X = np.array(p_e_X_list)
    p_e_Y = np.array(p_e_Y_list)
    cos_phi = 1 - 2 * p_e_X
    sin_phi = 2 * p_e_Y - 1
    phi_meas = np.unwrap(np.arctan2(sin_phi, cos_phi))
    
    return Phi_test, tau_list, phi_meas, p_e_X, p_e_Y
```

### 4.3 `analysis.py` 新增函数

```python
def build_calibration_curve(self, h_list, phi_cal):
    """
    构建 Cryoscope 标定曲线的插值函数
    
    参数:
        h_list:  z-bias 高度数组 (Φ0)
        phi_cal: 对应积累相位数组 (rad)，已 unwrap
    返回:
        phi_of_h:  正向插值 φ(h)
        h_of_phi:  反向插值 h(φ)
    """
    from scipy.interpolate import interp1d
    
    phi_of_h = interp1d(h_list, phi_cal, kind='cubic', fill_value='extrapolate')
    
    sorted_idx = np.argsort(phi_cal)
    phi_sorted = phi_cal[sorted_idx]
    h_sorted = np.array(h_list)[sorted_idx]
    unique_mask = np.concatenate(([True], np.diff(phi_sorted) > 1e-12))
    h_of_phi = interp1d(phi_sorted[unique_mask], h_sorted[unique_mask],
                         kind='cubic', fill_value='extrapolate')
    
    return phi_of_h, h_of_phi

def reconstruct_waveform_cryoscope(self, qubit, tau_list, phi_meas,
                                     h_of_phi=None, T_bias=None,
                                     method='sg', sg_window=5, sg_poly=2):
    """
    从 Cryoscope 测量结果还原波形
    
    参数:
        qubit:      TransmonQubit 对象
        tau_list:   截断延迟时间数组 (ns)
        phi_meas:   积累相位数组 (rad)，已 unwrap
        h_of_phi:   标定反函数（可选）
        T_bias:     标定步的 bias 持续时间 (ns)
        method:     'sg' 用 Savitzky-Golay 微分; 'diff' 用差分; 'calibration' 用标定曲线
        sg_window:  SG 窗口大小
        sg_poly:    SG 多项式阶数
    返回:
        tau_list, h_recon
    """
    from scipy.signal import savgol_filter
    dt = tau_list[1] - tau_list[0]
    
    if method == 'calibration' and h_of_phi is not None and T_bias is not None:
        # 差分相位 → 通过标定曲线反查 h
        delta_phi = np.gradient(phi_meas, tau_list) * T_bias
        h_recon = h_of_phi(delta_phi)
    
    elif method == 'sg':
        if sg_window >= len(tau_list):
            sg_window = max(3, len(tau_list) // 2 * 2 - 1)
        dphi_dt = savgol_filter(phi_meas, sg_window, sg_poly, deriv=1, delta=dt)
        delta_f = dphi_dt / (2 * np.pi)
        
        kappa = qubit.frequency_sensitivity(qubit.flux)
        if abs(kappa) > 1e-6:
            h_recon = delta_f / kappa
        else:
            f_Q = qubit.frequency - delta_f
            ratio = (f_Q + qubit.EC)**2 / (8 * qubit.EC * qubit.EJ_0)
            ratio = np.clip(ratio, 0, 1)
            h_recon = np.arccos(ratio) / np.pi
    
    elif method == 'diff':
        delta_phi = np.gradient(phi_meas, tau_list)
        delta_f = delta_phi / (2 * np.pi)
        kappa = qubit.frequency_sensitivity(qubit.flux)
        h_recon = delta_f / kappa if abs(kappa) > 1e-6 else np.zeros_like(delta_f)
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return tau_list, h_recon
```

---

## 5. 测试方案

### 5.1 测试 Notebook 代码

```python
import numpy as np
import matplotlib.pyplot as plt
from src.qubit import TransmonQubit
from src.signal import Signal
from src.protocal import Protocal
from src.analysis import Analysis

# ============ 初始化 qubit ============
qubit = TransmonQubit(
    EC  = 0.2 * 2 * np.pi,
    EJ  = 10.0 * 2 * np.pi,
    T1  = 100.0e3,
    T2  = 50.0e3,
    flux = np.arctan(np.sqrt(2)) / np.pi,  # 最优灵敏度工作点
    state = 0,
    n_levels = 2
)

print(f"Qubit frequency: {qubit.frequency/(2*np.pi):.4f} GHz")
print(f"Flux sensitivity κ: {qubit.frequency_sensitivity(qubit.flux)/(2*np.pi):.4f} GHz/Φ0")

# ============ Step 1: 标定 ============
T_bias = 20.0  # ns
h_list = np.linspace(-0.02, 0.02, 41)

cal_protocol = Protocal(type=6, T_bias=T_bias, h_list=h_list,
                         t_rabi=np.linspace(0, 10, 20))
cal_protocol.initialize(qubit)
h_arr, phi_cal, p_X, p_Y = cal_protocol.evolve(qubit)

# 绘制标定曲线
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].plot(h_arr, p_X, 'o-', label='$p_e^{(X)}$')
axes[0].plot(h_arr, p_Y, 's-', label='$p_e^{(Y)}$')
axes[0].set_xlabel('Bias height h (Φ₀)')
axes[0].set_ylabel('Excited state probability')
axes[0].legend()
axes[0].grid()

axes[1].plot(h_arr, phi_cal, 'o-')
axes[1].set_xlabel('Bias height h (Φ₀)')
axes[1].set_ylabel('Phase φ (rad)')
axes[1].set_title('Calibration: φ(h)')
axes[1].grid()

# 与解析值比较
kappa = qubit.frequency_sensitivity(qubit.flux)
phi_theory = 2 * np.pi * kappa * h_arr * T_bias
axes[2].plot(h_arr, phi_cal, 'o', label='Simulated')
axes[2].plot(h_arr, phi_theory, '-', label='Linear theory: 2πκhT')
axes[2].set_xlabel('Bias height h (Φ₀)')
axes[2].set_ylabel('Phase φ (rad)')
axes[2].legend()
axes[2].grid()
plt.tight_layout()
plt.show()

# 构建标定插值
analysis = Analysis()
phi_of_h, h_of_phi = analysis.build_calibration_curve(h_arr, phi_cal)

# ============ Step 2: 波形测量 ============
# 待测波形：高斯脉冲
t_sig = np.linspace(0, 200, 400)
Phi_test = Signal(type=3, t_list=t_sig, amplitude=0.01, center=100, width=30)
Phi_test.plot()

meas_protocol = Protocal(type=7, Phi_test=Phi_test,
                          t_rabi=np.linspace(0, 10, 20), T_sep_extra=100)
meas_protocol.initialize(qubit)
_, tau_list, phi_meas, p_X_meas, p_Y_meas = meas_protocol.evolve(qubit)

# ============ Step 3: 波形还原 ============
# 方法1: SG 微分 + 线性反演
tau_rec, h_rec_sg = analysis.reconstruct_waveform_cryoscope(
    qubit, tau_list, phi_meas, method='sg', sg_window=7, sg_poly=2
)

# 方法2: 标定曲线
tau_rec2, h_rec_cal = analysis.reconstruct_waveform_cryoscope(
    qubit, tau_list, phi_meas, h_of_phi=h_of_phi, T_bias=T_bias, method='calibration'
)

# ============ 绘图比较 ============
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# 原始测量数据
axes[0, 0].plot(tau_list, p_X_meas, 'o-', markersize=2, label='$p_e^{(X)}$')
axes[0, 0].plot(tau_list, p_Y_meas, 's-', markersize=2, label='$p_e^{(Y)}$')
axes[0, 0].set_xlabel('Truncation delay τ (ns)')
axes[0, 0].set_ylabel('$p_e$')
axes[0, 0].legend()
axes[0, 0].grid()
axes[0, 0].set_title('Raw measurement data')

# 相位
axes[0, 1].plot(tau_list, phi_meas, '-')
axes[0, 1].set_xlabel('Truncation delay τ (ns)')
axes[0, 1].set_ylabel('Phase φ (rad)')
axes[0, 1].set_title('Unwrapped phase φ(τ)')
axes[0, 1].grid()

# SG 微分方法重建
axes[1, 0].plot(t_sig, Phi_test.signal, 'k-', linewidth=2, label='Original')
axes[1, 0].plot(tau_rec, h_rec_sg, 'r--', label='SG reconstruction')
axes[1, 0].set_xlabel('Time (ns)')
axes[1, 0].set_ylabel('Flux (Φ₀)')
axes[1, 0].legend()
axes[1, 0].grid()
axes[1, 0].set_title('Waveform reconstruction (SG)')

# 标定方法重建
axes[1, 1].plot(t_sig, Phi_test.signal, 'k-', linewidth=2, label='Original')
axes[1, 1].plot(tau_rec2, h_rec_cal, 'b--', label='Calibration reconstruction')
axes[1, 1].set_xlabel('Time (ns)')
axes[1, 1].set_ylabel('Flux (Φ₀)')
axes[1, 1].legend()
axes[1, 1].grid()
axes[1, 1].set_title('Waveform reconstruction (Calibration)')

plt.tight_layout()
plt.show()
```

### 5.2 失真表征测试

```python
# ============ 阶跃响应测量（Cryoscope 的原始应用）============
# 模拟带有指数衰减失真的方波
t_step = np.linspace(0, 300, 600)
step_ideal = np.heaviside(t_step - 20, 0.5) * 0.01  # 理想方波

# 添加模拟失真：指数衰减过冲
tau_distort = 40  # ns
A_distort = 0.15  # 过冲幅度比
step_distorted = step_ideal * (1 + A_distort * np.exp(-(t_step - 20) / tau_distort))
step_distorted[t_step < 20] = 0

Phi_step = Signal(type=0, t_list=t_step)
Phi_step.signal = step_distorted

# 执行 Cryoscope 测量
meas_step = Protocal(type=7, Phi_test=Phi_step,
                      t_rabi=np.linspace(0, 10, 20), T_sep_extra=100)
meas_step.initialize(qubit)
_, tau_s, phi_s, _, _ = meas_step.evolve(qubit)

# 还原
_, h_step_rec = analysis.reconstruct_waveform_cryoscope(
    qubit, tau_s, phi_s, method='sg', sg_window=9, sg_poly=2
)

# 归一化为阶跃响应
h_steady = np.mean(h_step_rec[-50:])
s_reconstructed = h_step_rec / h_steady if abs(h_steady) > 1e-10 else h_step_rec

plt.figure(figsize=(10, 5))
plt.plot(t_step, step_distorted / 0.01, 'k-', label='True step response', linewidth=2)
plt.plot(tau_s, s_reconstructed, 'r--', label='Cryoscope reconstruction')
plt.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5)
plt.fill_between(tau_s, 0.999, 1.001, alpha=0.2, color='green', label='±0.1% band')
plt.xlabel('Time (ns)')
plt.ylabel('Normalized step response')
plt.legend()
plt.grid()
plt.title('Cryoscope: Step Response Characterization')
plt.show()
```

---

## 6. 参数推荐

### 6.1 标定步参数

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| $T_{\text{bias}}$ | 10–50 ns | 与待测波形的时间分辨率匹配 |
| $h$ 扫描范围 | $[-0.05, 0.05]$ $\Phi_0$ | 覆盖待测波形的幅度范围 |
| $h$ 采样点数 | 51–101 | 确保标定曲线插值精度 |
| $t_{\text{rabi}}$ | `np.linspace(0, 10, 20)` | π/2 脉冲 ~10 ns |
| 工作点 | $\Phi = \arctan(\sqrt{2})/\pi$ | 最大灵敏度点 |

### 6.2 测量步参数

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| $\Delta\tau$ (时间步长) | 0.4–2 ns | 由 `tau_list` 步长决定，对应 AWG 采样率 |
| $T_{\text{sep}}$ | $\tau_{\max} + 100$ ns | 比最大截断长 100 ns |
| 待测波形幅度 | $< 0.05$ $\Phi_0$ | 确保不超过标定范围 |

### 6.3 数据处理参数

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| SG 窗口 | 5–11 | 越大越平滑但分辨率越低 |
| SG 阶数 | 2 | 二阶多项式 |

---

## 7. 注意事项

1. **IQ 测量 vs 单通道测量**：本方案使用 IQ 测量（分别沿 X 和 Y 轴读出），可直接提取相位而无需模糊的 `arccos` 解缠。代价是每个 $\tau$ 点需要两次测量。如果相位变化不大，可以只用单通道（最后 $\pi/2$ 沿 $X$ 轴），通过 `arccos` + 解缠绕提取相位。

2. **甜点 vs 非甜点工作**：
   - 甜点处：一阶磁通灵敏度为零，Cryoscope 利用二次依赖关系，可以抑制关断瞬态误差。但标定曲线 $\varphi(h)$ 是偶函数，无法区分正负偏置。
   - 非甜点处：线性灵敏度高，标定曲线单调，反演简单。但对噪声更敏感。

3. **退相干**：长 $T_{\text{sep}}$ 会导致退相干衰减。建议 $T_{\text{sep}} \ll T_2$（默认 $T_2 = 50$ μs，远大于 ns 量级的测量时间，在仿真中不是问题）。

4. **计算效率**：对于 $N$ 个截断点，每个点需要 2 次（IQ）mesolve 演化，总共 $2N$ 次。如果 $N = 400$，使用预编译的 `QobjEvo` 列表格式可显著加速。

5. **Nyquist 问题**：当波形导致的瞬时频率偏移 $\Delta f_Q > f_{\text{Nyquist}} = 1/(2\Delta\tau)$ 时，会出现频率混叠。需要选择足够小的 $\Delta\tau$ 或降低波形幅度。
