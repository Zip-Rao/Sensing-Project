# 超导量子比特量子传感仿真平台 - 噪声与真实条件添加方案

> ⚠️ **时效声明(2026-07-17 补注)**：本文件是**针对旧 `src/` 架构的噪声集成技术方案**,非勾选式任务清单。部分已在 `sqc/` 落地(`sqc/simulation/noise.py` 的 1/f 噪声、qubit 的 T1/T2),但全量 Lindblad 集成到各 workflow 尚未完成(见 [`../RELEASE_TODO.md`](../RELEASE_TODO.md) P2)。本文件保留作为**设计参考**。

## 1. 项目概述

本仿真平台旨在模拟基于超导Transmon量子比特的量子传感过程。当前已实现基本框架，包括：
- `TransmonQubit`类：Transmon量子比特模型
- `Signal`类：信号生成（恒定、正弦、高斯脉冲等）
- `Pulse`类：脉冲定义与哈密顿量计算
- `Protocol`类：量子传感协议（Rabi振荡、Ramsey干涉、自旋回波等）
- `Analysis`类：数据分析与可视化

当前仿真处于理想条件，为提升仿真真实性，需要添加以下真实条件：
1. **退相干噪声**：T1弛豫和T2退相干
2. **磁通噪声**：1/f噪声及静态偏移
3. **脉冲畸变**：有限带宽、幅度/相位噪声、时间抖动
4. **测量噪声**：读取误差
5. **环境噪声**：温度漂移、电荷噪声等

## 2. 当前架构分析

### 2.1 现有噪声相关功能
- `TransmonQubit`类已包含`T1`、`T2`属性，但未在演化中使用
- `Signal`类有`back_signal()`方法（返回零信号），可用于添加背景噪声
- `Pulse`类有`back_pulse()`方法（返回零信号），可用于添加脉冲噪声
- `TransmonQubit`类有`flux_noise()`和`sweet_point()`方法占位符

### 2.2 演化机制
- 所有量子演化通过QuTiP的`mesolve`函数进行
- 当前未包含Lindblad算符（退相干）
- 脉冲哈密顿量通过`Pulse.get_hamiltonian()`生成，使用理想信号

## 3. 噪声源分类

### 3.1 退相干噪声
- **能量弛豫（T1）**：速率γ₁ = 1/T₁，Lindblad算符为湮灭算符`a`
- **纯退相位（T_φ）**：速率γ_φ = 1/T_φ，Lindblad算符为粒子数算符`n`
- **关系**：1/T₂ = 1/(2T₁) + 1/T_φ

### 3.2 磁通噪声
- **1/f噪声**：功率谱密度S(f) ∝ 1/f，影响Transmon频率
- **静态偏移**：实验间的随机磁通偏移
- **频率灵敏度**：∂ω/∂Φ，在甜点处为零

### 3.3 脉冲误差
- **幅度噪声**：Rabi频率的随机波动（白噪声）
- **相位噪声**：驱动相位的随机波动
- **时间抖动**：脉冲时序的不确定性
- **有限带宽**：脉冲边沿的非理想性
- **载波频率偏移**：驱动频率与qubit频率的失谐

### 3.4 测量噪声
- **散粒噪声**：有限测量次数导致的统计涨落
- **读取误差**：误判|0⟩为|1⟩或反之

### 3.5 环境噪声
- **电荷噪声**：影响Transmon的电容能量EC
- **温度漂移**：影响约瑟夫森能量EJ
- **交叉耦合**：多qubit系统中的串扰

## 4. 详细实现方案

### 4.1 退相干噪声实现

#### 4.1.1 TransmonQubit类扩展
在`qubit.py`的`TransmonQubit`类中添加方法：

```python
def get_collapse_operators(self):
    """
    返回Lindblad算符列表，用于mesolve
    返回: [sqrt(gamma1) * a, sqrt(gamma_phi) * n]
    """
    gamma1 = 1.0 / self.T1 if self.T1 > 0 else 0
    gamma_phi = 1.0 / self.T2 - 0.5 / self.T1 if self.T2 > 0 else 0
    gamma_phi = max(gamma_phi, 0)  # 确保非负

    c_ops = []
    if gamma1 > 0:
        c_ops.append(np.sqrt(gamma1) * self.a)
    if gamma_phi > 0:
        c_ops.append(np.sqrt(gamma_phi) * self.n)
    return c_ops
```

#### 4.1.2 Protocol类修改
在`protocol.py`的`evolve`方法中，使用collapse算符：

```python
# 修改mesolve调用，添加c_ops参数
c_ops = qubit.get_collapse_operators()
result = mesolve(H_rabi, qubit.state, t_rabi, c_ops=c_ops, e_ops=[psi_e * psi_e.dag()])
```

### 4.2 磁通噪声实现

#### 4.2.1 1/f噪声生成器
在`qubit.py`中添加辅助函数：

```python
def generate_1f_noise(t_list, amplitude, f_low=1e-3, f_high=1e3):
    """
    生成1/f噪声
    :param t_list: 时间序列 (ns)
    :param amplitude: 噪声幅度 (Φ₀)
    :param f_low: 低频截止 (GHz)
    :param f_high: 高频截止 (GHz)
    :return: 噪声序列
    """
    dt = t_list[1] - t_list[0]
    n = len(t_list)
    freqs = np.fft.fftfreq(n, d=dt)

    # 生成1/f频谱
    spectrum = np.zeros(n, dtype=complex)
    for i in range(1, n//2):
        f = abs(freqs[i])
        if f_low <= f <= f_high:
            spectrum[i] = amplitude / np.sqrt(f) * (np.random.normal() + 1j*np.random.normal())

    # 对称化
    spectrum[n//2+1:] = np.conj(spectrum[1:n//2][::-1])

    # 逆傅里叶变换
    noise = np.real(np.fft.ifft(spectrum))
    return noise
```

#### 4.2.2 TransmonQubit类扩展
添加磁通噪声相关方法：

```python
def add_flux_noise(self, noise_type='1f', **kwargs):
    """
    添加磁通噪声
    :param noise_type: '1f' 或 'static'
    :param kwargs: 噪声参数
    """
    if noise_type == 'static':
        self.flux += np.random.normal(0, kwargs.get('std', 0.01))
        # 重新计算受噪声影响的参数
        self.EJ = self.EJ_0 * abs(math.cos(math.pi * self.flux))
        self.frequency = self.calculate_frequency()
        self.anharmonicity = self.calculate_anharmonicity()
        self.hamiltonian = self.get_hamiltonian()
    elif noise_type == '1f':
        # 生成时变噪声，在仿真过程中动态应用
        pass
```

#### 4.2.3 甜点计算
实现`sweet_point`方法：

```python
def sweet_point(self):
    """
    计算甜点位置（磁通噪声灵敏度最低的点）
    返回: flux值列表，可能为0或0.5
    """
    # 频率对磁通的一阶导数
    # ω = √(8EJEC) - EC, EJ = EJ₀|cos(πΦ)|
    # dω/dΦ = (4EC/√(8EJEC)) * EJ₀ * π * sin(πΦ) * sign(cos(πΦ))
    # 在Φ=0和Φ=0.5处，sin(πΦ)=0，一阶导数为零
    return [0.0, 0.5]
```

### 4.3 脉冲畸变实现

#### 4.3.1 Signal类扩展
在`signal.py`中增强`back_signal`方法：

```python
def back_signal(self, t_list=None, noise_type='gaussian', noise_level=0.01):
    """
    背景信号（可加噪声）
    :param noise_type: 'gaussian', '1f', 'uniform'
    :param noise_level: 噪声水平（相对于信号幅度）
    :return: 噪声序列
    """
    if t_list is None:
        t_list = self.t_list

    n = len(t_list)
    noise = np.zeros_like(t_list)

    if noise_type == 'gaussian':
        noise = np.random.normal(0, noise_level * self.params.get('amplitude', 1.0), n)
    elif noise_type == 'uniform':
        noise = np.random.uniform(-noise_level, noise_level, n) * self.params.get('amplitude', 1.0)
    elif noise_type == '1f':
        noise = generate_1f_noise(t_list, noise_level * self.params.get('amplitude', 1.0))

    return noise
```

#### 4.3.2 Pulse类扩展
在`pulse.py`中增强脉冲噪声处理：

```python
class Pulse:
    def __init__(self, frame=0, omega_d=0.0, phase=0.0, Omega=None,
                 is_rwa=True, qubit=None, noise_params=None):
        # ... 现有参数 ...
        self.noise_params = noise_params or {}

    def get_Rabi_frequency(self, t, include_noise=True):
        """
        获取脉冲在时间t处的Rabi频率（含噪声）
        """
        base_freq = self._get_base_Rabi_frequency(t)  # 原始方法

        if not include_noise:
            return base_freq

        # 添加幅度噪声
        if 'amplitude_noise' in self.noise_params:
            noise_level = self.noise_params['amplitude_noise']
            base_freq += np.random.normal(0, noise_level * abs(base_freq))

        # 添加时间抖动
        if 'timing_jitter' in self.noise_params:
            jitter = np.random.normal(0, self.noise_params['timing_jitter'])
            t_effective = t + jitter
            base_freq = self._get_base_Rabi_frequency(t_effective)

        return base_freq

    def apply_bandwidth_limit(self, cutoff_freq):
        """
        应用有限带宽滤波（模拟控制电子设备带宽限制）
        :param cutoff_freq: 截止频率 (GHz)
        """
        # 实现低通滤波器
        pass
```

#### 4.3.3 DRAG脉冲误差
在`qubit.py`的`simulate_gate`方法中，考虑非谐性误差：

```python
def simulate_gate(self, theta, phi, T, sigma, alpha_error=0.0):
    """
    模拟单比特门操作，考虑非谐性估计误差
    :param alpha_error: 非谐性相对误差
    """
    # 实际非谐性
    alpha_actual = self.anharmonicity * (1 + alpha_error)
    beta = -1.0 / alpha_actual  # DRAG参数

    # 后续计算使用实际非谐性
    # ...
```

### 4.4 测量噪声实现

#### 4.4.1 读取误差模型
在`analysis.py`中添加：

```python
def add_measurement_noise(self, probabilities, readout_error=(0.01, 0.02)):
    """
    添加测量噪声（读取误差）
    :param probabilities: 原始概率数组
    :param readout_error: (p(1|0), p(0|1)) 误判概率
    :return: 带噪声的概率
    """
    p0, p1 = probabilities, 1 - probabilities
    p01, p10 = readout_error

    # 误判模型：P_meas(1) = P(1)*(1-p10) + P(0)*p01
    p_meas = p1 * (1 - p10) + p0 * p01
    return p_meas

def add_shot_noise(self, probabilities, n_shots=1000):
    """
    添加散粒噪声（有限测量次数）
    :param probabilities: 真实概率
    :param n_shots: 测量次数
    :return: 带噪声的概率估计
    """
    noisy_probs = np.random.binomial(n_shots, probabilities) / n_shots
    return noisy_probs
```

### 4.5 环境噪声实现

#### 4.5.1 电荷噪声
在`TransmonQubit`类中添加：

```python
def add_charge_noise(self, noise_level=0.01):
    """
    添加电荷噪声，影响EC
    :param noise_level: EC的相对波动
    """
    EC_noisy = self.EC * (1 + np.random.normal(0, noise_level))
    self.EC = EC_noisy
    # 重新计算频率和非谐性
    self.frequency = self.calculate_frequency()
    self.anharmonicity = self.calculate_anharmonicity()
    self.hamiltonian = self.get_hamiltonian()
```

#### 4.5.2 温度漂移
添加温度依赖的EJ：

```python
def set_temperature(self, T, dEJ_dT=-0.01):
    """
    设置温度，影响约瑟夫森能量
    :param T: 温度 (K)
    :param dEJ_dT: EJ对温度的灵敏度 (%/K)
    """
    T0 = 0.010  # 基准温度 10mK
    EJ_shift = self.EJ_0 * dEJ_dT * (T - T0)
    self.EJ = self.EJ_0 + EJ_shift
    # 重新计算相关参数
    self.frequency = self.calculate_frequency()
    self.hamiltonian = self.get_hamiltonian()
```

## 5. 集成到现有代码的具体修改

### 5.1 修改文件清单

1. **qubit.py**
   - 添加`get_collapse_operators()`方法
   - 实现`flux_noise()`方法
   - 实现`sweet_point()`方法
   - 添加`generate_1f_noise()`函数
   - 添加电荷噪声和温度漂移方法

2. **signal.py**
   - 增强`back_signal()`方法，支持多种噪声类型
   - 添加噪声生成辅助函数

3. **pulse.py**
   - 在`Pulse.__init__()`中添加`noise_params`参数
   - 修改`get_Rabi_frequency()`以包含噪声
   - 添加带宽限制方法

4. **protocol.py**
   - 在`evolve()`方法中使用collapse算符
   - 在适当位置添加磁通噪声

5. **analysis.py**
   - 添加测量噪声函数
   - 添加散粒噪声函数

6. **新增noise.py模块**（可选）
   - 集中管理噪声模型
   - 提供配置接口

### 5.2 向后兼容性
所有噪声特性默认关闭，确保现有代码不受影响。通过可选参数启用噪声：

```python
# 启用退相干
qubit = TransmonQubit(..., T1=100e3, T2=50e3)

# 启用脉冲噪声
pulse = Pulse(..., noise_params={'amplitude_noise': 0.01, 'timing_jitter': 0.1})

# 启用测量噪声
analysis = Analysis()
p_noisy = analysis.add_measurement_noise(p, readout_error=(0.01, 0.02))
```

## 6. 注意事项和最佳实践

### 6.1 计算效率
1. **噪声预生成**：对于长时间仿真，预生成噪声序列避免重复计算
2. **向量化操作**：使用NumPy向量化代替循环
3. **频率截断**：1/f噪声需要合理的频率截断范围

### 6.2 物理合理性
1. **噪声相关性**：不同噪声源可能相关，需考虑联合分布
2. **参数范围**：确保噪声添加后参数仍在物理合理范围内
3. **单位一致性**：所有参数使用一致单位（GHz、ns、Φ₀）

### 6.3 数值稳定性
1. **小时间步长**：对于快速变化的噪声，需要足够小的时间步长
2. **随机数种子**：设置随机种子以保证结果可重复
3. **收敛性检查**：验证仿真结果对噪声样本的收敛性

### 6.4 实验校准
1. **噪声参数测量**：从实验数据中提取噪声参数
2. **模型验证**：与实验数据对比验证噪声模型
3. **参数扫描**：研究不同噪声水平对传感性能的影响

## 7. 测试计划

### 7.1 单元测试
1. **退相干测试**：验证T1、T2对Ramsey条纹衰减的影响
2. **磁通噪声测试**：验证甜点处的噪声抑制
3. **脉冲噪声测试**：验证门保真度随噪声水平的变化
4. **测量噪声测试**：验证反卷积算法的鲁棒性

### 7.2 集成测试
1. **完整协议测试**：运行所有传感协议，检查噪声影响
2. **参数扫描测试**：系统研究不同噪声参数的影响
3. **性能基准测试**：对比不同噪声模型的仿真速度

### 7.3 验证测试
1. **与文献对比**：验证噪声模型与已发表结果的一致性
2. **极限情况测试**：测试极端噪声条件下的仿真稳定性
3. **长期稳定性测试**：长时间仿真测试数值稳定性

## 8. 未来扩展

### 8.1 高级噪声模型
1. **非马尔可夫噪声**：超出Lindblad主方程的范围
2. **时空相关噪声**：多qubit系统中的相关噪声
3. **非线性噪声响应**：大噪声幅度的非线性效应

### 8.2 优化工具
1. **噪声识别算法**：从实验数据中自动提取噪声参数
2. **最优传感协议设计**：针对特定噪声环境的协议优化
3. **误差缓解策略**：仿真误差缓解技术的效果

### 8.3 硬件接口
1. **实验数据导入**：直接使用实验测量的噪声数据
2. **实时仿真**：与实验控制系统集成
3. **数字孪生**：创建实际量子处理器的虚拟副本

## 9. 总结

本方案系统性地提出了在超导量子比特传感仿真平台中添加真实条件的实现方法。通过分层、模块化的设计，可以逐步集成各类噪声模型，同时保持代码的清晰性和可维护性。

**核心建议**：
1. **渐进式实现**：从退相干噪声开始，逐步添加其他噪声源
2. **配置驱动**：通过配置文件或参数控制噪声特性
3. **验证导向**：每个噪声模型添加后都与物理预期对比验证
4. **文档完善**：记录每个噪声模型的物理背景和实现细节

通过实施本方案，仿真平台将能更真实地模拟实际量子传感实验，为协议设计、性能评估和误差分析提供有力工具。