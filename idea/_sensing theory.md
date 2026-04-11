# 量子传感协议理论补充

对于磁通可调Transmon qubit，其频率随刺痛的响应函数可以写为
$$
\omega_{01}(\Phi) = \sqrt{8E_C E_J(\Phi)} - E_C
,/ ,/ ,/ ,
E_D(\Phi) = E_{J0} \cos(\pi \Phi / \Phi_0)
$$
在工作点$\Phi_w$处对其做展开，得到
$$
\omega_{01}(\Phi) \approx \omega_{01}(\Phi_w) + \left. \frac{d\omega_{01}}{d\Phi} \right|_{\Phi_w} (\Phi - \Phi_w) + \frac{1}{2} \left. \frac{d^2\omega_{01}}{d\Phi^2} \right|_{\Phi_w} (\Phi - \Phi_w)^2 + ...
$$
取一阶近似，得到
$$
\delta \omega = \omega_{01}(\Phi) - \omega_{01}(\Phi_w) \approx \kappa_{\Phi} \delta \Phi
$$
其中$\kappa_{\Phi} = \left. \frac{d\omega_{01}}{d\Phi} \right|_{\Phi_w}$是频率对磁通的灵敏度。

在旋转系下，考虑RWA近似，则系统哈密顿量可以写为
$$
H(t) = \frac{1}{2} \delta \omega(t) \sigma_z + H_{ctrl}(t)
$$
其中$H_{ctrl}(t)$是控制哈密顿量，包含了我们施加的脉冲，在自由演化时，控制哈密顿量为零。对于不同的协议，$H_{ctrl}(t)$的形式不同。控制脉冲包含$\pi/2$脉冲和$\pi$脉冲，其中$\pi/2$脉冲的作用是将量子态从$|0\rangle$旋转到$|+\rangle$（或其他方向），其通常作用在脉冲序列的开始或结尾，用于将初态制备为叠加态，或将末态投影到测量态，得到积累的相位信息，$\pi$脉冲的作用是将量子态从$|+\rangle$旋转到$|-\rangle$（或其他方向），其会将已经积累的相位变换为负值，并在后续演化中以变换后的相反方向继续积累相位。

通过以上的描述，可以得到自由演化过程积累的相位可以写为
$$
\phi = \int_0^{T}y(t)\delta \omega(t) dt
$$
其中$y(t)$是调制函数，取决于我们施加的控制脉冲序列。

演化末态可以写为
$$
\ket{\psi} = \frac{1}{\sqrt{2}} (|0\rangle + e^{i\phi} |1\rangle)
$$
施加$y-\pi/2$脉冲，测量得到的激发态概率为
$$
P_e = \frac{1}{2} (1 - \cos\phi)
$$


## Ramsey干涉
在Ramsey干涉协议中，$y(t) = 1$，时间范围为0到T，对其作fourier变换，得到频域的调制函数为因此积累的相位为
$$
Y(\omega) = \int_0^T y(t) e^{-i\omega t} dt = e^{-i\omega T/2} \frac{\sin(\omega T/2)}{\omega/2}
$$
根据频域函数，可以看到$Y(0) = T$最大，且第一个零点出现在$\omega = 2\pi/T$，因此Ramsey干涉协议中，测量时间越长， 带宽越窄，但是对低频信号响应强度更大。

积累的相位为
$$
\phi_{Ramsey} = \int_0^{T} \delta \omega(t) dt
$$
通过改变演化时间$T$，得到激发态概率随演化时间的关系，从而可以拟合得到相位关于演化时间的关系，随后对相位进行数值微分，可以得到给定时刻的频率偏移$\delta \omega$，进而得到磁场的瞬时值。

如果磁场信号是一个恒定信号，则可以直接对概率曲线进行拟合。但是如果磁场随时间变化，则需要先还原出相位曲线，但是由于概率曲线是一个周期函数，因此理论上的还原出的相位曲线是一个多值函数，因此需要保证相位足够小，使其能落入主值区域，例如缩短演化时间。或者对相位进行unwrap。

另外，对phi进行数值微分会引入额外噪声，因此在实际操作中，可能需要对phi进行平滑处理，或者直接对概率曲线进行拟合，得到一个关于时间的函数，然后对这个函数进行微分，得到频率偏移。

## Spin Echo
在Spin Echo协议中，$y(t)$是一个分段函数，在前半段为1，在后半段为-1，作fourier变换，得到频域的调制函数为
$$
Y(\omega) = \int_0^T y(t) e^{-i\omega t} dt = e^{i\omega T/2} \frac{i\sin^2(\omega T/4)}{\omega/4}
$$
根据频域函数，可以看到$Y(0) = 0, |Y(2\pi/T)| = \frac{T}{4\pi}$，因此Spin Echo协议中，低频分量被抑制，而对频率为$2\pi/T$的分量响应较强。
因此积累的相位为
$$
\phi(T = 2\tau) = \int_0^{\tau} \delta \omega(t) dt - \int_{\tau}^{2\tau} \delta \omega(t) dt
$$

根据相位定义，可以看到当磁场为0时， 相位为0，因此Spin Echo协议可以抵消掉慢变的噪声，从而提高对快速变化的信号的灵敏度。

spin echo 主要用于测量ac磁场，如$\kappa\Phi = A\cos(\omega t + \varphi)$，代入相位表达式，得到
$$
\phi = \int_0^{\tau} A\cos(\omega t + \varphi) dt - \int_{\tau}^{2\tau} A\cos(\omega t + \varphi) dt = \frac{4A}{\omega} \sin^2(\frac{\omega \tau}{2}) \sin(\omega \tau + \varphi)
$$
概率函数满足
$$
p(T) = \frac{1}{2} (1 - C\cos\phi + \phi_0) \
$$
若磁场信号频率已知，则可以将取$\tau\approx \frac{\pi}{\omega}$，此时echo对于频率的响应最大，随后固定$\tau$，通过扫描分析脉冲的相位$\phi_0$，可以得到概率关于$\phi_0$的函数。从而还原出相位$\phi$。

若磁场频率未知，则需要扫描$\tau$，随后扫描$\phi_0$，还原得到$\phi$关于$\tau$的函数，取峰值点未知，得到对应的最大响应$\tau$，从而重复上述过程，得到磁场频率的估计值。

## CPMG协议
在CPMG协议中，$y(t)$是一个分段函数，设第$j$个$\pi$脉冲的时间为$t_j = (j - \frac{1}{2})\tau, t_0 = 0, j = 1,2,...,N, N \geq 2$，则在$t_{2k}$和$t_{2k+1}$区间内，$y(t) = 1$，作fourier变换，得到频域的调制函数为
$$
\begin{aligned}
Y(\omega) &= \int_0^T y(t) e^{-i\omega t} dt \\  &= \sum_{k = 0}^{[(N+1)/2]-1}\int_{t_{2k}}^{t_{2k+1}} e^{-i\omega t} dt - \sum_{k = 0}^{[N/2]-1}\int_{t_{2k+1}}^{t_{2k+2}} e^{-i\omega t} dt \\ &= \frac{1}{i\omega} \left( 1 + (-1)^{N+1} e^{-i\omega T} + 2\sum_{k = 1}^{[N/2]} (-1)^k e^{-i\omega t_k} \right)
\end{aligned}
$$
则主响应峰在$\omega = \frac{\pi}{\tau}$，缝宽为$\frac{2\pi}{\tau}$，$\omega = 0$时，频域函数为0。则CPMG通过增加N，使得主响应峰更尖锐，从而提高对特定频率信号的灵敏度，同时抑制掉更多的低频噪声。相较于spinecho，其可以实现在中心频率不变的情况下，尽可能地缩短带宽。

若待测磁场为AC信号，同样可以计算计算出积累的相位为
$$
\phi = \int_0^T y(t) A\cos(\omega t + \varphi) dt 
$$

在实际测量时，若磁场频率已知，则取$\tau \approx \frac{\pi}{\omega}$，$N$取足够大，随后扫描分析脉冲的相位$\phi_0$，得到概率关于$\phi_0$的函数，从而还原出相位$\phi$。

若磁场为任意信号，对其做fourier展开，得到
$$
B(t) = \sum_{n = -\infty}^{+\infty} A_n \cos(\omega_n t + \varphi_n)
$$


基于上述傅里叶展开，CPMG协议可通过改变脉冲间隔 $\tau$ 实现对不同频率分量的选择性测量，从而重建磁场时域信号。

### 1. CPMG作为可调谐带通滤波器

CPMG序列的调制函数 $y_\tau(t)$ 在频域的响应 $Y_\tau(\omega)$ 具有带通特性，其中心频率 $\omega_c$ 与脉冲间隔 $\tau$ 满足：
$$
\omega_c = \frac{\pi}{\tau}
$$
当 $\tau$ 变化时，带通滤波器的中心频率随之移动。

对于包含多频率分量的磁场 $B(t) = \sum_n A_n \cos(\omega_n t + \varphi_n)$，积累的相位为：
$$
\phi(\tau) = \int_0^T y_\tau(t) B(t) dt = \sum_n A_n \int_0^T y_\tau(t) \cos(\omega_n t + \varphi_n) dt
$$
定义复响应函数：
$$
\tilde{Y}_\tau(\omega) = \int_0^T y_\tau(t) e^{-i\omega t} dt
$$
其实部对应余弦分量的响应。相位可表达为：
$$
\phi(\tau) = \sum_n A_n |\tilde{Y}_\tau(\omega_n)| \cos\left(\varphi_n + \arg\tilde{Y}_\tau(\omega_n)\right)
$$
因此，对于每个固定的 $\tau$，测量得到的相位 $\phi(\tau)$ 是各频率分量经 $Y_\tau(\omega)$ 加权滤波后的叠加。

### 2. 频率扫描采样

通过系统性地改变 $\tau$，可以获得磁场频谱的采样值。设扫描的 $\tau$ 序列为 $\{\tau_1, \tau_2, \dots, \tau_M\}$，对应的中心频率为 $\{\omega_{c1}, \omega_{c2}, \dots, \omega_{cM}\}$，其中 $\omega_{cj} = \pi/\tau_j$。

测量得到相位序列 $\phi(\tau_j)$，近似满足：
$$
\phi(\tau_j) \approx \sum_n A_n Y_{\tau_j}(\omega_n) \cos(\varphi_n + \delta_{jn})
$$
其中 $\delta_{jn} = \arg Y_{\tau_j}(\omega_n)$。

### 3. 磁场重建算法

重建磁场时域信号可视为解卷积问题。下面给出详细的数学框架和算法步骤，可直接用于代码实现。

#### 3.1 离散化与频率网格

将连续频谱离散化为 $N_f$ 个频率分量。设频率网格为 $\{\omega_1, \omega_2, \dots, \omega_{N_f}\}$，其范围由 $\tau$ 扫描范围决定：

- 最低频率：$\omega_{\min} = \pi / \tau_{\max}$
- 最高频率：$\omega_{\max} = \pi / \tau_{\min}$

频率点数 $N_f$ 的选择需权衡分辨率与计算复杂度，通常 $N_f \geq M$（$M$ 为 $\tau$ 采样点数）。可采用线性或对数间隔：

**线性网格**（适用于均匀带宽）：
$$
\omega_n = \omega_{\min} + \frac{n-1}{N_f-1}(\omega_{\max} - \omega_{\min}), \quad n=1,\dots,N_f
$$

**对数网格**（适用于宽频带测量）：
$$
\omega_n = \exp\left(\ln\omega_{\min} + \frac{n-1}{N_f-1}(\ln\omega_{\max} - \ln\omega_{\min})\right)
$$

#### 3.2 响应矩阵构建

对于每个 $\tau_j$ 和 $\omega_n$，计算 CPMG 调制函数的傅里叶变换值。根据前文推导，CPMG 序列的频域响应函数为：

$$
Y_\tau(\omega) = \frac{1}{i\omega} \left( 1 + (-1)^{N+1} e^{-i\omega T} + 2\sum_{k=1}^{[N/2]} (-1)^k e^{-i\omega t_k} \right)
$$

其中 $t_k = (k-1/2)\tau$，$T = N\tau$，$N$ 为 $\pi$ 脉冲数。

**响应矩阵 $\mathbf{Y} \in \mathbb{C}^{M \times N_f}$** 的元素定义为：
$$
Y_{jn} = Y_{\tau_j}(\omega_n), \quad j=1,\dots,M; \ n=1,\dots,N_f
$$

由于测量相位 $\phi(\tau_j)$ 是实数，实际参与运算的是 $\mathbf{Y}$ 的实部。定义 **实响应矩阵 $\mathbf{Y}_r = \Re\{\mathbf{Y}\}$**。

#### 3.3 线性系统建模

设复频谱向量 $\mathbf{a} \in \mathbb{C}^{N_f}$，其中 $a_n = A_n e^{i\varphi_n}$ 包含幅值 $A_n$ 和相位 $\varphi_n$。测量模型为：

$$
\boldsymbol{\phi} = \Re\{\mathbf{Y} \mathbf{a}\} + \boldsymbol{\epsilon}
$$

其中 $\boldsymbol{\phi} = [\phi(\tau_1), \dots, \phi(\tau_M)]^T \in \mathbb{R}^M$，$\boldsymbol{\epsilon} \in \mathbb{R}^M$ 为测量噪声。

将复数运算转换为实数运算：令 $\mathbf{a} = \mathbf{x} + i\mathbf{y}$，其中 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^{N_f}$ 分别为实部和虚部。则：

$$
\Re\{\mathbf{Y} \mathbf{a}\} = \Re\{\mathbf{Y}\}\mathbf{x} - \Im\{\mathbf{Y}\}\mathbf{y} = \mathbf{Y}_r \mathbf{x} - \mathbf{Y}_i \mathbf{y}
$$

其中 $\mathbf{Y}_i = \Im\{\mathbf{Y}\}$。定义增广向量 $\mathbf{z} = [\mathbf{x}^T, \mathbf{y}^T]^T \in \mathbb{R}^{2N_f}$ 和增广矩阵 $\mathbf{\tilde{Y}} = [\mathbf{Y}_r, -\mathbf{Y}_i] \in \mathbb{R}^{M \times 2N_f}$，则线性系统简化为：

$$
\boldsymbol{\phi} = \mathbf{\tilde{Y}} \mathbf{z} + \boldsymbol{\epsilon}
$$

#### 3.4 正则化求解

直接求解 $\mathbf{z}$ 通常病态（尤其当 $M < 2N_f$ 时），需引入正则化。采用 **Tikhonov 正则化**：

$$
\hat{\mathbf{z}} = \arg\min_{\mathbf{z}} \left\{ \|\mathbf{\tilde{Y}} \mathbf{z} - \boldsymbol{\phi}\|_2^2 + \lambda \|\mathbf{z}\|_2^2 \right\}
$$

其中 $\lambda > 0$ 为正则化参数。该优化问题的解析解为：

$$
\hat{\mathbf{z}} = (\mathbf{\tilde{Y}}^T \mathbf{\tilde{Y}} + \lambda \mathbf{I})^{-1} \mathbf{\tilde{Y}}^T \boldsymbol{\phi}
$$

**正则化参数选择**：
- **交叉验证**：将数据分为训练集和验证集，选择使验证误差最小的 $\lambda$
- **L-曲线法**：绘制 $\|\mathbf{\tilde{Y}}\mathbf{z}-\boldsymbol{\phi}\|$ 与 $\|\mathbf{z}\|$ 的关系曲线，选择拐点处的 $\lambda$
- **经验公式**：$\lambda = \sigma^2 / \|\mathbf{a}_{\text{prior}}\|^2$，其中 $\sigma^2$ 为噪声方差，$\|\mathbf{a}_{\text{prior}}\|^2$ 为先验信号能量估计

#### 3.5 频谱与时域重建

从 $\hat{\mathbf{z}} = [\hat{\mathbf{x}}^T, \hat{\mathbf{y}}^T]^T$ 恢复复频谱：
$$
\hat{a}_n = \hat{x}_n + i\hat{y}_n, \quad n=1,\dots,N_f
$$

幅值和相位分别为：
$$
\hat{A}_n = |\hat{a}_n|, \quad \hat{\varphi}_n = \arg(\hat{a}_n)
$$

重建的时域磁场信号为：
$$
\hat{B}(t) = \sum_{n=1}^{N_f} \hat{A}_n \cos(\omega_n t + \hat{\varphi}_n)
$$

#### 3.6 算法步骤总结

1. **数据采集**：
   - 扫描 $\tau$ 序列 $\{\tau_1, \dots, \tau_M\}$
   - 对每个 $\tau_j$，通过双点测量（分析脉冲相位 $\phi_0=0$ 和 $\pi/2$）得到 $\phi(\tau_j)$
   - 收集测量向量 $\boldsymbol{\phi} \in \mathbb{R}^M$

2. **参数设置**：
   - 选择频率网格 $\{\omega_1, \dots, \omega_{N_f}\}$（线性或对数）
   - 设定正则化参数 $\lambda$（通过交叉验证等）
   - 确定 CPMG 脉冲数 $N$（通常固定）

3. **矩阵计算**：
   - 对每个 $(\tau_j, \omega_n)$ 计算 $Y_{\tau_j}(\omega_n)$
   - 构建复数矩阵 $\mathbf{Y} \in \mathbb{C}^{M \times N_f}$
   - 分离实部 $\mathbf{Y}_r$ 和虚部 $\mathbf{Y}_i$
   - 构造增广矩阵 $\mathbf{\tilde{Y}} = [\mathbf{Y}_r, -\mathbf{Y}_i]$

4. **求解频谱**：
   - 解正则化方程：$\hat{\mathbf{z}} = (\mathbf{\tilde{Y}}^T \mathbf{\tilde{Y}} + \lambda \mathbf{I})^{-1} \mathbf{\tilde{Y}}^T \boldsymbol{\phi}$
   - 提取 $\hat{\mathbf{x}}, \hat{\mathbf{y}}$，组合成 $\hat{\mathbf{a}} = \hat{\mathbf{x}} + i\hat{\mathbf{y}}$

5. **信号重建**：
   - 计算幅值 $\hat{A}_n = |\hat{a}_n|$ 和相位 $\hat{\varphi}_n = \arg(\hat{a}_n)$
   - 合成时域信号：$\hat{B}(t) = \sum_n \hat{A}_n \cos(\omega_n t + \hat{\varphi}_n)$

#### 3.7 伪代码示例（Python）

```python
import numpy as np
from scipy.linalg import solve

def reconstruct_magnetic_field(tau_list, phi_list, N_pulses, lambda_reg=1e-3):
    """
    从CPMG扫描数据重建磁场
    
    参数:
        tau_list: τ扫描序列 (M,)
        phi_list: 测量的相位序列 (M,)
        N_pulses: CPMG中的π脉冲数
        lambda_reg: 正则化参数
    
    返回:
        omega_grid: 频率网格 (N_f,)
        A_est: 幅值估计 (N_f,)
        phi_est: 相位估计 (N_f,)
        B_recon: 重建的时域信号函数 B(t)
    """
    M = len(tau_list)
    
    # 1. 设置频率网格
    omega_min = np.pi / max(tau_list)
    omega_max = np.pi / min(tau_list)
    N_f = min(2*M, 200)  # 频率点数，经验值
    omega_grid = np.linspace(omega_min, omega_max, N_f)
    
    # 2. 计算CPMG响应函数
    def Y_tau_omega(tau, omega, N):
        """计算Y_tau(omega)"""
        T = N * tau
        t_k = np.array([(k - 0.5) * tau for k in range(1, N//2 + 1)])
        sum_term = 2 * np.sum([(-1)**k * np.exp(-1j * omega * t_k[k-1]) 
                              for k in range(1, N//2 + 1)])
        return (1 + (-1)**(N+1) * np.exp(-1j * omega * T) + sum_term) / (1j * omega)
    
    # 3. 构建响应矩阵
    Y_complex = np.zeros((M, N_f), dtype=complex)
    for j, tau in enumerate(tau_list):
        for n, omega in enumerate(omega_grid):
            Y_complex[j, n] = Y_tau_omega(tau, omega, N_pulses)
    
    Y_real = np.real(Y_complex)
    Y_imag = np.imag(Y_complex)
    
    # 4. 构建增广系统
    Y_tilde = np.hstack([Y_real, -Y_imag])  # (M, 2N_f)
    phi_vec = np.array(phi_list).reshape(-1, 1)  # (M, 1)
    
    # 5. 正则化求解
    I = np.eye(2 * N_f)
    z_est = solve(Y_tilde.T @ Y_tilde + lambda_reg * I, Y_tilde.T @ phi_vec)
    z_est = z_est.flatten()
    
    # 6. 提取频谱
    x_est = z_est[:N_f]  # 实部
    y_est = z_est[N_f:]  # 虚部
    a_est = x_est + 1j * y_est
    
    A_est = np.abs(a_est)
    phi_est = np.angle(a_est)
    
    # 7. 定义重建函数
    def B_reconstructed(t):
        """重建的磁场时域信号"""
        return np.sum(A_est * np.cos(omega_grid * t[:, np.newaxis] + phi_est), axis=1)
    
    return omega_grid, A_est, phi_est, B_reconstructed

# 使用示例
# tau_scan = np.linspace(10, 200, 50)  # 50个τ点，10-200 ns
# phi_measured = [...]  # 实验测量的相位
# omega, A, phi, B_func = reconstruct_magnetic_field(tau_scan, phi_measured, N_pulses=4)
# t_eval = np.linspace(0, 1000, 1000)
# B_recon = B_func(t_eval)
```

#### 3.8 注意事项

1. **矩阵条件数**：响应矩阵 $\mathbf{\tilde{Y}}$ 通常病态，需确保 $\lambda$ 足够大以稳定求解
2. **频率外推**：重建仅覆盖 $[\omega_{\min}, \omega_{\max}]$ 频带，该范围外的分量无法恢复
3. **计算效率**：直接求解 $(M \times 2N_f)$ 系统复杂度为 $O((2N_f)^3)$，对于大 $N_f$ 可考虑迭代法
4. **噪声模型**：上述假设高斯白噪声，若噪声有色需采用加权正则化

### 4. 实际考虑因素

#### 4.1 频率分辨率
- 最小可分辨频率 $\Delta f_{\min} \approx 1/T$，其中 $T = N\tau$ 为总演化时间。
- 最大可测频率 $f_{\max} \approx 1/(2\tau_{\min})$，受 $\tau$ 下限限制。

#### 4.2 噪声与正则化
- 测量噪声会通过解卷积放大，需采用 Tikhonov 正则化或 Wiener 滤波。
- 正则化参数可通过交叉验证确定。

#### 4.3 采样策略
- $\tau$ 采样应覆盖感兴趣频带，可按 $\omega_c = \pi/\tau$ 对数或线性分布。
- 对于宽带信号，建议先粗扫定位主要频率成分，再细扫提高分辨率。

#### 4.4 相位提取精度
- 需通过扫描分析脉冲相位 $\phi_0$ 精确提取 $\phi(\tau)$，方法同前：
  $$
  P_e(\phi_0) = \frac{1}{2}\left[1 - \cos(\phi(\tau) + \phi_0)\right]
  $$
  测量 $\phi_0 = 0$ 和 $\phi_0 = \pi/2$ 两点即可解出 $\phi(\tau)$。

### 5. 总结

CPMG频率扫描方案本质上是利用脉冲序列的**可调谐滤波特性**，通过改变 $\tau$ 移动带通滤波器的中心频率，从而采样磁场的频率成分。该方法适用于：
1. **未知频率的AC磁场测量**
2. **宽带磁场信号的频谱分析**
3. **时变磁场的非侵入式监测**

主要优势包括频率选择性好、抗低频噪声能力强；主要限制是频率分辨率受总演化时间 $T$ 约束，且对直流磁场不敏感。

该方案为量子磁强计提供了一种灵活的频域测量手段，可结合压缩感知等先进算法进一步提高重建效率。


## 差分回波协议
该协议采用等效采样计技术，下图描述了4中具体脉冲序列：
![alt text](image.png)
其中，纵坐标$B(t)$表示一个可重复的任意磁场信号，横坐标$t$表示相对于触发点的时间，$t_{rep}$表示两次采样之间的间隔时间，即重复时间。

b展示了标准Ramsey协议的脉冲序列，通过改变延迟时间可以得到不同的磁场从0到$\tau$的积分，对最终结果进行数值微分可以得到磁场的瞬时值。但是由于概率函数是周期函数，通过这种方法求相位会导致相位wrap，并且求导会引入额外的噪声。

c展示了一种Ramsey干涉的变形，在$t$时刻施加一个脉冲，演化一个小的时间间隔$t_{int}$后，施加另一个脉冲，在这种脉冲序列下，由于演化时间很短，于是相位可以直接写为
$$
\phi = \int_t^{t+t_{int}} \delta \omega(t') dt' \approx \kappa_{\Phi} t_{int} B(t)
$$
由此可以直接得到磁场的瞬时值，而不需要进行数值微分，同时由于相位足够小，因此不存在wrap的问题。但是由于间隔时间短，因此积累的相位较小，测量结果的信噪比较低，灵敏度较低。

d展示了差分回波协议的脉冲序列，在初始时施加一个$\pi/2$脉冲，随后在$t$时刻施加一个$\pi$脉冲，待信号结束，再次触发，产生信号，随后在$t+t_{int}$时刻施加另一个$\pi$脉冲，重复$k$上述过程，最后施加一个$\pi/2$脉冲。

以上的序列可以统一用一个调制函数$M(t, t')$描述，$t$表示采样延迟，即当前想要测量的磁场信号的时间点，相位可以写为
$$
\phi = k\int_{0}^{2t_{rep}} M(t, t') \delta \omega(t') dt' = k\int_{0}^{t_{rep}}( M_1(t, t') + M_2(t, t') ) \kappa_{\Phi} B(t') dt'
$$
其中考虑了磁场的可重复性。上式中，$M_1 + M_2$的作用等价于图c脉冲的作用。则
$$
\phi = -2k\kappa_{\Phi} t_{int} B(t)
$$

上式通过一个脉冲序列，自动实现了对磁场的微分，从而直接得到磁场的瞬时值，同时由于积累了$k$次相位，因此可以提高信噪比，提升灵敏度。


## Walsh协议
