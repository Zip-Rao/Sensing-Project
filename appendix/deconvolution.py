"""
Wiener反卷积模块
实现从观测信号中恢复原始磁场信号的功能
基于Wiener反卷积算法，用于量子传感信号处理
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve
from numpy.fft import fft, ifft


class WienerDeconvolution:
    """
    Wiener反卷积类
    实现频域Wiener滤波器，用于从含噪观测信号中恢复原始磁场信号
    """

    def __init__(self):
        """初始化Wiener反卷积对象"""
        pass

    def generate_mock_kernel(self, t_axis, center, width):
        """
        生成模拟核函数 k(t)
        使用高斯函数模拟典型的Ramsey/Spin-Echo类型核函数

        Parameters:
        -----------
        t_axis : array_like
            时间轴 (ns)
        center : float
            核函数中心位置 (ns)
        width : float
            核函数宽度 (ns, FWHM)

        Returns:
        --------
        k : ndarray
            核函数值
        """
        # 模拟一个类似文中Fig. 1c或2d的核函数（中心为峰值，两侧衰减）
        # 使用高斯函数近似
        sigma = width / 2.355  # FWHM转sigma
        k = np.exp(-(t_axis - center)**2 / (2 * sigma**2))

        # 归一化：通常核函数的峰值由控制脉冲决定，这里设为1用于演示
        return k

    def transient_field(self, t, t0, tau_rise, tau_decay, amplitude):
        """
        生成非对称瞬态磁场信号（例如：双指数衰减）

        Parameters:
        -----------
        t : array_like
            时间序列 (ns)
        t0 : float
            磁场起始时间 (ns)
        tau_rise : float
            上升时间常数 (ns)
        tau_decay : float
            衰减时间常数 (ns)
        amplitude : float
            磁场幅值

        Returns:
        --------
        B : ndarray
            瞬态磁场信号
        """
        val = np.zeros_like(t)
        mask = t >= t0

        # 简单的瞬态模型：(1 - exp(-t/rise)) * exp(-t/decay)
        t_rel = t[mask] - t0
        val[mask] = amplitude * (1 - np.exp(-t_rel/tau_rise)) * np.exp(-t_rel/tau_decay)

        return val

    def wiener_deconvolution(self, signal_p, kernel_k, lambda_reg, dt):
        """
        实现Wiener反卷积算法

        公式: B_hat(ω) = [K*(ω) / (|K(ω)|² + λ²)] * P(ω)

        Parameters:
        -----------
        signal_p : array_like
            观测信号 p(t)（含噪声）
        kernel_k : array_like
            传感核函数 k(t)
        lambda_reg : float
            正则化参数，控制平滑与噪声抑制的权衡
        dt : float
            时间步长 (ns)

        Returns:
        --------
        B_rec : ndarray
            恢复的磁场信号 B_hat(t)
        """
        # 1. 转换到频域
        # 注意：为了避免FFT引入的线性相位（时间平移），我们需要将核函数的峰值移到数组起点
        # 找到核函数最大值位置
        k_argmax = np.argmax(kernel_k)
        # 循环移位，使峰值位于index 0
        k_shifted = np.roll(kernel_k, -k_argmax)

        # 傅里叶变换，乘以dt对应连续傅里叶变换的数值近似
        K_w = fft(k_shifted) * dt
        P_w = fft(signal_p) * dt

        # 2. 构建Wiener滤波器
        # C(ω) = K*(ω) / (|K(ω)|² + λ²)
        # λ扮演正则化参数角色
        numerator = np.conj(K_w)
        denominator = (np.abs(K_w)**2 + lambda_reg**2)

        # 避免除零
        epsilon = 1e-12
        denominator = np.where(denominator < epsilon, epsilon, denominator)

        Wiener_filter = numerator / denominator

        # 3. 应用滤波器并逆变换
        B_w = Wiener_filter * P_w
        B_rec = ifft(B_w) / dt  # 除以dt还原

        # 4. 移位回原始时间参考
        B_rec_shifted = np.roll(np.real(B_rec), k_argmax)

        return B_rec_shifted

    def run_simulation(self, use_real_kernel=False, kernel_data=None):
        """
        运行完整的Wiener反卷积仿真

        Parameters:
        -----------
        use_real_kernel : bool
            是否使用真实核函数（从传感仿真获得）
        kernel_data : tuple or None
            如果use_real_kernel为True，提供(t_list, kernel)数据

        Returns:
        --------
        results : dict
            包含所有仿真结果的字典
        """
        # ==========================================
        # 0. 时间轴设置
        # ==========================================
        t_total = 100  # ns
        dt = 0.2       # ns
        t_list = np.arange(0, t_total, dt)
        n_points = len(t_list)

        # ==========================================
        # 1. 获取或生成核函数 k(t)
        # ==========================================
        if use_real_kernel and kernel_data is not None:
            # 使用真实核函数
            t_kernel, k_t = kernel_data
            # 确保时间轴匹配，可能需要插值
            if len(t_kernel) != len(t_list):
                # 简单线性插值（实际使用时可能需要更精确的插值）
                from scipy.interpolate import interp1d
                interp_func = interp1d(t_kernel, k_t, kind='linear',
                                       bounds_error=False, fill_value=0)
                k_t = interp_func(t_list)
            print("使用真实核函数")
        else:
            # 生成模拟核函数
            kernel_center = t_total / 2
            kernel_width = 8.0  # ns (FWHM)
            k_t = self.generate_mock_kernel(t_list, kernel_center, kernel_width)
            print("使用模拟核函数")

        # ==========================================
        # 2. 生成模拟瞬态磁场 B(t)
        # ==========================================
        # 模拟一个快速的磁场脉冲
        # 假设 gamma = 1 (归一化单位)，我们直接处理 gamma * B
        B_true = self.transient_field(t_list, t0=40, tau_rise=2.0,
                                      tau_decay=5.0, amplitude=1.0)

        # ==========================================
        # 3. 正向过程：卷积生成观测概率 p(t)
        # ==========================================
        # 理论公式: p(t) = k(t) * B(t)
        # discrete convolution需要乘dt来近似连续积分
        # mode='same'保证输出长度一致
        p_clean = convolve(B_true, k_t, mode='same') * dt

        # 添加测量噪声 n(t)
        noise_level = 0.02  # 模拟2%的读出噪声
        np.random.seed(42)  # 固定随机种子以确保可重复性
        noise = np.random.normal(0, noise_level, size=len(t_list))
        p_measured = p_clean + noise

        # ==========================================
        # 4. 逆向过程：Wiener反卷积重构
        # ==========================================
        # 测试不同正则化参数λ
        lambdas = [0.5, 0.1, 5e-4]
        colors = ['green', 'orange', 'red']
        reconstructed_signals = []

        for lam in lambdas:
            B_rec = self.wiener_deconvolution(p_measured, k_t, lam, dt)
            reconstructed_signals.append(B_rec)

        # ==========================================
        # 5. 计算重构误差
        # ==========================================
        reconstruction_errors = []
        for B_rec in reconstructed_signals:
            # 计算均方根误差 (RMSE)
            rmse = np.sqrt(np.mean((B_rec - B_true)**2))
            reconstruction_errors.append(rmse)

        # ==========================================
        # 6. 收集所有结果
        # ==========================================
        results = {
            't_list': t_list,
            'B_true': B_true,
            'k_t': k_t,
            'p_clean': p_clean,
            'p_measured': p_measured,
            'lambdas': lambdas,
            'reconstructed_signals': reconstructed_signals,
            'reconstruction_errors': reconstruction_errors,
            'colors': colors,
            'dt': dt
        }

        return results

    def visualize_results(self, results, save_path=None):
        """
        可视化Wiener反卷积结果

        Parameters:
        -----------
        results : dict
            run_simulation返回的结果字典
        save_path : str or None
            图像保存路径，如果为None则不保存
        """
        # 提取数据
        t_list = results['t_list']
        B_true = results['B_true']
        k_t = results['k_t']
        p_clean = results['p_clean']
        p_measured = results['p_measured']
        lambdas = results['lambdas']
        reconstructed_signals = results['reconstructed_signals']
        reconstruction_errors = results['reconstruction_errors']
        colors = results['colors']

        # 创建图形
        plt.figure(figsize=(12, 10))

        # 子图1: 核函数与原始磁场
        plt.subplot(3, 1, 1)
        plt.plot(t_list, k_t / np.max(k_t), 'b-', linewidth=2,
                 label='Sensing Kernel $k(t)$ (Normalized)')
        plt.plot(t_list, B_true, 'k-', linewidth=2,
                 label='True Field $B(t)$')
        plt.title('Input Signals: Kernel and True Field')
        plt.ylabel('Amplitude')
        plt.legend()
        plt.grid(True, alpha=0.3)

        # 子图2: 正向过程（卷积与噪声）
        plt.subplot(3, 1, 2)
        plt.plot(t_list, B_true, 'k-', linewidth=2, label='True Field $B(t)$')
        plt.plot(t_list, p_measured, 'b.', markersize=2, alpha=0.5,
                 label='Measured $p(t)$ (Noisy)')
        plt.plot(t_list, p_clean, 'c--', linewidth=1.5,
                 label='Ideal $p(t)$ (Smoothed)')
        plt.title('Forward Process: Convolution & Noise')
        plt.ylabel('Amplitude')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.text(5, 0.8, "Observation:\n$p(t)$ is broader/smoother\nthan $B(t)$ due to $k(t)$ width",
                 bbox=dict(facecolor='white', alpha=0.8))

        # 子图3: Wiener反卷积重构
        plt.subplot(3, 1, 3)
        plt.plot(t_list, B_true, 'k-', linewidth=2, alpha=0.6,
                 label='True Field $B(t)$')

        for i, lam in enumerate(lambdas):
            label_str = f'Wiener $\lambda={lam:.4f}$ (RMSE={reconstruction_errors[i]:.4f})'
            if i == 0:
                label_str += ' (Under-fit / Smooth)'
            elif i == 2:
                label_str += ' (Over-fit / Noisy)'
            else:
                label_str += ' (Good fit)'

            plt.plot(t_list, reconstructed_signals[i], color=colors[i],
                     linewidth=1.5, label=label_str)

        plt.title('Inverse Process: Wiener Deconvolution Reconstruction')
        plt.xlabel('Time (ns)')
        plt.ylabel('Reconstructed Field')
        plt.legend()
        plt.grid(True, alpha=0.3)

        plt.tight_layout()

        # 保存或显示图像
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"图像已保存至: {save_path}")

        plt.show()

        # 打印误差摘要
        print("\n" + "="*50)
        print("Wiener反卷积结果摘要")
        print("="*50)
        for i, lam in enumerate(lambdas):
            print(f"λ = {lam:.4f}: RMSE = {reconstruction_errors[i]:.6f}")
        print("="*50)


def run_deconvolution_demo():
    """
    运行Wiener反卷积演示
    """
    print("开始Wiener反卷积仿真...")

    # 创建反卷积对象
    deconv = WienerDeconvolution()

    # 运行仿真（使用模拟核函数）
    results = deconv.run_simulation(use_real_kernel=False)

    # 可视化结果
    deconv.visualize_results(results, save_path="wiener_deconvolution_results.png")

    return results


if __name__ == "__main__":
    run_deconvolution_demo()