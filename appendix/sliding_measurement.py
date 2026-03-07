"""
滑动测量仿真模块
实现图片中所示的协议：通过滑动磁场信号测量激发态概率
使用现有仿真平台的所有类，不修改原始代码
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# 添加src目录到Python路径
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./src'))

# 导入现有平台模块
from src.qubit import TransmonQubit
from src.signal import Signal
from src.pulse import Pulse, CompositePulse
from src.analysis import Analysis
# 导入之前创建的Wiener反卷积模块
try:
    from deconvolution import WienerDeconvolution
    WIENER_AVAILABLE = True
except ImportError:
    WIENER_AVAILABLE = False
    print("警告: Wiener反卷积模块未找到，将使用内置简化版本")


class SlidingMeasurement:
    """
    滑动测量类
    实现图片中所示的协议：通过滑动磁场信号B(t)相对于控制脉冲序列
    测量每个延迟时间τ下的激发态概率pe(τ)
    """

    def __init__(self, qubit, control_pulse, magnetic_signal, sensitivity=2*np.pi*1.0):
        """
        初始化滑动测量对象

        Parameters:
        -----------
        qubit : TransmonQubit
            量子比特对象
        control_pulse : CompositePulse or Pulse
            控制脉冲序列（使用CompositePulse或单个Pulse）
        magnetic_signal : Signal
            待测磁场信号B(t)
        sensitivity : float
            磁场灵敏度 γ (rad/ns/T)，默认 2π * 1.0
        """
        self.qubit = qubit
        self.control_pulse = control_pulse
        self.magnetic_signal = magnetic_signal
        self.sensitivity = sensitivity

        # 确保控制脉冲是CompositePulse类型
        if isinstance(control_pulse, Pulse):
            # 将单个Pulse转换为CompositePulse
            self.control_pulse = CompositePulse([control_pulse])

        # 获取控制脉冲的时间轴
        self.t_list = np.array(self.control_pulse.t_list)
        if len(self.t_list) < 2:
            raise ValueError("控制脉冲时间轴过短")

        self.dt = self.t_list[1] - self.t_list[0]
        self.total_time = self.t_list[-1]

        # 分析工具
        self.analysis = Analysis()

        # 结果存储
        self.delays = None
        self.pe_ideal = None
        self.pe_measured = None
        self.noise_level = 0.0

    def create_ramsey_pulse_sequence(self, pulse_duration=40, pi_over_2_angle=np.pi/2):
        """
        创建Ramsey脉冲序列（π/2 - 等待 - π/2）
        用于演示目的

        Parameters:
        -----------
        pulse_duration : float
            每个脉冲的持续时间 (ns)
        pi_over_2_angle : float
            π/2脉冲的角度 (rad)

        Returns:
        --------
        composite_pulse : CompositePulse
            Ramsey脉冲序列
        """
        # 第一个π/2脉冲
        t_pulse1 = np.linspace(0, pulse_duration, 100)
        # 使用恒定信号作为简单脉冲
        signal1 = Signal(
            type=1,  # 恒定信号
            t_list=t_pulse1,
            amplitude=pi_over_2_angle / pulse_duration  # 简单近似
        )

        pulse1 = Pulse(
            frame=1,  # 旋转系
            omega_d=0.0,  # 共振驱动
            phase=0.0,  # X轴旋转
            Omega=signal1,
            is_rwa=True,
            qubit=self.qubit
        )

        # 等待间隔（零脉冲）
        wait_time = 20  # ns
        t_wait = np.linspace(0, wait_time, 50)
        signal_wait = Signal(type=0, t_list=t_wait)  # 零信号

        pulse_wait = Pulse(
            frame=1,
            omega_d=0.0,
            phase=0.0,
            Omega=signal_wait,
            is_rwa=True,
            qubit=self.qubit
        )

        # 第二个π/2脉冲
        t_pulse2 = np.linspace(0, pulse_duration, 100)
        signal2 = Signal(
            type=1,
            t_list=t_pulse2,
            amplitude=pi_over_2_angle / pulse_duration
        )

        pulse2 = Pulse(
            frame=1,
            omega_d=0.0,
            phase=0.0,
            Omega=signal2,
            is_rwa=True,
            qubit=self.qubit
        )

        # 创建复合脉冲序列
        composite_pulse = CompositePulse([pulse1, pulse_wait, pulse2])

        # 更新控制脉冲
        self.control_pulse = composite_pulse
        self.t_list = np.array(composite_pulse.t_list)
        self.dt = self.t_list[1] - self.t_list[0]
        self.total_time = self.t_list[-1]

        return composite_pulse

    def run_single_measurement(self, delay, add_noise=False):
        """
        运行单次测量：给定延迟τ，计算激发态概率pe(τ)

        Parameters:
        -----------
        delay : float
            磁场信号相对于控制脉冲序列的延迟 (ns)
        add_noise : bool
            是否添加测量噪声

        Returns:
        --------
        p_e : float
            激发态概率
        """
        # 定义磁场信号的哈密顿量项
        # H_sense = 0.5 * γ * B(t - τ) * σ_z
        # 注意：使用sigmaz()作为两能级系统的Z算符

        # 获取量子比特的能级数
        n_levels = self.qubit.n_levels

        # 创建Z算符（适用于多能级系统）
        if n_levels == 2:
            sz = sigmaz()
        else:
            # 对于多能级系统，使用num算符的线性组合
            n = num(n_levels)
            sz = 2*n - qeye(n_levels)

        # 定义磁场信号函数（考虑延迟）
        def magnetic_signal_func(t, args):
            # t是仿真时间（从0开始）
            # args包含delay和magnetic_signal
            delay = args['delay']
            signal = args['magnetic_signal']
            sensitivity = args['sensitivity']

            # 计算全局时间：t_global = t - delay
            t_global = t - delay

            # 获取磁场信号值
            B_value = signal.value_at(t_global)

            # 返回哈密顿量系数：0.5 * γ * B(t-τ)
            return 0.5 * sensitivity * B_value

        # 控制脉冲的哈密顿量
        def control_hamiltonian(t, args):
            pulse = args['control_pulse']
            return pulse.get_hamiltonian(t)

        # 静态哈密顿量（在旋转系下）
        if self.control_pulse.frame == 0:
            # 实验系
            H0 = self.qubit.get_hamiltonian(self.qubit.frequency)
        else:
            # 旋转系（假设共振）
            H0 = self.qubit.get_hamiltonian_rwa(self.qubit.frequency)

        # 总哈密顿量
        H = [H0,
             [control_hamiltonian, {'control_pulse': self.control_pulse}],
             [sz, magnetic_signal_func]]

        # 参数
        args = {
            'delay': delay,
            'magnetic_signal': self.magnetic_signal,
            'sensitivity': self.sensitivity,
            'control_pulse': self.control_pulse
        }

        # 初始态（基态）
        psi0 = self.qubit.state

        # 演化算符（用于测量激发态概率）
        psi_e = basis(n_levels, 1)
        e_ops = [psi_e * psi_e.dag()]

        # 执行演化
        result = mesolve(H, psi0, self.t_list, [], e_ops, args=args,
                         options=Options(nsteps=5000, atol=1e-10, rtol=1e-10))

        # 获取最终时刻的激发态概率
        p_e = result.expect[0][-1]

        # 添加测量噪声（如果启用）
        if add_noise and self.noise_level > 0:
            noise = np.random.normal(0, self.noise_level)
            p_e = max(0, min(1, p_e + noise))  # 限制在[0,1]范围内

        return p_e

    def run_sliding_scan(self, delay_start=-20, delay_end=80, n_delays=141,
                         noise_level=0.02, verbose=True):
        """
        运行滑动扫描：在延迟范围内测量激发态概率

        Parameters:
        -----------
        delay_start : float
            起始延迟 (ns)
        delay_end : float
            结束延迟 (ns)
        n_delays : int
            延迟点数
        noise_level : float
            测量噪声水平（标准差）
        verbose : bool
            是否显示进度信息

        Returns:
        --------
        delays : ndarray
            延迟时间数组
        pe_ideal : ndarray
            理想激发态概率（无噪声）
        pe_measured : ndarray
            测量到的激发态概率（含噪声）
        """
        self.noise_level = noise_level

        # 生成延迟数组
        delays = np.linspace(delay_start, delay_end, n_delays)
        self.delays = delays

        pe_ideal = []
        pe_measured = []

        if verbose:
            print(f"开始滑动扫描，共 {n_delays} 个延迟点...")
            print(f"延迟范围: [{delay_start}, {delay_end}] ns")
            print(f"噪声水平: {noise_level}")

        # 固定随机种子以确保可重复性
        np.random.seed(42)

        for i, delay in enumerate(delays):
            if verbose and i % 20 == 0:
                print(f"  进度: {i+1}/{n_delays} (delay={delay:.1f} ns)")

            # 计算理想概率（无噪声）
            p_ideal = self.run_single_measurement(delay, add_noise=False)
            pe_ideal.append(p_ideal)

            # 计算含噪声的概率
            p_measured = self.run_single_measurement(delay, add_noise=True)
            pe_measured.append(p_measured)

        self.pe_ideal = np.array(pe_ideal)
        self.pe_measured = np.array(pe_measured)

        if verbose:
            print("滑动扫描完成!")
            print(f"理想概率范围: [{np.min(self.pe_ideal):.4f}, {np.max(self.pe_ideal):.4f}]")
            print(f"测量概率范围: [{np.min(self.pe_measured):.4f}, {np.max(self.pe_measured):.4f}]")

        return delays, self.pe_ideal, self.pe_measured

    def get_kernel(self):
        """
        获取控制脉冲的传感核函数

        Returns:
        --------
        t_kernel : ndarray
            核函数时间轴
        kernel : ndarray
            核函数值
        """
        # 使用Analysis类的get_kernel方法
        t_kernel, kernel = self.analysis.get_kernel(self.control_pulse, self.qubit)

        # 核函数可能需要重新对齐（使峰值在中心）
        # 这里简单返回原始结果
        return t_kernel, kernel

    def reconstruct_magnetic_field(self, lambdas=[0.01, 0.1, 1.0, 10.0],
                                   use_wiener_module=True):
        """
        使用Wiener反卷积从测量数据重建磁场信号

        Parameters:
        -----------
        lambdas : list
            正则化参数λ列表
        use_wiener_module : bool
            是否使用外部Wiener反卷积模块

        Returns:
        --------
        results : dict
            重建结果字典，包含不同λ的重建信号
        """
        if self.pe_measured is None:
            raise ValueError("请先运行滑动扫描 (run_sliding_scan)")

        # 获取核函数
        t_kernel, kernel = self.get_kernel()

        # 确保核函数和测量数据有相同的时间分辨率
        # 这里假设delays和t_kernel有相同的dt
        dt_delays = self.delays[1] - self.delays[0]
        dt_kernel = t_kernel[1] - t_kernel[0] if len(t_kernel) > 1 else dt_delays

        # 如果时间分辨率不同，需要插值
        if abs(dt_delays - dt_kernel) > 1e-10 and len(t_kernel) > 1:
            from scipy.interpolate import interp1d
            interp_func = interp1d(t_kernel, kernel, kind='cubic',
                                  bounds_error=False, fill_value=0)
            # 在delays时间轴上插值核函数
            kernel_interp = interp_func(self.delays)
        else:
            kernel_interp = kernel

        # 归一化核函数（可选）
        if np.max(np.abs(kernel_interp)) > 0:
            kernel_interp = kernel_interp / np.max(np.abs(kernel_interp))

        results = {
            'delays': self.delays,
            'pe_measured': self.pe_measured,
            'pe_ideal': self.pe_ideal,
            'kernel_original': (t_kernel, kernel),
            'kernel_interpolated': kernel_interp,
            'reconstructed_signals': {},
            'lambdas': lambdas
        }

        if use_wiener_module and WIENER_AVAILABLE:
            # 使用外部Wiener反卷积模块
            deconv = WienerDeconvolution()

            for lam in lambdas:
                # 使用deconvolution模块中的wiener_deconvolution方法
                B_rec = deconv.wiener_deconvolution(
                    self.pe_measured, kernel_interp, lam, dt_delays
                )
                results['reconstructed_signals'][lam] = B_rec

        else:
            # 使用内置简化版本的Wiener反卷积
            print("使用内置简化Wiener反卷积")

            for lam in lambdas:
                B_rec = self._simple_wiener_deconvolution(
                    self.pe_measured, kernel_interp, lam, dt_delays
                )
                results['reconstructed_signals'][lam] = B_rec

        return results

    def _simple_wiener_deconvolution(self, measured_sig, kernel_sig, lamb, dt):
        """
        简化版本的Wiener反卷积实现

        Parameters:
        -----------
        measured_sig : array_like
            测量信号 p(t)
        kernel_sig : array_like
            核函数 k(t)
        lamb : float
            正则化参数 λ
        dt : float
            时间步长

        Returns:
        --------
        B_rec : ndarray
            重建的磁场信号
        """
        from numpy.fft import fft, ifft

        # 1. 转换到频域
        P_omega = fft(measured_sig) * dt
        K_omega = fft(kernel_sig) * dt

        # 2. 应用Wiener滤波器公式
        # B(ω) = P(ω) * K*(ω) / (|K(ω)|² + λ²)
        numerator = P_omega * np.conj(K_omega)
        denominator = (np.abs(K_omega)**2 + lamb**2)

        # 避免除零
        epsilon = 1e-12
        denominator = np.where(denominator < epsilon, epsilon, denominator)

        B_omega = numerator / denominator

        # 3. 转换回时域
        B_rec = np.real(ifft(B_omega)) / dt

        return B_rec

    def visualize_results(self, reconstruction_results=None, save_path=None):
        """
        可视化滑动测量和重建结果

        Parameters:
        -----------
        reconstruction_results : dict or None
            reconstruct_magnetic_field返回的结果
        save_path : str or None
            图像保存路径
        """
        # 创建图形
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # 获取真实磁场信号（用于对比）
        # 注意：这里需要从magnetic_signal中提取
        B_true_values = []
        for t in self.delays:
            # 在t时刻的磁场值（无延迟）
            B_val = self.magnetic_signal.value_at(t)
            B_true_values.append(B_val)

        B_true = np.array(B_true_values)

        # 子图1: 控制脉冲序列
        ax = axes[0, 0]
        # 绘制Rabi频率包络线
        t_pulse = self.t_list
        Omega_values = [self.control_pulse.get_Omega(t) for t in t_pulse]

        ax.plot(t_pulse, Omega_values, 'b-', linewidth=1.5)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Rabi Frequency (GHz)')
        ax.set_title('Control Pulse Sequence')
        ax.grid(True, alpha=0.3)

        # 子图2: 磁场信号
        ax = axes[0, 1]
        # 绘制原始磁场信号
        t_signal = np.linspace(self.delays[0], self.delays[-1], 500)
        B_signal = [self.magnetic_signal.value_at(t) for t in t_signal]

        ax.plot(t_signal, B_signal, 'r-', linewidth=2, label='True B(t)')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Magnetic Field (arb.)')
        ax.set_title('Magnetic Field Signal')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 子图3: 滑动测量结果
        ax = axes[1, 0]
        ax.plot(self.delays, self.pe_ideal, 'k--', linewidth=1.5,
                label='Ideal P_e (no noise)', alpha=0.7)
        ax.plot(self.delays, self.pe_measured, 'b.', markersize=3,
                label='Measured P_e (with noise)')
        ax.set_xlabel('Delay τ (ns)')
        ax.set_ylabel('Excited State Probability')
        ax.set_title('Sliding Measurement Results')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 子图4: Wiener反卷积重建结果
        ax = axes[1, 1]

        if reconstruction_results is not None:
            # 绘制真实磁场信号
            ax.plot(self.delays, B_true / np.max(np.abs(B_true)),
                   'k--', linewidth=2, label='True B(t) (Normalized)', alpha=0.6)

            # 绘制不同λ的重建结果
            lambdas = reconstruction_results['lambdas']
            colors = ['green', 'orange', 'red', 'purple']

            for i, lam in enumerate(lambdas):
                if lam in reconstruction_results['reconstructed_signals']:
                    B_rec = reconstruction_results['reconstructed_signals'][lam]
                    # 归一化以便对比
                    if np.max(np.abs(B_rec)) > 0:
                        B_rec_norm = B_rec / np.max(np.abs(B_rec))
                    else:
                        B_rec_norm = B_rec

                    color_idx = i % len(colors)
                    ax.plot(self.delays, B_rec_norm, color=colors[color_idx],
                           linewidth=1.5, label=f'Reconstructed (λ={lam})')

        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Normalized Amplitude')
        ax.set_title('Wiener Deconvolution Reconstruction')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # 保存或显示图像
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"图像已保存至: {save_path}")

        plt.show()

        # 如果有核函数，也单独绘制
        if reconstruction_results is not None:
            t_kernel, kernel = reconstruction_results['kernel_original']

            plt.figure(figsize=(8, 4))
            plt.plot(t_kernel, kernel, 'b-', linewidth=2)
            plt.xlabel('Time (ns)')
            plt.ylabel('Kernel Amplitude')
            plt.title('Sensing Kernel k(t)')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()

            if save_path:
                kernel_save_path = save_path.replace('.png', '_kernel.png')
                plt.savefig(kernel_save_path, dpi=150, bbox_inches='tight')

            plt.show()


def run_sliding_measurement_demo():
    """
    运行滑动测量演示
    """
    print("="*60)
    print("滑动测量协议演示")
    print("="*60)

    # ==========================================
    # 1. 创建量子比特对象
    # ==========================================
    print("\n创建量子比特对象...")
    qubit = TransmonQubit(
        EC=0.2,      # 电容能量 (GHz)
        EJ=10.0,     # 约瑟夫森能量 (GHz)
        T1=100.0e3,  # 100 μs
        T2=50.0e3,   # 50 μs
        flux=0.0,    # 外加磁通
        state=0,     # 初始态 |0>
        n_levels=2   # 两能级系统
    )

    # ==========================================
    # 2. 创建控制脉冲序列 (Ramsey序列)
    # ==========================================
    print("\n创建控制脉冲序列 (Ramsey π/2 - wait - π/2)...")
    sliding_meas = SlidingMeasurement(
        qubit=qubit,
        control_pulse=None,  # 稍后设置
        magnetic_signal=None,  # 稍后设置
        sensitivity=2*np.pi*1.0  # γ = 2π * 1.0 rad/ns/T
    )

    # 创建Ramsey脉冲序列
    control_pulse = sliding_meas.create_ramsey_pulse_sequence(
        pulse_duration=20,  # 每个脉冲20ns
        pi_over_2_angle=np.pi/2
    )

    # ==========================================
    # 3. 创建磁场信号 (高斯脉冲)
    # ==========================================
    print("\n创建磁场信号 (高斯脉冲)...")
    # 创建时间轴（用于磁场信号）
    t_total = 100  # ns
    t_magnetic = np.linspace(-50, 150, 1000)  # 更宽的时间范围

    magnetic_signal = Signal(
        type=3,  # 高斯脉冲
        t_list=t_magnetic,
        amplitude=1.0,  # 幅度
        center=50.0,    # 中心在50ns
        width=8.0       # 宽度8ns
    )

    # ==========================================
    # 4. 更新滑动测量对象
    # ==========================================
    sliding_meas.magnetic_signal = magnetic_signal

    # ==========================================
    # 5. 运行滑动扫描
    # ==========================================
    print("\n运行滑动扫描...")
    delays, pe_ideal, pe_measured = sliding_meas.run_sliding_scan(
        delay_start=-20,   # 起始延迟
        delay_end=80,      # 结束延迟
        n_delays=141,      # 点数
        noise_level=0.02,  # 2%噪声
        verbose=True
    )

    # ==========================================
    # 6. 获取核函数
    # ==========================================
    print("\n获取传感核函数...")
    t_kernel, kernel = sliding_meas.get_kernel()
    print(f"核函数点数: {len(kernel)}")
    print(f"核函数最大值: {np.max(kernel):.6f}")

    # ==========================================
    # 7. Wiener反卷积重建
    # ==========================================
    print("\n运行Wiener反卷积重建...")
    reconstruction_results = sliding_meas.reconstruct_magnetic_field(
        lambdas=[0.01, 0.1, 1.0, 10.0],  # 测试不同λ值
        use_wiener_module=WIENER_AVAILABLE
    )

    # ==========================================
    # 8. 可视化结果
    # ==========================================
    print("\n生成可视化结果...")
    sliding_meas.visualize_results(
        reconstruction_results=reconstruction_results,
        save_path="sliding_measurement_results.png"
    )

    print("\n" + "="*60)
    print("滑动测量协议演示完成!")
    print("="*60)

    return sliding_meas, reconstruction_results


if __name__ == "__main__":
    run_sliding_measurement_demo()