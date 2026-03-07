"""
瞬态磁场测量协议演示
实现图片中所示的完整量子传感协议流程
基于现有仿真平台的所有模块，不修改原始代码
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# 添加src目录到Python路径
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./src'))

# 导入现有平台模块
from src.qubit import TransmonQubit
from src.signal import Signal
from src.pulse import Pulse, CompositePulse
from src.analysis import Analysis
from src.protocal import Protocal

# 导入自定义模块
try:
    from appendix.sliding_measurement import SlidingMeasurement
    SLIDING_AVAILABLE = True
except ImportError:
    SLIDING_AVAILABLE = False
    print("警告: 滑动测量模块未找到")

try:
    from deconvolution import WienerDeconvolution
    WIENER_AVAILABLE = True
except ImportError:
    WIENER_AVAILABLE = False
    print("警告: Wiener反卷积模块未找到")


class TransientFieldProtocol:
    """
    瞬态磁场测量协议类
    实现图片中所示的完整协议流程
    """

    def __init__(self, qubit_params=None, protocol_params=None):
        """
        初始化瞬态磁场测量协议

        Parameters:
        -----------
        qubit_params : dict
            量子比特参数
        protocol_params : dict
            协议参数
        """
        # 默认量子比特参数
        default_qubit_params = {
            'EC': 0.2,      # 电容能量 (GHz)
            'EJ': 10.0,     # 约瑟夫森能量 (GHz)
            'T1': 100.0e3,  # 100 μs
            'T2': 50.0e3,   # 50 μs
            'flux': 0.0,    # 外加磁通
            'state': 0,     # 初始态 |0>
            'n_levels': 2   # 两能级系统
        }

        # 默认协议参数
        default_protocol_params = {
            'sensitivity': 2*np.pi*1.0,  # 磁场灵敏度 γ (rad/ns/T)
            'pulse_duration': 20.0,      # 脉冲持续时间 (ns)
            'wait_time': 20.0,           # 等待时间 (ns)
            'delay_start': -20.0,        # 起始延迟 (ns)
            'delay_end': 80.0,           # 结束延迟 (ns)
            'n_delays': 141,             # 延迟点数
            'noise_level': 0.02,         # 测量噪声水平
            'magnetic_signal_center': 50.0,  # 磁场信号中心 (ns)
            'magnetic_signal_width': 8.0,    # 磁场信号宽度 (ns)
            'magnetic_signal_amplitude': 1.0,# 磁场信号幅度
            'wiener_lambdas': [0.01, 0.1, 1.0, 10.0]  # Wiener反卷积λ值
        }

        # 合并参数
        if qubit_params is not None:
            default_qubit_params.update(qubit_params)
        self.qubit_params = default_qubit_params

        if protocol_params is not None:
            default_protocol_params.update(protocol_params)
        self.protocol_params = default_protocol_params

        # 创建量子比特对象
        self.qubit = self._create_qubit()

        # 创建控制脉冲序列
        self.control_pulse = self._create_control_pulse()

        # 创建磁场信号
        self.magnetic_signal = self._create_magnetic_signal()

        # 创建分析工具
        self.analysis = Analysis()

        # 结果存储
        self.sliding_results = None
        self.kernel_results = None
        self.reconstruction_results = None

    def _create_qubit(self):
        """创建量子比特对象"""
        qubit = TransmonQubit(
            EC=self.qubit_params['EC'],
            EJ=self.qubit_params['EJ'],
            T1=self.qubit_params['T1'],
            T2=self.qubit_params['T2'],
            flux=self.qubit_params['flux'],
            state=self.qubit_params['state'],
            n_levels=self.qubit_params['n_levels']
        )
        return qubit

    def _create_control_pulse(self):
        """
        创建控制脉冲序列（Ramsey序列：π/2 - 等待 - π/2）

        Returns:
        --------
        composite_pulse : CompositePulse
            复合脉冲序列
        """
        pulse_duration = self.protocol_params['pulse_duration']
        wait_time = self.protocol_params['wait_time']
        pi_over_2_angle = np.pi/2

        # 第一个π/2脉冲
        t_pulse1 = np.linspace(0, pulse_duration, 100)
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

        return composite_pulse

    def _create_magnetic_signal(self):
        """
        创建待测磁场信号（高斯脉冲）

        Returns:
        --------
        magnetic_signal : Signal
            磁场信号
        """
        # 创建时间轴（更宽的范围以支持滑动）
        center = self.protocol_params['magnetic_signal_center']
        width = self.protocol_params['magnetic_signal_width']
        amplitude = self.protocol_params['magnetic_signal_amplitude']

        # 时间范围覆盖滑动延迟范围
        delay_start = self.protocol_params['delay_start']
        delay_end = self.protocol_params['delay_end']

        # 扩展时间范围以确保滑动时信号完整
        t_start = delay_start - 3*width
        t_end = delay_end + 3*width
        t_magnetic = np.linspace(t_start, t_end, 1000)

        magnetic_signal = Signal(
            type=3,  # 高斯脉冲
            t_list=t_magnetic,
            amplitude=amplitude,
            center=center,
            width=width
        )

        return magnetic_signal

    def run_sliding_measurement(self, use_sliding_module=True):
        """
        运行滑动测量

        Parameters:
        -----------
        use_sliding_module : bool
            是否使用滑动测量模块

        Returns:
        --------
        results : dict
            滑动测量结果
        """
        print("="*60)
        print("步骤1: 滑动测量")
        print("="*60)

        if use_sliding_module and SLIDING_AVAILABLE:
            # 使用滑动测量模块
            print("使用滑动测量模块...")

            sliding_meas = SlidingMeasurement(
                qubit=self.qubit,
                control_pulse=self.control_pulse,
                magnetic_signal=self.magnetic_signal,
                sensitivity=self.protocol_params['sensitivity']
            )

            delays, pe_ideal, pe_measured = sliding_meas.run_sliding_scan(
                delay_start=self.protocol_params['delay_start'],
                delay_end=self.protocol_params['delay_end'],
                n_delays=self.protocol_params['n_delays'],
                noise_level=self.protocol_params['noise_level'],
                verbose=True
            )

            self.sliding_results = {
                'delays': delays,
                'pe_ideal': pe_ideal,
                'pe_measured': pe_measured,
                'sliding_object': sliding_meas
            }

        else:
            # 手动实现滑动测量（简化版）
            print("手动实现滑动测量...")

            # 这里可以使用之前创建的SlidingMeasurement类的逻辑
            # 为了简化，这里直接调用SlidingMeasurement类
            from appendix.sliding_measurement import SlidingMeasurement

            sliding_meas = SlidingMeasurement(
                qubit=self.qubit,
                control_pulse=self.control_pulse,
                magnetic_signal=self.magnetic_signal,
                sensitivity=self.protocol_params['sensitivity']
            )

            delays, pe_ideal, pe_measured = sliding_meas.run_sliding_scan(
                delay_start=self.protocol_params['delay_start'],
                delay_end=self.protocol_params['delay_end'],
                n_delays=self.protocol_params['n_delays'],
                noise_level=self.protocol_params['noise_level'],
                verbose=True
            )

            self.sliding_results = {
                'delays': delays,
                'pe_ideal': pe_ideal,
                'pe_measured': pe_measured,
                'sliding_object': sliding_meas
            }

        print(f"滑动测量完成! 测量点数: {len(delays)}")
        return self.sliding_results

    def extract_kernel(self):
        """
        提取传感核函数

        Returns:
        --------
        results : dict
            核函数结果
        """
        print("\n" + "="*60)
        print("步骤2: 提取传感核函数")
        print("="*60)

        # 使用Analysis类的get_kernel方法
        t_kernel, kernel = self.analysis.get_kernel(self.control_pulse, self.qubit)

        print(f"核函数提取完成! 点数: {len(kernel)}")
        print(f"核函数最大值: {np.max(kernel):.6f}")
        print(f"核函数最小值: {np.min(kernel):.6f}")

        self.kernel_results = {
            't_kernel': t_kernel,
            'kernel': kernel
        }

        return self.kernel_results

    def run_wiener_deconvolution(self):
        """
        运行Wiener反卷积重建磁场信号

        Returns:
        --------
        results : dict
            重建结果
        """
        print("\n" + "="*60)
        print("步骤3: Wiener反卷积重建")
        print("="*60)

        if self.sliding_results is None:
            raise ValueError("请先运行滑动测量 (run_sliding_measurement)")

        if self.kernel_results is None:
            raise ValueError("请先提取核函数 (extract_kernel)")

        delays = self.sliding_results['delays']
        pe_measured = self.sliding_results['pe_measured']
        t_kernel = self.kernel_results['t_kernel']
        kernel = self.kernel_results['kernel']

        # 确保核函数和测量数据有相同的时间分辨率
        dt_delays = delays[1] - delays[0]
        if len(t_kernel) > 1:
            dt_kernel = t_kernel[1] - t_kernel[0]
        else:
            dt_kernel = dt_delays

        # 如果时间分辨率不同，需要插值
        if abs(dt_delays - dt_kernel) > 1e-10 and len(t_kernel) > 1:
            from scipy.interpolate import interp1d
            interp_func = interp1d(t_kernel, kernel, kind='cubic',
                                  bounds_error=False, fill_value=0)
            kernel_interp = interp_func(delays)
        else:
            kernel_interp = kernel

        # 归一化核函数（可选）
        if np.max(np.abs(kernel_interp)) > 0:
            kernel_interp = kernel_interp / np.max(np.abs(kernel_interp))

        # 获取真实磁场信号（用于对比）
        B_true = np.array([self.magnetic_signal.value_at(t) for t in delays])

        # 运行Wiener反卷积
        lambdas = self.protocol_params['wiener_lambdas']
        reconstructed_signals = {}

        if WIENER_AVAILABLE:
            # 使用Wiener反卷积模块
            deconv = WienerDeconvolution()

            for lam in lambdas:
                print(f"  运行λ={lam}的反卷积...")
                B_rec = deconv.wiener_deconvolution(
                    pe_measured, kernel_interp, lam, dt_delays
                )
                reconstructed_signals[lam] = B_rec

        else:
            # 使用简化版本
            from numpy.fft import fft, ifft

            for lam in lambdas:
                print(f"  运行λ={lam}的反卷积...")

                # 简化Wiener反卷积
                P_omega = fft(pe_measured) * dt_delays
                K_omega = fft(kernel_interp) * dt_delays

                numerator = P_omega * np.conj(K_omega)
                denominator = (np.abs(K_omega)**2 + lam**2)

                epsilon = 1e-12
                denominator = np.where(denominator < epsilon, epsilon, denominator)

                B_omega = numerator / denominator
                B_rec = np.real(ifft(B_omega)) / dt_delays

                reconstructed_signals[lam] = B_rec

        # 计算重建误差
        reconstruction_errors = {}
        for lam, B_rec in reconstructed_signals.items():
            # 计算归一化相关系数
            if np.std(B_rec) > 0 and np.std(B_true) > 0:
                correlation = np.corrcoef(B_rec, B_true)[0, 1]
            else:
                correlation = 0

            # 计算归一化均方根误差
            if np.max(np.abs(B_true)) > 0:
                B_true_norm = B_true / np.max(np.abs(B_true))
                B_rec_norm = B_rec / np.max(np.abs(B_rec))
                rmse = np.sqrt(np.mean((B_rec_norm - B_true_norm)**2))
            else:
                rmse = 0

            reconstruction_errors[lam] = {
                'correlation': correlation,
                'rmse': rmse
            }

        self.reconstruction_results = {
            'delays': delays,
            'pe_measured': pe_measured,
            'B_true': B_true,
            't_kernel': t_kernel,
            'kernel': kernel,
            'kernel_interpolated': kernel_interp,
            'reconstructed_signals': reconstructed_signals,
            'reconstruction_errors': reconstruction_errors,
            'lambdas': lambdas
        }

        print("Wiener反卷积完成!")
        print("重建误差:")
        for lam in lambdas:
            err = reconstruction_errors[lam]
            print(f"  λ={lam}: 相关系数={err['correlation']:.4f}, RMSE={err['rmse']:.6f}")

        return self.reconstruction_results

    def visualize_full_protocol(self, save_path=None):
        """
        可视化完整协议流程（类似图片中的图示）

        Parameters:
        -----------
        save_path : str or None
            图像保存路径
        """
        print("\n" + "="*60)
        print("步骤4: 可视化完整协议流程")
        print("="*60)

        if (self.sliding_results is None or
            self.kernel_results is None or
            self.reconstruction_results is None):
            raise ValueError("请先运行完整协议流程")

        # 创建图形（4个子图，模拟图片中的布局）
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))

        # 子图1: 控制脉冲序列
        ax = axes[0, 0]
        t_pulse = np.array(self.control_pulse.t_list)
        Omega_values = [self.control_pulse.get_Omega(t) for t in t_pulse]

        ax.plot(t_pulse, Omega_values, 'b-', linewidth=1.5)
        ax.fill_between(t_pulse, Omega_values, alpha=0.3, color='blue')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Rabi Frequency (GHz)')
        ax.set_title('(a) Control Pulse Sequence')
        ax.grid(True, alpha=0.3)

        # 子图2: 磁场信号
        ax = axes[0, 1]
        delays = self.sliding_results['delays']
        B_true = self.reconstruction_results['B_true']

        ax.plot(delays, B_true, 'r-', linewidth=2)
        ax.fill_between(delays, B_true, alpha=0.3, color='red')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Magnetic Field (arb.)')
        ax.set_title('(b) True Magnetic Field B(t)')
        ax.grid(True, alpha=0.3)

        # 子图3: 传感核函数
        ax = axes[0, 2]
        t_kernel = self.kernel_results['t_kernel']
        kernel = self.kernel_results['kernel']

        ax.plot(t_kernel, kernel, 'g-', linewidth=2)
        ax.fill_between(t_kernel, kernel, alpha=0.3, color='green')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Kernel Amplitude')
        ax.set_title('(c) Sensing Kernel k(t)')
        ax.grid(True, alpha=0.3)

        # 子图4: 滑动测量结果
        ax = axes[1, 0]
        pe_ideal = self.sliding_results['pe_ideal']
        pe_measured = self.sliding_results['pe_measured']

        ax.plot(delays, pe_ideal, 'k--', linewidth=1.5,
                label='Ideal P_e', alpha=0.7)
        ax.plot(delays, pe_measured, 'b.', markersize=3,
                label='Measured P_e (with noise)')
        ax.set_xlabel('Delay τ (ns)')
        ax.set_ylabel('Excited State Probability')
        ax.set_title('(d) Sliding Measurement P_e(τ)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 子图5: Wiener反卷积重建结果（最佳λ）
        ax = axes[1, 1]
        # 找到最佳λ（相关系数最高）
        best_lambda = None
        best_correlation = -1

        for lam, errors in self.reconstruction_results['reconstruction_errors'].items():
            if errors['correlation'] > best_correlation:
                best_correlation = errors['correlation']
                best_lambda = lam

        if best_lambda is not None:
            B_rec_best = self.reconstruction_results['reconstructed_signals'][best_lambda]

            # 归一化以便对比
            if np.max(np.abs(B_true)) > 0:
                B_true_norm = B_true / np.max(np.abs(B_true))
            else:
                B_true_norm = B_true

            if np.max(np.abs(B_rec_best)) > 0:
                B_rec_norm = B_rec_best / np.max(np.abs(B_rec_best))
            else:
                B_rec_norm = B_rec_best

            ax.plot(delays, B_true_norm, 'k--', linewidth=2,
                   label='True B(t)', alpha=0.6)
            ax.plot(delays, B_rec_norm, 'r-', linewidth=1.5,
                   label=f'Reconstructed (λ={best_lambda})')

            errors = self.reconstruction_results['reconstruction_errors'][best_lambda]
            ax.text(0.05, 0.95,
                   f"Correlation: {errors['correlation']:.4f}\nRMSE: {errors['rmse']:.6f}",
                   transform=ax.transAxes, fontsize=10,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Normalized Amplitude')
        ax.set_title('(e) Wiener Deconvolution Reconstruction')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 子图6: 不同λ的重建误差对比
        ax = axes[1, 2]
        lambdas = self.reconstruction_results['lambdas']
        correlations = [self.reconstruction_results['reconstruction_errors'][lam]['correlation']
                       for lam in lambdas]
        rmses = [self.reconstruction_results['reconstruction_errors'][lam]['rmse']
                for lam in lambdas]

        # 双y轴图
        ax1 = ax
        color1 = 'tab:blue'
        ax1.set_xlabel('Regularization Parameter λ')
        ax1.set_ylabel('Correlation', color=color1)
        ax1.plot(lambdas, correlations, 'o-', color=color1, linewidth=2)
        ax1.tick_params(axis='y', labelcolor=color1)
        ax1.set_xscale('log')

        ax2 = ax1.twinx()
        color2 = 'tab:red'
        ax2.set_ylabel('RMSE', color=color2)
        ax2.plot(lambdas, rmses, 's-', color=color2, linewidth=2)
        ax2.tick_params(axis='y', labelcolor=color2)

        ax.set_title('(f) Reconstruction Error vs λ')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # 保存或显示图像
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"完整协议流程图已保存至: {save_path}")

        plt.show()

        # 单独保存最佳重建结果
        if save_path and best_lambda is not None:
            fig2, ax2 = plt.subplots(figsize=(8, 5))

            # 归一化对比
            B_rec_best = self.reconstruction_results['reconstructed_signals'][best_lambda]

            if np.max(np.abs(B_true)) > 0:
                B_true_norm = B_true / np.max(np.abs(B_true))
            else:
                B_true_norm = B_true

            if np.max(np.abs(B_rec_best)) > 0:
                B_rec_norm = B_rec_best / np.max(np.abs(B_rec_best))
            else:
                B_rec_norm = B_rec_best

            ax2.plot(delays, B_true_norm, 'k-', linewidth=2,
                    label='True Magnetic Field', alpha=0.7)
            ax2.plot(delays, B_rec_norm, 'r-', linewidth=1.5,
                    label=f'Reconstructed (λ={best_lambda})')

            errors = self.reconstruction_results['reconstruction_errors'][best_lambda]
            ax2.text(0.05, 0.95,
                    f"Correlation: {errors['correlation']:.4f}\nRMSE: {errors['rmse']:.6f}",
                    transform=ax2.transAxes, fontsize=10,
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            ax2.set_xlabel('Time (ns)')
            ax2.set_ylabel('Normalized Amplitude')
            ax2.set_title(f'Best Reconstruction (λ={best_lambda})')
            ax2.legend()
            ax2.grid(True, alpha=0.3)

            best_save_path = save_path.replace('.png', '_best_reconstruction.png')
            plt.savefig(best_save_path, dpi=150, bbox_inches='tight')
            print(f"最佳重建结果已保存至: {best_save_path}")

            plt.show()

    def run_full_protocol(self, save_results=True):
        """
        运行完整协议流程

        Parameters:
        -----------
        save_results : bool
            是否保存结果

        Returns:
        --------
        results : dict
            所有结果
        """
        print("="*70)
        print("瞬态磁场测量协议 - 完整流程")
        print("="*70)

        # 步骤1: 滑动测量
        sliding_results = self.run_sliding_measurement()

        # 步骤2: 提取核函数
        kernel_results = self.extract_kernel()

        # 步骤3: Wiener反卷积
        reconstruction_results = self.run_wiener_deconvolution()

        # 步骤4: 可视化
        if save_results:
            save_path = "transient_field_protocol_results.png"
        else:
            save_path = None

        self.visualize_full_protocol(save_path=save_path)

        # 汇总结果
        all_results = {
            'sliding_results': sliding_results,
            'kernel_results': kernel_results,
            'reconstruction_results': reconstruction_results,
            'protocol_params': self.protocol_params,
            'qubit_params': self.qubit_params
        }

        print("\n" + "="*70)
        print("协议完成摘要")
        print("="*70)
        print(f"量子比特频率: {self.qubit.frequency:.3f} GHz")
        print(f"磁场灵敏度: {self.protocol_params['sensitivity']/(2*np.pi):.3f} × 2π rad/ns/T")
        print(f"滑动测量点数: {len(sliding_results['delays'])}")
        print(f"核函数点数: {len(kernel_results['kernel'])}")

        # 最佳重建结果
        best_lambda = None
        best_correlation = -1
        for lam, errors in reconstruction_results['reconstruction_errors'].items():
            if errors['correlation'] > best_correlation:
                best_correlation = errors['correlation']
                best_lambda = lam

        if best_lambda is not None:
            print(f"最佳正则化参数: λ={best_lambda}")
            print(f"最佳相关系数: {best_correlation:.4f}")
            print(f"最佳RMSE: {reconstruction_results['reconstruction_errors'][best_lambda]['rmse']:.6f}")

        print("="*70)

        return all_results


def run_protocol_demo():
    """
    运行协议演示
    """
    # 自定义协议参数（可选）
    protocol_params = {
        'sensitivity': 2*np.pi*0.5,  # 较低灵敏度
        'pulse_duration': 15.0,      # 较短脉冲
        'wait_time': 25.0,           # 较长等待
        'delay_start': -30.0,        # 更宽延迟范围
        'delay_end': 90.0,
        'n_delays': 161,
        'noise_level': 0.015,        # 较低噪声
        'magnetic_signal_center': 50.0,
        'magnetic_signal_width': 6.0,  # 较窄信号
        'magnetic_signal_amplitude': 1.2,
        'wiener_lambdas': [0.001, 0.01, 0.1, 1.0, 10.0]  # 更多λ值
    }

    # 创建协议对象
    protocol = TransientFieldProtocol(protocol_params=protocol_params)

    # 运行完整协议
    results = protocol.run_full_protocol(save_results=True)

    return protocol, results


if __name__ == "__main__":
    run_protocol_demo()