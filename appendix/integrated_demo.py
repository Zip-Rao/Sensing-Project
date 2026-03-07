"""
集成演示脚本：连接传感仿真与Wiener反卷积
演示完整的量子传感信号处理流程
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# 添加src目录到Python路径
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./src'))

# 导入自定义模块
try:
    from sensing_simulation import SensingSimulation, run_simulation as run_sensing_sim
    from deconvolution import WienerDeconvolution, run_deconvolution_demo
    print("模块导入成功")
except ImportError as e:
    print(f"模块导入失败: {e}")
    print("请确保sensing_simulation.py和deconvolution.py在同一目录下")
    sys.exit(1)


def run_integrated_demo():
    """
    运行集成演示：传感仿真 + Wiener反卷积
    """
    print("="*60)
    print("集成演示：量子传感信号处理完整流程")
    print("="*60)

    # ==========================================
    # 第一阶段：传感仿真获取真实核函数
    # ==========================================
    print("\n" + "-"*40)
    print("第一阶段：传感仿真获取真实核函数")
    print("-"*40)

    # 创建量子比特对象
    from src.qubit import TransmonQubit
    qubit = TransmonQubit(
        EC=0.2,      # 电容能量 (GHz)
        EJ=10.0,     # 约瑟夫森能量 (GHz)
        T1=100.0e3,  # 100 μs
        T2=50.0e3,   # 50 μs
        flux=0.0,    # 外加磁通
        state=0,     # 初始态 |0>
        n_levels=2   # 两能级系统
    )

    # 创建传感仿真对象
    sim = SensingSimulation(qubit, gamma=2*np.pi*1.0)

    # 时间轴设置
    t_total = 100  # ns
    steps = 500
    t_list = np.linspace(0, t_total, steps)

    # 生成理想控制脉冲并提取核函数
    pulse_center = 50.0  # ns
    pulse_width = 20.0   # ns

    print("生成理想控制脉冲...")
    ctrl_ideal = sim.generate_control_pulse(t_list, pulse_width, pulse_center,
                                           shape='gaussian', distortion=False)

    print("提取真实传感核函数...")
    t_scan_ideal, kernel_ideal, pe_base_ideal = sim.extract_kernel(
        t_list, ctrl_ideal, stim_strength=0.1, stim_width=1.0)

    # 由于提取的核函数是降采样后的，我们需要将其插值到完整时间轴
    from scipy.interpolate import interp1d
    interp_func = interp1d(t_scan_ideal, kernel_ideal, kind='cubic',
                           bounds_error=False, fill_value=0)
    kernel_full = interp_func(t_list)

    # 归一化核函数
    if np.max(np.abs(kernel_full)) > 0:
        kernel_full = kernel_full / np.max(np.abs(kernel_full))

    print(f"核函数提取完成，时间点: {len(t_scan_ideal)} -> 插值到: {len(t_list)}")
    print(f"核函数最大值: {np.max(kernel_full):.6f}")

    # ==========================================
    # 第二阶段：Wiener反卷积（使用真实核函数）
    # ==========================================
    print("\n" + "-"*40)
    print("第二阶段：Wiener反卷积（使用真实核函数）")
    print("-"*40)

    # 创建Wiener反卷积对象
    deconv = WienerDeconvolution()

    # 为反卷积准备时间轴（更高分辨率）
    dt = 0.2  # ns
    t_list_deconv = np.arange(0, t_total, dt)

    # 将核函数插值到反卷积时间轴
    interp_func_deconv = interp1d(t_list, kernel_full, kind='cubic',
                                  bounds_error=False, fill_value=0)
    kernel_deconv = interp_func_deconv(t_list_deconv)

    # 运行反卷积仿真
    results_real_kernel = deconv.run_simulation(
        use_real_kernel=True,
        kernel_data=(t_list_deconv, kernel_deconv)
    )

    # 可视化结果
    print("\n生成可视化结果（真实核函数）...")
    deconv.visualize_results(results_real_kernel,
                             save_path="wiener_deconv_real_kernel.png")

    # ==========================================
    # 第三阶段：对比模拟核函数与真实核函数
    # ==========================================
    print("\n" + "-"*40)
    print("第三阶段：模拟核函数与真实核函数对比")
    print("-"*40)

    # 运行模拟核函数的反卷积
    print("运行模拟核函数的反卷积...")
    results_mock_kernel = deconv.run_simulation(use_real_kernel=False)

    # 对比可视化
    compare_kernels_and_results(t_list, kernel_full, t_list_deconv,
                                kernel_deconv, results_real_kernel,
                                results_mock_kernel)

    print("\n" + "="*60)
    print("集成演示完成！")
    print("="*60)


def compare_kernels_and_results(t_list_sim, kernel_real, t_list_deconv,
                                kernel_deconv, results_real, results_mock):
    """
    对比模拟核函数与真实核函数的效果

    Parameters:
    -----------
    t_list_sim : array_like
        传感仿真的时间轴
    kernel_real : array_like
        真实核函数（传感仿真得到）
    t_list_deconv : array_like
        反卷积时间轴
    kernel_deconv : array_like
        插值到反卷积时间轴的核函数
    results_real : dict
        真实核函数的反卷积结果
    results_mock : dict
        模拟核函数的反卷积结果
    """
    # 创建对比图形
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 子图1: 核函数对比
    ax = axes[0, 0]
    ax.plot(t_list_sim, kernel_real / np.max(kernel_real),
            'b-', linewidth=2, label='Real Kernel (from Sensing)')
    ax.plot(t_list_deconv, kernel_deconv / np.max(kernel_deconv),
            'b--', linewidth=1.5, alpha=0.7, label='Interpolated for Deconv')

    # 生成模拟核函数用于对比
    deconv = WienerDeconvolution()
    kernel_mock = deconv.generate_mock_kernel(t_list_deconv,
                                              center=50, width=8.0)
    kernel_mock = kernel_mock / np.max(kernel_mock)
    ax.plot(t_list_deconv, kernel_mock, 'r-', linewidth=1.5,
            label='Mock Kernel (Gaussian)')

    ax.set_title('Kernel Comparison: Real vs Mock')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Normalized Amplitude')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 子图2: 重构误差对比
    ax = axes[0, 1]
    lambdas = results_real['lambdas']
    errors_real = results_real['reconstruction_errors']
    errors_mock = results_mock['reconstruction_errors']

    width = 0.35
    x = np.arange(len(lambdas))
    ax.bar(x - width/2, errors_real, width, label='Real Kernel', color='blue', alpha=0.7)
    ax.bar(x + width/2, errors_mock, width, label='Mock Kernel', color='red', alpha=0.7)

    ax.set_xlabel('Regularization Parameter λ')
    ax.set_ylabel('Reconstruction RMSE')
    ax.set_title('Reconstruction Error Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{lam:.4f}' for lam in lambdas])
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 子图3: 真实核函数的最佳重构结果（λ=0.1）
    ax = axes[1, 0]
    best_idx_real = np.argmin(errors_real)
    best_lam_real = lambdas[best_idx_real]
    B_rec_best_real = results_real['reconstructed_signals'][best_idx_real]

    ax.plot(results_real['t_list'], results_real['B_true'],
            'k-', linewidth=2, label='True Field $B(t)$', alpha=0.6)
    ax.plot(results_real['t_list'], B_rec_best_real,
            'b-', linewidth=1.5, label=f'Reconstructed (λ={best_lam_real:.4f})')

    ax.set_title(f'Best Reconstruction with Real Kernel\nλ={best_lam_real:.4f}, RMSE={errors_real[best_idx_real]:.6f}')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Field Amplitude')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 子图4: 模拟核函数的最佳重构结果
    ax = axes[1, 1]
    best_idx_mock = np.argmin(errors_mock)
    best_lam_mock = lambdas[best_idx_mock]
    B_rec_best_mock = results_mock['reconstructed_signals'][best_idx_mock]

    ax.plot(results_mock['t_list'], results_mock['B_true'],
            'k-', linewidth=2, label='True Field $B(t)$', alpha=0.6)
    ax.plot(results_mock['t_list'], B_rec_best_mock,
            'r-', linewidth=1.5, label=f'Reconstructed (λ={best_lam_mock:.4f})')

    ax.set_title(f'Best Reconstruction with Mock Kernel\nλ={best_lam_mock:.4f}, RMSE={errors_mock[best_idx_mock]:.6f}')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Field Amplitude')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # 保存对比图像
    save_path = "kernel_comparison_results.png"
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"对比图像已保存至: {save_path}")

    plt.show()

    # 打印对比摘要
    print("\n" + "="*50)
    print("核函数对比摘要")
    print("="*50)
    print(f"真实核函数形状参数:")
    print(f"  峰值位置: {t_list_sim[np.argmax(kernel_real)]:.1f} ns")
    print(f"  半高全宽: {calculate_fwhm(t_list_sim, kernel_real):.2f} ns")
    print(f"  面积: {np.trapz(kernel_real, t_list_sim):.4f}")

    print(f"\n模拟核函数形状参数:")
    print(f"  峰值位置: {t_list_deconv[np.argmax(kernel_mock)]:.1f} ns")
    print(f"  半高全宽: {calculate_fwhm(t_list_deconv, kernel_mock):.2f} ns")
    print(f"  面积: {np.trapz(kernel_mock, t_list_deconv):.4f}")

    print(f"\n最佳重构误差对比:")
    print(f"  真实核函数: λ={best_lam_real:.4f}, RMSE={errors_real[best_idx_real]:.6f}")
    print(f"  模拟核函数: λ={best_lam_mock:.4f}, RMSE={errors_mock[best_idx_mock]:.6f}")
    print(f"  误差比 (模拟/真实): {errors_mock[best_idx_mock]/errors_real[best_idx_real]:.3f}")
    print("="*50)


def calculate_fwhm(t, signal):
    """
    计算信号的半高全宽 (FWHM)

    Parameters:
    -----------
    t : array_like
        时间轴
    signal : array_like
        信号值

    Returns:
    --------
    fwhm : float
        半高全宽，如果无法计算则返回0
    """
    signal = np.abs(signal)
    max_val = np.max(signal)
    half_max = max_val / 2.0

    if max_val == 0:
        return 0.0

    # 找到超过半高全宽的点
    above_half = signal >= half_max

    if np.sum(above_half) == 0:
        return 0.0

    # 找到第一个和最后一个超过半高全宽的点
    indices = np.where(above_half)[0]
    if len(indices) == 0:
        return 0.0

    first_idx = indices[0]
    last_idx = indices[-1]

    if first_idx >= len(t) - 1 or last_idx >= len(t) - 1:
        return 0.0

    # 线性插值获取精确的半高全宽点
    if first_idx > 0:
        t1 = t[first_idx - 1] + (t[first_idx] - t[first_idx - 1]) * \
             (half_max - signal[first_idx - 1]) / (signal[first_idx] - signal[first_idx - 1])
    else:
        t1 = t[first_idx]

    if last_idx < len(t) - 1:
        t2 = t[last_idx] + (t[last_idx + 1] - t[last_idx]) * \
             (half_max - signal[last_idx]) / (signal[last_idx + 1] - signal[last_idx])
    else:
        t2 = t[last_idx]

    return t2 - t1


if __name__ == "__main__":
    run_integrated_demo()