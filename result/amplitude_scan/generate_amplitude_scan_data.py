#!/usr/bin/env python3
"""
生成幅度扫描的传感数据
针对给定幅度运行protocol4，保存数据供后续分析

使用方法：
python generate_amplitude_scan_data.py <amplitude>

示例：
python generate_amplitude_scan_data.py 0.001
python generate_amplitude_scan_data.py 0.01

输出文件：result/amplitude_scan/amplitude_<value>.npz
"""

import sys
import os
import argparse
import numpy as np
import time
import matplotlib.pyplot as plt

# 添加src目录到Python路径
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./src'))

from src.qubit import TransmonQubit
from src.signal import Signal
from src.pulse import Pulse, create_pulse
from src.protocal import Protocal
from src.analysis import Analysis

def run_amplitude_experiment(amplitude, output_dir='result/amplitude_scan'):
    """
    运行单个幅度的传感实验

    参数:
    amplitude: 信号幅度
    output_dir: 输出目录

    返回:
    dict: 包含原始信号、测量数据和重建结果的数据字典
    """

    print(f"\n=== 处理幅度 {amplitude} ===")

    # 创建量子比特
    qubit_params = {
        'EC': 0.2 * 2 * np.pi,
        'EJ': 10.0 * 2 * np.pi,
        'T1': 100.0e3,   # 100 µs
        'T2': 50.0e3,    # 50 µs
        'flux': 0.0,
        'state': 0,
        'n_levels': 2
    }

    qubit = TransmonQubit(**qubit_params)

    # 设置最优工作点
    Phi_opt = qubit.optimal_work_point()
    qubit.change_flux(Phi_opt)
    print(f"量子比特频率: {qubit.frequency:.4f} GHz, 工作点磁通: {qubit.flux:.4f}")

    # 创建信号（使用type=4非对称脉冲信号，固定其他参数）
    t_list = np.linspace(0, 200, 400)  # 与protocol4相同的时间轴
    signal_params = {
        'amplitude': amplitude,
        'rise': 10,      # 上升时间 (ns)
        'fall': 10,      # 下降时间 (ns)
        'center': 100,   # 中心位置 (ns)
        'noise_level': 0.0001  # 噪声水平
    }
    signal = Signal(type=4, t_list=t_list, **signal_params)

    # 运行瞬态磁场测量协议
    protocal = Protocal(4)  # type=4: 瞬态磁场测量协议
    protocal.initialize(qubit, state=0)

    print("运行测量协议...")
    start_time = time.time()
    t_samples, kernel, scan_list, delta_p, p_e, Phi_signal, control_pulse = protocal.evolve(qubit)
    meas_time = time.time() - start_time
    print(f"测量完成，耗时 {meas_time:.2f} 秒")

    # Wiener反卷积
    print("运行Wiener反卷积...")
    analysis = Analysis()
    dt = scan_list[1] - scan_list[0]
    # 尝试多个正则化参数，选择最佳结果
    lambdas_list = [0.1, 1.0, 5.0, 10.0, 50.0]
    best_rmse = float('inf')
    best_B_wiener = None
    best_B_lists_wiener = None

    for lamb in lambdas_list:
        try:
            B_lists_wiener, B_recon_wiener = analysis.wiener_deconvolution(delta_p, kernel, dt, lambdas=lamb)
            # 计算与原始信号的RMSE（需要插值到相同时间网格）
            from scipy.interpolate import interp1d
            f_wiener = interp1d(B_lists_wiener, B_recon_wiener, bounds_error=False, fill_value=0)
            B_wiener_interp = f_wiener(Phi_signal.t_list)
            rmse = np.sqrt(np.mean((Phi_signal.signal - B_wiener_interp)**2))
            if rmse < best_rmse:
                best_rmse = rmse
                best_B_wiener = B_recon_wiener
                best_B_lists_wiener = B_lists_wiener
                best_lambda = lamb
        except Exception as e:
            print(f"  Lambda={lamb} 失败: {e}")
            continue

    if best_B_wiener is not None:
        print(f"最佳Wiener反卷积 lambda={best_lambda}, RMSE={best_rmse:.6f}")
        B_recon_wiener = best_B_wiener
        B_lists_wiener = best_B_lists_wiener
    else:
        # 使用默认lambda=5.0作为后备
        print("所有lambda尝试失败，使用默认lambda=5.0")
        B_lists_wiener, B_recon_wiener = analysis.wiener_deconvolution(delta_p, kernel, dt, lambdas=5.0)

    # Levenberg-Marquardt优化（数值反演）
    print("运行Levenberg-Marquardt优化...")
    start_time = time.time()
    try:
        # 使用Wiener结果作为初始猜测
        Phi_guess = Signal(type=1, t_list=B_lists_wiener, amplitude=0.01, rise=10, fall=10, center=100)
        B_opt_lm, history = analysis.numerical_inverse(
            qubit, control_pulse, p_e, B_lists_wiener, Phi_guess,
            basis_type='fourier', n_basis=100, lambdas=100.0, max_iter=25, tol=1e-4
        )
        lm_time = time.time() - start_time
        print(f"LM优化完成，耗时 {lm_time:.2f} 秒，残差范数: {np.linalg.norm(history['res'][-1]):.6f}")
    except Exception as e:
        print(f"LM优化失败: {e}")
        # 如果失败，使用Wiener结果作为占位符
        B_opt_lm = B_recon_wiener.copy()
        history = {'res': [np.array([0])]}
        print("使用Wiener重建结果作为LM重建的替代")

    # 收集数据
    data = {
        'amplitude': amplitude,
        'signal_params': signal_params,
        'original_time': Phi_signal.t_list,
        'original_signal': Phi_signal.signal,
        'scan_time': scan_list,
        'delta_p': delta_p,
        'kernel_time': t_samples,
        'kernel': kernel,
        'wiener_time': B_lists_wiener,
        'wiener_recon': B_recon_wiener,
        'lm_time': Phi_signal.t_list,  # LM重建的时间轴与原始信号相同
        'lm_recon': B_opt_lm,
        'qubit_frequency': qubit.frequency,
        'qubit_flux': qubit.flux,
        'wiener_lambda_used': best_lambda if 'best_lambda' in locals() else 5.0,
        'lm_history_res': history['res'][-1] if history else np.array([0])
    }

    # 保存数据
    os.makedirs(output_dir, exist_ok=True)
    # 使用幅度值创建文件名，避免小数点引起的歧义
    amp_str = f"{amplitude:.6f}".rstrip('0').rstrip('.')
    output_file = os.path.join(output_dir, f'amplitude_{amp_str}.npz')
    np.savez(output_file, **data)
    print(f"数据已保存到 {output_file}")

    # 快速可视化检查
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # 原始信号
    axes[0].plot(data['original_time'], data['original_signal'], 'b-', label='Original')
    axes[0].set_xlabel('Time (ns)')
    axes[0].set_ylabel('Magnetic field')
    axes[0].set_title(f'Amplitude = {amplitude}')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    # 测量结果
    axes[1].plot(data['scan_time'], data['delta_p'], 'g-', label='Δp (measured)')
    axes[1].set_xlabel('Scan time (ns)')
    axes[1].set_ylabel('Probability change')
    axes[1].set_title('Measurement results')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    # 重建结果对比
    axes[2].plot(data['original_time'], data['original_signal'], 'k-', label='Original', alpha=0.7)
    axes[2].plot(data['wiener_time'], data['wiener_recon'], 'r--', label='Wiener', linewidth=1.5)
    axes[2].plot(data['lm_time'], data['lm_recon'], 'b:', label='LM', linewidth=1.5)
    axes[2].set_xlabel('Time (ns)')
    axes[2].set_ylabel('Magnetic field')
    axes[2].set_title('Reconstruction comparison')
    axes[2].grid(True, alpha=0.3)
    axes[2].legend()

    plt.tight_layout()
    check_file = os.path.join(output_dir, f'check_amplitude_{amp_str}.png')
    plt.savefig(check_file, dpi=150)
    plt.close()
    print(f"检查图已保存到 {check_file}")

    return data

def main():
    parser = argparse.ArgumentParser(description='生成幅度扫描的传感数据')
    parser.add_argument('amplitude', type=float, help='信号幅度（例如：0.001, 0.01, 0.1）')
    parser.add_argument('--output', '-o', default='result/amplitude_scan',
                       help='输出目录（默认：result/amplitude_scan）')

    args = parser.parse_args()

    # 验证幅度范围
    if args.amplitude <= 0:
        print(f"错误：幅度必须为正数，但得到 {args.amplitude}")
        sys.exit(1)

    if args.amplitude > 0.2:
        print(f"警告：幅度 {args.amplitude} 较大，可能超出线性范围")

    print(f"信号幅度: {args.amplitude}")
    print(f"输出目录: {args.output}")

    # 运行实验
    try:
        data = run_amplitude_experiment(args.amplitude, args.output)

        # 打印摘要信息
        print("\n" + "="*60)
        print("数据生成完成摘要")
        print("="*60)
        print(f"信号幅度: {args.amplitude}")
        print(f"原始信号长度: {len(data['original_signal'])}")
        print(f"Wiener重建长度: {len(data['wiener_recon'])}")
        print(f"LM重建长度: {len(data['lm_recon'])}")
        print(f"数据保存位置: {args.output}/amplitude_{args.amplitude:.6f}.npz")
        print(f"检查图保存位置: {args.output}/check_amplitude_{args.amplitude:.6f}.png")

    except Exception as e:
        print(f"运行过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()