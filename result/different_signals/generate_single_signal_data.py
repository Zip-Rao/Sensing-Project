#!/usr/bin/env python3
"""
生成单个信号类型的传感数据
运行Wiener反卷积和Levenberg-Marquardt优化算法
保存数据供绘图使用

使用方法：
python generate_single_signal_data.py <signal_type>
其中<signal_type>可以是：
  2 或 sine      : 正弦信号
  4 或 step      : 阶跃信号（非对称脉冲模拟）
  5 或 double    : 双峰高斯信号
  7 或 complex   : 复杂多峰信号

示例：
python generate_single_signal_data.py sine
python generate_single_signal_data.py 2
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

def parse_signal_type(signal_str):
    """将字符串或数字转换为信号类型整数和名称"""
    signal_map = {
        '2': (2, 'sine'),
        'sine': (2, 'sine'),
        '4': (4, 'step'),
        'step': (4, 'step'),
        '5': (5, 'double_peak'),
        'double': (5, 'double_peak'),
        'double_peak': (5, 'double_peak'),
        '7': (7, 'complex'),
        'complex': (7, 'complex')
    }

    if signal_str.lower() in signal_map:
        return signal_map[signal_str.lower()]
    else:
        try:
            signal_int = int(signal_str)
            for key, (sig_type, name) in signal_map.items():
                if sig_type == signal_int:
                    return (sig_type, name)
        except ValueError:
            pass
        raise ValueError(f"未知的信号类型: {signal_str}。可用选项: 2/sine, 4/step, 5/double, 7/complex")

def get_signal_params(signal_type):
    """根据信号类型返回默认参数"""
    if signal_type == 2:  # 正弦信号
        return {'amplitude': 0.01, 'frequency': 0.01}
    elif signal_type == 4:  # 阶跃信号（非对称脉冲）
        return {'amplitude': 0.01, 'center': 100, 'rise': 10, 'fall': 10}
    elif signal_type == 5:  # 双峰高斯信号
        return {'amplitude': 0.01, 'center': 100, 'width': 40}
    elif signal_type == 7:  # 复杂多峰信号
        return {'amplitude': 0.01}
    else:
        return {'amplitude': 0.01}

def run_signal_reconstruction(signal_type, signal_params, output_dir='result/different_signals'):
    """
    运行单个信号类型的传感和重建

    参数:
    signal_type: 信号类型（整数，对应Signal类的type参数）
    signal_params: 信号参数字典
    output_dir: 输出目录

    返回:
    dict: 包含原始信号、测量数据和重建结果的数据字典
    """

    print(f"\n=== 处理信号类型 {signal_type} ===")
    print(f"信号参数: {signal_params}")

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

    # 创建信号
    t_list = np.linspace(0, 200, 400)  # 与test5相同的时间轴
    signal = Signal(type=signal_type, t_list=t_list, **signal_params)

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
    B_lists_wiener, B_recon_wiener = analysis.wiener_deconvolution(delta_p, kernel, dt, lambdas=5.0)

    # Levenberg-Marquardt优化（数值反演）
    print("运行Levenberg-Marquardt优化...")
    start_time = time.time()
    try:
        # 使用Wiener结果作为初始猜测
        Phi = Signal(type = 1, t_list = B_lists_wiener, amplitude = 0.01, rise = 10, fall = 10, center = 100)
        B_opt_lm, history = analysis.numerical_inverse(
            qubit, control_pulse, p_e, B_lists_wiener, Phi,
            basis_type='fourier', n_basis=100, lambdas=100.0, max_iter=20, tol=1e-4
        )
        lm_time = time.time() - start_time
        print(f"LM优化完成，耗时 {lm_time:.2f} 秒，残差范数: {np.linalg.norm(history['res'][-1]):.6f}")
    except Exception as e:
        print(f"LM优化失败: {e}")
        # 如果失败，使用Wiener结果作为占位符
        B_opt_lm = B_recon_wiener.copy()
        print("使用Wiener重建结果作为LM重建的替代")

    # 收集数据
    data = {
        'signal_type': signal_type,
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
    }

    # 保存数据
    os.makedirs(output_dir, exist_ok=True)
    signal_name = get_signal_name(signal_type)
    output_file = os.path.join(output_dir, f'signal_{signal_name}.npz')
    np.savez(output_file, **data)
    print(f"数据已保存到 {output_file}")

    # 快速可视化检查
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # 原始信号
    axes[0].plot(data['original_time'], data['original_signal'], 'b-', label='Original')
    axes[0].set_xlabel('Time (ns)')
    axes[0].set_ylabel('Magnetic field')
    axes[0].set_title(f'{signal_name} Signal')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    # 重建结果对比
    axes[1].plot(data['original_time'], data['original_signal'], 'k-', label='Original', alpha=0.7)
    axes[1].plot(data['wiener_time'], data['wiener_recon'], 'r--', label='Wiener', linewidth=1.5)
    axes[1].plot(data['lm_time'], data['lm_recon'], 'b:', label='LM', linewidth=1.5)
    axes[1].set_xlabel('Time (ns)')
    axes[1].set_ylabel('Magnetic field')
    axes[1].set_title(f'{signal_name} Reconstruction')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    plt.tight_layout()
    check_file = os.path.join(output_dir, f'check_{signal_name}.png')
    plt.savefig(check_file, dpi=150)
    plt.close()
    print(f"检查图已保存到 {check_file}")

    return data

def get_signal_name(signal_type):
    """将信号类型转换为可读名称"""
    signal_names = {
        2: 'sine',
        4: 'step',
        5: 'double_peak',
        7: 'complex'
    }
    return signal_names.get(signal_type, f'type{signal_type}')

def main():
    parser = argparse.ArgumentParser(description='生成单个信号类型的传感数据')
    parser.add_argument('signal', help='信号类型：2/sine, 4/step, 5/double, 7/complex')
    parser.add_argument('--output', '-o', default='result/different_signals',
                       help='输出目录（默认：result/different_signals）')
    parser.add_argument('--amplitude', '-a', type=float, default=None,
                       help='信号幅度（覆盖默认值）')
    parser.add_argument('--center', '-c', type=float, default=None,
                       help='信号中心位置（ns，仅对step和double信号有效）')
    parser.add_argument('--frequency', '-f', type=float, default=None,
                       help='信号频率（GHz，仅对sine信号有效）')

    args = parser.parse_args()

    # 解析信号类型
    try:
        signal_type, signal_name = parse_signal_type(args.signal)
    except ValueError as e:
        print(f"错误: {e}")
        sys.exit(1)

    print(f"信号类型: {signal_type} ({signal_name})")

    # 获取默认参数
    signal_params = get_signal_params(signal_type)

    # 覆盖用户提供的参数
    if args.amplitude is not None:
        signal_params['amplitude'] = args.amplitude
    if args.center is not None:
        signal_params['center'] = args.center
    if args.frequency is not None and signal_type == 2:
        signal_params['frequency'] = args.frequency

    # 运行重建
    try:
        data = run_signal_reconstruction(signal_type, signal_params, args.output)

        # 打印摘要信息
        print("\n" + "="*60)
        print("数据生成完成摘要")
        print("="*60)
        print(f"信号类型: {signal_name}")
        print(f"原始信号长度: {len(data['original_signal'])}")
        print(f"Wiener重建长度: {len(data['wiener_recon'])}")
        print(f"LM重建长度: {len(data['lm_recon'])}")
        print(f"数据保存位置: {args.output}/signal_{signal_name}.npz")
        print(f"检查图保存位置: {args.output}/check_{signal_name}.png")

    except Exception as e:
        print(f"运行过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()