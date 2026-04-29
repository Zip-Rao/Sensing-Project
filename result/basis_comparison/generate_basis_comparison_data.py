#!/usr/bin/env python3
"""
基函数对比数据生成脚本
运行B-spline、Fourier和Legendre基函数的数值反演，收集收敛历史和RMSE数据
符合科研制图规范的数据收集

使用方法：
python generate_basis_comparison_data.py [输出目录] [参数文件]

示例：
python generate_basis_comparison_data.py result/basis_comparison config.json
"""

import numpy as np
import sys
import os
import json
import time
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# 添加src目录到Python路径
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./src'))

try:
    from src.qubit import TransmonQubit
    from src.signal import Signal
    from src.pulse import Pulse, CompositePulse, create_ramsey_pulse
    from src.analysis import Analysis
    from src.protocal import Protocal
    print("成功导入所有模块")
except ImportError as e:
    print(f"导入模块时出错: {e}")
    print("请确保src目录中存在所有必要模块")
    sys.exit(1)

def set_default_parameters():
    """设置默认仿真参数"""
    params = {
        # 量子比特参数
        'qubit': {
            'EC': 0.2,          # 电容能量 (GHz)
            'EJ': 10.0,         # 约瑟夫森能量 (GHz)
            'T1': 100.0e3,      # 100 μs
            'T2': 50.0e3,       # 50 μs
            'flux': 0.0,        # 外加磁通
            'state': 0,         # 初始态 |0>
            'n_levels': 2       # 能级数
        },

        # 原始磁场信号参数
        'magnetic_signal': {
            'type': 4,          # 瞬态信号类型
            't_min': 0,         # 时间最小值 (ns)
            't_max': 200,       # 时间最大值 (ns)
            'n_points': 400,    # 时间点数
            'amplitude': 0.1,   # 信号幅度
            'rise': 10,         # 上升时间 (ns)
            'fall': 10,         # 下降时间 (ns)
            'center': 100,      # 中心时间 (ns)
            'noise_level': 0.0001  # 噪声水平
        },

        # 控制脉冲参数
        'control_pulse': {
            't_min': 0,         # 脉冲时间最小值 (ns)
            't_max': 20,        # 脉冲时间最大值 (ns)
            'n_points': 40,     # 脉冲时间点数
            'omega_d': None,    # 驱动频率，设为None时使用量子比特频率
            'phase1': 0.0,      # 第一个脉冲相位
            'phase2': 0.0       # 第二个脉冲相位
        },

        # 数值反演参数
        'inversion': {
            'basis_types': ['bspline', 'fourier', 'legendre'],  # 要测试的基函数类型
            'n_basis': 20,      # 基函数数量
            'lambdas': 100.0,   # 正则化参数
            'max_iter': 30,     # 最大迭代次数
            'tol': 1e-6,        # 收敛容忍度
            'B_guess_scale': 0.5  # 初始猜测信号的缩放因子
        },

        # 输出设置
        'output': {
            'data_dir': 'result/basis_comparison',
            'save_raw_data': True,      # 是否保存原始数据
            'save_plots': False,        # 是否保存中间图（数据生成时不保存）
            'verbose': True             # 是否输出详细信息
        }
    }
    return params

def create_qubit(params):
    """创建量子比特对象"""
    qubit_params = params['qubit']
    qubit = TransmonQubit(
        EC=qubit_params['EC'],
        EJ=qubit_params['EJ'],
        T1=qubit_params['T1'],
        T2=qubit_params['T2'],
        flux=qubit_params['flux'],
        state=qubit_params['state'],
        n_levels=qubit_params['n_levels']
    )
    return qubit

def create_magnetic_signal(params):
    """创建原始磁场信号"""
    sig_params = params['magnetic_signal']
    t_list = np.linspace(sig_params['t_min'], sig_params['t_max'], sig_params['n_points'])

    signal = Signal(
        type=sig_params['type'],
        t_list=t_list,
        amplitude=sig_params['amplitude'],
        rise=sig_params['rise'],
        fall=sig_params['fall'],
        center=sig_params['center'],
        noise_level=sig_params['noise_level']
    )
    return signal

def create_control_pulse(params, qubit):
    """创建控制脉冲"""
    pulse_params = params['control_pulse']
    t_pulse = np.linspace(pulse_params['t_min'], pulse_params['t_max'], pulse_params['n_points'])

    omega_d = pulse_params['omega_d']
    if omega_d is None:
        omega_d = qubit.frequency

    control_pulse = create_ramsey_pulse(
        t_pulse,
        tau=0.0,
        omega_d=omega_d,
        phase1=pulse_params['phase1'],
        phase2=pulse_params['phase2']
    )
    return control_pulse

def run_sliding_measurement(qubit, magnetic_signal, control_pulse):
    """运行滑动测量，获取测量数据"""
    print("运行滑动测量...")
    start_time = time.time()

    protocol = Protocal(type=4)
    protocol.initialize(qubit, state=0)

    # 运行瞬态磁场测量协议
    t_samples, kernel, scan_list, delta_p, p_e, Phi, control_pulse = protocol.evolve(qubit)

    elapsed = time.time() - start_time
    print(f"滑动测量完成，耗时 {elapsed:.2f} 秒")
    print(f"测量点数: {len(delta_p)}, 核函数点数: {len(kernel)}")

    return {
        't_samples': t_samples,
        'kernel': kernel,
        'scan_list': scan_list,
        'delta_p': delta_p,
        'p_e': p_e,
        'Phi': Phi,
        'control_pulse': control_pulse
    }

def calculate_rmse(original_time, original_signal, recon_time, recon_signal):
    """计算重建信号的RMSE"""
    # 插值重建信号到原始时间网格
    from scipy.interpolate import interp1d
    f_recon = interp1d(recon_time, recon_signal, bounds_error=False, fill_value=0)
    recon_interp = f_recon(original_time)

    # 计算RMSE
    rmse = np.sqrt(np.mean((original_signal - recon_interp)**2))
    return rmse

def run_numerical_inversion(qubit, measurement_data, basis_type, params):
    """运行数值反演算法"""
    print(f"运行 {basis_type} 基函数的数值反演...")
    start_time = time.time()

    inversion_params = params['inversion']
    analysis = Analysis()

    # 获取测量数据
    t_list = measurement_data['Phi'].t_list
    p_meas = measurement_data['delta_p']
    control_pulse = measurement_data['control_pulse']

    # 创建初始猜测信号（缩放后的原始信号）
    B_guess = measurement_data['Phi'].signal * inversion_params['B_guess_scale']

    # 运行数值反演
    B_opt, history = analysis.numerical_inverse(
        qubit=qubit,
        control_pulse=control_pulse,
        p_meas=p_meas,
        t_list=t_list,
        B_guess=B_guess,
        basis_type=basis_type,
        n_basis=inversion_params['n_basis'],
        lambdas=inversion_params['lambdas'],
        max_iter=inversion_params['max_iter'],
        tol=inversion_params['tol']
    )

    # 计算RMSE
    original_signal = measurement_data['Phi'].signal
    rmse = calculate_rmse(t_list, original_signal, t_list, B_opt)

    # 提取收敛历史
    convergence_data = {
        'iterations': list(range(len(history['res']))),
        'residual_norms': [np.linalg.norm(res) for res in history['res']],
        'residuals': [res.tolist() if isinstance(res, np.ndarray) else res for res in history['res']],
        'parameters': [b.tolist() if isinstance(b, np.ndarray) else b for b in history['b']],
        'mu_values': history.get('mu', [])
    }

    elapsed = time.time() - start_time
    print(f"{basis_type} 反演完成，耗时 {elapsed:.2f} 秒")
    print(f"  迭代次数: {len(convergence_data['iterations'])}, 最终RMSE: {rmse:.6f}")

    return {
        'basis_type': basis_type,
        'reconstructed_signal': B_opt.tolist(),
        'reconstructed_time': t_list.tolist(),
        'rmse': rmse,
        'convergence': convergence_data,
        'final_parameters': history['b'][-1].tolist() if len(history['b']) > 0 else [],
        'computation_time': elapsed
    }

def save_results(results, original_data, params, output_dir):
    """保存结果数据"""
    os.makedirs(output_dir, exist_ok=True)

    # 准备保存的数据
    save_data = {
        'metadata': {
            'generation_date': datetime.now().isoformat(),
            'parameters': params,
            'description': '基函数对比数据 - B-spline vs Fourier vs Legendre'
        },
        'original_signal': {
            'time': original_data['Phi'].t_list.tolist(),
            'signal': original_data['Phi'].signal.tolist(),
            'parameters': params['magnetic_signal']
        },
        'measurement_data': {
            'delta_p': original_data['delta_p'].tolist(),
            'scan_list': original_data['scan_list'].tolist(),
            't_samples': original_data['t_samples'].tolist(),
            'kernel': original_data['kernel'].tolist()
        },
        'basis_results': {}
    }

    # 添加每种基函数的结果
    for result in results:
        basis_type = result['basis_type']
        save_data['basis_results'][basis_type] = result

    # 保存为JSON文件
    json_file = os.path.join(output_dir, 'basis_comparison_data.json')
    with open(json_file, 'w', encoding='utf-8') as f:
        json.dump(save_data, f, indent=2, ensure_ascii=False)

    print(f"结果已保存到: {json_file}")

    # 同时保存为NPZ文件（便于数值计算）
    npz_data = {}
    npz_data['original_time'] = original_data['Phi'].t_list
    npz_data['original_signal'] = original_data['Phi'].signal
    npz_data['delta_p'] = original_data['delta_p']
    npz_data['scan_list'] = original_data['scan_list']

    for result in results:
        basis_type = result['basis_type']
        npz_data[f'{basis_type}_recon'] = np.array(result['reconstructed_signal'])
        npz_data[f'{basis_type}_rmse'] = result['rmse']
        npz_data[f'{basis_type}_residual_norms'] = np.array(result['convergence']['residual_norms'])

    npz_file = os.path.join(output_dir, 'basis_comparison_data.npz')
    np.savez_compressed(npz_file, **npz_data)
    print(f"数值数据已保存到: {npz_file}")

    return json_file, npz_file

def print_summary(results):
    """打印结果摘要"""
    print("\n" + "="*80)
    print("基函数对比结果摘要")
    print("="*80)
    print(f"{'基函数类型':<12} {'RMSE':<15} {'迭代次数':<12} {'计算时间(s)':<15}")
    print("-"*80)

    for result in results:
        basis_type = result['basis_type']
        rmse = result['rmse']
        n_iter = len(result['convergence']['iterations'])
        comp_time = result['computation_time']
        print(f"{basis_type:<12} {rmse:<15.6f} {n_iter:<12} {comp_time:<15.2f}")

    print("-"*80)

    # 找出最佳基函数
    best_result = min(results, key=lambda x: x['rmse'])
    print(f"最佳基函数: {best_result['basis_type']} (RMSE = {best_result['rmse']:.6f})")

def main():
    """主函数"""
    import argparse
    parser = argparse.ArgumentParser(description='生成基函数对比数据')
    parser.add_argument('--output-dir', '-o', default='result/basis_comparison',
                       help='输出目录（默认：result/basis_comparison）')
    parser.add_argument('--config', '-c', default=None,
                       help='配置文件路径（JSON格式）')
    parser.add_argument('--quick-test', action='store_true',
                       help='快速测试模式（减少参数）')

    args = parser.parse_args()

    # 加载参数
    params = set_default_parameters()

    if args.config and os.path.exists(args.config):
        with open(args.config, 'r') as f:
            user_params = json.load(f)
            # 深度更新参数
            import copy
            def deep_update(target, src):
                for k, v in src.items():
                    if isinstance(v, dict) and k in target and isinstance(target[k], dict):
                        deep_update(target[k], v)
                    else:
                        target[k] = v
            deep_update(params, user_params)
        print(f"从 {args.config} 加载配置")

    if args.quick_test:
        print("启用快速测试模式")
        params['magnetic_signal']['n_points'] = 200
        params['control_pulse']['n_points'] = 20
        params['inversion']['max_iter'] = 10
        params['inversion']['n_basis'] = 10

    # 更新输出目录
    params['output']['data_dir'] = args.output_dir

    print("="*80)
    print("基函数对比数据生成")
    print("="*80)
    print(f"输出目录: {args.output_dir}")
    print(f"基函数类型: {params['inversion']['basis_types']}")
    print(f"基函数数量: {params['inversion']['n_basis']}")
    print(f"最大迭代次数: {params['inversion']['max_iter']}")

    # 创建量子比特、信号和脉冲
    print("\n1. 创建仿真对象...")
    qubit = create_qubit(params)
    magnetic_signal = create_magnetic_signal(params)
    control_pulse = create_control_pulse(params, qubit)

    print(f"   量子比特频率: {qubit.frequency:.3f} GHz")
    print(f"   磁场信号范围: {magnetic_signal.t_list[0]:.1f} - {magnetic_signal.t_list[-1]:.1f} ns")
    print(f"   控制脉冲时长: {control_pulse.t_list[-1]:.1f} ns")

    # 运行滑动测量
    print("\n2. 生成测量数据...")
    measurement_data = run_sliding_measurement(qubit, magnetic_signal, control_pulse)

    # 运行不同基函数的数值反演
    print("\n3. 运行数值反演...")
    results = []
    basis_types = params['inversion']['basis_types']

    for basis_type in basis_types:
        try:
            result = run_numerical_inversion(qubit, measurement_data, basis_type, params)
            results.append(result)
        except Exception as e:
            print(f"  {basis_type} 基函数反演失败: {e}")
            import traceback
            traceback.print_exc()

    if not results:
        print("错误：所有基函数反演都失败")
        sys.exit(1)

    # 保存结果
    print("\n4. 保存结果数据...")
    json_file, npz_file = save_results(results, measurement_data, params, args.output_dir)

    # 打印摘要
    print_summary(results)

    print("\n" + "="*80)
    print("数据生成完成!")
    print(f"- JSON数据: {json_file}")
    print(f"- NPZ数据: {npz_file}")
    print(f"- 包含 {len(results)} 种基函数的结果")
    print("\n下一步：使用 plot_basis_comparison.py 绘制图表")
    print("="*80)

if __name__ == '__main__':
    main()