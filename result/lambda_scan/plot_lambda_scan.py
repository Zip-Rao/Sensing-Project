#!/usr/bin/env python3
"""
绘制Wiener反卷积和LM优化在不同正则化参数lambda下的磁场信号重建结果
使用protocol4（瞬态磁场测量协议）获取测量数据p_meas
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import time

# 添加src目录到Python路径
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./src'))

from src.qubit import TransmonQubit
from src.signal import Signal
from src.protocal import Protocal
from src.analysis import Analysis

def run_protocol4():
    """
    运行protocol4（瞬态磁场测量协议）获取测量数据
    返回：原始信号、测量数据、核函数等
    """
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

    # 创建瞬态磁场信号（与protocol4中相同）
    t_list = np.linspace(0, 200, 400)
    Phi = Signal(type=4, t_list=t_list, amplitude=0.08, rise=10, fall=10, center=100, noise_level=0.0001)

    # 运行瞬态磁场测量协议
    protocol = Protocal(4)  # type=4: 瞬态磁场测量协议
    protocol.initialize(qubit, state=0)

    print("运行测量协议...")
    start_time = time.time()
    t_samples, kernel, scan_list, delta_p, p_e, Phi_signal, control_pulse = protocol.evolve(qubit)
    meas_time = time.time() - start_time
    print(f"测量完成，耗时 {meas_time:.2f} 秒")

    return {
        'qubit': qubit,
        'original_signal': Phi_signal,
        'control_pulse': control_pulse,
        't_samples': t_samples,
        'kernel': kernel,
        'scan_list': scan_list,
        'delta_p': delta_p,
        'p_e': p_e,
        'dt': scan_list[1] - scan_list[0]
    }

def wiener_lambda_scan(data, lambdas):
    """
    对一组lambda值运行Wiener反卷积
    返回：每个lambda的重建结果
    """
    delta_p = data['delta_p']
    kernel = data['kernel']
    dt = data['dt']

    wiener_results = []
    analysis = Analysis()

    for lam in lambdas:
        print(f"运行Wiener反卷积，lambda={lam}...")
        try:
            B_lists, B_recon = analysis.wiener_deconvolution(delta_p, kernel, dt, lambdas=lam)
            wiener_results.append({
                'lambda': lam,
                'time': B_lists,
                'reconstruction': B_recon
            })
        except Exception as e:
            print(f"  lambda={lam} 失败: {e}")
            wiener_results.append({
                'lambda': lam,
                'time': None,
                'reconstruction': None
            })

    return wiener_results

def lm_lambda_scan(data, lambdas, basis_type='fourier', n_basis=100, max_iter=20, tol=1e-4):
    """
    对一组lambda值运行LM优化
    注意：LM优化计算成本高，建议使用较少的lambda值
    """
    qubit = data['qubit']
    control_pulse = data['control_pulse']
    # 使用p_e作为测量结果，与generate_single_signal_data.py保持一致
    p_meas = data['p_e']

    # 获取Wiener反卷积结果作为初始猜测
    analysis = Analysis()
    delta_p = data['delta_p']
    kernel = data['kernel']
    dt = data['dt']
    # 使用中等lambda值进行Wiener反卷积，获取初始猜测
    B_lists_wiener, B_recon_wiener = analysis.wiener_deconvolution(delta_p, kernel, dt, lambdas=5.0)
    # Wiener反卷积返回的时间轴就是磁场信号的时间轴
    t_list = B_lists_wiener

    lm_results = []

    for lam in lambdas:
        print(f"运行LM优化，lambda={lam}...")
        start_time = time.time()
        try:
            # 使用Wiener重建结果作为初始猜测（数组）
            B_guess = Signal(type=1, t_list=B_lists_wiener, amplitude=0.01, rise=10, fall=10, center=100) # 生成与Wiener时间轴相同的初始猜测
            B_opt, history = analysis.numerical_inverse(
                qubit, control_pulse, p_meas, t_list, B_guess,  # 直接使用Wiener重建结果的数值作为初始猜测
                basis_type=basis_type, n_basis=n_basis, lambdas=lam,
                max_iter=max_iter, tol=tol
            )
            lm_time = time.time() - start_time
            print(f"  LM优化完成，耗时 {lm_time:.2f} 秒，残差范数: {np.linalg.norm(history['res'][-1]):.6f}")
            lm_results.append({
                'lambda': lam,
                'reconstruction': B_opt,
                'history': history,
                'time': t_list  # 重建信号的时间轴与Wiener反卷积相同
            })
        except Exception as e:
            print(f"  lambda={lam} 失败: {e}")
            # 如果失败，尝试使用零数组作为初始猜测
            try:
                B_guess_zero = np.zeros_like(B_recon_wiener)
                B_opt, history = analysis.numerical_inverse(
                    qubit, control_pulse, p_meas, t_list, B_guess_zero,
                    basis_type=basis_type, n_basis=n_basis, lambdas=lam,
                    max_iter=max_iter, tol=tol
                )
                lm_time = time.time() - start_time
                print(f"  LM优化完成（使用零初始猜测），耗时 {lm_time:.2f} 秒，残差范数: {np.linalg.norm(history['res'][-1]):.6f}")
                lm_results.append({
                    'lambda': lam,
                    'reconstruction': B_opt,
                    'history': history,
                    'time': t_list
                })
            except Exception as e2:
                print(f"  使用零初始猜测也失败: {e2}")
                lm_results.append({
                    'lambda': lam,
                    'reconstruction': None,
                    'history': None,
                    'time': None
                })

    return lm_results

def plot_wiener_results(original_signal, wiener_results, save_path=None):
    """
    绘制Wiener反卷积结果
    """
    plt.figure(figsize=(10, 6))

    # 绘制原始信号
    plt.plot(original_signal.t_list, original_signal.signal, 'k--', linewidth=1, label='Original signal')

    # 绘制每个lambda的重建结果
    colors = plt.cm.viridis(np.linspace(0, 1, len(wiener_results)))
    for i, res in enumerate(wiener_results):
        if res['reconstruction'] is not None:
            plt.plot(res['time'], res['reconstruction'], '-', color=colors[i],
                    linewidth=1.5, alpha=0.8, label=f'λ={res["lambda"]:.2f}')

    plt.xlabel('Time (ns)', fontsize=12)
    plt.ylabel('Magnetic field (a.u.)', fontsize=12)
    plt.title('Wiener Deconvolution Reconstruction with Different λ', fontsize=14)
    plt.legend(fontsize=10, loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Wiener反卷积图已保存到 {save_path}")
    else:
        plt.show()

def plot_lm_results(original_signal, lm_results, save_path=None):
    """
    绘制LM优化结果
    """
    plt.figure(figsize=(10, 6))

    # 绘制原始信号
    plt.plot(original_signal.t_list, original_signal.signal, 'k--', linewidth=1, label='Original signal')

    # 绘制每个lambda的重建结果
    colors = plt.cm.plasma(np.linspace(0, 1, len(lm_results)))
    for i, res in enumerate(lm_results):
        if res['reconstruction'] is not None:
            plt.plot(res['time'], res['reconstruction'], '-', color=colors[i],
                    linewidth=1.5, alpha=0.8, label=f'λ={res["lambda"]:.2f}')

    plt.xlabel('Time (ns)', fontsize=12)
    plt.ylabel('Magnetic field (a.u.)', fontsize=12)
    plt.title('LM Optimization Reconstruction with Different λ', fontsize=14)
    plt.legend(fontsize=10, loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"LM优化图已保存到 {save_path}")
    else:
        plt.show()

def main():
    """主函数"""
    print("=" * 60)
    print("绘制Wiener反卷积和LM优化在不同λ下的重建结果")
    print("=" * 60)

    # 运行protocol4获取测量数据
    data = run_protocol4()

    # 定义lambda值列表（Wiener反卷积使用较多值，LM优化使用较少值）
    wiener_lambdas = [2.0, 3.0, 4.0, 5.0, 10.0, 50.0, 100.0]
    lm_lambdas = [0.1, 1.0, 10.0, 100.0]  # LM优化计算成本高，减少数量

    # 运行Wiener反卷积扫描
    print("\n--- Wiener反卷积λ扫描 ---")
    wiener_results = wiener_lambda_scan(data, wiener_lambdas)

    # 运行LM优化扫描
    print("\n--- LM优化λ扫描 ---")
    lm_results = lm_lambda_scan(data, lm_lambdas, max_iter=10)  # 减少迭代次数以加快速度

    # 创建输出目录
    output_dir = "result/lambda_scan"
    os.makedirs(output_dir, exist_ok=True)

    # 绘制结果
    wiener_plot_path = os.path.join(output_dir, "wiener_lambda_scan.png")
    lm_plot_path = os.path.join(output_dir, "lm_lambda_scan.png")

    plot_wiener_results(data['original_signal'], wiener_results, save_path=wiener_plot_path)
    plot_lm_results(data['original_signal'], lm_results, save_path=lm_plot_path)

    # 保存数据供后续分析
    data_path = os.path.join(output_dir, "lambda_scan_data.npz")

    # 提取Wiener重建结果（只保存成功的）
    wiener_success = [res for res in wiener_results if res['reconstruction'] is not None]
    n_wiener = len(wiener_success)
    if n_wiener > 0:
        # 假设所有成功结果的时间轴相同
        wiener_time = wiener_success[0]['time']
        wiener_lambdas_arr = np.array([res['lambda'] for res in wiener_success])
        wiener_recon_arr = np.array([res['reconstruction'] for res in wiener_success])
    else:
        wiener_time = np.array([])
        wiener_lambdas_arr = np.array([])
        wiener_recon_arr = np.array([])

    # 提取LM重建结果
    lm_success = [res for res in lm_results if res['reconstruction'] is not None]
    n_lm = len(lm_success)
    if n_lm > 0:
        lm_time = lm_success[0]['time']
        lm_lambdas_arr = np.array([res['lambda'] for res in lm_success])
        lm_recon_arr = np.array([res['reconstruction'] for res in lm_success])
    else:
        lm_time = np.array([])
        lm_lambdas_arr = np.array([])
        lm_recon_arr = np.array([])

    np.savez(data_path,
             original_time=data['original_signal'].t_list,
             original_signal=data['original_signal'].signal,
             wiener_time=wiener_time,
             wiener_lambdas=wiener_lambdas_arr,
             wiener_reconstructions=wiener_recon_arr,
             lm_time=lm_time,
             lm_lambdas=lm_lambdas_arr,
             lm_reconstructions=lm_recon_arr,
             )
    print(f"\n数据已保存到 {data_path}")

    print("\n" + "=" * 60)
    print("完成！")
    print(f"Wiener反卷积图: {wiener_plot_path}")
    print(f"LM优化图: {lm_plot_path}")
    print("=" * 60)

if __name__ == '__main__':
    main()