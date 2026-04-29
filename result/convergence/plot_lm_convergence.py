#!/usr/bin/env python3
"""
绘制Levenberg-Marquardt算法残差随迭代次数的变化图，展示收敛性。
首先运行protocol4获取数据，然后运行LM算法，保存所有相关参数以便后续绘图。
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pickle

# 添加src目录到Python路径
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./src'))

# 导入自定义模块
from src.qubit import TransmonQubit
from src.protocal import Protocal
from src.analysis import Analysis
from src.signal import Signal

# 设置科研制图规范
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'mathtext.fontset': 'stix',  # 使用STIX字体，类似于Times New Roman
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 16,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'figure.autolayout': False,
})

def run_protocol4():
    """
    运行protocol4（瞬态磁场测量协议）获取数据。
    返回qubit, control_pulse, delta_p, Phi, scan_list等。
    """
    print("="*60)
    print("步骤1: 运行protocol4获取数据")
    print("="*60)

    # 创建量子比特对象（使用默认参数）
    qubit = TransmonQubit(
        EC=0.2 * 2 * np.pi,      # 电容能量 (GHz)
        EJ=10.0 * 2 * np.pi,     # 约瑟夫森能量 (GHz)
        T1=100.0e3,  # 100 μs
        T2=50.0e3,   # 50 μs
        flux=0.0,    # 外加磁通
        state=0,     # 初始态 |0>
        n_levels=2   # 两能级系统
    )
    Phi_opt = qubit.optimal_work_point()
    qubit.change_flux(Phi_opt)
    print(f"创建量子比特: 频率 = {qubit.frequency:.3f} GHz, 工作点磁通 = {qubit.flux:.4f}")

    # 创建protocol对象，类型为4（瞬态磁场测量协议）
    protocol = Protocal(type=4)

    # 初始化protocol（设置量子比特初始态）
    protocol.initialize(qubit, state=0)

    # 运行protocol4的evolve方法获取数据
    # 注意: evolve方法返回多个值，具体参见protocal.py中的case 4
    print("运行protocol4演化...")
    results = protocol.evolve(qubit)

    # 解析返回结果
    # 根据protocal.py第152行，返回值为:
    # t_samples, kernel, scan_list, delta_p, p_e, Phi, control_pulse
    t_samples, kernel, scan_list, delta_p, p_e, Phi, control_pulse = results

    print(f"protocol4完成:")
    print(f"  时间采样点: {len(t_samples)}")
    print(f"  核函数长度: {len(kernel)}")
    print(f"  扫描列表长度: {len(scan_list)}")
    print(f"  概率变化长度: {len(delta_p)}")
    print(f"  磁场信号长度: {len(Phi.t_list)}")
    print(f"  控制脉冲时间长度: {len(control_pulse.t_list)}")

    return qubit, control_pulse, delta_p, Phi, scan_list, t_samples, kernel, p_e

def run_lm_algorithm(qubit, control_pulse, p_e, Phi, scan_list,
                     basis_type='fourier', n_basis=20, reg=100.0,
                     max_iter=30, tol=1e-6):
    """
    运行Levenberg-Marquardt算法进行数值反演。

    参数:
        qubit: 量子比特对象
        control_pulse: 控制脉冲对象
        delta_p: 测量的概率变化
        Phi: 磁场信号对象（真实信号，用于获取时间轴）
        scan_list: 滑动测量时间轴
        basis_type: 基函数类型 ('fourier', 'bspline', 'legendre')
        n_basis: 基函数数量
        reg: 正则化参数
        max_iter: 最大迭代次数
        tol: 收敛容忍度

    返回:
        B_opt: 优化后的磁场信号
        history: 包含残差、参数等历史的字典
        results: 包含所有相关数据的字典
    """
    print("\n" + "="*60)
    print("步骤2: 运行Levenberg-Marquardt算法")
    print("="*60)

    # 获取磁场信号的时间轴
    t_list = Phi.t_list

    # 初始猜测：零磁场信号（或使用真实信号的平滑版本）
    B_guess = np.zeros_like(t_list)

    # 调用numerical_inverse函数（该函数内部调用levenberg_marquardt）
    print(f"使用基函数: {basis_type}, 数量: {n_basis}")
    print(f"正则化参数: {reg}, 最大迭代次数: {max_iter}")

    analysis = Analysis()
    B_opt, history = analysis.numerical_inverse(
        qubit=qubit,
        control_pulse=control_pulse,
        p_meas=p_e,
        t_list=t_list,
        B_guess=B_guess,
        basis_type=basis_type,
        n_basis=n_basis,
        lambdas=reg,
        max_iter=max_iter,
        tol=tol
    )

    print(f"LM算法完成，迭代次数: {len(history['res'])}")
    print(f"最终残差范数: {np.linalg.norm(history['res'][-1]):.6e}")

    # 收集所有结果
    results = {
        'B_opt': B_opt,
        'history': history,
        't_list': t_list,
        'delta_p': p_e,
        'scan_list': scan_list,
        'basis_type': basis_type,
        'n_basis': n_basis,
        'reg': reg,
        'max_iter': max_iter,
        'tol': tol,
        'qubit_params': {
            'EC': qubit.EC,
            'EJ': qubit.EJ,
            'frequency': qubit.frequency,
            'anharmonicity': qubit.anharmonicity,
            'T1': qubit.T1,
            'T2': qubit.T2,
        }
    }

    return B_opt, history, results

def plot_residual_convergence(history, save_path='lm_convergence.png'):
    """
    绘制残差随迭代次数的变化图。

    参数:
        history: LM算法的历史记录字典
        save_path: 图像保存路径
    """
    print("\n" + "="*60)
    print("步骤3: 绘制残差收敛图")
    print("="*60)

    # 提取残差范数
    residuals = history['res']
    iter_indices = np.arange(1, len(residuals) + 1)
    res_norms = [np.linalg.norm(r) for r in residuals]

    # 提取阻尼参数mu（如果存在）
    mu_history = history.get('mu', [])

    # 创建图形
    fig, axes = plt.subplots(1, 1, figsize=(14, 5))

    # 子图1: 残差范数随迭代次数的变化（线性坐标）
    ax1 = axes
    ax1.plot(iter_indices, res_norms, 'bo-', linewidth=2, markersize=6,
             markerfacecolor='white', markeredgewidth=1.5)
    ax1.set_xlabel('Iteration Number', fontsize=14)
    ax1.set_ylabel(r'Residual Norm $\|r\|$', fontsize=14)
    ax1.set_title('LM Algorithm Convergence', fontsize=16)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_yscale('log')  # 使用对数坐标更好显示收敛
    ax1.set_xlim(0, len(iter_indices) + 1)

    # 标记最终残差
    final_res = res_norms[-1]
    ax1.annotate(f'Final: {final_res:.2e}',
                xy=(iter_indices[-1], final_res),
                xytext=(iter_indices[-1] - 2, final_res * 10),
                arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
                fontsize=12,
                color='red')

    # 子图2: 阻尼参数mu的变化（如果存在）
    # ax2 = axes[1]
    # if len(mu_history) > 0:
    #     ax2.plot(iter_indices[:len(mu_history)], mu_history, 'rs-',
    #              linewidth=2, markersize=6, markerfacecolor='white',
    #              markeredgewidth=1.5)
    #     ax2.set_xlabel('Iteration Number', fontsize=14)
    #     ax2.set_ylabel('Damping Parameter $\mu$', fontsize=14)
    #     ax2.set_title('Damping Parameter Evolution', fontsize=16)
    #     ax2.grid(True, alpha=0.3, linestyle='--')
    #     ax2.set_yscale('log')
    #     ax2.set_xlim(0, len(iter_indices) + 1)
    # else:
    #     ax2.text(0.5, 0.5, 'Damping parameter data\nnot available',
    #             horizontalalignment='center',
    #             verticalalignment='center',
    #             transform=ax2.transAxes,
    #             fontsize=14)
    #     ax2.set_title('Damping Parameter Evolution', fontsize=16)

    plt.tight_layout()

    # 保存图像
    plt.savefig(save_path, dpi=300)
    print(f"收敛图已保存至: {save_path}")

    # 显示图像
    plt.show()

    return fig

def save_all_parameters(results, save_path='lm_parameters.pkl'):
    """
    保存所有相关参数以便后续使用。

    参数:
        results: 包含所有结果的字典
        save_path: 参数保存路径
    """
    print("\n" + "="*60)
    print("步骤4: 保存所有参数")
    print("="*60)

    # 确保目录存在
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)

    # 保存为pickle文件
    with open(save_path, 'wb') as f:
        pickle.dump(results, f)

    print(f"所有参数已保存至: {save_path}")
    print(f"保存的内容包括:")
    print(f"  - 优化后的磁场信号 B_opt (长度: {len(results['B_opt'])})")
    print(f"  - 历史记录 (迭代次数: {len(results['history']['res'])})")
    print(f"  - 基函数设置: {results['basis_type']}, {results['n_basis']}")
    print(f"  - 正则化参数: {results['reg']}")
    print(f"  - 量子比特参数")

    # 同时保存为文本文件（部分信息）
    txt_path = save_path.replace('.pkl', '.txt')
    with open(txt_path, 'w') as f:
        f.write("LM算法参数摘要\n")
        f.write("="*50 + "\n")
        f.write(f"迭代次数: {len(results['history']['res'])}\n")
        f.write(f"最终残差范数: {np.linalg.norm(results['history']['res'][-1]):.6e}\n")
        f.write(f"基函数类型: {results['basis_type']}\n")
        f.write(f"基函数数量: {results['n_basis']}\n")
        f.write(f"正则化参数: {results['reg']}\n")
        f.write(f"量子比特频率: {results['qubit_params']['frequency']:.3f} GHz\n")
        f.write(f"量子比特非谐性: {results['qubit_params']['anharmonicity']:.3f} GHz\n")

    print(f"参数摘要已保存至: {txt_path}")

def main():
    """主函数"""
    print("="*70)
    print("LM算法收敛性分析")
    print("="*70)

    try:
        # 步骤1: 运行protocol4获取数据
        qubit, control_pulse, delta_p, Phi, scan_list, t_samples, kernel, p_e = run_protocol4()

        # 步骤2: 运行LM算法
        B_opt, history, results = run_lm_algorithm(
            qubit=qubit,
            control_pulse=control_pulse,
            p_e=p_e,
            Phi=Phi,
            scan_list=scan_list,
            basis_type='fourier',  # 可以更改为'bspline'或'legendre'
            n_basis=100,
            reg=100.0,
            max_iter=100,
            tol=1e-6
        )
        save_all_parameters(results, save_path='lm_parameters.pkl')
        # 步骤3: 绘制残差收敛图
        fig = plot_residual_convergence(history, save_path='lm_convergence.png')

        # 步骤4: 保存所有参数
        #save_all_parameters(results, save_path='lm_parameters.pkl')

        print("\n" + "="*70)
        print("分析完成!")
        print("="*70)

    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())