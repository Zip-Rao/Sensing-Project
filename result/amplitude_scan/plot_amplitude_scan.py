#!/usr/bin/env python3
"""
绘制幅度扫描结果：RMSE vs 幅度曲线和线性失效边界图
符合科研制图规范

使用方法：
python plot_amplitude_scan.py [数据目录] [输出文件]

示例：
python plot_amplitude_scan.py result/amplitude_scan result/amplitude_scan_results.png
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import argparse
import glob
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

def set_plot_style(use_latex=False):
    """
    设置科研论文绘图样式

    参数:
    use_latex: 是否使用LaTeX渲染文本（需要系统安装LaTeX）
    """
    plt.style.use('seaborn-v0_8-whitegrid')

    # 基础样式设置
    params = {
        # 字体设置
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'DejaVu Serif', 'STIXGeneral'],
        'mathtext.fontset': 'stix',

        # 字体大小
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'legend.fontsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,

        # 图形设置
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,

        # 线宽和标记大小
        'lines.linewidth': 2.0,
        'lines.markersize': 8,

        # 坐标轴设置
        'axes.linewidth': 1.0,
        'axes.grid': True,
        'grid.linewidth': 0.5,
        'grid.alpha': 0.3,
        'grid.linestyle': ':',

        # 图例设置
        'legend.frameon': True,
        'legend.framealpha': 0.9,
        'legend.edgecolor': 'black',
        'legend.fancybox': False,
        'legend.loc': 'best',
    }

    if use_latex:
        try:
            params['text.usetex'] = True
            params['text.latex.preamble'] = r'\usepackage{amsmath}'
            print("使用LaTeX渲染文本")
        except:
            print("警告: LaTeX不可用，使用默认字体")
            use_latex = False

    plt.rcParams.update(params)
    return use_latex

def load_amplitude_data(data_dir):
    """
    加载所有幅度数据文件

    参数:
    data_dir: 数据目录

    返回:
    list: 数据字典列表，按幅度排序
    """
    data_files = glob.glob(os.path.join(data_dir, 'amplitude_*.npz'))
    if not data_files:
        raise FileNotFoundError(f"在目录 {data_dir} 中未找到幅度数据文件")

    all_data = []
    for file_path in data_files:
        try:
            data = dict(np.load(file_path, allow_pickle=True))
            # 从文件名提取幅度（如果数据中没有）
            if 'amplitude' not in data:
                filename = os.path.basename(file_path)
                # 提取幅度数值
                amp_str = filename.replace('amplitude_', '').replace('.npz', '')
                data['amplitude'] = float(amp_str)
            all_data.append(data)
        except Exception as e:
            print(f"警告: 无法加载文件 {file_path}: {e}")

    # 按幅度排序
    all_data.sort(key=lambda x: x['amplitude'])
    return all_data

def calculate_rmse(original_time, original_signal, recon_time, recon_signal):
    """
    计算重建信号的RMSE

    参数:
    original_time: 原始信号时间轴
    original_signal: 原始信号
    recon_time: 重建信号时间轴
    recon_signal: 重建信号

    返回:
    float: RMSE值
    """
    # 插值重建信号到原始时间网格
    f_recon = interp1d(recon_time, recon_signal, bounds_error=False, fill_value=0)
    recon_interp = f_recon(original_time)

    # 计算RMSE
    rmse = np.sqrt(np.mean((original_signal - recon_interp)**2))
    return rmse

def calculate_performance_metrics(all_data):
    """
    计算所有幅度的性能指标

    参数:
    all_data: 所有数据字典列表

    返回:
    tuple: (amplitudes, rmse_wiener, rmse_lm, nrmse_wiener, nrmse_lm, corr_wiener, corr_lm)
    """
    amplitudes = []
    rmse_wiener = []
    rmse_lm = []
    nrmse_wiener = []
    nrmse_lm = []
    corr_wiener = []
    corr_lm = []

    for data in all_data:
        amp = data['amplitude']
        original_time = data['original_time']
        original_signal = data['original_signal']

        # 计算Wiener RMSE
        rmse_w = calculate_rmse(original_time, original_signal,
                                data['wiener_time'], data['wiener_recon'])
        # 计算LM RMSE
        rmse_l = calculate_rmse(original_time, original_signal,
                                data['lm_time'], data['lm_recon'])

        # 计算归一化RMSE (NRMSE)
        signal_range = original_signal.max() - original_signal.min()
        if signal_range > 0:
            nrmse_w = rmse_w / signal_range
            nrmse_l = rmse_l / signal_range
        else:
            nrmse_w = nrmse_l = 0.0

        # 计算相关系数（需要插值）
        f_wiener = interp1d(data['wiener_time'], data['wiener_recon'],
                           bounds_error=False, fill_value=0)
        wiener_interp = f_wiener(original_time)
        corr_w = np.corrcoef(original_signal, wiener_interp)[0, 1]

        # LM已经在原始时间网格上
        corr_l = np.corrcoef(original_signal, data['lm_recon'])[0, 1]

        amplitudes.append(amp)
        rmse_wiener.append(rmse_w)
        rmse_lm.append(rmse_l)
        nrmse_wiener.append(nrmse_w)
        nrmse_lm.append(nrmse_l)
        corr_wiener.append(corr_w)
        corr_lm.append(corr_l)

    return (np.array(amplitudes), np.array(rmse_wiener), np.array(rmse_lm),
            np.array(nrmse_wiener), np.array(nrmse_lm),
            np.array(corr_wiener), np.array(corr_lm))

def find_crossover_point(amplitudes, rmse_wiener, rmse_lm):
    """
    找到Wiener开始不如LM的交叉点（RMSE_Wiener > RMSE_LM）

    参数:
    amplitudes: 幅度数组
    rmse_wiener: Wiener RMSE数组
    rmse_lm: LM RMSE数组

    返回:
    tuple: (crossover_amplitude, crossover_index)
           如果未找到交叉点，返回(None, None)
    """
    # 计算RMSE比值
    ratio = rmse_wiener / rmse_lm

    # 找到第一个比值大于1的点
    for i in range(len(ratio)):
        if ratio[i] > 1.0:
            return amplitudes[i], i

    # 如果未找到，尝试插值
    if len(amplitudes) > 1:
        # 找到比值最接近1的点
        idx_min = np.argmin(np.abs(ratio - 1.0))
        if 0 < idx_min < len(amplitudes) - 1:
            # 使用三点插值更精确地估计交叉点
            x = amplitudes[idx_min-1:idx_min+2]
            y = ratio[idx_min-1:idx_min+2] - 1.0
            # 线性插值求根
            if y[0] * y[2] < 0:  # 根在区间内
                # 线性插值
                root = x[0] - y[0] * (x[1] - x[0]) / (y[1] - y[0])
                return root, idx_min

    return None, None

def create_amplitude_scan_figure(amplitudes, rmse_wiener, rmse_lm,
                                 nrmse_wiener, nrmse_lm,
                                 corr_wiener, corr_lm,
                                 output_path=None, use_latex=False):
    """
    创建幅度扫描结果图

    参数:
    amplitudes: 幅度数组
    rmse_wiener, rmse_lm: RMSE数组
    nrmse_wiener, nrmse_lm: 归一化RMSE数组
    corr_wiener, corr_lm: 相关系数数组
    output_path: 输出文件路径，如果为None则显示而不保存
    use_latex: 是否使用LaTeX渲染文本
    """
    # 创建图形，1行2列
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    # 颜色设置（遵循ColorBrewer Set1）
    wiener_color = '#e41a1c'  # 红色
    lm_color = '#377eb8'      # 蓝色
    crossover_color = '#984ea3'  # 紫色
    threshold_color = '#ff7f00'  # 橙色

    # 标记样式
    wiener_marker = 'o'
    lm_marker = 's'

    # 子图1：RMSE vs 幅度（对数坐标）
    ax1 = axes[0]
    ax1.semilogx(amplitudes, rmse_wiener, color=wiener_color,
                 marker=wiener_marker, label='Wiener', linewidth=2)
    ax1.semilogx(amplitudes, rmse_lm, color=lm_color,
                 marker=lm_marker, label='LM', linewidth=2)

    # 添加平滑曲线（可选）
    if len(amplitudes) >= 5:
        try:
            # 使用Savitzky-Golay滤波器平滑
            window = min(5, len(amplitudes) // 2 * 2 + 1)  # 奇数
            rmse_w_smooth = savgol_filter(rmse_wiener, window, 2)
            rmse_l_smooth = savgol_filter(rmse_lm, window, 2)
            ax1.semilogx(amplitudes, rmse_w_smooth, color=wiener_color,
                        linestyle='--', alpha=0.5, linewidth=1.5)
            ax1.semilogx(amplitudes, rmse_l_smooth, color=lm_color,
                        linestyle='--', alpha=0.5, linewidth=1.5)
        except:
            pass  # 平滑失败，使用原始数据

    ax1.set_xlabel('Signal Amplitude (a.u.)', fontsize=11)
    ax1.set_ylabel('RMSE (a.u.)', fontsize=11)
    ax1.set_title('(a) Reconstruction Error vs Amplitude', fontsize=12, pad=10)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(fontsize=10)

    # 找到交叉点并标记
    crossover_amp, crossover_idx = find_crossover_point(amplitudes, rmse_wiener, rmse_lm)
    if crossover_amp is not None:
        ax1.axvline(x=crossover_amp, color=crossover_color, linestyle='--',
                   alpha=0.7, linewidth=1.5, label=f'Crossover: {crossover_amp:.4f}')
        # 在交叉点添加标记
        ax1.plot(crossover_amp, rmse_wiener[crossover_idx], 'x',
                color=crossover_color, markersize=12, markeredgewidth=2)
        ax1.legend(fontsize=10)  # 重新绘制图例以包含交叉点线

    # 子图2：RMSE比值和相关系数
    ax2 = axes[1]

    # 计算RMSE比值
    rmse_ratio = rmse_wiener / rmse_lm

    # 绘制RMSE比值（左侧y轴）
    ax2.semilogx(amplitudes, rmse_ratio, color=crossover_color,
                marker='^', linewidth=2, label='RMSE Ratio (Wiener/LM)')

    # 添加水平线y=1表示相等性能
    ax2.axhline(y=1.0, color=threshold_color, linestyle='--',
               alpha=0.7, linewidth=1.5, label='Equal Performance')

    ax2.set_xlabel('Signal Amplitude (a.u.)', fontsize=11)
    ax2.set_ylabel('RMSE Ratio (Wiener / LM)', fontsize=11, color=crossover_color)
    ax2.tick_params(axis='y', labelcolor=crossover_color)

    # 创建第二个y轴用于相关系数
    ax2b = ax2.twinx()
    ax2b.semilogx(amplitudes, corr_wiener, color=wiener_color,
                 marker=wiener_marker, linestyle=':', alpha=0.7, label='Wiener Corr')
    ax2b.semilogx(amplitudes, corr_lm, color=lm_color,
                 marker=lm_marker, linestyle=':', alpha=0.7, label='LM Corr')

    ax2b.set_ylabel('Correlation Coefficient', fontsize=11, color='black')
    ax2b.tick_params(axis='y', labelcolor='black')
    ax2b.set_ylim(0, 1.05)  # 相关系数范围

    ax2.set_title('(b) Performance Metrics vs Amplitude', fontsize=12, pad=10)
    ax2.grid(True, which='both', alpha=0.3)

    # 合并图例（需要特殊处理）
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='lower left')

    # 标记交叉点
    if crossover_amp is not None:
        ax2.axvline(x=crossover_amp, color=crossover_color, linestyle='--',
                   alpha=0.7, linewidth=1.5)
        ax2.text(crossover_amp, ax2.get_ylim()[1] * 0.95,
                f'Crossover\n{crossover_amp:.4f}',
                color=crossover_color, fontsize=9,
                ha='center', va='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # 保存或显示
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=300)
        print(f"图形已保存到: {output_path}")

        # 同时保存PDF版本（矢量图，适合出版物）
        pdf_path = output_path.replace('.png', '.pdf').replace('.jpg', '.pdf')
        plt.savefig(pdf_path, dpi=300, format='pdf')
        print(f"PDF版本已保存到: {pdf_path}")
    else:
        plt.show()

    return fig, axes, crossover_amp

def print_summary_table(amplitudes, rmse_wiener, rmse_lm, crossover_amp, output_file=None):
    """
    打印和保存摘要表格

    参数:
    amplitudes: 幅度数组
    rmse_wiener, rmse_lm: RMSE数组
    crossover_amp: 交叉点幅度
    output_file: 输出文本文件路径
    """
    print("\n" + "="*80)
    print("幅度扫描结果摘要")
    print("="*80)
    print(f"{'Amplitude':<12} {'RMSE Wiener':<15} {'RMSE LM':<15} {'Ratio (W/L)':<12} {'Better Method':<15}")
    print("-"*80)

    for i, amp in enumerate(amplitudes):
        ratio = rmse_wiener[i] / rmse_lm[i]
        better = "Wiener" if ratio < 1 else "LM" if ratio > 1 else "Equal"
        print(f"{amp:<12.6f} {rmse_wiener[i]:<15.6f} {rmse_lm[i]:<15.6f} {ratio:<12.4f} {better:<15}")

    print("-"*80)
    if crossover_amp:
        print(f"线性失效边界（Wiener开始不如LM的幅度）: {crossover_amp:.6f}")
        # 找到最接近的幅度点
        idx = np.argmin(np.abs(amplitudes - crossover_amp))
        print(f"  在幅度 {amplitudes[idx]:.6f} 处:")
        print(f"    RMSE Wiener: {rmse_wiener[idx]:.6f}")
        print(f"    RMSE LM: {rmse_lm[idx]:.6f}")
        print(f"    Ratio: {rmse_wiener[idx]/rmse_lm[idx]:.4f}")
    else:
        print("未找到线性失效边界（Wiener在所有测试幅度上都优于或等于LM）")

    # 保存为文本文件
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("幅度扫描结果摘要\n")
            f.write("="*80 + "\n")
            f.write(f"{'Amplitude':<12} {'RMSE Wiener':<15} {'RMSE LM':<15} {'Ratio (W/L)':<12} {'Better Method':<15}\n")
            f.write("-"*80 + "\n")

            for i, amp in enumerate(amplitudes):
                ratio = rmse_wiener[i] / rmse_lm[i]
                better = "Wiener" if ratio < 1 else "LM" if ratio > 1 else "Equal"
                f.write(f"{amp:<12.6f} {rmse_wiener[i]:<15.6f} {rmse_lm[i]:<15.6f} {ratio:<12.4f} {better:<15}\n")

            f.write("-"*80 + "\n")
            if crossover_amp:
                f.write(f"线性失效边界（Wiener开始不如LM的幅度）: {crossover_amp:.6f}\n")
                idx = np.argmin(np.abs(amplitudes - crossover_amp))
                f.write(f"  在幅度 {amplitudes[idx]:.6f} 处:\n")
                f.write(f"    RMSE Wiener: {rmse_wiener[idx]:.6f}\n")
                f.write(f"    RMSE LM: {rmse_lm[idx]:.6f}\n")
                f.write(f"    Ratio: {rmse_wiener[idx]/rmse_lm[idx]:.4f}\n")
            else:
                f.write("未找到线性失效边界（Wiener在所有测试幅度上都优于或等于LM）\n")
        print(f"摘要表格已保存到: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='绘制幅度扫描结果图')
    parser.add_argument('data_dir', nargs='?', default='result/amplitude_scan',
                       help='数据目录（默认：result/amplitude_scan）')
    parser.add_argument('--output', '-o', default='result/amplitude_scan_results.png',
                       help='输出文件路径（默认：result/amplitude_scan_results.png）')
    parser.add_argument('--latex', action='store_true',
                       help='使用LaTeX渲染文本（需要系统安装LaTeX）')
    parser.add_argument('--min-amplitude', type=float, default=0.001,
                       help='最小幅度（用于过滤数据，默认：0.001）')
    parser.add_argument('--max-amplitude', type=float, default=0.1,
                       help='最大幅度（用于过滤数据，默认：0.1）')

    args = parser.parse_args()

    # 设置绘图样式
    use_latex = set_plot_style(args.latex)

    # 检查数据目录
    if not os.path.exists(args.data_dir):
        print(f"错误：数据目录 {args.data_dir} 不存在")
        print("请先运行 generate_amplitude_scan_data.py 生成数据")
        sys.exit(1)

    # 加载数据
    print("加载幅度扫描数据...")
    try:
        all_data = load_amplitude_data(args.data_dir)
    except Exception as e:
        print(f"加载数据时出错: {e}")
        sys.exit(1)

    print(f"成功加载 {len(all_data)} 个幅度数据文件")

    # 计算性能指标
    print("计算性能指标...")
    amplitudes, rmse_wiener, rmse_lm, nrmse_wiener, nrmse_lm, corr_wiener, corr_lm = \
        calculate_performance_metrics(all_data)

    # 过滤幅度范围
    mask = (amplitudes >= args.min_amplitude) & (amplitudes <= args.max_amplitude)
    if not np.any(mask):
        print(f"错误：在范围 [{args.min_amplitude}, {args.max_amplitude}] 内没有数据")
        sys.exit(1)

    amplitudes = amplitudes[mask]
    rmse_wiener = rmse_wiener[mask]
    rmse_lm = rmse_lm[mask]
    nrmse_wiener = nrmse_wiener[mask]
    nrmse_lm = nrmse_lm[mask]
    corr_wiener = corr_wiener[mask]
    corr_lm = corr_lm[mask]

    print(f"分析范围: {len(amplitudes)} 个幅度点，从 {amplitudes[0]:.6f} 到 {amplitudes[-1]:.6f}")

    # 创建图形
    print("创建幅度扫描结果图...")
    try:
        fig, axes, crossover_amp = create_amplitude_scan_figure(
            amplitudes, rmse_wiener, rmse_lm,
            nrmse_wiener, nrmse_lm,
            corr_wiener, corr_lm,
            args.output, use_latex
        )
    except Exception as e:
        print(f"创建图形时出错: {e}")
        sys.exit(1)

    # 打印和保存摘要
    print("\n计算摘要信息...")
    try:
        summary_file = args.output.replace('.png', '_summary.txt').replace('.jpg', '_summary.txt')
        print_summary_table(amplitudes, rmse_wiener, rmse_lm, crossover_amp, summary_file)
    except Exception as e:
        print(f"生成摘要时出错: {e}")
        # 继续执行，不影响主图生成

    print(f"\n处理完成！")
    print(f"- 结果图: {args.output}")
    print(f"- PDF版本: {args.output.replace('.png', '.pdf').replace('.jpg', '.pdf')}")
    if 'summary_file' in locals():
        print(f"- 摘要文件: {summary_file}")

    # 显示提示信息
    print("\n提示：")
    print("1. 如需调整图形，可以修改脚本中的绘图参数")
    print("2. 对于出版物，建议使用PDF格式（矢量图）")
    print("3. 可以尝试 --latex 选项获得更好的排版效果（需要系统安装LaTeX）")
    print("4. 使用 --min-amplitude 和 --max-amplitude 参数可以限制分析范围")

if __name__ == '__main__':
    main()