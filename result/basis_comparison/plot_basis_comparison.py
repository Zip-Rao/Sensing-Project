#!/usr/bin/env python3
"""
基函数对比绘图脚本
绘制B-spline、Fourier和Legendre基函数的性能对比图，符合科研制图规范

包括：
1. 收敛曲线对比（残差 vs 迭代次数）
2. 最终RMSE对比
3. 信号重建对比（可选）

使用方法：
python plot_basis_comparison.py [数据文件] [输出文件]

示例：
python plot_basis_comparison.py result/basis_comparison/basis_comparison_data.npz result/basis_comparison/basis_comparison_figure.png
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import json
import argparse
from scipy.interpolate import interp1d
from matplotlib import rcParams

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

def load_data(data_path):
    """
    加载数据文件，支持NPZ和JSON格式

    参数:
    data_path: 数据文件路径

    返回:
    dict: 加载的数据字典
    """
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"数据文件不存在: {data_path}")

    if data_path.endswith('.npz'):
        # 加载NPZ文件
        data = np.load(data_path, allow_pickle=True)
        data_dict = {}
        for key in data.files:
            data_dict[key] = data[key]

        # 尝试加载对应的JSON文件获取元数据
        json_path = data_path.replace('.npz', '.json')
        if os.path.exists(json_path):
            with open(json_path, 'r', encoding='utf-8') as f:
                metadata = json.load(f)
            data_dict['metadata'] = metadata
        else:
            data_dict['metadata'] = {'source': 'NPZ file only'}

    elif data_path.endswith('.json'):
        # 加载JSON文件
        with open(data_path, 'r', encoding='utf-8') as f:
            data_dict = json.load(f)
    else:
        raise ValueError("不支持的文件格式，请使用.npz或.json文件")

    return data_dict

def extract_basis_results(data_dict):
    """
    从数据字典中提取基函数结果

    参数:
    data_dict: 数据字典

    返回:
    dict: 基函数结果字典
    """
    basis_results = {}

    # 检查数据结构
    if 'metadata' in data_dict and 'basis_results' in data_dict.get('metadata', {}):
        # JSON格式数据
        for basis_type, result in data_dict['metadata']['basis_results'].items():
            basis_results[basis_type] = result
    else:
        # NPZ格式数据
        for key in data_dict.keys():
            if key.endswith('_residual_norms'):
                basis_type = key.replace('_residual_norms', '')
                if basis_type in ['bspline', 'fourier', 'legendre']:
                    basis_results[basis_type] = {
                        'residual_norms': data_dict[key],
                        'rmse': data_dict.get(f'{basis_type}_rmse', np.nan),
                        'reconstructed_signal': data_dict.get(f'{basis_type}_recon', None)
                    }

    # 确保我们有原始信号数据
    if 'original_time' in data_dict and 'original_signal' in data_dict:
        basis_results['original'] = {
            'time': data_dict['original_time'],
            'signal': data_dict['original_signal']
        }

    return basis_results

def create_convergence_plot(basis_results, ax=None):
    """
    创建收敛曲线对比图

    参数:
    basis_results: 基函数结果字典
    ax: 可选的Axes对象

    返回:
    matplotlib.axes.Axes: 绘图的Axes对象
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    # 定义颜色和标记
    colors = {
        'bspline': '#e41a1c',   # 红色
        'fourier': '#377eb8',   # 蓝色
        'legendre': '#4daf4a',  # 绿色
    }

    markers = {
        'bspline': 'o',
        'fourier': 's',
        'legendre': '^',
    }

    labels = {
        'bspline': 'B-spline',
        'fourier': 'Fourier',
        'legendre': 'Legendre',
    }

    # 绘制每种基函数的收敛曲线
    for basis_type in ['bspline', 'fourier', 'legendre']:
        if basis_type in basis_results:
            result = basis_results[basis_type]
            if 'residual_norms' in result:
                residuals = result['residual_norms']
                iterations = np.arange(1, len(residuals) + 1)

                # 使用对数坐标显示残差
                ax.semilogy(iterations, residuals,
                           color=colors[basis_type],
                           marker=markers[basis_type],
                           markevery=max(1, len(iterations)//10),
                           label=labels[basis_type],
                           linewidth=2,
                           markersize=8)

    # 设置坐标轴标签
    ax.set_xlabel('Iteration Number', fontsize=11)
    ax.set_ylabel('Residual Norm (log scale)', fontsize=11)
    ax.set_title('(a) Convergence History: Residual Norm vs Iteration', fontsize=12, pad=10)

    # 添加网格
    ax.grid(True, which='both', alpha=0.3)

    # 添加图例
    ax.legend(fontsize=10, loc='best')

    # 设置坐标轴范围
    ax.set_xlim(left=1)

    return ax

def create_rmse_comparison_plot(basis_results, ax=None):
    """
    创建RMSE对比柱状图

    参数:
    basis_results: 基函数结果字典
    ax: 可选的Axes对象

    返回:
    matplotlib.axes.Axes: 绘图的Axes对象
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    # 定义颜色
    colors = {
        'bspline': '#e41a1c',   # 红色
        'fourier': '#377eb8',   # 蓝色
        'legendre': '#4daf4a',  # 绿色
    }

    labels = {
        'bspline': 'B-spline',
        'fourier': 'Fourier',
        'legendre': 'Legendre',
    }

    # 收集RMSE数据
    basis_types = []
    rmse_values = []
    bar_colors = []

    for basis_type in ['bspline', 'fourier', 'legendre']:
        if basis_type in basis_results:
            result = basis_results[basis_type]
            if 'rmse' in result:
                basis_types.append(labels[basis_type])
                rmse_values.append(result['rmse'])
                bar_colors.append(colors[basis_type])

    if not rmse_values:
        ax.text(0.5, 0.5, 'No RMSE data available',
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=12)
        return ax

    # 创建柱状图
    x_pos = np.arange(len(basis_types))
    bars = ax.bar(x_pos, rmse_values, color=bar_colors, edgecolor='black', alpha=0.8)

    # 添加数值标签
    for i, (bar, rmse) in enumerate(zip(bars, rmse_values)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.01 * max(rmse_values),
                f'{rmse:.4f}', ha='center', va='bottom', fontsize=9)

    # 设置坐标轴标签
    ax.set_xlabel('Basis Function Type', fontsize=11)
    ax.set_ylabel('RMSE', fontsize=11)
    ax.set_title('(b) Reconstruction Error Comparison', fontsize=12, pad=10)

    # 设置x轴刻度
    ax.set_xticks(x_pos)
    ax.set_xticklabels(basis_types, rotation=0)

    # 添加网格
    ax.grid(True, axis='y', alpha=0.3, linestyle='--')

    # 找出最佳基函数（最小RMSE）
    if len(rmse_values) > 0:
        best_idx = np.argmin(rmse_values)
        best_rmse = rmse_values[best_idx]
        best_type = basis_types[best_idx]

        # 添加标注
        ax.annotate(f'Best: {best_type}\nRMSE = {best_rmse:.4f}',
                   xy=(best_idx, best_rmse),
                   xytext=(best_idx, best_rmse * 1.2),
                   arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                   fontsize=10,
                   ha='center')

    return ax

def create_signal_comparison_plot(basis_results, ax=None):
    """
    创建信号重建对比图

    参数:
    basis_results: 基函数结果字典
    ax: 可选的Axes对象

    返回:
    matplotlib.axes.Axes: 绘图的Axes对象
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    # 检查是否有原始信号
    if 'original' not in basis_results:
        ax.text(0.5, 0.5, 'No original signal data available',
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=12)
        return ax

    # 获取原始信号
    original_time = basis_results['original']['time']
    original_signal = basis_results['original']['signal']

    # 定义颜色和线型
    colors = {
        'original': 'black',
        'bspline': '#e41a1c',   # 红色
        'fourier': '#377eb8',   # 蓝色
        'legendre': '#4daf4a',  # 绿色
    }

    linestyles = {
        'original': '--',
        'bspline': '-',
        'fourier': '-',
        'legendre': '-',
    }

    labels = {
        'original': 'Original Signal',
        'bspline': 'B-spline Reconstruction',
        'fourier': 'Fourier Reconstruction',
        'legendre': 'Legendre Reconstruction',
    }

    # 绘制原始信号
    ax.plot(original_time, original_signal,
            color=colors['original'],
            linestyle=linestyles['original'],
            linewidth=2.5,
            label=labels['original'],
            alpha=0.8)

    # 绘制重建信号
    for basis_type in ['bspline', 'fourier', 'legendre']:
        if basis_type in basis_results:
            result = basis_results[basis_type]
            if 'reconstructed_signal' in result and result['reconstructed_signal'] is not None:
                recon_signal = result['reconstructed_signal']
                # 如果时间轴不同，需要插值
                if 'reconstructed_time' in result:
                    recon_time = result['reconstructed_time']
                else:
                    recon_time = original_time  # 假设时间轴相同

                # 确保信号长度匹配
                if len(recon_time) != len(recon_signal):
                    # 简单处理：截断或填充
                    min_len = min(len(recon_time), len(recon_signal))
                    recon_time = recon_time[:min_len]
                    recon_signal = recon_signal[:min_len]

                ax.plot(recon_time, recon_signal,
                       color=colors[basis_type],
                       linestyle=linestyles[basis_type],
                       linewidth=1.5,
                       label=labels[basis_type],
                       alpha=0.7)

    # 设置坐标轴标签
    ax.set_xlabel('Time (ns)', fontsize=11)
    ax.set_ylabel('Signal Amplitude (a.u.)', fontsize=11)
    ax.set_title('(c) Signal Reconstruction Comparison', fontsize=12, pad=10)

    # 添加网格
    ax.grid(True, alpha=0.3)

    # 添加图例
    ax.legend(fontsize=10, loc='best')

    return ax

def create_lm_convergence_detail_plot(basis_results, selected_basis='bspline', ax=None):
    """
    创建LM算法收敛细节图（残差 vs 迭代次数）

    参数:
    basis_results: 基函数结果字典
    selected_basis: 要详细显示的基函数类型
    ax: 可选的Axes对象

    返回:
    matplotlib.axes.Axes: 绘图的Axes对象
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    if selected_basis not in basis_results:
        ax.text(0.5, 0.5, f'No data for {selected_basis} basis',
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=12)
        return ax

    result = basis_results[selected_basis]
    if 'residual_norms' not in result:
        ax.text(0.5, 0.5, f'No convergence data for {selected_basis} basis',
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=12)
        return ax

    # 获取收敛数据
    residuals = result['residual_norms']
    iterations = np.arange(1, len(residuals) + 1)

    # 颜色映射
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(residuals)))

    # 绘制收敛曲线
    ax.semilogy(iterations, residuals, 'o-', color='#984ea3',
               linewidth=2, markersize=8, markevery=max(1, len(iterations)//10))

    # 标记收敛行为
    if len(residuals) > 1:
        # 计算收敛率（最后5次迭代）
        if len(residuals) >= 5:
            last_residuals = residuals[-5:]
            convergence_rate = np.mean(np.diff(np.log10(last_residuals)))
            convergence_text = f'Convergence rate: {convergence_rate:.3f}'
        else:
            convergence_rate = np.diff(np.log10(residuals))[-1] if len(residuals) > 1 else 0
            convergence_text = f'Final convergence: {convergence_rate:.3f}'

        # 添加收敛信息
        ax.text(0.02, 0.98, convergence_text,
               transform=ax.transAxes, fontsize=10,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # 设置坐标轴标签
    ax.set_xlabel('Iteration Number', fontsize=11)
    ax.set_ylabel('Residual Norm (log scale)', fontsize=11)

    basis_labels = {
        'bspline': 'B-spline',
        'fourier': 'Fourier',
        'legendre': 'Legendre'
    }
    basis_label = basis_labels.get(selected_basis, selected_basis)
    ax.set_title(f'LM Algorithm Convergence: {basis_label} Basis', fontsize=12, pad=10)

    # 添加网格
    ax.grid(True, which='both', alpha=0.3)

    # 添加收敛阈值线（如果知道）
    if 'rmse' in result:
        ax.axhline(y=result['rmse'], color='r', linestyle='--', alpha=0.5,
                  label=f'Final RMSE: {result["rmse"]:.4f}')
        ax.legend(fontsize=10)

    return ax

def create_comprehensive_figure(basis_results, output_path=None, use_latex=False):
    """
    创建综合对比图

    参数:
    basis_results: 基函数结果字典
    output_path: 输出文件路径
    use_latex: 是否使用LaTeX渲染文本

    返回:
    matplotlib.figure.Figure: 创建的图形对象
    """
    # 设置绘图样式
    set_plot_style(use_latex)

    # 创建图形，2行2列
    fig = plt.figure(figsize=(14, 10), constrained_layout=True)

    # 创建子图网格
    gs = fig.add_gridspec(2, 2, hspace=0.15, wspace=0.15)

    # 子图1: 收敛曲线对比
    ax1 = fig.add_subplot(gs[0, 0])
    create_convergence_plot(basis_results, ax1)

    # 子图2: RMSE对比
    ax2 = fig.add_subplot(gs[0, 1])
    create_rmse_comparison_plot(basis_results, ax2)

    # 子图3: 信号重建对比
    ax3 = fig.add_subplot(gs[1, :])
    create_signal_comparison_plot(basis_results, ax3)

    # 添加总体标题
    fig.suptitle('Basis Function Comparison: B-spline vs Fourier vs Legendre',
                fontsize=14, y=1.02)

    # 保存或显示
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"综合对比图已保存到: {output_path}")

        # 同时保存PDF版本
        pdf_path = output_path.replace('.png', '.pdf').replace('.jpg', '.pdf')
        plt.savefig(pdf_path, dpi=300, format='pdf', bbox_inches='tight')
        print(f"PDF版本已保存到: {pdf_path}")
    else:
        plt.show()

    return fig

def create_convergence_figure(basis_results, output_path=None, use_latex=False):
    """
    创建收敛历史细节图

    参数:
    basis_results: 基函数结果字典
    output_path: 输出文件路径
    use_latex: 是否使用LaTeX渲染文本

    返回:
    matplotlib.figure.Figure: 创建的图形对象
    """
    # 设置绘图样式
    set_plot_style(use_latex)

    # 创建图形，1行3列
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    # 为每种基函数创建收敛细节图
    basis_types = ['bspline', 'fourier', 'legendre']
    titles = ['(a) B-spline Basis', '(b) Fourier Basis', '(c) Legendre Basis']

    for idx, (basis_type, title) in enumerate(zip(basis_types, titles)):
        ax = axes[idx]
        if basis_type in basis_results:
            create_lm_convergence_detail_plot(basis_results, basis_type, ax)
            ax.set_title(title, fontsize=12, pad=10)
        else:
            ax.text(0.5, 0.5, f'No data for {basis_type}',
                   horizontalalignment='center', verticalalignment='center',
                   transform=ax.transAxes, fontsize=12)
            ax.set_title(title, fontsize=12, pad=10)

    # 添加总体标题
    fig.suptitle('LM Algorithm Convergence History for Different Basis Functions',
                fontsize=14, y=1.02)

    # 保存或显示
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"收敛历史图已保存到: {output_path}")

        # 同时保存PDF版本
        pdf_path = output_path.replace('.png', '.pdf').replace('.jpg', '.pdf')
        plt.savefig(pdf_path, dpi=300, format='pdf', bbox_inches='tight')
        print(f"PDF版本已保存到: {pdf_path}")
    else:
        plt.show()

    return fig

def print_data_summary(basis_results):
    """打印数据摘要"""
    print("\n" + "="*80)
    print("数据摘要")
    print("="*80)

    # 检查原始信号
    if 'original' in basis_results:
        orig_time = basis_results['original']['time']
        orig_signal = basis_results['original']['signal']
        if hasattr(orig_time, '__len__'):
            print(f"原始信号: {len(orig_time)} 个时间点")
            if hasattr(orig_signal, '__len__'):
                print(f"  时间范围: {orig_time[0]:.1f} - {orig_time[-1]:.1f} ns")
                print(f"  信号范围: {np.min(orig_signal):.4f} - {np.max(orig_signal):.4f}")

    # 检查基函数结果
    basis_types = ['bspline', 'fourier', 'legendre']
    print("\n基函数结果:")
    for basis_type in basis_types:
        if basis_type in basis_results:
            result = basis_results[basis_type]
            print(f"\n{basis_type.upper()}:")

            if 'rmse' in result:
                print(f"  RMSE: {result['rmse']:.6f}")

            if 'residual_norms' in result:
                residuals = result['residual_norms']
                print(f"  收敛迭代次数: {len(residuals)}")
                if len(residuals) > 1:
                    print(f"  初始残差: {residuals[0]:.6f}")
                    print(f"  最终残差: {residuals[-1]:.6f}")
                    improvement = residuals[0] / residuals[-1] if residuals[-1] > 0 else float('inf')
                    print(f"  改进倍数: {improvement:.2f}")

    print("="*80)

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='绘制基函数对比图')
    parser.add_argument('data_file', nargs='?', default='result/basis_comparison/basis_comparison_data.npz',
                       help='数据文件路径（默认：result/basis_comparison/basis_comparison_data.npz）')
    parser.add_argument('--output', '-o', default='result/basis_comparison/basis_comparison_figure.png',
                       help='输出文件路径（默认：result/basis_comparison/basis_comparison_figure.png）')
    parser.add_argument('--convergence-output', '-c', default='result/basis_comparison/convergence_history.png',
                       help='收敛历史图输出路径（默认：result/basis_comparison/convergence_history.png）')
    parser.add_argument('--latex', action='store_true',
                       help='使用LaTeX渲染文本（需要系统安装LaTeX）')
    parser.add_argument('--show', action='store_true',
                       help='显示图形而不保存')
    parser.add_argument('--summary', action='store_true',
                       help='打印数据摘要')

    args = parser.parse_args()

    print("="*80)
    print("基函数对比绘图")
    print("="*80)
    print(f"数据文件: {args.data_file}")
    print(f"输出文件: {args.output}")
    print(f"收敛历史图: {args.convergence_output}")

    # 检查数据文件
    if not os.path.exists(args.data_file):
        print(f"错误：数据文件 {args.data_file} 不存在")
        print("请先运行 generate_basis_comparison_data.py 生成数据")
        sys.exit(1)

    # 加载数据
    print("\n加载数据...")
    try:
        data_dict = load_data(args.data_file)
        print(f"成功加载数据文件: {args.data_file}")
    except Exception as e:
        print(f"加载数据时出错: {e}")
        sys.exit(1)

    # 提取基函数结果
    print("提取基函数结果...")
    basis_results = extract_basis_results(data_dict)

    if not any(basis_type in basis_results for basis_type in ['bspline', 'fourier', 'legendre']):
        print("错误：数据中未找到基函数结果")
        print("请确保数据文件包含 bspline、fourier 或 legendre 基函数的结果")
        sys.exit(1)

    # 打印摘要
    if args.summary:
        print_data_summary(basis_results)

    # 创建图形
    print("\n创建图形...")

    if not args.show:
        # 创建综合对比图
        try:
            print("创建综合对比图...")
            fig1 = create_comprehensive_figure(
                basis_results,
                output_path=args.output if not args.show else None,
                use_latex=args.latex
            )
        except Exception as e:
            print(f"创建综合对比图时出错: {e}")
            import traceback
            traceback.print_exc()

        # 创建收敛历史图
        try:
            print("创建收敛历史图...")
            fig2 = create_convergence_figure(
                basis_results,
                output_path=args.convergence_output if not args.show else None,
                use_latex=args.latex
            )
        except Exception as e:
            print(f"创建收敛历史图时出错: {e}")
            import traceback
            traceback.print_exc()

        print("\n" + "="*80)
        print("绘图完成!")
        if not args.show:
            print(f"- 综合对比图: {args.output}")
            print(f"- PDF版本: {args.output.replace('.png', '.pdf').replace('.jpg', '.pdf')}")
            print(f"- 收敛历史图: {args.convergence_output}")
            print(f"- 收敛历史PDF: {args.convergence_output.replace('.png', '.pdf').replace('.jpg', '.pdf')}")
        print("="*80)
    else:
        # 显示图形
        print("显示图形...")
        create_comprehensive_figure(basis_results, use_latex=args.latex)
        create_convergence_figure(basis_results, use_latex=args.latex)
        plt.show()

if __name__ == '__main__':
    main()