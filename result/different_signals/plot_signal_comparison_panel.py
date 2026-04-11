#!/usr/bin/env python3
"""
绘制四种信号类型的Wiener和LM重建结果对比图
2行×4列布局，符合科研制图规范

使用方法：
python plot_signal_comparison_panel.py [数据目录] [输出文件]

示例：
python plot_signal_comparison_panel.py result/different_signals result/comparison.png
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import argparse

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
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'legend.fontsize': 9,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,

        # 图形设置
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,

        # 线宽和标记大小
        'lines.linewidth': 1.5,
        'lines.markersize': 4,

        # 坐标轴设置
        'axes.linewidth': 0.8,
        'axes.grid': True,
        'grid.linewidth': 0.5,
        'grid.alpha': 0.3,
        'grid.linestyle': ':',

        # 图例设置
        'legend.frameon': True,
        'legend.framealpha': 0.9,
        'legend.edgecolor': 'black',
        'legend.fancybox': False,
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

def load_signal_data(signal_name, data_dir):
    """加载单个信号的数据"""
    data_file = os.path.join(data_dir, f'signal_{signal_name}.npz')
    if not os.path.exists(data_file):
        # 尝试从合并文件中加载
        combined_file = os.path.join(data_dir, 'all_signals_data.npz')
        if os.path.exists(combined_file):
            data = np.load(combined_file, allow_pickle=True)
            # 提取该信号的数据
            signal_data = {}
            prefix = f'{signal_name}_'
            for key in data.files:
                if key.startswith(prefix):
                    new_key = key[len(prefix):]
                    signal_data[new_key] = data[key]
            if signal_data:
                return signal_data
        raise FileNotFoundError(f"数据文件 {data_file} 不存在")

    return dict(np.load(data_file, allow_pickle=True))

def create_panel_figure(signal_names, data_dir, output_path=None, use_latex=False):
    """
    创建2×4 panel对比图，符合科研规范

    参数:
    signal_names: 信号名称列表，按顺序对应各列
    data_dir: 数据目录
    output_path: 输出文件路径，如果为None则显示而不保存
    use_latex: 是否使用LaTeX渲染文本
    """

    if len(signal_names) != 4:
        raise ValueError("需要恰好4种信号类型")

    # 加载所有信号数据
    signals_data = {}
    for name in signal_names:
        signals_data[name] = load_signal_data(name, data_dir)

    # 创建图形，2行4列，增加高度为底部图例留出空间
    # 各列有自己的坐标轴（不共享y轴）
    fig, axes = plt.subplots(2, 4, figsize=(12, 6.0),
                             sharex='col', sharey=False,
                             constrained_layout=True)
    # 调整布局，为底部图例留出空间
    fig.tight_layout(rect=[0, 0.15, 1, 0.92])

    # 信号类型描述（用于标题）
    signal_descriptions = {
        'sine': 'Sinusoidal',
        'step': 'Step-like',
        'double_peak': 'Double-peak',
        'complex': 'Complex'
    }

    # 颜色和线型设置（遵循ColorBrewer Set1，避免红绿色盲问题）
    original_color = 'black'           # 原始信号：黑色
    wiener_color = '#e41a1c'           # Wiener重建：红色
    lm_color = '#377eb8'               # LM重建：蓝色

    original_style = '--'              # 原始信号：虚线
    wiener_style = '-'                 # Wiener重建：实线
    lm_style = '-'                    # LM重建：点划线

    line_widths = [1.5, 1.8, 2.0]      # 线宽：原始，Wiener，LM

    # 遍历所有子图
    for row in range(2):  # 行：重建方法 (0: Wiener, 1: LM)
        for col in range(4):  # 列：信号类型
            ax = axes[row, col]
            signal_name = signal_names[col]
            data = signals_data[signal_name]

            # 提取数据
            original_time = data['original_time']
            original_signal = data['original_signal']

            if row == 0:  # Wiener重建
                recon_time = data['wiener_time']
                recon_signal = data['wiener_recon']
                recon_color = wiener_color
                recon_style = wiener_style
                recon_label = 'Wiener' if col == 0 else ''
            else:  # LM重建
                recon_time = data['lm_time']
                recon_signal = data['lm_recon']
                recon_color = lm_color
                recon_style = lm_style
                recon_label = 'LM' if col == 0 else ''

            # 绘制原始信号（所有子图都绘制，但只在第一列图例中显示）
            ax.plot(original_time, original_signal,
                   color=original_color, linestyle=original_style,
                   linewidth=line_widths[0],
                   label='Original' if col == 0 else '')

            # 绘制重建信号
            ax.plot(recon_time, recon_signal,
                   color=recon_color, linestyle=recon_style,
                   linewidth=line_widths[1],
                   label=recon_label)

            # 添加子图标签 (a), (b), (c), (d) 等
            label_x = 0.02
            label_y = 0.95
            label_text = f'({chr(97 + col + row*4)})'
            ax.text(label_x, label_y, label_text,
                   transform=ax.transAxes,
                   fontsize=11, fontweight='bold',
                   verticalalignment='top')

            # 设置子图标题（仅第一行）
            if row == 0:
                title = signal_descriptions.get(signal_name, signal_name)
                ax.set_title(title, fontsize=10, pad=10)

            # 设置轴标签（仅最外侧）
            if row == 1:  # 最后一行显示x轴标签
                ax.set_xlabel('Time (ns)', fontsize=10)

            if col == 0:  # 第一列显示y轴标签
                if row == 0:
                    ax.set_ylabel('Wiener\nField (arb. units)', fontsize=10)
                else:
                    ax.set_ylabel('LM\nField (arb. units)', fontsize=10)

            # 设置坐标轴范围，添加5%边距
            x_min = min(original_time.min(), recon_time.min())
            x_max = max(original_time.max(), recon_time.max())
            x_range = x_max - x_min
            if x_range > 0:
                x_margin = 0.05 * x_range
                ax.set_xlim(x_min - x_margin, x_max + x_margin)
            else:
                ax.set_xlim(x_min - 0.1, x_max + 0.1)

            # 根据列设置特定的y轴范围
            # ae (列0): ±0.012, bf (列1): -0.002到0.012, cg (列2): -0.002到0.012, dh (列3): ±0.003
            y_limits = {
                0: (-0.012, 0.012),      # 第一列: ae
                1: (-0.002, 0.012),      # 第二列: bf
                2: (-0.002, 0.012),      # 第三列: cg
                3: (-0.002, 0.0027)       # 第四列: dh
            }
            y_min, y_max = y_limits[col]
            ax.set_ylim(y_min, y_max)

            # 添加细网格
            ax.grid(True, linestyle=':', alpha=0.4, linewidth=0.5)

    # 添加全局图例（放置在图形底部）
    # 创建自定义图例句柄
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=original_color, linestyle=original_style,
               linewidth=line_widths[0], label='Original'),
        Line2D([0], [0], color=wiener_color, linestyle=wiener_style,
               linewidth=line_widths[1], label='Wiener'),
        Line2D([0], [0], color=lm_color, linestyle=lm_style,
               linewidth=line_widths[2], label='LM')
    ]

    # 在图形底部添加图例，调整位置避免遮挡横坐标
    fig.legend(handles=legend_elements, loc='lower center',
               ncol=3, fontsize=9, frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, 0.08))

    # 调整子图间距
    plt.subplots_adjust(wspace=0.15, hspace=0.25)

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

    return fig, axes

def calculate_performance_metrics(signals_data):
    """计算重建性能指标"""
    metrics = {}

    for signal_name, data in signals_data.items():
        original = data['original_signal']
        wiener = data['wiener_recon']
        lm = data['lm_recon']

        # 确保长度一致（插值到相同网格）
        from scipy.interpolate import interp1d

        # 插值Wiener结果到原始时间网格
        if len(data['wiener_time']) != len(original):
            f_wiener = interp1d(data['wiener_time'], wiener,
                               bounds_error=False, fill_value='extrapolate')
            wiener_interp = f_wiener(data['original_time'])
        else:
            wiener_interp = wiener

        # LM结果已经在原始时间网格上

        # 计算均方根误差（RMSE）
        rmse_wiener = np.sqrt(np.mean((original - wiener_interp)**2))
        rmse_lm = np.sqrt(np.mean((original - lm)**2))

        # 计算归一化均方根误差（NRMSE）
        range_original = original.max() - original.min()
        if range_original > 0:
            nrmse_wiener = rmse_wiener / range_original
            nrmse_lm = rmse_lm / range_original
        else:
            nrmse_wiener = nrmse_lm = 0.0

        # 计算相关系数
        corr_wiener = np.corrcoef(original, wiener_interp)[0, 1]
        corr_lm = np.corrcoef(original, lm)[0, 1]

        metrics[signal_name] = {
            'rmse_wiener': rmse_wiener,
            'rmse_lm': rmse_lm,
            'nrmse_wiener': nrmse_wiener,
            'nrmse_lm': nrmse_lm,
            'corr_wiener': corr_wiener,
            'corr_lm': corr_lm,
        }

    return metrics

def print_metrics_table(metrics, output_file=None):
    """打印性能指标表格"""
    print("\n" + "="*70)
    print("重建性能指标")
    print("="*70)
    print(f"{'信号类型':<15} {'方法':<10} {'RMSE':<12} {'NRMSE':<12} {'相关系数':<12}")
    print("-"*70)

    for signal_name, metric in metrics.items():
        print(f"{signal_name:<15} {'Wiener':<10} {metric['rmse_wiener']:<12.4f} "
              f"{metric['nrmse_wiener']:<12.4f} {metric['corr_wiener']:<12.4f}")
        print(f"{'':<15} {'LM':<10} {metric['rmse_lm']:<12.4f} "
              f"{metric['nrmse_lm']:<12.4f} {metric['corr_lm']:<12.4f}")
        print("-"*70)

    # 保存为文本文件
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("信号重建性能指标\n")
            f.write("="*70 + "\n")
            f.write(f"{'信号类型':<15} {'方法':<10} {'RMSE':<12} {'NRMSE':<12} {'相关系数':<12}\n")
            f.write("-"*70 + "\n")

            for signal_name, metric in metrics.items():
                f.write(f"{signal_name:<15} {'Wiener':<10} {metric['rmse_wiener']:<12.4f} "
                       f"{metric['nrmse_wiener']:<12.4f} {metric['corr_wiener']:<12.4f}\n")
                f.write(f"{'':<15} {'LM':<10} {metric['rmse_lm']:<12.4f} "
                       f"{metric['nrmse_lm']:<12.4f} {metric['corr_lm']:<12.4f}\n")
                f.write("-"*70 + "\n")
        print(f"指标表格已保存到: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='绘制2×4 panel信号重建对比图')
    parser.add_argument('data_dir', nargs='?', default='result/different_signals',
                       help='数据目录（默认：result/different_signals）')
    parser.add_argument('--output', '-o', default='result/comparison.png',
                       help='输出文件路径（默认：result/comparison.png）')
    parser.add_argument('--latex', action='store_true',
                       help='使用LaTeX渲染文本（需要系统安装LaTeX）')
    parser.add_argument('--signals', '-s', nargs=4,
                       default=['sine', 'step', 'double_peak', 'complex'],
                       help='信号名称列表，按顺序对应各列（默认：sine step double_peak complex）')

    args = parser.parse_args()

    # 设置绘图样式
    use_latex = set_plot_style(args.latex)

    # 检查数据目录
    if not os.path.exists(args.data_dir):
        print(f"错误：数据目录 {args.data_dir} 不存在")
        print("请先运行 generate_single_signal_data.py 生成数据")
        sys.exit(1)

    # 检查所有信号数据是否存在
    missing_signals = []
    for signal_name in args.signals:
        data_file = os.path.join(args.data_dir, f'signal_{signal_name}.npz')
        if not os.path.exists(data_file):
            missing_signals.append(signal_name)

    if missing_signals:
        print(f"错误：以下信号数据文件缺失: {missing_signals}")
        print(f"请确保已在 {args.data_dir} 目录中生成所有信号数据")
        sys.exit(1)

    # 创建对比图
    print("创建2×4 panel对比图...")
    try:
        fig, axes = create_panel_figure(args.signals, args.data_dir, args.output, use_latex)
    except Exception as e:
        print(f"创建图形时出错: {e}")
        sys.exit(1)

    # 计算性能指标
    print("\n计算性能指标...")
    try:
        # 加载数据
        signals_data = {}
        for name in args.signals:
            signals_data[name] = load_signal_data(name, args.data_dir)

        metrics = calculate_performance_metrics(signals_data)

        # 保存指标表格
        metrics_file = args.output.replace('.png', '_metrics.txt').replace('.jpg', '_metrics.txt')
        print_metrics_table(metrics, metrics_file)

    except Exception as e:
        print(f"计算性能指标时出错: {e}")
        # 继续执行，不影响主图生成

    print(f"\n处理完成！")
    print(f"- 对比图: {args.output}")
    print(f"- PDF版本: {args.output.replace('.png', '.pdf').replace('.jpg', '.pdf')}")

    # 显示提示信息
    print("\n提示：")
    print("1. 如需调整图形，可以修改脚本中的绘图参数")
    print("2. 对于出版物，建议使用PDF格式（矢量图）")
    print("3. 可以尝试 --latex 选项获得更好的排版效果（需要系统安装LaTeX）")

if __name__ == '__main__':
    main()