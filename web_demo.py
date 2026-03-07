#!/usr/bin/env python3
"""
量子传感仿真平台网页交互Demo
基于现有src模块创建的可交互网页界面
使用Gradio构建

注意: 本demo不修改原代码，仅使用现有功能
未实现的功能会在界面上给出提示
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import traceback

# ============================================================================
# 路径设置 - 确保能够导入src模块
# ============================================================================

# 获取当前脚本所在目录
current_dir = Path(__file__).parent.absolute()
project_root = current_dir

# 添加项目根目录和src目录到Python路径
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "src"))

print(f"项目根目录: {project_root}")
print(f"Python路径: {sys.path}")

# ============================================================================
# 检查核心依赖
# ============================================================================

required_modules = ['qutip', 'numpy', 'matplotlib']
missing_modules = []

for module in required_modules:
    try:
        __import__(module)
    except ImportError:
        missing_modules.append(module)

if missing_modules:
    print("错误: 缺少必需的Python模块:")
    for module in missing_modules:
        print(f"  - {module}")
    print("\n请安装所有依赖:")
    print("  pip install qutip numpy matplotlib scipy")
    print("\n或者使用项目中的requirements.txt:")
    print("  pip install -r requirements.txt")
    sys.exit(1)

# ============================================================================
# 导入项目模块
# ============================================================================

try:
    # 注意: 原代码中使用 "protocal" 而不是 "protocol"
    from src.qubit import TransmonQubit
    from src.signal import Signal
    from src.pulse import Pulse, CompositePulse, create_pulse, create_ramsey_pulse
    from src.protocal import Protocal  # 注意拼写: protocal 不是 protocol
    from src.analysis import Analysis
    print("✓ 成功导入所有项目模块")
except ImportError as e:
    print(f"✗ 导入错误: {e}")
    print("\n可能的原因:")
    print("1. 确保在项目根目录下运行此脚本")
    print("2. 检查src目录是否包含所有必要的.py文件")
    print("3. 原代码中使用 'protocal' 而不是 'protocol'，请确认文件名")
    print("\n必需的文件:")
    print("  - src/qubit.py")
    print("  - src/signal.py")
    print("  - src/pulse.py")
    print("  - src/protocal.py")
    print("  - src/analysis.py")
    sys.exit(1)

# ============================================================================
# 检查Gradio
# ============================================================================

try:
    import gradio as gr
    print("✓ Gradio已安装")
except ImportError:
    print("✗ Gradio未安装")
    print("请使用以下命令安装:")
    print("  pip install gradio")
    sys.exit(1)

# ============================================================================
# 仿真函数
# ============================================================================

def run_simulation(ec, ej, t1, t2, flux, n_levels, initial_state, protocol_type,
                   t_max, tau_max, t_rabi, signal_amplitude, signal_frequency,
                   signal_width, lambda_param):
    """
    运行量子传感仿真

    参数说明:
    ec: 电容能量 (GHz)
    ej: 约瑟夫森能量 (GHz)
    t1: 弛豫时间 (ns)
    t2: 退相干时间 (ns)
    flux: 外加磁通 (单位：Φ0)
    n_levels: 能级数
    initial_state: 初始能级状态 (0或1)
    protocol_type: 协议类型 (0-拉比振荡, 1-Ramsey, 4-瞬态磁场测量)
    t_max: 总仿真时间 (ns)
    tau_max: Ramsey延迟最大时间 (ns)
    t_rabi: 拉比脉冲时间 (ns)
    signal_amplitude: 信号幅度 (瞬态磁场测量用)
    signal_frequency: 信号频率 (GHz)
    signal_width: 信号宽度 (ns)
    lambda_param: 维纳反卷积参数
    """
    try:
        # ====================================================================
        # 1. 创建qubit对象
        # ====================================================================
        print(f"创建Qubit: EC={ec} GHz, EJ={ej} GHz, flux={flux}, n_levels={n_levels}")

        # 注意: 原代码中EC和EJ需要乘以2π转换为角频率
        qubit = TransmonQubit(
            EC=ec * 2 * np.pi,  # 转换为角频率
            EJ=ej * 2 * np.pi,
            T1=t1,
            T2=t2,
            flux=flux,
            state=initial_state,
            n_levels=n_levels
        )

        print(f"Qubit创建成功: 频率={qubit.frequency/(2*np.pi):.3f} GHz, "
              f"非谐性={qubit.anharmonicity/(2*np.pi):.3f} GHz")

        # ====================================================================
        # 2. 根据协议类型设置参数
        # ====================================================================
        protocol_params = {}

        if protocol_type == 0:  # 拉比振荡
            print(f"协议0: 拉比振荡测量, t_max={t_max} ns")
            t_list = np.linspace(0, t_max, 1000)
            protocol_params['t_list'] = t_list

        elif protocol_type == 1:  # Ramsey测量
            print(f"协议1: Ramsey干涉测量, tau_max={tau_max} ns, t_rabi={t_rabi} ns")
            tau_list = np.linspace(0, tau_max, 100)
            t_rabi_list = np.linspace(0, t_rabi, 100)
            protocol_params['tau_list'] = tau_list
            protocol_params['t_rabi'] = t_rabi_list

        elif protocol_type == 4:  # 瞬态磁场测量
            print(f"协议4: 瞬态磁场测量, t_max={t_max} ns")
            t_list = np.linspace(0, t_max, 400)
            protocol_params['t_list'] = t_list

        else:
            # 其他协议类型（2, 3, 5等）在原代码中未完全实现
            error_msg = f"协议类型 {protocol_type} 在原代码中未完全实现\n\n"
            error_msg += "已实现的协议类型:\n"
            error_msg += "  0: 拉比振荡测量\n"
            error_msg += "  1: Ramsey干涉测量\n"
            error_msg += "  4: 瞬态磁场测量\n"
            error_msg += "\n请查看TODO.md了解需要完善的功能"
            raise NotImplementedError(error_msg)

        # ====================================================================
        # 3. 创建协议对象并运行
        # ====================================================================
        print(f"创建协议对象: type={protocol_type}")
        protocol = Protocal(type=protocol_type, params=protocol_params)
        protocol.initialize(qubit, state=initial_state)

        print("运行协议演化...")
        results = protocol.evolve(qubit)
        print("协议演化完成")

        # ====================================================================
        # 4. 处理和分析结果
        # ====================================================================
        analysis = Analysis()

        if protocol_type == 0:  # 拉比振荡
            # results应该是mesolve结果对象
            p_e = analysis.get_expectation_values(results, e_ops_index=0)
            t_list = protocol_params['t_list']

            # 创建图表
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(t_list, p_e, 'b-', linewidth=2, label='激发态概率')
            ax.set_xlabel('时间 (ns)', fontsize=12)
            ax.set_ylabel('激发态概率', fontsize=12)
            ax.set_title('拉比振荡测量', fontsize=14)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=11)
            ax.set_ylim(0, 1)

            # 显示qubit信息
            qubit_info = f"""
            **Qubit参数**:
            - 频率: {qubit.frequency/(2*np.pi):.3f} GHz
            - 非谐性: {qubit.anharmonicity/(2*np.pi):.3f} GHz
            - T1: {t1:,} ns
            - T2: {t2:,} ns
            - 磁通: {flux} Φ₀
            - 能级数: {n_levels}
            - 初始态: |{initial_state}⟩

            **协议参数**:
            - 协议类型: 0 (拉比振荡)
            - 总时间: {t_max} ns
            - 采样点数: 1000
            """

        elif protocol_type == 1:  # Ramsey测量
            # results应该是概率列表
            tau_list = protocol_params['tau_list']

            # 创建图表
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(tau_list, results, 'g-', linewidth=2, label='激发态概率')
            ax.set_xlabel('延迟时间 τ (ns)', fontsize=12)
            ax.set_ylabel('激发态概率', fontsize=12)
            ax.set_title('Ramsey干涉条纹', fontsize=14)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=11)
            ax.set_ylim(0, 1)

            # 显示qubit信息
            qubit_info = f"""
            **Qubit参数**:
            - 频率: {qubit.frequency/(2*np.pi):.3f} GHz
            - 非谐性: {qubit.anharmonicity/(2*np.pi):.3f} GHz
            - T1: {t1:,} ns
            - T2: {t2:,} ns
            - 磁通: {flux} Φ₀
            - 能级数: {n_levels}
            - 初始态: |{initial_state}⟩

            **协议参数**:
            - 协议类型: 1 (Ramsey干涉)
            - 最大延迟: {tau_max} ns
            - 拉比脉冲时间: {t_rabi} ns
            """

        elif protocol_type == 4:  # 瞬态磁场测量
            # results应该是元组: (t_samples, kernel, scan_list, delta_p)
            t_samples, kernel, scan_list, delta_p = results

            # 执行维纳反卷积
            dt = scan_list[1] - scan_list[0]
            B_lists, B_recon = analysis.wiener_deconvolution(
                delta_p, kernel, dt, lambdas=lambda_param
            )

            # 创建多个子图
            fig, axes = plt.subplots(2, 2, figsize=(14, 10))

            # 子图1: 核函数
            axes[0, 0].plot(t_samples, kernel, 'b-', linewidth=2)
            axes[0, 0].set_xlabel('时间 (ns)', fontsize=11)
            axes[0, 0].set_ylabel('核函数值', fontsize=11)
            axes[0, 0].set_title('控制脉冲核函数', fontsize=12)
            axes[0, 0].grid(True, alpha=0.3)

            # 子图2: 测量概率变化
            axes[0, 1].plot(scan_list, delta_p, 'g-', linewidth=2)
            axes[0, 1].set_xlabel('扫描时间 (ns)', fontsize=11)
            axes[0, 1].set_ylabel('ΔP (概率变化)', fontsize=11)
            axes[0, 1].set_title('测量到的概率变化', fontsize=12)
            axes[0, 1].grid(True, alpha=0.3)

            # 子图3: 重建的磁场信号
            axes[1, 0].plot(B_lists, B_recon, 'r-', linewidth=2)
            axes[1, 0].set_xlabel('时间 (ns)', fontsize=11)
            axes[1, 0].set_ylabel('重建磁场 (任意单位)', fontsize=11)
            axes[1, 0].set_title('维纳反卷积重建的瞬态磁场', fontsize=12)
            axes[1, 0].grid(True, alpha=0.3)

            # 子图4: 激发态概率原始测量
            # 注意: delta_p是变化量，需要加上基线概率(0.5)
            axes[1, 1].plot(scan_list, np.array(delta_p) + 0.5, 'm-', linewidth=2)
            axes[1, 1].set_xlabel('扫描时间 (ns)', fontsize=11)
            axes[1, 1].set_ylabel('激发态概率', fontsize=11)
            axes[1, 1].set_title('原始激发态概率测量', fontsize=12)
            axes[1, 1].grid(True, alpha=0.3)
            axes[1, 1].set_ylim(0, 1)

            plt.tight_layout()

            # 显示重建信息
            qubit_info = f"""
            **Qubit参数**:
            - 频率: {qubit.frequency/(2*np.pi):.3f} GHz
            - 非谐性: {qubit.anharmonicity/(2*np.pi):.3f} GHz
            - 磁通: {flux} Φ₀
            - 能级数: {n_levels}
            - 初始态: |{initial_state}⟩

            **重建参数**:
            - 维纳参数 λ: {lambda_param}
            - 重建信号长度: {len(B_recon)} 点
            - 最大重建值: {np.max(np.abs(B_recon)):.3e}
            - 信号幅度: {signal_amplitude}
            - 信号频率: {signal_frequency} GHz
            - 信号宽度: {signal_width} ns
            """

        plt.tight_layout()
        return fig, qubit_info

    except NotImplementedError as e:
        # 协议未实现的特定错误
        fig, ax = plt.subplots(figsize=(10, 6))
        error_msg = str(e)
        ax.text(0.5, 0.5, error_msg,
               ha='center', va='center', fontsize=12, wrap=True)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        return fig, f"**错误**: {error_msg}"

    except Exception as e:
        # 其他错误
        error_trace = traceback.format_exc()
        print(f"仿真错误: {e}")
        print(f"详细错误信息:\n{error_trace}")

        fig, ax = plt.subplots(figsize=(10, 6))
        error_msg = f"仿真过程中发生错误:\n\n{str(e)}\n\n请检查参数是否合理。"
        ax.text(0.5, 0.5, error_msg,
               ha='center', va='center', fontsize=12, wrap=True)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')

        error_info = f"""
        **仿真错误**:
        {str(e)}

        **可能的原因**:
        1. 参数值超出物理合理范围
        2. 原代码中的某些限制
        3. 计算过程中出现数值问题

        **建议**:
        - 检查输入参数
        - 查看控制台输出获取更多信息
        - 参考Simulation.ipynb中的示例参数
        """
        return fig, error_info


# ============================================================================
# UI更新函数
# ============================================================================

def update_ui(protocol_type):
    """根据选择的协议类型更新UI显示"""
    # 将protocol_type转换为整数
    try:
        protocol_type = int(protocol_type)
    except:
        protocol_type = 0

    # 显示/隐藏相关参数控件
    visible_t_max = protocol_type in [0, 4]
    visible_tau_max = protocol_type == 1
    visible_t_rabi = protocol_type in [0, 1]
    visible_signal_params = protocol_type == 4
    visible_lambda_param = protocol_type == 4

    return (
        gr.update(visible=visible_t_max),
        gr.update(visible=visible_tau_max),
        gr.update(visible=visible_t_rabi),
        gr.update(visible=visible_signal_params),  # signal_amplitude
        gr.update(visible=visible_signal_params),  # signal_frequency
        gr.update(visible=visible_signal_params),  # signal_width
        gr.update(visible=visible_lambda_param)
    )


# ============================================================================
# 创建Gradio界面
# ============================================================================

def create_interface():
    """创建Gradio界面"""
    with gr.Blocks(title="量子传感仿真平台", theme=gr.themes.Soft()) as demo:
        gr.Markdown("# 🧪 量子传感仿真平台 - 交互式Demo")
        gr.Markdown("""
        本交互式演示基于超导量子比特传感仿真平台。
        调整Qubit参数和协议参数，运行仿真并查看结果。

        **注意**: 本demo使用原代码功能，不修改任何源代码。
        某些功能可能受限，详情请查看项目根目录下的TODO.md。
        """)

        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("## 📊 Qubit参数")

                with gr.Group():
                    ec = gr.Slider(
                        minimum=0.01, maximum=2.0, value=0.2,
                        step=0.01, label="电容能量 EC (GHz)",
                        info="Transmon电容能量，典型值: 0.1-0.5 GHz"
                    )
                    ej = gr.Slider(
                        minimum=1.0, maximum=50.0, value=10.0,
                        step=0.1, label="约瑟夫森能量 EJ (GHz)",
                        info="Transmon约瑟夫森能量，典型值: 5-20 GHz"
                    )
                    t1 = gr.Slider(
                        minimum=1e3, maximum=1e6, value=100e3,
                        step=1e3, label="弛豫时间 T1 (ns)",
                        info="能量弛豫时间，典型值: 10⁴-10⁶ ns"
                    )
                    t2 = gr.Slider(
                        minimum=1e3, maximum=1e6, value=50e3,
                        step=1e3, label="退相干时间 T2 (ns)",
                        info="相位退相干时间，通常小于T1"
                    )
                    flux = gr.Slider(
                        minimum=0.0, maximum=0.5, value=0.0,
                        step=0.01, label="外加磁通 Φ/Φ₀",
                        info="归一化磁通，范围: 0-0.5"
                    )
                    n_levels = gr.Dropdown(
                        choices=[2, 3, 4, 5], value=2,
                        label="能级数",
                        info="计算中考虑的能级数量"
                    )
                    initial_state = gr.Dropdown(
                        choices=[0, 1], value=0,
                        label="初始能级状态",
                        info="Qubit的初始状态"
                    )

                gr.Markdown("## ⚙️ 协议参数")

                protocol_type = gr.Dropdown(
                    choices=[
                        ("0", "0 - 拉比振荡测量"),
                        ("1", "1 - Ramsey干涉测量"),
                        ("4", "4 - 瞬态磁场测量")
                    ],
                    value="0",
                    label="协议类型",
                    info="选择要执行的量子传感协议"
                )

                t_max = gr.Slider(
                    minimum=10, maximum=500, value=100,
                    step=10, label="总仿真时间 (ns)",
                    visible=True,
                    info="仿真的总时间长度"
                )
                tau_max = gr.Slider(
                    minimum=10, maximum=500, value=100,
                    step=10, label="最大延迟时间 τ (ns)",
                    visible=False,
                    info="Ramsey协议中两个脉冲之间的最大延迟"
                )
                t_rabi = gr.Slider(
                    minimum=1, maximum=100, value=40,
                    step=1, label="拉比脉冲时间 (ns)",
                    visible=True,
                    info="π/2脉冲的持续时间"
                )

                with gr.Group(visible=False) as signal_group:
                    gr.Markdown("### 瞬态磁场参数")
                    signal_amplitude = gr.Slider(
                        minimum=0.001, maximum=0.1, value=0.01,
                        step=0.001, label="信号幅度",
                        info="瞬态磁场信号的幅度"
                    )
                    signal_frequency = gr.Slider(
                        minimum=0.1, maximum=10.0, value=1.0,
                        step=0.1, label="信号频率 (GHz)",
                        info="瞬态磁场信号的频率"
                    )
                    signal_width = gr.Slider(
                        minimum=1, maximum=50, value=10,
                        step=1, label="信号宽度 (ns)",
                        info="瞬态磁场信号的宽度"
                    )

                lambda_param = gr.Slider(
                    minimum=0.0, maximum=10.0, value=1.0,
                    step=0.1, label="维纳反卷积参数 λ",
                    visible=False,
                    info="维纳反卷积中的正则化参数"
                )

                run_btn = gr.Button("🚀 运行仿真", variant="primary", size="lg")

                gr.Markdown("""
                ### 💡 使用提示
                1. 从Simulation.ipynb中选择已验证的参数组合
                2. 复杂计算可能需要几秒钟时间
                3. 查看控制台输出了解详细进度
                """)

            with gr.Column(scale=2):
                gr.Markdown("## 📈 仿真结果")
                plot_output = gr.Plot(label="结果图表", show_label=True)
                info_output = gr.Markdown(
                    label="仿真信息",
                    value="调整参数后点击'运行仿真'查看结果..."
                )

        # ====================================================================
        # 事件绑定
        # ====================================================================
        protocol_type.change(
            update_ui,
            inputs=[protocol_type],
            outputs=[t_max, tau_max, t_rabi, signal_amplitude,
                    signal_frequency, signal_width, lambda_param]
        )

        run_btn.click(
            run_simulation,
            inputs=[
                ec, ej, t1, t2, flux, n_levels, initial_state,
                protocol_type, t_max, tau_max, t_rabi,
                signal_amplitude, signal_frequency, signal_width,
                lambda_param
            ],
            outputs=[plot_output, info_output]
        )

        gr.Markdown("""
        ## 📋 协议说明

        ### 协议0: 拉比振荡测量
        - **目的**: 测量Qubit在共振微波驱动下的激发态概率振荡
        - **物理**: 施加共振驱动场，观测|0⟩和|1⟩态之间的Rabi振荡
        - **输出**: 激发态概率随时间的变化

        ### 协议1: Ramsey干涉测量
        - **目的**: 测量Qubit的频率，检测微小频率变化
        - **物理**: 两个π/2脉冲之间自由演化，观测干涉条纹
        - **输出**: 激发态概率随延迟时间的变化

        ### 协议4: 瞬态磁场测量
        - **目的**: 重建随时间变化的瞬态磁场信号
        - **物理**: 滑动测量 + 维纳反卷积
        - **输出**: 核函数、测量结果、重建的磁场信号

        ## ⚠️ 已知限制

        1. **未实现的协议**: 协议2(自旋回波)和3(CPMG)在原代码中未完全实现
        2. **计算时间**: 某些计算可能需要较长时间
        3. **参数范围**: 某些参数组合可能导致数值问题

        详细问题列表请查看项目根目录下的TODO.md文件。

        ## 🔧 技术信息

        - **仿真平台**: 基于QuTiP的量子光学工具箱
        - **Qubit模型**: Transmon qubit (非线性谐振子)
        - **单位制**: 自然单位制 (ħ=1)，频率单位GHz，时间单位ns
        - **代码状态**: 使用原始代码，未经修改

        ---

        *最后更新: 2026-03-07*
        *量子传感仿真平台 Web Demo*
        """)

    return demo


# ============================================================================
# 主程序
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("量子传感仿真平台 - Web交互式演示")
    print("=" * 70)

    # 检查是否在项目根目录下运行
    src_dir = Path(__file__).parent / "src"
    required_files = ["qubit.py", "signal.py", "pulse.py", "protocal.py", "analysis.py"]
    missing_files = []

    for file in required_files:
        if not (src_dir / file).exists():
            missing_files.append(file)

    if missing_files:
        print("错误: 找不到必需的源文件:")
        for file in missing_files:
            print(f"  - {file}")
        print(f"\n请确保在项目根目录下运行此脚本。")
        print(f"当前目录: {Path.cwd()}")
        print(f"期望的src目录: {src_dir}")
        sys.exit(1)

    print("✓ 所有必需文件都存在")
    print("\n启动Web服务器...")
    print("请稍候，初始化可能需要几秒钟")
    print("\n服务器启动后，请在浏览器中访问: http://localhost:7860")
    print("按Ctrl+C停止服务器")
    print("-" * 70)

    # 创建并启动界面
    try:
        demo = create_interface()
        demo.launch(
            server_name="0.0.0.0",
            server_port=7861,
            share=True,  # 生成公网链接
            show_error=True,
            quiet=False,  # 显示启动信息
            show_api=False
        )
    except Exception as e:
        print(f"启动失败: {e}")
        print("\n可能的原因:")
        print("1. 端口7860已被占用")
        print("2. 网络权限问题")
        print("3. Gradio安装问题")
        print("\n尝试以下解决方案:")
        print("1. 关闭占用7860端口的程序")
        print("2. 使用其他端口: 修改server_port参数")
        print("3. 检查防火墙设置")
        sys.exit(1)