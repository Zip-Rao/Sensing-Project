"""
量子传感协议演示入口点
提供多种演示选项，展示图片中所示的协议实现
"""

import sys
import os

# 添加src目录到Python路径
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./src'))

def main():
    """主函数"""
    print("="*70)
    print("量子传感协议演示平台")
    print("="*70)
    print("基于现有仿真平台实现图片中所示的完整协议流程")
    print("")
    print("可用的演示选项:")
    print("  1. 滑动测量演示 (Sliding Measurement)")
    print("  2. Wiener反卷积演示 (Wiener Deconvolution)")
    print("  3. 瞬态磁场测量完整协议 (Transient Field Protocol)")
    print("  4. 集成演示 (所有功能)")
    print("  5. 退出")
    print("")

    try:
        choice = input("请选择演示选项 (1-5): ").strip()

        if choice == '1':
            run_sliding_measurement_demo()
        elif choice == '2':
            run_wiener_deconvolution_demo()
        elif choice == '3':
            run_transient_field_protocol_demo()
        elif choice == '4':
            run_integrated_demo()
        elif choice == '5':
            print("退出演示")
            return
        else:
            print("无效选择，请重新运行")
            return

    except KeyboardInterrupt:
        print("\n用户中断")
    except Exception as e:
        print(f"\n运行出错: {e}")
        import traceback
        traceback.print_exc()


def run_sliding_measurement_demo():
    """运行滑动测量演示"""
    print("\n" + "="*70)
    print("滑动测量演示")
    print("="*70)

    try:
        from sliding_measurement import run_sliding_measurement_demo as run_demo
        run_demo()
    except ImportError as e:
        print(f"导入错误: {e}")
        print("请确保 sliding_measurement.py 文件存在")


def run_wiener_deconvolution_demo():
    """运行Wiener反卷积演示"""
    print("\n" + "="*70)
    print("Wiener反卷积演示")
    print("="*70)

    try:
        from deconvolution import run_deconvolution_demo as run_demo
        run_demo()
    except ImportError as e:
        print(f"导入错误: {e}")
        print("请确保 deconvolution.py 文件存在")


def run_transient_field_protocol_demo():
    """运行瞬态磁场测量协议演示"""
    print("\n" + "="*70)
    print("瞬态磁场测量协议演示")
    print("="*70)

    try:
        from appendix.transient_field_protocol import run_protocol_demo as run_demo
        run_demo()
    except ImportError as e:
        print(f"导入错误: {e}")
        print("请确保 transient_field_protocol.py 文件存在")


def run_integrated_demo():
    """运行集成演示"""
    print("\n" + "="*70)
    print("集成演示 - 完整量子传感流程")
    print("="*70)

    try:
        # 导入所需模块
        from src.qubit import TransmonQubit
        from src.signal import Signal
        from src.pulse import Pulse, CompositePulse
        from src.analysis import Analysis

        print("1. 检查模块导入... [✓]")

        # 创建量子比特
        qubit = TransmonQubit(
            EC=0.2, EJ=10.0, T1=100.0e3, T2=50.0e3,
            flux=0.0, state=0, n_levels=2
        )
        print(f"2. 创建量子比特 [✓] 频率: {qubit.frequency:.3f} GHz")

        # 创建控制脉冲
        t_pulse = np.linspace(0, 20, 100)
        signal = Signal(type=1, t_list=t_pulse, amplitude=0.1)
        pulse = Pulse(frame=1, omega_d=0.0, phase=0.0,
                     Omega=signal, is_rwa=True, qubit=qubit)
        control_pulse = CompositePulse([pulse])
        print("3. 创建控制脉冲 [✓]")

        # 创建磁场信号
        t_magnetic = np.linspace(-50, 150, 1000)
        magnetic_signal = Signal(
            type=3, t_list=t_magnetic,
            amplitude=1.0, center=50.0, width=8.0
        )
        print("4. 创建磁场信号 [✓]")

        # 导入并运行滑动测量
        try:
            from sliding_measurement import SlidingMeasurement
            sliding_meas = SlidingMeasurement(
                qubit=qubit,
                control_pulse=control_pulse,
                magnetic_signal=magnetic_signal,
                sensitivity=2*np.pi*1.0
            )
            print("5. 创建滑动测量对象 [✓]")

            # 运行简化演示（避免长时间计算）
            print("\n运行简化演示...")
            print("注意: 完整滑动扫描需要较长时间")
            print("这里运行一个小规模演示")

            # 创建简化Ramsey序列
            sliding_meas.create_ramsey_pulse_sequence(pulse_duration=10)

            # 运行少量延迟点的扫描
            delays, pe_ideal, pe_measured = sliding_meas.run_sliding_scan(
                delay_start=-10, delay_end=30, n_delays=21,
                noise_level=0.02, verbose=True
            )
            print(f"6. 滑动测量完成 [✓] 点数: {len(delays)}")

            # 获取核函数
            t_kernel, kernel = sliding_meas.get_kernel()
            print(f"7. 提取核函数 [✓] 点数: {len(kernel)}")

            # 运行Wiener反卷积
            try:
                from deconvolution import WienerDeconvolution
                deconv = WienerDeconvolution()

                # 简化反卷积
                dt = delays[1] - delays[0]
                kernel_interp = kernel  # 简化处理
                if len(kernel) != len(delays):
                    from scipy.interpolate import interp1d
                    interp_func = interp1d(t_kernel, kernel, kind='linear',
                                          bounds_error=False, fill_value=0)
                    kernel_interp = interp_func(delays)

                B_rec = deconv.wiener_deconvolution(
                    pe_measured, kernel_interp, 0.1, dt
                )
                print("8. Wiener反卷积完成 [✓]")

                # 简单可视化
                import matplotlib.pyplot as plt
                plt.figure(figsize=(10, 6))

                plt.subplot(2, 2, 1)
                t_pulse_full = np.array(control_pulse.t_list)
                Omega_vals = [control_pulse.get_Omega(t) for t in t_pulse_full]
                plt.plot(t_pulse_full, Omega_vals, 'b-')
                plt.title('Control Pulse')
                plt.xlabel('Time (ns)')
                plt.ylabel('Rabi Freq (GHz)')
                plt.grid(True, alpha=0.3)

                plt.subplot(2, 2, 2)
                B_true = [magnetic_signal.value_at(t) for t in delays]
                plt.plot(delays, B_true, 'r-')
                plt.title('Magnetic Field')
                plt.xlabel('Time (ns)')
                plt.ylabel('B(t)')
                plt.grid(True, alpha=0.3)

                plt.subplot(2, 2, 3)
                plt.plot(delays, pe_measured, 'b.', markersize=3)
                plt.title('Measured P_e(τ)')
                plt.xlabel('Delay τ (ns)')
                plt.ylabel('P_e')
                plt.grid(True, alpha=0.3)

                plt.subplot(2, 2, 4)
                # 归一化
                if np.max(np.abs(B_true)) > 0:
                    B_true_norm = np.array(B_true) / np.max(np.abs(B_true))
                else:
                    B_true_norm = np.array(B_true)

                if np.max(np.abs(B_rec)) > 0:
                    B_rec_norm = B_rec / np.max(np.abs(B_rec))
                else:
                    B_rec_norm = B_rec

                plt.plot(delays, B_true_norm, 'k--', label='True', alpha=0.6)
                plt.plot(delays, B_rec_norm, 'r-', label='Reconstructed')
                plt.title('Wiener Deconvolution')
                plt.xlabel('Time (ns)')
                plt.ylabel('Normalized')
                plt.legend()
                plt.grid(True, alpha=0.3)

                plt.tight_layout()
                plt.savefig("integrated_demo_result.png", dpi=150, bbox_inches='tight')
                plt.show()

                print("\n集成演示完成!")
                print(f"结果已保存至: integrated_demo_result.png")
                print(f"滑动测量点数: {len(delays)}")
                print(f"核函数点数: {len(kernel)}")

            except ImportError as e:
                print(f"Wiener反卷积模块导入失败: {e}")
                print("继续进行其他演示...")

        except ImportError as e:
            print(f"滑动测量模块导入失败: {e}")
            print("请确保所有模块文件存在")

    except ImportError as e:
        print(f"基本模块导入失败: {e}")
        print("请确保src目录中的模块可用")
    except Exception as e:
        print(f"运行出错: {e}")
        import traceback
        traceback.print_exc()


# 添加必要的导入（在函数外部）
import numpy as np

if __name__ == "__main__":
    main()