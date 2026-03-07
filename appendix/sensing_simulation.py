"""
传感核函数仿真脚本
基于现有仿真平台 (src模块) 实现传感核函数提取
功能对应于示例代码中的 TransmonSensingSim 类
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from qutip import *

# 添加src目录到Python路径
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./src'))

from src.qubit import TransmonQubit
from src.signal import Signal
from src.pulse import Pulse, create_pulse


class SensingSimulation:
    """
    传感仿真类，利用现有平台模块实现传感核函数提取
    """
    def __init__(self, qubit, gamma=2*np.pi*1.0):
        """
        初始化传感仿真对象

        Parameters:
        -----------
        qubit : TransmonQubit
            量子比特对象
        gamma : float
            旋磁比/磁场灵敏度 (rad/ns/T)，默认值 2π * 1.0
        """
        self.qubit = qubit
        self.gamma = gamma

        # 获取量子比特频率 (GHz)
        self.w_01 = qubit.frequency  # 单位: GHz

        # 定义两能级系统的Pauli算符（使用qubit的能级数）
        self.n_levels = qubit.n_levels
        if self.n_levels < 2:
            raise ValueError("Qubit must have at least 2 levels")

        # 使用qubit的能级数创建Pauli算符
        self.sx = self._pauli_x()
        self.sy = self._pauli_y()
        self.sz = self._pauli_z()
        self.id = qeye(self.n_levels)

        # 初始态 (基态 |0><0|)
        self.rho0 = basis(self.n_levels, 0) * basis(self.n_levels, 0).dag()
        self.psi_e_proj = basis(self.n_levels, 1) * basis(self.n_levels, 1).dag()  # 用于测量 Pe

    def _pauli_x(self):
        """创建X Pauli算符（适用于n_levels能级系统）"""
        a = destroy(self.n_levels)
        return a + a.dag()

    def _pauli_y(self):
        """创建Y Pauli算符（适用于n_levels能级系统）"""
        a = destroy(self.n_levels)
        return -1j * (a - a.dag())

    def _pauli_z(self):
        """创建Z Pauli算符（适用于n_levels能级系统）"""
        n = num(self.n_levels)
        return 2*n - self.id

    def generate_control_pulse(self, t_list, width, center, shape='gaussian', distortion=False):
        """
        生成控制脉冲包络 Omega(t)，使用现有平台的Signal类

        Parameters:
        -----------
        t_list : array_like
            时间序列 (ns)
        width : float
            脉冲宽度 (ns)
        center : float
            脉冲中心时间 (ns)
        shape : str
            脉冲形状，'gaussian' 或 'square'
        distortion : bool
            是否添加畸变（模拟带宽限制）

        Returns:
        --------
        pulse_signal : Signal
            控制脉冲信号对象
        """
        if shape == 'gaussian':
            # 使用Signal类的高斯脉冲类型 (type=3)
            pulse_signal = Signal(
                type=3,  # 高斯脉冲
                t_list=t_list,
                amplitude=1.0,  # 初始幅度，后续会归一化
                center=center,
                width=width
            )
            # 获取脉冲包络
            pulse_envelope = pulse_signal.signal

        elif shape == 'square':
            # 方波脉冲：使用恒定信号加窗函数
            # 先创建零信号，然后手动设置方波
            pulse_envelope = np.zeros_like(t_list)
            # 方波区域
            pulse_envelope[(t_list > center - width/2) & (t_list < center + width/2)] = 1.0
            # 创建恒定信号并替换其信号值
            pulse_signal = Signal(
                type=1,  # 恒定信号
                t_list=t_list,
                amplitude=1.0
            )
            pulse_signal.signal = pulse_envelope
        else:
            raise ValueError(f"Unsupported pulse shape: {shape}")

        # 归一化幅度，使其积分为 Pi (构成一个 Pi 脉冲)
        # 注意：这里只是粗略归一化，精确需要 Rabi 校准
        dt = t_list[1] - t_list[0]
        pulse_area = np.sum(pulse_envelope) * dt

        if pulse_area > 0:
            scale_factor = np.pi / pulse_area
        else:
            scale_factor = 1.0

        pulse_envelope = pulse_envelope * scale_factor

        # 模拟畸变 (Distortion)：使用高斯滤波模拟有限带宽
        if distortion:
            # sigma_filter 越大，低通滤波越强，波形越"圆"
            pulse_envelope = gaussian_filter1d(pulse_envelope, sigma=5)
            # 重新归一化以保持相同的脉冲面积
            pulse_area_new = np.sum(pulse_envelope) * dt
            if pulse_area_new > 0:
                pulse_envelope = pulse_envelope * (np.pi / pulse_area_new)

        # 更新信号值
        pulse_signal.signal = pulse_envelope

        return pulse_signal

    def generate_stimulus_pulse(self, t_list, t_center, t_width, strength):
        """
        生成极短的高斯刺激信号 B_stim(t)，用于探测核函数

        Parameters:
        -----------
        t_list : array_like
            时间序列 (ns)
        t_center : float
            刺激信号中心时间 (ns)
        t_width : float
            刺激信号宽度 (ns, FWHM)
        strength : float
            刺激信号强度

        Returns:
        --------
        stim_signal : Signal
            刺激信号对象
        """
        # 使用Signal类的高斯脉冲类型 (type=3)
        stim_signal = Signal(
            type=3,  # 高斯脉冲
            t_list=t_list,
            amplitude=strength,
            center=t_center,
            width=t_width  # Signal类中width参数对应sigma*4
        )

        return stim_signal

    def run_density_matrix_sim(self, t_list, control_signal, stim_signal):
        """
        执行密度矩阵主方程求解

        Hamiltonian: H = 0.5 * Omega(t) * sx + 0.5 * gamma * B_stim(t) * sz
        假设驱动与Qubit共振 (Delta=0)，在旋转参考系下。

        Parameters:
        -----------
        t_list : array_like
            时间序列 (ns)
        control_signal : Signal
            控制脉冲信号对象
        stim_signal : Signal
            刺激信号对象

        Returns:
        --------
        p_e : float
            最终时刻的激发态布居数 Pe
        """
        # 将信号转换为数组形式
        if isinstance(control_signal, Signal):
            omega_t = control_signal.signal
        else:
            omega_t = control_signal

        if isinstance(stim_signal, Signal):
            stim_t = stim_signal.signal
        else:
            stim_t = stim_signal

        # 构建哈密顿量项
        # 控制项 (X轴旋转): 0.5 * Omega(t) * sx
        H_control = [0.5 * self.sx, omega_t]

        # 传感/刺激项 (Z轴失谐/旋转): 0.5 * gamma * B_stim(t) * sz
        H_stim = [0.5 * self.sz, self.gamma * stim_t]

        # 总哈密顿量
        H = [H_control, H_stim]

        # 求解 Lindblad 主方程 (这里暂设无耗散 c_ops=[])
        result = mesolve(H, self.rho0, t_list, c_ops=[],
                         options=Options(nsteps=5000, atol=1e-10, rtol=1e-10))

        # 返回最终时刻的激发态布居数 Pe
        rho_final = result.states[-1]
        p_e = (rho_final * self.psi_e_proj).tr().real

        return p_e

    def extract_kernel(self, t_list, control_signal, stim_strength=0.1, stim_width=1.0):
        """
        扫描提取传感核函数 k(t)

        原理: delta_p(tau) = P_e(with_stim at tau) - P_e(no_stim)

        Parameters:
        -----------
        t_list : array_like
            时间序列 (ns)
        control_signal : Signal
            控制脉冲信号对象
        stim_strength : float
            刺激信号强度
        stim_width : float
            刺激信号宽度 (ns, FWHM)

        Returns:
        --------
        scan_times : ndarray
            扫描时间点
        kernel : ndarray
            核函数值
        p_e_baseline : float
            无刺激时的基准激发态概率
        """
        print("正在计算基准演化 (无刺激)...")

        # 1. 计算基准值 (无刺激信号)
        zeros_signal = Signal(type=0, t_list=t_list)  # 零信号
        p_e_baseline = self.run_density_matrix_sim(t_list, control_signal, zeros_signal)

        kernel = []
        scan_indices = []

        # 降采样扫描，避免计算量过大 (每隔几个点算一次)
        step = max(1, len(t_list) // 250)  # 大约250个扫描点
        scan_indices = range(0, len(t_list), step)
        scan_times = t_list[scan_indices]

        print(f"开始扫描提取核函数，共 {len(scan_times)} 个时间点...")

        for i, idx in enumerate(scan_indices):
            if i % 10 == 0:
                print(f"  进度: {i+1}/{len(scan_indices)}")

            tau = t_list[idx]

            # 生成在 tau 时刻的微小刺激
            stim_pulse = self.generate_stimulus_pulse(t_list, tau, stim_width, stim_strength)

            # 计算有刺激时的结果
            p_e_stim = self.run_density_matrix_sim(t_list, control_signal, stim_pulse)

            # 计算响应差值 (即核函数的值)
            # 根据线性响应理论，delta_p ~ k(tau) * Area_stim
            delta_p = p_e_stim - p_e_baseline
            kernel.append(delta_p)

        return scan_times, np.array(kernel), p_e_baseline


def run_simulation():
    """
    执行完整的传感仿真流程
    """
    # ==========================================
    # 1. 创建量子比特对象
    # ==========================================
    print("创建量子比特对象...")
    qubit = TransmonQubit(
        EC=0.2,      # 电容能量 (GHz)
        EJ=10.0,     # 约瑟夫森能量 (GHz)
        T1=100.0e3,  # 100 μs
        T2=50.0e3,   # 50 μs
        flux=0.0,    # 外加磁通
        state=0,     # 初始态 |0>
        n_levels=2   # 两能级系统
    )

    # ==========================================
    # 2. 实例化传感仿真对象
    # ==========================================
    print("初始化传感仿真对象...")
    sim = SensingSimulation(qubit, gamma=2*np.pi*1.0)  # gamma = 2π * 1.0 rad/ns/T

    # ==========================================
    # 3. 时间轴设置
    # ==========================================
    t_total = 100  # ns
    steps = 500
    t_list = np.linspace(0, t_total, steps)

    # ==========================================
    # 4. 定义 Pi 脉冲参数
    # ==========================================
    pulse_center = 50.0  # ns
    pulse_width = 20.0   # ns

    # --- 场景 A: 理想高斯脉冲 ---
    print("\n=== 场景 A: 理想高斯控制脉冲 ===")
    ctrl_ideal = sim.generate_control_pulse(t_list, pulse_width, pulse_center,
                                           shape='gaussian', distortion=False)
    t_scan_ideal, kernel_ideal, pe_base_ideal = sim.extract_kernel(
        t_list, ctrl_ideal, stim_strength=0.1, stim_width=1.0)

    # --- 场景 B: 畸变/有噪声的脉冲 ---
    print("\n=== 场景 B: 畸变控制脉冲 (模拟带宽限制) ===")
    ctrl_distorted = sim.generate_control_pulse(t_list, pulse_width, pulse_center,
                                               shape='gaussian', distortion=True)

    # 重新校准幅度 (畸变会改变脉冲面积，为了公平对比，简单归一化使其面积一致)
    # 在实际实验中，这对应于重新校准 Rabi 频率
    dt = t_list[1] - t_list[0]
    ideal_area = np.sum(ctrl_ideal.signal) * dt
    distorted_area = np.sum(ctrl_distorted.signal) * dt

    if distorted_area > 0:
        scale_factor = ideal_area / distorted_area
        ctrl_distorted.signal = ctrl_distorted.signal * scale_factor

    t_scan_dist, kernel_dist, pe_base_dist = sim.extract_kernel(
        t_list, ctrl_distorted, stim_strength=0.1, stim_width=1.0)

    # ==========================================
    # 5. 结果可视化
    # ==========================================
    print("\n生成可视化结果...")
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), sharex='col')

    # --- 左列：理想脉冲 ---
    # 1. 脉冲波形
    ax[0, 0].plot(t_list, ctrl_ideal.signal / (2*np.pi), 'b-',
                  label='Ideal Gaussian Control $\Omega(t)$')
    ax[0, 0].set_ylabel('Amplitude (GHz)')
    ax[0, 0].set_title('(a) Ideal Control Pulse')
    ax[0, 0].legend()
    ax[0, 0].grid(True, alpha=0.3)

    # 2. 传感核函数
    ax[1, 0].plot(t_scan_ideal, kernel_ideal, 'k-', lw=2,
                  label='Sensing Kernel $k(t)$')
    ax[1, 0].fill_between(t_scan_ideal, kernel_ideal, color='gray', alpha=0.3)
    ax[1, 0].set_ylabel('Kernel Amplitude (arb.)')
    ax[1, 0].set_xlabel('Time (ns)')
    ax[1, 0].set_title('(c) Kernel Response (Ideal)')
    ax[1, 0].legend()
    ax[1, 0].grid(True, alpha=0.3)

    # --- 右列：畸变脉冲 ---
    # 1. 脉冲波形
    ax[0, 1].plot(t_list, ctrl_ideal.signal / (2*np.pi), 'b--', alpha=0.3,
                  label='Reference')
    ax[0, 1].plot(t_list, ctrl_distorted.signal / (2*np.pi), 'r-',
                  label='Distorted Control $\Omega_{dist}(t)$')
    ax[0, 1].set_title('(b) Distorted Control Pulse')
    ax[0, 1].legend()
    ax[0, 1].grid(True, alpha=0.3)

    # 2. 传感核函数
    ax[1, 1].plot(t_scan_ideal, kernel_ideal, 'k--', alpha=0.3,
                  label='Ideal Ref')
    ax[1, 1].plot(t_scan_dist, kernel_dist, 'r-', lw=2,
                  label='Distorted Kernel $k_{dist}(t)$')
    ax[1, 1].fill_between(t_scan_dist, kernel_dist, color='red', alpha=0.1)
    ax[1, 1].set_xlabel('Time (ns)')
    ax[1, 1].set_title('(d) Kernel Response (Distorted)')
    ax[1, 1].legend()
    ax[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()

    # 保存图像
    output_file = "sensing_kernel_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"结果已保存至: {output_file}")

    plt.show()

    # ==========================================
    # 6. 打印结果摘要
    # ==========================================
    print("\n" + "="*50)
    print("仿真结果摘要")
    print("="*50)
    print(f"量子比特频率: {qubit.frequency:.3f} GHz")
    print(f"无刺激基准概率 (理想脉冲): {pe_base_ideal:.6f}")
    print(f"无刺激基准概率 (畸变脉冲): {pe_base_dist:.6f}")
    print(f"核函数最大值 (理想脉冲): {np.max(np.abs(kernel_ideal)):.6f}")
    print(f"核函数最大值 (畸变脉冲): {np.max(np.abs(kernel_dist)):.6f}")
    print("="*50)


if __name__ == "__main__":
    run_simulation()