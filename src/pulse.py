'''
创建脉冲类，定义传感中的虚拟脉冲，可用于构建脉冲序列及其哈密顿量
支持旋转系和实验系的脉冲定义，支持控制旋转角度
TODO:给脉冲类加一个kernel属性，表示该脉冲的控制核函数，未来可以直接用来做卷积和反卷积
'''

import numpy as np
from qutip import *
from src.qubit import TransmonQubit
from src.signal import Signal

class Pulse:
    def __init__(self, frame = 0, omega_d = 0.0, phase = 0.0, Omega = None, is_rwa = True, qubit=None):
        '''
        初始化Pulse对象
        :param frame: 参考系，0-实验系，1-旋转系
        :param omega_d: 驱动频率（GHz）
        :param phase: 旋转轴相位（rad），旋转轴在xy平面的分量与Phase有关，在z轴分量由失谐调节
        :param Omega: Rabi 频率包络线（GHz），属于Signal类对象或常数
        :param is_rwa: 是否使用旋转波近似
        :param qubit: 作用的Qubit对象，用于获取能级数（可选，如果未提供则默认使用2能级）
        '''
        self.frame = frame
        self.omega_d = omega_d
        self.phase = phase
        self.Omega = Omega
        self.is_rwa = is_rwa
        self.qubit = qubit
        if qubit is not None:
            self.n_levels = qubit.n_levels
        else:
            self.n_levels = 2  # 默认两能级，向后兼容

    def get_Rabi_frequency(self, t):
        '''
        获取脉冲在时间t处的Rabi频率
        :param t: 时间点（ns）
        '''
        if isinstance(self.Omega, Signal):
            return self.Omega.value_at(t)
        
        elif isinstance(self.Omega, (int, float)):
            return self.Omega
    
    def get_hamiltonian(self, t):
        '''
        获取脉冲在时间t处的哈密顿量，使用湮灭算符a和产生算符a^†表示
        :param t: 时间点（ns）
        '''
        n = self.n_levels
        a = destroy(n)  # 湮灭算符
        adag = create(n)  # 产生算符
        Omega_t = self.get_Rabi_frequency(t)
        if self.frame == 0:  # 实验系
            phase_term = self.omega_d * t + self.phase
            H_pulse = Omega_t * np.cos(phase_term) * (a + adag)
            return H_pulse
        else:  # 旋转系
            if self.is_rwa:
                H_rwa = Omega_t / 2 * (a * np.exp(1j * self.phase) + adag * np.exp(-1j * self.phase))
                return H_rwa
            else:
                H_rwa = Omega_t / 2 * (a * np.exp(1j * self.phase) + adag * np.exp(-1j * self.phase))
                phase_cr = 2 * self.omega_d * t + self.phase
                H_cr = Omega_t / 2 * (a * np.exp(-1j * phase_cr) + adag * np.exp(1j * phase_cr))
                return H_rwa + H_cr
    
    def get_angle(self, qubit:TransmonQubit=None):
        '''
        获取脉冲的旋转角度，与qubit的性质有关
        :param qubit: 作用的Qubit对象（可选），如果未提供则使用初始化时传入的qubit
        '''
        if qubit is not None:
            self.qubit = qubit
            self.n_levels = qubit.n_levels

        n_levels = self.n_levels
        psi_0 = qubit.state if self.qubit is not None else basis(n_levels, 0)
        t_list = self.Omega.t_list
        dt = t_list[1] - t_list[0]
        U = qeye(n_levels)
        for i in range(len(t_list)-1):
            t_i = ( t_list[i] + t_list[i+1] ) / 2
            H_i = self.get_hamiltonian(t_i)
            U_i = (-1j * H_i * dt).expm()
            U = U_i * U
        # 计算旋转角度
        psi_final = U * psi_0  
        rho_sub = psi_final.ptrace([0]) if n_levels > 2 else psi_final * psi_final.dag()
        a = U[0, 0]
        a = np.real(a)
        angle = 2 * np.arccos(a)
        c0 = psi_final[0,0]
        c1 = psi_final[1,0]
        norm = np.sqrt(abs(c0)**2 + abs(c1)**2)
    
    # 4. 计算极角 Theta (Rabi Angle)
    # P1 = |c1|^2 (处于激发态的概率)
    # P1 = sin^2(theta/2)  =>  theta = 2 * arcsin(sqrt(P1))
        p1 = abs(c1)**2 / norm**2
        theta = 2 * np.arcsin(np.sqrt(p1))
    
    # 5. 计算方位角 Phi (Phase)
    # 相对相位：c1 / c0 = tan(theta/2) * exp(i*phi)
    # phi = arg(c1) - arg(c0)
        phi = np.angle(c1) - np.angle(c0)
    
        return theta, phi  # 返回两个角度更完整
    
    def get_angle_simple(self):
        '''
        获取脉冲的旋转角度（rad），考虑共振情况，角度等于Rabi频率积分
        '''
        num_points = 1000
        t_list = np.linspace(self.Omega.t_list[0], self.Omega.t_list[-1], num_points)  # 假设脉冲持续时间为1ns
        Omegas = [self.get_Rabi_frequency(t) for t in t_list]
        angle = np.trapezoid(Omegas, t_list)  # 计算积分
        return angle
    
    def back_pulse(self):
        '''
        背景脉冲（可加噪声）
        '''
        return np.zeros_like(self.Omega.t_list)

    def get_kernel(self, qubit):
        '''
        获取脉冲的控制核函数，表示该脉冲对测量结果的影响，未来可以直接用来做卷积和反卷积
        :param qubit: 作用的Qubit对象，用于计算控制核
        '''
        kernel = []
        t_list = self.Omega.t_list
        samples = range(0,len(t_list),1)
        t_samples = t_list[::1]
        psi_e = basis(qubit.n_levels, 1)

        if self.frame == 0:
            H_0 = qubit.get_hamiltonian(qubit.frequency)
        else:
            H_0 = qubit.get_hamiltonian_rwa(qubit.frequency)
        H_pulse = lambda t: self.get_hamiltonian(t)
        H_base = lambda t, args: H_0 + H_pulse(t)
        result_base = mesolve(H_base, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()])
        p_e_base = result_base.expect[0][-1]
        for i in samples:
            t_i = t_list[i]
            stim_pulse = Signal(
                type = 3,
                t_list = t_list,
                amplitude = 1,  # 待完善：自动调整幅度，远大于脉冲幅度
                center = t_i,
                width = 2  # 待完善：自动调整宽度，足够窄
            )
            H_stim = [sigmaz(), stim_pulse]  # 待完善，支持更多能级
            stim_area = stim_pulse.params['amplitude'] * stim_pulse.params['width'] / 4 * np.sqrt(2 * np.pi)
            H_total = lambda t, args: H_0 + H_pulse(t) + H_stim[0] * H_stim[1].value_at(t)
            result_stim = mesolve(H_total, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()])
            p_e_stim = result_stim.expect[0][-1]
            kernel.append((p_e_stim-p_e_base) / stim_area)
        self.t_samples = t_samples
        self.kernel = kernel



class CompositePulse:
    def __init__(self, pulses):
        '''
        初始化CompositePulse对象，包含多个Pulse对象
        :param pulses: Pulse对象列表
        '''
        self.pulses = pulses
        self.t_list = self.get_t_list()
        self.frame = pulses[0].frame 
        self.omega_d = pulses[0].omega_d
    def get_t_list(self):
        '''
        获取复合脉冲的时间列表，返回list形式
        '''
        t_list = []
        curr = 0.0
        for pulse in self.pulses:
            pulse_list = [t + curr for t in pulse.Omega.t_list]
            t_list.extend(pulse_list)
            if pulse_list:
                curr = pulse_list[-1]
        return t_list
    
    def get_Omega(self,t):
        '''
        获取复合脉冲的Rabi频率包络线，返回t时刻的幅值
        '''
        curr = 0.0
        for pulse in self.pulses:
            if t >= curr and t <= curr + pulse.Omega.t_list[-1]:
                return pulse.get_Rabi_frequency(t - curr)
            curr += pulse.Omega.t_list[-1]


    def get_hamiltonian(self, t):
        '''
        获取复合脉冲在时间t处的哈密顿量
        '''
        curr = 0.0
        for pulse in self.pulses:
            if t >= curr and t <= curr + pulse.Omega.t_list[-1]:
                return pulse.get_hamiltonian(t - curr)
            curr += pulse.Omega.t_list[-1]

        return Qobj(np.zeros((pulse.n_levels, pulse.n_levels)))  # 超出范围返回零哈密顿量
    def plot(self):
        '''
        绘制复合脉冲的Rabi频率包络线
        '''
        import matplotlib.pyplot as plt
        t_list = self.t_list
        Omega_values = [self.get_Omega(t) for t in t_list]
        plt.figure(figsize=(8, 4))
        plt.plot(t_list, Omega_values)
        plt.xlabel('Time (ns)')
        plt.ylabel('Rabi Frequency (GHz)')
        plt.title('Composite Pulse Rabi Frequency Envelope')
        plt.grid()
        plt.show()
    
    def get_kernel(self, qubit):
        '''
        获取复合脉冲的控制核函数，表示该复合脉冲对测量结果的影响，未来可以直接用来做卷积和反卷积
        :param qubit: 作用的Qubit对象，用于计算控制核
        '''
        kernel = []
        t_list = self.t_list
        samples = range(0,len(t_list),1)
        t_samples = t_list[::1]
        psi_e = basis(qubit.n_levels, 1)

        if self.frame == 0:
            H_0 = qubit.get_hamiltonian(qubit.frequency)
        else:
            H_0 = qubit.get_hamiltonian_rwa(qubit.frequency)
        H_pulse = lambda t: self.get_hamiltonian(t)
        H_base = lambda t, args: H_0 + H_pulse(t)
        result_base = mesolve(H_base, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()])
        p_e_base = result_base.expect[0][-1]
        for i in samples:
            t_i = t_list[i]
            stim_pulse = Signal(
                type = 3,
                t_list = t_list,
                amplitude = 1,  # 待完善：自动调整幅度，远大于脉冲幅度
                center = t_i,
                width = 2  # 待完善：自动调整宽度，足够窄
            )
            H_stim = [sigmaz(), stim_pulse]  # 待完善，支持更多能级
            stim_area = stim_pulse.params['amplitude'] * stim_pulse.params['width'] / 4 * np.sqrt(2 * np.pi)
            H_total = lambda t, args: H_0 + H_pulse(t) + H_stim[0] * H_stim[1].value_at(t)
            result_stim = mesolve(H_total, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()])
            p_e_stim = result_stim.expect[0][-1]
            kernel.append((p_e_stim-p_e_base) / stim_area)
        self.t_samples = t_samples
        self.kernel = kernel
    

#定义一些常用脉冲序列，返回含时哈密顿量

def create_pulse(qubit:TransmonQubit, frame, type, t_list, omega_d, phase, angle = None, **kwargs):
    '''
    创建Pulse对象，给定参考系，波形，驱动频率，旋转角度（可选），旋转轴和其他参数，自动计算Rabi频率包络线，默认使用RWA，未来研究BC频移可设置is_rwa=False
    :param type: 脉冲类型，与signal类一致
    :param t_list: 时间列表（ns），用于生成Rabi频率序列
    :param kwargs: 脉冲参数，可以包括幅度，频率，相位，中心位置，宽度，直流偏置等
    '''
    Omega_signal = Signal(type=type, t_list=t_list, **kwargs)
    Omega_pulse = Pulse(frame, omega_d, phase, Omega=Omega_signal, is_rwa=True, qubit=qubit)
    if qubit.frequency == omega_d:  # 共振情况
        current_angle = Omega_pulse.get_angle_simple()
    else:
        current_angle = Omega_pulse.get_angle_simple()
    # 调整幅度以实现目标旋转角度
    if angle is not None:
        Omega_signal.params['amplitude'] *= angle / current_angle
    # 重新生成脉冲对象
    kwargs['amplitude'] = Omega_signal.params['amplitude']
    signal = Signal(type=type, t_list=t_list, **kwargs)
    pulse = Pulse(frame, omega_d, phase, Omega=signal, is_rwa=True, qubit=qubit)
    H_t = lambda t: pulse.get_hamiltonian(t)
    return H_t



def create_ramsey_pulse(t_rabi, tau, omega_d=0.0):
    '''
    创建Ramsey序列的复合脉冲对象，包含两个π/2脉冲和一个等待时间tau
    :param t_rabi: π/2脉冲的时间列表（ns）
    :param tau: 等待时间（ns）
    '''
    if tau != 0.0:
        Omega_0 = Signal(
            type = 0,
            t_list = np.linspace(0, tau, 100)  # ns
        )
    

    Omega_1 = Signal(
        type = 1,
        t_list = t_rabi,
        amplitude = ( np.pi / 2.0 ) / ( t_rabi[-1] - t_rabi[0] )  # GHz
    )
    pulses = []
    pulses.append(Pulse(
        frame = 1,
        omega_d = omega_d,
        phase = np.pi / 2.0,
        Omega = Omega_1,
        is_rwa = True
    ))
    if tau != 0.0:
        pulses.append(Pulse(
            frame = 1,
            omega_d = omega_d,
            phase = 0.0,
            Omega = Omega_0,
            is_rwa = True
        ))
    pulses.append(Pulse(
        frame = 1,
        omega_d = omega_d,
        phase = 0.0,
        Omega = Omega_1,
        is_rwa = True
    ))
    composite_pulse = CompositePulse(pulses)
    return composite_pulse

