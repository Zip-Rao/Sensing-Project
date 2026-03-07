'''
创建分析类，返回各种数值分析结果，以及图表
包括对mesolve结果的分析
'''

import numpy as np
from qutip import *
from src.qubit import TransmonQubit
from src.signal import Signal
from src.pulse import Pulse,CompositePulse
class Analysis:
    def __init__(self):
        pass
    
    def get_expectation_values(self, result, e_ops_index):
        '''
        从mesolve结果中提取期望值，返回该期望值随时间变化的数组
        :param result: mesolve结果对象
        :param e_ops_index: 期望值操作符的索引
        '''
        if hasattr(result, 'expect'):
            if isinstance(result.expect, list) and len(result.expect) > e_ops_index:
                return np.array(result.expect[e_ops_index])
            else:
                return np.array(result.expect)


    def get_population(self, result, level):
        '''
        从mesolve结果中提取特定能级的占据概率，返回该概率随时间变化的数组
        :param result: mesolve结果对象
        :param level: 能级索引
        '''
        if hasattr(result, 'state'):
            populations = []
            for state in result.states:
                proj = basis(state.dims[0][0], level) * basis(state.dims[0][0], level).dag()
                pop = expect(proj, state)
                populations.append(pop)
            return np.array(populations)
        
    def get_kernel(self, control_pulse : CompositePulse, qubit:TransmonQubit):
        '''
        获取脉冲的控制核函数
        :param control_pulse: 脉冲对象
        '''
        kernel = []
        t_list = control_pulse.t_list
        samples = range(0,len(t_list),1)
        t_samples = t_list[::1]
        psi_e = basis(qubit.n_levels, 1)

        if control_pulse.frame == 0:
            H_0 = qubit.get_hamiltonian(qubit.frequency)
        else:
            H_0 = qubit.get_hamiltonian_rwa(qubit.frequency)
        H_pulse = lambda t: control_pulse.get_hamiltonian(t)
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
        return t_samples,np.array(kernel)

    def wiener_deconvolution(self, delta_p, kernel, dt, lambdas):
        '''
        Wiener 反卷积
        :param delta_p: 测量的概率变化 (signal_len)
        :param kernel: 控制核 (kernel_len)
        '''
        N_del = len(delta_p)
        N_ker = len(kernel)
        N = N_del - N_ker + 1  # 原始信号的长度
        
        N_fft = N_del
        
        # 补零到完整卷积长度
        p_pad = np.zeros(N_fft)
        p_pad[:N_del] = delta_p
        
        k_pad = np.zeros(N_fft)
        k_pad[:N_ker] = kernel
        
        Y = np.fft.fft(p_pad)
        H = np.fft.fft(k_pad)
        
        H_conj = np.conj(H)
        G = H_conj / (np.abs(H)**2 + lambdas**2)
        
        X_w = Y * G / dt
        x_rec = np.real(np.fft.ifft(X_w))
        
        x_final = x_rec[:N]
        x_lists = np.linspace(0, (N-1)*dt, N)  # 原始信号的时间轴
        return x_lists, x_final