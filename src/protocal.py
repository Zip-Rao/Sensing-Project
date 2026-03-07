'''
创建协议类，定义各种量子传感协议
'''

import numpy as np
from qutip import *
from src.qubit import TransmonQubit
from src.signal import Signal
from src.pulse import Pulse, CompositePulse
from src.pulse import create_pulse, create_ramsey_pulse
from src.analysis import Analysis
class Protocal:
    def __init__(self, type = 0, **kwargs):
        '''
        初始化Protocal对象
        :param type: 协议类型，0-拉比振荡测量，1-Ramsey测量，2-自旋回波测量，3-CPMG测量 4-瞬态磁场测量协议
        :param kwargs: 协议参数，例如脉冲序列，时间间隔等
        '''
        self.type = type
        self.params = kwargs


    def initialize(self, qubit:TransmonQubit, state = 0):
        '''
        根据协议类型初始化量子态
        :param qubit: 作用的Qubit对象
        :param state: 初始能级状态（整数或Qobj）
        :param n_levels: 能级数
        '''
        if isinstance(state, int) and 0 <= state < qubit.n_levels:
            qubit.state = basis(qubit.n_levels, state)
        elif isinstance(state, Qobj) and state.dims == [[qubit.n_levels], [1]]:
            qubit.state = state.unit()

        default_params = {
            't_list': np.linspace(0, 100, 1000),  # 时间列表（ns）
            'tau_list': np.linspace(0, 100, 100),  # 时间间隔列表（ns）
            't_rabi': np.linspace(0, 40, 100),  # 拉比振荡时间列表（ns）
        } 
        # 补全参数，用默认值代替
        for key, value in default_params.items():
            if key not in self.params:
                self.params[key] = value

    def evolve(self, qubit):
        '''
        根据协议类型演化量子态
        :param qubit: 作用的Qubit对象
        '''
        match self.type:
            case 0: # 拉比振荡测量
                psi_e = basis(qubit.n_levels, 1)
                t_rabi = np.linspace(0, 40, 1000)
                H_0 = qubit.get_hamiltonian_rwa(qubit.frequency)
                H_pulse = create_pulse(qubit, 1, 1, t_rabi, qubit.frequency, phase=0.0)
                H_rabi = lambda t, args: H_0 + H_pulse(t)
                result = mesolve(H_rabi, qubit.state, t_rabi, [], e_ops = [psi_e * psi_e.dag()])
                return result

            case 1: # Ramsey测量
                psi_e = basis(qubit.n_levels, 1)
                detuning = qubit.frequency * 0.01  # 失谐频率（GHz）
                omega_d = qubit.frequency - detuning
                H_0 = qubit.get_hamiltonian_rwa(omega_d)
                H_1 = create_pulse(qubit, 1, 1, self.params['t_rabi'], omega_d, phase = 0.0, angle = np.pi/2)
                U_1 = (-1j * H_1(0) * self.params['t_rabi'][-1]).expm()
                result = []
                for tau in self.params['tau_list']:
                    U_0 = (-1j * H_0 * tau).expm()
                    psi_1 = U_1 * qubit.state
                    psi_2 = U_0 * psi_1
                    psi_3 = U_1 * psi_2
                    p_e = expect(psi_e * psi_e.dag(), psi_3)
                    result.append(p_e)
                return result
                    # 处理result，计算Ramsey fringes等
            case 2: # 自旋回波测量
                pass
            case 3: # CPMG测量
                pass
            #TODO：区分正信号和负信号
            case 4: # 瞬态磁场测量协议
                #Phi = Signal(type = 1, t_list = np.linspace(0, 200, 400), amplitude = 0.01)
                Phi = Signal(type = 4, t_list = np.linspace(0, 200, 400), amplitude = 0.01, rise = 10, fall = 10, center = 100)
                Phi.plot()
                Phi_0 = Signal(type = 1, t_list = np.linspace(0, 200, 400), amplitude = 0.0)
                t_rabi = np.linspace(0, 10, 20)
                control_pulse =create_ramsey_pulse(t_rabi, tau = 0.0, omega_d = qubit.frequency)
                scan_list, p_e = self.sliding_measrement(qubit, Phi, control_pulse)
                scan_list_base, p_e_base = self.sliding_measrement(qubit, Phi_0, control_pulse)
                analysis = Analysis()
                t_samples, kernel = analysis.get_kernel(control_pulse, qubit)
                delta_p = np.array(p_e) - np.array(p_e_base)
                return t_samples, kernel, scan_list, delta_p
            
    def single_measurement(self, qubit:TransmonQubit, Phi_signal:Signal, control_pulse:CompositePulse, t_delay):
        '''
        在给定时间延迟下，进行单次测量
        :param t_delay: 磁场信号和脉冲信号的时间延迟（ns）
        待完善：这里不考虑退相干，因此将时间窗口缩小到只包含脉冲，即以脉冲时间为时间轴
        '''
    
        # 计算Qubit在延迟后的磁场信号下的频率变化
        qubit_t = qubit.qubit_under_mag(Phi_signal)

        # 计算脉冲序列下的哈密顿量
        H_pulse = lambda t: control_pulse.get_hamiltonian(t)
        
        delta = t_delay - 0.5 * control_pulse.t_list[-1]
        # 定义总哈密顿量
        def H_total(t, args):
            if t + delta < 0:
                H_0 = qubit.get_hamiltonian() if control_pulse.frame == 0 else qubit.get_hamiltonian_rwa(control_pulse.omega_d)
            elif 0 <= t + delta < Phi_signal.t_list[-1]:
                index = np.searchsorted(Phi_signal.t_list, t + delta)
                qubit_current = qubit_t[index]
                H_0 = qubit_current.get_hamiltonian() if control_pulse.frame == 0 else qubit_current.get_hamiltonian_rwa(control_pulse.omega_d)
            else:
                H_0 = qubit.get_hamiltonian() if control_pulse.frame == 0 else qubit.get_hamiltonian_rwa(control_pulse.omega_d)
            return H_0 + H_pulse(t)
            

        # 演化量子态
        result = mesolve(H_total, qubit.state, control_pulse.t_list, [], e_ops=[])
        final_state = result.states[-1]
        p_e = expect(basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag(), final_state)
        return p_e
    
    def sliding_measrement(self, qubit:TransmonQubit, Phi_signal:Signal, control_pulse:CompositePulse):
        '''
        通过改变磁场信号和脉冲信号的时间延迟，实现滑动测量
        用于支持瞬态磁场测量协议
        :param start: 瞬态磁场信号开始时间（ns）
        '''
        delay_start = -0.5 * control_pulse.t_list[-1]  # 最小延迟，确保脉冲完全覆盖在信号上
        delay_end = Phi_signal.t_list[-1] + 0.5 * control_pulse.t_list[-1]  # 最大延迟，确保脉冲完全覆盖在信号上
        n_samples = len(Phi_signal.t_list) + len(control_pulse.t_list) - 1  # 滑动测量的样本数量，等于信号和脉冲的卷积
        scan_list = np.linspace(delay_start, delay_end, n_samples)  # 滑动扫描时间点列表（ns），其长度与B与k卷积的长度一致
        p_e = []
        for i, t_delay in enumerate(scan_list):
            p_e.append(self.single_measurement(qubit, Phi_signal, control_pulse, t_delay))
        return scan_list, p_e


        

        
