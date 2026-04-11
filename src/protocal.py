'''
创建协议类，定义各种量子传感协议
'''

import numpy as np
from qutip import *
from src.qubit import TransmonQubit
from src.signal import Signal
from src.pulse import Pulse, CompositePulse
from src.pulse import create_pulse, create_ramsey_pulse, create_echo_pulse
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
                H_rabi = QobjEvo(H_0) + QobjEvo(H_pulse, tlist = t_rabi)
                result = mesolve(H_rabi, qubit.state, t_rabi, [], e_ops = [psi_e * psi_e.dag()])
                return result

            case 1: # Ramsey测量
                Phi = Signal(type = 1, t_list = np.linspace(0, 250, 500), amplitude = 0.01)
                omega_d = qubit.frequency
                #detuning = qubit.frequency * 0.01  # 失谐频率（GHz）
                qubit.qubit_in_mag(Phi, frame = 1, omega_d = omega_d)
                t_rabi = np.linspace(0, 20, 40)
                tau_list = np.linspace(0, 200, 400)
                p_e_list = []
                psi_e = basis(qubit.n_levels, 1)
                for tau in tau_list:
                    control_pulse =create_ramsey_pulse(t_rabi, tau, omega_d = omega_d, phase1 = 0.0, phase2 = 0.0)
                    
                    H = QobjEvo(qubit.H_list, tlist = qubit.mag_signal.t_list, order = 1) + QobjEvo(control_pulse.hamiltonian, tlist = control_pulse.t_list, order = 1)
                    result = mesolve(H, qubit.state, control_pulse.t_list, [], e_ops = [psi_e * psi_e.dag()])
                    p_e_list.append(result.expect[0][-1])
                def pe(t, C, Delta, phase):
                    return 0.5 * (1 + C * np.cos(Delta * t + phase))
                import scipy.optimize as optimize
                popt, pcov = optimize.curve_fit(pe, tau_list, p_e_list, p0=[1, 2 * np.pi * 0.01, 0])
                Delta_opt = popt[1]
                return tau_list, p_e_list, popt, pcov

                # psi_e = basis(qubit.n_levels, 1)
                # detuning = qubit.frequency * 0.01  # 失谐频率（GHz）
                # omega_d = qubit.frequency - detuning
                # H_0 = qubit.get_hamiltonian_rwa(omega_d)
                # H_1 = create_pulse(qubit, 1, 1, self.params['t_rabi'], omega_d, phase = 0.0, angle = np.pi/2)
                # U_1 = (-1j * H_1(0) * self.params['t_rabi'][-1]).expm()
                # result = []
                # for tau in self.params['tau_list']:
                #     U_0 = (-1j * H_0 * tau).expm()
                #     psi_1 = U_1 * qubit.state
                #     psi_2 = U_0 * psi_1
                #     psi_3 = U_1 * psi_2
                #     p_e = expect(psi_e * psi_e.dag(), psi_3)
                #     result.append(p_e)
                # return result

                    # 处理result，计算Ramsey fringes等
            case 2: # spin echo测量
                Phi = Signal(type = 2, t_list = np.linspace(0, 250, 500), amplitude = 0.01, noise_level = 0.001)
                omega_d = qubit.frequency
                qubit.qubit_in_mag(Phi, frame = 1, omega_d = omega_d)
                t_rabi = np.linspace(0, 20, 40)
                tau_list = np.linspace(0, 200, 400)
                psi_e = basis(qubit.n_levels, 1)
                p_e_list = [[], []]  # 分别存储两种不同相位的自旋回波测量结果
                for tau in tau_list:
                    control_pulse1 = create_echo_pulse(t_rabi, tau, omega_d = omega_d)
                    control_pulse2 = create_echo_pulse(t_rabi, tau, omega_d = omega_d, phase3 = np.pi/2)
                    H_1 = QobjEvo(qubit.H_list, tlist = qubit.mag_signal.t_list, order = 1) + QobjEvo(control_pulse1.hamiltonian, tlist = control_pulse1.t_list, order = 1)
                    H_2 = QobjEvo(qubit.H_list, tlist = qubit.mag_signal.t_list, order = 1) + QobjEvo(control_pulse2.hamiltonian, tlist = control_pulse2.t_list, order = 1)
                    result1 = mesolve(H_1, qubit.state, control_pulse1.t_list, [], e_ops = [psi_e * psi_e.dag()])
                    result2 = mesolve(H_2, qubit.state, control_pulse2.t_list, [], e_ops = [psi_e * psi_e.dag()])
                    p_e_list[0].append(result1.expect[0][-1])
                    p_e_list[1].append(result2.expect[0][-1])
                return tau_list, p_e_list
            case 3: # CPMG测量
                pass
            #TODO：区分正信号和负信号
            case 4: # 瞬态磁场测量协议
                #Phi = Signal(type = 1, t_list = np.linspace(0, 200, 400), amplitude = 0.01)
                t_list = np.linspace(0,200, 400)
                Phi = Signal(type = 4, t_list = t_list, amplitude = 0.06, rise = 10, fall = 10, center = 100, noise_level = 0.0001)
                Phi.plot()
                Phi_0 = Signal(type = 1, t_list = t_list, amplitude = 0.0)
                t_rabi = np.linspace(0, 10, 20)
                print("create control pulse")
                control_pulse =create_ramsey_pulse(t_rabi, tau = 0.0, omega_d = qubit.frequency)
                print("start measurement")
                scan_list, p_e = self.sliding_measrement(qubit, Phi, control_pulse)
                
                scan_list_base, p_e_base = self.sliding_measrement(qubit, Phi_0, control_pulse)
                control_pulse.get_kernel(qubit)
                t_samples, kernel = control_pulse.t_samples, control_pulse.kernel
                #analysis = Analysis()
                #t_samples, kernel = analysis.get_kernel(control_pulse, qubit)
                delta_p = np.array(p_e) - np.array(p_e_base)
                return t_samples, kernel, scan_list, delta_p, p_e, Phi, control_pulse
            
    def single_measurement(self, qubit:TransmonQubit, Phi_signal:Signal, control_pulse:CompositePulse, t_delay, qubit_t=None, H=None, t_evole=None, index = None):
        '''
        在给定时间延迟下，进行单次测量
        :param t_delay: 磁场信号和脉冲信号的时间延迟（ns）
        :param qubit_t: 预计算的qubit时间序列，由qubit.qubit_under_mag(Phi_signal)得到，避免重复计算
        待完善：这里不考虑退相干，因此将时间窗口缩小到只包含脉冲，即以脉冲时间为时间轴
        '''

        # 计算Qubit在延迟后的磁场信号下的频率变化
        import time
        
        start = time.time()
        if qubit_t is None:
            qubit_t = qubit.qubit_under_mag(Phi_signal)

        # 计算脉冲序列下的哈密顿量
        
        delta = t_delay - 0.5 * control_pulse.t_list[-1]
        if t_evole is None:
            t_start = min(delta, 0)
            t_end = max(Phi_signal.t_list[-1], t_delay + 0.5 * control_pulse.t_list[-1])
            t_list = np.linspace(t_start, t_end, len(Phi_signal.t_list) + len(control_pulse.t_list) - 1)
        else:
            t_list = t_evole
        if index is None:
            delay_start = Phi_signal.t_list[0] - 0.5 * control_pulse.t_list[-1]  # 最小延迟，确保脉冲完全覆盖在信号上
            delay_end = Phi_signal.t_list[-1] + 0.5 * control_pulse.t_list[-1]  # 最大延迟，确保脉冲完全覆盖在信号上
            n_samples = len(Phi_signal.t_list) + len(control_pulse.t_list) - 1  # 滑动测量的样本数量，等于信号和脉冲的卷积
            scan_list = np.linspace(delay_start, delay_end, n_samples)  # 滑动扫描时间点列表（ns），其长度与B与k卷积的长度一致
            index = np.searchsorted(scan_list, t_delay)
        if H is None:
        # 定义总哈密顿量
            H_pulse = lambda t: control_pulse.get_hamiltonian_at(t)

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
        else:
            
            H_total = QobjEvo(H, tlist = t_list, order = 1)
        
        # 演化量子态
        dt_coeff = t_list[1] - t_list[0]
        options = {
            "atol" :1e-10,
            "rtol" :1e-8,
            "max_step":dt_coeff / 2,
            "nsteps":10000
        }
        #result = mesolve(H_total, qubit.state, control_pulse.t_list, [], e_ops=[], options = options)
        result = mesolve(H_total, qubit.state, t_list, [], e_ops=[], options = options)
        #result = mesolve(H_total, qubit.state, t_list, [], e_ops=[])
        final_state = result.states[-1]
        p_e = expect(basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag(), final_state)
        end = time.time()
        #print(f"Single measurement with delay {t_delay:.2f} ns took {end - start:.2f} seconds, p_e = {p_e:.4f}")
        return p_e
    
    def sliding_measrement(self, qubit:TransmonQubit, Phi_signal:Signal, control_pulse:CompositePulse):
        '''
        通过改变磁场信号和脉冲信号的时间延迟，实现滑动测量
        用于支持瞬态磁场测量协议
        :param start: 瞬态磁场信号开始时间（ns）
        '''
        import time
        start = time.time()
        delay_start = Phi_signal.t_list[0] - 0.5 * control_pulse.t_list[-1]  # 最小延迟，确保脉冲完全覆盖在信号上
        delay_end = Phi_signal.t_list[-1] + 0.5 * control_pulse.t_list[-1]  # 最大延迟，确保脉冲完全覆盖在信号上
        n_samples = len(Phi_signal.t_list) + len(control_pulse.t_list) - 1  # 滑动测量的样本数量，等于信号和脉冲的卷积
        scan_list = np.linspace(delay_start, delay_end, n_samples)  # 滑动扫描时间点列表（ns），其长度与B与k卷积的长度一致
        # 预计算qubit，哈密顿量，演化时间轴
        qubit_t = qubit.qubit_under_mag(Phi_signal)  # 预计算，所有t_delay共用同一个Phi_signal，只需计算一次
        n_levels = qubit.n_levels
        n_op = qubit.n
        def H(t_delay):
            delta = t_delay - 0.5 * control_pulse.t_list[-1]
            t_start = min(delta, 0)
            t_end = max(Phi_signal.t_list[-1], t_delay + 0.5 * control_pulse.t_list[-1])
            N = len(Phi_signal.t_list) + len(control_pulse.t_list) - 1
            t_evole = np.linspace(t_start, t_end, N)
            freq_coeffs = np.zeros(N)       #随时间变化的项主要是频率项，理论上非谐项也会变化，但后续在考虑（TODO）
            EJ_coeffs = np.zeros(N)      
            for i, t in enumerate(t_evole):
                if 0 <= t <= Phi_signal.t_list[-1]:
                    index = np.searchsorted(Phi_signal.t_list, t)
                    index = min(index, len(qubit_t) - 1)  # 确保索引不越界
                    qubit_current = qubit_t[index]
                else:
                    qubit_current = qubit

                if control_pulse.frame == 0:
                    freq_coeffs[i] = qubit_current.frequency
                else:
                    freq_coeffs[i] = qubit_current.frequency - control_pulse.omega_d
                EJ_coeffs[i] = -qubit_current.EJ + 0.25 * qubit_current.EC 
            H_list = []
            id = np.ones(N, dtype = complex)
            if control_pulse.frame == 0:
                H_list.append([qubit.anharmonicity * 0.5 * (qubit.n * qubit.n - qubit.n), id])
                H_list.append([qubit.n + 0.5 * qeye(n_levels), freq_coeffs])
                #H_list.append([qeye(n_levels), EJ_coeffs])
            else:
                H_list.append([qubit.anharmonicity * 0.5 * (qubit.n * qubit.n - qubit.n), id])
                H_list.append([qubit.n, freq_coeffs])
                #H_list.append([qeye(n_levels), EJ_coeffs])
            t_pulse = control_pulse.t_list
            for op, coeffs in control_pulse.hamiltonian:
                coeff_global = np.zeros(N, dtype=complex)
                for i, t in enumerate(t_evole):
                     t_loc = t - delta
                     if 0 <= t_loc <= t_pulse[-1]:
                        index = np.clip(np.searchsorted(t_pulse, t_loc), 0, len(t_pulse) - 1)  # 确保索引在有效范围内
                        index = min(index, len(control_pulse.hamiltonian[0][1]) - 1)  # 确保索引不越界
                        coeff_global[i] = coeffs[index]
                H_list.append([op, coeff_global])
            return H_list, t_evole


        #from joblib import Parallel, delayed
        #p_e = []
        #p_e = Parallel(n_jobs=6, backend='threading')(delayed(self.single_measurement)(qubit, Phi_signal, control_pulse, t_delay, qubit_t, H_list) for t_delay in scan_list)
        p_e = []
        for i, t_delay in enumerate(scan_list):
             H_total, t_evole = H(t_delay)
             p_e.append(self.single_measurement(qubit, Phi_signal, control_pulse, t_delay, qubit_t, H_total, t_evole))
        end = time.time()
        print(f"Sliding measurement with {len(scan_list)} samples took {end - start:.2f} seconds")
        return scan_list, p_e


        

        
