'''
创建协议类，定义各种量子传感协议
'''

import numpy as np
from qutip import *
from src.qubit import TransmonQubit
from src.signal import Signal, CompositeSignal
from src.pulse import Pulse, CompositePulse
from src.pulse import create_pulse, create_ramsey_pulse, create_diff_echo_pulse, create_cpmg_pulse, create_cryoscope_pulse
from src.analysis import Analysis
class Protocal:
    def __init__(self, type = 0, **kwargs):
        '''
        初始化Protocal对象
        :param type: 协议类型，0-拉比振荡测量，1-Ramsey测量，2-差分回波测量，3-CPMG作为带通滤波器测量 4-瞬态磁场测量协议 5-cryoscope协议
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
            't_global': np.linspace(-50, 300, 700),  # 全局时间列表（ns），用于演化和测量
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
                # 创建测试信号
                Phi = Signal(type = 2, t_list = np.linspace(0, 250, 500), amplitude = 0.001, frequency = 0.01, rise = 10, fall = 10, center = 100, noise_level = 0.000)
                #Phi = Signal(type = 4, t_list = np.linspace(0, 250, 500), amplitude = 0.01, frequency = 0.004, rise = 10, fall = 10, center = 100, noise_level = 0.0001)
                Phi.plot()
                omega_d = qubit.frequency
                #detuning = qubit.frequency * 0.01  # 失谐频率（GHz）
                qubit.qubit_in_mag(Phi, frame = 1, omega_d = omega_d)
                t_rabi = np.linspace(0, 20, 40)
                tau_list = np.linspace(0, 250, 500)
                p_e_list = []
                psi_e = basis(qubit.n_levels, 1)
                for tau in tau_list:
                    # 执行IQ调制，生成Ramsey脉冲序列
                    control_pulse =create_ramsey_pulse(t_rabi, tau, omega_d = omega_d, phase1 = 0.0, phase2 = 0.0)
                    # 对其时间轴，时间轴0点为第一个pi/2脉冲的结尾
                    control_pulse.t_list = control_pulse.t_list - t_rabi[-1]

                    H = QobjEvo(qubit.H_list, tlist = qubit.mag_signal.t_list, order = 1) + QobjEvo(control_pulse.hamiltonian, tlist = control_pulse.t_list, order = 1)
                    result = mesolve(H, qubit.state, self.params['t_global'], [], e_ops = [psi_e * psi_e.dag()])
                    p_e_list.append(result.expect[0][-1])
                return Phi, tau_list, p_e_list
                varphi = np.unwrap(np.arccos(2 * np.array(p_e_list) - 1))
                B = np.zeros_like(varphi)
                for i, tau in enumerate(tau_list):
                    if 0 < i < len(varphi) - 1:
                        B[i] = (varphi[i + 1] - varphi[i - 1]) / (2 * (tau_list[1] - tau_list[0]))
                    elif i == 0:
                        B[i] = (varphi[i + 1] - varphi[i]) / (tau_list[1] - tau_list[0])
                    else:
                        B[i] = (varphi[i] - varphi[i - 1]) / (tau_list[1] - tau_list[0])
                return tau_list, p_e_list, varphi, B
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
                # 创建测试信号
                k = 5
                t_list = np.linspace(0, 100, 200)
                Phi = Signal(type = 3, t_list = t_list, amplitude = 0.01, frequency = 0.01, rise = 10, fall = 10, center = 50, noise_level = 0.0001)
                Phi_list = [Phi.copy() for _ in range(2*k)]
                composite_phi = CompositeSignal(Phi_list)
                composite_phi.plot()
                omega_d = qubit.frequency
                qubit.qubit_in_mag(composite_phi, frame = 1, omega_d = omega_d)
                t_rabi = np.linspace(0, 10, 20)
                t_int = (t_rabi[-1] - t_rabi[0]) * 0.5
                t_rep = Phi.t_list[-1] - Phi.t_list[0]
                tau_list = Phi.t_list.copy()
                psi_e = basis(qubit.n_levels, 1)
                p_e_list = []
                t_global = np.linspace(-10, 1010, 2020)
                for tau in tau_list:
                    control_pulse = create_diff_echo_pulse(t_rabi, tau, t_int, t_rep, k = k, omega_d = omega_d)
                    control_pulse.t_list = control_pulse.t_list - t_rabi[-1]
                    H = QobjEvo(qubit.H_list, tlist = qubit.mag_signal.t_list, order = 1) + QobjEvo(control_pulse.hamiltonian, tlist = control_pulse.t_list, order = 1)
                    result = mesolve(H, qubit.state, t_global, [], e_ops = [psi_e * psi_e.dag()])
                    p_e = result.expect[0][-1]
                    p_e_list.append(p_e)
                return Phi, tau_list, p_e_list, k, t_int
            case 3: # CPMG测量
                pass
            case 4: # 瞬态磁场测量协议
                # 创建测试信号
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
            case 5: # cryoscope协议
                # 测试信号
                Phi = Signal(type = 2, t_list = np.linspace(0, 80, 160), amplitude = 0.01)
                Phi.plot()
                trunc_list = Phi.t_list[140:20:-1]  # 从后往前截取，模拟不同的时间延迟
                t_rabi = np.linspace(0, 10, 20)
                dt = trunc_list[0] - trunc_list[1]
                tau = 100
                psi_e = basis(qubit.n_levels, 1)
                p_e_list = [[], []]  # 分别存储I和Q分量的结果
                # 从后往前截取Phi，从而避免每次重构Phi
                for trunc in trunc_list:
                    Phi.truncate(0, trunc)
                    qubit.qubit_in_mag(Phi, frame = 1, omega_d = qubit.frequency)
                    p_e_I, p_e_Q = IQ_readout(qubit, type = 3, tau = tau)
                    p_e_list[0].append(p_e_I)
                    p_e_list[1].append(p_e_Q)
                varphi_list = np.arctan2(np.array(p_e_list[1]) - 0.5, np.array(p_e_list[0]) - 0.5)
                varphi_list = varphi_list[::-1]  # 对齐时间轴，去掉最后一个点
                return trunc_list, varphi_list, Phi, p_e_list


            
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


        

class Calibration:
    def __init__(self, qubit, type = 0, **kwargs):
        '''
        初始化Calibration对象
        :param type: 校准类型， 频率f_01标定：0-Ramsey标定  频率-磁通f(Phi)标定：1-ramsey标定 2-瞬态磁场测量协议校准 3-cryoscope协议校准
        :param kwargs: 校准参数，例如扫描范围，扫描步长等
        '''
        self.qubit = qubit
        self.type = type
        self.params = kwargs


    def calibrate(self):
        match self.type:
            case 0: # 频率f_01标定：Ramsey标定
                pass
            case 1: # 频率-磁通f(Phi)标定：Ramsey标定
                pass
            case 2: # 瞬态磁场测量协议校准
                pass
            case 3: # 相位-磁通varphi(h)标定：cryoscope协议校准
                qubit = self.qubit
                h_list = np.linspace(-0.03, 0.03, 21)
                varphi_list = []
                tau = 100
                def signals(h):
                    # 设计一个信号，确保在两个脉冲之间有一个持续时间为tau的平坦区域，且该区域内磁场强度为h
                    t_list = np.linspace(0, tau + 20, 240)
                    signal = np.zeros_like(t_list)
                    signal[(t_list >= 10) & (t_list <= tau + 10)] = h
                    return signal
                # z轴磁场应该在两个脉冲之间
                Phi = Signal(type = 8, t_list = np.linspace(0, tau + 20, 240), signal = signals(0))
                for h in h_list:
                    Phi.update_signal(signal = signals(h))
                    qubit.qubit_in_mag(Phi, frame = 1, omega_d = qubit.frequency)
                    p_e_I, p_e_Q = IQ_readout(qubit, type = 2, tau = tau, h = h)
                    varphi = np.arctan2(p_e_Q - 0.5, p_e_I - 0.5)
                    varphi_list.append(varphi)
                varphi_list = np.unwrap(varphi_list, period = np.pi)
                return h_list, varphi_list, tau


def IQ_readout(qubit, type, **kwargs):
    '''
    根据协议类型进行IQ读出
    :param qubit: 作用的Qubit对象，要求已经包含磁场信号
    :param Phi_signal: 作用的磁场信号
    :param type: 协议类型，0-rasmey 测量 1-echo测量，2-cryoscope标定 3-cryoscope测量
    '''
    # 创建脉冲
    match type:
        case 0: # Ramsey测量
            pass
        case 1: # echo测量
            pass
        case 2: # cryoscope标定
            tau = 20 if 'tau' not in kwargs else kwargs['tau']
            h = 0.01 if 'h' not in kwargs else kwargs['h']
            t_rabi = np.linspace(0, 10, 20)
            control_pulse_I = create_ramsey_pulse(t_rabi, tau, omega_d = qubit.frequency, phase1 = np.pi/2, phase2 = 0.0)
            control_pulse_Q = create_ramsey_pulse(t_rabi, tau, omega_d = qubit.frequency, phase1 = np.pi/2, phase2 = np.pi/2)
            t_evolve = control_pulse_I.t_list
        case 3: # cryoscope测量
            tau = 20 if 'tau' not in kwargs else kwargs['tau']
            t_rabi = np.linspace(0, 10, 20) if 't_rabi' not in kwargs else kwargs['t_rabi']
            control_pulse_I = create_ramsey_pulse(t_rabi, tau, omega_d = qubit.frequency, phase1 = np.pi/2, phase2 = 0.0)
            control_pulse_Q = create_ramsey_pulse(t_rabi, tau, omega_d = qubit.frequency, phase1 = np.pi/2, phase2 = np.pi/2)
            t_evolve = control_pulse_I.t_list
    H_I = QobjEvo(control_pulse_I.hamiltonian, tlist = control_pulse_I.t_list, order = 1) + QobjEvo(qubit.H_list, tlist = qubit.mag_signal.t_list, order = 1)
    H_Q = QobjEvo(control_pulse_Q.hamiltonian, tlist = control_pulse_Q.t_list, order = 1) + QobjEvo(qubit.H_list, tlist = qubit.mag_signal.t_list, order = 1)
        
    result_I = mesolve(H_I, qubit.state, t_evolve, [], e_ops = [basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()],options={"store_states": True})
    result_Q = mesolve(H_Q, qubit.state, t_evolve, [], e_ops = [basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()],options={"store_states": True})
    index = np.searchsorted(t_evolve, 10)

    print(f"State at t={t_evolve[index]:.2f} ns: {result_I.states[index]}")
    return result_I.expect[0][-1], result_Q.expect[0][-1]