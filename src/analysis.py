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
            # TODO:刺激信号的t_list与脉冲信号一致，是否会导致其带宽过大，是否要做截断
            stim_pulse = Signal(
                type = 3,
                t_list = t_list,
                amplitude = 0.0215,  # 待完善：自动调整幅度，远大于脉冲幅度
                center = t_i,
                width = 3  # 待完善：自动调整宽度，足够窄
            )
            stim_area = np.trapezoid(stim_pulse.signal, stim_pulse.t_list)
            qubit_t = qubit.qubit_under_mag(stim_pulse)
            # H_stim = [0.5 * sigmaz(), stim_pulse]  # 待完善，支持更多能级
            # H_total = lambda t, args: H_0 + H_pulse(t) + H_stim[0] * H_stim[1].value_at(t)
            def H_total(t, args):
                index = np.searchsorted(stim_pulse.t_list, t)
                index = min(index, len(stim_pulse.t_list)-1)  # 确保索引不越界
                qubit_current = qubit_t[index]
                H_0_current = qubit_current.get_hamiltonian(qubit.frequency) if control_pulse.frame == 0 else qubit_current.get_hamiltonian_rwa(qubit.frequency)
                return H_0_current + H_pulse(t)
            result_stim = mesolve(H_total, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()])
            p_e_stim = result_stim.expect[0][-1]
            kernel.append((p_e_stim-p_e_base) / stim_area)  # TODO：调整系数，考虑刺激信号的幅度和宽度
        return t_samples,np.array(kernel)

    def get_volterra_kernel(self, control_pulse : CompositePulse, qubit:TransmonQubit):
        '''
        获取二阶Volterra核函数，考虑非线性响应
        '''
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

    def hammerstein_wiener_deconvolution(self, qubit, delta_p, kernel, dt, lambdas):
        '''
        Hammerstein-Wiener 反卷积，考虑非线性响应
        块结构为：磁场B->响应函数omega(B)->相位phi->测量结果delta_p
        :param delta_p: 测量的概率变化 (signal_len)
        :param kernel: 控制核 (kernel_len)
        '''
        omega_lists, omega = self.wiener_deconvolution(delta_p, kernel, dt, lambdas)  # 先进行线性反卷积，得到频率响应
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 5))
        plt.plot(omega_lists, omega, label='Deconvolved Frequency Response')
        plt.xlabel('Time (ns)')
        plt.ylabel('Frequency Shift (GHz)')
        plt.legend()
        plt.show()
        B = (1 / np.pi) * np.arccos((omega + qubit.frequency + qubit.EC) ** 2 / (8 * qubit.EC * qubit.EJ) )   # 根据Transmon的频率-磁场关系，得到磁场响应
        B_list = omega_lists  # 磁场的时间轴与频率响应的时间轴相同
        return B_list, B
    
    def numerical_inverse(self, qubit, control_pulse, p_meas, t_meas, basis_type = 'bspline', n_basis=10, lambdas=1.0, max_iter=100, tol=1e-6):
        '''
        基于全密度矩阵模拟的数值反演算法，考虑非线性响应和退相干
        通过优化算法调整输入磁场信号，使得模拟的测量结果与实际测量结果p_meas尽可能接近
        :param p_meas: 实际测量的概率变化 (signal_len)
        :param t_meas: 实际测量的时间轴 (signal_len)
        :param basis_type: 用于表示输入磁场信号的基函数类型，支持'bspline','fourier','legendre'等
        :param n_basis: 基函数的数量
        :param lambdas: 正则化参数，控制解的平滑程度
        :param max_iter: 最大迭代次数
        :param tol: 收敛容忍度
        '''

        

# 下面定义一些辅助函数

def generate_basis_functions(basis_type, n_basis, t_min, t_max):
    '''
    生成基函数
    :param basis_type: 基函数类型，支持'bspline','fourier','legendre'等
    :param n_basis: 基函数数量
    :param t_min: 时间轴最小值
    :param t_max: 时间轴最大值
    '''
    if basis_type == 'bspline':
        from scipy.interpolate import BSpline
        knots = np.linspace(t_min, t_max, n_basis - 2)  # degree=3的B样条需要n_basis-2个内部节点
        knots = np.r_[[t_min] * 3, knots, [t_max] * 3]  # 添加边界节点
        basis_functions = []
        for i in range(n_basis):
            coeffs = np.zeros(n_basis)
            coeffs[i] = 1.0
            basis_functions.append(BSpline(knots, coeffs, 3))
        return basis_functions
    elif basis_type == 'fourier':
        T = t_max - t_min
        basis_functions = []
        basis_functions.append(lambda t: np.ones_like(t) / np.sqrt(T))  # 常数项
        n_max = (n_basis - 1) // 2
        for n in range(1, n_max + 1):
            basis_functions.append(lambda t, n=n: np.sqrt(2 / T) * np.sin(2 * np.pi * n * (t - t_min) / T))
            basis_functions.append(lambda t, n=n: np.sqrt(2 / T) * np.cos(2 * np.pi * n * (t - t_min) / T))
        return basis_functions[:n_basis]  # 只返回前n_basis个函数
    elif basis_type == 'legendre':
        from scipy.special import legendre
        basis_functions = [lambda t, n=n: legendre(n)((2 * (t - t_min) / (t_max - t_min)) - 1) for n in range(n_basis)]
        return basis_functions
    else:
        raise ValueError("Unsupported basis type")

def basis_function_decomposition(signal, t_array, basis_functions, n_basis):
    '''
    将信号分解到基函数上，得到基函数系数
    :param signal: 待分解的信号 (signal_len)
    :param t_array: 信号的时间轴 (signal_len)
    :param basis_functions: 基函数列表
    '''
    from scipy.integrate import simpson
    from scipy.interpolate import interp1d
    f = interp1d(t_array, signal, kind='cubic', fill_value="extrapolate")
    b = np.zeros(n_basis)
    for i in range(n_basis):
        basis_i = basis_functions[i](t_array)
        integrand = f(t_array) * basis_i
        b[i] = simpson(integrand, t_array)
    return b

def D(n):
    '''
    二阶差分矩阵，用于正则化
    '''
    D = np.zeros((n-2, n))
    for i in range(n-2):
        D[i, i] = 1
        D[i, i+1] = -2
        D[i, i+2] = 1
    return D

def forward_simulation(qubit, control_pulse, b, t_list):
    '''
    正向模拟，根据输入的基函数系数b生成磁场信号，并模拟延迟t时刻施加脉冲信号的测量结果
    '''
    B = np.zeros_like(t_list)
    for k, phi_k in enumerate(b):
        B += b[k] * phi_k(t_list)
    result = []
    qubit_t = qubit.qubit_under_mag(B)
    def H_total(t, t_i, args):
        index = np.searchsorted(t_list, t)
        index = min(index, len(t_list)-1)  # 确保索引不越界
        qubit_current = qubit_t[index]
        H_0_current = qubit_current.get_hamiltonian(qubit.frequency) if control_pulse.frame == 0 else qubit_current.get_hamiltonian_rwa(qubit.frequency)
        if t_i - 0.5 * control_pulse.t_list[-1] <= t <= t_i + 0.5 * control_pulse.t_list[-1]:
            return H_0_current + control_pulse.get_hamiltonian(t - t_i + 0.5 * control_pulse.t_list[-1])
        else:
            return H_0_current
        
    for t_i in t_list:
        H = lambda t, args: H_total(t, t_i, args)
        result.append(mesolve(H, qubit.state, t_list, [], e_ops=[basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()]))
    return result

def compute_jacobian(qubit, control_pulse, b, n_basis, t_lists):
    '''
    伴随方法计算雅可比矩阵
    '''
    N = len(t_lists)
    M = n_basis
    J = np.zeros((N, M))
    dim = qubit.state.shape[0]
    
    B = np.zeros_like(t_lists)
    for k, phi_k in enumerate(b):
        B += b[k] * phi_k(t_lists)
    qubit_t = qubit.qubit_under_mag(B)
    sensitivity = [qubit_t[n].calculate_sensitivity() for n in range(len(t_lists))]
    def H_total(t, t_i, args):
        index = np.searchsorted(t_lists, t)
        index = min(index, len(t_lists)-1)  # 确保索引不越界
        qubit_current = qubit_t[index]
        H_0_current = qubit_current.get_hamiltonian(qubit.frequency) if control_pulse.frame == 0 else qubit_current.get_hamiltonian_rwa(qubit.frequency)
        if t_i - 0.5 * control_pulse.t_list[-1] <= t <= t_i + 0.5 * control_pulse.t_list[-1]:
            return H_0_current + control_pulse.get_hamiltonian(t - t_i + 0.5 * control_pulse.t_list[-1])
        else:    
            return H_0_current
    # 前向传播
    result = forward_simulation(qubit, control_pulse, b, t_lists)
    states = result.states
    p_sim = np.array([res.expect[0] for res in result])  # 模拟的测量结果

    for i in range(N):
        t_i = t_lists[i]
        from scipy.integrate import solve_ivp

        def adjoint(s, mu, *args):
            t_curr = t_i - s
            H = H_total(t_curr, t_i, args)
            rhs = 1j * (H @ mu - mu @ H)
            for k in range(len(qubit.c_ops)):
                c_k = qubit.c_ops[k]
                rhs += c_k.dag() @ mu @ c_k - 0.5 * (c_k.dag() @ c_k @ mu + mu @ c_k.dag() @ c_k)
            return rhs
        
        mu_0 = basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()  # 期望值操作符
        s_max = t_i
        s_list = t_i - t_lists[:i+1][::-1]
        sol = solve_ivp(
                adjoint,
                t_span=(0, s_max),
                y0=mu_0,
                method='RK45',
                t_eval=s_list,
                rtol=1e-9,
                atol=1e-11,
                max_step=1e-10
            )
        
        N_s = len(sol.t)
        lambda_t = np.zeros((i+1, dim, dim), dtype=complex)
        for n in range(N_s):
            j = i - n
            lambda_t[j] = sol.y[:, n].reshape((dim, dim))

        trace_values = np.zeros(i + 1, dtype=complex)
        for n in range(i + 1):
            rho_n = states[n]       # ρ(t_n), 前向轨迹
            lam_n = lambda_t[n]    # λ(t_n), 伴随轨迹
            G_mat = qubit.n 
            trace_values[n] = np.trace(lam_n @ (G_mat @ rho_n - rho_n @ G_mat))

        g_values = -1j * sensitivity[:i+1] * trace_values
        t_int = t_lists[:i+1]  # 积分的时间轴
        for k in range(M):
            integrand_k = g_values * phi_k(t_lists[:i+1])
            integral = np.trapezoid(integrand_k, t_int)
            J[i, k] = np.real(integral)  # 取实部，虚部应为数值零
    return J




    
def levenberg_marquardt(qubit, p_meas, t_meas, control_pulse, basis_functions, n_basis, reg, b_init = None, max_iter = 50, tol = 1e-6, mu_init = 1e-3):
    '''
    Levenberg-Marquardt优化算法，用于最小化模拟测量结果与实际测量结果之间的差异
    :param qubit: 量子比特对象
    :param p_meas: 实际测量的概率变化 (signal_len)
    :param t_meas: 实际测量的时间轴 (signal_len)
    :param control_pulse: 控制脉冲对象，用于模拟测量结果
    :param basis_functions: 用于表示输入磁场信号的基函数列表
    :param b_init: 初始的基函数系数，可以是随机的或基于先验知识的
    :param max_iter: 最大迭代次数
    :param tol: 收敛容忍度
    :param mu_init: 初始的阻尼参数
    '''
    if b_init is None:
        # 用wiener反卷积的结果作为初始值
        b_init = basis_function_decomposition(p_meas, t_meas, basis_functions, n_basis)
    b = b_init.copy()
    mu = mu_init
    history = {
        'b': [b],
        'res': [],
        'mu': [mu]
    }    
    
    for iter in range(max_iter):
        # 正向模拟
        result = forward_simulation(qubit, control_pulse, b, t_meas)
        p_sim = np.array([res.expect[0] for res in result])  # 模拟的测量结果
        # 计算残差
        res = np.linalg.norm(p_meas - p_sim)
        history['res'].append(res)
        # 计算雅可比矩阵
        J = compute_jacobian(qubit, control_pulse, b, n_basis, t_meas)

        A = J.T @ J + mu * np.eye(n_basis) + reg * D(n_basis).T @ D(n_basis) # 正则化的Hessian矩阵

        delta_b = np.linalg.solve(A, J.T @ res)  # 计算参数更新
        b_trial = b + delta_b
        result_trial = forward_simulation(qubit, control_pulse, b_trial, t_meas)
        p_sim_trial = np.array([res.expect[0] for res in result_trial])
        res_trial = np.linalg.norm(p_meas - p_sim_trial)
        if res > res_trial:
            b = b_trial
            mu = max(mu / 2, 1e-8)  # 减小阻尼参数
            history['b'].append(b.copy())
            history['mu'].append(mu)
            if np.linalg.norm(delta_b) / (np.linalg.norm(b) + 1e-8) < tol:
                print(f'Converged at iteration {iter}')
                break
        else:
            mu = min(mu * 2, 1e8)  # 增加阻尼参数
            history['b'].append(b.copy())  # 保证b和mu的历史记录长度一致
            history['mu'].append(mu)

    return b, history
