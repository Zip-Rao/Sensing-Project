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
        import time 
        start = time.time()
        kernel = []
        t_list = control_pulse.t_list
        samples = range(0,len(t_list),1)
        t_samples = t_list[::1]
        psi_e = basis(qubit.n_levels, 1)

        if control_pulse.frame == 0:
            H_0 = qubit.get_hamiltonian(qubit.frequency)
        else:
            H_0 = qubit.get_hamiltonian_rwa(qubit.frequency)
        H_pulse = lambda t: control_pulse.get_hamiltonian_at(t)
        H_base = lambda t, args: H_0 + H_pulse(t)
        result_base = mesolve(H_base, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()])
        p_e_base = result_base.expect[0][-1]
        '''
        from joblib import Parallel, delayed
        def compute_kernel_at(i, t_list, qubit, control_pulse, p_e_base):
            t_i = t_list[i]
            stim_pulse = Signal(
                type = 3,
                t_list = t_list,
                amplitude = 0.0215,  # 待完善：自动调整幅度，远大于脉冲幅度
                center = t_i,
                width = 3  # 待完善：自动调整宽度，足够窄
            )
            stim_area = np.trapezoid(stim_pulse.signal, stim_pulse.t_list)

            qubit_t_curr = qubit.qubit_under_mag(stim_pulse)
            def H_total(t, args):
                index = np.searchsorted(stim_pulse.t_list, t)
                index = min(index, len(stim_pulse.t_list)-1)  # 确保索引不越界
                qubit_current = qubit_t_curr[index]
                H_0_current = qubit_current.get_hamiltonian(qubit.frequency) if control_pulse.frame == 0 else qubit_current.get_hamiltonian_rwa(qubit.frequency)
                return H_0_current + H_pulse(t)
            result_stim = mesolve(H_total, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()])
            p_e_stim = result_stim.expect[0][-1]
            return (p_e_stim-p_e_base) / stim_area  # TODO：调整系数，考虑刺激信号的幅度和宽度   
        kernel = Parallel(n_jobs=-1, backend = 'threading')(delayed(compute_kernel_at)(i, t_list, qubit, control_pulse, p_e_base) for i in samples)         
        '''
        
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
        
        end = time.time()
        print(f"Kernel computation took {end - start:.2f} seconds")
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
    
    def numerical_inverse(self, qubit, control_pulse, p_meas, t_list, B_guess, basis_type = 'fourier', n_basis=100, lambdas=100.0, max_iter=10, tol=1e-6):
        '''
        基于全密度矩阵模拟的数值反演算法，考虑非线性响应和退相干
        通过优化算法调整输入磁场信号，使得模拟的测量结果与实际测量结果p_meas尽可能接近
        :param p_meas: 实际测量的概率变化 (signal_len)
        :param t_list: 磁场信号的时间轴，测量时间轴应为磁场时间轴与脉冲时间轴的卷积
        :param basis_type: 用于表示输入磁场信号的基函数类型，支持'bspline','fourier','legendre'等
        :param n_basis: 基函数的数量
        :param lambdas: 正则化参数，控制解的平滑程度
        :param max_iter: 最大迭代次数
        :param tol: 收敛容忍度
        '''
        #basis_funcs = generate_basis_functions(basis_type, n_basis, t_meas[0], t_meas[-1])
        B_init = Signal(type = 6, t_list = t_list, n_basis = n_basis, basis_type = basis_type)
        # b_init = basis_function_decomposition(B_guess, t_list, B_init.basis_functions)
        #Phi = Signal(type = 1, t_list = np.linspace(0, 200, 400), amplitude = 0.01, rise = 10, fall = 10, center = 100)
        b_init = basis_function_decomposition(B_guess, t_list, B_init.basis_functions)
        B_init.update_signal(b=b_init)  # 将初始猜测的磁场信号分解到基函数上，得到初始的基函数系数，并更新B_init的信号
        print("figure of the initial guess signal")
        B_init.plot()
        b_opt, history = levenberg_marquardt(qubit, p_meas, t_list, control_pulse, b_init, B_init, lambdas, max_iter, tol)
        B_init.update_signal(b=b_opt)  # 更新B_init的参数b，得到优化后的磁场信号
        B_opt = B_init.signal
        #for k, phi_k in enumerate(B_init.basis_functions):
         #   B_opt += b_opt[k] * phi_k(t_list)
        return B_opt, history
    
    
    

        

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

def basis_function_decomposition(sig, t_array, basis_functions):
    '''
    将信号分解到基函数上，得到基函数系数
    :param signal: 待分解的信号 (signal_len)
    :param t_array: 信号的时间轴 (signal_len)
    :param basis_functions: 基函数列表
    '''
    from scipy.integrate import simpson
    from scipy.interpolate import interp1d
    if isinstance(sig, Signal):
        signal = sig.signal
    else:
        signal = sig
    n_basis = len(basis_functions)
    n_points = len(t_array)
    # 构建设计矩阵 A (n_points × n_basis)
    A = np.zeros((n_points, n_basis))
    for i, basis_func in enumerate(basis_functions):
        A[:, i] = basis_func(t_array)
    
    # 求解最小二乘问题 A @ b = signal
    b, _, _, _ = np.linalg.lstsq(A, signal, rcond=None)
    
    return b

def R(n, basis_type):
    '''
    正则化矩阵，鼓励解的平滑性
    '''
    if basis_type == 'fourier':
        freq_order = np.array([(i + 1)//2 for i in range(n)])  # 频率顺序，常数项频率为0
        R = np.diag(freq_order**2)  # 二阶导数正则化，频率越高惩罚越大
    else:
        D = np.zeros((n-2, n))
        for i in range(n-2):
            D[i, i] = 1
            D[i, i+1] = -2
            D[i, i+2] = 1
        R = D.T @ D  # 二阶差分正则化，鼓励解的平滑性
    return R

def forward_simulation(qubit, control_pulse, B_curr, t_meas, H_list, t_evolve_list):
    '''
    正向模拟，根据输入的基函数系数b生成磁场信号，并模拟延迟t时刻施加脉冲信号的测量结果
    t_meas表示滑动测量的时间轴，表示脉冲信号的延迟时间，应与等效时间采样的时间轴一致
    H_list表示各个延迟时间的哈密顿量列表，t_evolve表示演化的时间轴，起点是测量开始时间，终点是测量结束时间
    '''
    import time 
    start = time.time()
    B = B_curr
    result = []
    #qubit_t = qubit.qubit_under_mag(B)
    
    n_levels = qubit.n_levels
    n_op = qubit.n
    
    if False:
        def H_total(t, t_i, args):
            index = np.searchsorted(t_meas, t)
            index = max(index, 0)  # 确保索引非负
            index = min(index, len(B.t_list)-1)  # 确保索引不越界
            qubit_current = qubit_t[index]
            H_0_current = qubit_current.get_hamiltonian(qubit.frequency) if control_pulse.frame == 0 else qubit_current.get_hamiltonian_rwa(qubit.frequency)
            if t_i - 0.5 * control_pulse.t_list[-1] <= t <= t_i + 0.5 * control_pulse.t_list[-1]:
                return H_0_current + control_pulse.get_hamiltonian_at(t - t_i + 0.5 * control_pulse.t_list[-1])
            else:
                return H_0_current
    
    for i, t_i in enumerate(t_meas):

        H = H_list[i]
        t_evolve = t_evolve_list[i]
        dt_evolve = t_evolve[1] - t_evolve[0]
        options = {"store_states": True, "atol": 1e-10, "rtol": 1e-8, "max_step": dt_evolve / 2, "nsteps": 10000}  # 启用存储态

        H_total = QobjEvo(H, tlist=t_evolve, order = 1)

        #H = lambda t, args: H_total(t, t_i, args)
        #t_m = min(B.t_list[0], t_i - 0.5 * control_pulse.t_list[-1])  # 测量开始时间，考虑脉冲的持续时间
        #t_M = max(B.t_list[-1], t_i + 0.5 * control_pulse.t_list[-1])   # 测量结束时间
        #t_evolve = np.linspace(t_m, t_M, len(B.t_list) + len(control_pulse.t_list))  # 演化时间轴，分辨率足够高
        result.append(mesolve(H_total, qubit.state, t_evolve, [], e_ops=[basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()], options=options))
    end = time.time()
    print(f"Forward simulation time: {end - start:.2f} s")

    return result
# TODO:预计算H_total，优化adjoint
def compute_jacobian(qubit, control_pulse, B_curr, t_meas, results, H_list, t_evolve_list):
    '''
    伴随方法计算雅可比矩阵
    t_lists表示滑动测量的时间轴，表示脉冲信号的延迟时间，应与等效时间采样的时间轴一致
    t_evolve表示演化的时间轴，起点是
    '''
    import time
    start = time.time()
    N = len(t_meas)
    M = len(B_curr.params['b'])
    J = np.zeros((N, M))
    dim = qubit.state.shape[0]
    
    B = B_curr
    tB_list = B.t_list
    
    sensitivity = [qubit.frequency_sensitivity(flux) for flux in [qubit.flux + B.value_at(t) for t in tB_list]]
    # 扩展sensitivity的长度至整个测量时间轴
    #indices = np.searchsorted(tB_list, t_meas, side='right') - 1
    #indices = np.where(t_meas < tB_list[0], 0, indices)
    #indices = np.where(t_meas > tB_list[-1], len(sensitivity0) - 1, indices)
    #indices = np.clip(indices, 0, len(sensitivity0) - 1)
    #sensitivity = sensitivity0[indices]
    if False:
        def H_total(t, t_i, args):
            index = np.searchsorted(t_meas, t)
            index = max(index, 0)  # 确保索引非负
            index = min(index, len(tB_list)-1)  # 确保索引不越界
            qubit_current = qubit_t[index]
            H_0_current = qubit_current.get_hamiltonian(qubit.frequency) if control_pulse.frame == 0 else qubit_current.get_hamiltonian_rwa(qubit.frequency)
            if t_i - 0.5 * control_pulse.t_list[-1] <= t <= t_i + 0.5 * control_pulse.t_list[-1]:
                return H_0_current + control_pulse.get_hamiltonian_at(t - t_i + 0.5 * control_pulse.t_list[-1])
            else:    
                return H_0_current
        # 前向传播
    #result = forward_simulation(qubit, control_pulse, B_curr, t_lists)
    result = results
    states = [[] for _ in range(len(result))]
    
    for i, res in enumerate(result):
        states[i] = (np.array([(s * s.dag()).full() for s in res.states]))  # 存储每个时间点的密度矩阵
    p_sim = np.array([res.expect[0][-1] for res in result])  # 模拟的测量结果
    # 预计算
    #H_list = [lambda t, args, t_i=t_meas[i]: H_total(t, t_i, args).full() for i in range(N)]
    c_ops_list = [qubit.c_ops[k].full() for k in range(len(qubit.c_ops))]  # 假设c_ops不随时间变化
    c_ops_dag_list = [c_ops_list[k].conj().T for k in range(len(qubit.c_ops))]
    ops = [op.full() for op, _ in H_list[0]]
    mu_0 = basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()  # 期望值操作符
    mu_0 = mu_0.full().flatten()  # 转换为一维数组
    
    G_mat = qubit.n.full()  # G(t)，假设不随时间变化
    


    for i in range(N):
        start_i = time.time()
        t_i = t_meas[i]
        # 创建演化时间轴
        #t_m = min(tB_list[0], t_i - 0.5 * control_pulse.t_list[-1])  # 测量开始时间，考虑脉冲的持续时间
        #t_M = max(tB_list[-1], t_i + 0.5 * control_pulse.t_list[-1])   # 测量结束时间
        #t_evolve = np.linspace(t_m, t_M, len(tB_list) + len(control_pulse.t_list))  # 演化时间轴，分辨率足够高
        t_evolve = t_evolve_list[i]
        t_m = t_evolve[0]
        t_M = t_evolve[-1]
        N_e = len(t_evolve)
        # 创建基函数的时间映射，由于t_evolve的范围可能比tB_list更大，因此需要将s_list映射到tB_list的范围内
        t_mapped = np.empty_like(t_evolve)
        for s in range(N_e):
            if t_evolve[s] < tB_list[0]:
                t_mapped[s] = tB_list[0]
            elif t_evolve[s] >= tB_list[-1]:
                t_mapped[s] = tB_list[-1]
            else:
                index = np.searchsorted(tB_list, t_evolve[s])
                index = min(index, len(tB_list)-1)  # 确保索引不越界
                t_mapped[s] = tB_list[index]

        H_evolve = np.zeros((N_e, dim, dim), dtype=complex)
        for j in range(N_e):
            for m, (op, coeffs) in enumerate(H_list[i]):
                coeff = coeffs[j] if isinstance(coeffs, np.ndarray) else coeffs
                H_evolve[j] += ops[m] * coeff
        from scipy.integrate import solve_ivp

        def adjoint(s, mu, *args):
            mu = mu.reshape((dim, dim))
            t_curr = t_M - s
            j_float = (t_curr - t_m) / (t_M - t_m) * (N_e - 1)
            j_lo = max(int(np.floor(j_float)), 0)
            j_hi = min(j_lo + 1, N_e - 1)
            alpha = j_float - j_lo
            H = (1 - alpha) * H_evolve[j_lo] + alpha * H_evolve[j_hi]  # 线性插值计算当前时间点的哈密顿量
            
            rhs = 1j * (H @ mu - mu @ H)
            # 暂时不考虑耗散
            # for k in range(len(qubit.c_ops)):
            #     rhs += c_ops_dag_list[k] @ mu @ c_ops_list[k] - 0.5 * (c_ops_dag_list[k] @ c_ops_list[k] @ mu + mu @ c_ops_dag_list[k] @ c_ops_list[k])
            return rhs.flatten()
        
        #mu_0 = basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()  # 期望值操作符
        #mu_0 = mu_0.full().flatten()  # 转换为一维数组
        s_list = t_M - t_evolve[::-1]  # 伴随方程的时间点，逆序对应积分时间轴
        ds = s_list[1] - s_list[0]  # 步长限制，防止自适应步长跳过脉冲区域
        sol = solve_ivp(
                adjoint,
                t_span=(s_list[0], s_list[-1]),
                y0=mu_0,
                method='RK45',
                t_eval=s_list,
                rtol=1e-6,
                atol=1e-8,
                max_step=ds * 10,
            )
        
        N_s = len(sol.t)
        lambda_t = np.zeros((N_s, dim, dim), dtype=complex)
        for n in range(N_s):
            lambda_t[N_s - 1 - n] = sol.y[:, n].reshape((dim, dim))  # 伴随轨迹，逆序对应测量时间轴

        trace_values = np.zeros(N_s, dtype=complex)
        # for n in range(N_s):
        #     rho_n = states[i][n]      # ρ^i(t_n), 前向轨迹
        #     lam_n = lambda_t[n]    # λ(t_n), 伴随轨迹
        #     #G_mat =qubit.n.full()  # G(t_n)，
        #     trace_values[n] = np.trace(lam_n @ (G_mat @ rho_n - rho_n @ G_mat))
        rho_stack = states[i][:N_s]
        commutator = G_mat @ rho_stack - rho_stack @ G_mat  # 计算[G, ρ(t_n)]，得到一个形状为(N_s, dim, dim)的数组
        trace_values = np.trace(lambda_t @ commutator, axis1 = -2, axis2 = -1)  # 计算Tr[λ(t_n) @ (G @ ρ(t_n) - ρ(t_n) @ G)]，得到一个长度为N_s的数组
        g_values = (-1j * trace_values)  
        for s in range(N_s):
            if t_evolve[s] < tB_list[0]:
                g_values[s] *= 0.0  # 信号范围外的时间点对梯度没有贡献
            elif t_evolve[s] >= tB_list[-1]:
                g_values[s] *= 0.0  # 信号范围外的时间点对梯度没有贡献
            else:
                index = np.searchsorted(tB_list, t_evolve[s])
                index = min(index, len(sensitivity)-1)  # 确保索引不越界
                g_values[s] *= sensitivity[index]
        for k, phi_k in enumerate(B_curr.basis_functions):
            integrand_k = g_values * phi_k(t_mapped)
            integral = np.trapezoid(integrand_k, t_evolve)
            
            J[i, k] = np.real(integral)  # 取实部，虚部应为数值零
        end_i = time.time()
        #print(f"Jacobian row {i} computation time: {end_i - start_i:.2f} s")
    
    end = time.time()
    print(f"Jacobian computation time: {end - start:.2f} s")

    return J

def compute_jacobian_finite_difference(qubit, control_pulse, B_curr, t_meas, p_sim, H_list, t_evolve_list, epsilon=1e-6):
    '''
    有限差分计算雅可比矩阵，作为伴随方法的对照
    '''
    N = len(t_meas)
    M = len(B_curr.params['b'])
    J_fd = np.zeros((N, M))
    t_signal = B_curr.t_list

    def build_h_for_current_qubit():
        H_local = [[] for _ in range(N)]
        t_evolve_local = [[] for _ in range(N)]
        for i, t_delay in enumerate(t_meas):
            delta = t_delay - 0.5 * control_pulse.t_list[-1]
            t_start = min(t_signal[0], delta)
            t_end = max(t_signal[-1], delta + control_pulse.t_list[-1])
            N_e = len(t_signal) + len(control_pulse.t_list) - 1
            t_evolve = np.linspace(t_start, t_end, N_e)
            t_evolve_local[i] = t_evolve

            freq_coeffs = np.zeros(N_e)
            for j, t in enumerate(t_evolve):
                if t_signal[0] <= t <= t_signal[-1]:
                    index = np.searchsorted(t_signal, t)
                    index = min(index, len(t_signal) - 1)
                    freq_coeffs[j] = qubit.freq_coeffs[index] if control_pulse.frame == 0 else qubit.freq_coeffs[index] - control_pulse.omega_d
                else:
                    freq_coeffs[j] = qubit.frequency if control_pulse.frame == 0 else qubit.frequency - control_pulse.omega_d  # 信号范围外使用未受扰动的qubit频率

            id_coeffs = np.ones(N_e, dtype=complex)
            if control_pulse.frame == 0:
                H_local[i].append([qubit.anharmonicity * qubit.n * (qubit.n - 1) * 0.5, id_coeffs])
                H_local[i].append([qubit.n + 0.5 * qeye(qubit.n_levels), freq_coeffs])
            else:
                H_local[i].append([qubit.anharmonicity * qubit.n * (qubit.n - 1) * 0.5, id_coeffs])
                H_local[i].append([qubit.n, freq_coeffs])

            for op, coeffs in control_pulse.hamiltonian:
                coeff_global = np.zeros(N_e, dtype=complex)
                for j, t in enumerate(t_evolve):
                    t_loc = t - delta
                    if 0 <= t_loc <= control_pulse.t_list[-1]:
                        index = np.searchsorted(control_pulse.t_list, t_loc)
                        index = min(index, len(control_pulse.t_list) - 1)
                        coeff_global[j] = coeffs[index]
                H_local[i].append([op, coeff_global])
        return H_local, t_evolve_local

    for k in range(M):
        b_perturb = B_curr.params['b'].copy()
        b_perturb[k] += epsilon
        B_perturb = B_curr.copy()
        B_perturb.update_signal(b=b_perturb)

        # 有限差分每一列都必须重新更新qubit频率并构建对应H，否则结果会退化为0
        qubit.qubit_in_mag(B_perturb)
        H_perturb, t_evolve_perturb = build_h_for_current_qubit()
        result_perturb = forward_simulation(qubit, control_pulse, B_perturb, t_meas, H_perturb, t_evolve_perturb)
        p_sim_perturb = np.array([res.expect[0][-1] for res in result_perturb])
        J_fd[:, k] = (p_sim_perturb - p_sim) / epsilon

    # 恢复基准点，避免影响后续迭代状态
    qubit.qubit_in_mag(B_curr)
    return J_fd


    
def levenberg_marquardt(qubit, p_meas, t_list, control_pulse, b_init, B_init, reg, max_iter = 50, tol = 1e-6, mu_init = 1e-3):
    '''
    Levenberg-Marquardt优化算法，用于最小化模拟测量结果与实际测量结果之间的差异
    :param qubit: 量子比特对象
    :param p_meas: 实际测量的概率变化 (signal_len)
    :param t_list: 磁场信号的时间轴，测量时间轴应为磁场时间轴与脉冲时间轴的卷积
    :param control_pulse: 控制脉冲对象，用于模拟测量结果
    :param b_init: 初始的基函数系数
    :param B_init: 初始的磁场信号，类型为Signal
    :param reg: 正则化参数
    :param max_iter: 最大迭代次数
    :param tol: 收敛容忍度
    :param mu_init: 初始的阻尼参数
    '''
    n_basis = len(b_init) 
    b = b_init.copy()
    B_curr = B_init.copy()
    basis_type = B_curr.params['basis_type']
    mu = mu_init
    qubit.qubit_in_mag(B_curr)
    history = {
        'b': [b],
        'res': [],
        'mu': [mu]
    }   
    #预计算
     
    # 标定测量时间轴
    meas_start = t_list[0] - 0.5 * control_pulse.t_list[-1]  # 测量开始时间，考虑脉冲的持续时间
    meas_end = t_list[-1] + 0.5 * control_pulse.t_list[-1]   # 测量结束时间
    t_meas = np.linspace(meas_start, meas_end, len(t_list) + len(control_pulse.t_list) - 1)  # 测量时间轴

    # 给定延迟时间的哈密顿量列表
    def H(qubit):
        H_list = [[] for _ in range(len(t_meas))]  # 每个延迟时间对应一个哈密顿量列表
        t_evolve_list = [[] for _ in range(len(t_meas))]  # 每个延迟时间对应一个演化时间轴
        for i, t_delay in enumerate(t_meas):
            delta = t_delay - 0.5 * control_pulse.t_list[-1]  # 调整时间，使得t_delay表示脉冲中心的时间
            t_start = min(t_list[0], delta)
            t_end = max(t_list[-1], delta + control_pulse.t_list[-1])
            N_e = len(t_list) + len(control_pulse.t_list) - 1
            t_evolve = np.linspace(t_start, t_end, N_e)  # 演化时间轴，分辨率足够高
            t_evolve_list[i] = t_evolve
            freq_coeffs = np.zeros(N_e)
            for j, t in enumerate(t_evolve):
                if t_list[0] <= t <= t_list[-1]:
                    index = np.searchsorted(t_list, t)
                    index = min(index, len(t_list)-1)  # 确保索引不越界
                    freq_coeffs[j] = qubit.freq_coeffs[index] if control_pulse.frame == 0 else qubit.freq_coeffs[index] - control_pulse.omega_d
                else:
                    freq_coeffs[j] = qubit.frequency if control_pulse.frame == 0 else qubit.frequency - control_pulse.omega_d  # 信号范围外使用未受扰动的qubit频率，与protocal.py一致
            id = np.ones(N_e, dtype = complex)
            if control_pulse.frame == 0:
                H_list[i].append([qubit.anharmonicity * qubit.n * (qubit.n - 1) * 0.5, id])
                H_list[i].append([qubit.n + 0.5 * qeye(qubit.n_levels), freq_coeffs])
            else:
                H_list[i].append([qubit.anharmonicity * qubit.n * (qubit.n - 1) * 0.5, id])
                H_list[i].append([qubit.n, freq_coeffs])
            for op, coeffs in control_pulse.hamiltonian:
                coeff_global = np.zeros(N_e, dtype = complex)
                for j, t in enumerate(t_evolve):
                    t_loc = t - delta
                    if 0 <= t_loc <= control_pulse.t_list[-1]:
                        index = np.searchsorted(control_pulse.t_list, t_loc)
                        index = min(index, len(control_pulse.t_list)-1)  # 确保索引不越界
                        coeff_global[j] = coeffs[index]
                H_list[i].append([op, coeff_global])

        return H_list, t_evolve_list
                    
    H_curr_list, t_evolve_list = H(qubit)  # 给定延迟时间的哈密顿量列表和演化时间轴
    for iter in range(max_iter):
        #B_curr.plot()
        # 正向模拟

        results = forward_simulation(qubit, control_pulse, B_curr, t_meas, H_curr_list, t_evolve_list)
        p_sim = np.array([result.expect[0][-1] for result in results])  # 模拟的测量结果
        # 计算残差
        res = p_meas - p_sim
        history['res'].append(res)
        # 计算雅可比矩阵
        J = compute_jacobian(qubit, control_pulse, B_curr, t_meas, results, H_curr_list, t_evolve_list)
        #J = compute_jacobian_finite_difference(qubit, control_pulse, B_curr, t_meas, p_sim, H_curr_list, t_evolve_list)
                # ===== 分离调试：检查 sensitivity * phi_k 是否匹配频率变化的有限差分 =====
        eps = 1e-5
        k_test = 0
        phi_k = B_curr.basis_functions[k_test]

        # 当前频率
        freq_orig = qubit.freq_coeffs.copy()

        # 扰动后频率
        b_plus = b.copy(); b_plus[k_test] += eps
        B_plus = B_curr.copy(); B_plus.update_signal(b=b_plus)
        qubit.qubit_in_mag(B_plus)
        freq_plus = qubit.freq_coeffs.copy()
        qubit.qubit_in_mag(B_curr)  # 恢复

        # 有限差分：频率对 b_k 的导数
        dfreq_fd = (freq_plus - freq_orig) / eps  # 长度 = len(t_list)

        # adjoint chain rule：sensitivity * phi_k(t_list)
        tB = B_curr.t_list
        sensitivity_arr = np.array([
            qubit.frequency_sensitivity(qubit.flux + B_curr.value_at(t)) for t in tB
        ])
        dfreq_adj = sensitivity_arr * phi_k(np.array(tB))

        print("=== Chain Rule 验证 ===")
        for idx in [0, len(tB)//4, len(tB)//2, 3*len(tB)//4, len(tB)-1]:
            print(f"  t={tB[idx]:.2f}: fd={dfreq_fd[idx]:.6f}, adj={dfreq_adj[idx]:.6f}, ratio={dfreq_adj[idx]/(dfreq_fd[idx]+1e-30):.4f}")
                
        # ===== 有限差分验证（调试完删掉）=====
        eps = 1e-5
        k_test = 0  # 测试第 0 列
        b_plus = b.copy(); b_plus[k_test] += eps
        B_plus = B_curr.copy(); 
        B_plus.update_signal(b=b_plus)
        qubit.qubit_in_mag(B_plus)
        H_plus, t_ev_plus = H(qubit)
        res_plus = forward_simulation(qubit, control_pulse, B_plus, t_meas, H_plus, t_ev_plus)
        p_plus = np.array([r.expect[0][-1] for r in res_plus])
        J_fd = (p_plus - p_sim) / eps
        qubit.qubit_in_mag(B_curr)  # 恢复
        print("dimension check: ", J.shape, J_fd.shape)
        print(f"Adjoint J[:,{k_test}] = {J[:5, k_test]}")
        print(f"FinDiff J[:,{k_test}] = {J_fd[:5]}")
        print(f"Ratio: {J[:5, k_test] / (J_fd[:5] + 1e-30)}")
        # ===== 验证结束 =====
        A = J.T @ J + mu * np.eye(n_basis) + reg * R(n_basis, basis_type)  # 正则化的Hessian矩阵

        delta_b = np.linalg.solve(A, J.T @ res)  # 计算参数更新
        b_trial = b + delta_b
        B_trial = B_curr.copy()  # 更新磁场信号
        B_trial.update_signal(b=b_trial)    
        qubit.qubit_in_mag(B_trial)  # 更新量子比特状态
        H_curr_list, _ = H(qubit)  # 计算更新后的哈密顿量列表
        result_trial = forward_simulation(qubit, control_pulse, B_trial, t_meas, H_curr_list, t_evolve_list)  # 计算更新后的模拟结果
        p_sim_trial = np.array([result.expect[0][-1] for result in result_trial])
        res_trial = p_meas - p_sim_trial
        if np.linalg.norm(res) > np.linalg.norm(res_trial):
            b = b_trial
            B_curr = B_trial
            mu = max(mu / 2, 1e-8)  # 减小阻尼参数
            history['b'].append(b.copy())
            history['mu'].append(mu)
            if np.linalg.norm(delta_b) / (np.linalg.norm(b) + 1e-8) < tol:
                print(f'Converged at iteration {iter}')
                break
        else:
            mu = min(mu * 2, 1e8)  # 增加阻尼参数
            qubit.qubit_in_mag(B_curr)  # 恢复原始状态
            H_curr_list, _ = H(qubit)  # 恢复原始哈密顿量列表
            history['b'].append(b.copy())  # 保证b和mu的历史记录长度一致
            history['mu'].append(mu)
        
        # 打印迭代信息
        print(f'Iter {iter}: res={np.linalg.norm(res):.6f}, res_trial={np.linalg.norm(res_trial):.6f}, mu={mu:.6e}')

    return b, history
