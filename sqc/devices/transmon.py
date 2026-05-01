import numpy as np
import math
from qutip import *
from dataclasses import dataclass

from sqc.devices.base import Device
from sqc.control.flux_signal import Signal
'''
qubit类
创建Transmon qubit对象，包括实验相关参数：电容，电感，频率，非谐性；弛豫时间，噪声模型等
实现态投影计算，单比特门等功能
本仿真平台全部基于Fock基底展开，使用自然单位制（ħ=1）
同时定义了一些两比特门函数，如CiSWAP门，CZ门等
'''


@dataclass(frozen=True)
class QubitSpec:
    """Pure parameter description of a flux-tunable Transmon qubit.

    Units follow the legacy code: EC and EJ are in rad*GHz, time is in ns,
    and flux is in units of Phi_0.
    """

    name: str
    EC: float
    EJ: float
    T1: float
    T2: float
    flux_bias: float = 0.0
    n_levels: int = 3

    def EJ_at(self, flux: float | None = None) -> float:
        """Return Josephson energy at the given flux bias."""
        f = self.flux_bias if flux is None else flux
        return self.EJ * abs(math.cos(math.pi * f))

    def frequency(self, flux: float | None = None) -> float:
        """Return f01(Phi) = sqrt(8 EJ(Phi) EC) - EC."""
        return float(np.sqrt(8 * self.EJ_at(flux) * self.EC) - self.EC)

    def anharmonicity(self) -> float:
        """Return the transmon anharmonicity with legacy sign convention."""
        return -self.EC

    def sensitivity(self, flux: float | None = None, delta: float = 1e-6) -> float:
        """Return df01/dPhi by centered finite difference."""
        f = self.flux_bias if flux is None else flux
        return (self.frequency(f + delta) - self.frequency(f - delta)) / (2 * delta)

    def with_flux(self, new_flux: float) -> "QubitSpec":
        """Return a new spec with a changed static flux bias."""
        return QubitSpec(
            name=self.name,
            EC=self.EC,
            EJ=self.EJ,
            T1=self.T1,
            T2=self.T2,
            flux_bias=new_flux,
            n_levels=self.n_levels,
        )

    @staticmethod
    def optimal_work_point() -> float:
        """Flux with maximal first-order sensitivity, in Phi_0 units."""
        return float(np.arctan(np.sqrt(2)) / np.pi)


class TransmonQubit(Device):
    def __init__(self, EC, EJ, T1, T2, flux = 0, state=0,n_levels=3, name="Q"):
        '''
        初始化Transmon Qubit对象
        :param EC: 电容能量（rad·GHz）
        :param EJ: 约瑟夫森能量（rad·GHz）
        :param T1: 弛豫时间（ns）
        :param T2: 退相干时间（ns）
        :param flux: 外加磁通（单位：Φ0）
        :param state: 初始能级状态（整数或Qobj）
        :param n_levels: 能级数
        '''
        self._name = name
        self._spec = QubitSpec(
            name=name,
            EC=EC,
            EJ=EJ,
            T1=T1,
            T2=T2,
            flux_bias=flux,
            n_levels=n_levels,
        )
        self.EC = EC  # 电容能量
        self.EJ = EJ * abs(math.cos(math.pi * flux)) # 约瑟夫森能量
        self.EJ_0 = EJ  # 零磁通时的约瑟夫森能量
        self.flux = flux  # 外加磁通
        self.n_levels = n_levels  # 能级数

        # 初始化态
        if(isinstance(state, int) and 0 <= state < n_levels):
            self.state = basis(n_levels, state)
        elif(isinstance(state, Qobj) and state.dims == [[n_levels], [1]]):
            self.state = state.unit()

        # 计算频率和非谐性
        self.frequency = self.calculate_frequency()
        self.anharmonicity = self.calculate_anharmonicity()

        # 弛豫时间和退相干时间
        self.T1 = T1
        self.T2 = T2

        self.a = destroy(n_levels)
        self.a_dag = self.a.dag()
        self.n = num(n_levels)
        self.I = qeye(n_levels)

        # 计算哈密顿量
        self.hamiltonian = self.get_hamiltonian()
        self.c_ops = self.get_collapse_operators()
        
    @property
    def name(self):
        """Device name."""
        return self._name

    def spec(self):
        """Return a pure parameter snapshot for new sqc APIs."""
        return self._spec

    def hilbert_dim(self):
        """Return the Hilbert-space truncation dimension."""
        return self.n_levels

    def hamiltonian_static(self):
        """Return the static Hamiltonian."""
        return self.hamiltonian

    def collapse_operators(self):
        """Return Lindblad collapse operators."""
        return self.c_ops

        
    def calculate_frequency(self):
        '''
        计算Transmon Qubit的频率（GHz）
        '''
        return np.sqrt(8 * self.EJ * self.EC) - self.EC
    def calculate_anharmonicity(self):
        '''
        计算Transmon Qubit的非谐性（GHz）
        '''
        return -self.EC

    def frequency_sensitivity(self, flux, delta_flux=1e-6):
        '''
        计算Transmon Qubit频率对磁通的敏感度（GHz/Φ0）
        使用中心差分法计算数值导数
        :param delta_flux: 磁通变化量（Φ0），用于有限差分计算
        '''
        # 计算磁通增加delta_flux时的频率
        flux_plus = flux + delta_flux
        EJ_plus = self.EJ_0 * abs(math.cos(math.pi * flux_plus))
        f_plus = np.sqrt(8 * EJ_plus * self.EC) - self.EC

        # 计算磁通减少delta_flux时的频率
        flux_minus = flux - delta_flux
        EJ_minus = self.EJ_0 * abs(math.cos(math.pi * flux_minus))
        f_minus = np.sqrt(8 * EJ_minus * self.EC) - self.EC

        # 中心差分计算导数
        df_dphi = (f_plus - f_minus) / (2 * delta_flux)
        return df_dphi

    def get_hamiltonian(self):
        '''
        计算Transmon Qubit的哈密顿量 ，
        '''
        n_levels = self.n_levels
        a = self.a
        a_dag = self.a_dag
        n = self.n
        H_0 = (-self.EJ + 0.25 * self.EC) * qeye(n_levels)
        H_1 = self.frequency * (n + 0.5 * qeye(n_levels))
        H_2 = (self.anharmonicity / 2) * (n * n - n)
        H = H_0 + H_1 + H_2
        return H

    def get_hamiltonian_rwa(self, omega_d):
        '''
        计算Transmon Qubit在旋转波近似下的哈密顿量
        :param omega_d: 驱动频率（GHz）
        '''
        n_levels = self.n_levels
        a = self.a
        a_dag = self.a_dag
        n = self.n
        Delta = self.frequency - omega_d
        H_0 = Delta * n
        H_1 = (self.anharmonicity / 2) * (n * n - n)
        H_rwa = H_0 + H_1
        return H_rwa
    
    def get_collapse_operators(self):
        '''
        返回Lindblad算符列表，包含退相干噪声算符
        '''
        gamma_1 = 1 / self.T1  # 弛豫率
        gamma_phi = 1 / self.T2 - 0.5 * gamma_1  # 纯退相干率
        c_ops = []
        # 添加退相干噪声算符
        c_ops.append(np.sqrt(gamma_1) * self.a)
        c_ops.append(np.sqrt(gamma_phi) * self.n)
        return c_ops
    
    def generate_1f_noise(self, t_lists, amplitude, f_min, f_max):
        '''
        生成1/f噪声时间序列
        :param amplitude: 噪声幅度
        :param f_min: 最小频率（Hz）
        :param f_max: 最大频率（Hz）
        :param t_lists: 时间列表（ns）  
        '''
        dt = t_lists[1] - t_lists[0]  # 时间步长
        n = len(t_lists)
        freqs = np.fft.fftfreq(n, dt)
        spectrum = np.zeros(n, dtype=complex)
        for i in range(1, n // 2):
            f = abs(freqs[i])
            if f_min <= f <= f_max:
                spectrum[i] = amplitude / np.sqrt(f) * (np.random.normal() + 1j * np.random.normal())
        # 对称化，保证实数
        spectrum[n // 2 + 1:] = np.conj(spectrum[1:n // 2][::-1])
        noise = np.fft.ifft(spectrum).real
        return noise



    def calculate_state_projection(self, target_state):
        '''
        计算Transmon Qubit处在目标态的投影率
        :param target_state: 目标态（整数或Qobj）
        '''
        if(isinstance(target_state, int) and 0 <= target_state < self.n_levels):
            target = basis(self.n_levels, target_state)
        elif(isinstance(target_state, Qobj) and target_state.dims == [[self.n_levels], [1]]):
            target = target_state.unit()
        proj = target.dag() * target
        prob = expect(proj, self.state)
        return float(prob.real)
    
    def qubit_under_mag(self,Phi_signal:Signal, is_noise = False):
        '''
        计算Transmon Qubit在外加磁通信号下的频率变化
        :param Phi_signal: 外加磁通信号（Signal对象）
        '''
        qubit = []
        if is_noise:
            noise = self.generate_1f_noise(Phi_signal.t_list, amplitude=0.001, f_min=1e-3, f_max=1e3)
        for t in Phi_signal.t_list:
            qubit.append(TransmonQubit(
                EC=self.EC,
                EJ=self.EJ_0,
                T1=self.T1,
                T2=self.T2,
                flux=(self.flux + Phi_signal.value_at(t) + noise[t]) if is_noise else (self.flux + Phi_signal.value_at(t)),
                state=self.state,
                n_levels=self.n_levels
            ))
        
        return qubit
    
    def qubit_under_mag_hamiltonian(self, qubit_t, t_list, frame = 0, omega_d = None):
        '''
        当qubit处于时变磁场下时，计算qubit的哈密顿量随时间的变化
        '''
        H_list = []
        freq_coeffs = np.zeros(len(t_list))

        for i, t in enumerate(t_list):
            qubit_current = qubit_t[i]
            if frame == 0:
                freq_coeffs[i] = qubit_current.frequency
            else:
                freq_coeffs[i] = qubit_current.frequency - omega_d

        H_list.append(qubit_current.anharmonicity * 0.5 * (self.n * self.n - self.n))
        H_list.append([self.n + 0.5 * qeye(self.n_levels), freq_coeffs])
        return H_list
    
    def qubit_in_mag(self, Phi_signal, frame = 0, omega_d = None):
        '''
        优化后的计算Transmon Qubit在外加磁通信号下的频率变化的方法，直接计算哈密顿量随时间的变化，而不是每个时间点都构建一个新的TransmonQubit对象
        '''
        from sqc.simulation.hamiltonian import HamiltonianBuilder

        self.isinmag = True
        self.mag_signal = Phi_signal
        self.H_list, _ = HamiltonianBuilder.build(
            self,
            flux_signal=Phi_signal,
            frame="lab" if frame == 0 else "rotating",
            omega_d=omega_d,
        )
        self.freq_coeffs = np.asarray(self.H_list[1][1])



    def change_flux(self, flux):
        '''
        改变Transmon Qubit的外加磁通，更新频率和哈密顿量
        :param flux: 新的外加磁通（单位：Φ0）
        '''
        self.flux = flux
        self.EJ = self.EJ_0 * abs(math.cos(math.pi * flux))
        self._spec = self._spec.with_flux(flux)
        self.frequency = self.calculate_frequency()
        self.anharmonicity = self.calculate_anharmonicity()
        self.hamiltonian = self.get_hamiltonian()

    def sensitivity(self):
        '''
        计算Transmon Qubit的灵敏度
        '''

    
    def optimal_work_point(self):
        '''
        计算Transmon Qubit的最优工作位置
        此处的工作位置意为频率对磁通变化斜率最大的磁通点
        '''
        Phi = np.arctan(np.sqrt(2))
        return Phi

    def flux_noise(self):
        '''
        计算Transmon Qubit的磁通噪声影响
        '''

        pass
    def ideal_gate(self, theta, phi):
        '''
        计算理想单比特门操作的演化算符
        :param theta: 旋转角度（弧度）
        :param phi: 旋转轴的相位（弧度）
        '''
        a = destroy(2)
        U = (-1j * theta / 2 * ( (a + a.dag()) * np.cos(phi) - 1j * (-a + a.dag()) * np.sin(phi))).expm()
        return U
    def simulate_gate(self, theta, phi, T, sigma):
        '''
        模拟单比特门操作，使用DRAG脉冲，计算门的有效矩阵，平均泄露和保真度
        :param theta: 旋转角度（弧度）
        :param phi: 旋转轴的相位（弧度）
        :param T: 门操作时间（ns）
        :param sigma: 门操作的宽度（ns）
        :param qubit: 进行门操作的TransmonQubit对象
        '''
        def Omega_I(t, args):
                A = args["A"]
                t0 = args["t0"]
                sigma = args["sigma"]
                beta = args["beta"]
                phi = args["phi"]

                env = A * np.exp(-0.5 * ((t - t0) / sigma) ** 2)
                drag = - A * (t - t0) / (sigma ** 2) * np.exp(-0.5 * ((t - t0) / sigma) ** 2)
                return env * np.cos(phi) - beta * drag * np.sin(phi)
        def Omega_Q(t, args):
                A = args["A"]
                t0 = args["t0"]
                sigma = args["sigma"]
                beta = args["beta"]
                phi = args["phi"]

                env = A * np.exp(-0.5 * ((t - t0) / sigma) ** 2)
                drag = - A * (t - t0) / (sigma ** 2) * np.exp(-0.5 * ((t - t0) / sigma) ** 2)
                return env * np.sin(phi) + beta * drag * np.cos(phi)
        

        t_0 = T / 2
        t_list = np.linspace(0, T, 1000)
        Omega = (theta / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((t_list - t_0) / sigma) ** 2) # 归一化的高斯脉冲

        beta = -1.0/(self.anharmonicity)
        args = {"A": theta / (sigma * np.sqrt(2 * np.pi)), "t0": t_0, "sigma": sigma, "beta": beta, "phi": phi}
        ket_0 = basis(self.n_levels, 0)
        ket_1 = basis(self.n_levels, 1)
        a = destroy(self.n_levels)
        H = [self.get_hamiltonian_rwa(self.frequency), [0.5 * (a + a.dag()) , Omega_I] ,[0.5 * (-1j) * (a - a.dag()) , Omega_Q]]
        result_0= sesolve(H, ket_0, t_list, args=args)
        result_1= sesolve(H, ket_1, t_list, args=args)
        final_state_0 = result_0.states[-1]
        final_state_1 = result_1.states[-1]

        amp00 = ket_0.overlap(final_state_0)  # <0|psi0>
        amp01 = ket_0.overlap(final_state_1)  # <0|psi1>
        amp10 = ket_1.overlap(final_state_0)  # <1|psi0>
        amp11 = ket_1.overlap(final_state_1)  # <1|psi1>

        U_eff = Qobj(np.array([[amp00, amp01],
                       [amp10, amp11]], dtype=complex),
             dims=[[2], [2]])
        U_ideal = self.ideal_gate(theta, phi)
        leakage = 1 - 0.5 * (np.sum(np.abs(U_eff.full())**2))
        F = (abs((U_ideal.dag() @ U_eff).tr())**2 + 2)/ (2 * 3)
        return U_eff, leakage, F
    
class Cavity:
    def __init__(self, frequencies, n_levels):
        '''
        初始化Cavity对象
        :param frequencies: 不同模式的腔体频率列表（GHz）
        :param M: 模式数
        :param n_levels: 不同模式的截断能级数
        '''
        self.frequencies = frequencies
        self.M = len(frequencies)
        self.n_levels = n_levels
        self.build_operators()
        self.hamiltonian = self.get_hamiltonian()

    def build_operators(self):
        '''
        构建Cavity不同模式的算符列表
        '''
        self.a = []
        self.adag = []
        self.n = []
        self.I = tensor([qeye(n) for n in self.n_levels])

        for M in range(self.M):
            ops = []
            for m in range(self.M):
                if m == M:
                    a_m = destroy(self.n_levels[m])
                    ops.append(a_m)
                else:
                    ops.append(qeye(self.n_levels[m]))
            a_M = tensor(ops)
            self.a.append(a_M)
            self.adag.append(a_M.dag())
            self.n.append(a_M.dag() * a_M)
    
    def get_hamiltonian(self):
        '''
        计算Cavity的哈密顿量
        暂时只考虑线性项
        '''
        H = 0 * self.I
        for i in range(self.M):
            H += self.frequencies[i] * self.n[i]
        return H
    def get_hamiltonian_rot(self, omega_d):
        '''
        计算Cavity在旋转框架下的哈密顿量
        :param omega_d: 驱动频率（GHz）
        '''
        H = 0 * self.I
        for i in range(self.M):
            H += (self.frequencies[i] - omega_d) * self.n[i]
        return H
# 两个Transmon qubit通过一个coupler与多模腔耦合系统
class Coupled_System:
    def __init__(self, qubit1, qubit2, cavity):
        '''
        初始化Coupled_System对象，为了简化模型，不显式构建coupler qubit，而是等效为一个可调的交换耦合
        :param qubit1: 第一个TransmonQubit对象
        :param qubit2: 第二个TransmonQubit对象
        :param qubitc: Coupler对象
        :param cavity: Cavity对象
        :param g1c: qubit1与coupler的耦合强度（GHz）
        :param g2c: qubit2与coupler的耦合强度（GHz）
        :param gccav: coupler与cavity的耦合强度（GHz）
        '''
        self.qubit1 = qubit1
        self.qubit2 = qubit2
        self.cavity = cavity
        self.a1 = tensor(self.qubit1.a, self.qubit2.I, self.cavity.I)
        self.a2 = tensor(self.qubit1.I, self.qubit2.a, self.cavity.I)
        self.acav = [tensor(self.qubit1.I, self.qubit2.I, self.cavity.a[m]) for m in range(self.cavity.M)]
        self.g = self.initialize_g()
        

    def initialize_g(self):
        '''
        初始化qubit-qubit交换耦合强度，等效于coupler的调节，返回2*M的列表，分别对应coupler与每个腔模式的耦合强度
        '''
        g = [[],[]]
        for m in range(self.cavity.M):
            g[0].append(2 * np.pi * 0.01)  # qubit1与coupler的耦合强度（GHz）
            g[1].append(2 * np.pi * 0.01)  # qubit2与coupler的耦合强度（GHz）
        return g
        
    def control_g(self, g, is_time_dependent = False, f = None):
        '''
        控制qubit-qubit交换耦合强度，等效于coupler的调节
        :param g: 2*M的列表，分别对应coupler与每个腔模式的耦合强度
        :param is_time_dependent: 是否为时间依赖的耦合
        :param f: 如果是时间依赖的耦合，f应该是一个函数，输入时间t，输出2*M的列表，分别对应coupler与每个腔模式的耦合强度
        '''
        self.g = g
        if not is_time_dependent:
            H_qc1 = sum([self.g[0][m] * (self.a1.dag() * self.acav[m] + self.a1 * self.acav[m].dag()) for m in range(self.cavity.M)])
            H_qc2 = sum([self.g[1][m] * (self.a2.dag() * self.acav[m] + self.a2 * self.acav[m].dag()) for m in range(self.cavity.M)])
            self.H = self.H_0 + H_qc1 + H_qc2
        else:
            self.H = [self.H_0]
            for M in range(self.cavity.M):
                H_qc1 = self.g[0][M] * (self.a1.dag() * self.acav[M] + self.a1 * self.acav[M].dag())
                H_qc2 = self.g[1][M] * (self.a2.dag() * self.acav[M] + self.a2 * self.acav[M].dag())
                f1 = lambda t, args: f(t, args)[0][M]
                f2 = lambda t, args: f(t, args)[1][M]
                self.H.append([H_qc1, f1])
                self.H.append([H_qc2, f2])
                
                
    def build_system(self, is_rwa=False, omega_d = None):
        '''
        构建两个Transmon qubit通过一个coupler与多模腔耦合的系统
        :param qubit1: 第一个TransmonQubit对象
        :param qubit2: 第二个TransmonQubit对象
        :param qubitc: Coupler对象
        :param cavity: Cavity对象
        :param g1c: qubit1与coupler的耦合强度（GHz）
        :param g2c: qubit2与coupler的耦合强度（GHz）
        :param gccav: coupler与cavity的耦合强度（GHz）
        '''
        # 孤立系统的哈密顿量
        if is_rwa:
            self.H_q1 = tensor(self.qubit1.get_hamiltonian_rwa(omega_d), self.qubit2.I, self.cavity.I)
            self.H_q2 = tensor(self.qubit1.I, self.qubit2.get_hamiltonian_rwa(omega_d), self.cavity.I)
            self.H_cav = tensor(self.qubit1.I, self.qubit2.I, self.cavity.get_hamiltonian_rot(omega_d))
        else:
            self.H_q1 = tensor(self.qubit1.get_hamiltonian(), self.qubit2.I, self.cavity.I)
            self.H_q2 = tensor(self.qubit1.I, self.qubit2.get_hamiltonian(), self.cavity.I)
            self.H_cav = tensor(self.qubit1.I, self.qubit2.I, self.cavity.get_hamiltonian())

        self.H_0 = self.H_q1 + self.H_q2 + self.H_cav
        # 等效qubit-cavity耦合
        H_qc1 = sum([self.g[0][m] * (self.a1.dag() * self.acav[m] + self.a1 * self.acav[m].dag()) for m in range(self.cavity.M)])
        H_qc2 = sum([self.g[1][m] * (self.a2.dag() * self.acav[m] + self.a2 * self.acav[m].dag()) for m in range(self.cavity.M)])

        self.H = self.H_q1 + self.H_q2 + self.H_cav + H_qc1 + H_qc2


    def prepare_ket11(self, T, sigma):
        '''
        制备高保真度的11态，关闭与腔的耦合，对两个qubit施加DRAG pi脉冲
        '''
        t_lists = np.linspace(0, T, 1000)
        def Omega_I(t, args):
            return args["omega_I"] * np.exp(-0.5 * ((t - args["t0"]) / args["sigma"]) ** 2)
        def Omega_Q(t, args):
            return args["omega_Q"] * np.exp(-0.5 * ((t - args["t0"]) / args["sigma"]) ** 2) * (-(t - args["t0"]) / (args["sigma"] ** 2 ))
        args1 = {"omega_I": np.pi / (sigma * np.sqrt(2 * np.pi)), "omega_Q": -np.pi / (sigma * np.sqrt(2 * np.pi) * self.qubit1.anharmonicity), "t0": T/2, "sigma": sigma}
        args2 = {"omega_I": np.pi / (sigma * np.sqrt(2 * np.pi)), "omega_Q": -np.pi / (sigma * np.sqrt(2 * np.pi) * self.qubit2.anharmonicity), "t0": T/2, "sigma": sigma}
        H_d = [[0.5 * (self.a1 + self.a1.dag()) , Omega_I] , [0.5 * (-1j) * (self.a1 - self.a1.dag()) , Omega_Q],
         [0.5 * (self.a2 + self.a2.dag()) , Omega_I] , [0.5 * (-1j) * (self.a2 - self.a2.dag()) , Omega_Q]]
        # 关闭耦合
        g_shut = [[0 for m in range(self.cavity.M)], [0 for m in range(self.cavity.M)]]
        self.control_g(g_shut)

        H = [self.H] + H_d
        psi0 = tensor(self.qubit1.state, self.qubit2.state, *[basis(self.cavity.n_levels[m], 0) for m in range(self.cavity.M)])
        result = sesolve(H, psi0, t_lists, args=args1)
        final_state = result.states[-1]
        target_state = tensor(basis(self.qubit1.n_levels, 1), basis(self.qubit2.n_levels, 1), *[basis(self.cavity.n_levels[m], 0) for m in range(self.cavity.M)])
        fidelity = np.abs(target_state.overlap(final_state)) ** 2
        return final_state, fidelity
    
    def simulate_iSWAP(self, g):
        '''
        模拟iSWAP门操作，计算门的有效矩阵和保真度
        两qubit共振，腔模远失谐
        '''
        def f(t,args):
            T, T1, T2 = args["T"], args["T1"], args["T2"]
            if t < 0 or t > T:
                return [[0 for m in range(self.cavity.M)], [0 for m in range(self.cavity.M)]]
            elif t < T1:
                return [[np.sin(0.5 * np.pi * t / T1) ** 2 for m in range(self.cavity.M)], [np.sin(0.5 * np.pi * t / T1) ** 2 for m in range(self.cavity.M)]]
            elif t < T - T2:
                return [[1 for m in range(self.cavity.M)], [1 for m in range(self.cavity.M)]]
            elif t <= T:
                return [[np.sin(0.5 * np.pi * (T - t) / T2) ** 2 for m in range(self.cavity.M)], [np.sin(0.5 * np.pi * (T - t) / T2) ** 2 for m in range(self.cavity.M)]]
        g_lists = [[g for m in range(self.cavity.M)], [g for m in range(self.cavity.M)]]
        self.control_g(g_lists, is_time_dependent=True, f=f)
        J = sum([self.g[0][m] * self.g[1][m] / np.abs(self.qubit1.frequency - self.cavity.frequencies[m]) for m in range(self.cavity.M)])
        T_ideal = np.pi / (2 * J) 
        args = {"T": T_ideal, "T1": T_ideal/10, "T2": T_ideal/10}
        U = propagator(self.H, T_ideal, args=args, options = {"nsteps":1000000})
        # 计算子空间的态矢量
        kets = [tensor(basis(self.qubit1.n_levels, i), basis(self.qubit2.n_levels, j), *[basis(self.cavity.n_levels[m], 0) for m in range(self.cavity.M)]) for i in range(2) for j in range(2)]
        U_c = np.zeros((4, 4), dtype=complex)
        for i, ket_i in enumerate(kets):
            for j, ket_j in enumerate(kets):
                U_c[i, j] = (ket_i.dag() * U * ket_j)

        U_eff = Qobj(U_c, dims=[[2,2], [2,2]])
        U_ideal = ideal_iSWAP()
        leakage = 1 - np.mean(np.sum(np.abs(U_eff.full())**2, axis=0))
        F = (abs((U_ideal.dag() @ U_eff).tr())**2 + 4) / (4 * 5)
        import scipy.linalg as la

        Uu, _ = la.polar(U_eff.full())   # U_eff ≈ Uu @ P
        Uu = Qobj(Uu, dims=U_eff.dims)
        
        def Rz(phi):
            return Qobj([[np.exp(-1j*phi/2), 0],
                        [0, np.exp(1j*phi/2)]])

        def avg_gate_fidelity(U, V):
            d = U.shape[0]
            tr = (V.dag() * U).tr()
            return (abs(tr)**2 + d) / (d*(d+1))

        U_id = ideal_iSWAP()
        from scipy.optimize import minimize

        def cost(x):
            a,b,c,d = x
            L = tensor(Rz(a), Rz(b))   # post Z
            R = tensor(Rz(c), Rz(d))   # pre Z
            Uc = L * Uu * R
            return -avg_gate_fidelity(Uc, U_id)

        res = minimize(cost, x0=np.zeros(4))
        F_loc = -res.fun
        a, b, c, d = res.x
        L_opt = tensor(Rz(a), Rz(b))
        R_opt = tensor(Rz(c), Rz(d))
        U_cal = L_opt * Uu * R_opt
        #print("Fidelity up to local Z:", F_loc)
        #print("best angles:", res.x)
        return U_cal, leakage, F_loc

# 下面定义一些两比特门函数，如CiSWAP门，CZ门等
def ideal_iSWAP():
    '''
    计算理想iSWAP门的矩阵
    '''
    iSWAP = Qobj([[1, 0, 0, 0],
                  [0, 0, 1j, 0],
                  [0, 1j, 0, 0],
                  [0, 0, 0, 1]], dims=[[2,2], [2,2]])
    return iSWAP
    
def simulate_iSWAP(qubit1, qubit2, g):
    '''
    模拟iSWAP门操作，计算门的有效矩阵和保真度
    :param qubit1: 进行门操作的第一个TransmonQubit对象
    :param qubit2: 进行门操作的第二个TransmonQubit对象
    :param g: 两个qubit之间的耦合强度（GHz）
    :param T: 门操作时间（ns）
    '''
    H_0 = tensor(qubit1.get_hamiltonian_rwa(qubit1.frequency),qeye(qubit2.n_levels)) + tensor(qeye(qubit1.n_levels),qubit2.get_hamiltonian_rwa(qubit2.frequency))
    H_int = g * (tensor(destroy(qubit1.n_levels), destroy(qubit2.n_levels).dag()) + tensor(destroy(qubit1.n_levels).dag(), destroy(qubit2.n_levels)))
    H = H_0+ H_int
    T = np.pi / (2 * g)  # iSWAP门的操作时间
    U = propagator(H, T)
    idx = [0, 1, qubit2.n_levels, qubit2.n_levels + 1]  # 选择前4个能级对应的子空间
    U_eff = Qobj(U.full()[np.ix_(idx, idx)], dims=[[2,2], [2,2]])
    U_eff = U_eff.dag()  # TODO：为什么
    U_ideal = ideal_iSWAP()
    leakage = 1 - np.mean(np.sum(np.abs(U_eff.full())**2, axis=0))
    F = (abs((U_ideal.dag() @ U_eff).tr())**2 + 4) / (4 * 5)
    return U_eff, leakage, F

def ideal_CZ():
    '''
    计算理想CZ门的矩阵
    '''
    CZ = Qobj([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, -1]], dims=[[2,2], [2,2]])
    return CZ

def simulate_CZ(qubit1, qubit2, g):
    '''
    模拟CZ门操作，计算门的有效矩阵和保真度
    :param qubit1: 进行门操作的第一个TransmonQubit对象
    :param qubit2: 进行门操作的第二个TransmonQubit对象
    :param g: 两个qubit之间的耦合强度（GHz）
    :param T: 门操作时间（ns）
    '''
    H_0 = tensor(qubit1.get_hamiltonian_rwa(qubit1.frequency),qeye(qubit2.n_levels)) + tensor(qeye(qubit1.n_levels),qubit2.get_hamiltonian_rwa(qubit2.frequency))
    H_int = g * (tensor(destroy(qubit1.n_levels), destroy(qubit2.n_levels).dag()) + tensor(destroy(qubit1.n_levels).dag(), destroy(qubit2.n_levels)))
    H_1 = H_0+ H_int
    def Delta_1(t, args):
        T, T1, T2 = args["T"], args["T1"], args["T2"]
        A = args["A"]
        if t < 0 or t > T:
            return 0
        elif t < T1:
            return A * np.sin(0.5 * np.pi * t / T1) ** 2
        elif t < T - T2:
            return A
        elif t <= T:
            return A * np.sin(0.5 * np.pi * (T - t) / T2) ** 2
    T = np.pi / (np.sqrt(2.0) * g)  # CZ门的操作时间
    T1 = T / 4
    T2 = T / 4
    A = -qubit2.anharmonicity
    args = {"T": T, "T1": T1, "T2": T2, "A": A}
    n_1 = tensor(qubit1.n, qeye(qubit2.n_levels))
    H = [H_1, [n_1, Delta_1]]
    U = propagator(H, T, args=args)
    idx = [0, 1, qubit2.n_levels, qubit2.n_levels + 1]  # 选择前4个能级对应的子空间
    U_eff = Qobj(U.full()[np.ix_(idx, idx)], dims=[[2,2], [2,2]])
    phi00 = np.angle(U_eff[0, 0])
    phi10 = np.angle(U_eff[1, 1])
    phi01 = np.angle(U_eff[2, 2])
    phi1 = phi10 - phi00
    phi2 = phi01 - phi00
    D = Qobj(np.diag([np.exp(-1j * phi00), np.exp(-1j * phi10), np.exp(-1j * phi01), np.exp(-1j * (phi10 + phi01 - phi00))]), dims=[[2,2], [2,2]])
    U_eff_corr = D @ U_eff
    U_ideal = ideal_CZ()
    leakage = 1 - 0.25 * (np.sum(np.abs(U_eff_corr.full())**2))
    F = (abs((U_ideal.dag() @ U_eff_corr).tr())**2 + 4) / (4 * 5)
    return U_eff_corr, leakage, F




