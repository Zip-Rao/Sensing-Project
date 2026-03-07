'''
该文件只作为参考，不参与仿真计算
'''

import numpy as np
import qutip as qt
from qutip import Qobj, basis, destroy, num, qeye
from scipy.linalg import eigh
from scipy.special import factorial

class TransmonBasisRepresentation:
    """Transmon在不同基下的表示"""
    
    def __init__(self, EC, EJ, n_levels=5, n_charge_range=10):
        self.EC = EC
        self.EJ = EJ
        self.n_levels = n_levels
        self.n_charge_range = n_charge_range
        
        # 基本参数
        self.omega_p = np.sqrt(8 * EJ * EC)  # 等离子体频率
        self.omega_q = self.omega_p - EC      # 量子比特频率
        self.alpha = -EC                     # 非谐性
        
    # ==================== 1. FOCK基 ====================
    def hamiltonian_fock_basis(self):
        """
        Fock基 {|n⟩, n=0,1,2,...}
        
        H = ω_p (a†a + 1/2) - EJ - (EC/2) a†a†aa
          ≈ ω_q a†a + (α/2) a†a†aa + 常数
        
        性质: 
        - a|n⟩ = √n|n-1⟩
        - a†|n⟩ = √(n+1)|n+1⟩
        - ⟨n|m⟩ = δ_{nm}
        - 无限维，计算时截断
        """
        a = destroy(self.n_levels)
        n = a.dag() * a
        
        # 完整哈密顿量
        H = (-self.EJ * qeye(self.n_levels) + 
             self.omega_p * (n + 0.5 * qeye(self.n_levels)) + 
             (self.alpha / 2) * (n * n - n))
        
        return H
    
    # ==================== 2. 电荷基 ====================
    def hamiltonian_charge_basis(self):
        """
        电荷基 {|n⟩, n = -N,...,N}
        
        H = 4EC Σ (n - ng)² |n⟩⟨n| - (EJ/2) Σ (|n⟩⟨n+1| + |n+1⟩⟨n|)
        
        性质:
        - n|n⟩ = n|n⟩ (Cooper对数算符本征态)
        - ⟨n|m⟩ = δ_{nm}
        - 周期边界: |n+N⟩ = |n⟩
        - 完全基，无限维但可截断
        """
        N = self.n_charge_range
        dim = 2 * N + 1
        n_values = np.arange(-N, N + 1)
        
        # 充电能量项
        H_charge = Qobj(np.zeros((dim, dim)))
        ng = 0.0  # 门电荷，通常设为0
        for i, n in enumerate(n_values):
            H_charge[i, i] = 4 * self.EC * (n - ng)**2
        
        # 约瑟夫森能量项
        H_josephson = Qobj(np.zeros((dim, dim)))
        for i in range(dim - 1):
            H_josephson[i, i+1] = -self.EJ / 2
            H_josephson[i+1, i] = -self.EJ / 2
        
        H = H_charge + H_josephson
        
        return H, n_values
    
    # ==================== 3. 能量本征基 ====================
    def hamiltonian_eigen_basis(self):
        """
        能量本征基 {|E_k⟩, k=0,1,2,...}
        
        H = Σ E_k |E_k⟩⟨E_k|
        
        性质:
        - H|E_k⟩ = E_k|E_k⟩
        - ⟨E_k|E_l⟩ = δ_{kl}
        - 物理测量基
        - 与电荷基通过U矩阵变换
        """
        # 从电荷基对角化得到
        H_charge, _ = self.hamiltonian_charge_basis()
        evals, evecs = H_charge.eigenstates()
        
        # 截断到低能级
        evals = evals[:self.n_levels] - evals[0]  # 以基态为0
        evecs = evecs[:self.n_levels]
        
        # 能量本征基下的哈密顿量（对角）
        H = Qobj(np.diag(evals))
        
        return H, evals, evecs
    
    # ==================== 4. 相位基 ====================
    def hamiltonian_phase_basis(self, n_phi=100):
        """
        相位基 {|φ⟩, φ ∈ [0, 2π)}
        
        |φ⟩ = (1/√(2π)) Σ_{n=-∞}^{∞} e^{inφ} |n⟩
        H = -4EC ∂²/∂φ² - EJ cos φ
        
        性质:
        - 连续基，⟨φ|φ'⟩ = δ(φ-φ')
        - 相位算符本征态: e^{iφ̂}|φ⟩ = e^{iφ}|φ⟩
        - 周期边界: |φ+2π⟩ = |φ⟩
        - 适用于大EJ/EC极限
        """
        # 离散化相位空间
        phi = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
        dphi = phi[1] - phi[0]
        
        # 动能项: -4EC ∂²/∂φ² (二阶导数矩阵)
        T = np.zeros((n_phi, n_phi))
        for i in range(n_phi):
            T[i, i] = -2
            T[i, (i+1)%n_phi] = 1
            T[i, (i-1)%n_phi] = 1
        T = -4 * self.EC * T / dphi**2
        
        # 势能项: -EJ cos φ
        V = np.diag(-self.EJ * np.cos(phi))
        
        H = Qobj(T + V)
        
        return H, phi
    
    # ==================== 5. Bloch基（两能级子空间） ====================
    def hamiltonian_bloch_basis(self):
        """
        Bloch基 {|0⟩, |1⟩} (两能级近似)
        
        H = (ω_q/2) σ_z
        
        性质:
        - σ_z|0⟩ = +|0⟩, σ_z|1⟩ = -|1⟩
        - 泡利算符完备
        - 适用于低能动力学
        - 忽略泄漏到|2⟩
        """
        # 标准Bloch表示
        H = (self.omega_q / 2) * qt.sigmaz()
        
        # 包含非谐修正的有效哈密顿量
        # H_eff = (ω_q/2) σ_z - (χ/2) σ_z ⊗ 其他自由度
        return H
    
    # ==================== 6. 坐标/动量基 ====================
    def hamiltonian_position_basis(self, x_range=(-5, 5), n_x=200):
        """
        坐标基 {|x⟩, x ∈ ℝ}
        
        ψ_n(x) = ⟨x|n⟩ = (1/√(2^n n!)) (mω/πħ)^{1/4} e^{-mωx²/2} H_n(√(mω/ħ)x)
        
        性质:
        - 连续基，⟨x|x'⟩ = δ(x-x')
        - 位置算符本征态: x̂|x⟩ = x|x⟩
        - 非紧致，需要截断
        """
        x = np.linspace(x_range[0], x_range[1], n_x)
        dx = x[1] - x[0]
        
        # 动能项: -ħ²/2m ∂²/∂x² (设ħ=1, m=1)
        T = np.zeros((n_x, n_x))
        for i in range(n_x):
            T[i, i] = -2
            if i > 0:
                T[i, i-1] = 1
            if i < n_x-1:
                T[i, i+1] = 1
        T = -0.5 * T / dx**2
        
        # 势能项: 1/2 ω_p² x² - EJ cos(φ_zpf x) 
        # φ_zpf = (2EC/EJ)^{1/4}
        phi_zpf = (2 * self.EC / self.EJ)**0.25
        V = np.diag(0.5 * self.omega_p**2 * x**2 - self.EJ * np.cos(phi_zpf * x))
        
        H = Qobj(T + V)
        
        return H, x
    


class BasisTransformation:
    """Transmon基变换矩阵"""
    
    def __init__(self, transmon):
        self.t = transmon
        
    # ==================== 变换矩阵 ====================
    
    def charge_to_fock_matrix(self):
        """
        电荷基 → Fock基 变换矩阵
        
        U_fc = ⟨n_fock|n_charge⟩
        
        Fock基是电荷基在EJ=0时的本征态
        对于大EJ/EC，两者通过Bogoliubov变换联系
        """
        N = self.t.n_charge_range
        dim_charge = 2 * N + 1
        dim_fock = self.t.n_levels
        
        U = np.zeros((dim_fock, dim_charge), dtype=complex)
        
        # 近似变换：无相互作用电荷基与Fock基的关系
        # |n_fock⟩ ≈ Σ_k c_{nk} |k_charge⟩
        for n in range(dim_fock):
            for k in range(-N, N+1):
                # 高斯波包近似
                idx = k + N
                sigma = (self.t.EJ / (32 * self.t.EC))**0.25
                U[n, idx] = np.exp(-(k - n)**2 / (2 * sigma**2))
        
        # 正交归一化
        U, _ = np.linalg.qr(U)
        
        return Qobj(U)
    
    def fock_to_eigen_matrix(self):
        """
        Fock基 → 能量本征基 变换矩阵
        
        U_fe = ⟨E_k|n⟩
        
        对于Transmon，Fock基不是能量本征态
        需要通过对角化H_fock得到
        """
        H_fock = self.t.hamiltonian_fock_basis()
        evals, evecs = H_fock.eigenstates()
        
        # U_fe的列是能量本征态在Fock基下的表示
        dim = self.t.n_levels
        U = Qobj(np.zeros((dim, dim), dtype=complex))
        
        for i, evec in enumerate(evecs[:dim]):
            U[:, i] = evec.full().flatten()
        
        return U, evals[:dim]
    
    def charge_to_eigen_matrix(self):
        """
        电荷基 → 能量本征基 变换矩阵
        
        U_ce = ⟨E_k|n_charge⟩
        
        最直接的变换：对角化电荷基哈密顿量
        """
        H_charge, _ = self.t.hamiltonian_charge_basis()
        evals, evecs = H_charge.eigenstates()
        
        # 截断到低能级
        dim_charge = 2 * self.t.n_charge_range + 1
        dim_fock = self.t.n_levels
        
        U = Qobj(np.zeros((dim_fock, dim_charge), dtype=complex))
        
        for i in range(dim_fock):
            U[i, :] = evecs[i].full().flatten()
        
        return U, evals[:dim_fock]
    
    def bloch_to_fock_matrix(self):
        """
        Bloch基 → Fock基 变换矩阵
        
        U_bf = ⟨n|0/1⟩
        
        |0⟩ ≈ |n=0⟩
        |1⟩ ≈ |n=1⟩
        """
        dim = self.t.n_levels
        
        U = Qobj(np.zeros((dim, 2), dtype=complex))
        
        # 零级近似：Bloch基就是低能Fock态
        U[0, 0] = 1.0  # |0⟩_B = |0⟩_F
        U[1, 1] = 1.0  # |1⟩_B = |1⟩_F
        
        # 更精确的变换包含高阶修正
        # |0⟩ ≈ |0⟩ + ε₂|2⟩ + ...
        # |1⟩ ≈ |1⟩ + ε₃|3⟩ + ...
        if dim > 2:
            # 非谐性引起的修正
            epsilon_2 = -self.t.alpha / (2 * np.sqrt(2) * self.t.omega_q)
            epsilon_3 = -np.sqrt(3) * self.t.alpha / (2 * self.t.omega_q)
            
            U[2, 0] = epsilon_2
            U[3, 1] = epsilon_3
        
        # 归一化
        U = U.unit()
        
        return U
    
    # ==================== 变换函数 ====================
    
    def transform_operator(self, operator, from_basis, to_basis):
        """
        在不同基之间变换算符
        
        A' = U† A U
        """
        # 获取变换矩阵
        U = self.get_transformation_matrix(from_basis, to_basis)
        
        if U is None:
            # 尝试通过中间基变换
            intermediate = 'fock'
            U1 = self.get_transformation_matrix(from_basis, intermediate)
            U2 = self.get_transformation_matrix(intermediate, to_basis)
            
            if U1 is not None and U2 is not None:
                U = U2 * U1
            else:
                raise ValueError(f"Cannot transform from {from_basis} to {to_basis}")
        
        return U.dag() * operator * U
    
    def get_transformation_matrix(self, from_basis, to_basis):
        """获取两个基之间的变换矩阵"""
        
        # 直接变换
        key = f"{from_basis}_to_{to_basis}"
        
        transformations = {
            'charge_to_fock': self.charge_to_fock_matrix(),
            'fock_to_eigen': self.fock_to_eigen_matrix()[0],
            'charge_to_eigen': self.charge_to_eigen_matrix()[0],
            'fock_to_bloch': self.bloch_to_fock_matrix().dag(),
            'bloch_to_fock': self.bloch_to_fock_matrix(),
        }
        
        if key in transformations:
            return transformations[key]
        
        # 逆变换
        inv_key = f"{to_basis}_to_{from_basis}"
        if inv_key in transformations:
            return transformations[inv_key].dag()
        
        return None