# 基于全密度矩阵模拟的数值反演方案

## 概述

本方案旨在实现一个基于全密度矩阵模拟的数值反演算法，用于从量子测量数据中高精度重建瞬态磁场信号。该算法针对非线性响应较强的Transmon量子比特系统，克服传统线性反卷积方法的局限性，提供更准确的重建结果。

## 算法原理

### 问题描述

给定实验测量数据：激发态概率曲线 $p_{\text{meas}}(t_i), i=1,...,N_t$（$t_i$表示等效时间采样时测量的延迟时间），需要重建原始磁场信号 $B(t)$，使得该磁场信号的模拟结果与测量数据尽可能接近。系统响应由量子比特的密度矩阵演化描述，包含非线性特性。

### 数学模型

将磁场信号参数化为基函数的线性组合：

$$
B(t) = \sum_{k=1}^{M} b_k \phi_k(t)
$$

其中 $\{\phi_k(t)\}$ 为预先选定的基函数集（如B样条、傅里叶基等），$\mathbf{b} = [b_1, ..., b_M]^T$ 为待求系数。

### 前向模型

对于给定的系数向量 $\mathbf{b}$，通过求解Lindblad主方程计算模拟测量概率：

$$
\frac{d\rho}{dt} = -i[H(t;\mathbf{b}), \rho] + \mathcal{L}[\rho]
$$

其中哈密顿量 $H(t;\mathbf{b}) = H_0 + H_{\text{pulse}}(t) + H_{\text{field}}(B(t;\mathbf{b}))$，$\mathcal{L}$ 为耗散超算符。模拟概率 $p_{\text{sim}}(t_i;\mathbf{b}) = \text{Tr}[\rho(t_i) \cdot |e\rangle\langle e|]$。

### 优化问题

最小化残差：

$$
\min_{\mathbf{b}} \left\| \mathbf{p}_{\text{meas}} - \mathbf{p}_{\text{sim}}(\mathbf{b}) \right\|^2 + \lambda \mathbf{b}^T \mathbf{D}^T \mathbf{D} \mathbf{b}
$$

其中 $\lambda$ 为正则化参数，$\mathbf{D}$ 为平滑算子（通常为二阶差分矩阵）。

### Levenberg-Marquardt算法

采用带正则化的Levenberg-Marquardt算法迭代求解：

1. **初始化**：使用Hammerstein-Wiener逐块反演结果作为初始估计 $\mathbf{b}^{(0)}$

2. **迭代步骤**（$n=0,1,2,...$）：

   a. 前向模拟：计算 $p_{\text{sim}}(t_i;\mathbf{b}^{(n)})$ 和残差 $\mathbf{r}^{(n)} = \mathbf{p}_{\text{meas}} - \mathbf{p}_{\text{sim}}(\mathbf{b}^{(n)})$

   b. 计算Jacobian矩阵：$J_{ik}^{(n)} = \frac{\partial p_{\text{sim}}(t_i)}{\partial b_k} \bigg|_{\mathbf{b}=\mathbf{b}^{(n)}}$

   c. 更新方程：$(\mathbf{J}^{(n)T}\mathbf{J}^{(n)} + \mu^{(n)}\mathbf{I} + \lambda\mathbf{D}^T\mathbf{D}) \delta\mathbf{b} = \mathbf{J}^{(n)T} \mathbf{r}^{(n)}$

   d. 参数更新：$\mathbf{b}^{(n+1)} = \mathbf{b}^{(n)} + \delta\mathbf{b}$
   
   e. 阻尼调整：若 $\|\mathbf{r}^{(n+1)}\| < \|\mathbf{r}^{(n)}\|$，则 $\mu^{(n+1)} = \mu^{(n)}/2$；否则 $\mu^{(n+1)} = 2\mu^{(n)}$ 并拒绝更新

3. **收敛判断**：当 $\|\delta\mathbf{b}\|/\|\mathbf{b}\| < \epsilon_{\text{tol}}$ 或 $\|\mathbf{r}\| < \epsilon_{\text{data}}$ 时停止

### Jacobian计算

采用伴随法（adjoint method）高效计算Jacobian矩阵$$J_{ik} = \frac{\partial p_{\text{sim}}(t_i)}{\partial b_k}$$

1. 前向传播：密度矩阵 $\rho(t)$满足$$\frac{d\rho}{dt} = \mathcal{L}[\rho]=-i[H(t;\mathbf{b}), \rho] + \sum_k \gamma_k(L_k\rho L_k^\dagger - \frac{1}{2} \{L_k^\dagger L_k, \rho\})$$
初始条件为 $\rho(0) = |0\rangle\langle 0|$，观测量为$p(t_i) = \text{Tr}[\rho(t_i) \cdot |e\rangle\langle e|]$。
2. 反向传播：为了得到雅可比矩阵，需要计算方程$$\frac{\partial p(t_i)}{\partial b_k} = \text{Tr} [\frac{\partial \rho(t)}{\partial b_k} \ket{1}\bra{1}] = \text{Tr}[\sigma_{k}(t_i)\ket{1}\bra{1}]$$

$\sigma_k(t)$满足方程$$ \frac{d\sigma_k}{dt} = \frac{\partial }{\partial b_k}\frac{\partial }{\partial t}\rho = \frac{\partial \mathcal{L}}{\partial b_k}[\rho] + \mathcal{L}[\sigma_k] = \mathcal{L}[\sigma_k] + \mathcal{S_k}[\rho] $$
其中，$\mathcal{S_k}[\rho]$为源项，当耗散超算符与$b_k$无关时，$\mathcal{S_k}[\rho] = -i[\frac{\partial H(t)}{\partial b_k},\rho]$。

$\sigma_k$满足$\sigma_{k}(0) = 0$，利用这个性质，
引入函数$\lambda(t)$，有$$\begin{aligned}\text{Tr}[\ket{1}\bra{1}\sigma_{k}(t_i)] &= \text{Tr}[\lambda(t_i)\sigma_k(t_i)] - \text{Tr}[\lambda(0)\sigma_k(0)]\\ &= \int_0^{t_i} \frac{d}{dt}\bigg\{\text{Tr}\left[\lambda(t) \sigma_k(t)\right]\bigg\} dt \\ &= \int_0^{t_i} \text{Tr}\bigg[\dot\lambda\sigma_k + \lambda \dot\sigma_k \bigg]dt \\ &= \int_0^{t_i} \text{Tr}(\dot\lambda\sigma_{k} ) + \text{Tr}(\lambda\mathcal{L[\sigma_k]}) -i \ \text{Tr}(\lambda[\frac{\partial H}{\partial b_k},\rho])dt \\ &= \int_0^{t_i} \text{Tr}\left[(\dot\lambda(t) + \mathcal{L}^\dagger[\lambda])\sigma_k\right] -i\  \text{Tr}(\lambda[\frac{\partial H}{\partial b_k},\rho])dt \\ &=-i\  \text{Tr}(\lambda[\frac{\partial H}{\partial b_k},\rho])dt \end{aligned}$$
其中，$\lambda(t) $满足：$$\frac{d\lambda}{dt} = -\mathcal{L}^\dagger[\lambda], \quad \lambda(t_i) = |e\rangle\langle e|$$
上述推导利用了超算符的性质：$$\text{Tr}[\lambda \mathcal{L}[\sigma_k]] = \text{Tr}[\mathcal{L}^\dagger[\lambda]\sigma_k]$$

于是得到了完整的伴随方程：$$\dot\lambda = -i[H, \lambda] + \sum_k \gamma_k(L_k^\dagger \lambda L_k - \frac{1}{2} \{L_k^\dagger L_k, \lambda\})$$

作变量替换$s = t_i - t, \mu(s) = \lambda(t) =  \lambda(t_i - s) $，则$$\frac{d\mu}{ds} = i[H(t_i - s), \mu] + \sum_k \gamma_k(L_k^\dagger \mu L_k - \frac{1}{2} \{L_k^\dagger L_k, \mu\}), \quad \mu(0) = |e\rangle\langle e|$$


3. 梯度计算：$\frac{\partial p_{\text{sim}}(t_i)}{\partial b_k} = \text{Re} \bigg\{-i\int_0^{t_i} \text{Tr}\left[\lambda(t) [\frac{\partial H(t)}{\partial b_k}, \rho(t)]\right] dt\bigg\}$

旋转坐标系下，Transmon qubit的哈密顿量为$$ H(t, b) = \frac{\Delta(B)}{2}G + H_{drive}(t)$$
其中，$G$为z方向的旋转算符，二能级时为$\sigma_z$。

代入$B(t) = \sum_k b_k \phi_k(t)$，则$$\frac{\partial H}{\partial b_k} = \frac{\partial \Delta(B)}{\partial B}\phi_k(t)\frac{G}{2}$$
若脉冲频率与磁场无关，则上式第一项即Transmon qubit对磁场的灵敏度。

## 集成方案

### 扩展现有Analysis类

在 `src/analysis.py` 中添加新的方法类，与现有Wiener反卷积和Hammerstein-Wiener反卷积方法保持兼容。

#### 新增类：NumericalInversion

```python
class NumericalInversion:
    """基于全密度矩阵模拟的数值反演算法"""

    def __init__(self, qubit, control_pulse, basis_functions, reg_param=1e-3):
        """
        初始化数值反演器

        参数：
            qubit: TransmonQubit对象
            control_pulse: CompositePulse对象，控制脉冲序列
            basis_functions: 基函数列表，每个元素为可调用函数 phi_k(t)
            reg_param: 正则化参数 λ
        """
        self.qubit = qubit
        self.control_pulse = control_pulse
        self.basis_functions = basis_functions
        self.M = len(basis_functions)
        self.reg_param = reg_param

        # 构建平滑算子 D（二阶差分）
        self.D = self._build_smoothing_operator()

    def _build_smoothing_operator(self):
        """构建二阶差分平滑算子"""
        D = np.zeros((self.M-2, self.M))
        for i in range(self.M-2):
            D[i, i] = 1
            D[i, i+1] = -2
            D[i, i+2] = 1
        return D

    def b_to_B(self, b):
        """将系数向量转换为磁场信号"""
        t_array = np.linspace(0, self.control_pulse.t_list[-1], 1000)
        B = np.zeros_like(t_array)
        for k, phi_k in enumerate(self.basis_functions):
            B += b[k] * phi_k(t_array)
        return t_array, B

    def forward_simulation(self, b, t_points=None):
        """
        前向模拟：计算给定系数b下的激发态概率

        参数：
            b: 系数向量 (M,)
            t_points: 时间点数组，None则使用控制脉冲时间点

        返回：
            t_sim: 时间点数组
            p_sim: 模拟概率数组
        """
        if t_points is None:
            t_points = self.control_pulse.t_list

        # 将b转换为磁场信号B(t)
        t_array, B_array = self.b_to_B(b)

        # 构建含磁场信号的哈密顿量
        H_total = self._build_hamiltonian_with_field(B_array, t_array)

        # 求解主方程
        result = mesolve(H_total, self.qubit.state, t_points,
                         self.qubit.get_collapse_operators(),
                         e_ops=[self.qubit.excited_state_projector()])

        return t_points, np.array(result.expect[0])

    def _build_hamiltonian_with_field(self, B_array, t_array):
        """构建包含磁场信号的哈密顿量"""
        # 插值函数将离散B(t)转换为连续函数
        from scipy.interpolate import interp1d
        B_func = interp1d(t_array, B_array, kind='cubic',
                          bounds_error=False, fill_value=0.0)

        def H_total(t, args):
            # 基础哈密顿量
            if self.control_pulse.frame == 0:
                H_0 = self.qubit.get_hamiltonian(self.qubit.frequency)
            else:
                H_0 = self.qubit.get_hamiltonian_rwa(self.control_pulse.omega_d)

            # 脉冲哈密顿量
            H_pulse = self.control_pulse.get_hamiltonian(t)

            # 磁场项：通过改变磁通实现
            B_t = B_func(t)
            flux_t = self.qubit.flux + B_t  # 假设B(t)直接加在磁通上
            qubit_t = TransmonQubit(self.qubit.EC, self.qubit.EJ_0,
                                    self.qubit.T1, self.qubit.T2,
                                    flux=flux_t,
                                    state=self.qubit.state,
                                    n_levels=self.qubit.n_levels)

            if self.control_pulse.frame == 0:
                H_field = qubit_t.get_hamiltonian(self.qubit.frequency) - H_0
            else:
                H_field = qubit_t.get_hamiltonian_rwa(self.control_pulse.omega_d) - H_0

            return H_0 + H_pulse + H_field

        return H_total

    def compute_jacobian_adjoint(self, b, t_points):
        """
        使用伴随法计算Jacobian矩阵

        参数：
            b: 当前系数向量
            t_points: 时间点数组

        返回：
            J: Jacobian矩阵 (N_t × M)
        """
        N_t = len(t_points)
        M = self.M
        J = np.zeros((N_t, M))

        # 前向模拟得到密度矩阵轨迹
        t_array, B_array = self.b_to_B(b)
        H_total = self._build_hamiltonian_with_field(B_array, t_array)

        # 求解主方程并保存状态轨迹
        result = mesolve(H_total, self.qubit.state, t_points,
                         self.qubit.get_collapse_operators(),
                         options=dict(store_states=True))

        states = result.states

        # 对每个时间点计算梯度
        for i, t_i in enumerate(t_points):
            # 构造伴随方程的初始条件
            rho_i = states[i]
            P_e = self.qubit.excited_state_projector()

            # 简化计算：使用有限差分近似Jacobian
            # 实际实现应使用真正的伴随法或自动微分

        return J

    def levenberg_marquardt(self, p_meas, t_meas, b_init=None, max_iter=50,
                           tol=1e-6, mu_init=1.0):
        """
        Levenberg-Marquardt优化主循环

        参数：
            p_meas: 测量概率数组 (N_t,)
            t_meas: 测量时间点数组 (N_t,)
            b_init: 初始系数向量，None则使用Hammerstein-Wiener结果
            max_iter: 最大迭代次数
            tol: 收敛容差
            mu_init: 初始阻尼参数

        返回：
            b_opt: 优化后的系数向量
            B_rec: 重建的磁场信号 (t_array, B_array)
            history: 优化历史记录
        """
        if b_init is None:
            # 使用Hammerstein-Wiener反卷积作为初始估计
            from .analysis import Analysis
            analysis = Analysis()
            dt = t_meas[1] - t_meas[0]
            t_kernel, kernel = analysis.get_kernel(self.control_pulse, self.qubit)
            omega_list, omega = analysis.wiener_deconvolution(p_meas, kernel, dt, self.reg_param)
            b_init = self._frequency_to_coefficients(omega_list, omega)

        b = b_init.copy()
        mu = mu_init
        N_t = len(p_meas)

        history = {
            'b': [b.copy()],
            'residual': [],
            'mu': [mu]
        }

        for iter in range(max_iter):
            # 前向模拟
            t_sim, p_sim = self.forward_simulation(b, t_meas)

            # 计算残差
            r = p_meas - p_sim
            residual_norm = np.linalg.norm(r)
            history['residual'].append(residual_norm)

            # 计算Jacobian（简化版：使用有限差分）
            J = self._finite_difference_jacobian(b, t_meas, p_sim)

            # 构建正规方程
            JtJ = J.T @ J
            A = JtJ + mu * np.eye(self.M) + self.reg_param * self.D.T @ self.D

            # 求解更新量
            delta_b = np.linalg.solve(A, J.T @ r)

            # 试验更新
            b_trial = b + delta_b
            _, p_trial = self.forward_simulation(b_trial, t_meas)
            r_trial = p_meas - p_trial
            residual_trial = np.linalg.norm(r_trial)

            # 判断是否接受更新
            if residual_trial < residual_norm:
                # 接受更新
                b = b_trial
                mu = max(mu/2, 1e-8)
                history['b'].append(b.copy())
                history['mu'].append(mu)

                # 收敛检查
                if np.linalg.norm(delta_b) / (np.linalg.norm(b) + 1e-8) < tol:
                    print(f"收敛于迭代 {iter+1}")
                    break
            else:
                # 拒绝更新，增加阻尼
                mu = min(mu * 2, 1e8)
                history['mu'].append(mu)

        # 最终重建
        t_array, B_array = self.b_to_B(b)

        return b, (t_array, B_array), history

    def _finite_difference_jacobian(self, b, t_points, p_sim, epsilon=1e-6):
        """使用有限差分法计算Jacobian矩阵（简化实现）"""
        N_t = len(t_points)
        M = self.M
        J = np.zeros((N_t, M))

        for k in range(M):
            b_plus = b.copy()
            b_plus[k] += epsilon
            _, p_plus = self.forward_simulation(b_plus, t_points)

            J[:, k] = (p_plus - p_sim) / epsilon

        return J

    def _frequency_to_coefficients(self, t_omega, omega):
        """将频率响应转换为基函数系数"""
        # 将频率响应投影到基函数空间
        from scipy.interpolate import interp1d
        from scipy.integrate import simpson

        # 创建插值函数
        omega_func = interp1d(t_omega, omega, kind='cubic',
                             bounds_error=False, fill_value=0.0)

        # 时间网格
        t_min, t_max = t_omega[0], t_omega[-1]
        t_grid = np.linspace(t_min, t_max, 1000)

        # 计算投影系数
        b = np.zeros(self.M)
        for k, phi_k in enumerate(self.basis_functions):
            integrand = omega_func(t_grid) * phi_k(t_grid)
            b[k] = simpson(integrand, t_grid)

        return b
```

#### 在Analysis类中添加接口方法

在现有 `Analysis` 类中添加一个包装方法，便于使用：

```python
class Analysis:
    # 现有方法...

    def numerical_field_inversion(self, p_meas, t_meas, qubit, control_pulse,
                                  basis_type='bspline', n_basis=20,
                                  reg_param=1e-3, max_iter=50, tol=1e-6):
        """
        数值场反演主接口

        参数：
            p_meas: 测量概率数组
            t_meas: 测量时间点数组
            qubit: TransmonQubit对象
            control_pulse: CompositePulse对象
            basis_type: 基函数类型 ('bspline', 'fourier', 'legendre')
            n_basis: 基函数数量
            reg_param: 正则化参数
            max_iter: 最大迭代次数
            tol: 收敛容差

        返回：
            t_rec: 重建时间数组
            B_rec: 重建磁场数组
            info: 反演信息字典
        """
        # 生成基函数
        if basis_type == 'bspline':
            basis_funcs = self._generate_bspline_basis(t_meas[0], t_meas[-1], n_basis, degree=3)
        elif basis_type == 'fourier':
            basis_funcs = self._generate_fourier_basis(t_meas[0], t_meas[-1], n_basis)
        elif basis_type == 'legendre':
            basis_funcs = self._generate_legendre_basis(t_meas[0], t_meas[-1], n_basis)
        else:
            raise ValueError(f"未知的基函数类型: {basis_type}")

        # 创建数值反演器
        inverter = NumericalInversion(qubit, control_pulse, basis_funcs, reg_param)

        # 运行优化
        b_opt, (t_rec, B_rec), history = inverter.levenberg_marquardt(
            p_meas, t_meas, max_iter=max_iter, tol=tol)

        # 收集反演信息
        info = {
            'coefficients': b_opt,
            'residual_history': history['residual'],
            'damping_history': history['mu'],
            'basis_type': basis_type,
            'n_basis': n_basis,
            'final_residual': history['residual'][-1] if history['residual'] else None
        }

        return t_rec, B_rec, info

    def _generate_bspline_basis(self, t_min, t_max, n_basis, degree=3):
        """生成B样条基函数"""
        from scipy.interpolate import BSpline
        import numpy as np

        # 均匀节点
        knots = np.linspace(t_min, t_max, n_basis - degree + 1)
        knots = np.r_[[t_min]*degree, knots, [t_max]*degree]  # 添加边界重复节点

        basis_funcs = []
        for i in range(n_basis):
            coeffs = np.zeros(n_basis)
            coeffs[i] = 1.0
            spline = BSpline(knots, coeffs, degree)
            basis_funcs.append(spline)

        return basis_funcs

    def _generate_fourier_basis(self, t_min, t_max, n_basis):
        """生成傅里叶基函数"""
        T = t_max - t_min

        basis_funcs = []
        # 常数项
        basis_funcs.append(lambda t: np.ones_like(t) / np.sqrt(T))

        # 正弦和余弦项
        n_max = (n_basis - 1) // 2
        for n in range(1, n_max + 1):
            omega_n = 2 * np.pi * n / T
            basis_funcs.append(lambda t, n=n: np.sqrt(2/T) * np.sin(omega_n * (t - t_min)))
            basis_funcs.append(lambda t, n=n: np.sqrt(2/T) * np.cos(omega_n * (t - t_min)))

        return basis_funcs[:n_basis]  # 截断到指定数量
```

### 基函数选择策略

1. **B样条基**（默认推荐）：
   - 优点：局部支持，数值稳定，适合非周期信号
   - 参数：3次样条，节点数 = n_basis - degree + 1

2. **傅里叶基**：
   - 优点：频域稀疏，适合周期性或准周期性信号
   - 缺点：全局支持，可能产生吉布斯现象

3. **勒让德多项式基**：
   - 优点：正交性良好，数值稳定
   - 缺点：全局支持，不适合局部特征

## 使用示例

### 基本使用流程

```python
import numpy as np
import matplotlib.pyplot as plt
from src.qubit import TransmonQubit
from src.pulse import create_ramsey_pulse
from src.analysis import Analysis

# 1. 准备数据
qubit = TransmonQubit(EC=0.2*2*np.pi, EJ=10.0*2*np.pi,
                      T1=100e3, T2=50e3, flux=0.0)
t_rabi = np.linspace(0, 40, 100)
control_pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

# 2. 获取测量数据（示例：模拟测量）
analysis = Analysis()
t_kernel, kernel = analysis.get_kernel(control_pulse, qubit)

# 模拟一个测试磁场信号
t_meas = np.linspace(0, 200, 400)
B_true = 0.01 * np.exp(-((t_meas - 100) / 20)**2) * np.sin(2*np.pi*0.01*t_meas)

# 通过前向模型生成测量数据（模拟实验）
# 这里简化处理，实际应从实验获取p_meas

# 3. 运行数值反演
t_rec, B_rec, info = analysis.numerical_field_inversion(
    p_meas=p_meas,  # 实验测量概率
    t_meas=t_meas,  # 测量时间点
    qubit=qubit,
    control_pulse=control_pulse,
    basis_type='bspline',
    n_basis=30,
    reg_param=1e-3,
    max_iter=100,
    tol=1e-6
)

# 4. 可视化结果
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(t_meas, p_meas, 'b-', label='Measured probability')
plt.xlabel('Time (ns)')
plt.ylabel('$p_e$')
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(t_rec, B_rec, 'r-', linewidth=2, label='Reconstructed')
plt.plot(t_meas, B_true, 'b--', alpha=0.7, label='True')
plt.xlabel('Time (ns)')
plt.ylabel('$B(t)$')
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 3)
plt.semilogy(info['residual_history'], 'o-')
plt.xlabel('Iteration')
plt.ylabel('Residual norm')
plt.title('Convergence history')
plt.grid(True)

plt.subplot(2, 2, 4)
basis_coeff = info['coefficients']
plt.stem(range(len(basis_coeff)), basis_coeff)
plt.xlabel('Basis index')
plt.ylabel('Coefficient value')
plt.title('Basis coefficients')
plt.grid(True)

plt.tight_layout()
plt.show()

# 5. 评估重建质量
from scipy.interpolate import interp1d

# 插值到相同时间网格
B_rec_interp = interp1d(t_rec, B_rec, kind='cubic')(t_meas)

# 计算相对误差
relative_error = np.linalg.norm(B_rec_interp - B_true) / np.linalg.norm(B_true)
print(f"Relative reconstruction error: {relative_error:.2%}")

# 计算信噪比改进
SNR_improvement = 10 * np.log10(np.var(B_true) / np.var(B_rec_interp - B_true))
print(f"SNR improvement: {SNR_improvement:.2f} dB")
```

### 与现有方法的对比使用

```python
# 比较不同反演方法
methods = ['wiener', 'hammerstein_wiener', 'numerical']
results = {}

# Wiener反卷积
dt = t_meas[1] - t_meas[0]
t_wiener, B_wiener = analysis.wiener_deconvolution(p_meas, kernel, dt, lambdas=0.1)
results['wiener'] = (t_wiener, B_wiener)

# Hammerstein-Wiener反卷积
t_hw, B_hw = analysis.hammerstein_wiener_deconvolution(qubit, p_meas, kernel, dt, lambdas=0.1)
results['hammerstein_wiener'] = (t_hw, B_hw)

# 数值反演
t_num, B_num, info = analysis.numerical_field_inversion(
    p_meas, t_meas, qubit, control_pulse,
    basis_type='bspline', n_basis=25, reg_param=1e-3
)
results['numerical'] = (t_num, B_num)

# 可视化对比
plt.figure(figsize=(10, 6))
colors = {'wiener': 'blue', 'hammerstein_wiener': 'green', 'numerical': 'red'}

for method, (t_method, B_method) in results.items():
    plt.plot(t_method, B_method, color=colors[method],
             label=f'{method}', linewidth=2, alpha=0.8)

plt.plot(t_meas, B_true, 'k--', linewidth=3, alpha=0.5, label='True')
plt.xlabel('Time (ns)')
plt.ylabel('$B(t)$')
plt.title('Comparison of inversion methods')
plt.legend()
plt.grid(True)
plt.show()
```

## 性能优化建议

### 计算加速

1. **Jacobian计算优化**：
   - 实现真正的伴随法，避免有限差分的数值误差和计算成本
   - 使用自动微分（如JAX）计算精确梯度

2. **并行化**：
   - 不同时间点的前向模拟可并行计算
   - 基函数维度的Jacobian列可并行计算

3. **稀疏性利用**：
   - 对于局部基函数（如B样条），Jacobian矩阵是稀疏的
   - 使用稀疏矩阵运算减少内存和计算量

### 内存管理

1. **状态轨迹存储**：
   - 对于长时间演化，只存储必要的时间点状态
   - 使用检查点技术减少内存占用

2. **矩阵预计算**：
   - 预先计算基函数在时间网格上的值
   - 缓存哈密顿量构建中的不变部分

## 误差分析与调试

### 常见问题及解决方案

1. **不收敛**：
   - 检查初始估计质量：使用Hammerstein-Wiener结果作为warm start
   - 调整阻尼参数 μ：初始值设为1.0，自适应调整
   - 增加正则化强度 λ：抑制过拟合

2. **过拟合**：
   - 增加正则化参数 λ
   - 减少基函数数量 n_basis
   - 使用交叉验证选择超参数

3. **计算时间过长**：
   - 减少时间点数量（适当降采样）
   - 使用较少的基函数
   - 启用并行计算

### 诊断工具

在 `NumericalInversion` 类中添加诊断方法：

```python
def diagnose_convergence(self, history):
    """诊断收敛问题"""
    residuals = history['residual']
    damping = history['mu']

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # 残差历史
    axes[0, 0].semilogy(residuals, 'o-')
    axes[0, 0].set_xlabel('Iteration')
    axes[0, 0].set_ylabel('Residual norm')
    axes[0, 0].grid(True)

    # 阻尼参数历史
    axes[0, 1].semilogy(damping, 's-')
    axes[0, 1].set_xlabel('Iteration')
    axes[0, 1].set_ylabel('Damping μ')
    axes[0, 1].grid(True)

    # 系数变化
    b_history = history['b']
    n_iter = len(b_history)
    b_array = np.array(b_history)

    for k in range(min(5, self.M)):
        axes[1, 0].plot(range(n_iter), b_array[:, k], label=f'b_{k}')
    axes[1, 0].set_xlabel('Iteration')
    axes[1, 0].set_ylabel('Coefficient value')
    axes[1, 0].legend()
    axes[1, 0].grid(True)

    # 条件数估计
    # 计算最后一次迭代的Jacobian条件数
    if 'J' in history:
        J = history['J'][-1]
        cond_number = np.linalg.cond(J.T @ J + self.reg_param * self.D.T @ self.D)
        axes[1, 1].text(0.1, 0.5, f'Condition number: {cond_number:.2e}')

    plt.tight_layout()
    return fig
```

## 实验集成指南

### 与现有实验流程对接

1. **数据预处理**：
   ```python
   # 从实验文件加载数据
   experimental_data = np.loadtxt('experiment_p_meas.csv')
   t_meas = experimental_data[:, 0]
   p_meas = experimental_data[:, 1]

   # 去噪和归一化
   from scipy.signal import savgol_filter
   p_meas_smooth = savgol_filter(p_meas, window_length=11, polyorder=3)

   # 移除基线漂移
   p_meas_corrected = p_meas_smooth - np.mean(p_meas_smooth[:10])
   ```

2. **参数校准**：
   ```python
   # 使用已知测试信号校准正则化参数
   calibration_signals = [
       ('gaussian', 0.01, 20, 100),  # 类型, 幅度, 宽度, 中心
       ('sine', 0.005, 0.01, 0),     # 类型, 幅度, 频率, 相位
   ]

   # 网格搜索最优参数
   param_grid = {
       'n_basis': [15, 20, 25, 30],
       'reg_param': [1e-4, 1e-3, 1e-2, 1e-1],
       'basis_type': ['bspline', 'fourier']
   }
   ```

3. **批量处理**：
   ```python
   def batch_inversion(data_files, qubit_params, pulse_params):
       """批量处理多个实验数据文件"""
       results = {}

       for file in data_files:
           # 加载数据
           data = load_experimental_data(file)

           # 创建量子比特和脉冲对象
           qubit = TransmonQubit(**qubit_params)
           control_pulse = create_ramsey_pulse(**pulse_params)

           # 运行反演
           t_rec, B_rec, info = analysis.numerical_field_inversion(
               data['p_meas'], data['t_meas'],
               qubit, control_pulse,
               basis_type='bspline', n_basis=25, reg_param=1e-3
           )

           results[file] = {
               't': t_rec,
               'B': B_rec,
               'info': info
           }

       return results
   ```

## 扩展与未来发展

### 算法扩展方向

1. **非线性基函数**：
   - 引入非线性参数化，如神经网络表示 $B(t)$
   - 使用物理信息神经网络（PINN）结合物理约束

2. **不确定性量化**：
   - 贝叶斯反演框架，提供后验分布
   - 马尔可夫链蒙特卡洛（MCMC）采样

3. **实时反演**：
   - 递归最小二乘实现在线更新
   - 滑动窗口处理流数据

### 硬件加速

1. **GPU加速**：
   - 使用CuPy或PyTorch实现矩阵运算
   - 并行化多个时间点的密度矩阵演化

2. **量子硬件辅助**：
   - 在量子处理器上直接模拟前向模型
   - 量子梯度计算加速优化

## 结论

本方案提出的基于全密度矩阵模拟的数值反演方法，为Transmon量子比特系统的瞬态磁场测量提供了高精度重建工具。通过与现有代码模块的无缝集成，用户可以：

1. 在处理强非线性响应时获得比传统方法更准确的结果
2. 通过基函数参数化灵活适应不同信号特征
3. 利用正则化和优化技术稳定反演过程
4. 通过诊断工具调试和优化算法性能

该方案为量子传感领域的信号处理提供了先进的数值工具，有望在精密测量、量子计量等领域发挥重要作用。