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
   Jacobian矩阵表示第k个参数对第i个测量点的敏感度

   c. 更新方程：$(\mathbf{J}^{(n)T}\mathbf{J}^{(n)} + \mu^{(n)}\mathbf{I} + \lambda\mathbf{D}^T\mathbf{D}) \delta\mathbf{b} = \mathbf{J}^{(n)T} \mathbf{r}^{(n)}$

   d. 参数更新：$\mathbf{b}^{(n+1)} = \mathbf{b}^{(n)} + \delta\mathbf{b}$
   
   e. 阻尼调整：若 $\|\mathbf{r}^{(n+1)}\| < \|\mathbf{r}^{(n)}\|$，则 $\mu^{(n+1)} = \mu^{(n)}/2$；否则 $\mu^{(n+1)} = 2\mu^{(n)}$ 并拒绝更新

3. **收敛判断**：当 $\|\delta\mathbf{b}\|/\|\mathbf{b}\| < \epsilon_{\text{tol}}$ 或 $\|\mathbf{r}\| < \epsilon_{\text{data}}$ 时停止


### Jacobian计算

采用伴随法（adjoint method）高效计算Jacobian矩阵$$J_{ik} = \frac{\partial p_{\text{sim}}(t_i)}{\partial b_k}$$

1. 前向传播：密度矩阵 $\rho(t)$满足$$\frac{d\rho_i}{dt} = \mathcal{L^i}[\rho]=-i[H_i(t;\mathbf{b}), \rho] + \sum_k \gamma_k(L_k\rho L_k^\dagger - \frac{1}{2} \{L_k^\dagger L_k, \rho\})$$
设脉冲的起点为$\delta_i = t_i - T_{pulse}/2$，演化时间范围为$t_m = min(0, \delta_i), t_M = max(T_{field}, t_i + T_{pulse}/2)$，初始条件为 $\rho_i(t_m) = \rho_0$，观测量为$p(t_i) = \text{Tr}[\rho_i(t_M) \cdot |1\rangle\langle 1|]$。
2. 反向传播：为了得到雅可比矩阵，需要计算方程$$\frac{\partial p(t_i, t)}{\partial b_k} = \text{Tr} [\frac{\partial \rho(t_i, t)}{\partial b_k} \ket{1}\bra{1}] = \text{Tr}[\sigma_{k}^i(t)\ket{1}\bra{1}]$$

$\sigma_k^i(t)$满足方程$$ \frac{d\sigma_k^i}{dt} = \frac{\partial }{\partial b_k}\frac{\partial }{\partial t}\rho = \frac{\partial \mathcal{L}}{\partial b_k}[\rho] + \mathcal{L}[\sigma_k] = \mathcal{L}[\sigma_k] + \mathcal{S_k}[\rho] $$
其中，$\mathcal{S_k}[\rho]$为源项，当耗散超算符与$b_k$无关时，$\mathcal{S_k}[\rho] = -i[\frac{\partial H(t)}{\partial b_k},\rho]$。

$\sigma_k$满足$\sigma_{k}^i(t_m) = 0$，利用这个性质，
引入函数$\lambda(t)$，有$$\begin{aligned}\text{Tr}[\ket{1}\bra{1}\sigma_{k}^i(t_M)] &= \text{Tr}[\lambda(t_M)\sigma_k(t_M)] - \text{Tr}[\lambda(t_m)\sigma_k(t_m)]\\ &= \int_{t_m}^{t_M} \frac{d}{dt}\bigg\{\text{Tr}\left[\lambda(t) \sigma_k(t)\right]\bigg\} dt \\ &= \int_{t_m}^{t_M} \text{Tr}\bigg[\dot\lambda\sigma_k + \lambda \dot\sigma_k \bigg]dt \\ &= \int_{t_m}^{t_M} \text{Tr}(\dot\lambda\sigma_{k} ) + \text{Tr}(\lambda\mathcal{L[\sigma_k]}) -i \ \text{Tr}(\lambda[\frac{\partial H}{\partial b_k},\rho])dt \\ &= \int_{t_m}^{t_M} \text{Tr}\left[(\dot\lambda(t) + \mathcal{L}^\dagger[\lambda])\sigma_k\right] -i\  \text{Tr}(\lambda[\frac{\partial H}{\partial b_k},\rho])dt \\ &=-i\int_{t_m}^{t_M}  \text{Tr}(\lambda[\frac{\partial H}{\partial b_k},\rho])dt \end{aligned}$$
其中，$\lambda(t) $满足：$$\frac{d\lambda}{dt} = -\mathcal{L}^\dagger[\lambda], \quad \lambda(t_M) = |e\rangle\langle e|$$
上述推导利用了超算符的性质：$$\text{Tr}[\lambda \mathcal{L}[\sigma_k]] = \text{Tr}[\mathcal{L}^\dagger[\lambda]\sigma_k]$$

于是得到了完整的伴随方程：$$\dot\lambda = -i[H, \lambda] + \sum_k \gamma_k(L_k^\dagger \lambda L_k - \frac{1}{2} \{L_k^\dagger L_k, \lambda\})$$

作变量替换$s = t_M - t, \mu(s) = \lambda(t) =  \lambda(t_M - s) $，则$$\frac{d\mu}{ds} = i[H(t_M - s), \mu] + \sum_k \gamma_k(L_k^\dagger \mu L_k - \frac{1}{2} \{L_k^\dagger L_k, \mu\}), \quad \mu(0) = |e\rangle\langle e|$$
上式将终值条件问题转换为初值问题，数值求解更为稳定。
求得$\lambda(t)$后，雅可比矩阵为$$J_{ik} = \int_{t_m}^{t_M} \text{Tr}\left[\lambda^i(t) \mathcal{S}_k^i[\rho_i(t)]\right] dt$$

3. 梯度计算：$\frac{\partial p_{\text{sim}}(t_i)}{\partial b_k} = \text{Re} \bigg\{-i\int_{t_m}^{t_M} \text{Tr}\left[\lambda(t) [\frac{\partial H(t)}{\partial b_k}, \rho(t)]\right] dt\bigg\}$

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



    # ================================================================
    #  主方法: 伴随法Jacobian
    # ================================================================
    def compute_jacobian_adjoint(self, b, t_points):
        """
        使用伴随法计算Jacobian矩阵 J_{ik} = ∂p(t_i)/∂b_k

        数学基础
        ========
        正向方程:  dρ/dt = L[ρ] = -i[H,ρ] + D[ρ]
        伴随方程:  dλ/dt = -L†[λ] = -i[H,λ] - D†[λ]

        其中 D†[λ] = Σ_k γ_k(L_k† λ L_k - ½{L_k†L_k, λ})

        梯度公式:
          ∂p(t_i)/∂b_k = Re{ -i ∫₀^{t_i} Tr[λ(t) [∂H/∂b_k, ρ(t)]] dt }

        其中: ∂H/∂b_k = (∂Δω/∂B · φ_k(t)) · Ĝ
              Ĝ = σ_z/2 (二能级) 或 diag(0,1,...) (多能级)

        算法流程
        ========
        1. 前向传播: ρ(0) → ρ(T), 保存 ρ(t), B(t) 在所有时间步
        2. 对每个观测点 t_i:
           a. 后向传播: λ(t_i)=Ô → λ(0)
           b. 预计算: f(t) = Tr[λ(t)·[Ĝ, ρ(t)]] (不依赖k)
           c. 梯度: J[i,k] = Re{-i ∫ (∂Δω/∂B · φ_k(t)) · f(t) dt}

        计算量: 2×N_t 次ODE积分 + N_t×M 次标量积分 (不含额外ODE)

        参数
        ====
          b : array (M,)
              基函数系数向量
          t_points : array (N_t,)
              观测时间点

        返回
        ====
          J : array (N_t, M)
              Jacobian矩阵
        """
        N_t = len(t_points)
        M = self.M
        J = np.zeros((N_t, M))

        # ==============================================================
        # 第一步: 构建细密时间网格 + 前向传播
        # ==============================================================

        # 时间范围 (在观测窗口前后留余量)
        t_margin = getattr(self, 'tau_pulse', 5e-9) * 1.5
        t_start = t_points[0] - t_margin
        t_end = t_points[-1] + t_margin
        dt_fine = getattr(self, 'dt_sim', 0.01e-9)  # 默认10ps
        t_fine = np.arange(t_start, t_end + dt_fine / 2, dt_fine)
        N_fine = len(t_fine)

        # 计算 B(t) 和基函数值
        B_fine, phi_vals = self._evaluate_field_and_basis(b, t_fine)
        # B_fine:  shape (N_fine,)  — 各时间步的磁场值
        # phi_vals: shape (M, N_fine) — 各基函数在各时间步的值

        # 前向传播 (使用QuTiP的mesolve)
        H_total = self._build_hamiltonian_with_field(B_fine, t_fine)
        c_ops = self.qubit.get_collapse_operators()
        rho_0 = self.qubit.state

        opts = dict(store_states=True, nsteps=50000, max_step=dt_fine * 10)
        result = mesolve(H_total, rho_0, t_fine, c_ops, options=opts)

        # 转换为numpy数组轨迹
        dim = rho_0.shape[0]
        rho_traj = np.zeros((N_fine, dim, dim), dtype=complex)
        for n in range(N_fine):
            rho_traj[n] = result.states[n].full()

        # 测量算符 (numpy矩阵)
        O_mat = self.qubit.excited_state_projector().full()

        # ==============================================================
        # 第二步: 预计算不依赖于k和i的量
        # ==============================================================

        # ∂Δω/∂B 在每个时间步的值
        dw_dB_array = np.array([self._compute_dw_dB(B_fine[n])
                                for n in range(N_fine)])

        # Ĝ算符: ∂H/∂(Δω) 中的算符部分
        G_mat = self._get_dH_dbk_operator(dim)

        # collapse算符的numpy矩阵形式 (后向传播用)
        c_mats = [c.full() for c in c_ops]
        c_dag_mats = [c.dag().full() for c in c_ops]     # L†
        cdc_mats = [cd @ c for c, cd in zip(c_mats, c_dag_mats)]  # L†L

        # B(t) 插值器 (后向传播中在任意时间获取B值)
        B_interp = interp1d(t_fine, B_fine, kind='linear',
                            fill_value=(B_fine[0], B_fine[-1]),
                            bounds_error=False)

        # ==============================================================
        # 第三步: 对每个观测点——后向传播 + 梯度积分
        # ==============================================================

        for i in range(N_t):
            t_obs = t_points[i]

            # 找到 t_obs 在 t_fine 中的最近索引
            idx_obs = np.argmin(np.abs(t_fine - t_obs))
            # 后向传播覆盖 t_fine[0] 到 t_fine[idx_obs]
            N_bwd = idx_obs + 1

            if N_bwd < 2:
                # 观测点过早,跳过
                continue

            # ──────────────────────────────────────────────
            #  3a. 后向传播伴随态
            # ──────────────────────────────────────────────
            #
            # 变量替换: s = t_obs - t,  μ(s) = λ(t_obs - s)
            #
            # dμ/ds = +i[H(t_obs-s), μ] + D†(t_obs-s)[μ]
            #
            # 其中 D†[μ] = Σ_k γ_k (L_k† μ L_k - ½{L_k†L_k, μ})
            #
            # 初始条件: μ(0) = λ(t_obs) = Ô
            #
            # s从0积分到 (t_obs - t_fine[0])

            t_obs_actual = t_fine[idx_obs]  # 对齐到网格
            s_max = t_obs_actual - t_fine[0]

            # s网格: 对应原始时间从 t_obs 到 t_fine[0]
            # s_n = t_obs - t_fine[idx_obs - n], n = 0,1,...,idx_obs
            s_eval = t_obs_actual - t_fine[:N_bwd][::-1]
            # s_eval[0] = 0 (对应t=t_obs)
            # s_eval[-1] = t_obs - t_fine[0] (对应t=t_fine[0])

            def adjoint_rhs(s, mu_flat, *args):
                """
                伴随方程右端项 (s变量):
                  dμ/ds = +i[H(t_obs-s), μ]
                        + Σ_k γ_k (L_k† μ L_k - ½{L_k†L_k, μ})
                """
                mu = mu_flat.reshape(dim, dim)
                t_curr = t_obs_actual - s

                # 获取当前时刻的Hamiltonian
                B_curr = float(B_interp(t_curr))
                H_curr = self._get_hamiltonian_matrix(t_curr, B_curr)

                # ── 酉部分: +i[H, μ] ──
                dmu = 1j * (H_curr @ mu - mu @ H_curr)

                # ── Lindblad伴随部分: Σ γ_k(L†μL - ½{L†L, μ}) ──
                for c_mat, cd_mat, cdc_mat in zip(c_mats, c_dag_mats, cdc_mats):
                    dmu += (cd_mat @ mu @ c_mat
                            - 0.5 * (cdc_mat @ mu + mu @ cdc_mat))

                return dmu.flatten()

            # 初始条件
            mu_0 = O_mat.copy().flatten()

            # 使用scipy积分
            sol = solve_ivp(
                adjoint_rhs,
                t_span=(0, s_max),
                y0=mu_0,
                method='RK45',
                t_eval=s_eval,
                rtol=1e-9,
                atol=1e-11,
                max_step=dt_fine * 5
            )

            # 检查积分是否成功
            if not sol.success:
                print(f"  警告: 观测点{i} (t={t_obs*1e9:.2f}ns) "
                      f"后向积分失败: {sol.message}")
                continue

            # ──────────────────────────────────────────────
            #  3b. 将μ(s)转换回λ(t), 对齐到t_fine网格
            # ──────────────────────────────────────────────
            #
            # sol.t 中的 s_n 对应 t = t_obs - s_n
            # s_n 从小到大 → t 从大到小
            # 需要反转以对齐 t_fine[0:N_bwd]

            N_sol = len(sol.t)
            lambda_traj = np.zeros((N_bwd, dim, dim), dtype=complex)

            for n in range(N_sol):
                s_n = sol.t[n]
                # s_n 对应原始时间 t_obs - s_n
                # 在t_fine中的索引 (从后往前)
                j = N_bwd - 1 - n
                if 0 <= j < N_bwd:
                    lambda_traj[j] = sol.y[:, n].reshape(dim, dim)

            # 对缺失点做插值 (如果sol.t和s_eval不完全对齐)
            # 通常如果t_eval给得正确, 应该没有缺失

            # ──────────────────────────────────────────────
            #  3c. 计算梯度积分
            # ──────────────────────────────────────────────
            #
            # ∂p(t_i)/∂b_k = Re{ -i ∫₀^{t_obs} Tr[λ(t)·[∂H/∂b_k, ρ(t)]] dt }
            #
            # ∂H/∂b_k = (∂Δω/∂B)(t) · φ_k(t) · Ĝ
            #
            # 所以:
            # ∂p/∂b_k = Re{ -i ∫ (∂Δω/∂B)·φ_k(t) · Tr[λ(t)·[Ĝ, ρ(t)]] dt }
            #
            # 关键优化: Tr[λ·[Ĝ,ρ]] 不依赖于k, 可预计算一次

            # 预计算 f(t) = Tr[λ(t) · [Ĝ, ρ(t)]]
            trace_values = np.zeros(N_bwd, dtype=complex)

            for n in range(N_bwd):
                rho_n = rho_traj[n]       # ρ(t_n), 前向轨迹
                lam_n = lambda_traj[n]    # λ(t_n), 伴随轨迹

                # [Ĝ, ρ] = Ĝρ - ρĜ
                comm_G_rho = G_mat @ rho_n - rho_n @ G_mat

                # Tr[λ · [Ĝ, ρ]]
                trace_values[n] = np.trace(lam_n @ comm_G_rho)

            # 被积函数 (不含φ_k): g(t) = -i · (∂Δω/∂B)(t) · f(t)
            g_values = -1j * dw_dB_array[:N_bwd] * trace_values

            # 对每个参数k, 乘以φ_k(t)并积分
            t_bwd = t_fine[:N_bwd]  # 积分的时间网格

            for k in range(M):
                # 被积函数 = g(t) · φ_k(t)
                integrand_k = g_values * phi_vals[k, :N_bwd]

                # 梯形法则积分
                integral = np.trapz(integrand_k, t_bwd)

                # 取实部 (虚部应为数值零)
                J[i, k] = np.real(integral)

            # 进度报告
            if (i + 1) % max(1, N_t // 10) == 0 or i == N_t - 1:
                p_i = np.real(np.trace(O_mat @ rho_traj[idx_obs]))
                print(f"  伴随法进度: {i+1}/{N_t}, "
                      f"p({t_obs*1e9:.1f}ns) = {p_i:.4f}")

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




