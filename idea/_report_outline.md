# 基于超导量子比特的瞬态磁场量子传感仿真平台 — 报告大纲

---

## 一、引言

### 1.1 研究背景与动机
- 量子传感的基本概念与优势（灵敏度、空间分辨率）
- 超导量子比特（Transmon）作为传感器的独特优势：可调谐频率、与磁通的强耦合、cQED 架构成熟
- 瞬态磁场测量的应用场景（材料科学、生物医学、基础物理）
- 现有方案的局限性：时间分辨率不足、仅支持线性响应、无法区分磁场极性

### 1.2 研究目标
- 建立完整的基于 Transmon 量子比特的瞬态磁场传感数值仿真平台
- 实现从量子动力学模拟、测量协议、到信号反演的全链路仿真
- 探索非线性响应、双极性磁场传感、微分传感等新方案
- 优化数值反演算法性能（伴随态方法、Levenberg-Marquardt）

### 1.3 报告结构概览

---

## 二、理论基础

### 2.1 超导量子电路基础
- LC 谐振器的量子化：从经典 LC 振荡器到量子谐振子
- 共面波导（CPW）谐振器：准 TEM 模式、品质因子 $Q_{int}$/$Q_{ext}$、零点电压涨落
- Josephson 结物理：
  - 两个基本关系：$I = I_c \sin\varphi$，$\dot{\varphi} = \frac{2\pi}{\Phi_0} V$
  - Josephson 电感 $L_J = L_{J0}/\cos\varphi$
  - 余弦势能 $E = -E_J \cos\varphi$

### 2.2 Transmon 量子比特
- Cooper pair box 哈密顿量：$H_T = 4E_C(\hat{n} - n_g)^2 - E_J \cos\hat{\varphi}$
- Transmon 近似（$E_J/E_C \gg 1$）：电荷色散指数压制 $\sim e^{-\sqrt{8E_J/E_C}}$
- 等离子体频率 $\omega_p = \sqrt{8E_JE_C}/\hbar$
- 量子比特频率 $\omega_{01} = \sqrt{8E_JE_C} - E_C$
- 非谐性 $\alpha = -E_C$
- Fock 基下的 Kerr 哈密顿量：$H \approx \hbar\omega_q \hat{b}^\dagger\hat{b} - \frac{E_C}{2}\hat{b}^\dagger\hat{b}^\dagger\hat{b}\hat{b}$
- 多种基矢表示（电荷基、相位基、Fock 基、能量本征基、Bloch 基）及其变换关系

### 2.3 磁通可调 Transmon（SQUID 结构）
- SQUID 等效 Josephson 能：$E_J(\Phi_x) = E_{J\Sigma}|\cos(\pi\Phi_x/\Phi_0)|$
- 频率-磁通关系：$\omega_{01}(\Phi) = \sqrt{8E_J(\Phi)E_C} - E_C$
- Sweet spot（$\Phi = 0, 0.5\Phi_0$）处 $d\omega/d\Phi = 0$ 的物理意义
- 最优工作点：$\Phi_{opt} = \arctan(\sqrt{2})$，最大 $|d\omega/d\Phi|$

### 2.4 光-物质相互作用与色散读出
- 交换相互作用与 Jaynes-Cummings 模型
- 色散区间：$H_{disp} = \hbar\omega_r' a^\dagger a + \frac{\hbar\omega_q'}{2}\sigma_z + \hbar\chi a^\dagger a \sigma_z$
- 色散位移 $\chi$ 与读出原理
- 测量链路：衰减器 → 环形器 → JPA/TWPA → HEMT → IQ 混频 → ADC/FPGA

### 2.5 退相干与噪声
- 能量弛豫（$T_1$）：坍缩算符 $\sqrt{\gamma_1}\hat{a}$，$\gamma_1 = 1/T_1$
- 纯退相（$T_\varphi$）：坍缩算符 $\sqrt{\gamma_\varphi}\hat{n}$
- 关系：$1/T_2 = 1/(2T_1) + 1/T_\varphi$
- Lindblad 主方程：$\dot{\rho} = -i[H,\rho] + \sum_k \gamma_k \mathcal{D}[L_k]\rho$
- 噪声源分类：电磁环境（Purcell 衰变）、材料/界面损耗（TLS）、准粒子、磁通 $1/f$ 噪声、电荷噪声、光子噪声

---

## 三、传感协议设计

### 3.1 基本量子操控协议
- **Rabi 振荡（协议 0）**：连续驱动下 $|1\rangle$ 布居数随时间的振荡，用于标定驱动强度
- **Ramsey 干涉（协议 1）**：$\pi/2 - \tau - \pi/2$ 脉冲序列，通过自由演化相位积累测量频率失谐
  - 理论：$P_e(\tau) = \frac{1}{2}[1 - \cos(\Delta\omega \cdot \tau)]$
- **自旋回波 / CPMG（协议 2, 3）**：通过重聚脉冲消除低频噪声（框架已建立，待实现）

### 3.2 瞬态磁场传感协议（协议 4）— 核心方案
- **基本原理**：滑动测量（sliding measurement）实现量子卷积
  $$p_e(\tau) = \int k(t) B(\tau - t) \, dt$$
- **传感核函数 $k(t)$**：
  - 物理含义：量子态对磁场脉冲扰动的线性响应函数
  - 提取方法：逐点施加窄高斯磁通扰动，记录 $\Delta p_e / \text{area}$
  - 理想 Ramsey 核：$k(t) = \sin\left(\frac{\pi}{2T}(T-t)\right)$
- **测量流程**：
  1. 制备 $|0\rangle$ 态
  2. 施加 Ramsey 控制脉冲序列
  3. 在不同时间延迟 $\tau$ 下施加磁场信号
  4. 测量 $p_e(\tau)$ 得到滑动扫描曲线
  5. 通过反卷积恢复 $B(t)$

### 3.3 双极性磁场传感方案（新方案）
- **问题**：Sweet spot 处 $\omega(\Phi)$ 为偶函数，$\arccos$ 反演无法区分 $\pm B$
- **解决方案**：磁通偏置工作点
  - 偏置展开：$\omega(\Phi_{bias} + \delta\Phi) \approx \omega_0 + \frac{d\omega}{d\Phi}\bigg|_{\Phi_{bias}} \cdot \delta\Phi + \cdots$
  - 一阶项携带磁场符号信息
- **品质因数优化**：
  $$\text{FoM}(\Phi_{bias}) = \left|\frac{d\omega}{d\Phi}\right|_{\Phi_{bias}} \cdot T_2^*(\Phi_{bias})$$
  - 灵敏度与相干时间的权衡
  - $1/T_2^* = 1/T_2 + \pi|d\omega/d\Phi| \cdot S_\Phi^{1/2}$
- **反演修正**：
  - 解析法：$\Phi_{total} = \frac{1}{\pi}\arccos(\cos\_val)$，$B = \Phi_{total} - \Phi_{bias}$
  - 查找表法：预计算 $\omega \to B$ 插值表

### 3.4 微分传感方案（新方案）
- **动机**：传统积分型传感 $\phi(t) = \int_0^t \omega(\tau)d\tau$ 时间分辨率受限
- **微分原理**：$d\phi/dt = \omega(t)$，直接测量相位变化率获取瞬时频率
- **四种实现途径**：
  1. 短时微分 Ramsey：灵敏度 $\delta\omega_{min} \approx 1/(\Delta t \cdot \sqrt{N})$
  2. 连续弱测量 + Kalman 滤波
  3. 微分核函数设计：频域优化 $\min \int |\mathcal{F}\{k_d\}(\omega) - i\omega|^2 d\omega$
  4. 双量子比特梯度计
- **噪声特性**：$S_{out}(\omega) = |\omega|^2 S_{in}(\omega)$，低频噪声被抑制，高频噪声被放大
- **正则化**：Tikhonov 微分 $\min \|Dx - y\|^2 + \alpha\|Lx\|^2$

---

## 四、信号反演理论与算法

### 4.1 线性反演：Wiener 反卷积
- 频域 Wiener 滤波器：
  $$\hat{B}(\omega) = \frac{K^*(\omega)}{|K(\omega)|^2 + \lambda^2} \hat{Y}(\omega)$$
- 正则化参数 $\lambda$ 的作用：噪声放大（过小）vs 过平滑（过大）
- 计算复杂度：$O(N\log N)$（FFT）

### 4.2 非线性反演：Hammerstein-Wiener 模型
- 两阶段反演：
  1. Wiener 反卷积得到 $\omega(t)$
  2. 反转 Transmon 频率-磁通关系：$B = \frac{1}{\pi}\arccos\left(\frac{(\omega + \omega_0 + E_C)^2}{8E_CE_J}\right)$
- 适用条件与局限性

### 4.3 全密度矩阵数值反演：Levenberg-Marquardt 算法
- **信号参数化**：$B(t) = \sum_{k=1}^M b_k \phi_k(t)$
  - 基函数选择：B-spline（局部支撑，推荐）、Fourier、Legendre
- **正向模型**：Lindblad 主方程
  $$\frac{d\rho}{dt} = -i[H(t;\mathbf{b}), \rho] + \mathcal{L}[\rho]$$
- **优化目标**：
  $$\min_{\mathbf{b}} \|\mathbf{p}_{meas} - \mathbf{p}_{sim}(\mathbf{b})\|^2 + \lambda \mathbf{b}^T D^T D \mathbf{b}$$
- **Levenberg-Marquardt 更新方程**：
  $$(J^TJ + \mu I + \lambda D^TD)\delta\mathbf{b} = J^T \mathbf{r}$$
- **自适应阻尼策略**：残差减小则接受并减半 $\mu$，否则拒绝并加倍 $\mu$

### 4.4 伴随态方法计算 Jacobian（核心理论推导）
- **正向传播**：$\rho_i(t)$ 满足 Lindblad 方程
- **伴随变量** $\lambda(t)$ 的定义与方程（从 $t_M$ 反向传播至 $t_m$）：
  $$\frac{d\lambda}{dt} = -i[H, \lambda] + \sum_k \gamma_k\left(L_k^\dagger \lambda L_k - \frac{1}{2}\{L_k^\dagger L_k, \lambda\}\right)$$
  $$\lambda(t_M) = |e\rangle\langle e|$$
- **Jacobian 表达式推导**：
  $$J_{ik} = \text{Re}\left\{-i\int_{t_m}^{t_M} \left(\frac{\partial\Delta\omega}{\partial B}\right)\phi_k(t)\,\text{Tr}[\lambda(t)[G, \rho(t)]]\,dt\right\}$$
- **关键优化**：$\text{Tr}[\lambda \cdot [G, \rho]]$ 与基函数索引 $k$ 无关，预计算一次后投影到各基函数
- **变量替换**：$s = t_M - t$，将终值问题转化为初值问题（数值稳定）
- **正则化矩阵构造**：
  - Fourier 基：$R = \text{diag}(n^2)$（频率惩罚）
  - 其他基：$R = D_2^T D_2$（二阶差分，惩罚曲率）

---

## 五、非线性响应扩展理论

### 5.1 线性假设的局限
- 现有核提取和 Wiener 反卷积均假设线性响应
- Transmon 频率响应本质非线性：$\omega(\Phi) = \sqrt{8E_J(\Phi)E_C} - E_C$

### 5.2 Volterra 级数表示
$$P_e(t) = \sum_{n=1}^N \int \cdots \int h_n(\tau_1, \ldots, \tau_n) \Phi(t-\tau_1) \cdots \Phi(t-\tau_n) d\tau_1 \cdots d\tau_n$$
- $h_1$：线性核（一阶响应）
- $h_2$：二阶核（非线性校正）
- 工作点线性化展开

### 5.3 非线性系统辨识
- Toeplitz 矩阵最小二乘法提取一阶核
- 残差法提取二阶核
- 高斯白噪声激励下的互相关方法

### 5.4 非线性反演
- 迭代 L-BFGS-B 最小化
- Tikhonov 非线性反演：$(H^TH + \alpha D^TD)\mathbf{x} = H^T\mathbf{y}$

---

## 六、仿真平台实现

### 6.1 软件架构设计
- 模块化结构：
  ```
  src/
  ├── qubit.py      — 量子比特模型（TransmonQubit, Cavity, Coupled_System）
  ├── signal.py     — 信号生成（7 种波形类型 + 基函数展开）
  ├── pulse.py      — 微波控制脉冲（Pulse, CompositePulse, DRAG 脉冲）
  ├── protocal.py   — 传感协议（Rabi, Ramsey, 瞬态场传感）
  └── analysis.py   — 信号分析与反演（Wiener, Hammerstein-Wiener, LM + 伴随态）
  ```
- 数据流：$\text{Signal} \to \text{Qubit} \to \text{Pulse} \to \text{Protocol} \to \text{Analysis}$

### 6.2 量子比特模块（qubit.py）
- Transmon 量子比特：哈密顿量构造（实验室系 / 旋转系）、坍缩算符、DRAG 门模拟
- 磁通响应：`qubit_in_mag()` 方法——预计算时变频率系数数组，返回 QuTiP list 格式哈密顿量
- 灵敏度计算：`frequency_sensitivity()` 中心差分法
- 1/f 噪声生成：逆 FFT 方法
- 多体系统：`Cavity`（多模谐振腔）、`Coupled_System`（双比特耦合系统、iSWAP/CZ 门模拟）

### 6.3 信号与脉冲模块（signal.py, pulse.py）
- 信号类型：零信号、常数、正弦、高斯、非对称脉冲、Slepian（待实现）、基函数展开
- 脉冲哈密顿量：实验室系、旋转系（RWA / 含反旋项）
- 复合脉冲序列：`CompositePulse` 拼接、全局时间轴映射
- Ramsey 脉冲工厂函数：`create_ramsey_pulse()`

### 6.4 协议与分析模块（protocal.py, analysis.py）
- 滑动测量实现：遍历所有时间延迟 $\tau$，构造时变哈密顿量并求解 `mesolve`
- 核函数提取：刺激-响应法
- 三级反演管线：Wiener → Hammerstein-Wiener → 数值 LM 反演
- 伴随态 Jacobian 计算的工程实现

### 6.5 交互演示系统
- Gradio Web 演示（`web_demo.py`）：参数可调的交互式传感仿真
- 命令行演示集（`appendix/`）：滑动测量、反卷积、完整协议流程、集成演示
- Jupyter Notebook（`Simulation.ipynb`）：开发与测试环境

---

## 七、性能优化

### 7.1 性能瓶颈分析
- 核心瓶颈：`mesolve` 的 Python 回调函数 `H_total` 在每个 ODE 步中创建 ~3 个 `Qobj` 对象
  - 420 扫描点 × 200 ODE 步 × 6 RK45 阶段 × 3 Qobj ≈ **150 万次 Qobj 创建**
  - GIL 锁竞争：多线程下比串行更慢
- LM 迭代中 `forward_simulation` 的重复调用（每迭代 ~24s，50 迭代 ~20 分钟）
- `store_states=True` 默认存储全部中间态（内存膨胀 ~420 倍）
- Jacobian 计算中 `g_values` 索引 bug 导致收敛缓慢

### 7.2 关键优化方案
1. **消除 Python 回调**：将哈密顿量转换为 QuTiP list/array 格式 `[H_op, coeff_array]`
   - 预计算系数数组，`mesolve` 在 C 层插值——GIL 全程释放
   - 预期加速：5~10 倍
2. **消除冗余 `forward_simulation`**：将正向模拟结果传递给 `compute_jacobian`
   - 预期节省：每迭代 ~25%
3. **修复 Jacobian bug**：`g_values[i]` → `g_values[s]`
   - 预期效果：收敛迭代数减半（50 → ~20 次）
4. **按需存储中间态**：默认 `store_states=False`
   - 内存节省：~500 倍
5. **伴随方程预计算**：构建 3D numpy 数组 `H_evolve_mats[j]`，消除 `adjoint` 中的 Python 回调
6. **向量化 trace 计算**：numpy broadcasting 替代双重循环

### 7.3 并行化策略
- GIL 分析：消除 Python 回调后，`mesolve` 全程在 C 层执行，GIL 释放
- 推荐方案：`threading` + `n_jobs=6`（仅使用 P-core，避免 E-core 拖慢）
- 备选方案：`loky` 多进程（模块级 worker 函数，numpy 序列化）
- **预期总加速：20~50 倍**（从 ~20 分钟降至 ~1.5 分钟）

### 7.4 多进程改造
- 三条铁律：模块级函数、可 pickle 参数、worker 内重建 QuTiP 对象
- Windows `spawn` 模式限制与解决方案
- 数据流改造：`QuTiP → .full() → numpy → pickle → subprocess → Qobj() → mesolve`

---

## 八、仿真结果与验证

### 8.1 基础协议验证
- Rabi 振荡仿真结果与理论对比
- Ramsey 干涉条纹：不同磁通偏置下的频率偏移验证

### 8.2 核函数提取与验证
- 数值核 vs 解析理想核 $k(t) = \sin(\frac{\pi}{2T}(T-t))$ 的对比
- 脉冲畸变对核函数形状的影响

### 8.3 瞬态磁场重建
- Wiener 反卷积结果：不同正则化参数 $\lambda$ 下的重建质量（RMSE、Pearson 相关系数）
- Hammerstein-Wiener 非线性校正效果
- 全密度矩阵数值反演结果

### 8.4 新方案初步验证
- 双极性磁场传感：偏置工作点下的符号区分能力
- 1/f 噪声对传感性能的影响

---

## 九、讨论与展望

### 9.1 当前工作总结
- 建立了完整的 Transmon 量子传感仿真平台
- 实现了从线性到非线性的多级反演算法
- 提出了双极性传感、微分传感等新方案
- 完成了详细的性能瓶颈分析与优化路线

### 9.2 待完善工作
- Spin Echo / CPMG 协议实现
- 完整噪声模型集成（退相干、磁通噪声、脉冲误差、读出噪声）
- 二阶 Volterra 核的实际提取与验证
- GPU 加速与 JAX 自动微分集成

### 9.3 未来方向
- 非马尔科夫噪声模型
- 基于机器学习的脉冲序列优化
- 量子增强微分传感（纠缠辅助）
- 贝叶斯反演与 MCMC 不确定性量化
- PINN（物理信息神经网络）加速正向模拟
- 实时递归最小二乘在线传感
- 硬件接口与数字孪生

---

## 附录

### A. Transmon 量子比特多基矢表示详细推导
（对应 `appendix/appendix.py`：Fock 基、电荷基、相位基、能量本征基、Bloch 基的完整构造与变换矩阵）

### B. 伴随态方法的完整数学推导
（对应 `idea/numerical_inversion_scheme.md`：从灵敏度方程到伴随方程的详细推导，变量替换，Jacobian 投影公式）

### C. 代码结构与 API 文档
（各模块类和方法的完整说明）

### D. 演示系统使用指南
（对应 `RUN_DEMO.md`：Web 演示安装、运行、参数说明）

---

## 参考文献

- QuTiP 文档：时间依赖哈密顿量、`parallel_map`
- cQED 笔记（项目内 `note/cQED note/`）
- SQC 笔记（项目内 `note/SQC note/`）
- Intel Thread Director 白皮书（混合架构 CPU 并行优化）
- Tikhonov 正则化、Levenberg-Marquardt 算法相关文献
- Volterra 级数与非线性系统辨识相关文献
