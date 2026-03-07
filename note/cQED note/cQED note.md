# 第一部分：基础与 Transmon 模型 (对应任务 1 & 2)

  

## II. Superconducting quantum circuits

  

### A. The quantum LC resonator

  

#### 1. 核心思想：从经典电路到量子电路

  

*   **集总元件 (Lumped Element) 近似**:

    *   当电路尺寸 $d \ll \lambda$ (波长) 时，电路可以用分立的 L (电感) 和 C (电容) 来描述。

    *   **LC 振荡器 = 谐振子**:

        *   能量公式 $H = \frac{Q^2}{2C} + \frac{\Phi^2}{2L}$ 完美对应机械谐振子 $H = \frac{p^2}{2m} + \frac{1}{2}kx^2$。

        *   **映射关系**:

            *   电荷 $Q \leftrightarrow$ 动量 $p$

            *   磁通 $\Phi \leftrightarrow$ 位置 $x$

            *   电容 $C \leftrightarrow$ 质量 $m$

            *   电感 $L^{-1} \leftrightarrow$ 弹簧系数 $k$

  

*   **正则量子化 (Canonical Quantization)**:

    *   将经典变量 $Q, \Phi$ 变为算符 $\hat{Q}, \hat{\Phi}$。

    *   施加对易关系 $[\hat{\Phi}, \hat{Q}] = i\hbar$。

    *   引入产生/湮灭算符 ($\hat{a}^\dagger, \hat{a}$)，哈密顿量变为 $\hbar\omega_r(\hat{a}^\dagger\hat{a} + 1/2)$。

    *   **物理意义**: 电路中的电压和电流不再是连续值，而是存在最小的能量包（微波光子）。

  

*   **实现量子行为的两个条件**:

    1.  **低损耗 (High Q)**: $Q = \omega/\kappa \gg 1$。光子在消失前能振荡很多次。这需要使用超导体（如铝 Al，铌 Nb）来消除电阻。

    2.  **低温度 (Low T)**: $k_B T \ll \hbar \omega$。防止热涨落激发光子。

        *   典型频率: 5-10 GHz ($\sim 240-480$ mK)

        *   工作温度: 10-20 mK (稀释制冷机)

        *   结果: 电路自然处于基态 $|0\rangle$。

  

*   **零点涨落 (Zero Point Fluctuations, ZPF)**:

    *   即使在基态（真空态），电压和电流也不为零，而是存在涨落。

    *   特征阻抗 $Z_r = \sqrt{L/C}$ 决定了涨落的分配：

        *   $Z_r$ 高 $\rightarrow$ $\Phi_{zpf}$ 大 (磁通噪声大)

        *   $Z_r$ 低 $\rightarrow$ $Q_{zpf}$ 大 (电荷噪声大)

    *   **关键点**: 真空电压涨落 $\Delta V_0 \sim 1 \mu V$ 在微观尺度下产生的电场极大，这是实现强光-物质耦合（Strong Coupling）的基础。

  

#### 2. 谐振器类型 (Resonators)

  

##### A.2D 共面波导 (Coplanar Waveguide, CPW)



## 1. 目标与场景：为什么讨论 2D 谐振腔？

1) **电路 QED 的基本器件**：微波谐振器是量子谐振子的一种实现（另一类是集总参数 LC 振荡器）。在电路量子电动力学中，谐振器扮演"光腔"角色——要么存储微波光子作为量子存储器，要么作为读出量子比特状态的探测器。

2) **追求高品质因数**：实验室可实现很宽范围的品质因数：
   $$Q \sim 10^3 \text{ 到 } 10^8$$

3) **材料选择与损耗意识**：
   - 为获得更大超导能隙（以及相关性能），有时用 **Nb 替代 Al**（Nb 的能隙对应更高频率上限）
   - 不仅要关注金属/介质的内损耗，还要控制与外部电路耦合带来的外损耗——因为耦合既是驱动/读出所必需的，也是退相干的主要来源之一

---

## 2. 进入"量子工作区"的两条硬条件

### 2.1 热激发要小：$\hbar\omega \gg k_BT$
- 相邻本征态能级差 $\hbar\omega$ 必须显著大于热能 $k_BT$，否则热涨落会随机激发光子，破坏量子相干性。
- **换算记忆点**：**1 GHz 对应 $\hbar\omega/k_B \approx 50$ mK**（推导：$\hbar\omega/k_B = h\times 10^9 / k_B \approx 6.626\times10^{-25} / 1.38\times10^{-23} \approx 0.048$ K）。
- 稀释制冷机典型工作温度 **$T\sim 10$ mK** 下，GHz 量级（5-15 GHz）电路容易满足量子条件，但几十 MHz 的电路即使在此低温下热占据数仍可能大于1。

### 2.2 模必须"窄线宽"：$\kappa_m \ll \omega_m$
- 线宽（耗散率）$\kappa_m$ 要远小于模频率 $\omega_m$，否则光子还没完成一次振荡就耗散掉了，无法形成清晰的能级。
- 用品质因数表示：
  $$\kappa_m=\frac{\omega_m}{Q_m},\qquad \frac{\kappa_m}{2\pi}=\frac{f_m}{Q_m}$$
- 物理意义：$Q$ 越大，光子在谐振器内振荡的次数越多（$Q \sim$ 振荡次数），量子相干性保持得越好。

---

## 3. 2D 谐振腔的定义：2D vs 3D

- **2D（平面）谐振腔**：电磁场主要限制在平面结构附近（几何上"近二维"），通过光刻技术在芯片上制备。
- **3D 谐振腔**：电磁场限制在三维腔体体积中（如挖空的金属块），通常由整块超导材料加工而成。
- **关键区别**：3D 腔的模体积大，场主要分布在真空中，表面参与率低，因此内损耗极小（$Q$ 可达 $10^6-10^8$）；2D 腔模体积小，便于集成，但表面损耗更大。
- 无论哪种几何，边界条件使场离散成一组本征模；**每个模都可看作独立谐波振子**。

---

## 4. 核心器件：共面波导（CPW）谐振器长什么样？

### 4.1 几何结构（对应图 2a/2b）
CPW 谐振器 = **有限长度 $d$** 的 CPW 传输线段：
- 中心导体：宽度 $w$，厚度 $t$
- 两侧接地平面：与中心导体同厚度 $t$
- 中心导体与地之间的间隙：每侧间距 $s$
- 典型制作：沉积在**低损耗介质衬底**上，衬底介电常数 $\epsilon$，衬底厚度通常 **远大于 $w,s,t$**（避免衬底中的电磁场模式与 CPW 模式耦合）
![[Pasted image 20260218163822.png]]
### 4.2 场分布与传播模式：quasi‑TEM
- CPW 在平面上充当传输线，传播方式类似同轴线：把电磁场限制在中心导体与地之间的小体积。
- **quasi‑TEM 模式**：场一部分在衬底介质中，一部分在衬底上方真空/空气中；但**最强集中在 gap（缝隙）区域**（这是量子比特与谐振器耦合的关键位置——量子比特通常放在此处以最大化耦合）。
- **关键设计思想**：通过选择 $w,s,t$ 及介质相关几何，使场尽量集中在中心导体与地之间，从而**减少向其他方向辐射**，并最大化与量子比特的耦合强度。

---

## 5. 传输线参数：用 $l_0,c_0$ 把几何"压缩"为可计算模型

### 5.1 两个单位长度参数
- 单位长度电感 $l_0$（由几何和穿透深度决定）
- 单位长度对地电容 $c_0$（由几何和介电常数决定）

### 5.2 两个最常用派生量（必须会背）
$$Z_r=\sqrt{\frac{l_0}{c_0}},\qquad v_0=\frac{1}{\sqrt{l_0c_0}}$$
- 典型量级（图中给出）：
  - $Z_r\sim 50~\Omega$（与标准微波阻抗匹配，便于连接外部设备）
  - $v_0\sim 1.3\times10^8~\text{m/s}$（约真空光速的 1/3，由有效介电常数决定）

### 5.3 几何如何影响阻抗（工程要点）
- 改变 $w,s,t$ 会改变 $l_0,c_0$，进而改变 $Z_r$。
- **常见技巧**：让 $w/s$ 近似保持恒定，可在改变线宽的同时保持近似恒阻抗——便于端口做宽（易于 bonding 和连接）、腔体内部做窄（减小模体积以增强与量子比特的耦合）。

---

## 6. "谐振腔"从哪里来：边界条件 + 有限长度

把长度为 $d$ 的传输线段两端施加边界条件即可形成谐振器——相当于在传输线上形成驻波。

### 6.1 两类边界条件
- **开路端（open）**：端点电流为零（可用中心导体做 gap 实现，gap 形成高阻抗）
- **短路端（short/grounded）**：端点电压为零（直接接地）

### 6.2 $\lambda/2$ 与 $\lambda/4$ 谐振器（频率公式必须熟）
- **两端开路（open-open）** → $\lambda/2$ 谐振器：
  $$f_0=\frac{v_0}{2d},\qquad f_m=(m+1)f_0$$
  （电压驻波：两端为波腹，中间为波节）

- **一端开路一端短路（open-short）** → $\lambda/4$ 谐振器：
  $$f_0=\frac{v_0}{4d}$$
  （电压驻波：开路端为波腹，短路端为波节）

### 6.3 典型频段与极端例子（图中文字给出的信息）
- **典型设计**：$d=1$ cm、$v_0=1.3\times10^8$ m/s → $f_0\approx 6.5$ GHz（处于 circuit QED 常用频段）。
- **高频上限**：受超导能隙限制（超过能隙会产生准粒子损耗）；对 **Al**，图中给出约 **82 GHz** 的量级。
- **低频极端**：可通过超长蛇形线实现低频率，例如 $d=0.68$ m，$f_0=92$ MHz；但在 $10$ mK 下该频率的热占据数可能仍不小（图中给出基模热占据约 $\sim 1.8$ 的量级），因此不一定在真空态。
- **circuit QED 常用频率范围**：**5–15 GHz**（此频段兼顾了低温条件、微波电子学成熟度、以及量子比特频率设计的灵活性）。

---

## 7. 耦合与损耗：$Q_{\text{int}}, Q_{\text{ext}}, Q_L$ 一套逻辑

### 7.1 外耦合怎么做（图 2a & 文字）
- 输入/输出端口通常通过**电容耦合**到外部传输线（这样才能把信号送进/取出谐振器，又不至于严重破坏其量子性）。
- 耦合电容可由：
  - 简单 gap 形成的电容
  - 或叉指电容（interdigitated capacitors，可精确控制电容值）

### 7.2 两类损耗的物理起源
- **内损耗（internal）**：与不可控自由度耦合导致能量耗散
  - 介质与导体表面/界面损耗（主要是非晶态介质中的两能级系统 TLS）
  - 衬底介质损耗（如声子激发）
  - 非平衡准粒子（超导能隙未被完全打开的激发）
  - 磁通涡旋运动（若存在 trapped flux）
  - **TLS（two-level systems）** 尤为重要——它们存在于各种无序介质中（如氧化物层、表面吸附物），尤其在电场强的界面区域影响最大
  
- **外损耗（external）**：与输入/输出端口耦合导致能量泄露（读出/驱动需要，但会降低 Q——这是必要的代价）

### 7.3 总耗散率与加载 Q（必须会用）
$$\kappa_m=\kappa_{\text{ext},m}+\kappa_{\text{int},m}$$
$$Q_{L,m}=\left(Q_{\text{ext},m}^{-1}+Q_{\text{int},m}^{-1}\right)^{-1}$$
- 工艺进步后，图中文字给出：**$Q_{\text{int}}\sim 10^5$** 已较常见。
- **过耦合（overcoupled）**：$Q_{\text{ext}}<Q_{\text{int}}$（线宽主要由外耦合决定，利于"快测量"——信号能快速进出谐振器）
- **欠耦合（undercoupled）**：$Q_{\text{ext}}>Q_{\text{int}}$（线宽主要由内损耗决定；若内损耗也很小，可用于长时间存储微波光子——量子存储器）

---

## 8. 为什么 2D CPW 容易实现强耦合：零点涨落更大

### 8.1 电压零点涨落（图中文字的典型数值链）
给一组典型参数（示例）：$L\sim0.8$ nH、$C\sim0.4$ pF、$\omega_r/2\pi\sim 8$ GHz、$Z_r\sim50\Omega$，可得基态电压零点涨落量级：
$$\Delta V_0\simeq \sqrt{\frac{\hbar\omega_r}{2C}}\sim 1~\mu\text{V}$$
（推导：谐振器基态能量 $\hbar\omega_r/2 = C(\Delta V_0)^2/2$）

### 8.2 横向尺寸可做得很小 → 电场零点涨落更大
- 常用横向几何：$w\sim 10~\mu\text{m}$，$s\sim 5~\mu\text{m}$。
- 进一步缩小到亚微米时，会受超导薄膜**穿透深度**限制（图中给 **100–200 nm** 量级——若线宽小于穿透深度，电感会增加，阻抗变化）。
- 用 $\Delta V_0\sim 1~\mu\text{V}$、$s\sim 5~\mu\text{m}$ 估算电场零点涨落：
  $$\Delta E_0\approx \frac{\Delta V_0}{s}\sim 0.2~\text{V/m}$$
- **图中文字结论**：这比 3D 腔 QED 的典型零点电场强**至少大两个数量级**；配合超导人工原子的大电偶极矩（约瑟夫森结的相位涨落对应大有效电荷），形成 circuit QED 的强光-物质耦合优势（$g/2\pi$ 可达几十到几百 MHz，而 3D 腔 QED 通常只有几十 kHz）。

---

## 9. 量子化推导主线（从"电报模型"到"模态谐振子"）

这一部分的逻辑是：**离散电路 → 连续场 → 波动方程 → 正常模 → 每模一个谐振子 → 升格为量子算符**。这是从经典传输线到量子谐振子的标准路径。

### 9.1 为什么要走这条路？
虽然集总 LC 谐振子可以直接量子化，但传输线谐振器是分布参数系统，有多个模式。必须通过正常模分解，才能正确写出每个模式的量子化形式。

### 9.2 离散电报模型（图 3）与哈密顿量（式 6）
把传输线离散成 $N$ 段（每段 $\delta x$）：
- 串联电感 $L_0=\delta x\,l_0$（每段储存磁能）
- 并联对地电容 $C_0=\delta x\,c_0$（每段储存电能）
- 节点变量：磁通 $\Phi_n$（类比"位置"）、电荷 $Q_n$（类比"动量"）

经典哈密顿量（每个节点的电能 + 相邻节点间的磁能）：
$$H=\sum_{n=0}^{N-1}\left[\frac{Q_n^2}{2C_0}+\frac{(\Phi_{n+1}-\Phi_n)^2}{2L_0}\right]$$

### 9.3 连续极限（式 7）与变量定义
取 $\delta x\to 0$，定义连续场：
$$\Phi(x_n)=\Phi_n,\qquad Q(x_n)=\frac{Q_n}{\delta x}$$
得到连续哈密顿量：
$$H=\int_0^d dx\left[\frac{Q(x,t)^2}{2c_0}+\frac{(\partial_x\Phi(x,t))^2}{2l_0}\right]$$
并有共轭关系：
$$Q(x,t)=c_0\,\partial_t\Phi(x,t)$$
广义磁通与电压的关系：
$$\Phi(x,t)=\int_{-\infty}^{t}dt'\,V(x,t')$$
（磁通是电压的时间积分，类比于力学中的位移是速度的时间积分）

### 9.4 波动方程（式 8）
由哈密顿方程可得：
$$v_0^2\partial_x^2\Phi-\partial_t^2\Phi=0,\qquad v_0=\frac{1}{\sqrt{l_0c_0}}$$
这是一维波动方程，$v_0$ 是波速。

### 9.5 正常模展开（式 9–10）
解可写为正常模叠加（时间部分简谐振动 × 空间分布）：
$$\Phi(x,t)=\sum_{m=0}^{\infty}u_m(x)\Phi_m(t)$$
其中时间部分满足 $\ddot{\Phi}_m=-\omega_m^2\Phi_m$，空间模函数为：
$$u_m(x)=A_m\cos(k_m x+\varphi_m),\quad k_m=\omega_m/v_0$$
（这是波动方程的通解形式，具体由边界条件决定 $k_m$ 和 $\varphi_m$）

### 9.6 开路端边界条件（式 11）与归一化（式 12）
开路端电流为零：
$$I(x)= -\frac{1}{l_0}\partial_x\Phi(x,t),\quad I(0)=I(d)=0$$
这要求 $\partial_x\Phi(0)=\partial_x\Phi(d)=0$，即 $u_m'(0)=u_m'(d)=0$。

模函数正交归一化条件（能量正交）：
$$\int_0^d dx\,u_m(x)u_{m'}(x)=\delta_{mm'}$$
图中文字给出一种结果：可取 $A_m=\sqrt{2}$（并指出这意味着模幅随长度按 $\sqrt{d}$ 标度变化——腔越长，单位长度的场强越小）。

### 9.7 模态哈密顿量（式 13）→ 一组独立谐振子
将展开代入原哈密顿量，利用正交性可得：
$$H=\sum_{m=0}^{\infty}\left[\frac{Q_m^2}{2C_r}+\frac{1}{2}C_r\omega_m^2\Phi_m^2\right]$$
其中：
- $C_r=c_0d$（谐振器总电容）
- $Q_m=C_r\dot{\Phi}_m$（第 $m$ 模的广义电荷）
- **物理意义**：每个模式都是一个独立的谐振子，有自己有效的"质量"($C_r$)和"弹簧常数"($C_r\omega_m^2$)

### 9.8 量子化（式 14–16）
将 $\Phi_m,Q_m$ 升格为算符，满足 $[\hat{\Phi}_m,\hat{Q}_{m'}]=i\hbar\delta_{mm'}$，引入产生湮灭算符：
$$\hat{\Phi}_m=\sqrt{\frac{\hbar Z_m}{2}}(\hat a_m^\dagger+\hat a_m)$$
$$\hat{Q}_m=i\sqrt{\frac{\hbar}{2Z_m}}(\hat a_m^\dagger-\hat a_m)$$
其中模阻抗 $Z_m=\sqrt{L_m/C_r}$，$L_m^{-1}=C_r\omega_m^2$。

最终量子化哈密顿量（忽略零点能的紧凑形式）：
$$\hat H=\sum_{m=0}^{\infty}\hbar\omega_m\,\hat a_m^\dagger\hat a_m$$

### 9.9 文中强调的"简化与扩展"
- **简化假设**：常为简化忽略端口耦合电容对边界条件/模函数/频率的影响；但实际中耦合电容会：
  - 改变边界条件（从理想开路变为电容负载）
  - 改变模函数形状和频率（产生频率偏移）
  - 引入外部耗散（决定 $Q_{\text{ext}}$）
  
- **非线性扩展**：在中心导体中嵌入约瑟夫森结可引入 **Kerr 非线性**（$\propto (\hat a^\dagger\hat a)^2$ 项），用于：
  - 近量子极限参量放大
  - 分岔放大（bifurcation amplification）
  - 研究量子加热（quantum heating）
  - 产生非经典态（压缩态、猫态）
  
- 更一般的含非均匀/嵌入结的处理有专门文献（如 Bourassa 等），涉及如何将结作为集总非线性元件耦合到多模场中。

---

## 逻辑链总结

**为什么需要** → 谐振器是 circuit QED 的核心元件（量子存储器/读出接口）

**怎么构成** → CPW 结构（中心导体 + 间隙 + 地），通过边界条件形成驻波

**用什么描述** → 单位长度参数 $l_0,c_0$ → 特性阻抗 $Z_r$、波速 $v_0$ → 频率 $f_m$

**如何工作** → 边界条件决定模式频率 → 品质因数 $Q$ 描述能量保持能力 → $Q_{\text{int}}$（材料/工艺限制） vs $Q_{\text{ext}}$（耦合设计）

**为什么强耦合** → 小横向尺寸 → 大零点电场涨落 → 强耦合强度 $g$

**如何量子化** → 离散电报方程 → 连续场 → 波动方程 → 正常模分解 → 每模独立谐振子 → 升格为量子算符 → 得到光子数态
  
  

#### B. 3D 腔 (3D Cavities)

  

*   **结构**: 挖空的金属块（如矩形腔、同轴腔）。

*   **优势**: 极高的 Q 值 ($10^6 - 10^8$)。

*   **原因**: 电磁场主要分布在真空中，表面参与率极低，受表面 TLS 影响小。

*   **应用**: 作为长寿命的量子存储器。

  

#### 3. 引入非线性：约瑟夫森结 (Josephson Junction)

  

*   **为什么需要非线性？**

    *   LC 振荡器是**线性**的，能级间隔是等距的 ($\hbar \omega, 2\hbar \omega, \dots$)。

    *   如果你想从 $|0\rangle$ 激发到 $|1\rangle$，同样的脉冲会把你激发到 $|2\rangle, |3\rangle \dots$。你无法单独控制 $|0\rangle$ 和 $|1\rangle$ 作为一个**量子比特**。

    *   我们需要**非谐性**，使得 $\omega_{01} \neq \omega_{12}$。

  

*   **约瑟夫森结**:

    *   **结构**: 超导体-绝缘体-超导体。

    *   **物理**: 非耗散的非线性电感。

    *   **方程**: $I = I_c \sin \varphi$ (电流与相位的正弦关系，而非线性的 $I \propto \Phi$)。

    *   **势能**: $U(\varphi) = -E_J \cos \varphi$ (余弦势阱，而非 LC 的抛物线势阱)。

    *   正是这个余弦势阱提供了能级的不等间距，使得我们可以制造出 Transmon 等量子比特。

  

#### C. The transmon artificial atom
# Transmon Qubit（Transmon 人工原子）逻辑笔记

## 1. 动机：为什么需要 Transmon（非线性人工原子）？

### 1.1 线性谐振子不利于编码量子信息
即使上一节讨论的 LC 谐振器或 CPW 谐振器能被制备到量子基态，要在这类**线性系统**里实现量子信息处理仍然困难。线性谐振子的能级是等间距的：
$$E_n = \hbar\omega\left(n+\frac{1}{2}\right)$$
这意味着驱动 $|0\rangle \rightarrow |1\rangle$ 的脉冲同样会驱动 $|1\rangle \rightarrow |2\rangle$、$|2\rangle \rightarrow |3\rangle$……无法单独控制特定的两个能级作为一个量子比特。

### 1.2 量子信息需要非线性
为了在电路中编码与操控量子信息，需要**非谐性（anharmonicity）**——即能级间隔不等，使得 $\omega_{01} \neq \omega_{12} \neq \omega_{23} \dots$。这样我们可以用频率选择性地驱动特定跃迁，把系统近似看作一个二能级系统（量子比特）。

### 1.3 超导的优势：引入非线性但保持低损耗
超导材料本身是无损耗的（直流电阻为零），但线性电感（几何电感）只能产生等间距能级。要引入非线性同时保持低损耗，关键元件是**约瑟夫森结（Josephson junction, JJ）**。它是唯一能与高 Q 超导谐振器、毫开尔文温度兼容的非线性无源元件。

### 1.4 历史点
约瑟夫森结的物理最早由 Brian Josephson 于 1962 年理论预言，随后在实验中被证实，他也因此获得诺贝尔奖。如今，约瑟夫森结是几乎所有超导量子比特的核心元件。

---

## 2. Josephson 结的两条基本关系

### 2.1 无耗散超电流关系（DC Josephson 效应）
两个超导电极之间隔着薄绝缘层（~1-2 nm）时，库珀对可以通过量子隧穿效应穿过绝缘层，形成无耗散的超电流：
$$I = I_c \sin\varphi$$
其中：
- $\varphi = \varphi_2 - \varphi_1$ 是两个超导电极的宏观波函数相位差
- $I_c$ 是**临界电流**——结在库珀对被破坏前可承受的最大超电流，由结的尺寸、材料和温度决定

**物理意义**：即使结两端电压为零，也可以有持续的超电流流动（直流约瑟夫森效应）。电流与相位差的正弦关系是**非线性的根源**。

### 2.2 相位与电压关系（AC Josephson 效应）
当结两端有电压 $V$ 时，相位差随时间演化：
$$\frac{d\varphi}{dt} = \frac{2\pi}{\Phi_0}V$$
其中磁通量子（超导通用常数）：
$$\Phi_0 = \frac{h}{2e} \approx 2.07 \times 10^{-15} \text{ Wb}$$

**物理意义**：恒定电压会产生线性增长的相位，导致交流电流（交流约瑟夫森效应）。这类似于电感两端电压与电流的关系，但约瑟夫森结是非线性的。

常把相位与"广义磁通"$\Phi(t)$ 联系起来：
$$\varphi(t) = \frac{2\pi}{\Phi_0}\Phi(t) \quad (\text{mod } 2\pi),\qquad \Phi(t) = \int^t dt'\,V(t')$$
但需要注意：$\varphi$ 是**紧致变量**（$\varphi \equiv \varphi + 2\pi$，因为波函数的相位周期性），而 $\Phi$ 可取任意实数。

---

## 3. Josephson 电感：把 JJ 看作"非线性电感"

在小信号极限下（$I \ll I_c$），约瑟夫森结可以看作一个非线性电感。定义**约瑟夫森电感**（式 18）：
$$L_J(\Phi) \equiv \left(\frac{\partial I}{\partial \Phi}\right)^{-1} = \frac{\Phi_0}{2\pi I_c \cos(2\pi\Phi/\Phi_0)} = \frac{L_{J0}}{\cos\varphi}$$
其中 $L_{J0} = \Phi_0/(2\pi I_c)$ 是零偏置时的电感。

**与几何电感的关键区别**：
- 几何电感 $L$ 是常数，与电流/磁通无关
- 约瑟夫森电感 $L_J$ **随工作点（相位/磁通/电流）变化**——当 $\varphi \rightarrow \pi/2$ 时，$L_J \rightarrow \infty$（电流饱和）；当 $\varphi > \pi/2$ 时，$L_J$ 变为负值（对应"失稳"区域）

因此，约瑟夫森结在低于临界电流时可视为**非线性电感元件**——这正是产生非谐性的物理根源。

---

## 4. 能量视角：JJ 产生余弦势阱

### 4.1 线性电感能量 vs 约瑟夫森能量
- **线性电感**的储能公式：$E = \int I\,d\Phi = \frac{\Phi^2}{2L}$，对应**抛物线势阱**。
- **约瑟夫森结**的储能（式 19）：
  $$E = -\int I\,d\Phi = -\int I_c\sin\varphi \cdot \frac{\Phi_0}{2\pi}d\varphi = -E_J\cos\varphi + \text{常数}$$
  其中**约瑟夫森能量**定义为：
  $$E_J = \frac{\Phi_0 I_c}{2\pi}$$

**结论**：约瑟夫森结对应的是**余弦势阱**，与 LC 谐振子的二次势完全不同。图 5(a) 对比了余弦势与二次势的形状差异——余弦势在底部稍平坦，但两侧更陡，导致能级不等间距。

### 4.2 "人工原子"的由来
将 LC 谐振子的几何电感 $L$ 替换为约瑟夫森结（图 5(b)/(c) 的思路）后，电路变为**非线性**，能级不再等间隔。如果非谐性和品质因数足够好，能谱类似原子——能级可分辨且非均匀分布，因此称为**人工原子**。

实际用作量子比特时，通常只取最低两个能级（基态 $|g\rangle$ 与第一激发态 $|e\rangle$）编码量子信息，更高能级（$|f\rangle, |h\rangle, \dots$）可能会在门操作中引入泄漏错误，需要额外处理。

---

## 5. Transmon 的电路与哈密顿量（核心公式链）

### 5.1 电路结构：电容分流的约瑟夫森结
固定频率 transmon（图 5b）的核心结构：
- 一个约瑟夫森结（参数 $E_J$、结电容 $C_J$）
- 并联一个较大的分流电容 $C_S$

总电容：
$$C_\Sigma = C_J + C_S$$

**设计思想**：通过增加 $C_S$ 增大总电容，从而降低充电能 $E_C$，使系统进入 $E_J/E_C \gg 1$ 的 transmon 区。

### 5.2 变量定义
- **电荷数算符**（单位：库珀对数）：
  $$\hat n \equiv \frac{\hat Q}{2e}$$
  （$\hat Q$ 是结上的总电荷算符）

- **相位算符**：
  $$\hat\varphi \equiv \frac{2\pi}{\Phi_0}\hat\Phi \quad (\text{mod } 2\pi)$$
  （$\hat\Phi$ 是广义磁通算符）

- **充电能**（charging energy）：
  $$E_C \equiv \frac{e^2}{2C_\Sigma}$$
  这是单个电子（不是库珀对）的充电能——注意因子 2 的来源。

- **偏置电荷** $n_g$：来自与外界电荷源的电容耦合（杂散电荷或外加栅压）。通常写为 $n_g = Q_g/(2e)$，其中 $Q_g$ 是栅电容 $C_g$ 上的感应电荷。

### 5.3 Transmon 哈密顿量（式 20）
$$\hat H_T = \frac{(\hat Q - Q_g)^2}{2C_\Sigma} - E_J\cos\left(\frac{2\pi\hat\Phi}{\Phi_0}\right) = 4E_C(\hat n - n_g)^2 - E_J\cos\hat\varphi$$

**物理意义**：
- 第一项 $4E_C(\hat n - n_g)^2$：电容的充电能，类似于动能项
- 第二项 $-E_J\cos\hat\varphi$：约瑟夫森结的势能项
- 两者对易关系：$[\hat\varphi, \hat n] = i$（类似于位置与动量）

### 5.4 能谱的一般形式
无论参数区间如何，总能写成对角形式：
$$\hat H = \sum_j \hbar\omega_j |j\rangle\langle j|$$
常用本征态标记：$|g\rangle,|e\rangle,|f\rangle,|h\rangle,\dots$（有时也用 $|0\rangle,|1\rangle,|2\rangle,\dots$；需与谐振器 Fock 态区分时会特别说明）。

---

## 6. 关键控制参量：比值 $E_J/E_C$ 决定量子比特类型

$E_J/E_C$ 是区分不同类型超导电荷量子比特的关键参数。

### 6.1 小 $E_J/E_C$ 区（充电能主导，$E_J/E_C < 1$）
- 本征态近似是电荷数算符的本征态 $|n\rangle$（局域在特定电荷数）。
- 此时能级对 $n_g$ **极度敏感**：环境中不可避免的电荷涨落会引起跃迁频率大幅波动 → 严重的电荷噪声退相干（dephasing）。
- 这是早期 Cooper Pair Box 的工作区间，相干时间很短。

### 6.2 Transmon 区：大 $E_J/E_C$（典型值 20–80）
通过增大 $E_J/E_C$（即增大 $E_J$ 或减小 $E_C$）进入 transmon 区，电荷自由度在余弦势阱中变得**强烈离域**——波函数在相位空间扩展，电荷数不再是好量子数。

**结果**：最低几个能级的能量对 $n_g$ 变得几乎不敏感。图 6 展示了 $E_J/E_C = 2,10,50$ 时能级随 $n_g$ 的起伏变化——比值越大，曲线越"平坦"。

**重要补充**：即使能级对静态 $n_g$ 不敏感，我们仍然可以用外部交流电压源（通过电容耦合）驱动 transmon 的态间跃迁——这相当于在 $n_g$ 上施加一个交流调制。

### 6.3 代价与收益的权衡
增大 $E_J/E_C$ 带来两个相反的效果：
- **收益**：电荷色散（charge dispersion）随 $E_J/E_C$ **指数级降低**（图中给出 $\sim e^{-\sqrt{8E_J/E_C}}$ 量级），相干时间大幅延长
- **代价**：非谐性（anharmonicity）随 $E_J/E_C$ 增大而**缓慢降低**，约为 $\sim (E_J/E_C)^{-1/2}$ 量级

**权衡结果**：选择 $E_J/E_C \sim 20-80$ 可在相干时间（指数增益）和非谐性（幂律损失）之间取得最佳平衡。

---

## 7. 在 Transmon 区的近似：从余弦势到弱非简谐振子

当 $E_J/E_C \gg 1$ 时，相位在势阱底部较局域（但不像 Cooper Pair Box 那样局域在电荷数空间），可以把余弦势在最小值附近做泰勒展开。

### 7.1 去掉 $n_g$ 并重写哈密顿量
在 transmon 区，相关低能级的频率对 $n_g$ 几乎不敏感，因此常将 $n_g$ 略去，直接写为：
$$\hat H_q \approx 4E_C\hat n^2 - E_J\cos\hat\varphi$$

### 7.2 余弦展开后的哈密顿量（式 22）
将 $\cos\hat\varphi$ 在 $\hat\varphi=0$ 附近展开到四阶：
$$\cos\hat\varphi \approx 1 - \frac{\hat\varphi^2}{2} + \frac{\hat\varphi^4}{24} + \mathcal{O}(\hat\varphi^6)$$
代入并忽略常数项，得到：
$$\hat H_q \approx 4E_C\hat n^2 + \frac{1}{2}E_J\hat\varphi^2 - \frac{1}{24}E_J\hat\varphi^4$$

**物理意义**：transmon 是一个**弱非简谐振子**——主导项是二次型（简谐振子），非线性来自 $-\varphi^4$ 项（负号表示非谐性为负）。

> **重要警告**（文中脚注）：这种截断近似会导致近似哈密顿量在数学上"下无界"（当 $\hat\varphi$ 很大时，$-\varphi^4$ 项会使能量趋于 $-\infty$），这是截断余弦级数导致的伪影。使用时应在原希尔伯特空间的**截断子空间**内谨慎应用，不能外推到高激发态。

### 7.3 引入产生湮灭算符（式 23-24）
选取算符 $\hat b, \hat b^\dagger$ 使二次部分对角化：
$$\hat\varphi = \left(\frac{2E_C}{E_J}\right)^{1/4}(\hat b^\dagger + \hat b)$$
$$\hat n = \frac{i}{2}\left(\frac{E_J}{2E_C}\right)^{1/4}(\hat b^\dagger - \hat b)$$

**重要结论**：当 $E_J/E_C$ 增大时：
- $\hat\varphi$ 的零点涨落 $\propto (E_C/E_J)^{1/4}$ **减小**（相位更局域）
- 共轭变量 $\hat n$ 的零点涨落 $\propto (E_J/E_C)^{1/4}$ **增大**（电荷涨落更大）

这正是 transmon 对电荷噪声不敏感的原因——电荷涨落大意味着电荷数不确定，因此外部电荷扰动影响小。

### 7.4 得到 Kerr 形式的有效哈密顿量（式 25）
代入并在旋转波近似下保留"数算符守恒"项（即保留 $\hat b^\dagger\hat b^\dagger\hat b\hat b$ 形式的项，舍弃 $\hat b^\dagger\hat b^\dagger + \hat b\hat b$ 等非能量守恒项）：
$$\hat H_q \approx \sqrt{8E_CE_J}\,\hat b^\dagger \hat b - \frac{E_C}{12}(\hat b^\dagger+\hat b)^4 \approx \hbar\omega_q\,\hat b^\dagger \hat b - \frac{E_C}{2}\hat b^\dagger\hat b^\dagger\hat b\hat b$$

其中：
$$\hbar\omega_q \approx \sqrt{8E_CE_J} - E_C$$

**旋转波近似有效的条件**：$\hbar\omega_q \gg E_C$，这在 transmon 区（$E_J/E_C \gg 1$）很容易满足。

### 7.5 三个必须记住的物理量（transmon 区）
1. **约瑟夫森等离子体频率**（plasma frequency）：
   $$\omega_p = \frac{\sqrt{8E_CE_J}}{\hbar}$$
   对应余弦势阱底部"小振动"的角频率——即如果把余弦近似为抛物线得到的简谐振子频率。

2. **基态到第一激发态的跃迁频率**：
   $$\omega_{ge} \approx \omega_p - \frac{E_C}{\hbar}$$
   比等离子体频率略低，偏移量正好是 $E_C/\hbar$。

3. **非谐性**（anharmonicity）：
   $$\alpha \equiv \omega_{ef} - \omega_{ge} \approx -\frac{E_C}{\hbar}$$
   典型量级：$E_C/h \sim 100\text{–}400$ MHz，所以 $\alpha/2\pi \sim -100$ 到 $-400$ MHz。

### 7.6 为什么"非线性虽小但仍够用"
尽管 $E_C \ll \hbar\omega_q$（即非谐性远小于跃迁频率），但 $|\alpha|$ 通常仍在 $100$ MHz 量级，而量子比特的线宽 $\gamma$ 通常在 kHz-MHz 量级。因此：
$$|\alpha| \gg \gamma$$
能级不等间距仍远大于谱线宽度，从而可以用频率选择性地驱动 $|g\rangle \leftrightarrow |e\rangle$ 跃迁而不激发到 $|f\rangle$。

**同时要牢记**：transmon 本质上是多能级系统，在一些高保真度门操作或特定实验中，$|f\rangle$ 态的布居可能会成为限制因素，需要额外考虑。

---

## 8. Flux-tunable Transmon：用 SQUID 实现频率可调

### 8.1 电路思想
用一个 SQUID（两个并联的约瑟夫森结）替代单个约瑟夫森结，使等效 $E_J$ 受外加磁通 $\Phi_x$ 控制，从而量子比特频率可调。图 5(c) 展示了这种结构。

### 8.2 SQUID 情况的哈密顿量（式 26）
$$\hat H_T = 4E_C\hat n^2 - E_{J1}\cos\hat\varphi_1 - E_{J2}\cos\hat\varphi_2$$

外加磁通 $\Phi_x$ 穿过 SQUID 环路，且忽略环路几何电感时（即环路电感远小于约瑟夫森电感），有磁通量子化条件：
$$\hat\varphi_1 - \hat\varphi_2 = 2\pi\frac{\Phi_x}{\Phi_0} \quad (\text{mod } 2\pi)$$

### 8.3 约化为单结形式
定义平均相位 $\hat\varphi = (\hat\varphi_1 + \hat\varphi_2)/2$，通过三角恒等式可化为（式 27）：
$$\hat H_T = 4E_C\hat n^2 - E_J(\Phi_x)\cos(\hat\varphi - \varphi_0)$$

### 8.4 有效约瑟夫森能量（式 28）
$$E_J(\Phi_x) = E_{J\Sigma} \cos\left(\pi\frac{\Phi_x}{\Phi_0}\right) \sqrt{1 + d^2 \tan^2\left(\pi\frac{\Phi_x}{\Phi_0}\right)}$$
其中：
- $E_{J\Sigma} = E_{J1} + E_{J2}$（总约瑟夫森能量）
- $d = \frac{E_{J2} - E_{J1}}{E_{J\Sigma}}$（结的不对称度）

相位偏置 $\varphi_0 = \arctan\!\left(d\tan(\pi\Phi_x/\Phi_0)\right)$ 对**静态磁通**通常可忽略（因仅影响波函数相位，不改变能谱）。

**对称结特例**（$E_{J1} = E_{J2} = E_J$，$d=0$）：
$$E_J(\Phi_x) = 2E_J \left|\cos\left(\pi\frac{\Phi_x}{\Phi_0}\right)\right|$$

### 8.5 可调频率与调谐速度
量子比特频率随磁通变化：
$$\omega_q(\Phi_x) \approx \frac{\sqrt{8E_C E_J(\Phi_x)} - E_C}{\hbar}$$

**实验参数**：
- 调谐范围：约 **1-2 GHz**（典型值）
- 调谐速度：可快到 **10-20 ns**（通过快速磁通脉冲）
- 动态范围可通过增加磁通线带宽扩展

这种可调谐性用于多种应用：
- 快速将量子比特调到与谐振器共振（实现 iSWAP 门）
- 调到特定避免交叉点（实现 CZ 门）
- 调到远离谐振器（减少耦合，保护相干性）

### 8.6 代价：磁通噪声退相干
额外的调谐旋钮会引入磁通噪声导致的退相干风险。磁通噪声来源包括：
- 外部环境中的随机磁通涨落
- 临近涡旋的运动
- 电流源的噪声

**不对称性的权衡**：
- 更大的结不对称 $d$ → 调谐范围减小（$E_J(\Phi_x)$ 的最小值变大）
- 但也会使器件在磁通偏置点对磁通噪声**更不敏感**（因为 $dE_J/d\Phi_x$ 更小）

这是一个典型的"可调谐性 vs 噪声敏感性"的权衡。

### 8.7 拓展方向
文中提到实现"电压可调"的 transmon：例如用**半导体纳米线**替代 SQUID 环路，通过栅压控制结的传输特性，实现电学可调的约瑟夫森能量。

---

## 9. 图像与直观总结（对应图 5 与图 6）

### 图 5 的核心信息
- **图 5(a)**：transmon 的余弦势阱（实线）与 LC 的二次势（虚线）对比。余弦势底部稍平坦，但两侧更陡，导致能级 $|g\rangle,|e\rangle,|f\rangle,|h\rangle,\dots$ 不等间隔，非谐性量级约为 $-E_C$。
- **图 5(b)**：固定频率 transmon 的电路图——单个约瑟夫森结并联分流电容。
- **图 5(c)**：频率可调 transmon——用 SQUID（两个结）替代单结，通过外磁通控制等效 $E_J$。

### 图 6 的核心信息
展示不同 $E_J/E_C$ 比值下，最低三个能级对 offset charge $n_g$ 的敏感性：
- $E_J/E_C = 2$：能级剧烈振荡，对电荷噪声极度敏感（早期 Cooper Pair Box）
- $E_J/E_C = 10$：振荡幅度显著减小
- $E_J/E_C = 50$：能级几乎平坦，对 $n_g$ 不敏感（典型 transmon 区）

**结论**：比值越大，电荷色散被压制得越彻底——这正是 transmon 获得长相干时间的关键。

---

## 10. Transmon 设计参数总结表

| 参数 | 符号 | 典型值 | 作用 |
|------|------|--------|------|
| 充电能 | $E_C/h$ | 100-400 MHz | 决定非谐性 $\alpha \approx -E_C/\hbar$ |
| 约瑟夫森能 | $E_J/h$ | 10-30 GHz | 决定跃迁频率 $\omega_{ge} \approx \sqrt{8E_CE_J}/\hbar$ |
| 比值 | $E_J/E_C$ | 20-80 | 决定电荷噪声敏感性（指数级影响） |
| 跃迁频率 | $\omega_{ge}/2\pi$ | 4-8 GHz | 量子比特的工作频率 |
| 非谐性 | $\alpha/2\pi$ | -100 到 -400 MHz | 确保能级可分辨，$\alpha \gg \gamma$ |
| 电荷色散 | $\epsilon$ | kHz-MHz 量级 | 随 $E_J/E_C$ 指数减小 |
| 调谐范围（可调版）| $\Delta\omega/2\pi$ | 1-2 GHz | 通过磁通调节频率的范围 |

---

## 逻辑链总结

**为什么需要** → 线性谐振子不能做量子比特（能级等距）→ 需要非线性 → 约瑟夫森结提供非线性（余弦势）

**如何工作** → 电荷能 $4E_C(\hat n - n_g)^2$ + 约瑟夫森能 $-E_J\cos\hat\varphi$ → 哈密顿量

**关键参数** → $E_J/E_C$ 比值决定量子比特类型：小比值→电荷敏感（相干差），大比值→电荷不敏感（transmon）

**transmon 区近似** → $E_J/E_C \gg 1$ → 弱非简谐振子近似 → $\hat H \approx \hbar\omega_q \hat b^\dagger \hat b - (E_C/2)\hat b^\dagger\hat b^\dagger\hat b\hat b$

**三个关键频率** → $\omega_p = \sqrt{8E_CE_J}/\hbar$（等离子体频率），$\omega_{ge} \approx \omega_p - E_C/\hbar$（工作频率），$\alpha \approx -E_C/\hbar$（非谐性）

**可调版本** → SQUID 替代单结 → $E_J(\Phi_x)$ 受磁通控制 → 频率可调，但引入磁通噪声风险

**最终结果** → transmon 实现了"长相干时间"（通过电荷噪声指数压制）与"足够非谐性"（$|\alpha| \gg \gamma$）的平衡，成为当前超导量子计算的主流量子比特类型






## **III. LIGHT–MATTER INTERACTION IN CIRCUIT QED**  
（电路 QED 中的光—物质相互作用）

---

# **A. Exchange interaction between a transmon and an oscillator**  
（transmon 与谐振子之间的交换相互作用）

## **A1. 场景与耦合方式（从“经典门电压”到“量子谐振子”）**
- 前文已引入两大主角：  
  1) **量子简谐振子**（微波谐振器的模式）  
  2) **transmon 人工原子**
- transmon 由于需要 **较大的电荷能（即较大电容）**，因此很自然地可通过**电容耦合**接入微波谐振器（见 Fig.7 的示意）。
- 当谐振器在电路中扮演“经典电压源 $V_g$”的角色时，可以在 transmon 哈密顿量中把**经典门电压**对应的 **$n_g$** 替换成一个**量子化的门电荷偏置**：
  $$
  n_g \ \rightarrow\ -\hat n_r
  $$
  其中 $-\hat n_r$ 表示由谐振器引入、作用在 transmon 上的量子化电荷偏置（符号正负是文献中的常用约定；文中说明采用该约定，见 Appendix A）。

## **A2. 多模情况下的总哈密顿量（Eq. 29）**
耦合后的总系统哈密顿量（Blais et al., 2004）写成：
$$
\hat H
=4E_C(\hat n+\hat n_r)^2
- E_J\cos\hat\phi
+\sum_m \hbar\omega_m \hat a_m^\dagger \hat a_m,
\tag{29}
$$
并给出：
- $\hat n$：transmon 的库珀对数（电荷数）算符  
- $\hat\phi$：结两端超导相位差算符  
- $E_C$：充电能；$E_J$：约瑟夫森能  
- 第 $m$ 个谐振器模：频率 $\omega_m$，产生/湮灭算符 $\hat a_m^\dagger,\hat a_m$
- 谐振器引入的“门电荷”项：
  $$
  \hat n_r=\sum_m \hat n_m,\qquad
  \hat n_m=\Big(\frac{C_g}{C_m}\Big)\frac{\hat q_m}{2e},
  $$
  其中  
  - $C_g$：耦合电容（gate capacitance）  
  - $C_m$：第 $m$ 个模式对应的“模式电容”  
  - $\hat q_m$：第 $m$ 模对 transmon 产生门偏置的电荷坐标（以电荷形式出现）  
- 文中为简化使用了假设：  
  $$
  C_g \ll C_\Sigma,\ C_m
  $$
  其中 $C_\Sigma$ 是 transmon 的总电容（有效电容）。  
- 文中说明：对 Eq.(29) 的更完整推导（不仅是简单替换 $n_g\to -\hat n_r$，且不依赖上述小电容近似）可见 Appendix A（单个 LC 振子耦合 transmon 的情形）。

## **A3. 单模近似与等效电路解释（Fig.7 与“black-box quantization”）**
- 若 transmon 频率更接近某一个谐振器模式而远离其他模式，满足近似条件（原文表述）：
  $$
  |\omega_0-\omega_q|\ll|\omega_m-\omega_q|\quad (m\ge 1),
  $$
  则可以将 Eq.(29) 中对 $m$ 的求和**截断**为单一项（single-mode approximation）。
- 在单模近似下，总系统可视为：**频率为 $\omega_r$ 的单一简谐振子 + transmon**。
- 文中特别强调：不论这个“振子模式”在物理上来自  
  - 2D/3D 分布参数谐振器的某一模，或  
  - 集总元件 LC  
  都可以用**等效集总元件电路**表示：transmon 电容耦合到 LC（Fig.7(b)）。
- 这种把复杂几何结构“等效成集总元件”的形式化量子化思路，称为 **black-box quantization**（Nigg et al., 2012），文中指出会在 Sec. III.D 更详细讨论。
- 文中同时提醒：在很多实验相关情形里，忽略谐振器的**多模性**会导致不准确预测（将在 Sec. IV.E 讨论）。

## **A4. 单模下、用产生湮灭算符写出的耦合哈密顿量（Eq. 30）**
在单模近似并引入前文的产生/湮灭算符后（谐振器用 $\hat a$，transmon 的弱非谐性用 $\hat b$ 表示），得到：
$$
\hat H \approx
\hbar\omega_r \hat a^\dagger\hat a
+\hbar\omega_q \hat b^\dagger\hat b
-\frac{E_C}{2}\hat b^\dagger\hat b^\dagger \hat b\hat b
-\hbar g(\hat b^\dagger-\hat b)(\hat a^\dagger-\hat a).
\tag{30}
$$
- $\omega_r$：所选谐振器模式频率  
- $\omega_q$：transmon 的（近似）本征频率  
- $-\frac{E_C}{2}\hat b^\dagger\hat b^\dagger \hat b\hat b$：transmon 的 Duffing 型非线性（弱非谐性）
- 最后项是“电偶极—电场”型的耦合（以 $(\hat b^\dagger-\hat b)(\hat a^\dagger-\hat a)$ 形式出现）

**脚注 3（原文要点）**  
- 可能担心 Eq.(29) 展开后出现 $\hat n_r^2$ 项；文中说明该项可并入谐振器模式的“充电能”从而导致谐振器频率的重整化，因此为简洁起见省略；更多细节见文中引用的附录公式（A9）、（A10）。

## **A5. 旋波近似（RWA）与交换型耦合（Eq. 32）**
在实验常见的弱耦合条件下：
$$
|g|\ll \omega_r,\ \omega_q,
$$
采用 rotating-wave approximation（旋波近似）后，Eq.(30) 化为：
$$
\hat H \approx
\hbar\omega_r \hat a^\dagger\hat a
+\hbar\omega_q \hat b^\dagger\hat b
-\frac{E_C}{2}\hat b^\dagger\hat b^\dagger \hat b\hat b
+\hbar g(\hat b^\dagger \hat a+\hat b \hat a^\dagger).
\tag{32}
$$
- 最后一项 $\hbar g(\hat b^\dagger \hat a+\hat b \hat a^\dagger)$ 明确对应“交换（exchange）相互作用”：在两系统间相干交换一个量子。

## **A6. 耦合强度 $g$ 的表达式（Eq. 31）及参数含义**
文中给出振子—transmon（光—物质）耦合强度：
$$
g=\omega_r\frac{C_g}{C_\Sigma}\Big(\frac{E_J}{2E_C}\Big)^{1/4}\sqrt{\frac{\pi Z_r}{R_K}}.
\tag{31}
$$
各量定义：
- $Z_r$：谐振器模式的特征阻抗（characteristic impedance）
- $R_K=h/e^2 \approx 25.8~\text{k}\Omega$：电阻量子（resistance quantum）
- $\frac{C_g}{C_\Sigma}$：耦合“分压/分电荷”比例因子
- $\Big(\frac{E_J}{2E_C}\Big)^{1/4}$：与 transmon 的电荷涨落尺度相关（文中指出可由 Eq.(24) 看出这一点）

## **A7. 用“偶极矩 × 零点电场”的物理图像理解 $g$**
- 文中引入一个长度尺度 $l$：表示库珀对穿越 transmon 结隧穿时的等效距离。由此可把 Eq.(31) 诱导性地解释为：
  $$
  \hbar g = d_0 E_0,
  $$
  其中  
  - $d_0=2el\,(E_J/32E_C)^{1/4}$：transmon 的电偶极矩  
  - $E_0=(\omega_r/l)(C_g/C_\Sigma)\sqrt{\hbar Z_r/2}$：谐振器的零点电场（以 transmon “看到”的方式定义）
- 结论：因为 $d_0$ 与 $E_0$ 都可以做得很大（尤其在 transmon 区域 $E_J/E_C$ 大时），所以电路 QED 中的电偶极耦合强度可远大于天然原子在腔 QED 中的耦合。

## **A8. 进一步写成与精细结构常数相关的形式（Eq. 33）与工程提升路径**
文中进一步把 Eq.(31) 改写成与精细结构常数有关的形式（Devoret et al., 2007）：
$$
g=\omega_r\frac{C_g}{C_\Sigma}\Big(\frac{E_J}{2E_C}\Big)^{1/4}
\sqrt{\frac{Z_r}{Z_{\rm vac}}}\,\sqrt{2\pi\alpha},
\tag{33}
$$
其中
- $\alpha=Z_{\rm vac}/(2R_K)$：精细结构常数
- $Z_{\rm vac}=\sqrt{\mu_0/\varepsilon_0}\approx 377~\Omega$：真空阻抗
- 文中指出：这里“相互作用强度”被 $\frac{Z_r}{Z_{\rm vac}}$ 与 $\frac{C_g}{C_\Sigma}$ 这类小于 1 的因子削弱；但仍可通过增大 $E_J/E_C$（transmon 区域）实现很大的耦合。
- 代价：更大的 $g$ 往往以减小 transmon 的相对非谐性为代价；文中给出相对非谐性标度（原文表达）：
  $$
  -\frac{E_C}{\hbar\omega_q}\sim \sqrt{\frac{E_C}{8E_J}}.
  $$
- 工程上提升 $g$ 的方式之一：**提高谐振器阻抗 $Z_r$**，例如用结阵列替代谐振器中心导体（Andersen and Blais, 2017; Stockklauser et al., 2017）。
- 从 2D 到 3D：3D 谐振腔使模式体积显著变大，从而降低真空电场涨落；但可通过把 transmon 做大（增加偶极矩）在不改变 $g$ 数量级的情况下补偿（Paik et al., 2011）。文中形象描述：transmon 在 3D 腔中相当于“天线”，放置在合适位置可强耦合到某一腔模（Fig.7(c)）。

## **A9. 两能级近似得到 Jaynes–Cummings 哈密顿量（Eq. 34）**
为了进一步强化与腔 QED 的类比，文中将 transmon 限制到前两能级 $\{|g\rangle,|e\rangle\}$，做替换：
$$
\hat b^\dagger\rightarrow \hat\sigma_+ = |e\rangle\langle g|,\qquad
\hat b\rightarrow \hat\sigma_- = |g\rangle\langle e|.
$$
于是从 Eq.(30) 得到著名的 Jaynes–Cummings（JC）模型（Blais et al., 2004; Haroche and Raimond, 2006）：
$$
\hat H_{\rm JC}
=\hbar\omega_r\hat a^\dagger\hat a
+\frac{\hbar\omega_q}{2}\hat\sigma_z
+\hbar g(\hat a^\dagger\hat\sigma_-+\hat a\hat\sigma_+),
\tag{34}
$$
其中采用
$$
\hat\sigma_z=|e\rangle\langle e|-|g\rangle\langle g|.
$$
- 最后一项描述“单量子”的相干交换：一个谐振器光子 ↔ transmon 的一个激发。

---

# **B. The Jaynes–Cummings spectrum**  
（Jaynes–Cummings 能谱）

## **B1. 模型地位与适用性（文中评价）**
- JC 哈密顿量是**可精确求解**的模型，能准确描述许多情形：把（天然或人工）原子当作两能级系统，与电磁场的单一模式相互作用。
- 在只关心 transmon 前两能级 $|\sigma\rangle=\{|g\rangle,|e\rangle\}$ 的实验里，JC 模型常能给出**定性一致**的结果。
- 但文中也强调：要达到**定量一致**，通常必须考虑  
  1) transmon 的更高能级  
  2) 场的多模性  
  这点在电路 QED 中经常很关键。尽管如此，JC 模型仍提供大量重要物理直觉，因此接下来聚焦其谱结构。

## **B2. 裸态（bare states）与“总激发数”守恒结构**
- 当 $g=0$（无耦合）时，qubit–field 系统的裸态标记为 $|\sigma,n\rangle$，其中  
  - $\sigma\in\{g,e\}$ 表示两能级系统状态  
  - $n$ 表示谐振器的光子数
- JC 模型具有按“总激发数”分块对角的结构（下文引入算符 $\hat N_T$）。

## **B3. 用 Bogoliubov-like 幺正变换对角化（Eq. 35–37）**
文中给出 dressed 态可由幺正变换从裸态获得：
$$
|\overline{\sigma,n}\rangle=\hat U^\dagger|\sigma,n\rangle,
$$
其中变换取（Boissonneault et al., 2009; Carbonaro et al., 1979）：
$$
\hat U=\exp\!\Big[\Lambda(\hat N_T)\big(\hat a^\dagger\hat\sigma_- - \hat a\hat\sigma_+\big)\Big],
\tag{35}
$$
并定义
$$
\Lambda(\hat N_T)=\frac{\arctan\!\big(2\lambda\sqrt{\hat N_T}\big)}{2\sqrt{\hat N_T}}.
\tag{36}
$$
这里
- $$
  \hat N_T=\hat a^\dagger\hat a+\hat\sigma_+\hat\sigma_-
  $$
  是**总激发数算符**；文中指出它与 $\hat H_{\rm JC}$ 对易。
- $$
  \lambda=\frac{g}{\Delta},\qquad \Delta=\omega_q-\omega_r
  $$
  其中 $\Delta$ 是 qubit–resonator 失谐（detuning）。

对角化后：
$$
\hat H_D=\hat U^\dagger \hat H_{\rm JC}\hat U
=\hbar\omega_r\hat a^\dagger\hat a
+\frac{\hbar\omega_q}{2}\hat\sigma_z
-\frac{\hbar\Delta}{2}\Big(1-\sqrt{1+4\lambda^2\hat N_T}\Big)\hat\sigma_z.
\tag{37}
$$

## **B4. JC 能谱的“双态（doublets）”结构与本征能量（Eq. 38）**
由 Eq.(37) 可直接读出 dressed 态能量；JC 谱由固定激发数的双重态（doublet）
$\{|g,n\rangle,\ |e,n-1\rangle\}$ 构成（文中脚注 4 说明：推导中给 $\hat H_D$ 加了 $\hbar\omega_r/2$ 的整体能量平移，这不影响物理结论）。

能量为：
$$
E_{g,n}=n\hbar\omega_r-\frac{\hbar}{2}\sqrt{\Delta^2+4g^2n},
\qquad
E_{e,n-1}=n\hbar\omega_r+\frac{\hbar}{2}\sqrt{\Delta^2+4g^2n}.
\tag{38}
$$
并给出基态 $|g,0\rangle$ 的能量：
$$
E_{g,0}=-\frac{\hbar\omega_q}{2}.
$$

## **B5. dressed 本征态（混合角形式，Eq. 39）**
对应的激发态 dressed 态为：
$$
|\overline{g,n}\rangle=\cos(\theta_n/2)|g,n\rangle-\sin(\theta_n/2)|e,n-1\rangle,
$$
$$
|\overline{e,n-1}\rangle=\sin(\theta_n/2)|g,n\rangle+\cos(\theta_n/2)|e,n-1\rangle,
\tag{39}
$$
其中混合角
$$
\theta_n=\arctan\!\Big(\frac{2g\sqrt n}{\Delta}\Big).
$$

## **B6. 共振时的 $\sqrt n$ 分裂与实验意义（结合 Fig.8 的说明）**
- 能谱的关键特征是对光子数 $n$ 的标度：在**零失谐**（共振）$\Delta=0$ 时，
  - doublet 能级分裂为
    $$
    2g\sqrt n.
    $$
- 文中强调：这不同于两个耦合的**谐振子**（harmonic oscillators），后者的分裂与 $n$ 无关；因此 **$2g\sqrt n$** 是 JC 量子特征之一。
- 实验上测量该谱可用于评估耦合系统的量子性质（Carmichael et al., 1996; Fink et al., 2008）；文中指出将于 Sec. VI.A 回到相关讨论。
- Fig.8（文中描述要点）：
  - 灰线：uncoupled（未耦合）能级  
  - 蓝线：dressed（耦合后）能级  
  - 在 $\Delta=0$ 时，$|g\rangle,|e\rangle$ 与光子数 $n=0,1,2,\dots$ 形成的简并被 $2g\sqrt n$ 劈裂
  - 图外浅蓝线标出 transmon 的第三激发态 $|f\rangle$（虽未纳入 JC 两能级模型）；文中指出即便未显式画出，该能级的存在会在 $n\ge 2$ 的 manifold 中引起 dressed 态能量位移（Fink et al., 2008）——这是“必须考虑 transmon 多能级性”的直接提示。

---

### **整体逻辑总结**

**III.A 部分** 从物理图像和形式推导两个层面，建立了 transmon 与谐振器之间交换相互作用的完整描述。  
- 首先，将经典门电压替换为谐振器量子化电荷偏置，得到多模耦合哈密顿量 (29)。  
- 然后通过单模近似和旋波近似，简化得到交换型哈密顿量 (32)，并给出耦合强度 $g$ 的表达式 (31)，展示了其与电路参数（阻抗、电容比、$E_J/E_C$ 等）的关系。  
- 进一步通过偶极矩-零点电场图像 (A7) 和精细结构常数形式 (33)，揭示了电路 QED 可实现超强耦合的物理根源及工程提升途径。  
- 最后，将 transmon 截断为两能级系统，得到 Jaynes-Cummings 模型 (34)，为后续能谱分析奠定基础。

**III.B 部分** 则聚焦于 JC 模型的能谱。  
- 首先指出 JC 模型的适用性和局限性，强调在电路 QED 中需考虑多能级和多模效应才能达到定量一致。  
- 接着利用总激发数守恒，通过幺正变换 (35)-(37) 对角化哈密顿量，得到 dressed 态能量 (38) 和本征态 (39)。  
- 重点分析了共振时的 $\sqrt{n}$ 依赖的能级分裂，这是 JC 模型区别于耦合谐振子的关键量子特征，也是实验上检验系统量子行为的重要依据。  
- 最后结合 Fig.8 的说明，指出了更高能级（如 $|f\rangle$）存在所带来的影响，为后续讨论更精确的模型埋下伏笔。

整体上，III.A 完成了从实际电路到抽象 JC 模型的映射，III.B 则深入剖析了 JC 模型的能谱结构，两者共同构成了电路 QED 中光-物质相互作用的理论基础。


## **III.C Dispersive regime（色散/非共振工作区）**

> 核心思想：当 qubit–resonator **失谐足够大**、交换过程不再共振时，光子与人工原子主要通过**虚跃迁**相互作用，表现为“频率拉拽/状态依赖频移”、以及更高阶的 Kerr 非线性。

### **C0. 进入色散区的条件与物理图像（承接前文）**
- 色散条件（前文已给出）：  
  $$
  |\lambda|=\left|\frac{g}{\Delta}\right|\ll 1,\qquad \Delta=\omega_q-\omega_r
  $$
- 此时 Eq.(32) 的交换项（$\hbar g(\hat b^\dagger\hat a+\hat b\hat a^\dagger)$）不再导致真实能量交换，耦合主要通过**虚光子过程**产生能级位移（Lamb shift）与腔频拉拽（dispersive shift）。

---

### **C1. Schrieffer–Wolff（SW）方法得到的一般色散哈密顿量（Eq. 40–41，承接你上一张图的末尾）**
从 Duffing 形式的 transmon–resonator 哈密顿量 Eq.(32) 出发，在色散极限做二阶 SW 变换，得到（文中给出的近似形式）：
$$
\hat H_{\rm disp}\approx 
\hbar\omega_r\hat a^\dagger\hat a
+\hbar\omega_q\hat b^\dagger\hat b
-\frac{E_C}{2}\hat b^\dagger\hat b^\dagger\hat b\hat b
+\sum_{j=0}^{\infty}\hbar\left(\Lambda_j+\chi_j\hat a^\dagger\hat a\right)|j\rangle\langle j|,
\tag{40}
$$
其中 $|j\rangle$ 为 transmon 的能级本征态标记（此处在 Eq.(25) 的近似下，可视作 $\hat b^\dagger\hat b$ 的数态本征态）。

系数关系（原文 Eq. 41a,b）：
$$
\Lambda_j=\chi_{j-1,j},
\qquad
\chi_j=\chi_{j-1,j}-\chi_{j,j+1},
\tag{41a}
$$
$$
\chi_{j-1,j}=\frac{jg^2}{\Delta-(j-1)E_C/\hbar}.
\tag{41b}
$$

补充说明（你给的 C 开头那段文字）：
- 对 $j>0$：$\Lambda_0=0$，且 $\chi_0=-g^2/\Delta$。
- $\chi_j$ 通常称为 **dispersive shifts（色散位移）**。
- $\Lambda_j$ 称为 **Lamb shifts**，体现真空涨落效应（Bethe 1947；Fragner 2008；Lamb & Retherford 1947）。

---

### **C2. 截断到两能级 → 标准两能级色散哈密顿量（Eq. 42）**
把 Eq.(40) 截断到 transmon 前两能级（$|g\rangle,|e\rangle$）得到更常见的形式（Blais et al., 2004）：
$$
\hat H_{\rm disp}\approx
\hbar\omega_r'\hat a^\dagger\hat a
+\frac{\hbar\omega_q'}{2}\hat\sigma_z
+\hbar\chi\,\hat a^\dagger\hat a\,\hat\sigma_z.
\tag{42}
$$
- $\chi$：qubit 状态依赖的腔频移动（cavity pull），决定读出时腔响应随 qubit 态变化的差异。

---

### **C3. 色散区中“可观测的”重整化频率与 $\chi$（Eq. 43）**
Eq.(42) 中参数由（Koch et al., 2007）：
$$
\omega_r'=\omega_r-\frac{g^2}{\Delta-E_C/\hbar},
\qquad
\omega_q'=\omega_q+\frac{g^2}{\Delta},
\tag{43a}
$$
$$
\chi=-\frac{g^2(E_C/\hbar)}{\Delta(\Delta-E_C/\hbar)}.
\tag{43b}
$$

原文强调的实验含义：
- **色散区实验测到的是 dressed（重整化）频率**（即 $\omega_r',\omega_q'$ 以及由 $\chi$ 导致的拉拽）。
- 但 Eq.(43) 右端出现的 $\omega_r,\omega_q$ 是 **bare（裸）频率**。

> **Fig.9（原文图注要点）**
> - 灰线：未耦合能级；蓝线：色散区 dressed 能级。
> - qubit 的 $|g\rangle\leftrightarrow|e\rangle$ 过渡在腔内光子数 $n$ 下被 Lamb shift 并进一步被 $\omega_q+\chi n$ “拉拽”。
> - 腔频率也被 qubit 态拉拽为 $\omega_r\pm\chi$（态依赖）。

---

### **C4. 更高阶效应：自 Kerr（self-Kerr）与四阶非线性**
- 原文指出：SW 变换不仅给出 Eq.(42)/(43) 的频移，还会在**四阶**产生  
  - resonator self-Kerr  
  - qubit self-Kerr  
  等非线性（Zhu et al., 2013）。

---

### **C5. 近似有效性的更精确判据：临界光子数 $n_{\rm crit}$（Eq. 44）**
虽然常用条件是 $|g/\Delta|\ll 1$，但由于耦合矩阵元会随腔内光子数、qubit 激发数增长，文中给出更精确（更保守）的有效性界：
$$
\bar n \ll n_{\rm crit}
=\frac{1}{2j+1}\left(\frac{|\Delta-jE_C/\hbar|^2}{4g^2}-j\right),
\tag{44}
$$
- 这里 $j=0,1,\dots$ 按原文“指代 qubit 状态”（在该推导中用于得到不同的界）。
- 当 $j=0$ 时，恢复 JC 模型常见的结果：
  $$
  n_{\rm crit}=\left(\frac{\Delta}{2g}\right)^2.
  $$
- 文中说 $j=1$ 会给出更保守的估计。
- 同时强调：该条件只是粗略估计，用于判断何时高阶效应开始重要。

---

### **C6. 与两能级 JC 色散结果的对比：$\chi=g^2/\Delta$ 何时成立？**
- 对 JC 两能级模型做色散近似会得到：
  $$
  \chi=\frac{g^2}{\Delta}.
  $$
- 文中指出：
  - 在某些极限下该结果与上面对 transmon 的结果一致；
  - 并且在很多 transmon 实验里，由于 $E_C/h$ 往往相对 $\Delta$ 较小，JC 两能级的 $\chi=g^2/\Delta$ 在实践中常能给出不错预测。
- 关键物理直觉：  
  **$E_C$** 决定 transmon 的非谐性。若 $E_C\to 0$，系统趋向两个耦合谐振子，只能产生**与状态无关**的频移，因此 **色散的“态依赖频移”必须消失**。

---

### **C7. Bogoliubov 方法（文中 2. Bogoliubov approach）：更“省事”的推导路线（Eq.45–51）**
文中给出一种在电路 QED 文献中常用、比 SW 更直接的推导方式。

#### **C7.1 把 Eq.(32) 拆成线性 + 非线性（Eq.45–46）**
把哈密顿量写成 $\hat H=\hat H_L+\hat H_{NL}$：
$$
\hat H_L=\hbar\omega_r\hat a^\dagger\hat a+\hbar\omega_q\hat b^\dagger\hat b
+\hbar g(\hat b^\dagger\hat a+\hat b\hat a^\dagger),
\tag{45}
$$
$$
\hat H_{NL}=-\frac{E_C}{2}\hat b^\dagger\hat b^\dagger\hat b\hat b.
\tag{46}
$$

#### **C7.2 用 Bogoliubov 变换精确对角化线性部分（Eq.47–49）**
幺正变换：
$$
\hat U_{\rm disp}=\exp\!\left[\Lambda(\hat a^\dagger\hat b-\hat a\hat b^\dagger)\right],
\tag{47}
$$
在该变换下（原文给出）湮灭算符混合：
- $\hat a$ 与 $\hat b$ 发生 $\cos\Lambda,\sin\Lambda$ 的线性旋转。
- 取
  $$
  \Lambda=\frac12\arctan(2\lambda),\qquad \lambda=\frac{g}{\Delta},
  $$
可使 $\hat H_L$ 完全对角化：
$$
\hat U_{\rm disp}^\dagger \hat H_L \hat U_{\rm disp}
=\hbar\bar\omega_r\hat a^\dagger\hat a+\hbar\bar\omega_q\hat b^\dagger\hat b.
\tag{48}
$$
对应 dressed（混合后）线性本征频率：
$$
\bar\omega_r=\frac12\left(\omega_r+\omega_q-\sqrt{\Delta^2+4g^2}\right),
\tag{49a}
$$
$$
\bar\omega_q=\frac12\left(\omega_r+\omega_q+\sqrt{\Delta^2+4g^2}\right).
\tag{49b}
$$

#### **C7.3 作用到非线性项并在色散区展开 → Kerr 与 cross-Kerr（Eq.50–51）**
对 $\hat H_{NL}$ 做同样变换，在色散区按 $\Lambda$ 展开（细节见 Appendix B.3），得到：
$$
\hat H_{\rm disp}\approx
\hbar\bar\omega_r\hat a^\dagger\hat a+\hbar\bar\omega_q\hat b^\dagger\hat b
+\frac{\hbar K_a}{2}\hat a^\dagger\hat a^\dagger\hat a\hat a
+\frac{\hbar K_b}{2}\hat b^\dagger\hat b^\dagger\hat b\hat b
+\hbar\chi_{ab}\hat a^\dagger\hat a\,\hat b^\dagger\hat b.
\tag{50}
$$
并给出近似表达式：
$$
K_a\approx -\frac{E_C}{2\hbar}\left(\frac{g}{\Delta}\right)^4,\qquad
K_b\approx -\frac{E_C}{\hbar},
\tag{51a}
$$
$$
\chi_{ab}\approx -2\frac{g^2(E_C/\hbar)}{\Delta(\Delta-E_C/\hbar)}.
\tag{51b}
$$

解释与备注（按原文要点）：
- $K_a$：resonator self-Kerr；$K_b$：transmon self-Kerr；$\chi_{ab}$：cross-Kerr。
- 在通常色散区这三者为负。
- $\chi_{ab}$ 的表达式来自进一步 SW 消去 $\propto \hat b^\dagger\hat b^\dagger\hat a\hat b$ 的项；以及把随频率 $\sim\Delta$ 或更快振荡的高阶项丢弃（文中指出这些项写在 Eq.(B32)）。
- 把 Eq.(50) 截断到 transmon 两能级，能正确回到 Eq.(42) 与 Eq.(43)。

#### **C7.4 这些色散公式何时失效？“straddling regime”**
文中强调 Eq.(50)/(51) **不适用**的典型情形：
- resonator 或 transmon 的激发数过大；
- 或 $|\Delta|\sim E_C/\hbar$ 等导致近似不可靠。
- 特别地，若
  $$
  0<\Delta<E_C
  $$
  称为 **straddling regime**，与通常色散区性质不同：
  - $K_a$、$\chi_{ab}$ 可为正；
  - 更适合用 Eq.(29) 的**精确数值对角化**处理（Koch et al., 2007）。

#### **C7.5 强色散与大非线性：电路 QED 的“优势区”**
文中给出的结论性描述：
- 电路 QED 在色散区可实现非常大的非线性与频移：
  - $\chi$ 可大于腔或 qubit 的线宽（$\chi>\kappa,\gamma$）→ **strong dispersive coupling**（Gambetta 2006；Schuster 2007）。
  - 甚至可以实现 resonator 的大 self-Kerr：$K_a>\kappa$。
- 非线性可通过在谐振器中心导体嵌入 Josephson 结增强（Bourassa 2012；Ong 2013），用于量子极限参量放大器、制备微波场量子态等应用。

---

## **III.D Josephson junctions embedded in multimode electromagnetic environments**  
（多模电磁环境中的 Josephson 结）

> 目标：从“单模耦合”推广到更一般的情形：一个（或多个）结处在复杂的线性电磁环境中（如 3D 腔），需要系统地得到多模频移与 Kerr / cross-Kerr。

### **D1. 问题设置：结 + 任意线性环境（Fig.10a）**
- 之前主要讨论 transmon 与单模谐振子耦合。
- 许多实验中，需要考虑 transmon（甚至多个 transmon）嵌入复杂几何的电磁环境（例如 3D cavity），可用一个等效的 **阻抗 $Z(\omega)$** 来表征环境（Fig.10a）。

### **D2. 结的线性化分解（Fig.10b）与总哈密顿量拆分**
- 一个电容分流（capacitively shunted）的 Josephson 结可以分解为：
  - 线性电容（结电容 + 分流电容的并联）：$C_\Sigma=C_S+C_J$
  - 线性电感：  
    $$
    L_J=E_J^{-1}\left(\frac{\Phi_0}{2\pi}\right)^2
    $$
  - 以及一个“纯非线性元件”（Fig.10(b) 的蜘蛛符号）
- 假设外部电磁环境是 **线性、非磁性、无自由电荷与电流**。
- 将线性部分（含 $C_\Sigma,L_J$ 与环境）合并为 $\hat H_L$，剩余非线性为 $\hat H_{NL}$。非线性项写成：
$$
\hat H_{NL}=-E_J\left(\cos\hat\phi+\frac12\hat\phi^2\right),
\tag{52}
$$
其中 $\frac12\hat\phi^2$ 表示把二次项从 $\cos\hat\phi$ 中剥离出去并已归入线性部分（与前文对 $\cos$ 展开/线性化的思想一致）。

### **D3. 模展开量子化：从电磁场本征模得到 $\hat H_L$（Eq.53–54）**
选用电磁环境的正则场变量（文中选取）：
- 电位移场 $\hat{\mathbf D}(\mathbf x)$ 与磁场 $\hat{\mathbf B}(\mathbf x)$：
$$
\hat{\mathbf D}(\mathbf x)=\sum_m\left[\mathbf D_m(\mathbf x)\hat a_m+H.c.\right],
\tag{53a}
$$
$$
\hat{\mathbf B}(\mathbf x)=\sum_m\left[\mathbf B_m(\mathbf x)\hat a_m+H.c.\right],
\tag{53b}
$$
且 $[\hat a_m,\hat a_n^\dagger]=\delta_{mn}$。
- 通过选择满足正交与归一化条件的模函数，可使线性哈密顿量成为：
$$
\hat H_L=\sum_m \hbar\omega_m \hat a_m^\dagger \hat a_m.
\tag{54}
$$
- 默认假设本征模是离散谱；若存在连续谱（如开放波导），求和需替换为对频率范围的积分。

> 原文强调：求模函数 $\{\mathbf D_m,\mathbf B_m\}$ 本质是**经典电磁本征模问题**，可用有限元等数值软件求解（Minev 2019）。文中后续将只讨论离散谱。

### **D4. 结相位算符的模分解（Eq.55）**
线性部分对角化后，结两端相位差 $\hat\phi$ 可写成线性叠加：
$$
\hat\phi=\sum_m\left[\varphi_m \hat a_m+H.c.\right],
\tag{55}
$$
- $\varphi_m$：第 $m$ 个模在结处“看到的”**零点涨落幅度**（dimensionless magnitude of zero-point fluctuations），由模函数与边界条件决定（Fig.10(a) 给出结的位置与积分路径标记）。

### **D5. 多模色散区的 Kerr / cross-Kerr 结构（Eq.56–58）**
- 将 Eq.(55) 代入 Eq.(52)，在 transmon 区（小非谐性）把 $\cos$ 展开到四阶，类比前文做法。
- 在“各本征频率充分分离”的色散条件下，忽略快旋项（类比 Sec. III.C），得到：
$$
\hat H_{NL}\approx
\sum_m \hbar\Delta_m \hat a_m^\dagger\hat a_m
+\frac12\sum_m \hbar K_m(\hat a_m^\dagger\hat a_m)^2
+\sum_{m>n}\hbar\chi_{mn}\hat a_m^\dagger\hat a_m\,\hat a_n^\dagger\hat a_n,
\tag{56}
$$
并定义（原文紧随其后）：
- $\Delta_m=\frac12\sum_n \chi_{mn}$（频率位移）
- $K_m=\chi_{mm}$（self-Kerr）
- 以及
$$
\hbar\chi_{mn}=-E_J\varphi_m^2\varphi_n^2.
\tag{57}
$$

进一步引入 **能量参与比（energy participation ratio）** $p_m$：第 $m$ 模的总电感能中存储在结中的比例。定义满足
$$
p_m=\left(\frac{2E_J}{\hbar\omega_m}\right)\varphi_m^2 \quad \text{(Minev 2019)},
$$
从而得到更“工程化”的 Kerr/cross-Kerr：
$$
\chi_{mn}=-\frac{\hbar\omega_m\omega_n}{4E_J}p_m p_n.
\tag{58}
$$

> 关键结论（原文总结语气）：求非线性哈密顿量可归结为求系统本征模与结处零点涨落（或 $p_m$）。这在复杂几何中可能麻烦，但本质仍是经典电磁学问题。

### **D6. 用阻抗 $Z(\omega)$ 表示线性环境（Fig.10c）：等效 LC 串联模展开**
- 文中给出替代路线：把纯非线性元件“看到的”线性环境用一个阻抗 $Z(\omega)$ 表示（Fig.10c）。
- 忽略损耗时，任意这样的阻抗可等效为**（可能无限多个）串联的 LC 振子链**。
- 本征频率 $\omega_m$ 可由导纳
  $$
  Y(\omega)=Z^{-1}(\omega)
  $$
  的零点（更准确地说与其零/极结构相联系）得到；每个模的有效阻抗也可由 $Y(\omega)$ 在零点附近的性质确定。
- 文中还给出与零点涨落相关的“有效阻抗”表达式（你截图里出现的等式）：
  $$
  Z_m^{\rm eff}=2\left(\frac{\Phi_0}{2\pi}\right)^2\frac{\varphi_m^2}{\hbar}
  =\frac{R_K}{4\pi}\varphi_m^2.
  $$
- 结论：量子化步骤在这种表述下，归结为确定 $Z(\omega)$ 随频率的函数形式。

---

## **III.E Beyond the transmon: multilevel artificial atom**  
（超越 transmon 近似：一般多能级人工原子）

> 动机：前文多次提到两点——(i) transmon 实际是多能级；(ii) $\cos$ 势能的微扰展开依赖 $E_J/E_C\gg 1$。本节给出不依赖“近谐振子近似”的更一般写法。

### **E1. 从“$\cos$”势的微扰展开转向精确对角化**
- 之前在 transmon 区（$E_J/E_C\gg 1$）用 $\cos$ 势的微扰展开得到 Duffing 模型。
- 若要超越该近似，需要对 transmon 哈密顿量做**精确数值对角化**。

### **E2. 把 full transmon–resonator 哈密顿量写成“多能级原子 + 腔模 + 一般耦合”（Eq.59–60）**
回到完整的 Eq.(29)（Koch et al., 2007），可写为（文中写法）：
$$
\hat H
=4E_C\hat n^2-E_J\cos\hat\phi+\hbar\omega_r\hat a^\dagger\hat a+8E_C\hat n\hat n_r
$$
并在 bare transmon 本征态 $|j\rangle$（由
$\hat H_T=4E_C\hat n^2-E_J\cos\hat\phi$ 数值对角化得到）基底下写成：
$$
\hat H
=\sum_j \hbar\omega_j |j\rangle\langle j|
+\hbar\omega_r\hat a^\dagger\hat a
+i\sum_{ij}\hbar g_{ij}|i\rangle\langle j|(\hat a^\dagger-\hat a),
\tag{59}
$$
其中耦合矩阵元定义为：
$$
\hbar g_{ij}=2e\frac{C_g}{C_\Sigma}Q_{\rm zpf}\langle i|\hat n|j\rangle.
\tag{60}
$$
- $\omega_j$、$\langle i|\hat n|j\rangle$：
  - 可在电荷基底中数值计算；
  - 也可利用相位表象下出现的 Mathieu 方程结构（Eq.(20) 的形式）得到（文中提及：Mathieu 方程有已知精确解，见 Cotlet 2002；Koch 2007）。

> 原文强调：Eq.(59) 这种“按能级写出”的形式非常一般，可描述**任意多能级人工原子**与腔模的电容耦合。

### **E3. 一般多能级系统的色散近似（Eq.61–62）**
在色散区，若对所有相关跃迁 $i\neq j$ 满足 $|g_{ij}|\sqrt{n+1}\ll |\omega_{ij}-\omega_r|$（文中口头描述：对所有相关原子跃迁与光子数都应远失谐），可做 SW 变换二阶近似对角化（Zhu et al., 2013）得到：
$$
\hat H\approx
\sum_j \hbar(\omega_j+\Lambda_j)|j\rangle\langle j|
+\hbar\omega_r\hat a^\dagger\hat a
+\sum_j \hbar\chi_j \hat a^\dagger\hat a\,|j\rangle\langle j|.
\tag{61}
$$
其中
$$
\Lambda_j=\sum_i\frac{|g_{ij}|^2}{\omega_j-\omega_i-\omega_r},
\tag{62a}
$$
$$
\chi_j=\sum_i |g_{ij}|^2\left(
\frac{1}{\omega_j-\omega_i-\omega_r}
-\frac{1}{\omega_i-\omega_j-\omega_r}
\right).
\tag{62b}
$$
- $\Lambda_j$：多能级的一般 Lamb shift。
- $\chi_j$：处于 $|j\rangle$ 时腔频的状态依赖位移（读出所用的“cavity pull”在多能级下的推广）。

原文补充：该结果相当一般，可用于多种人工原子在色散极限下与谐振器耦合的情形；更高阶表达式可见 Boissonneault 2010；Zhu 2013。

---

## **III.F Alternative coupling schemes**  
（其它耦合方案）

> 本节强调：电路 QED 中“光—物质耦合”不只限于电偶极（电容）耦合，还可以通过磁耦合、以及纵向耦合等实现新的物理与应用。

### **F1. 磁耦合：利用互感将 flux qubit 与谐振器耦合**
- 常见耦合：qubit 电偶极矩与腔的零点电场耦合（前文主线）。
- 另一类：利用 **flux qubit** 与谐振器中心导体之间的**互感**，使 qubit 的磁偶极与谐振器零点磁场耦合。
- 更强耦合可通过把 flux qubit **galvanically**（直流电连接）到传输线谐振器中心导体来实现（Bourassa et al., 2009）。
- 这种情况下耦合强度可以做到与系统频率同量级甚至更大，从而进入 **ultrastrong coupling regime**（文中指向 Sec. VI.C）。

### **F2. 纵向（longitudinal）耦合：只移频、不翻转（Eq.63）**
另一种思路：设计耦合使得腔场不引起 qubit 的跃迁，只让 qubit 频率随腔场位移而变——即纵向耦合（Billangeon 2015；Didier 2015；Kerman 2013；Richer & DiVincenzo 2016 等）。哈密顿量写为：
$$
\hat H_z=
\hbar\omega_r\hat a^\dagger\hat a
+\frac{\hbar\omega_q}{2}\hat\sigma_z
+\hbar g_z(\hat a^\dagger+\hat a)\hat\sigma_z.
\tag{63}
$$
- 与 JC（横向耦合，$\propto \sigma_x$）不同，这里相互作用 $\propto \sigma_z$。
- 原文指出的重要差异：由于耦合与 $\sigma_z$ 成正比，纵向相互作用**不会像 Sec. III.B 那样导致腔场对 qubit 的 dressing**（以及相应的 JC 混合结构）。
- 与读出相关的后果（尤其是 qubit readout）在 Sec. V.C.3 讨论。

---

### **整体逻辑总结**

**III.C 部分** 系统阐述了电路 QED 中最重要的色散区物理。首先给出色散条件，然后通过两种方法（Schrieffer–Wolff 和 Bogoliubov 变换）推导出色散哈密顿量，得到重整化的频率、态依赖的腔频位移 $\chi$、以及高阶 Kerr 非线性。特别强调了临界光子数 $n_{\rm crit}$ 作为近似有效性的判据，并指出 straddling regime 等失效情形。最后展示了电路 QED 在色散区可实现强非线性和强色散耦合的优势。

**III.D 部分** 将分析推广到多模电磁环境。通过将 Josephson 结的线性部分与外部环境合并，将非线性部分展开，并引入能量参与比 $p_m$，得到了多模 Kerr 和 cross-Kerr 的简洁表达式。这为处理复杂几何（如 3D 腔）中的多模效应提供了系统方法。

**III.E 部分** 超越 transmon 的 Duffing 近似，回到更一般的多能级人工原子描述。通过数值对角化 transmon 哈密顿量得到本征态和耦合矩阵元，将系统哈密顿量写成一般形式，并给出了多能级色散近似的 Lamb shift 和 cavity pull 表达式，为处理任意多能级系统提供了统一框架。

**III.F 部分** 简要介绍了除电容耦合以外的两种耦合方案：磁耦合（可实现超强耦合）和纵向耦合（只移频不翻转），为后续讨论更丰富的物理和应用埋下伏笔。

整体而言，第三章从最简单的交换相互作用开始，逐步深入到色散区、多模环境、多能级系统以及替代耦合方案，构建了电路 QED 中光-物质相互作用的完整理论体系。


# IV. 与外界耦合：环境的作用（Coupling to the outside world: the role of the environment）

  

到目前为止讨论的是**孤立量子系统**。但对量子电路的完整描述必须包含系统如何与环境耦合（测量装置、控制线路等）。环境在量子技术中具有双重角色：  

- **不可避免性**：完全隔离不现实，因为总会存在非期望的环境自由度耦合。  

- **必要性**：完全隔离也“没用”，因为我们将无法控制、驱动或测量系统。  

  

因此本节考虑：量子系统与**外部传输线**耦合，并引入对电路QED读出至关重要的**输入-输出理论（input-output theory）**。

  

---

  

## A. 用传输线连接量子系统（Wiring up quantum systems with transmission lines）

  

### A.1 传输线作为环境：连续谱极限

考虑单个量子系统（示例：LC振子）通过电容耦合到一根**半无限**传输线（如共面波导，CPW）。  

- 传输线长度 $d\to \infty$ 时，模频率变得**极其稠密**，需视为**连续谱（continuum）**。

- 传输线哈密顿量写为连续模积分形式：

  

$$

\hat H_{\mathrm{tml}}=\int_0^\infty d\omega\,\hbar\omega\,\hat b_\omega^\dagger \hat b_\omega .

\tag{64}

$$

  

其中 $\hat b_\omega$ 为频率为 $\omega$ 的传输线模湮灭算符，并满足连续模对易关系（文中随后用到）

$$

[\hat b_\omega,\hat b_{\omega'}^\dagger]=\delta(\omega-\omega') .

$$

  

---

  

### A.2 传输线的场变量（通量与电荷）的连续模展开

传输线在位置 $x$ 处的通量与电荷算符（与离散模展开类似，在连续极限下）为：

  

$$

\hat\Phi_{\mathrm{tml}}(x)=\int_0^\infty d\omega\,

\sqrt{\frac{\hbar}{\pi\omega c v}}\,

\cos\!\left(\frac{\omega x}{v}\right)\,(\hat b_\omega^\dagger+\hat b_\omega),

\tag{65a}

$$

  

$$

\hat Q_{\mathrm{tml}}(x)= i\int_0^\infty d\omega\,

\sqrt{\frac{\hbar\omega c}{\pi v}}\,

\cos\!\left(\frac{\omega x}{v}\right)\,(\hat b_\omega^\dagger-\hat b_\omega).

\tag{65b}

$$

  

- 这里 $v=1/\sqrt{lc}$ 为传输线相速度；$l,c$ 分别是单位长度电感、电容（注意与后文的耦合电容符号区分）。  

- 在海森堡绘景下，$\hat\Phi_{\mathrm{tml}}(x,t)=e^{i\hat H_{\mathrm{tml}}t/\hbar}\hat\Phi_{\mathrm{tml}}(x)e^{-i\hat H_{\mathrm{tml}}t/\hbar}$（$\hat Q$ 同理）。

  

---

  

### A.3 通过电容将传输线耦合到谐振子（LC）

考虑 LC 振子（系统）在 $x=0$ 端与传输线通过耦合电容 $C_k$ 相连（图11）。系统哈密顿量为

$$

\hat H_S=\hbar\omega_r \hat a^\dagger \hat a ,

$$

其中 $\hat a$ 为谐振子湮灭算符，$\omega_r$ 为共振频率。

  

在小耦合近似下，耦合后的总哈密顿量为（系统 + 传输线 + 相互作用）：

  

$$

\hat H=\hat H_S+\hat H_{\mathrm{tml}}

-\hbar\int_0^\infty d\omega\,\lambda(\omega)\,(\hat b_\omega^\dagger-\hat b_\omega)(\hat a^\dagger-\hat a).

\tag{66}

$$

  

耦合强度的频率依赖形式为

$$

\lambda(\omega)=\left(\frac{C_k}{\sqrt{C_r C_t}}\right)\sqrt{\frac{\omega_r}{2\pi}\,\omega},

$$

其中：

- $C_r$ 为谐振子电容；

- $C_t$ 为传输线等效电容（文中以“归一化/有效电容”方式出现）；

- 这些表达式可视作因耦合导致的电容小重整化（附录A讨论）。

  

---

  

### A.4 窄带/高Q与旋转波近似：得到标准系统-浴耦合形式

假设 $\lambda(\omega)$ 在 $\omega_r$ 附近变化缓慢且相对小，可当作微扰；系统 $Q$ 因子大，主要只响应 $\omega\approx\omega_r$ 的窄带频率。于是取近似 $\lambda(\omega)\simeq \lambda(\omega_r)$ 并丢弃快速振荡项，得到（Gardiner & Zoller 1999）：

  

$$

\hat H \approx \hat H_S+\hat H_{\mathrm{tml}}

+\hbar\int_0^\infty d\omega\,\lambda(\omega_r)\,(\hat a\hat b_\omega^\dagger-\hat a^\dagger\hat b_\omega).

\tag{67}

$$

  

这就是典型的“能量守恒型”（RWA）耦合：系统放光到浴 / 从浴吸收。

  

---

  

### A.5 Born–Markov 近似：Lindblad 主方程（谐振子衰减与热激发）

在标准 Born–Markov 近似下，由式(67)可得系统密度矩阵 $\rho$ 的 Markovian 主方程（Lindblad 形式）：

  

$$

\dot\rho=-i[\hat H_S,\rho]+\kappa(\bar n_\kappa+1)\mathcal D[\hat a]\rho

+\kappa\bar n_\kappa\,\mathcal D[\hat a^\dagger]\rho.

\tag{68}

$$

  

- $\kappa$：谐振子的能量衰减率（线宽）。文中给出

$$

\kappa=2\pi\lambda(\omega_r)^2

= Z_{\mathrm{tml}}\omega_r^2\frac{C_k^2}{C_r},

$$

其中 $Z_{\mathrm{tml}}$ 为传输线特性阻抗。  

- $\bar n_\kappa=\bar n_\kappa(\omega_r)$：传输线在频率 $\omega_r$ 处的热光子数（Bose–Einstein 分布给出），并满足浴相关函数

$$

\langle \hat b_\omega^\dagger \hat b_{\omega'}\rangle=\bar n_\kappa(\omega)\delta(\omega-\omega').

$$

  

耗散超算符定义为

  

$$

\mathcal D[\hat O]\rho=\hat O\rho\hat O^\dagger-\frac12\{\hat O^\dagger\hat O,\rho\},

\tag{69}

$$

  

其中 $\{\cdot,\cdot\}$ 为反对易子。

  

**物理解释（对应式(68)两项耗散）：**

- $\mathcal D[\hat a]$ 项：对应**光子损失**（衰减）。当作用在数态 $\lvert n\rangle$ 上，$\hat a\lvert n\rangle=\sqrt n\lvert n-1\rangle$，体现 $n\to n-1$ 的跃迁。

- $\mathcal D[\hat a^\dagger]$ 项：对应从热浴**吸收光子**导致的激发。

- $T\to 0$（稀释制冷极低温）常取 $\bar n_\kappa\to 0$，但实际可能因室温线缆热辐射等导致残余热占据，需要工程抑制（吸收器、滤波、隔离等）。

  

---

  

## B. 电网络中的输入-输出理论（Input-output theory in electrical networks）

  

### B.1 为什么需要输入-输出理论

主方程（如式(68)）描述系统的阻尼动力学，但**不直接给出系统向外辐射的场**。实验中可测量的正是这些辐射/反射信号，因此需要输入-输出理论来连接：

- 传输线上的输入场（噪声/驱动）

- 系统动力学

- 输出场（可测读出的信号）

  

本节介绍一种常用推导路线：将传输线模分解为**左行（输入）**与**右行（输出）**辐射场，并用 $x=0$ 处边界条件联系二者。

  

---

  

### B.2 左/右行场与端口电压算符

在 $x=0$（系统端口位置）处，将电压分为输入与输出贡献：

$$

\hat V(x=0,t)=\hat V_{\mathrm{in}}(t)+\hat V_{\mathrm{out}}(t).

$$

  

文中给出端口输入/输出电压的模展开（用左右行模算符 $\hat b_{L/R,\omega}$ 表示）：

  

$$

\hat V_{\mathrm{in/out}}(t)=

i\int_0^\infty d\omega\,

\sqrt{\frac{\hbar\omega}{4\pi c v}}\,

e^{-i\omega t}\,\hat b_{L/R,\omega}+ \mathrm{H.c.}

\tag{70}

$$

  

---

  

### B.3 边界条件（Kirchhoff 定律）与端口电流

端口处满足 Kirchhoff 电流定律给出的边界条件：

  

$$

\hat I(t)=\frac{\hat V_{\mathrm{out}}(t)-\hat V_{\mathrm{in}}(t)}{Z_{\mathrm{tml}}},

\tag{71}

$$

  

其中左侧 $\hat I(t)$ 是由样品注入传输线的电流。文中说明该电流可写为

$$

\hat I(t)=\left(\frac{C_r}{C_k}\right)\dot{\hat Q}_r(t),

$$

$\hat Q_r$ 为谐振子电荷（细节见附录C的推导）。

  

---

  

### B.4 标准输入-输出关系与输入/输出场定义

通过对相关算符做模展开并采用前述近似（系统仅响应 $\omega\approx\omega_r$，忽略快速旋转项，可将积分范围扩展到 $(-\infty,\infty)$），得到**标准输入-输出关系**：

  

$$

\hat b_{\mathrm{out}}(t)-\hat b_{\mathrm{in}}(t)=\sqrt{\kappa}\,\hat a(t).

\tag{72}

$$

  

输入/输出场算符定义为（以频率相对系统频率的失谐表示）：

  

$$

\hat b_{\mathrm{in}}(t)=

\frac{-i}{\sqrt{2\pi}}\int_{-\infty}^{\infty} d\omega\,\hat b_{L,\omega}\,e^{-i(\omega-\omega_r)t},

\tag{73}

$$

  

$$

\hat b_{\mathrm{out}}(t)=

\frac{-i}{\sqrt{2\pi}}\int_{-\infty}^{\infty} d\omega\,\hat b_{R,\omega}\,e^{-i(\omega-\omega_r)t},

\tag{74}

$$

  

并满足白噪声型对易关系

$$

[\hat b_{\mathrm{in}}(t),\hat b_{\mathrm{in}}^\dagger(t')]=

[\hat b_{\mathrm{out}}(t),\hat b_{\mathrm{out}}^\dagger(t')]=\delta(t-t').

$$

  

> 关键近似：丢弃 $\pm\omega_r$ 的快速旋转项（RWA/窄带），并假设相关频率都靠近 $\omega_r$，从而得到 Markov 极限的 $\delta(t-t')$ 关联。

  

---

  

### B.5 海森堡-朗之万方程：谐振子由输入场驱动并以 $\kappa$ 衰减

使用同样的近似，可得到谐振子场算符的运动方程（海森堡绘景）：

  

$$

\dot{\hat a}(t)=i[\hat H_S,\hat a(t)]-\frac{\kappa}{2}\hat a(t)+\sqrt{\kappa}\,\hat b_{\mathrm{in}}(t).

\tag{75}

$$

  

含义：

- $-\kappa\hat a/2$：辐射阻尼（能量泄露到线中）

- $\sqrt{\kappa}\hat b_{\mathrm{in}}$：输入场（噪声/外加驱动）对系统的驱动

- 结合(72)可由输入与系统动力学预测输出，从而与实验可测信号建立联系。

  

---

  

### B.6 传输线中某位置的电压（可测信号表达式）

文中指出：可在 $x>0$ 处测量电压以探测输出场。在上述近似下，电压算符可写为

  

$$

\hat V(x,t)\simeq

\sqrt{\frac{\hbar\omega_r Z_{\mathrm{tml}}}{2}}

\Big[

e^{i\omega_r x/v-i\omega_r t}\hat b_{\mathrm{out}}(t)

+

e^{-i\omega_r x/v-i\omega_r t}\hat b_{\mathrm{in}}(t)

+\mathrm{H.c.}

\Big].

\tag{76}

$$

  

并强调：该表达式假设相关频率都在 $\omega_r$ 附近，同时忽略非Markov时间延迟效应。

  

---

  

### B.7 扩展性：从单端口到量子网络

本节主要讨论单个系统连接半无限传输线端点的简单情形。更一般地：

- 多个系统可耦合到同一条线形成网络；

- 多条线可作为多个端口。

这些更复杂结构可用 **SLH 形式化**系统化处理（文中提及 Combes、Gough & James 等）。

  

---

  

## C. 量子比特的弛豫与退相干（Qubit relaxation and dephasing）

  

### C.1 从谐振子到一般量子系统：主方程的普适性

式(68)虽由“谐振子-传输线”模型推得，但式(66)本质是非常通用的系统-浴哈密顿量，可用来描述多种噪声源导致的耗散（如 Caldeira–Leggett 思路）。  

对任意量子系统，只需将式(66)中的 $\hat a$ 替换为与浴耦合的系统算符即可。

  

---

  

### C.2 Transmon 的能量弛豫（连接传输线导致的衰减）

对 transmon（图12）：  

- 将式(68)中 $\hat H_S$ 替换为 transmon 哈密顿量 $\hat H_q$（文中指向前面 Eq.(25)）。  

- 并将耗散算符替换为 transmon 的跃迁算符：

  - $\mathcal D[\hat a]\to \mathcal D[\hat b]$

  - $\mathcal D[\hat a^\dagger]\to \mathcal D[\hat b^\dagger]$

- 同时将 $\kappa\to\gamma$，其中

$$

\gamma=2\pi\lambda(\omega_q)^2

$$

为在量子比特频率 $\omega_q$ 处评估的弛豫率（与耦合强度相关）。

  

得到量子比特（与环境交换能量）主方程：

  

$$

\dot\rho=-i[\hat H_q,\rho]+\gamma(\bar n_\gamma+1)\mathcal D[\hat b]\rho

+\gamma\bar n_\gamma\,\mathcal D[\hat b^\dagger]\rho,

\tag{77}

$$

  

其中：

- $\rho$ 现在是 transmon 的密度矩阵；

- $\bar n_\gamma$ 是量子比特环境在 $\omega_q$ 处的热占据数。常取 $\bar n_\gamma\to 0$，但实验上也常观察到残余热激发。

  

---

  

### C.3 纯退相干（pure dephasing）的表型模型

超导量子电路除了能量弛豫，还会因参数涨落（如通量噪声、色散耦合导致的频率抖动等）产生退相干。对 transmon，可在主方程中加入表型纯退相干项：

  

$$

2\gamma_\phi\,\mathcal D[\hat b^\dagger \hat b]\rho,

\tag{78}

$$

  

其中 $\gamma_\phi$ 为纯退相干率。文中讨论要点：

- 对 0–1 跃迁的 transmon，因其对电荷噪声不敏感，$\gamma_\phi$ 往往较小；但高能级因电荷色散随能级数增长而变得更明显，可能出现更强退相干。

- 另一个退相干来源：与读出谐振子色散耦合时，谐振子中残余热光子数的涨落会引起量子比特频率抖动（“光子数涨落退相干”）。在这种情况下，也可在谐振子主方程中加入类似式(78)但以 $\hat a^\dagger \hat a$ 替代的项；通常该贡献较小而常被忽略。

- 还存在材料/界面双能级系统（TLS）、准粒子、红外与电离辐射等引起的弛豫与退相干来源。

  

---

  

### C.4 合并后（取 $\bar n_\gamma\to 0$）的 transmon 主方程

将能量弛豫（零温极限）与纯退相干合并：

  

$$

\dot\rho=-i[\hat H_q,\rho]+\gamma\mathcal D[\hat b]\rho

+2\gamma_\phi\,\mathcal D[\hat b^\dagger \hat b]\rho.

\tag{79}

$$

  

常将其化为两能级近似形式：令

$$

\hat H_q \to \hbar\omega_q\frac{\hat\sigma_z}{2},\quad

\hat b^\dagger \to \hat\sigma_+,\quad

\hat b \to \hat\sigma_-.

$$

  

---

  

### C.5 $T_1$ 与 $T_2$：弛豫与相干时间

文中给出与实验常用指标的关系：

- 纵向弛豫时间（能量弛豫）：

$$

T_1=\frac{1}{\gamma_1}=\frac{1}{\gamma_\downarrow+\gamma_\uparrow}\simeq \frac{1}{\gamma},

\tag{80a}

$$

其中在有限温度时 $\gamma_\downarrow=(\bar n_\gamma+1)\gamma$、$\gamma_\uparrow=\bar n_\gamma\gamma$；零温近似下 $T_1\simeq 1/\gamma$。

- 横向退相干时间（相干寿命）：

$$

T_2=\frac{1}{\gamma_2}=\left(\frac{\gamma_1}{2}+\gamma_\phi\right)^{-1}.

\tag{80b}

$$

  

并进一步指出实验典型量级与材料体系相关（文中列举了多种器件的大致范围），同时强调：

- **能量弛豫**过程往往可较好地用 Markov 主方程描述；

- **退相干**（尤其由低频噪声导致）常具有显著非Markov特征，Markov 处理对退相干只是某种“粗略近似”，但在预测稳态响应等场景仍常有效。

  

---

  
  

## D. 色散区的耗散（Dissipation in the dispersive regime）

  

本节讨论 **transmon–谐振子** 在色散耦合区（dispersive regime）下的耗散描述。关键点：在色散区里，系统的本征态是**穿衣态（dressed states）**，比特与腔并非真正独立；因此即使耗散通道看似“分别作用在比特/腔”，穿衣结构也会导致**交叉耗散效应**（例如 Purcell 衰减、dressed dephasing）。

  

### D.1 “独立浴”近似下的复合系统主方程（bare basis 的直觉式）

假设 transmon 与谐振子分别耦合到独立环境（图13示意），并取 $\bar n\to 0$（为简化），可写复合系统（比特+谐振子）的主方程为

  

$$

\dot\rho = -i[\hat H,\rho] + \kappa\mathcal D[\hat a]\rho + \gamma \mathcal D[\hat b]\rho + 2\gamma_\phi \mathcal D[\hat b^\dagger \hat b]\rho .

\tag{81}

$$

  

- $\rho$：**总系统**密度矩阵  

- $\hat H$：耦合的 transmon–谐振子哈密顿量（文中指向前面 Eq.(32)）  

- $\kappa$：腔衰减率（光子泄露）  

- $\gamma$：比特能量弛豫率  

- $\gamma_\phi$：比特纯退相干率  

- $\hat a$：谐振子湮灭算符；$\hat b$：transmon 的跃迁/降低算符（在弱非简谐下近似类振子）

  

**重要限制/提醒**：式(81)只对 $\omega_q/\omega_r$ 的小值区间有效。更本质的原因是：真实能量衰减发生在**系统本征态之间**的跃迁；而式(81)的耗散子更像是在未耦合的“裸态”（bare states）之间写跃迁，因而不是在任意耦合强度/任意频率比下都可靠。（文中提到一般推导可见 Beaudoin et al. 2011）

  

---

  

### D.2 穿衣态导致的交叉耗散：色散变换后的主方程

乍看式(81)似乎表示腔耗散与比特耗散彼此独立。但在色散区，由于穿衣态将 $|e,0\rangle$ 与 $|g,1\rangle$ 混合（图13），腔的阻尼可能诱导比特跃迁，反之亦然。

  

对系统做色散区的幺正变换（文中指向 Eq.(47) 的色散变换），不仅哈密顿量变换，主方程也会变换。忽略快振荡项并保留到 $\lambda$ 的二阶修正（其中 $\lambda\sim g/\Delta$，并与 $\kappa,\gamma,\gamma_\phi$ 同阶，即 $\kappa,\gamma,\gamma_\phi = O(E_C g^2/\Delta^2)$ 的量级一致），得到色散表象下主方程（Boissonneault et al. 2009）：

  

$$$$

\dot\rho_{\mathrm{disp}}

= -i[\hat H_{\mathrm{disp}},\rho_{\mathrm{disp}}]

+(\kappa+\kappa_\gamma)\mathcal D[\hat a]\rho_{\mathrm{disp}}

+(\gamma+\gamma_\kappa)\mathcal D[\hat b]\rho_{\mathrm{disp}}

+2\gamma_\phi \mathcal D[\hat b^\dagger\hat b]\rho_{\mathrm{disp}}

+\gamma_\Delta \mathcal D[\hat a^\dagger\hat a\,\hat b]\rho_{\mathrm{disp}}

+\gamma_\Delta \mathcal D[\hat a^\dagger\hat a\,\hat b^\dagger]\rho_{\mathrm{disp}} .

$$$$

$$

\tag{82}

$$

  

其中 $\rho_{\mathrm{disp}} = \hat U_{\mathrm{disp}}^\dagger\,\rho\,\hat U_{\mathrm{disp}}$。

  

#### D.2.1 新出现的三类修正速率（Purcell/反向Purcell/测量诱导跃迁）

式(82)中出现三种“交叉”速率（文中首先给出白噪声近似下的简单形式）：

  

$$

\kappa_\gamma=\left(\frac{g}{\Delta}\right)^2\kappa,\qquad

\gamma_\kappa=\left(\frac{g}{\Delta}\right)^2\gamma,\qquad

\gamma_\Delta=2\left(\frac{g}{\Delta}\right)^2\gamma_\phi .

\tag{83}

$$

  

并且文中强调式(82)“有三个新速率”，其中第一个是 **Purcell 衰减率**：

  

- **Purcell 衰减率 $\gamma_\kappa$**：腔的损耗 $\kappa$ 通过穿衣态混合导致比特弛豫。直觉图像：穿衣激发态近似为

  $$

  |e,0\rangle \ \to\ |e,0\rangle + \left(\frac{g}{\Delta}\right)|g,1\rangle,

  $$

  即主要是比特激发，但含有小概率 $\left(\frac{g}{\Delta}\right)^2$ 的“腔中有一个光子、比特在基态”的分量；该光子以速率 $\kappa$ 泄露，从而把系统带到 $|g,0\rangle$，表现为比特弛豫。

- **反向效应 $\kappa_\gamma$**：同理，比特的衰减通道也会导致腔的有效衰减修正（图13亦标注为 inverse Purcell）。

- **最后一行（含 $\gamma_\Delta$ 的两项）更微妙**：它们对应由比特退相干噪声在色散耦合下“上/下变频”，从而诱导比特在读出时发生非预期跃迁的过程（见下）。

  

---

  

### D.3 消去腔自由度得到“dressed dephasing”（读出引发的伪跃迁）

按 Boissonneault et al. (2008, 2009) 的处理，从式(82)可通过近似消去谐振子自由度得到只对 transmon 的有效主方程。其结果：比特的弛豫与激发速率会额外得到与腔内平均光子数 $\bar n$ 成正比的贡献，量级约为

  

- 弛豫/激发速率 $\sim \bar n\,\gamma_\Delta$

  

这在实验上常表现为量子比特测量期间的**非期望状态跃迁**，可解释为：退相干噪声在失谐频率 $\Delta$ 处被读出光子上/下变频，导致“dressed dephasing”并引发跃迁。

  

---

  

### D.4 更精细推导：耗散率的频率依赖与“噪声谱”视角

文中指出：上面对式(82)的“直接把色散变换作用到主方程”的捷径忽略了耗散率的频率依赖。更严格的做法是对**系统+浴的哈密顿量**先做色散变换，再推导主方程（Boissonneault et al. 2009）。这样得到的形式仍类似式(82)，但速率表达式会改变，并显示出“在不同频率探测环境”的物理：

  

- 记 $\kappa=\kappa(\omega_r)$，$\gamma=\gamma(\omega_q)$，强调：

  - 腔光子弛豫“探测”环境在 $\omega_r$ 处的噪声/阻抗

  - 比特弛豫“探测”环境在 $\omega_q$ 处的噪声/阻抗

  

更严格推导下，前两项修正速率变为

$$

\gamma_\kappa=\left(\frac{g}{\Delta}\right)^2\kappa(\omega_q),\qquad

\kappa_\gamma=\left(\frac{g}{\Delta}\right)^2\gamma(\omega_r).

$$

从而更清楚地看到：**Purcell 衰减是“在比特频率 $\omega_q$ 处发射一个光子”，而不是在腔频率 $\omega_r$ 处**（这也解释了白噪声近似可能造成的概念混淆）。

  

同样，对退相干率写成 $\gamma_\phi=\gamma_\phi(\omega\to 0)$，突出低频噪声的重要性。此时式(82)最后两项对应的速率可记为（文中给出）

$$

\gamma_{\Delta-}=2\left(\frac{g}{\Delta}\right)^2\gamma_\phi(-\Delta),\qquad

\gamma_{\Delta+}=2\left(\frac{g}{\Delta}\right)^2\gamma_\phi(\Delta),

$$

并总结为：

  

> **dressed dephasing 实际上在探测 transmon–谐振子失谐频率 $\Delta$ 处的噪声。**

  

文中提到该效应已被用于在 GHz 频段探测相关噪声（Slichter et al. 2012）。

  

---

  

### D.5 与纵向耦合（longitudinal coupling）的对比

本节结论依赖于 Jaynes–Cummings 型电偶极耦合导致的穿衣（dressing）。若将相互作用替换为纵向耦合（文中指向 Eq.(63) 形式），则不会产生同样的光-物质穿衣，因此：

  

- **无 Purcell 衰减**

- **无 dressed dephasing**

  

这也是纵向耦合方案的一个优势（Billangeon et al. 2015a; Kerman 2013）。

  

---

  

## E. 多模 Purcell 效应与 Purcell 滤波器（Multi-mode Purcell effect and Purcell filters）

  

### E.1 从单模到多模：Purcell 率的求和与发散问题

前面讨论的是比特色散耦合到**单模**谐振子。若谐振子是多模结构（多模腔/传输线谐振器），则“对各模的 dressing”会带来 Purcell 率的多模累加。朴素类比会得到修正速率形如

$$

\sum_m \left(\frac{g_m}{\Delta_m}\right)^2 \kappa_m,

$$

其中 $m$ 为模编号。

  

然而若进一步考虑 $f_m$、$g_m$、$\Delta_m$ 的频率依赖，该求和在某些模型下会出现**发散**（Houck et al. 2008）。因此需要更精细的模型（考虑 transmon 的有限尺寸、谐振子输入/输出电容的频率依赖阻抗等；Bourassa 2012; Malekakhlagh 2017）以得到物理上良好的结果。

  

---

  

### E.2 用“环境导纳”计算 Purcell：工程视角更直接

文中给出一个更简洁且工程上常用的表达：在量子电路中阻尼率由经典参数决定，可证明

  

$$

\gamma_\kappa=\mathrm{Re}[Y(\omega_q)]\,\frac{C_\Sigma}{C_s},

$$

  

其中

- $Y(\omega)=1/Z(\omega)$ 为 transmon “看到的”电磁环境的导纳

- $C_\Sigma$、$C_s$ 为与器件电容相关的参数（文中给出该形式以强调比例关系与“探测环境”的思想）

  

该表达再次强调：弛豫在系统频率处探测环境（导纳的实部决定耗散）。同时它也指向一个工程策略：

  

> 通过设计 $Y(\omega)$ 使其在 $\omega_q$ 附近尽可能“纯反应性”（$\mathrm{Re}\,Y(\omega_q)\approx 0$），可以抑制 Purcell 衰减。

  

---

  

### E.3 Purcell 滤波器（Purcell filter）的思想与实现

一种常见做法是在谐振子的输出端加入合适长度并在开路端终止的传输线支节（transmission-line stub），使得在比特频率处形成抑制辐射的响应。这种结构称为 **Purcell filter**（Reed et al. 2010b）。

  

其效果（与图13文字一致）：

- 在不显著牺牲读出需要的腔耗散/带宽的情况下，

- **降低比特频率处的“腔穿衣态密度”与环境耦合**，

- 从而抑制 Purcell 衰减。

  

并指出：由于可优化参数更多（一定程度上把“读出所需的 $\kappa$”与“比特弛豫”解耦），各种 Purcell filter 在实验中很常见（Bonn et al. 2015; Jeffrey et al. 2014; Walter et al. 2017）。

  

---

  

## F. 用微波驱动控制量子系统（Controlling quantum systems with microwave drives）

  

连接到传输线会引入损耗，但对控制与测量是必要的。本节用输入-输出框架描述经典微波驱动。

  

### F.1 连续波驱动的输入场建模：$b_{\mathrm{in}}(t)$ 加上经典部分

考虑对谐振子输入端施加频率 $\omega_d$、相位 $\phi_d$ 的连续微波。用 IV.B 的输入-输出观点，可通过将式(75)中的输入场替换为

  

$$

\hat b_{\mathrm{in}}(t)\ \to\ \hat b_{\mathrm{in}}(t)+\beta(t),

$$

  

其中经典驱动（输入场的经典部分）

$$

\beta(t)=A(t)\,e^{-i\omega_d t-i\phi_d},

$$

$A(t)$ 是（可随时间变化的）驱动幅度。

  

这样在朗之万方程中会出现 $\sqrt{\kappa}\beta(t)$ 的经典驱动项，并可等效吸收到系统哈密顿量中：$\hat H_S\to \hat H_S+\hat H_d$。

  

---

  

### F.2 驱动哈密顿量形式

驱动对应的哈密顿量为（文中给出）

  

$$

\hat H_d=\hbar\Big[\epsilon(t)\hat a^\dagger e^{-i\omega_d t-i\phi_d}

+\epsilon^*(t)\hat a\, e^{i\omega_d t+i\phi_d}\Big],

\tag{84}

$$

  

其中

$$

\epsilon(t)= i\sqrt{\kappa}\,A(t)

$$

是谐振子“看到的”驱动强度（可随时间变）。

  

- 推广到多个驱动、或直接驱动 transmon，在形式上是直接的（更换相应算符与端口耦合）。

  

---

  

### F.3 驱动产生相空间位移：相干态与位移算符

$\hat H_d$ 本质上是谐振子相空间中的**位移（displacement）**生成元。适当选择驱动参数可将谐振子态从真空制备到任意相干态（Carmichael 2002; Gardiner & Zoller 1999）：

  

$$

|\alpha\rangle=\hat D(\alpha)|0\rangle

= e^{-|\alpha|^2/2}\sum_{n=0}^\infty \frac{\alpha^n}{\sqrt{n!}}|n\rangle.

\tag{85}

$$

  

位移算符为

$$

\hat D(\alpha)=e^{\alpha \hat a^\dagger-\alpha^*\hat a}.

\tag{86}

$$

  

文中提示：相干态在下一节的电路 QED 读出中起关键作用（例如读出腔在驱动下建立不同的指针态）。

  

---

  

### F.4 近似来源与适用性：驱动项来自 RWA，大振幅时会偏离

文中强调：$\hat H_d$ 的写法来自旋转波近似（RWA），可从类似 Eq.(20) 的形式理解；在该近似下驱动常写成

  

$$

i\hbar\,\epsilon(t)\cos(\omega_d t+\phi_d)\,(\hat a^\dagger-\hat a).

$$

  

对多数常见实验参数足够好，但当驱动幅度很大时，可能观察到偏离式(84)预测的效应（Pietikäinen et al. 2017; Verney et al. 2019）。

  

---


好的，我已将您提供的笔记内容中的数学环境修正为正确的 LaTeX 格式，并统一了排版风格，使其更清晰易读。

---

## 测量章节信息笔记（Circuit QED）

> 覆盖范围：相干态/位移算符的承接内容 + **V. Measurements in circuit QED** 中的 **A. Microwave field detection** 与 **B. Phase-space representations...**（含图 14–16 与式 (87)–(102) 等）。

---


---

### V. CIRCUIT QED 中的测量（Measurements in circuit QED）

#### 1) circuit QED 出现前的测量困难与 circuit QED 的优势
- 早期：用外加的测量器件（例如单电子晶体管）靠近超导量子比特测量。
- 核心挑战：读出电路在测量时要**强耦合**以快速提取信息（时间尺度小于 $T_1$），但测量关闭时又要与量子比特**良好退耦**以避免回馈（back-action）。
- 测量通常涉及耗散（dissipation），同时满足上述两点困难。
- circuit QED 的优势（文中点出）：
  1. 读出通过测量与量子比特耦合的**谐振腔上探测信号（probe tone）的散射**实现；读出开/关比（on/off ratio）更好（只有在探测信号存在时才发生读出）。
  2. 必要耗散可发生在远离量子比特处（例如室温电压表/电子学），而不是靠近芯片处。
  3. 在**色散（dispersive）**工作区，量子比特有效不吸收探测光子；对量子比特的回馈主要趋向于不可避免的退相干 → 接近 **QND（quantum non-demolition）** 读出理想。

#### 2) 微波频段测量的现实条件
- 微波光子能量小，单光子探测器仍在发展中。
- 因此 circuit QED 中常用策略：先用**近量子极限放大器**对弱微波信号放大，再用**外差/同相（heterodyne/homodyne）**方式测量场正交分量。

---

### A. Microwave field detection（微波场探测）

#### A1. 典型测量链路（图 14）
**链路结构（从信号源到数字处理）**
- 室温微波源产生 RF 信号 → 经多级衰减器（Att.）进入低温腔输入端。
- 衰减器作用：吸收从室温沿线路向样品传播的**热噪声**（thermal noise），并实现热锚定（不同温区分段衰减）。
- 腔输出场先经过**环行器（circulator）**再进入放大链：
  - 环行器：只允许信号正向传输，强烈抑制反向传播（来自放大器的噪声回灌到腔/样品）。
  - 工程限制：环行器常依赖永磁体，体积大、难片上集成，还引入插入损耗与片外损耗；目前有研发超导片上环行器的努力。
- 第一极：接近量子极限的低温放大（如 JPA/TWPA 等）。
- 后续：HEMT 放大（低温高电子迁移率晶体管放大器）。
- 然后与本振 LO 在 **IQ mixer** 混频，下变频到 IF。
- IF 信号由 ADC 采样，再由 FPGA 实时处理/记录。

---

#### A2. 有限带宽下的“滤波输出场”（式 (87)）
- 实际测量链路与电缆具有有限带宽。为简化讨论，引入滤波后的输出场模式 $\hat{a}_f(t)$：
  $$
  \hat{a}_f(t)=(f*\hat{b}_{out})(t)
  =\int_{-\infty}^{\infty}d\tau\, f(t-\tau)\hat{b}_{out}(\tau)
  =\int_{-\infty}^{\infty}d\tau\, f(t-\tau)\big[\sqrt{\kappa}\hat{a}(\tau)+\hat{b}_{in}(\tau)\big].
  $$
  - 最后一行使用输入-输出边界条件（文中引用此前方程）。
  - 归一化条件（确保为规范玻色模）：
    $$
    \int_{-\infty}^{\infty}dt\,|f(t)|^2=1 \quad\Rightarrow\quad [\hat{a}_f(t),\hat{a}_f^\dagger(t)]=1 .
    $$
- 滤波函数的用途：
  - 表征测量带宽；
  - 在读出中可用来优化不同量子比特态的可区分性（distinguishability）。

---

#### A3. 放大：相位保持放大器的输入输出关系（式 (88)）
- 忽略环行器并假设第一极为**相位保持（phase-preserving）**放大器（两个正交分量等增益）。
- 放大变换（式 (88)）：
  $$
  \hat{a}_{amp}=\sqrt{G}\,\hat{a}_f+\sqrt{G-1}\,\hat{h}^\dagger ,
  $$
  - $G$：功率增益（power gain）。
  - $\hat{h}^\dagger$：放大器引入的附加噪声（用“idler mode”表征）。
- 物理含义：
  - 理想相位保持放大器是对“信号模 + idler 模”的 Bogoliubov 变换。
  - 为保持输出算符满足玻色对易关系 $[\hat{a}_{amp},\hat{a}_{amp}^\dagger]=1$，必须引入附加噪声项。
  - “无附加噪声地同时放大两个正交分量”会与海森堡不确定性矛盾。

**典型增益/噪声（文中给出量级）**
- 近量子极限放大器：常见 $\sim 20\ \mathrm{dB}$ 增益。
- 后接 HEMT：进一步放大 $\sim 30\text{–}40\ \mathrm{dB}$，但附加噪声更大。
- 低温 HEMT（4–8 GHz）噪声数可低至 $\langle \hat{h}^\dagger\hat{h}\rangle\sim 5\text{–}10$。
- 注意：前级到后级之间的衰减会显著恶化有效噪声表现。

---

#### A4. 含衰减与多级放大的总附加噪声（图 15 与式 (89)）
- 用透过率为 $\eta_1,\eta_2$ 的“分束器（beam splitter）”来建模链路中的衰减（等效为与真空端口耦合引入真空噪声），随后串接两级放大器增益 $G_1,G_2$，噪声模占据数 $N_i=\langle \hat{h}_i^\dagger\hat{h}_i\rangle$（$i=1,2$），总等效噪声数为 $N_T=\langle \hat{h}_T^\dagger\hat{h}_T\rangle$。
- 总附加噪声（式 (89)）：
  $$
  N_T=\frac{1}{G_T}\left[
  \frac{1}{\eta_1}(G_1-1)G_2(N_1+1)+(G_2-1)(N_2+1)
  \right]-1
  \approx \frac{1}{\eta_1}\left[1+N_1+\frac{N_2}{\eta_2G_1}\right]-1 .
  $$
  - $G_T=\eta_1\eta_2G_1G_2$（总增益，文中在式后给出这一含义）。
  - 近似式对应“大增益极限”下的主导项分析。
- 结论（文中强调）：
  - 若第一极增益 $G_1$ 足够大，总噪声主要由第一极噪声 $N_1$ 决定；
  - 因此必须在链路第一极使用近量子极限低噪声放大器。

**量子效率（quantum efficiency）定义**
- 常用定义：
  $$
  \eta=\frac{1}{N_T+1},
  $$
  理想情况 $N_T=0 \Rightarrow \eta=1$。

---

#### A5. 量子效率的另一种等效定义（图 15(b)，式 (90)–(92)）
- 另一种模型：把“有噪声放大器（增益 $G$）”等效成
  - 一个“无噪声放大器（增益 $G/\bar{\eta}$）”
  - 前接一个虚拟分束器（透过率 $\bar{\eta}$），通过向输入端注入真空噪声来等效附加噪声。
- 该网络输出（文中给出）可写为：
  $$\hat{a}_{amp}=\sqrt{G/\bar{\eta}}\big(\sqrt{\bar{\eta}}\hat{a}_f+\sqrt{1-\bar{\eta}}\hat{v}\big)$$
  （$\hat{v}$ 为真空模）。
- 输出幅度平方期望（式 (90)）：
  $$
  \langle |\hat{a}_{amp}|^2\rangle
  =\frac{G}{\bar{\eta}}\left[\left(1-\bar{\eta}^{-1}\right)\frac{1}{2}+\bar{\eta}\langle |\hat{a}_f|^2\rangle\right],
  $$
  其中
  $$
  \langle |O|^2\rangle=\frac{\langle \{O^\dagger,O\}\rangle}{2}
  $$
  为对称化涨落（symmetrized fluctuations）。
- 同一结果也可写为（式 (91)）：
  $$
  \langle |\hat{a}_{amp}|^2\rangle = G\big(A+\langle |\hat{a}_f|^2\rangle\big),
  $$
  并定义附加噪声（式 (92)）：
  $$
  A=\left(\frac{G-1}{G}\right)\left(\langle \hat{h}^\dagger\hat{h}\rangle+\frac{1}{2}\right).
  $$
- 低噪声且大增益极限：
  - 当 $\langle \hat{h}^\dagger\hat{h}\rangle\to 0$ 且 $G\gg 1$ 时，
    $$
    A\ge \frac{1-G^{-1}}{2}\approx \frac{1}{2},
    $$
    即相位保持放大器的附加噪声下界对应“半个光子”的噪声（文中引用 Caves 结果）。
- 由该定义得到量子效率（文中给出）：
  $$
  \bar{\eta}=\frac{1}{2A+1}\le \frac{1}{2},
  $$
  且理想情况下上界为 $1/2$。
- 文中强调：量子效率概念不仅适用于放大器，也可扩展到整条测量链路（如图 14）。

---

#### A6. 放大后的电压表达式（式 (93)）
- 结合此前关系（文中引用式 (76) 与 (88)），放大后电压可写为（式 (93)）：
  $$
  \hat{V}_{amp}(t)\simeq \sqrt{\frac{\hbar\omega_{RF}Z_{tm}}{2}}
  \left[e^{-i\omega_{RF}t}\hat{a}_{amp}+\mathrm{H.c.}\right].
  $$
  - $\omega_{RF}$：信号频率。
  - 文中说明：为简化，忽略有限电缆长度相关相位；并忽略朝向放大器反向传播输入场的贡献（该贡献未被放大，通常较小）。

---

#### A7. IQ 混频与外差检测（图 16，式 (94)–(96)）
- 下一阶段使用 IQ mixer 从放大后的 RF 信号提取信息。
- IQ mixer 结构要点（图 16）：
  - 先经功分（power divider；图中以分束器示意）→ 两路；
  - 两路分别与 LO 混频，其中一路 LO 相位移 $\pi/2$；
  - 同时获得两正交分量的测量（I/Q）。
- 本振（LO）电压（式 (94)）：
  $$
  V_{LO}(t)=A_{LO}\cos(\omega_{LO}t-\phi_{LO}).
  $$
- 单个混频器输出（式 (95)）：
  $$
  V_{mixer}(t)=K V_{RF}(t)V_{LO}(t)
  =\frac{1}{2}KA_{LO}A_{RF}\Big[
  \cos((\omega_{LO}-\omega_{RF})t-\phi_{LO})
  +\cos((\omega_{LO}+\omega_{RF})t-\phi_{LO})
  \Big],
  $$
  - $K$：电压转换损耗系数（conversion losses）。
  - 低通滤波去除高频项，保留中频 $\omega_{IF}=\omega_{LO}-\omega_{RF}$ → **heterodyne detection（外差）**。
- 选择 $\omega_{IF}$：
  - 实验上通常为几十到几百 MHz；
  - 数字化由 ADC 采样率与 IF 频带设置决定；带宽在 circuit QED 中常为 MHz–几十 MHz（由信号带宽、腔线宽 $\kappa/2\pi$ 等决定）。
  - 信号随后可做平均、实时 FPGA 分析或离线处理。

**量子形式的 I/Q 输出（式 (96a)(96b)）**
- 合并式 (93) 与 (95)，并令 I 路 $\phi_{LO}=0$、Q 路 $\phi_{LO}=\pi/2$，得到：
  $$
  \hat{V}_I(t)=V_{IF}\Big[\hat{X}_f(t)\cos(\omega_{IF}t)-\hat{P}_f(t)\sin(\omega_{IF}t)\Big]+\hat{V}_{noise,I}(t),
  $$
  $$
  \hat{V}_Q(t)= -V_{IF}\Big[\hat{P}_f(t)\cos(\omega_{IF}t)+\hat{X}_f(t)\sin(\omega_{IF}t)\Big]+\hat{V}_{noise,Q}(t),
  $$
  其中
  $$
  V_{IF}=KA_{LO}\sqrt{\frac{GZ_{tm}\hbar\omega_{RF}}{2}},
  $$
  而 $\hat{V}_{noise,I/Q}(t)$ 汇集放大器噪声及其它附加噪声贡献。

---

#### A8. 场正交分量与旋转参考系（式 (97)–(98)）
- 定义输出（滤波）模式的两正交分量（式 (97)）：
  $$
  \hat{X}_f=\frac{\hat{a}_f^\dagger+\hat{a}_f}{2},\qquad
  \hat{P}_f=i\frac{\hat{a}_f^\dagger-\hat{a}_f}{2},
  $$
  并满足对易关系 $[\hat{X}_f,\hat{P}_f]=i/2$。
- 物理含义：
  - $\hat{V}_I(t),\hat{V}_Q(t)$ 随时间描绘在 $x_f\text{–}p_f$ 平面上的轨迹并携带两正交分量信息。
  - 可通过数模/数字处理将信号变换到“静止（stationary）”参考系。
- 旋转矩阵（式 (98)）：
  $$
  R(t)=
  \begin{pmatrix}
  \cos(\omega_{IF}t) & -\sin(\omega_{IF}t)\\
  \sin(\omega_{IF}t) & \cos(\omega_{IF}t)
  \end{pmatrix}.
  $$

---

#### A9. 同相（homodyne）检测作为特例（式 (99)）
- 当 $\omega_{IF}=0$ 时称为 **homodyne detection（同相检测）**。
- 某一路 IQ 输出与下式成正比（式 (99)）：
  $$
  \hat{X}_{f,\phi_{LO}}=
  \frac{\hat{a}_f^\dagger e^{i\phi_{LO}}+\hat{a}_f e^{-i\phi_{LO}}}{2}
  =\hat{X}_f\cos\phi_{LO}+\hat{P}_f\sin\phi_{LO}.
  $$
- 文中指出的利弊：
  - 表面更简单、似乎更有优势；
  - 但因为同相信号在直流（DC），更易受 $1/f$ 噪声与漂移影响；
  - 与光学中可无噪声同相探测不同：此处在放大器与 IQ mixer 的噪声端口处引入噪声。

---

### B. Phase-space representations and their relation to field detection
（相空间表象及其与场探测的关系）

#### B1. 为什么引入相空间表象
- 在场探测语境中，用相空间表示电磁场量子态特别方便。
- 文中聚焦两类表象：**Wigner 函数** 与 **Husimi-Q 分布**（此处先展开 Wigner 的定义）。

#### B2. Wigner 函数与特征函数（式 (100)–(101)）
- Wigner 函数是准概率分布，可由特征函数的傅里叶变换给出（式 (100)）：
  $$
  W_\rho(x,p)=\frac{1}{\pi^2}\iint_{-\infty}^{\infty}dx'\,dp'\,
  C_\rho(x',p')\,e^{2i(px'-x p')}.
  $$
- 特征函数定义为（式 (101)）：
  $$
  C_\rho(x,p)=\mathrm{Tr}\Big\{\rho\,e^{2i(p\hat{X}-x\hat{P})}\Big\}.
  $$

#### B3. 与位移算符的联系（式 (102)）
- 对电磁场态 $\rho$，特征函数可理解为位移算符的期望值。
- 位移算符的相空间形式（式 (102)）：
  $$
  \hat{D}(\alpha)=e^{2i(p\hat{X}-x\hat{P})}=e^{\alpha\hat{a}^\dagger-\alpha^*\hat{a}},
  \qquad \alpha=x+ip.
  $$
- 文中连接点：
  - 这与前述相干态（式 (85)）一致：相干态具有特别简单的 Wigner 函数结构（后续图 17 与后文将继续展开）。

---

### 图示要点速记
- **图 14**：完整微波测量链路（衰减器、环行器、低温近量子极限放大、HEMT、IQ mixer、ADC/FPGA），并标注不同温区（室温/4K/10mK）工作位置。
- **图 15(a)**：两级放大 + 两段衰减的链路模型，用分束器透过率 $\eta_1,\eta_2$ 等效衰减并引入真空噪声端口。
- **图 15(b)**：把“有噪声放大器”改写成“无噪声放大器 + 等效输入损耗（虚拟分束器）”，引出另一种量子效率定义与 $\le 1/2$ 的上界。
- **图 16**：IQ mixer 原理：RF 经功分后与 LO 混频，其中一路 LO 移相 $\pi/2$，从而同时测得 I/Q 两路（两正交分量）。

---

好的，我已将您提供的续篇笔记内容中的数学环境修正为正确的 LaTeX 格式，并统一了排版风格。

---

## 测量章节信息笔记（续）

> 覆盖范围：图 17–19；式 (103)–(120)；小节 **B** 的延续（Wigner/Q 与检测的关系）、小节 **C. Dispersive qubit readout**（含稳态腔场、SNR/保真度、其它方案），以及 **VI. Qubit-resonator coupling regimes** 的开头。

---

### B. Phase-space representations and their relation to field detection（续）

#### B4. 相干态的 Wigner 函数（图 17，式 (103)）
- 图 17：示意相干态在相空间中的分布，以及沿旋转正交分量轴 $\hat{X}_\phi$ 的边缘分布（marginal）。
- 相干态 $\lvert \beta\rangle$ 的 Wigner 函数为以 $\beta$ 为中心的高斯（式 (103)）：
  $$
  W_{\lvert \beta\rangle}(\alpha)=\frac{2}{\pi}e^{-2|\alpha-\beta|^2}.
  $$
- 高斯宽度 $1/\sqrt{2}$ 是量子真空噪声的特征；相干态饱和海森堡不等式：
  $$
  \Delta X\,\Delta P=\frac{1}{2},\qquad
  \Delta O^2=\langle \hat{O}^2\rangle-\langle \hat{O}\rangle^2 .
  $$
- 与相干态不同：**非经典态**的 Wigner 函数可取负值（Wigner negativity）。

#### B5. Wigner 边缘分布与正交分量测量（式 (104a)(104b)）
- 在色散读出等场测量语境中，Wigner 函数与理想正交分量测量结果概率分布直接相关。
- 对应两条边缘分布（式 (104)）：
  $$
  P(x)=\int_{-\infty}^{\infty}dp\,W_\rho(x,p)=\langle x\lvert\rho\rvert x\rangle,
  \qquad (104a)
  $$
  $$
  P(p)=\int_{-\infty}^{\infty}dx\,W_\rho(x,p)=\langle p\lvert\rho\rvert p\rangle.
  \qquad (104b)
  $$
- 因而理想 homodyne 测量某个旋转正交分量 $\hat{X}_\phi$ 的结果分布 $P(X_\phi)$，可视为对 $W_\rho$ 沿正交轴 $\hat{X}_{\phi+\pi/2}$ 积分得到（图 17 直观展示）。

#### B6. Husimi-Q 分布与异频（heterodyne）检测（式 (105)–(106)）
- Husimi-Q 分布定义（式 (105)）：
  $$
  Q_\rho(\alpha)=\frac{1}{\pi}\langle \alpha\lvert\rho\rvert \alpha\rangle .
  $$
  - 物理含义：在相干态 POVM 下“得到结果为 $\alpha$”的概率分布。
  - 性质：与 Wigner 函数不同，$Q_\rho(\alpha)$ **总为非负**。
- Q 与 Wigner 的关系（式 (106)）：
  $$
  Q_\rho(\alpha)=\frac{2}{\pi}\int_{-\infty}^{\infty}d^2\beta\,
  W_\rho(\beta)e^{-2|\alpha-\beta|^2}
  =W_\rho(\alpha)*W_{\lvert 0\rangle}(\alpha).
  $$
  - 结论：Q 分布是 Wigner 与高斯核（真空态 Wigner）的卷积，因此更“平滑”；可理解为在 Wigner 的基础上加入真空噪声。
- 与 IQ 异频检测的关键联系：
  - IQ mixer 的 heterodyne 检测在测量前等效加入（理想）真空噪声；
  - 因此“同时测两条正交分量”的统计分布，与其说对应 Wigner 的边缘，不如说对应 **Husimi-Q 分布的边缘**（文中据此得出结论）。

---

### C. Dispersive qubit readout（色散量子比特读出）

#### C1. 稳态腔内场（Steady-state intra-cavity field）

##### C1.1 色散哈密顿量与“指针态”（式 (107)–(108)）
- 色散区（dispersive regime）下，transmon–resonator 哈密顿量可近似为（式 (107)）：
  $$
  \hat{H}_{disp}\approx \hbar(\omega_r+\chi\hat{\sigma}_z)\hat{a}^\dagger\hat{a}
  +\frac{\hbar\omega_q}{2}\hat{\sigma}_z .
  $$
  文中说明的近似/截断：
  - 将 transmon 截断为两能级；
  - 吸收了频率的 Lamb shift；
  - 忽略了腔的 transmon 诱导非线性（Kerr 非线性 $K$，此前在式 (50) 讨论过）。
- 物理结论：腔频率随量子比特态改变  
  - 若比特在 $\lvert g\rangle$，则 $\langle \hat{\sigma}_z\rangle=-1$，腔频率为 $\omega_r-\chi$；
  - 若在 $\lvert e\rangle$，则 $\langle \hat{\sigma}_z\rangle=+1$，腔频率为 $\omega_r+\chi$。
- 因此驱动腔会产生与比特态相关的相干态（指针态）$\lvert \alpha_g\rangle,\lvert\alpha_e\rangle$；若比特初态为叠加 $c_g\lvert g\rangle+c_e\lvert e\rangle$，系统演化为纠缠态（式 (108)）：
  $$
  c_g\lvert g,\alpha_g\rangle+c_e\lvert e,\alpha_e\rangle .
  $$
- 文中用 Stern–Gerlach 类比解释“测量如何发生”：通过与可观测量相关的场态分离，将“自旋信息”映射为可读出的连续变量分布。

##### C1.2 腔内相干幅度的动力学（式 (109a)(109b)，图 18）
- 在 Langevin 方程框架下（文中引用先前式 (75)），腔内复振幅 $\alpha_\sigma(t)$（$\sigma=g,e$）满足（式 (109)）：
  $$
  \dot{\alpha}_e(t)=-i\varepsilon(t)-(\delta_r+\chi)\alpha_e(t)-\frac{\kappa}{2}\alpha_e(t),
  \qquad (109a)
  $$
  $$
  \dot{\alpha}_g(t)=-i\varepsilon(t)-(\delta_r-\chi)\alpha_g(t)-\frac{\kappa}{2}\alpha_g(t),
  \qquad (109b)
  $$
  其中 $\delta_r=\omega_r-\omega_d$ 为测量驱动相对“裸腔频率”的失谐，$\kappa$ 为腔能量衰减率，$\varepsilon(t)$ 为驱动幅度（复包络）。
- 图 18（a–c）：展示在相空间中从初态到稳态的轨迹；并对比不同 $2\chi/\kappa$ 时两指针态的分离与边缘分布的重叠（圈半径 $1/\sqrt{2}$ 表示真空噪声尺度）。

##### C1.3 稳态解与稳态正交分量（式 (110)–(113)）
- 在稳态（$\dot{\alpha}_\sigma=0$）且驱动恒定时（式 (110)）：
  $$
  \alpha^{s}_{e/g}=\frac{-\varepsilon}{(\delta_r\pm\chi)-i\kappa/2},
  \qquad (110)
  $$
  其中对 $e$ 取 $+$、对 $g$ 取 $-$。
- 对应稳态腔内正交分量期望（式 (111)）：
  $$
  \langle \hat{X}\rangle^{s}_{e/g}=
  \frac{\varepsilon(\delta_r\pm\chi)}{(\delta_r\pm\chi)^2+(\kappa/2)^2},
  \qquad (111a)
  $$
  $$
  \langle \hat{P}\rangle^{s}_{e/g}=
  \frac{\varepsilon\,\kappa/2}{(\delta_r\pm\chi)^2+(\kappa/2)^2}.
  \qquad (111b)
  $$
- 当驱动在裸腔频率（$\delta_r=0$）时，文中指出：关于比特的信息只在 $X$ 正交分量中（图 18(a–c) 的结论）。
- 定义稳态幅度（式 (112)）：
  $$
  A^s_{e/g}=\sqrt{\langle \hat{X}\rangle_{e/g}^2+\langle \hat{P}\rangle_{e/g}^2}
  =\frac{2\varepsilon}{\sqrt{(\kappa/2)^2+(\delta_r\pm\chi)^2}} .
  $$
- 定义稳态相位（式 (113)）：
  $$
  \phi^s_{e/g}=\arctan\!\left(\frac{\langle \hat{X}\rangle^s_{e/g}}{\langle \hat{P}\rangle^s_{e/g}}\right)
  =\arctan\!\left(\frac{\delta_r\pm\chi}{\kappa/2}\right).
  $$
- 图 19：给出两比特态下的透射谱（虚线）与相移（实线）。驱动频率靠近“被拉移的共振频率”时，响应强烈依赖比特态；由此可通过透射/反射信号的幅度与/或相位实现读出。

##### C1.4 实验上测“外部场”与端口耦合结构的要点（文字总结段落）
- 以上用腔内场幅相表述；实际通常测量耦合到谐振腔的**外部传输线中的场**（透射或反射），其与腔内场的关系由 input–output 理论给出（文中指向 Sec. IV 与附录）。
- 非对称腔（一个端口强耦合、另一端口弱耦合）：
  - 若从弱端口驱动，关于比特态的大部分信息在强端口辐射出的场中；
  - 若从强端口驱动，输出为“直接反射的驱动 + 腔辐射”的叠加；扫频穿过共振时相位演化特性不同（文中指出一类情况相移为 $\pi$，另一类为 $2\pi$），影响对色散频移的敏感度。
- 对称腔：比特信息被“分到两个端口”，导致该构型读出效率不佳。

---

#### C2. 信噪比与测量保真度（Signal-to-noise ratio and measurement fidelity）

##### C2.1 积分测量算符与权重函数（式 (114)）
- IQ 异频得到 $\hat{X}_f(t),\hat{P}_f(t)$ 后，考虑对测量记录积分（积分时间 $\tau_m$），定义测量算符（式 (114)）：
  $$
  \hat{M}(\tau_m)=\int_0^{\tau_m}dt\,\Big\{
  w_X(t)\big[V_{IF}\hat{X}_f(t)+\hat{V}_{noise,X}(t)\big]
  +w_P(t)\big[V_{IF}\hat{P}_f(t)+\hat{V}_{noise,P}(t)\big]
  \Big\}.
  $$
- 权重函数选择（文中给出常用选法与物理理由）：
  $$
  w_X(t)=\langle \hat{X}_f\rangle_e-\langle \hat{X}_f\rangle_g,\qquad
  w_P(t)=\langle \hat{P}_f\rangle_e-\langle \hat{P}_f\rangle_g,
  $$
  用于增强两态可区分性。
  - 直观解释：由于能量弛豫，长时间处权重应更小（晚时刻更倾向于读到基态信息）。
  - 对图 18 所示情形：$P$ 正交分量不含比特信息，因此取 $w_P(t)=0$，避免把该正交分量的噪声积分进来。

##### C2.2 SNR 定义（式 (115)）与“高斯近似”
- 文中指出：结合前述 A、B 小节，$\hat{M}(\tau_m)$ 多次重复测量的结果分布可视为高斯分布，并可由腔场 Q 分布的边缘刻画。
- 定义信噪比（式 (115)）：
  $$
  \mathrm{SNR}^2(t)=
  \frac{\big|\langle \hat{M}(t)\rangle_e-\langle \hat{M}(t)\rangle_g\big|^2}
  {\langle \hat{M}_N^2(t)\rangle_e+\langle \hat{M}_N^2(t)\rangle_g},
  $$
  其中 $\langle \hat{M}\rangle_\sigma$ 表示量子比特在态 $\sigma$ 时的平均积分信号，
  $$
  \hat{M}_N=\hat{M}-\langle \hat{M}\rangle
  $$
  为噪声算符（包含附加噪声与量子真空噪声贡献）。

##### C2.3 测量保真度（式 (116)）与与 SNR 的关系
- 定义测量保真度（式 (116)）：
  $$
  F_m=1-\big[P(e\lvert g)+P(g\lvert e)\big]=1-E_m .
  $$
  其中 $P(\sigma\lvert \sigma')$ 表示真实态为 $\sigma'$ 却判为 $\sigma$ 的概率；$E_m$ 为测量误差。
- 文中说明：对图 18(a) 的几何解释，$E_m$ 可视为两态分布（对应两条 Q 分布边缘）的重叠面积（并提到可通过选择 LO 相位最小化）。
- 若边缘分布为高斯，则保真度与 SNR 的关系为（文中给出）：
  $$
  F_m=1-\mathrm{erfc}\!\left(\frac{\mathrm{SNR}}{2}\right),
  $$
  其中 $\mathrm{erfc}$ 为互补误差函数。
- 重要限定：该结论仅在分布近似高斯时成立；实际中量子比特弛豫、色散哈密顿量高阶项、以及 Kerr 非线性等会导致指针态畸变并产生非高斯边缘。文中指出 circuit QED 中 Kerr 非线性可将相干态拉伸成“香蕉形”（bananaization）。

> 旁注（脚注给出的另一种定义）：常见的 *assignment fidelity* 定义为  
> $$
> F_a = 1-\frac{1}{2}\big[P(e\lvert g)+P(g\lvert e)\big],
> $$
> 其取值在 $[0,1]$，而上面的 $F_m$ 形式上可在 $[-1,1]$。

##### C2.4 长时间极限的 SNR 标度与参数优化（式 (117) 及随后的讨论）
- 在 $\delta_r=0$ 且忽略增益/混频等前因子时，长时间极限的 SNR（式 (117)）：
  $$
  \mathrm{SNR}(\tau_m\rightarrow\infty)\propto
  \left(\frac{2\varepsilon}{\kappa}\right)\sqrt{2\kappa\tau_m}\,|\sin 2\phi|,
  $$
  其中 $\phi$ 由式 (113) 给出。
- 由此可见：
  - $\mathrm{SNR}\propto \sqrt{\tau_m}$ 是“相干态读出”的标准量子极限标度；
  - 选择 $\chi/\kappa=1/2$ 可最大化该长时间极限 SNR（文中指出实验上常用该比值；但对有限测量时间并非总最优，图 18(d) 展示了最优点会随积分时间偏移）。
- 小 $\chi$ 极限下，文中将前因子 $2\varepsilon/\kappa$ 解释为与稳态平均测量光子数（腔内测量光子数）相关的量，从而提示：增大测量光子数可提高 SNR。

##### C2.5 测量光子数的上限：色散近似失效与非 QND
- 文中强调：提高测量光子数 $\bar{n}$ 不能无限进行，因为会破坏用于推导色散哈密顿量式 (107) 的近似。
- 色散近似的关键小参数并非 $g/\Delta$，而是 $\bar{n}/n_{crit}$（此前定义过临界光子数 $n_{crit}$）。
- 在到达 $\bar{n}/n_{crit}\sim 1$ 之前，高阶项可能已显著，且实验上常见色散测量在 $\bar{n}\sim 1\text{–}10$ 就出现 QND 特性下降与“虚假翻转”（spurious flips）。
- 文中提到：由“dressed dephasing”对 $\bar{n}\ll n_{crit}$ 的非 QND 预期与实验并不总一致；另指出在 $\bar{n}\sim n_{crit}$ 的跃迁可能与系统中的偶然共振有关（引用相关工作）。

##### C2.6 读出速度与 Purcell 衰减的约束；最优参数的受限性
- 读出需尽量快于 $T_1$。
- 加快读出的一种策略：增大 $\kappa$（低 Q 腔），使测量光子更快泄露到可检测端口；但需避免 Purcell 效应导致的量子比特额外衰减。
- 解决手段：在谐振腔输出端加入 Purcell filter（文中引用相关实现）。
- 在固定 $\kappa$ 并取最优 $\chi/\kappa$ 时，文中指出达到稳态响应的时间尺度为 $1/\chi$。
- 增大 $\chi$ 的途径：更大耦合与更大非简谐性（如更大充电能 $E_C$）；但不能无限增大，因为色散与去相干会随 $E_C$ 增强而恶化。

##### C2.7 实验水平与放大器的重要性；量子效率与退相干的联系（式 (118)）
- 文中给出当代实验（示例）：
  - 在最小化读出时间时，$\tau_m=48\,\mathrm{ns}$ 可达 $F_m\approx 98.25\%$；
  - 在最大化保真度时，$\tau_m=88\,\mathrm{ns}$ 可达 $F_m\approx 99.2\%$；
  - 两者均使用约 $\bar{n}\approx 2.5$ 的腔内测量光子数（并指出主要限制来自相对较短的 $T_1$，如 $7.6\,\mu\mathrm{s}$ 的量级）。
- 强调近量子极限放大器的重要性：其出现使单次（single-shot / projective）读出成为可能，并可观测到 transmon 的量子跃迁（quantum jumps）。
- 整条测量链路的量子效率可由 SNR 提取（式 (118)）：
  $$
  \eta=\frac{\mathrm{SNR}^2}{4\beta_m},
  $$
  其中
  $$
  \beta_m=2\chi\int_0^{\tau_m}dt\,\mathrm{Im}\!\Big[\alpha_g(t)\alpha_e(t)^*\Big]
  $$
  与测量诱导退相干（measurement-induced dephasing）相关。
- 物理意义：体现“信息获取速率”与“不可避免的回馈（退相干）”之间的基本联系。

---

#### C3. 其它读出方案（Other approaches）

##### C3.1 Josephson Bifurcation Amplifier（JBA）
- 方案：使用含 Josephson 结的非线性 transmission-line resonator（Kerr 非线性，文中举例量级 $\sim -500\,\mathrm{kHz}$）。
- 在合适幅度/频率的相干驱动下系统发生分岔（bifurcation），从低光子数态跃迁到高光子数态。
- 通过色散耦合使分岔阈值/动力学依赖量子比特态，从而读出比特。

##### C3.2 高功率读出与量子比特 “punch out”
- 现象：在中等测量光子数下，读出非 QND 会损害保真度；但在非常大的测量功率极限，反而可恢复快速且高保真的单次读出。
- 直观解释（来自 Jaynes–Cummings 图景）：当 $\bar{n}\gg n_{crit}$ 时，腔对驱动的响应趋向于裸腔频率 $\omega_r$，仿佛“把量子比特打穿/打出去（punching out）”。
- 对多能级 transmon：该“量子到经典”的过渡功率依赖于 transmon 状态，可用于高保真测量；但代价是完全失去色散读出的 QND 性质。

##### C3.3 Squeezing（压缩）以提升 SNR 标度
- 文中指出 $\mathrm{SNR}\propto \sqrt{\tau_m}$ 的标度对应“相干态读出”的标准量子极限；若使用压缩输入（squeezed input），可在理想情况下达到更优的（Heisenberg-limited）标度，使 SNR 随测量光子数线性增长。
- 困难点：在高保真色散读出所需的强色散耦合下，指针态从相空间中心向稳态演化时会导致压缩角显著旋转；从而反压缩（anti-squeezed）正交分量的噪声会被混入，等效噪声增加。
- 文中提到可借助“两模压缩 / quantum-mechanics-free subsystem”等思路在特定方案中逼近 Heisenberg 标度；也提到“强驱动色散耦合腔”作为实现路径之一。

##### C3.4 Longitudinal readout（纵向读出）与合成纵向耦合（式 (119)–(120)）
- 纵向读出基于纵向耦合哈密顿量（文中引用此前的式 (63) 形式）：
  - 与色散耦合导致的“相空间旋转”不同，纵向耦合导致**与比特态条件相关的线性位移**，可更快产生指针态分离。
- 文中强调其优点：
  - $[\hat{H}_z,\hat{\sigma}_z]=0$，对应 QND；
  - 不受色散近似失效与 $n_{crit}$ 限制的同类问题影响，因此原则上可用更大的测量光子数。
- 但若 $\omega_r\gg g_z$，由纵向耦合产生的稳态位移（文中指向式 (110) 类比）很小，难以用于实际读出；因此需要“激活/增强”纵向耦合。
- 一种增强方式：调制耦合强度 $g_z(t)=\tilde{g}_z\cos(\omega_r t)$，在旋转参考系并略去快速振荡项后得到（式 (119)）：
  $$
  \hat{H}_z=\frac{\hbar\tilde{g}_z}{2}\,(\hat{a}^\dagger+\hat{a})\hat{\sigma}_z .
  $$
  此时稳态位移变为 $\pm \tilde{g}_z/\kappa$，可在中等调制幅度下变得显著。
- 另一种“合成纵向耦合”思路：强驱动色散耦合腔产生大位移 $\hat{a}\rightarrow \hat{a}+\alpha$。代入色散项并取 $\alpha$ 为实，得到（式 (120) 的结构）：
  $$
  \chi \hat{a}^\dagger\hat{a}\hat{\sigma}_z
  \rightarrow
  \chi \hat{a}^\dagger\hat{a}\hat{\sigma}_z
  +\chi\alpha(\hat{a}^\dagger+\hat{a})\hat{\sigma}_z
  +\chi\alpha^2\hat{\sigma}_z .
  $$
  - 中间项等效纵向耦合，幅度 $g_z=\alpha\chi$；
  - 在 $\chi$ 小但驱动强导致 $\alpha$ 大的极限，可使 $\alpha\chi$ 保持常数，从而得到有效纵向读出；
  - 文中给出直观解释：强驱动下指针态在相空间中旋转的圆半径很大，短时间内运动近似线性位移。
- 文中还提到一种更“微妙”的合成方式：对量子比特以 Rabi 频率 $\Omega_R$ 驱动，同时在谐振腔的边带频率 $\omega_r\pm\Omega_R$ 驱动，从而得到期望的纵向相互作用；但由于仍基于色散框架，仍会受到 Purcell 衰减与非 QND 的影响。

---
