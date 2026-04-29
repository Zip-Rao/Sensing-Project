# 差分回波协议（Differential Echo Protocol）实现方案

## 1. 文献核心内容总结

### 1.1 协议原理

文献 Zopes & Degen, Phys. Rev. Applied 12, 054028 (2019) 提出了一种**无需重建的量子传感方案**，通过差分自旋回波直接测量时变磁场波形。

**核心思想**：在两次连续的波形触发之间插入两个 $\pi$ 脉冲（分别位于时间 $t$ 和 $t + t_{\text{int}}$），选择性地积累时间窗口 $[t, t + t_{\text{int}}]$ 内的相位，同时取消其余时间段的相位贡献。

### 1.2 脉冲序列结构（Fig. 1(d)）

```
(π/2)_X — [波形通道1: π脉冲@t] — [波形通道2: π脉冲@(t+t_int)] — ... 重复k次 ... — (π/2)_Y
```

详细时序：
1. $(\pi/2)_X$ 脉冲：初始化叠加态
2. **重复 k 次**以下单元：
   - 第1次波形通道：在 $t$ 时刻施加 $\pi$ 脉冲
   - 第2次波形通道：在 $t + t_{\text{int}}$ 时刻施加 $\pi$ 脉冲
3. $(\pi/2)_Y$ 脉冲：读出

### 1.3 关键公式

**测量信号**（线性近似 $\phi \ll 1$）：

$$p(t) = 0.5 + 2k \cdot \gamma_e \cdot B(t) \cdot t_{\text{int}}$$

**磁场反演**：

$$B(t) = \frac{p(t) - 0.5}{2k \cdot \gamma_e \cdot t_{\text{int}}}$$

**调制函数**（理想短 $\pi$ 脉冲情况）：

$$M_1(t, t') = \begin{cases} +1 & t' \leq t \\ -1 & t' > t \end{cases}$$

$$M_2(t, t') = \begin{cases} -1 & t' \leq t + t_{\text{int}} \\ +1 & t' > t + t_{\text{int}} \end{cases}$$

**有限脉冲宽度的调制函数**（$t_{\text{int}} = t_\pi$ 时）：

$$M(t, t') = \begin{cases} 1 + \cos\left(\frac{\pi[t' - t - t_\pi/2]}{t_\pi}\right) & t - \frac{t_\pi}{2} < t' \leq t + \frac{3t_\pi}{2} \\ 0 & \text{otherwise} \end{cases}$$

此时调制函数等价于一个特征宽度为 $2t_\pi$ 的 Hann 窗函数。

**时间分辨率**：

$$\tau \approx \frac{2}{\pi}\sqrt{t_\pi^2 + \frac{\pi}{2}t_{\text{int}}^2}$$

**灵敏度**：

$$B_{\min,\text{Diff}} = \frac{\sqrt{t_m + 2kt_{\text{rep}}} \cdot \exp(2kt_{\text{rep}}/T_2)}{\gamma_e C \cdot 2k \cdot t_{\text{int}}}$$

其中 $t_m$ 为初始化和读出时间，$C$ 为读出效率，$T_2$ 为相干时间。

**最优重复次数**：$k_{\text{opt}} \approx T_2 / (4 t_{\text{rep}})$

### 1.4 关键实验参数（NV center）

| 参数 | 值 |
|------|------|
| $t_\pi$ | ~20 ns |
| $t_{\text{int}}$ | ~20 ns |
| $t_s$（采样步长） | 4-8 ns |
| $t_{\text{rep}}$（波形重复周期） | 344-1400 ns |
| $k$（重复次数） | 1-8 |
| 时间分辨率 $\tau$ | ~20 ns |
| 带宽 $f_{-3\text{dB}}$ | ~25 MHz |
| 灵敏度 | ~4 $\mu$T/$\sqrt{\text{Hz}}$ |

---

## 2. 与现有仿真平台的对应关系

### 2.1 物理系统映射

| 文献（NV center） | 仿真平台（Transmon qubit） |
|-------------------|--------------------------|
| 电子自旋 $S=1$，$\{m_s=0, m_s=-1\}$ | Transmon qubit，$\{\|0\rangle, \|1\rangle\}$ |
| 磁场 $B(t)$ 通过 $\gamma_e B$ 调制频率 | 磁通 $\Phi(t)$ 通过 $\kappa(\Phi_0) \cdot \Phi$ 调制频率 |
| 旋磁比 $\gamma_e = 2\pi \times 28$ GHz/T | 磁通灵敏度 $\kappa = \partial f / \partial \Phi$ |
| $T_2 \sim 14\ \mu s$ | $T_2 = 50\ \mu s$（可调） |
| Rabi频率 ~25 MHz | Rabi频率由脉冲幅度决定 |
| $\pi$ 脉冲 ~20 ns | $\pi$ 脉冲 ~10-40 ns（可调） |

### 2.2 代码模块映射

| 功能需求 | 对应模块 | 现有状态 |
|----------|---------|---------|
| 差分回波脉冲序列 | `pulse.py` → `create_differential_echo_pulse()` | **需新增** |
| 协议执行逻辑 | `protocal.py` → `case 5` | **需新增** |
| 波形采样（等效时间采样） | `protocal.py` → `sliding_measurement()` | **可复用**，需适配 |
| 磁场直接反演 | `analysis.py` → `get_signal_from_diff_echo()` | **需新增** |
| 信号生成 | `signal.py` | **可直接复用** |
| qubit 磁通调制 | `qubit.py` → `qubit_in_mag()` | **可直接复用** |

---

## 3. 可行性分析

### 3.1 完全可行的部分

1. **脉冲序列构建**：现有 `CompositePulse` 类完全支持构建任意复合脉冲序列，包括 $\pi/2$ 脉冲、$\pi$ 脉冲和自由演化段的组合。
2. **滑动测量框架**：`sliding_measurement()` 已实现通过改变脉冲与信号的时间延迟进行等效时间采样，这正是差分回波协议扫描采样时间 $t$ 的核心机制。
3. **信号直接反演**：差分回波协议的最大优势是无需反卷积，直接通过公式 $B(t) = (p(t) - 0.5) / (2k\kappa t_{\text{int}})$ 即可得到磁场波形。这比现有的 Wiener 反卷积和 LM 数值反演简单得多。
4. **qubit 磁通调制**：`qubit_in_mag()` 和 `qubit_under_mag()` 已经可以计算在时变磁通下的频率变化。

### 3.2 需要适配的部分

1. **波形重复触发机制**：文献中波形被触发 $2k$ 次，在仿真中需要在总演化时间内让波形重复 $2k$ 次，每次触发间隔 $t_{\text{rep}}$。这需要在 `Signal` 类中支持周期性重复信号，或在协议层面处理信号拼接。

2. **总演化时间**：差分回波协议的总时间为 $\sim 2k \cdot t_{\text{rep}} + t_{\pi/2}$，对于 $k=8$, $t_{\text{rep}}=350$ ns，总时间约 5600 ns，mesolve 的时间步数和计算量需要注意。

3. **$\pi/2$ 脉冲相位**：文献中第一个 $\pi/2$ 沿 $X$ 轴（phase=0），最后一个 $\pi/2$ 沿 $Y$ 轴（phase=$\pi/2$），这确保在零信号时测量概率恰好为 0.5。

### 3.3 物理限制

1. **线性近似**：公式 $p = 0.5 + 2k\kappa B t_{\text{int}}$ 要求 $\phi = 2k\kappa B t_{\text{int}} \ll 1$，对于 Transmon qubit，需要信号幅度足够小或 $k$ 不太大。超出线性范围时需要相位循环技术。
2. **相干时间约束**：要求 $2k \cdot t_{\text{rep}} \lesssim T_2$，即波形重复总时间不能超过相干时间。平台默认 $T_2 = 50\ \mu s$，足够支持多次重复。
3. **Transmon 非线性**：Transmon 频率对磁通的响应 $f(\Phi) = \sqrt{8E_J E_C |\cos(\pi\Phi)|} - E_C$ 是非线性的，在最佳工作点附近近似线性，但大信号会引入非线性误差。

### 3.4 与现有方法的优势比较

| 方面 | 现有 Ramsey 滑动测量（case 4） | 差分回波协议 |
|------|--------------------------|------------|
| 反演方法 | 需要反卷积/数值优化 | **直接读出，无需反演** |
| 噪声放大 | 反卷积放大噪声 | **无噪声放大** |
| 灵敏度 | 受限于单次 $t_{\text{int}}$ | **通过 $k$ 次重复提升** |
| 实现复杂度 | 较简单 | 较复杂（需要波形重复触发） |
| 适用条件 | 单次触发 | 需要波形可重复触发 |

---

## 4. 详细实现方案

### 4.1 Step 1: 新增差分回波脉冲序列（`pulse.py`）

新增函数 `create_differential_echo_pulse(t_rabi, t_int, t_rep, k, omega_d)`：

```python
def create_differential_echo_pulse(t_rabi, t_int, t_rep, k, omega_d, t_sample):
    """
    创建差分回波序列的复合脉冲对象
    
    脉冲序列结构：
    (π/2)_X — [free(t_sample) — π — free(t_rep - t_sample - t_pi)
              — free(t_sample + t_int) — π — free(t_rep - t_sample - t_int - t_pi)] × k — (π/2)_Y
    
    参数:
        t_rabi: π/2 脉冲的时间数组 (ns)
        t_int:  差分积分时间 (ns)，即两个 π 脉冲间隔对应的有效传感窗口
        t_rep:  波形重复周期 (ns)
        k:      重复次数
        omega_d: 驱动频率 (GHz)
        t_sample: 采样时间点 (ns)，即第一个 π 脉冲在波形中的位置
    """
```

**关键设计细节**：

- $\pi/2$ 脉冲的 Rabi 频率幅度：$\Omega = (\pi/2) / t_{\text{rabi,duration}}$
- $\pi$ 脉冲的 Rabi 频率幅度：$\Omega = \pi / t_{\text{rabi,duration}}$
- 第一个 $\pi/2$ 脉冲相位 `phase = 0`（沿 X 轴）
- 最后一个 $\pi/2$ 脉冲相位 `phase = π/2`（沿 Y 轴）
- 自由演化段使用 `Signal(type=0)` 零信号
- 重复单元内的时间安排确保 $\pi$ 脉冲准确落在 $t_{\text{sample}}$ 和 $t_{\text{sample}} + t_{\text{int}}$ 处

### 4.2 Step 2: 实现协议演化逻辑（`protocal.py`）

在 `Protocal.evolve()` 中新增 `case 5`（差分回波测量）：

```python
case 5:  # 差分回波测量
    # 1. 创建测试波形信号 Phi(t)
    # 2. 创建零信号基线 Phi_0(t)
    # 3. 对每个采样时间 t_sample 扫描:
    #    a. 构建差分回波脉冲序列 create_differential_echo_pulse(...)
    #    b. 在波形信号下执行 single_measurement()
    #    c. (可选) 在零信号下执行基线测量
    # 4. 收集 p_e(t) 数组
    # 5. 通过公式直接反演 B(t)
```

**具体实现方式有两种选择**：

#### 方式 A：复用 `sliding_measurement()` 框架（推荐）

核心思路：差分回波的采样过程本质上也是滑动测量——改变 $t_{\text{sample}}$ 等效于改变脉冲在波形中的采样位置。

不同之处在于：
- 滑动测量中脉冲序列不变，改变的是信号与脉冲的相对延迟
- 差分回波中，不同的 $t_{\text{sample}}$ 对应不同的脉冲序列（$\pi$ 脉冲位置不同）

因此需要在每个采样点重新构建脉冲序列。

#### 方式 B：直接循环实现

对每个 $t_{\text{sample}}$：
1. 构建包含 $2k$ 次波形通道的总磁场信号（将波形重复 $2k$ 次）
2. 构建对应的差分回波脉冲序列
3. 进行 mesolve 演化
4. 提取末态概率

**推荐方式 B**，因为差分回波的脉冲结构随采样点变化，直接循环更清晰。

### 4.3 Step 3: 磁场直接反演（`analysis.py`）

新增函数 `get_signal_from_diff_echo()`：

```python
def get_signal_from_diff_echo(self, qubit, p_e_list, t_int, k):
    """
    从差分回波测量结果直接反演磁场波形
    
    B(t) = (p(t) - 0.5) / (2k * kappa * t_int)
    
    参数:
        qubit:    TransmonQubit 对象
        p_e_list: 激发态概率数组, 对应各采样时间
        t_int:    积分时间 (ns)
        k:        重复次数
    返回:
        B:        反演磁场数组
    """
    kappa = qubit.frequency_sensitivity(qubit.flux)
    p_e = np.array(p_e_list)
    B = (p_e - 0.5) / (2 * k * kappa * t_int)
    return B
```

### 4.4 Step 4: 协议参数设计

**推荐的仿真参数**（适配现有平台的 Transmon qubit）：

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| $t_\pi$（$\pi$ 脉冲时长） | 10 ns | 对应 Rabi 频率 ~50 MHz |
| $t_{\text{int}}$（积分时间） | 10 ns | 取 $t_{\text{int}} = t_\pi$ |
| $t_s$（采样步长） | 4 ns | 横向分辨率 |
| $t_{\text{rep}}$（波形周期） | 200-500 ns | 需满足 $2k \cdot t_{\text{rep}} \ll T_2$ |
| $k$（重复次数） | 1, 2, 4, 8 | 可扫描比较 |
| 信号幅度 | 0.001-0.01 $\Phi_0$ | 确保线性近似成立 |
| 工作点 | $\Phi = \arctan(\sqrt{2})/\pi$ | 最佳灵敏度 |

### 4.5 Step 5: 验证方案

1. **方波验证**（对应文献 Fig. 2）：
   - 使用恒定信号（type=1）模拟方波
   - 比较 Ramsey 积分方法（case 1）和差分回波方法的输出
   - 验证差分回波可以直接还原方波，无噪声放大

2. **正弦波验证**（对应文献 Fig. 3）：
   - 使用正弦信号（type=2），频率在带宽内
   - 扫描 $k = 1, 2, 4, 8$，验证信号增益与 $k$ 的线性关系
   - 验证灵敏度随 $k$ 的改善

3. **复杂波形验证**（对应文献 Fig. 4）：
   - 使用复杂信号（type=7）
   - 验证差分回波可以直接还原复杂波形
   - 与现有 Wiener 反卷积方法比较

4. **时间分辨率验证**：
   - 使用阶跃信号，细扫 $t_s = 2$ ns
   - 验证上升时间 $\tau \approx (2/\pi)\sqrt{t_\pi^2 + (\pi/2)t_{\text{int}}^2}$
   - 计算传递函数，与 Hann 窗口的傅里叶变换比较

---

## 5. 实现顺序与工作量估计

### Phase 1: 基础实现

1. 在 `pulse.py` 中实现 `create_differential_echo_pulse()` 函数
2. 在 `protocal.py` 的 `case 5` 中实现差分回波协议的演化逻辑
3. 在 `analysis.py` 中实现 `get_signal_from_diff_echo()` 直接反演函数

### Phase 2: 验证与调试

4. 在 Notebook 中编写测试单元，使用简单信号（恒定、正弦）验证协议正确性
5. 调整参数，确保线性近似成立
6. 比较有限脉冲宽度与理想脉冲的差异

### Phase 3: 性能比较

7. 与现有 Ramsey 滑动测量（case 4）进行系统比较
8. 绘制灵敏度 vs $k$ 的曲线
9. 绘制灵敏度 vs $t_{\text{rep}}$ 的曲线，与文献 Fig. 5 比较

### Phase 4: 扩展（可选）

10. 增加退相干效应（$T_2$ 衰减）的模拟
11. 实现相位循环技术，扩展线性范围
12. 增加传递函数计算和反滤波功能

---

## 6. 注意事项

1. **信号重复机制**：仿真中需要让同一个波形在 $2k$ 个 $t_{\text{rep}}$ 窗口内重复出现。可以通过：
   - 构建一个长度为 $2k \cdot t_{\text{rep}}$ 的信号，其中包含 $2k$ 次相同波形的拷贝
   - 或者在 `Signal` 类中增加周期性重复功能

2. **时间轴对齐**：$\pi$ 脉冲必须精确落在目标时间点。由于 `CompositePulse` 使用连续拼接，需要精确计算自由演化段的时长。

3. **计算效率**：总演化时间 $\sim 2k \cdot t_{\text{rep}}$ 可能较长，mesolve 的时间步数需要足够密（至少与 $\pi$ 脉冲分辨率匹配）。建议使用预计算的列表哈密顿量而非回调函数。

4. **相位约定**：确保 $(\pi/2)_X$ 对应 `phase=0`，$(\pi/2)_Y$ 对应 `phase=π/2`，与 `Pulse` 类的相位定义一致。

5. **基线校准**：实际测量中 $p = 0.5$ 对应零信号，但仿真中可能存在数值偏差。建议同时进行零信号基线测量，用 $\Delta p = p_{\text{signal}} - p_{\text{baseline}}$ 进行校准。
