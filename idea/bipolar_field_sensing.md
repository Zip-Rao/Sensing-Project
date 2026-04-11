# 支持正负磁场的瞬态磁场传感方案

## 1. 问题分析

### 1.1 当前协议的限制

在现有的瞬态磁场测量协议（`protocal.py` case 4）中，Hammerstein-Wiener反卷积的非线性逆映射使用了：

```python
B = (1/π) * arccos((ω + ω_01 + E_C)² / (8·E_C·E_J))
```

这来源于Transmon频率与磁通的关系：

$$\omega_{01}(\Phi) = \sqrt{8 E_J E_C |\cos(\pi\Phi)|} - E_C$$

**核心问题：** `|cos(πΦ)|` 中的绝对值使得 `ω(Φ)` 关于 `Φ=0` 是**偶函数**。即：

$$\omega(+B) = \omega(-B)$$

因此，`arccos` 逆映射无法区分 `+B` 和 `-B`，正负磁场产生完全相同的频率响应。

### 1.2 物理图像

在sweet spot（`Φ=0`）处，频率-磁通曲线呈抛物线顶点形状：

```
ω(Φ)
  │    ╱‾‾‾╲
  │   ╱     ╲
  │  ╱       ╲
  │ ╱         ╲
  │╱           ╲
  ┼──────┬──────── Φ
      -Φ₀  0  +Φ₀
```

在sweet spot处 `dω/dΦ = 0`，不仅无法区分正负，而且对磁场的一阶灵敏度为零。

---

## 2. 方案：偏置磁通工作点（Flux-Biased Operating Point）

### 2.1 核心思想

将qubit的静态工作点从sweet spot（`Φ_bias = 0`）偏移到频率-磁通曲线的**斜坡区域**，使得：

- `ω(Φ)` 在工作点附近**单调**，正负磁场产生不同方向的频率偏移
- 一阶灵敏度 `dω/dΦ ≠ 0`，可以线性区分正负信号

```
ω(Φ)
  │    ╱‾‾╲
  │   ╱  ● ╲        ← 偏置工作点 Φ_bias
  │  ╱  ↙ ↘ ╲       ← +B使ω下降，-B使ω上升（或反之）
  │ ╱         ╲
  │╱           ╲
  ┼──────┬──────── Φ
      -Φ₀  0  +Φ₀
```

### 2.2 最优偏置点选择

代码中 `qubit.py` 已实现 `optimal_work_point()` 方法：

```python
def optimal_work_point(self):
    Phi = np.arctan(np.sqrt(2))  # ≈ 0.3041 Φ₀
    return Phi
```

这是频率对磁通灵敏度 `|dω/dΦ|` 最大的点。然而，实际偏置点的选择需要在以下因素间权衡：

| 因素 | 偏置越大 | 偏置越小 |
|------|---------|---------|
| 一阶灵敏度 `dω/dΦ` | 先增后减 | 接近零 |
| 退相干时间 T2 | 显著缩短（对磁通噪声敏感） | 最长（sweet spot保护） |
| 频率-磁通线性度 | 在最优点附近较好 | 二次展开主导 |
| 可测磁场动态范围 | 较小（易超出线性区） | 较大 |

**推荐偏置点：** `Φ_bias ≈ 0.25 ~ 0.35 Φ₀`，具体取值可通过数值优化确定（见第4节）。

### 2.3 偏置后的频率-磁通关系

在偏置点 `Φ_bias` 附近，对小信号 `δΦ` 展开：

$$\omega(\Phi_{bias} + \delta\Phi) \approx \omega_0 + \left.\frac{d\omega}{d\Phi}\right|_{\Phi_{bias}} \cdot \delta\Phi + \frac{1}{2}\left.\frac{d^2\omega}{d\Phi^2}\right|_{\Phi_{bias}} \cdot \delta\Phi^2 + \cdots$$

其中一阶项 `dω/dΦ` 携带了磁场方向信息——这正是sweet spot处所缺失的。

---

## 3. 对现有协议框架的修改方案

### 3.1 修改总览

保持现有协议的整体框架（sliding measurement + kernel extraction + Wiener deconvolution + 非线性逆映射），仅需修改以下模块：

```
┌─────────────────────────────────────────────────────┐
│           现有协议流程（保持不变）                      │
│                                                     │
│  磁场B(t) → 滑动测量 → Δp(t) → Wiener反卷积 → ω(t) │
│                                                     │
│         ↓ 仅修改此处 ↓                               │
│                                                     │
│  ω(t) → [新的非线性逆映射] → B(t)（含正负）           │
└─────────────────────────────────────────────────────┘
```

### 3.2 模块一：Qubit初始化——添加偏置磁通（`qubit.py`）

**修改位置：** `protocal.py` case 4 中创建qubit时

**当前代码：** qubit在sweet spot工作（`flux=0`，默认值）

**修改方案：** 创建qubit时指定偏置磁通

```python
# 在 protocal.py case 4 中
# 原来: qubit默认 flux=0
# 修改为:
Phi_bias = qubit.optimal_work_point()  # ≈ 0.3041 Φ₀
# 或者由用户指定:
# Phi_bias = 0.25  # 可调参数

qubit_biased = TransmonQubit(
    EC=qubit.EC, EJ=qubit.EJ_0, T1=qubit.T1, T2=qubit.T2,
    flux=Phi_bias,  # 关键修改：偏置磁通
    state=qubit.state, n_levels=qubit.n_levels
)
```

**注意：** 偏置后的qubit频率会改变，控制脉冲的驱动频率 `omega_d` 也需要相应更新。

### 3.3 模块二：控制脉冲频率更新（`protocal.py`）

Ramsey脉冲的驱动频率需匹配偏置后的qubit频率：

```python
# 原来：
control_pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

# 修改为：
control_pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit_biased.frequency)
```

### 3.4 模块三：核提取——在偏置点处提取（`analysis.py`）

**当前 `get_kernel()` 方法**已经接受qubit参数，只要传入偏置后的qubit，kernel会自然反映偏置点处的响应特性。

**关键变化：** 在偏置点处，kernel不再关于时间轴中心对称，因为 `dω/dΦ ≠ 0` 引入了奇对称分量。这意味着kernel能够区分正负频率偏移。

```python
# 无需修改 get_kernel 函数本身，只需传入偏置后的qubit
t_samples, kernel = analysis.get_kernel(control_pulse, qubit_biased)
```

但需要注意：当前kernel提取使用的刺激脉冲（`Signal(type=3, amplitude=0.0215)`）是单向的（正振幅），在偏置点处正负刺激的响应不完全对称。对于线性化近似，这在小信号下是可接受的；对于更高精度，可考虑双向刺激（见3.7节扩展）。

### 3.5 模块四：Hammerstein-Wiener非线性逆映射（`analysis.py`）——核心修改

这是最关键的修改。需要将 `arccos` 逆映射替换为支持符号的数值逆映射。

**方案A：解析逆映射（推荐用于小信号）**

在偏置点附近，频率-磁通关系可以精确求逆：

$$\Phi = \frac{1}{\pi}\arccos\left(\frac{(\omega + E_C)^2}{8 E_C E_{J,0}}\right)$$

由于工作在 `Φ_bias > 0` 的斜坡上，`arccos` 返回的值对应**正确的磁通分支**。磁场信号为：

$$B(t) = \Phi(t) - \Phi_{bias}$$

当 `B > 0` 时 `Φ > Φ_bias`，当 `B < 0` 时 `Φ < Φ_bias`，符号自然保留。

```python
def hammerstein_wiener_deconvolution_biased(self, qubit, delta_p, kernel, dt, lambdas):
    """
    支持正负磁场的 Hammerstein-Wiener 反卷积
    要求qubit工作在偏置磁通点（非sweet spot）
    """
    # Step 1: 线性Wiener反卷积，得到频率偏移 δω(t)
    omega_lists, delta_omega = self.wiener_deconvolution(delta_p, kernel, dt, lambdas)

    # Step 2: 将频率偏移转换为绝对频率
    omega_abs = delta_omega + qubit.frequency  # qubit.frequency 已经是偏置点的频率

    # Step 3: 数值求逆——从 ω 反推 Φ
    # ω = sqrt(8 * EJ_0 * |cos(πΦ)| * EC) - EC
    # => |cos(πΦ)| = (ω + EC)² / (8 * EC * EJ_0)
    cos_val = (omega_abs + qubit.EC) ** 2 / (8 * qubit.EC * qubit.EJ_0)
    cos_val = np.clip(cos_val, 0, 1)  # 数值安全

    # Step 4: 选择正确的磁通分支
    # 在偏置点 Φ_bias 附近，cos(πΦ) 的符号和单调性已知
    # 对于 0 < Φ_bias < 0.5，cos(πΦ) > 0 且递减
    # arccos 返回 [0, π]，对应 Φ ∈ [0, 1]
    Phi_total = (1 / np.pi) * np.arccos(cos_val)

    # Step 5: 磁场 = 总磁通 - 偏置磁通
    B = Phi_total - qubit.flux

    return omega_lists, B
```

**方案B：数值查表逆映射（推荐用于大信号或高精度需求）**

当磁场幅度较大，超出线性区时，使用预计算的查找表进行逆映射：

```python
def build_inverse_lookup(self, qubit, Phi_range=(-0.1, 0.1), n_points=10000):
    """
    构建频率→磁通的逆映射查找表
    Phi_range: 相对于偏置点的磁通范围
    """
    Phi_bias = qubit.flux
    Phi_array = np.linspace(Phi_bias + Phi_range[0], Phi_bias + Phi_range[1], n_points)
    omega_array = np.array([
        np.sqrt(8 * qubit.EJ_0 * abs(np.cos(np.pi * phi)) * qubit.EC) - qubit.EC
        for phi in Phi_array
    ])

    # 检查单调性（在偏置点附近应为单调的）
    # 对于 0 < Φ_bias < 0.5 且 Phi_range 不跨越 sweet spot
    # omega 应该是 Φ 的单调递减函数

    # 返回插值函数: omega -> B
    from scipy.interpolate import interp1d
    B_array = Phi_array - Phi_bias
    # 注意：omega_array 可能是递减的，需要翻转
    sort_idx = np.argsort(omega_array)
    interp_func = interp1d(omega_array[sort_idx], B_array[sort_idx],
                           kind='cubic', fill_value='extrapolate')
    return interp_func

def hammerstein_wiener_deconvolution_biased_v2(self, qubit, delta_p, kernel, dt, lambdas):
    """
    使用查找表的正负磁场 Hammerstein-Wiener 反卷积
    """
    omega_lists, delta_omega = self.wiener_deconvolution(delta_p, kernel, dt, lambdas)
    omega_abs = delta_omega + qubit.frequency

    inverse_map = self.build_inverse_lookup(qubit)
    B = inverse_map(omega_abs)

    return omega_lists, B
```

### 3.6 模块五：信号生成——添加双极性信号类型（`signal.py`）

为测试正负磁场，需要添加能产生正负交替信号的类型：

```python
# 建议在 Signal 类中添加新的信号类型（如 type=7）
case 7:  # 双极性瞬态脉冲（正负交替）
    for amp, ctr, r, f in zip(amplitudes, centers, rises, falls):
        # amp 可以为正或负
        signal += amp * np.where(
            t_list < ctr, 0,
            (1 - np.exp(-(t_list - ctr) / r)) * np.exp(-(t_list - ctr) / f)
        )

# 或者更简单：允许 type=4 的 amplitude 为负值
# 当前 type=4 本身并不限制 amplitude 的符号，
# 但可以组合多个 type=4 信号来构建正负交替的磁场
```

也可以直接使用 type=2（正弦信号）作为正负交替的测试信号。

### 3.7 扩展：双刺激Kernel提取（可选，提高精度）

在偏置点处，由于非线性效应，正向刺激和负向刺激的kernel可能略有不同。可扩展kernel提取为双向：

```python
def get_kernel_bipolar(self, control_pulse, qubit):
    """
    双向kernel提取：分别对正、负刺激提取kernel，
    取平均作为线性kernel，取差作为非线性校正项
    """
    kernel_pos = self._extract_kernel(control_pulse, qubit, stim_amplitude=+0.0215)
    kernel_neg = self._extract_kernel(control_pulse, qubit, stim_amplitude=-0.0215)

    kernel_linear = (kernel_pos + kernel_neg) / 2  # 线性部分
    kernel_nonlinear = (kernel_pos - kernel_neg) / 2  # 非线性校正

    return kernel_linear, kernel_nonlinear
```

---

## 4. 偏置点优化策略

### 4.1 灵敏度-退相干权衡

偏置点的选择本质上是灵敏度与相干时间的权衡。可以定义一个综合性能指标：

$$\text{FoM}(\Phi_{bias}) = \left|\frac{d\omega}{d\Phi}\right|_{\Phi_{bias}} \cdot T_2^*(\Phi_{bias})$$

其中 `T_2*(Φ_bias)` 是偏置点处的有效退相干时间（受磁通噪声限制）：

$$\frac{1}{T_2^*} = \frac{1}{T_2} + \pi \left|\frac{d\omega}{d\Phi}\right| \cdot S_\Phi^{1/2}$$

`S_Φ` 为磁通噪声功率谱密度。

### 4.2 数值优化实现

```python
def optimize_bias_point(qubit, flux_noise_amplitude=1e-5):
    """
    数值搜索最优偏置点
    """
    Phi_range = np.linspace(0.05, 0.45, 1000)
    fom = []
    for phi in Phi_range:
        q = TransmonQubit(EC=qubit.EC, EJ=qubit.EJ_0,
                          T1=qubit.T1, T2=qubit.T2, flux=phi)
        sensitivity = abs(q.frequency_sensitivity())
        # 估算 T2*
        T2_star = 1.0 / (1.0/qubit.T2 + np.pi * sensitivity * flux_noise_amplitude)
        fom.append(sensitivity * T2_star)

    optimal_idx = np.argmax(fom)
    return Phi_range[optimal_idx]
```

---

## 5. 完整修改后的协议流程

```python
# protocal.py case 4 修改后的完整流程

case 4:  # 瞬态磁场测量协议（支持正负磁场）
    # 1. 选择偏置工作点
    Phi_bias = qubit.optimal_work_point()  # 或 optimize_bias_point(qubit)
    qubit_biased = TransmonQubit(
        EC=qubit.EC, EJ=qubit.EJ_0, T1=qubit.T1, T2=qubit.T2,
        flux=Phi_bias, state=qubit.state, n_levels=qubit.n_levels
    )

    # 2. 生成测试磁场（可包含正负分量）
    Phi = Signal(type=2, t_list=np.linspace(0, 200, 400),
                 amplitude=0.01, frequency=0.02)  # 正弦信号，含正负
    Phi_0 = Signal(type=1, t_list=np.linspace(0, 200, 400), amplitude=0.0)

    # 3. 创建控制脉冲（使用偏置后的频率）
    t_rabi = np.linspace(0, 10, 20)
    control_pulse = create_ramsey_pulse(t_rabi, tau=0.0,
                                         omega_d=qubit_biased.frequency)

    # 4. 滑动测量（无需修改，自动使用偏置后的qubit）
    scan_list, p_e = self.sliding_measurement(qubit_biased, Phi, control_pulse)
    scan_list_base, p_e_base = self.sliding_measurement(qubit_biased, Phi_0, control_pulse)

    # 5. 核提取（在偏置点处）
    analysis = Analysis()
    t_samples, kernel = analysis.get_kernel(control_pulse, qubit_biased)

    # 6. 反卷积
    delta_p = np.array(p_e) - np.array(p_e_base)

    # 7. Hammerstein-Wiener反卷积（使用新的支持正负的逆映射）
    B_list, B = analysis.hammerstein_wiener_deconvolution_biased(
        qubit_biased, delta_p, kernel, dt, lambdas
    )

    return t_samples, kernel, scan_list, delta_p, p_e, Phi, control_pulse
```

---

## 6. 方案对比与局限性

### 6.1 与sweet spot方案的对比

| 特性 | Sweet Spot (Φ=0) | 偏置工作点 (Φ≈0.3) |
|------|------------------|-------------------|
| 正负磁场区分 | 不可能 | 可以 |
| 一阶灵敏度 | 0 | ~数GHz/Φ₀ |
| T2退相干 | 最长（受保护） | 缩短（对磁通噪声敏感） |
| 频率稳定性 | 最佳 | 需要精确磁通控制 |
| 可测信号动态范围 | 大（二阶响应） | 受限于线性区宽度 |
| 非线性逆映射 | arccos（有歧义） | arccos或查表（无歧义） |

### 6.2 局限性与缓解措施

1. **T2退相干缩短**
   - 缓解：使用自旋回波（spin echo）或CPMG序列消除低频磁通噪声
   - 缓解：选择较小的偏置（如 `Φ_bias ≈ 0.15`），牺牲部分灵敏度换取更长T2

2. **动态范围受限**
   - 当 `|B|` 超过线性区（约 `±0.1 Φ₀`），非线性失真增大
   - 缓解：使用方案B的查找表逆映射处理大信号
   - 缓解：对大信号使用 `numerical_inverse()` 全数值反演

3. **磁通偏置精度要求**
   - 偏置点需要精确设定和稳定维持
   - 缓解：通过校准程序定期测量和调整偏置点

4. **kernel的非对称性**
   - 偏置后kernel不再对称，可能增加反卷积的数值误差
   - 缓解：增加kernel提取的采样点数；使用双向kernel提取

---

## 7. 实施步骤

按优先级排序：

1. **[必需]** 在 `protocal.py` case 4 中添加偏置qubit的创建和使用
2. **[必需]** 在 `analysis.py` 中实现 `hammerstein_wiener_deconvolution_biased()` 方法
3. **[必需]** 更新控制脉冲驱动频率为偏置后的qubit频率
4. **[推荐]** 在 `signal.py` 中添加双极性测试信号（或直接使用type=2正弦）
5. **[推荐]** 实现偏置点数值优化函数
6. **[可选]** 实现双向kernel提取以提高精度
7. **[可选]** 在 `qubit.py` 中添加偏置点处的T2估算方法
