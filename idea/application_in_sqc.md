# 瞬态磁场协议在超导量子计算中的应用分析

## 1. 核心问题：测量传递矩阵 $\hat{\mathbf{H}}(\omega)$

超导量子处理器中，所有基于 flux pulse 的操作（CZ 门、iSWAP 门、频率停泊、参数门等）都依赖于准确控制 qubit 频率的时域轨迹。AWG 输出的理想波形经过控制线路后到达芯片的实际磁通可能与期望值不同。

对于 $N$ 个 qubit 的处理器，所有 Z 线的输入-输出关系统一为：

$$\hat{\mathbf{\Phi}}(\omega) = \hat{\mathbf{H}}(\omega) \hat{\mathbf{V}}(\omega)$$

其中 $\hat{\mathbf{H}}(\omega) \in \mathbb{C}^{N \times N}$ 是频率依赖的传递矩阵：

- **对角元** $H_{ii}(\omega)$：qubit $i$ 自身 Z 线的传递函数（包含线缆失真、滤波器响应、阻抗失配等）
- **非对角元** $H_{ji}(\omega)$，$j \neq i$：从 Z 线 $i$ 到 qubit $j$ 的串扰传递函数

理想系统中 $\hat{\mathbf{H}} = \text{diag}(\alpha_1, \ldots, \alpha_N)$，实际系统中对角元偏离理想值（自身失真），非对角元非零（串扰）。

**预失真、串扰补偿、频率轨迹监测、ring-down 测量**本质上都是同一个问题的不同侧面：**标定 $\hat{\mathbf{H}}(\omega)$ 并设计补偿滤波器 $\hat{\mathbf{H}}^{-1}(\omega)$**。

| 通常名称 | 本质 | 涉及的矩阵元素 |
|---------|------|--------------|
| 预失真 / 脉冲校正 | 补偿自身 Z 线的传递函数 | 对角元 $H_{ii}$ |
| 频率轨迹监测 | 验证 $H_{ii}$ 补偿后的残差 | 对角元 $H_{ii}$ |
| Ring-down 测量 | 测量 $H_{ii}$ 的长尾衰减分量 | 对角元 $H_{ii}$ |
| Z-crosstalk 补偿 | 补偿跨线串扰 | 非对角元 $H_{ji}$ |

全局补偿公式统一为：

$$\hat{\mathbf{V}}_{\text{comp}}(\omega) = \hat{\mathbf{H}}^{-1}(\omega) \hat{\mathbf{\Phi}}_{\text{target}}(\omega)$$

---

## 2. 两种测量协议的分工

### 2.1 Cryoscope：测量对角元 $H_{ii}$

Cryoscope（Rol et al., APL 2020）在 qubit 自身的 Z 线上施加阶跃脉冲，用截断 Ramsey 测量阶跃响应。

**核心优势**：在甜点（$\Phi = 0$, $\kappa = 0$）工作时，截断信号的关断瞬态误差被 $\Delta f_Q \propto \Phi^2$ 的二次非线性自然抑制，精度可达 0.1%。

**适用范围**：对角元 $H_{ii}$。自身 flux pulse 的幅度通常较大（$\Phi \sim 0.1\text{-}0.25\,\Phi_0$），即使在甜点的二次响应下，信号仍然足够强。

### 2.2 瞬态磁场协议：测量非对角元 $H_{ji}$

对 Z 线 $i$ 施加测试脉冲，在 qubit $j$（$j \neq i$）上使用滑动 Ramsey + 反卷积测量串扰波形。

**核心优势**：探测 qubit 不被 flux pulse 直接驱动，可以在最优灵敏度点工作（$\kappa$ 最大，线性响应）。

**适用范围**：非对角元 $H_{ji}$。串扰信号小（$\Phi_{\text{xtalk}} \sim 0.001\text{-}0.01\,\Phi_0$），需要高灵敏度的线性响应。

### 2.3 为什么不能反过来用

**Cryoscope 测非对角元的困难**：

Cryoscope 的精度依赖甜点工作。但在甜点处 $\kappa = 0$，频率对磁通是二次响应：

$$\delta\omega_{\text{sweet}} = \alpha \Phi^2$$

对于串扰信号 $\Phi_{\text{xtalk}} \sim 0.003\,\Phi_0$：

$$\delta\omega_{\text{sweet}} \sim \alpha \times (0.003)^2 = \alpha \times 9 \times 10^{-6}\,\Phi_0^2$$

二次响应把小串扰信号进一步压缩了约 300 倍。若改在非甜点工作以获得线性灵敏度，则 Cryoscope 的截断瞬态误差抑制优势消失。

| Cryoscope 工作点 | 对串扰的灵敏度 | 截断瞬态误差 | 综合 |
|---|---|---|---|
| 甜点（$\kappa = 0$） | 极低（$\propto \Phi^2$） | 被抑制 | 信号太弱 |
| 非甜点（$\kappa \neq 0$） | 正常（$\propto \Phi$） | 不被抑制 | 失去核心优势 |

**瞬态协议测对角元的困难**：

测量自身 Z 线的传递函数时，flux pulse 和探测脉冲作用在同一个 qubit 上。flux pulse 导致 qubit 频率大幅变化（如 CZ 门中 $\sim 200\text{-}500$ MHz），探测脉冲失谐，核函数变得**信号依赖**。

- 小信号情况（频率停泊等，$\delta\omega \lesssim 5$ MHz）：核函数近似不变，Wiener 反卷积可用
- 大信号情况（CZ 门等，$\delta\omega \sim 200\text{-}500$ MHz）：需要完整密度矩阵模拟或分段线性化近似，精度和计算量都不如 Cryoscope

### 2.4 分工总结

| | 对角元 $H_{ii}$（自身失真） | 非对角元 $H_{ji}$（串扰） |
|---|---|---|
| **最优方法** | **Cryoscope** | **瞬态磁场协议** |
| 原因 | 甜点 + 大信号 → 高精度 | 最优灵敏度点 + 小信号 → 高灵敏度 |
| 替代方法 | 瞬态协议（仅限小信号） | Cryoscope（仅限非甜点，但精度降低） |

两者互补，共同标定完整传递矩阵 $\hat{\mathbf{H}}(\omega)$。

---

## 3. 非对角元 $H_{ji}$ 的详细理论

这是瞬态磁场协议相对于 Cryoscope 有不可替代价值的部分。

### 3.1 串扰的物理来源

Z 线之间的串扰来源于多种寄生耦合路径：

**（1）互感耦合（dominant）**

Z 控制线之间通过空间电磁耦合产生互感 $M_{ij}$。当 Z 线 $i$ 中流过电流 $I_i(t)$ 时，在 qubit $j$ 的 SQUID 环路中感应出磁通 $\Phi_{j \leftarrow i} = M_{ij} I_i(t)$。对于典型多层布线结构，最近邻串扰系数 $M_{ij}/L_j \sim 1\%\text{-}5\%$。此路径的频率响应较平坦（宽带），串扰波形近似与源脉冲同步。

**（2）地平面回流**

flux pulse 的返回电流通过公共地平面流动，在其他 qubit 的 SQUID 环路中产生附加磁通。频率特性取决于地平面的 sheet impedance 和接地通孔分布，通常表现为低通特性。

**（3）键合线/封装耦合**

键合线之间的互感、PCB 走线串扰。在高频分量上可能显著。

**（4）衬底涡流**

磁通脉冲在衬底中感应涡流，产生长程慢衰减的串扰（$\mu$s 量级）。表现为 flux pulse 结束后的缓慢 ring-down。

总串扰为各路径贡献的叠加：

$$\Phi_{\text{xtalk}}^{(j \leftarrow i)}(t) = \underbrace{M_{ij} I_i(t)}_{\text{瞬时互感}} + \underbrace{\int_0^t h_{\text{ground}}(t-t') I_i(t') dt'}_{\text{地平面}} + \underbrace{\int_0^t h_{\text{eddy}}(t-t') I_i(t') dt'}_{\text{衬底涡流}} + \cdots$$

总传递函数 $H_{ji}(\omega) = M_{ij} + \hat{h}_{\text{ground}}(\omega) + \hat{h}_{\text{eddy}}(\omega) + \cdots$ 是频率依赖的，其零频极限 $H_{ji}(0)$ 即为 DC Ramsey 测量得到的静态串扰系数。

### 3.2 为什么静态补偿不够

静态串扰矩阵 $\mathbf{M} = \hat{\mathbf{H}}(0)$ 仅描述 DC 行为。对于 ns 量级的 flux pulse（带宽 $\sim$ GHz），高频分量的串扰系数可能与 DC 值显著不同：

- 互感路径是宽带的：$H_{ji}^{\text{mutual}}(\omega) \approx \text{const}$
- 地平面路径是低通的：$|H_{ji}^{\text{ground}}(\omega)|$ 随频率衰减
- 衬底涡流路径是极低频的：仅在 $\omega \to 0$ 时有贡献

因此同一对 qubit 之间的串扰，在 DC 和高频处可能幅度和相位都不同。仅用静态矩阵补偿会残留动态误差。

### 3.3 串扰对门保真度的影响

对 qubit A 施加 CZ 门 flux pulse $V_A(t)$，qubit B 受到的串扰磁通：

$$\Phi_B^{\text{xtalk}}(t) = h_{BA} * V_A(t)$$

qubit B 积累的寄生相位：

$$\phi_B^{\text{parasitic}} = \int_0^{T_{\text{gate}}} \kappa_B \cdot \Phi_B^{\text{xtalk}}(t) \, dt$$

对于 CZ 门（$T_{\text{gate}} \sim 40$ ns）、串扰 $\sim 3\%$、$\Phi_A \sim 0.1\,\Phi_0$、$\kappa_B \sim 10$ GHz/$\Phi_0$：

$$\phi_B \sim 10 \times 0.003 \times 40 = 1.2 \text{ rad}$$

这远超容错阈值（$\sim 10^{-3}$ rad），说明**不补偿的串扰会灾难性地破坏门保真度**。静态补偿可以消除 DC 分量，但残留的动态分量（过冲、ring-down 等）仍可能达到 $\sim 0.01\text{-}0.1$ rad 量级。

### 3.4 现有方法的局限

| 方法 | 得到什么 | 得不到什么 |
|------|---------|----------|
| DC Ramsey | 静态串扰系数 $M_{ji} = H_{ji}(0)$ | 频率依赖性、时域波形细节 |
| Chevron 实验 | 定性：有/无串扰 | 定量波形、传递函数 |
| Cryoscope | 自身传递函数 $h_{ii}(t)$ | 跨线传递函数 $h_{ji}(t)$（甜点灵敏度两难） |

**核心盲区**：$h_{ji}(t)$ 的时域波形和 $H_{ji}(\omega)$ 的频率依赖性。

### 3.5 瞬态磁场协议测量串扰

#### 测量方案

```
Qubit A Z 线:   ──────┤ V_A(t) ├──────────
                       │
                       │ h_BA (串扰传递函数)
                       ▼
Qubit B SQUID:  ──── Φ_xtalk(t) = h_BA * V_A(t) ────
                       │
                       │ κ_B (磁通灵敏度)
                       ▼
Qubit B 频率:   ──── δω_B(t) = κ_B · Φ_xtalk(t) ────
                       │
                       │ k(t) (探测核函数)
                       ▼
测量结果:       ──── δp_e(t_d) = [k * δω_B](t_d) ────
```

探测 qubit B 工作在最优灵敏度点（$\kappa_B$ 最大），不施加 flux pulse，核函数 $k(t)$ 良好定义且不依赖信号。

#### 测量模型

$$\delta p_e^{(B)}(t_d) = \kappa_B \int k(t' - t_d) \cdot \Phi_B^{\text{xtalk}}(t') \, dt' = \kappa_B \left[k * \Phi_B^{\text{xtalk}}\right](t_d)$$

#### 反卷积还原串扰波形

$$\Phi_B^{\text{xtalk}}(t) = \frac{1}{\kappa_B} \mathcal{F}^{-1}\left[\frac{\hat{K}^*(\omega)}{|\hat{K}(\omega)|^2 + \lambda^2} \delta\hat{P}_e(\omega)\right]$$

#### 提取传递函数

若 $V_A(t)$ 已知：

$$\hat{H}_{BA}(\omega) = \frac{\hat{\Phi}_B^{\text{xtalk}}(\omega)}{\hat{V}_A(\omega)}$$

使用阶跃信号 $V_A(t) = V_0 u(t)$ 可直接测量串扰阶跃响应 $s_{BA}(t) = \Phi_B^{\text{xtalk}}(t) / V_0$。

### 3.6 核函数选择

由于串扰信号可能跨越多个时间尺度（ns 级互感 + μs 级涡流），推荐多核函数联合策略：

| 串扰分量 | 特征时间 | 推荐核函数 | 理由 |
|----------|---------|-----------|------|
| 瞬时互感串扰 | $\sim 1\text{-}40$ ns | Ramsey $\tau = 0$ | 最高时间分辨率 |
| 地平面回流 | $\sim 30\text{-}500$ ns | Ramsey $\tau = 50$ ns | 兼顾灵敏度和分辨率 |
| 衬底涡流 | $\sim 0.5\text{-}10\ \mu$s | Ramsey $\tau = 500$ ns | 灵敏度优先 |
| 有 $1/f$ 噪声 | — | Echo 核 | 抑制低频噪声 |

多核函数测量后在频域拼接，得到宽频段的 $\hat{H}_{ji}(\omega)$。

### 3.7 实验可行性

**信号幅度**：串扰 3%、$\Phi_A = 0.1\,\Phi_0$、$\kappa_B = 10$ GHz/$\Phi_0$、Ramsey $\tau=0$ 核面积 $\sim 2/\Omega$：

$$\delta p_e \sim 10 \times 0.003 \times 20 \times 10^{-3} = 0.6$$

完全可测。0.1% 串扰也可通过 $10^6$ 次平均达到 $\sim 10^{-3}$ 精度。

**时间分辨率**：Rabi $\Omega/2\pi = 50$ MHz → $t_{\min} \approx 12$ ns（反卷积后 $\sim 5$ ns），足以分辨 40 ns CZ 门的波形结构。

**测量时间**：200 个延迟点 × $10^6$ 次平均 × $5\,\mu$s 周期 $\approx 15$ min，与 Cryoscope 相当。

---

## 4. 完整标定与补偿流程

| 步骤 | 内容 | 方法 | 矩阵元素 |
|------|------|------|---------|
| 1 | 标定对角元 $H_{ii}$ | Cryoscope（甜点） | $N$ 个 |
| 2 | 标定非对角元 $H_{ji}$ | 瞬态磁场协议（最优灵敏度点） | $\sim 2N\text{-}4N$ 个（近邻） |
| 3 | 组装传递矩阵 $\hat{\mathbf{H}}(\omega)$ | 数据处理 | — |
| 4 | 设计补偿滤波器 $\hat{\mathbf{H}}^{-1}(\omega)$ | 逐频率矩阵求逆 | — |
| 5 | 实现为 AWG 实时滤波 | FIR/IIR 参数化 | — |
| 6 | 验证 | 重新测量残余失真/串扰 | — |

---

## 5. 其他潜在应用

以下应用在物理上也归结为测量传递矩阵的某个侧面，但有其独立的实验意义。

### 5.1 可调耦合器动力学

可调耦合器的 flux pulse 通过色散耦合引起邻居 qubit 的频率移动 $\delta\omega_{\text{disp}}(t) \propto g^2(t)/\Delta$。瞬态协议可以测量 $\delta\omega(t)$，间接推断 $g(t)$。

局限：$g^2/\Delta$ 是二次依赖，对小 $g$ 不灵敏。

### 5.2 频率标定（$f_{01}$ 和 $f(\Phi)$）

静态标定不是该协议的强项（标准 Ramsey 灵敏度更优）。但施加已知斜坡信号 + 反卷积可以一次提取非线性系数 $\kappa, \kappa', \kappa''$（见 `frequency_calibration_via_transient.md`）。

---

## 6. 总结

### 一句话定位

> **该协议在超导量子计算中的核心价值是标定传递矩阵 $\hat{\mathbf{H}}(\omega)$ 的非对角元（Z-crosstalk 传递函数），与 Cryoscope（标定对角元）互补，共同实现完整的多 qubit flux 控制标定。**

### 方法对比

| | Cryoscope | 瞬态磁场协议 | DC Ramsey |
|---|---|---|---|
| **测量内容** | $H_{ii}(\omega)$（对角元） | $H_{ji}(\omega)$（非对角元） | $H_{ji}(0)$（零频） |
| **信息维度** | 完整传递函数 | 完整传递函数 | 仅 DC |
| **最优工作点** | 甜点 | 最优灵敏度点 | 任意 |
| **大信号** | ✓ 自然处理 | 需近似或数值模拟 | ✓ |
| **小信号** | ✗ 甜点处灵敏度极低 | ✓ 线性响应 | ✓ |
| **核函数灵活性** | 无 | Echo/CPMG/多 $\tau$ | 无 |
| **精度** | 0.1% | 受反卷积正则化限制 | 受拟合精度限制 |

---

## 7. 对仿真平台的启示

最有说服力的 demo：**模拟双 qubit 系统的 Z-crosstalk 标定**。

1. 建立双 qubit 模型（`Coupled_System` + Z-crosstalk 通道）
2. 对 qubit A 施加 CZ 门 flux pulse，在 qubit B 上用瞬态协议测量串扰波形
3. 展示 DC Ramsey 只给出平均串扰，瞬态协议还原完整时域波形
4. 设计动态补偿脉冲，验证补偿后串扰被消除
5. 对比补偿前/后 qubit B 的寄生相位
