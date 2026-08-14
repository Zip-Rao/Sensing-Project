# 理论

本页概述三条 v1 管道背后的物理:从 Transmon 器件模型,到磁通传感的相位积累,再到
**短 Ramsey 型脉冲**及其**响应核(response kernel)**,后者把测量结果与待测磁通波形
联系起来。本页只保留理解 v1 管道所需的最小推导;完整的多阶 Volterra 展开与核提取方法
见项目结题报告。

## Transmon 量子比特

Transmon 是一种基于约瑟夫森结的超导量子比特,其哈密顿量为

$$H = 4 E_C\, n^2 - E_J(\Phi)\, \cos\varphi,$$

其中 $E_C = e^2/2C_\Sigma$ 是充电能,$n$ 与 $\varphi$ 是 Cooper 对数与相位算符,
约瑟夫森能通过 SQUID 环上的外部磁通调制:

$$E_J(\Phi) = E_{J0}\,\lvert\cos(\pi \Phi/\Phi_0)\rvert.$$

在 Transmon 区域 $E_J/E_C \gg 1$,系统是弱非谐振子:

$$\omega_T(\Phi) = \sqrt{8 E_J(\Phi) E_C} - E_C, \qquad \alpha = -E_C,$$

其中 $\omega_T$ 是 $|0\rangle\!\to\!|1\rangle$ 跃迁频率,$\alpha$ 是使二能级子空间
可寻址的非谐性。

## 磁通传感

外部磁通 $\Phi(t)$ 改变 $E_J$,进而改变 qubit 频率 $\omega_T(\Phi)$。Ramsey 序列
$(\pi/2)\!-\!\tau\!-\!(\pi/2)$ 累积相干相位

$$\phi(\tau) = \int_0^\tau \Delta\omega(t')\,dt' = \int_0^\tau \kappa\,\Phi(t')\,dt',$$

其中 $\kappa = d\omega_T/d\Phi$ 是磁通灵敏度。激发态布居为

$$p_e(\tau) = \tfrac{1}{2}\bigl(1 - \cos[\phi(\tau)]\bigr),$$

由此反演 $\phi(\tau)$ —— 从而恢复 $\Phi(t)$。

### 甜点(sweet spot)

由于 $E_J(\Phi) \propto \lvert\cos(\pi\Phi/\Phi_0)\rvert$,频率在整数磁通量子处
平坦($d\omega_T/d\Phi = 0$)。这个*甜点*对磁通噪声一阶不敏感,但传感灵敏度为零;
传感协议会偏置到远离甜点、$\kappa$ 有限之处。

## 控制脉冲

单比特门使用共振微波驱动。在旋转坐标系 + 旋波近似下,

$$H_d \approx \tfrac{\Omega(t)}{2}\bigl(a\,e^{i\phi} + a^\dagger e^{-i\phi}\bigr),$$

旋转角 $\theta = \int \Omega(t)\,dt$。$\pi$ 脉冲翻转布居;$\pi/2$ 脉冲制备叠加态。
DRAG 整形(Motzoi 2009)抑制向第三能级的泄漏。

## 短 Ramsey 型脉冲

上一节的 Ramsey 相位公式默认自由演化时间 $\tau$ 远大于脉冲时长,脉冲期间的相位
积累可忽略。分辨快变磁通则需要缩短序列。受量子速度极限(quantum speed limit)启发的
一类**短 Ramsey 型脉冲**将两段相位正交的控制脉冲前后相接,不留独立的自由演化窗口
(Herb 和 Degen 2024)。两段脉冲可取一般转角$\alpha$,即
$R_y(\alpha)$和$R_{\pm x}(\alpha)$,而不要求固定为$\pi/2$。固定最大Rabi频率时，
减小$\alpha$可缩短控制时间并收窄响应核，但同时减弱布居响应。因此，转角用于权衡
时间分辨率、灵敏度和测量成本。该协议把动态信号的时间分辨率降到脉冲时长量级；Herb 等(2025)
在 NV 中心平台上实现了 1.1 ns。Transmon 的微波控制可在纳秒尺度上编程,适合这一区间。

此时脉冲期间的相位不可忽略,失谐 $\delta\omega(t)$ 与 Rabi 包络 $\Omega(t)$ 都随时间
变化,固定转轴的 Bloch 旋转图像不再适用,需要更一般的响应描述。

## 响应核

把旋转坐标系下的哈密顿量拆成**已知控制项**与**待测微扰项**:

$$H(t) = H_0(t) + V(t), \qquad H_0(t)=\tfrac{\Omega(t)}{2}\bigl(\cos\phi_d\,\sigma_x+\sin\phi_d\,\sigma_y\bigr), \quad V(t)=\tfrac{\delta\omega(t)}{2}\,\sigma_z,$$

其中 $H_0$ 是已知的控制序列,$V(t)$ 由待测失谐 $\delta\omega(t)$ 驱动。在 $H_0$ 定义的
相互作用绘景下,末态测量期望 $\langle M\rangle_T$ 是 $\delta\omega(t)$ 的泛函。对小信号
作一阶展开,测量相对本底的变化是一层**卷积**:

$$\Delta\langle M\rangle \;\approx\; \int_0^T k_1(t)\,\delta\omega(t)\,dt, \qquad k_1(t) = i\,\langle 0|\,[\,W(t),\,M_I(T)\,]\,|0\rangle,$$

其中 $W(t)=\tfrac12 U_0^\dagger(t)\,\sigma_z\,U_0(t)$,$U_0$ 是控制哈密顿量的传播子。
**一阶响应核** $k_1(t)$ 只由初态、控制序列 $H_0$ 与测量算符 $M$ 决定,是传感协议对
失谐的时间权重函数:其幅值与符号给出各时刻失谐对最终测量的贡献强度与方向。

在磁通传感中,小信号下失谐与归一化磁通扰动近似线性,$\delta\omega(t)\approx\kappa\,
\delta\Phi(t)/\Phi_0$($\kappa=d\omega_T/d\Phi$),于是测量与磁通波形之间同样是卷积
$\Delta p(t)\approx (k*\Phi)(t)$。该线性关系是
{doc}`波形重建 <examples/waveform_reconstruction>`反问题的基础。

信号幅值较大时,响应超出线性近似,展开延伸到高阶 Volterra 核 $k_n$。平台已实现这些核:
{py:class}`~sqc.reconstruction.KernelEstimator` 可估计任意阶核,含完整非对角
$k_n(t_1,\dots,t_n)$;瞬态频率标定默认使用三阶核 $\iiint k_3$ 修正主导非线性。
Volterra 展开的完整推导,以及频率核与磁通核的区分,见项目结题报告。

## 推荐参数范围

| 参数 | 范围 | 默认值 |
|---|---|---|
| $E_J/h$ | 10–25 GHz | 15 GHz |
| $E_C/h$ | 160–400 MHz | 200 MHz |
| $E_J/E_C$ | ~50 | ~75 |
| $f_{01}$ | 4–8 GHz | ~4.7 GHz |
| $\alpha/h$ | 200–300 MHz | 200 MHz |

这些使 $f_{01}$ 落在商用 AWG/HEMT 带宽内,且非谐性足够大,可在 ~10–20 ns 门时长下
避免泄漏。

## 参考文献

- Z. Gao, M. Rol, S. Touzard, C. Wang, *Practical Guide for Building
  Superconducting Quantum Devices*, PRX Quantum **2**, 040202 (2021).
- J. Koch et al., *Charge-insensitive qubit design derived from the Cooper pair
  box*, Phys. Rev. A **76**, 042319 (2007).
- F. Motzoi et al., *Simple Pulses for Elimination of Leakage in Weakly
  Nonlinear Qubits*, Phys. Rev. Lett. **103**, 110501 (2009).
- K. Herb, C. L. Degen, *Quantum speed limit in quantum sensing*, Phys. Rev.
  Lett. **133**, 210802 (2024).
- K. Herb et al., *Quantum magnetometry of transient signals with a time
  resolution of 1.1 nanoseconds*, Nat. Commun. **16**, 822 (2025).
