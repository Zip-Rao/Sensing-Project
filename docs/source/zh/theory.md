# 理论

本页概述三条 v1 管道背后的物理。这是一份自洽的概览;推导详见 Gao、Rol、Touzard
和 Wang(2021,*PRX Quantum* 2, 040202)与 Koch 等人(2007)。

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
