# 概览

**sqc** 是一个基于 QuTiP 的、可组合的全栈仿真框架,用于超导 Transmon 量子比特的
时变磁场传感。它在八层 cQED 栈上建模"磁通 → 比特频率"的传导链路,并从仿真测量中
重建被测波形。

## 一个框架,而非三个应用

平台内置三条端到端科研管道:

1. **波形重建** —— 传感未知的时变磁通并将其重建(Wiener 去卷积、
   Hammerstein–Wiener 或 Levenberg–Marquardt),含 Ramsey / 差分回波 / cryoscope
   协议。
2. **频率标定** —— 基于 Ramsey 的 $f(\Phi)$ / $f_{01}$ 标定。
3. **预失真** —— 建模 AWG → 芯片传递函数并设计补偿滤波器。

这些是**组合各层的范例** —— 而非平台能力的边界。每一层(devices、hardware、control、
simulation、experiments、reconstruction、calibration、workflows)都暴露扩展接口。
你可以继承其中任意一层来构建自己的传感应用:新协议、新反演算法、新器件类型。见
{doc}`extending`。

## 关键特性

- **纯 Python + QuTiP** —— 通过 `mesolve` 进行全密度矩阵时间演化。
- **自然单位制**(ħ = 1):rad·GHz、ns、Φ₀。
- **不可变器件参数**与**无副作用**的 Hamiltonian 构造。
- **中心化配置** —— 单一来源 {py:data}`sqc.config.CONFIG`。
- **八层架构**,遵循 Gao 等人(2021)的 cQED 工程标准。

## 下一步

- {doc}`install` —— 让它跑起来。
- {doc}`quickstart` —— 三行传感管道。
- {doc}`architecture` —— 八层栈。
- {doc}`building_blocks/index` —— 各层详解。
```{note}
此页由 README 派生作为初稿,供你按叙事需要修订。
```
