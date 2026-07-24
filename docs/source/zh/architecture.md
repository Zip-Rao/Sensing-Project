# 架构

`sqc` 按**八层全栈 cQED 架构**组织,遵循 Gao、Rol、Touzard 和 Wang(2021,
*PRX Quantum* 2, 040202)提出的六层栈 —— 这是电路 QED 工程的事实标准。每一层都
对应超导量子比特系统中真实的物理组件或控制职责。

## 八层

从底层(物理器件)到顶层(科研流程):

```
┌───────────────────────────────────────────────────────────────┐
│ workflows/      顶层科研管道                                    │  组合 ↓
├───────────────────────────────────────────────────────────────┤
│ calibration/    标定工作流 → CalibrationTable                   │
├───────────────────────────────────────────────────────────────┤
│ reconstruction/ 波形重建算法                                    │
├───────────────────────────────────────────────────────────────┤
│ experiments/    实验协议(器件+磁通+序列+读出)                │
├───────────────────────────────────────────────────────────────┤
│ simulation/     QuTiP 调用:Hamiltonian、mesolve、噪声          │
├───────────────────────────────────────────────────────────────┤
│ control/        控制脉冲:Waveform、FluxSignal、Pulse           │
├───────────────────────────────────────────────────────────────┤
│ hardware/       控制电子学/链路:失真、传递函数                 │
├───────────────────────────────────────────────────────────────┤
│ devices/        物理器件:Transmon、Resonator、Chip             │
└───────────────────────────────────────────────────────────────┘
```

**依赖规则。** 每一层只能依赖其下方的层。禁止反向依赖 —— `reconstruction`
绝不可 import `workflows`。这使每一层都可独立使用与测试。

每一层都暴露一个**抽象基类**作为扩展点({py:class}`~sqc.devices.Device`、
{py:class}`~sqc.experiments.Experiment`、{py:class}`~sqc.reconstruction.Reconstruction`
等)。继承它即可添加你自己的器件、协议或算法 —— 见 {doc}`extending`。

## 与 Gao 2021 的对应

| Gao 2021 层 | 涵盖 | `sqc/` 模块 |
|---|---|---|
| 量子算法 | 编译 | *(仅留接口)* |
| 控制软件 | 脉冲标定、序列 | `experiments/`、`calibration/`、`workflows/` |
| 控制电子学 | AWG、ADC、FPGA | `hardware/electronics.py` |
| 微波信号处理 | IQ mixer、LO、HEMT、滤波器 | `hardware/distortion.py`、`hardware/readout.py` |
| 制冷 + 互连 | 控制线、衰减器 | `hardware/control_line.py`、`hardware/transfer_matrix.py` |
| 器件 | Transmon、resonator、SQUID | `devices/` |

平台聚焦于工程循环的**软件仿真侧**(Hamiltonian 设计 → 仿真 → 表征 → 反馈),
**不涉及**芯片设计、加工或制冷。

## 设计不变量

- **自然单位制**(ħ = 1):频率/能量单位 rad·GHz,时间单位 ns,磁通单位 Φ₀。
- **不可变器件参数**:{py:class}`~sqc.devices.QubitSpec` 是 frozen dataclass ——
  物理参数不会被实验状态污染。
- **无副作用 Hamiltonian 构造**:{py:class}`~sqc.simulation.HamiltonianBuilder`
  是纯函数。
- **中心化配置**:所有时间网格、AWG 参数与 qubit 默认值都从单一来源
  {py:data}`sqc.config.CONFIG` 推导。
