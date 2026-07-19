# sqc — 超导量子比特量子传感仿真平台

[English](README.md) | **中文**

基于 QuTiP 的全栈仿真平台,用于超导 Transmon 量子比特的**时变磁场传感**。它对
磁通 → 量子比特频率的转导链路建模,并从仿真测量中重建被感知的波形,遵循
Gao、Rol、Touzard、Wang(2021)提出的六层 cQED 架构。

- **版本:** 1.0.0 &nbsp;·&nbsp; **许可证:** MIT &nbsp;·&nbsp; **Python:** ≥ 3.10
- 全程自然单位(ħ = 1):频率/能量用 GHz,时间用 ns,磁通用 Φ₀。

## v1 提供的能力

三条端到端科研主线,外加交互式 Web 演示:

1. **波形重建** — 通过滑动瞬态协议感知未知的时变磁通并重建它(Wiener 反卷积、
   Hammerstein–Wiener,或 Levenberg–Marquardt),另含 Ramsey / 差分回波 /
   Cryoscope 协议。
2. **量子比特频率标定** — 基于 Ramsey 的 `f(Φ)` / `f₀₁` 标定。
3. **波形预失真** — 测量控制线失真并设计 IIR/FIR 逆滤波器,附端到端验证流程。

分层 API(`sqc.devices → control → hardware → simulation → experiments →
reconstruction → calibration → workflows`)让你在所需的任意层次工作。

## 安装

需要 Python ≥ 3.10。可编辑安装让 `import sqc` 在任意位置生效:

```bash
pip install -e .            # 核心(numpy, scipy, qutip, matplotlib)
pip install -e ".[demo]"    # + Gradio Web 演示
pip install -e ".[test]"    # + pytest 测试工具
```

或纯依赖安装:`pip install -r requirements.txt`。

## 快速开始

统一入口是 `SensingWorkflow` —— 配置、运行、绘图:

```python
from sqc.workflows import SensingWorkflow

wf = SensingWorkflow()
wf.configure(
    protocol="transient",       # ramsey | echo | transient | cryoscope | delay_ramsey
    signal_amplitude=0.01,      # Φ₀
    reconstruction="wiener",    # wiener | hammerstein | lm
    t_rabi_duration=20,         # ns
)
result = wf.run()               # 测量 + 重建
wf.plot()                       # matplotlib 图
```

扫描参数或对比重建方法:

```python
sweep = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02, 0.05])

wf.run(measure=True, reconstruct=False)
cmp = wf.compare(methods=["wiener", "hammerstein", "lm"])
print("best:", cmp.best)
```

### 底层 API

每一层都可直接使用,例如:

```python
from sqc.devices.transmon import TransmonQubit          # EC, EJ, T1, T2, flux
from sqc.experiments import TransientSensingExperiment
from sqc.reconstruction import KernelEstimator, TransientReconstruction
```

完整模块参考见 [`docs/architecture.md`](docs/architecture.md)。

## Web 演示

交互式 Gradio 界面(器件配置、协议、重建、预失真):

```bash
pip install -e ".[demo]"
python web_demo_v2.py
```

## 文档

- **教程 notebook:** [`Simulation_sqc.ipynb`](Simulation_sqc.ipynb) —— Rabi → Ramsey →
  差分回波 → 瞬态 → Cryoscope,端到端演示。
- **架构 / 模块参考与路线图:** [`docs/architecture.md`](docs/architecture.md)
  (中文,六层栈、逐模块 API、扩展指南、§A2 post-v1 路线图)。

## 测试与开发

测试套件、冻结的旧实现 `src/` 及其 `src_mirror` facade 位于**开发树**中
(不随发行包分发)。在开发检出中:

```bash
pip install -e ".[test]"
pytest tests/ -v
pytest tests/regression -m regression      # 物理回归 baseline
```

## 目录结构(分发)

```
sqc/                  平台本体(devices, control, hardware, simulation,
                      experiments, reconstruction, calibration, workflows)
web_demo_v2.py        Gradio Web 演示(基于 sqc)
Simulation_sqc.ipynb  端到端教程 notebook
docs/                 架构 / 技术文档
```

开发树另含测试套件(`tests/`)、冻结的旧实现(`src/`、`src_mirror/`)及科研
资料 —— 这些不随发行包分发。

## v1 未包含(计划中)

以下能力已存在于代码库中,但在 v1 公开接口中**被隐藏**,计划在后续版本引入:
两比特 Z-crosstalk 提取、基于瞬态的频率标定、CPMG 协议、可调耦合器与电子学
(AWG/ADC)硬件层。路线图见 [`docs/architecture.md`](docs/architecture.md) §A2。

## 引用

若在学术工作中使用本平台,请引用 —— 见 [`CITATION.cff`](CITATION.cff)。

## 许可证

MIT © 2026 Zip。见 [`LICENSE`](LICENSE)。
