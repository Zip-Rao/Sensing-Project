# 扩展平台

三条内置管道是示例,不是边界。栈的每一层都暴露一个抽象基类;继承它,你的对象就能
与平台其余部分组合。本页给出两个最常见的配方。

## 添加新的实验协议

场景:一个新的标准实验(ALLXY、T1、CPMG),或自定义磁通传感协议。

1. 新建 `sqc/experiments/my_protocol.py`。
2. 继承 {py:class}`~sqc.experiments.Experiment`,实现 `build_sequence()` 与 `run()`:

```python
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

from sqc.config import CONFIG
from sqc.experiments.base import Experiment
from sqc.control.flux_signal import FluxSignal
from sqc.simulation.result import ExperimentResult


@dataclass
class MyProtocolExperiment(Experiment):
    """我的新传感协议。

    Parameters
    ----------
    qubit : TransmonQubit
    flux_signal : FluxSignal or None
    param1 : float
    """
    qubit: object
    flux_signal: FluxSignal | None = None
    param1: float = 1.0
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())

    def build_sequence(self):
        ...  # 返回一个 PulseSequence

    def run(self) -> ExperimentResult:
        ...  # 运行 mesolve,返回 ExperimentResult(data=..., axes=..., metadata=...)
```

{py:class}`~sqc.simulation.ExperimentResult` 携带 `data`、`axes`、`metadata`、
`config` 字典 —— 这是每个实验返回的统一契约。

## 添加新的重建算法

场景:一个新的反演算法(深度学习反演、压缩感知、Volterra 核展开)。

1. 新建 `sqc/reconstruction/my_algo.py`。
2. 继承 {py:class}`~sqc.reconstruction.Reconstruction`:

```python
from __future__ import annotations
import numpy as np
from sqc.reconstruction.base import Reconstruction


class MyAlgoReconstruction(Reconstruction):
    """用我的算法从测量重建 Phi(t)。"""

    def reconstruct(self, measurement, **kwargs) -> np.ndarray:
        p_e = measurement.data["p_e"]
        # ... 反演 p_e -> phi -> Phi(t) ...
        return reconstructed_flux
```

注册后,即可通过
{py:meth}`SensingWorkflow.configure(reconstruction=...) <sqc.workflows.SensingWorkflow.configure>`
按名字选择,并用 {py:meth}`~sqc.workflows.SensingWorkflow.compare` 与内置方法比较。

## 其他扩展点

| 层 | 基类 | 添加… |
|---|---|---|
| devices | {py:class}`~sqc.devices.Device` | 新 qubit 类型(如 fluxonium) |
| hardware | {py:class}`~sqc.hardware.DistortionModel` | 新的链路失真模型 |
| control | {py:class}`~sqc.control.PulseBase` | 新的脉冲形状 |
| simulation | {py:class}`~sqc.simulation.RunnerBase` | 新的求解策略 |
| calibration | {py:class}`~sqc.calibration.Calibration` | 新的标定流程 |
| workflows | {py:class}`~sqc.workflows.Workflow` | 新的端到端管道 |

遵循依赖规则({doc}`architecture`):你的子类可以 import 其下方的层,绝不可 import
上方的层。
