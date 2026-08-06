# 快速开始

运行端到端传感管道最快的方式是高层的 {py:class}`~sqc.workflows.SensingWorkflow`。
它把 qubit 器件、磁通信号、控制序列、QuTiP 仿真和重建算法组合在一个流式接口后面。

## 传感并重建磁通波形

```python
from sqc.workflows import SensingWorkflow

wf = SensingWorkflow().configure(
    protocol="ramsey",         # 传感协议
    signal_type=3,             # 磁通信号形状(高斯脉冲)
    signal_amplitude=0.01,     # 单位为 Phi_0
    reconstruction="unwrap",   # 重建算法(Ramsey 可选 "unwrap" 或 "iq")
)
result = wf.run(measure=True, reconstruct=True)

wf.plot()              # 画出「测量 Δp + 重建波形 B(t)」
# result 同样携带原始测量与重建波形
```

`configure()` 只修改你显式传入的参数,其余保持默认(由
{py:data}`sqc.config.CONFIG` 推导)。它返回 `self`,因此可以链式调用。`run()`
构建实验、调用 QuTiP `mesolve`,并完成重建。

## 比较重建算法

多方法比较适用于 `transient`(瞬态)协议,其基于核卷积的算法家族(`wiener`、
`hammerstein`)从同一份测量重建:

```python
wf = SensingWorkflow().configure(protocol="transient", signal_type=3)
wf.run(measure=True, reconstruct=False)      # 先测量一次
comparison = wf.compare(methods=["wiener", "lm"])
print(comparison.best)
```

每个协议只接受各自的重建方法:`transient` → `wiener` / `hammerstein` / `lm`;
`ramsey` → `unwrap` / `iq`。

## 扫描参数

```python
wf = SensingWorkflow().configure(protocol="ramsey", reconstruction="unwrap")
sweep = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02])
```

## 下一步

- {doc}`architecture`:workflow 所基于的八层栈。
- {doc}`building_blocks/index`:直接使用各层以获得完全控制。
- {doc}`examples/index`:深入三条内置管道。
- {doc}`extending`:构建新的传感应用。
