# 范例

这三条管道随平台附带,是把各**基础构件**组合成完整感知应用的**示例**。每篇都遵循
同一套四段结构——目标、物理原理、端到端代码、结果解读——且每段代码都可原样运行。

- {doc}`waveform_reconstruction` —— qubit 当传感器:从激发态布居恢复未知瞬态磁通
  $\Phi(t)$,并用一个 `SensingWorkflow` 入口对比多种重建算法。
- {doc}`frequency_calibration` —— 表征器件本身:测 $f_{01}(\Phi)$ 响应曲线、查目标
  频率对应的偏置、并闭环整定到该频率。
- {doc}`predistortion` —— 校正控制线:测其传递函数、设计逆滤波器,验证片上波形是否
  落回目标。

```{toctree}
:maxdepth: 1

waveform_reconstruction
frequency_calibration
predistortion
```
