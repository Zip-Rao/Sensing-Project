# sqc — 超导量子比特磁场传感仿真平台

```{note}
本文档站正在建设中。页面正在逐步填充;下方结构反映规划的“框架优先”布局。
```

**sqc** 是一个基于 QuTiP 的、可组合的全栈仿真框架,用于超导 Transmon 量子比特的
时变磁场传感。它在八层 cQED 栈上建模“磁通 → 比特频率”的传导链路,并从仿真测量中
重建被测波形。

三条内置科研主线(波形重建、频率标定、预失真)是组合各层的**范例**,而非平台能力的
边界。每一层都暴露扩展接口,你可以据此构建自己的传感应用。

```{toctree}
:maxdepth: 2
:caption: 开始使用

overview
install
quickstart
```

```{toctree}
:maxdepth: 2
:caption: 框架

architecture
building_blocks/index
theory
```

```{toctree}
:maxdepth: 2
:caption: 指南

examples/index
extending
```

```{toctree}
:maxdepth: 2
:caption: 完整示例

tutorial
```

```{toctree}
:maxdepth: 2
:caption: 参考

api/index
```

```{toctree}
:maxdepth: 2
:caption: 关于

roadmap
```
