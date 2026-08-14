# 路线图

v1 有意提供一个聚焦、稳定的公共 API:三条科研管道(波形重建、频率标定、预失真)
加上 web demo。下列能力在代码库中处于不同成熟度,但**不属于 v1 公共接口**。它们规划
于后续版本,在此之前对公共 API、前端与本参考隐藏。

## 规划中的能力

- **Z 串扰重建** —— 多比特磁通串扰的表征与补偿。在提升之前,正在解决瞬态重建核中的
  一个数值缩放问题。
- **瞬态磁通响应标定** —— 从未知瞬态信号拟合 $\Delta\omega(\Phi)$ 色散曲线,即
  `FluxResponseCalibration(method="transient")`。仅此一个方法抛
  `NotImplementedError`。瞬态单点频率**测量**
  (`FrequencyMeasurement`、`SinglePointFrequencyCalibration`,以及
  `FrequencyCalibrationWorkflow`)已实现,属于 v1 公开 API。
  `FrequencyCalibrationWorkflow` 支持通过 `stages=[CalibrationStage(...), ...]`
  编排任意多段流水线(每段独立选测量法、步进器、容差与测量时间轴),瞬态→Ramsey
  hybrid 为默认 2 段预设。gradient 段还支持 `advance_when` 谓词,可按残差阈值以外的
  运行时信号自适应切换阶段。
- **CPMG 协议** —— Carr–Purcell–Meiboom–Gill 动力学解耦传感。
- **可调耦合器** —— 用于双比特门方案的 `TunableCoupler` 器件。
- **电子学层** —— 显式的 AWG/ADC/LO 抽象(`hardware/electronics.py`)。
- **SensingWorkflow 规划方法** —— `save`、`load`、`diff`、`benchmark`、`pipeline`、
  `multi_qubit`、`crosstalk`、`find_optimal_work_point`、`detectability_limit`、
  `noise_characterize`、`cross_validate`。它们目前抛出 `NotImplementedError`,并从
  API 参考中略去。

## 已知限制(v1)

- 对 chip/resonator 器件,`collapse_operators()` 返回空列表,在这些配置下静默丢弃
  耗散。
- 短网格上的瞬态波形重建只恢复部分幅度;定量结果请使用文档给出的网格。

## 基础设施路线图

面向完全社区化发行的成熟度工作:公开测试套件 + CI、PyPI 分发、GitHub Release 标签、
README 徽章、贡献指南,以及可引用的 Zenodo DOI。

```{note}
有两项瞬态能力不要混淆。瞬态波形*重建*
({py:class}`~sqc.experiments.TransientSensingExperiment`)与瞬态单点频率*测量*都是
v1 核心特性,已完全支持。仅上面那条瞬态*磁通响应标定*(从未知信号拟合
$\Delta\omega(\Phi)$ 曲线)被推迟。
```
