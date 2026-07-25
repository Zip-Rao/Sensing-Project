# 波形重建

## 目标

给定一段作用在 qubit 上、但**形状未知**的瞬态磁通 $\Phi(t)$,只通过测量 qubit 的
激发态布居,把这段波形**反演**出来。这是平台的第一条产品主线:qubit 当传感器,
外界磁通当被测信号,用 {py:class}`~sqc.workflows.SensingWorkflow` 一个入口跑完
「配置 → 测量 → 重建」整条管道。

## 物理原理

外部磁通 $\Phi(t)$ 通过 SQUID 环调制约瑟夫森能,进而移动 qubit 频率
$\omega_T(\Phi)$。把一个 Ramsey 型控制脉冲**滑过**这段磁通信号:在每个延迟处,
qubit 累积的相干相位正比于该时刻附近磁通对频率的推动,于是激发态布居 $p_e$
随滑动延迟画出磁通信号的一个「模糊化」的像。这个模糊由**控制核**
$k(t)$ 描述——它是脉冲对瞬时磁通的响应函数(见 {doc}`../theory` 的磁通传感一节)。

于是测量与波形之间是一层卷积:

$$\Delta p(t) \;\approx\; (k * \Phi)(t),$$

其中 $\Delta p$ 是「有信号」与「零磁通参考」两次滑动测量之差,已扣掉与磁通无关的
本底。**重建就是这层卷积的反问题**:已知 $\Delta p$ 与核 $k$,解出 $\Phi(t)$。
{doc}`../building_blocks/reconstruction` 层提供多种解法——线性 Wiener 反卷积、
考虑 $\omega(\Phi)$ 非线性的 Hammerstein-Wiener、以及全密度矩阵的 Levenberg-Marquardt
反演——本例用 `SensingWorkflow` 把它们串起来并做 A/B 对比。

## 端到端代码

```python
from sqc.workflows import SensingWorkflow

# ── 1. 配置一次感知实验 ───────────────────────────────────────────────
wf = SensingWorkflow()
wf.configure(
    protocol="transient",     # 瞬态场感知(滑动 Ramsey)
    flux_bias=0.1,            # 偏置到远离甜点、κ 有限处
    signal_type=4,            # 待恢复波形:非对称冲激
    signal_amplitude=0.02,    # 峰值磁通(Φ₀)
    signal_center=100,        # 冲激中心(ns)
    reconstruction="wiener",  # 先用线性 Wiener 反卷积
    lambda_reg=5.0,           # Wiener 正则化强度
)

# ── 2. 跑完整管道:测量 + 重建 ────────────────────────────────────────
result = wf.run()                         # WorkflowResult
rec = result.reconstructed_signal         # FluxSignal:恢复出的 Φ(t)
print("恢复波形点数:", len(rec.signal))

# ── 3. 同一份测量上比较三种重建算法 ───────────────────────────────────
cmp = wf.compare(methods=["wiener", "hammerstein", "lm"])
print("最优方法:", cmp.best)             # 按 RMSE 自动选优
for m in cmp.methods:
    print(f"  {m}: rmse={cmp.metrics[m]['rmse']:.3e}  snr={cmp.metrics[m]['snr']:.2f}")

# ── 4. 自动出图:上图 Δp,下图恢复波形 vs 真值 ───────────────────────
wf.plot()
```

```{note}
`compare()` 复用 `run()` 缓存的那份测量数据,只换重建算法,不重跑
`mesolve`——三种方法比的是**同一批** $\Delta p$。`best` 按 RMSE 自动选优;因为这是
仿真,真值波形 `flux_samples` 已知,RMSE 才有意义。
```

## 结果解读

- `result.measurement` 里有四样东西:`p_e`(激发态布居)、`delta_p`(扣本底后的
  信号)、`kernel`(控制核 $k$)、`flux_samples`(仿真注入的真值波形,仅仿真可得)。
  轴信息在 `result.measurement.axes`:`scan`(滑动延迟)、`t_flux`(波形时间轴)。
- `result.reconstructed_signal` 是恢复出的 {py:class}`~sqc.control.FluxSignal`——
  即从纯布居数据倒推出的 $\Phi(t)$。`wf.plot()` 会把它叠在真值上,直观看重建质量。
- `compare()` 的 `metrics` 给每种方法一份 `rmse` / `snr` / `peak`。典型结论:
  Wiener 最快、对小信号够用;Hammerstein 补上 $\omega(\Phi)$ 的弱非线性,大幅值时
  RMSE 更低;LM 最准但最慢(全密度矩阵正向仿真 + 迭代)。`best` 字段替你挑出
  RMSE 最低者。
- 想看灵敏度极限,把这条管道接 `wf.sweep("signal.amplitude", [...])` 扫幅值,
  即得 SNR-幅值曲线——完整 API 见 {doc}`../building_blocks/workflows`,各重建器的
  公式与适用范围见 {doc}`../building_blocks/reconstruction`。
