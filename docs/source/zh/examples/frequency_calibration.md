# 频率标定

## 目标

测出 qubit 频率随磁通的响应曲线 $f_{01}(\Phi)$,并据此回答两个实用问题:**在某个
磁通偏置下 qubit 频率是多少**,以及**要把频率调到某个目标值该施加多少磁通偏置**。
这是平台的第二条产品主线:不重建外界信号,而是**表征器件本身**——投片后确定工作点、
建查表、必要时闭环整定到目标频率。本例直接用 {doc}`../building_blocks/calibration`
层的标定类。

## 物理原理

Transmon 频率由外部磁通经 SQUID 环调制约瑟夫森能而定(见 {doc}`../theory`):

$$f_{01}(\Phi) \approx \frac{1}{2\pi}\left(\sqrt{8 E_J(\Phi)\, E_C} - E_C\right),
\qquad E_J(\Phi) = E_{J0}\,\lvert\cos(\pi\Phi/\Phi_0)\rvert.$$

因为 $E_J\propto\lvert\cos(\pi\Phi/\Phi_0)\rvert$,$f_{01}(\Phi)$ 是关于 $\Phi=0$ 的
**偶函数**,在整数磁通量子处取极大——即**甜点(sweet spot)**,那里
$\mathrm{d}f/\mathrm{d}\Phi=0$,对磁通噪声一阶不敏感。

单点测频用 Ramsey 序列:施加已知的**人工失谐** $f_a$,让 qubit 在自由进动期
以 $\lvert f_a - \Delta\rvert$ 的频率振荡($\Delta$ 是真实失谐),对 $p_e(\tau)$ 做 FFT
取峰即得 $f_{01}$。逐点扫磁通、每点测一次频率,就画出 $f_{01}(\Phi)$ 曲线。

有了曲线,**正向**用插值($\Phi\to f$),**反向**用反插值($f\to\Phi$)。注意曲线是
偶函数、整体不单调,反插值只在**单调的一支**上有定义——所以查目标频率对应的磁通时,
要落在甜点一侧。

## 端到端代码

```python
import numpy as np
from sqc.devices.transmon import TransmonQubit
from sqc.calibration import FluxResponseCalibration, FrequencyMeasurement

# ── 器件:EC/EJ 用角频率(rad·GHz)传入 ────────────────────────────────
qubit = TransmonQubit(
    EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
    T1=10_000, T2=8_000, flux=0.0, n_levels=3,
)

# ── 1. 单点测频:甜点处的 f01 ─────────────────────────────────────────
fm = FrequencyMeasurement(qubit=qubit, flux=0.0, method="ramsey")
f01 = fm.measure()                       # 有符号角频率(rad·GHz)
print(f"甜点 f01 = {f01 / (2*np.pi):.4f} GHz")

# ── 2. 扫磁通建 f(Φ) 查表 ────────────────────────────────────────────
#     h_list 越密越准,但每点一次 Ramsey 扫描,耗时线性增长(见下方说明)。
cal = FluxResponseCalibration(
    qubit=qubit,
    method="ramsey",
    h_list=np.linspace(-0.03, 0.03, 5),  # 演示用粗网格
)
table = cal.calibrate()                  # CalibrationTable, kind="f_phi"

# ── 3. 正向 / 反向查表 ───────────────────────────────────────────────
f_at_bias = table.evaluate(np.array([0.015]))          # Φ → f
f_target = table.outputs.max() * 0.999                 # 略低于甜点(留在单调支)
phi_needed = table.inverse(np.array([f_target]))       # f → Φ
print(f"目标 f={f_target/(2*np.pi):.4f} GHz 需偏置 Φ={phi_needed[0]:.4f}")
```

反插值给的是**开环估计**。要更准,可把 $f(\Phi)$ 给出的磁通范围当括号,用
`SinglePointFrequencyCalibration` **闭环整定**到目标频率:

```python
from sqc.calibration import SinglePointFrequencyCalibration

tuner = SinglePointFrequencyCalibration(
    qubit=qubit,
    f_target=f_target,
    V_a=0.0, V_b=0.03,          # 括号边界,取自 f(Φ) 的单调支
    step_method="secant",       # 割线法,典型 1–3 次迭代收敛
)
result = tuner.calibrate()      # CalibrationTable, kind="f01"
print("整定后偏置:", result.fit_params["V_opt"],
      "收敛:", result.fit_params["converged"])
```

```{note}
`FluxResponseCalibration` 每个磁通点跑一次完整 Ramsey τ 扫描(几十次 `mesolve`),
默认 `h_list` 有 51 个点——完整标定是**分钟级**任务。上面的 5 点粗网格只为演示流程;
真跑请按精度需求加密,或只在关心的频段附近细扫。
```

## 结果解读

- `fm.measure()` 返回单个工作点的有符号角频率($\mathrm{rad\cdot GHz}$);除以 $2\pi$
  得 GHz。默认单扫模式假设 $\lvert\Delta\rvert<0.1$ GHz;若可能远离甜点,置
  `f_artificial=None` 走双扫模式,对任意失谐稳健且给符号。
- `table` 是 `kind="f_phi"` 的 {py:class}`~sqc.calibration.CalibrationTable`:
  `inputs` 是磁通点、`outputs` 是对应角频率。`evaluate` 做磁通→频率的三次样条插值,
  `inverse` 做频率→磁通的反插值。因 $f(\Phi)$ 是偶函数,`inverse` 自动取单调支——
  查询目标频率务必落在甜点一侧,否则解不唯一。
- 典型曲线在 $\Phi=0$ 取极大(甜点),两侧对称下降;峰值即器件的最高工作频率。
  想要一阶抗磁通噪声就工作在甜点,想要传感灵敏度($\kappa=\mathrm{d}\omega/\mathrm{d}\Phi$
  有限)就偏置到旁边——这正是 {doc}`waveform_reconstruction` 里 `flux_bias` 的由来。
- 闭环整定的 `result.fit_params` 里有完整 `history`(每次迭代的 $V$、$f$、残差),
  可用来画收敛曲线;`converged` 标记是否在 `epsilon_f` 容差内收敛。
- 各标定类的字段、方法选项与调度机制见 {doc}`../building_blocks/calibration`;
  标定表如何服务重建的查表反演见 {doc}`../building_blocks/reconstruction`。
