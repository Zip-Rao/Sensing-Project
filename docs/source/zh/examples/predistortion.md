# 预畸变

## 目标

AWG 发出的电压波形,经控制线(线缆、滤波器、bias-tee)到达芯片时会被**畸变**——
理想方波会拖出一条指数尾巴。预畸变(predistortion)就是**预先给 AWG 波形做逆变换**,
让它经过控制线后恰好还原成目标波形。这是平台的第三条产品主线:测控制线传函 →
设计逆滤波器 → 验证补偿效果。本例用
{py:class}`~sqc.workflows.PredistortionValidationWorkflow` 一键跑完整条验证管道。

## 物理原理

控制线是一个**线性时不变(LTI)系统**,用传递函数 $H(\omega)$ 或等价的阶跃响应
$s(t)$ 描述。片上波形是 AWG 波形与控制线冲激响应的卷积:

$$\Phi_\mathrm{chip}(t) = (h * V_\mathrm{AWG})(t) \quad\Longleftrightarrow\quad
\Phi_\mathrm{chip}(\omega) = H(\omega)\, V_\mathrm{AWG}(\omega).$$

最常见的畸变是**单指数尾巴**:一个理想阶跃过控制线后变成
$s(t) = 1 - A\,e^{-t/\tau}$,即到达终值前有一段幅度 $A$、时间常数 $\tau$ 的拖尾
(见 {doc}`../building_blocks/hardware` 的 `DistortionModel`)。

预畸变的思路是给 AWG 波形串一个**逆滤波器** $H^{-1}(\omega)$,使总传函变成
$H(\omega)\,H^{-1}(\omega) = 1$:

$$V_\mathrm{AWG} = H^{-1} * \Phi_\mathrm{target}
\;\Longrightarrow\;
\Phi_\mathrm{chip} = H * H^{-1} * \Phi_\mathrm{target} = \Phi_\mathrm{target}.$$

对指数型畸变,逆滤波器有解析的 IIR 形式;一般畸变则走频域反演
$H^{-1}=\bar H/(|H|^2+\varepsilon^2)$($\varepsilon$ 是正则化,压制 $|H|$ 近零处的噪声放大)。
两条路都由 {py:class}`~sqc.calibration.PredistortionDesigner` 统一提供。

## 端到端代码

一键跑完整验证管道(注入畸变 → 测传函 → 设计逆滤波 → 施加预畸变 → 复测 → 算改善):

```python
import numpy as np
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.workflows import PredistortionValidationWorkflow

# ── 目标片上波形:一段 0.05 Φ₀ 的平台 ────────────────────────────────
dt = 1.0                                  # ns
t = np.arange(0, 500, dt)
target = Waveform(t_list=t, samples=np.ones_like(t) * 0.05)

# ── 真值畸变:幅度 0.04、时间常数 200 ns 的单指数尾巴 ────────────────
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)

# ── 一键验证 ─────────────────────────────────────────────────────────
wf = PredistortionValidationWorkflow(
    target_waveform=target,
    true_distortion=distortion,
)
result = wf.run()

m = result["metrics"]
print(f"补偿前 RMSE = {m['rmse_uncorrected']:.3e}")
print(f"补偿后 RMSE = {m['rmse_corrected']:.3e}")
print(f"改善倍数    = {m['improvement_factor']:.1f}×")
```

也可以手动走底层三步,看清「测传函 → 设计逆 → 施加」各自在做什么:

```python
from sqc.calibration import WaveformCalibration, PredistortionDesigner

# 1) 测传函 + 拟合成畸变模型
wf_cal = WaveformCalibration(distortion=distortion, fit_type="single_exp")
fwd_model = wf_cal.to_distortion_model()      # 拟合出的**正向**畸变模型

# 2) 设计逆滤波器
designer = PredistortionDesigner(method="auto")
inverse = designer.design(fwd_model, dt=dt)   # 逆模型 H⁻¹

# 3) 把预畸变施加到目标波形,得到该送进 AWG 的波形
awg_waveform = inverse.apply_to_waveform(target)
# 或一步到位:designer.predistort(target, transfer=distortion)
```

```{note}
`to_distortion_model()` 返回的是拟合出的**正向**畸变模型,逆滤波器要再经
`PredistortionDesigner.design()` 得到——两者别弄混。`method="auto"` 对指数型模型走
解析 IIR 逆,其余走频域反演。
```

## 结果解读

- `result` 字典有六样东西:`target`(目标)、`on_chip_uncorrected`(不补偿时的片上
  波形,带尾巴)、`on_chip_corrected`(补偿后的片上波形)、`awg_predistorted`(送进
  AWG 的预畸变波形)、`inverse_model`(逆滤波器)、`measured_model`(测出的传函模型)。
- `metrics` 给出 `rmse_uncorrected` / `rmse_corrected` / `improvement_factor`,外加
  补偿前后的整定时间 `settling_uncorrected_ns` / `settling_corrected_ns`。
- **改善倍数会非常大**:本例走**解析路径**——直接拿真值畸变的阶跃响应拟合,再解析求逆,
  相当于知道 $H$ 的精确形式,补偿后 RMSE 可低到机器精度($\sim10^{-15}$),改善倍数因而
  达到天文数字。这是解析上限,**不代表真实系统**。
- 真实标定要给 `PredistortionValidationWorkflow` 传 `qubit` + `measurement_protocol`
  (如 `"cryoscope"`),走量子仿真测阶跃响应——拟合有限精度、逆滤波非理想,改善倍数会
  落到有限的量级(典型 10–100×),更能反映实际预畸变的收益。
- 想换更贴近真实的畸变,可用 `MultiExponentialDistortion`(多条尾巴)或
  `IIRDistortion`;设计器的 `regularization` 越大越稳但补偿越保守。各畸变/逆滤波类的
  字段见 {doc}`../building_blocks/hardware` 与 {doc}`../building_blocks/calibration`。
