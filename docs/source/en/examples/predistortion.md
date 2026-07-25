# Predistortion

## Goal

The voltage waveform an AWG emits is **distorted** by the time it reaches the
chip through the control line (cables, filters, bias-tee) — an ideal square pulse
develops an exponential tail. Predistortion **pre-applies the inverse transform**
to the AWG waveform so that, after passing through the control line, it lands as
exactly the target waveform. This is the platform's third product mainline:
measure the control-line transfer function → design the inverse filter → verify
the correction. This example runs the whole validation pipeline in one call with
{py:class}`~sqc.workflows.PredistortionValidationWorkflow`.

## Physics

The control line is a **linear time-invariant (LTI) system**, described by a
transfer function $H(\omega)$ or the equivalent step response $s(t)$. The on-chip
waveform is the convolution of the AWG waveform with the control-line impulse
response:

$$\Phi_\mathrm{chip}(t) = (h * V_\mathrm{AWG})(t) \quad\Longleftrightarrow\quad
\Phi_\mathrm{chip}(\omega) = H(\omega)\, V_\mathrm{AWG}(\omega).$$

The most common distortion is a **single-exponential tail**: an ideal step
through the control line becomes $s(t) = 1 - A\,e^{-t/\tau}$ — a tail of amplitude
$A$ and time constant $\tau$ before it reaches its final value (see
`DistortionModel` in {doc}`../building_blocks/hardware`).

Predistortion chains an **inverse filter** $H^{-1}(\omega)$ onto the AWG waveform
so the total transfer function becomes $H(\omega)\,H^{-1}(\omega) = 1$:

$$V_\mathrm{AWG} = H^{-1} * \Phi_\mathrm{target}
\;\Longrightarrow\;
\Phi_\mathrm{chip} = H * H^{-1} * \Phi_\mathrm{target} = \Phi_\mathrm{target}.$$

For exponential distortions the inverse filter has an analytical IIR form; general
distortions use frequency-domain inversion
$H^{-1}=\bar H/(|H|^2+\varepsilon^2)$ ($\varepsilon$ is a regularisation that
suppresses noise amplification where $|H|$ is near zero). Both routes are provided
uniformly by {py:class}`~sqc.calibration.PredistortionDesigner`.

## End-to-end code

Run the whole validation pipeline in one call (inject distortion → measure
transfer function → design inverse filter → apply predistortion → re-measure →
compute improvement):

```python
import numpy as np
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.workflows import PredistortionValidationWorkflow

# ── target on-chip waveform: a 0.05 Φ₀ plateau ───────────────────────
dt = 1.0                                  # ns
t = np.arange(0, 500, dt)
target = Waveform(t_list=t, samples=np.ones_like(t) * 0.05)

# ── ground-truth distortion: single-exp tail, amplitude 0.04, tau 200 ns ──
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)

# ── one-call validation ──────────────────────────────────────────────
wf = PredistortionValidationWorkflow(
    target_waveform=target,
    true_distortion=distortion,
)
result = wf.run()

m = result["metrics"]
print(f"uncorrected RMSE = {m['rmse_uncorrected']:.3e}")
print(f"corrected RMSE   = {m['rmse_corrected']:.3e}")
print(f"improvement      = {m['improvement_factor']:.1f}×")
```

You can also walk the low-level three steps manually, to see what "measure
transfer function → design inverse → apply" each do:

```python
from sqc.calibration import WaveformCalibration, PredistortionDesigner

# 1) measure transfer function + fit a distortion model
wf_cal = WaveformCalibration(distortion=distortion, fit_type="single_exp")
fwd_model = wf_cal.to_distortion_model()      # the fitted **forward** distortion model

# 2) design the inverse filter
designer = PredistortionDesigner(method="auto")
inverse = designer.design(fwd_model, dt=dt)   # inverse model H⁻¹

# 3) apply predistortion to the target → the waveform to send to the AWG
awg_waveform = inverse.apply_to_waveform(target)
# or in one step: designer.predistort(target, transfer=distortion)
```

```{note}
`to_distortion_model()` returns the fitted **forward** distortion model; the
inverse filter comes from `PredistortionDesigner.design()` — don't confuse the
two. `method="auto"` uses an analytical IIR inverse for exponential models and
frequency-domain inversion otherwise.
```

## Reading the results

- `result` is a dict with six things: `target` (the goal),
  `on_chip_uncorrected` (the on-chip waveform without correction, with its tail),
  `on_chip_corrected` (the corrected on-chip waveform), `awg_predistorted` (the
  predistorted waveform sent to the AWG), `inverse_model` (the inverse filter),
  and `measured_model` (the measured transfer-function model).
- `metrics` reports `rmse_uncorrected` / `rmse_corrected` /
  `improvement_factor`, plus the pre- and post-correction settling times
  `settling_uncorrected_ns` / `settling_corrected_ns`.
- **The improvement factor will be enormous**: this example takes the
  **analytical path** — it fits the ground-truth distortion's step response and
  inverts it analytically, effectively knowing the exact form of $H$, so the
  corrected RMSE can drop to machine precision ($\sim10^{-15}$) and the
  improvement factor becomes astronomical. This is the analytical ceiling and
  **does not represent a real system**.
- A real calibration passes `qubit` + `measurement_protocol` (e.g.
  `"cryoscope"`) to `PredistortionValidationWorkflow`, measuring the step response
  via quantum simulation — the fit has finite precision and the inverse filter is
  imperfect, so the improvement factor lands at a finite scale (typically
  10–100×) that better reflects the real gain from predistortion.
- For a more realistic distortion, use `MultiExponentialDistortion` (multiple
  tails) or `IIRDistortion`; a larger designer `regularization` is more stable but
  more conservative in its correction. Fields for each distortion / inverse-filter
  class are in {doc}`../building_blocks/hardware` and
  {doc}`../building_blocks/calibration`.
