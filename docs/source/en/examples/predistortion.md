# Predistortion

## Overview

Predistortion is the platform's third product mainline: the voltage waveform
emitted by an arbitrary waveform generator (AWG) is **distorted** as it
propagates through the control line (cables, filters, bias-tee, etc.) to the
superconducting chip—an ideal step develops an exponential tail. The core idea
of predistortion is to **pre-apply the inverse of the control-line transfer
function** at the AWG end, so that after the forward distortion of the control
line the waveform lands exactly as the target.

Physically, the control line is modelled as a linear time-invariant (LTI)
system, fully characterised by a transfer function $H(\omega)$ or equivalently
a step response $s(t)$. The actual on-chip flux waveform is the convolution of
the AWG output with the control-line impulse response:

$$\Phi_\text{chip}(t) = (h * V_\text{AWG})(t)
\quad\Longleftrightarrow\quad
\Phi_\text{chip}(\omega) = H(\omega) \, V_\text{AWG}(\omega).$$

Predistortion chains an **inverse filter** $H^{-1}(\omega)$ before the AWG
output, so the total transfer function reduces to the identity:

$$V_\text{AWG} = H^{-1} * \Phi_\text{target}
\;\Longrightarrow\;
\Phi_\text{chip} = H * H^{-1} * \Phi_\text{target} = \Phi_\text{target}.$$

## Pipeline Architecture

The predistortion pipeline executes in the following order: define target
waveform → model control-line distortion → measure transfer function → design
inverse filter → verify compensation. The sections below follow this logical
sequence, introducing each layer as it appears in the pipeline.

### Step 1: Define the Target Waveform — Control Layer (`Waveform`)

The starting point of any predistortion pipeline is the **target waveform**—the
flux signal you want to appear on the chip.
{py:class}`~sqc.control.waveform.Waveform` is a generic time-domain signal data
structure holding a time axis `t_list` and sample values `samples`. It provides
semantically neutral operations: `value_at(t)` (sample-and-hold),
`truncate(t_start, t_end)` (time-window truncation), `samples_on(t_global)`
(projection onto a global axis), and `copy()` (deep copy).

```{note}
In the predistortion pipeline both the target waveform and the AWG output are
passed as plain `Waveform` objects—at this stage the signal has not yet reached
the qubit and carries no flux semantics. The flux-specific subclass
{py:class}`~sqc.control.flux_signal.FluxSignal` only appears in the
quantum-simulation measurement path.
```

### Step 2: Model the Control-Line Distortion — Hardware Layer (`DistortionModel` + `ControlLine`)

After the target waveform leaves the AWG it must traverse a physical control
line to reach the chip. Two layers of abstraction model this transmission:
distortion models define the mathematical form; the control line wraps them in
physical context.

#### Distortion Models

{py:class}`~sqc.hardware.distortion.DistortionModel` is the abstract base class
for all distortion models, defining a uniform LTI interface:

- `apply(waveform, dt)`: apply distortion to uniformly-sampled data.
- `step_response(t)` / `impulse_response(t)` / `frequency_response(omega)`:
  the step response $s(t)$, impulse response $h(t)$, and complex frequency
  response $H(\omega)$.
- `apply_to_waveform(wf)`: convenience wrapper for `Waveform` objects.

Six concrete models span the range from simple to complex:

| Class | Mathematical Form | Typical Use |
|---|---|---|
| `SingleExponentialDistortion` | $s(t) = 1 - A e^{-t/\tau}$ | Most common flux-line tail |
| `MultiExponentialDistortion` | $s(t) = 1 - \sum_k A_k e^{-t/\tau_k}$ | Multi-stage filtering, reflections |
| `IIRDistortion` | $y[n] = \sum b_k x[n-k] - \sum_{k>0} a_k y[n-k]$ | General IIR filter, inverse filters |
| `FIRDistortion` | $y[n] = \sum b_k x[n-k]$ | General FIR filter, residual correction |
| `CustomTransferDistortion` | User-supplied $H(\omega)$ frequency grid | Arbitrary frequency-domain LTI system |
| `CascadeDistortion` | $H(\omega) = \prod_k H_k(\omega)$ | Multi-stage series compensation (IIR + FIR) |

All models support a `smooth` flag: `smooth=False` (default) preserves the
direct (delta-function) path, so a step input has a jump at $t=0$;
`smooth=True` gives pure-lowpass behaviour with no direct path.

#### Control Line

{py:class}`~sqc.hardware.control_line.ControlLine` places the distortion model
in a physical context, modelling a specific control line (xy/z/readout). It
holds a `transfer_function` (a `DistortionModel` instance) plus auxiliary
attributes: name, source/target, impedance, attenuation. Two key methods
correspond to the two directions of the pipeline:

- `apply(awg_waveform)` — **forward**: AWG waveform $\to$ control line $\to$
  on-chip waveform.
- `predistort(target, designer)` — **inverse**: target on-chip waveform $\to$
  inverse filter $\to$ required AWG waveform.

### Step 3: Measure the Transfer Function — Calibration Layer (`WaveformCalibration`)

With the distortion model and control line in place, the next step is to
**measure** the actual transfer function—inject a known signal, observe the
output, and fit the model parameters. This is the role of
`WaveformCalibration` in the calibration layer.

{py:class}`~sqc.calibration.waveform.WaveformCalibration` provides a unified
measurement and fitting interface with two paths:

- **Analytical path** (`distortion` + `fit_type`): directly fits the analytical
  step-response function of a known `DistortionModel`. `fit_type` may be
  `"single_exp"`, `"multi_exp"`, `"fir"`, or `"iir"`.
- **Quantum-simulation path** (`measurement_protocol` + `qubit` +
  `control_line`): measures the step response via a qubit-based sensing
  protocol (see below).

After fitting, `to_distortion_model()` exports the fitted parameters as a
`DistortionModel` subclass instance—this is the **forward** distortion model,
representing how the control line distorts signals.

### Step 4: Design the Inverse Filter — Calibration Layer (`PredistortionDesigner`)

With the forward transfer function measured, the next step is to invert it.
{py:class}`~sqc.calibration.waveform.PredistortionDesigner` is the standalone
inverse-filter designer:

| `method` | Algorithm | Applicability |
|---|---|---|
| `"auto"` | Auto-select based on model type | Default, recommended |
| `"iir_inverse"` | Analytical IIR inverse (bilinear transform or z-domain pole-zero) | Exponential distortion, closed-form |
| `"fir_inverse"` | Frequency-domain Wiener inverse + FIR truncation | General LTI distortion |
| `"frequency_inverse"` | $H^{-1} = \bar{H}/(|H|^2 + \varepsilon^2)$ | Arbitrary frequency-domain transfer function |

For the most common single-exponential distortion $s(t) = 1 - A e^{-t/\tau}$,
`"iir_inverse"` yields the IIR inverse coefficients analytically via bilinear
transform; `SingleExponentialDistortion` bundles this in its `design_inverse()`
method. For multi-exponential distortion, `design_inverse()` sorts components
by descending $\tau$, allocates an IIR inverse stage to each slow component,
and appends an FIR stage to correct the fast-component residual
(Rol 2020, §IV).

The `regularization` parameter controls how aggressive the inversion is: larger
values yield more conservative (stable) correction; smaller values risk noise
amplification where $|H|$ is near zero.

```{note}
`to_distortion_model()` returns the fitted **forward** distortion model; the
inverse filter must be obtained separately via
`PredistortionDesigner.design()`. Do not confuse the two.
```

### Step 5: End-to-End Orchestration — Workflow Layer (`PredistortionValidationWorkflow`)

{py:class}`~sqc.workflows.PredistortionValidationWorkflow` chains the four
steps above into a single validation pipeline. It introduces no new physics,
only orchestrates the existing layers:

1. Construct a `ControlLine` with the ground-truth distortion.
2. **Forward, uncorrected**: target $\to$ control line $\to$ on-chip waveform
   (with tail).
3. Measure transfer function and fit forward distortion model via
   `WaveformCalibration`.
4. Design inverse filter via `PredistortionDesigner`.
5. **Apply predistortion**: target $\to$ inverse filter $\to$ AWG waveform.
6. **Forward, corrected**: predistorted AWG waveform $\to$ control line $\to$
   on-chip waveform; compare with target; compute RMSE, improvement factor,
   and settling time.

Passing only `target_waveform` and `true_distortion` and calling `run()`
yields the full before/after comparison. The `measurement_protocol` field
switches between the analytical and quantum-simulation calibration paths.

### Quantum-Simulation Measurement Path

The pipeline above defaults to the analytical path—fitting the known
distortion's analytical step-response formula directly. To simulate the
**laboratory process of measuring the control-line transfer function via a
qubit**, pass `measurement_protocol` and `qubit` to
`PredistortionValidationWorkflow` or `WaveformCalibration`. The system then
invokes three additional layers:

**Device layer** — {py:class}`~sqc.devices.transmon.TransmonQubit` provides
the qubit physical model. Different protocols impose different flux-bias
requirements: Cryoscope works best at the sweet spot ($\Phi=0$); delay Ramsey
and transient protocols achieve highest sensitivity where
$\kappa = d\omega/d\Phi$ is maximal.

**Experiment layer** — the corresponding experiment class runs a quantum
simulation: injects a test signal (default: step) through the control line and
obtains the qubit response via `mesolve()`. **Reconstruction layer** — the
corresponding reconstruction class inverts the qubit response back to an
on-chip flux waveform; normalising by the step amplitude yields the step
response.

The core of the protocol-driven path is the internal
`_ProtocolDrivenMeasurement`: construct experiment $\to$ run `mesolve` $\to$
reconstruct on-chip waveform $\to$ normalise. Four protocols are supported:

| `measurement_protocol` | Experiment Class | Reconstruction Class | Principle |
|---|---|---|---|
| `"cryoscope"` | `CryoscopeExperiment` | `CryoscopeReconstruction` | Square flux + Ramsey phase demodulation |
| `"delay_ramsey"` | `DelayRamseyExperiment` | `DelayRamseyReconstruction` | Delayed Ramsey tail measurement |
| `"transient"` | `TransientSensingExperiment` | `TransientReconstruction` | Sliding Ramsey + Wiener deconvolution |
| `"pi_pulse"` | `PiPulseCompensationExperiment` | `PiPulseCompReconstruction` | $\pi$-pulse compensation tail recovery |

## Usage

### End-to-End Validation (Recommended)

```python
import numpy as np
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.workflows import PredistortionValidationWorkflow

# Target on-chip waveform: a 0.05 Φ₀ plateau
dt = 1.0  # ns
t = np.arange(0, 500, dt)
target = Waveform(t_list=t, samples=np.ones_like(t) * 0.05)

# Ground-truth distortion: single-exponential tail, A=0.04, τ=200 ns
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)

# One-call validation
wf = PredistortionValidationWorkflow(
    target_waveform=target,
    true_distortion=distortion,
)
result = wf.run()

m = result["metrics"]
print(f"uncorrected RMSE = {m['rmse_uncorrected']:.3e}")
print(f"corrected RMSE   = {m['rmse_corrected']:.3e}")
print(f"improvement      = {m['improvement_factor']:.1f}×")
print(f"uncorrected settling time = {m['settling_uncorrected_ns']:.2f} ns")
print(f"corrected settling time   = {m['settling_corrected_ns']:.2f} ns")
```

### Step-by-Step Pipeline

The following example follows the logical pipeline order: measure transfer
function → design inverse → apply predistortion.

```python
from sqc.calibration.waveform import WaveformCalibration, PredistortionDesigner
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion
import numpy as np

# Define target
dt = 1.0
t = np.arange(0, 500, dt)
target = Waveform(t_list=t, samples=np.ones_like(t) * 0.05)

# Known distortion
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)

# 1. Measure transfer function and fit distortion model
wf_cal = WaveformCalibration(
    distortion=distortion, fit_type="single_exp",
)
fwd_model = wf_cal.to_distortion_model()  # forward distortion model

# 2. Design inverse filter
designer = PredistortionDesigner(method="auto", regularization=1e-4)
inverse = designer.design(fwd_model, dt=dt)

# 3. Apply predistortion: target waveform → AWG waveform
awg_waveform = inverse.apply_to_waveform(target)

print(f"target points: {target.n_points}")
print(f"AWG points:    {awg_waveform.n_points}")
print(f"inverse type:  {type(inverse).__name__}")
```

### Quantum-Simulation Measurement Path

```python
from sqc.devices import TransmonQubit
from sqc.hardware.control_line import ControlLine
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.calibration.waveform import WaveformCalibration

# Qubit and distortion
qubit = TransmonQubit(
    EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
    T1=10_000, T2=8_000, flux=0.0,  # sweet spot for Cryoscope
)
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)

# Control line
line = ControlLine(
    name="Z0", kind="z", source="AWG0", target="Q0",
    transfer_function=distortion,
)

# Measure step response via Cryoscope quantum simulation
cal = WaveformCalibration(
    qubit=qubit,
    control_line=line,
    measurement_protocol="cryoscope",
    method="transfer_function",
    fit_type="single_exp",
)
fwd_model = cal.to_distortion_model()  # forward model from quantum simulation
```

## Reading the Results

- `result["target"]`: the target on-chip waveform (`Waveform`).
- `result["on_chip_uncorrected"]`: the on-chip waveform after forward
  distortion without predistortion (shows the tail).
- `result["on_chip_corrected"]`: the on-chip waveform after predistortion
  compensation.
- `result["awg_predistorted"]`: the waveform to be sent to the AWG after
  applying the inverse filter.
- `result["inverse_model"]`: the designed inverse filter object.
- `result["measured_model"]`: the forward distortion model obtained from
  calibration measurement and fitting.

`metrics` provides four key quantities:

- `rmse_uncorrected` / `rmse_corrected`: root-mean-square error between the
  on-chip waveform and the target, before and after correction.
- `improvement_factor`: $=$ `rmse_uncorrected / rmse_corrected`.
- `settling_uncorrected_ns` / `settling_corrected_ns`: time until the waveform
  settles within tolerance (default $10^{-3}$ of the target), before and after
  correction.

```{important}
The improvement factor on the analytical path can be extremely large. This is
because `WaveformCalibration` fits the analytical step-response function of the
known distortion directly—equivalent to knowing the exact form of $H$ before
inverting it—so the corrected RMSE can drop to machine precision
($\sim 10^{-15}$) and the improvement factor becomes commensurately large
(machine-precision reciprocal). This is the analytical ceiling and **does not
represent real-system performance**.

To obtain physically meaningful improvement factors, use the quantum-simulation
measurement path (pass `measurement_protocol="cryoscope"` and `qubit`), which
measures the step response via `mesolve` simulation. The fit is then limited by
measurement resolution and noise, and the improvement factor lands in a finite
range (typically $10$–$10^2\times$).
```

- `SingleExponentialDistortion` has a built-in `design_inverse()` method that
  analytically computes the IIR inverse filter via bilinear transform. This is
  the most commonly used analytical path.
- For more complex distortions (multi-exponential sums, arbitrary
  frequency-domain transfer functions), use `MultiExponentialDistortion`'s
  cascade design or `PredistortionDesigner`'s frequency-domain methods.
- See {doc}`../building_blocks/hardware` for the full field definitions and
  physical interpretation of each distortion model.
- See {doc}`../building_blocks/calibration` for the complete API of calibration
  and inverse filter design.
