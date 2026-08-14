# Hardware

## What this layer provides

The `hardware` layer models the analog link between the AWG output and the
qubits on chip. After leaving the arbitrary waveform generator, a signal travels
through coax lines, attenuators, filters, and bias-tees before reaching the
chip, distorted along the way by the transfer function of the cabling and filter
network. This layer describes that distortion as a linear time-invariant (LTI)
system and provides operations in both directions: forward (AWG waveform to real
on-chip waveform) and inverse/predistortion (desired on-chip waveform to the AWG
waveform to send).

This layer only describes the transfer-function distortion and readout
transduction of the link itself. It does not cover qubit physics (see
{doc}`devices`), pulse design (see {doc}`control`), or time evolution (see
{doc}`simulation`).

The inputs and outputs of this layer are uniformly the `Waveform` (a
semantically neutral time-domain container) and `FluxSignal` (its subclass, with
samples in units of $\Phi_0$) types defined in the {doc}`control` layer. Both
are defined in `control` and referenced here by `import`. Signal types are not
redefined in this layer because of the code dependency direction: `hardware`
depends on `control` but not the reverse, so placing the type definitions in
`control` avoids a circular import.

One convention runs through this layer: distortion models are pure LTI
descriptions. Given a sample spacing `dt`, a single model consistently yields
four self-consistent views, the time-domain `apply()`, step response, impulse
response, and frequency response, with DC gain fixed at 1 (a static bias is left
unchanged by the distortion; only transients are distorted).

## Class overview

| Class | Role | Notes |
|---|---|---|
| `DistortionModel` | Abstract base | Common contract (LTI) for all control-line distortions; the layer's main extension point |
| `SingleExponentialDistortion` | Distortion | Single exponential tail, the most common Z-line distortion model |
| `MultiExponentialDistortion` | Distortion | Sum of K single-exponential tails |
| `FIRDistortion` | Distortion | Finite impulse response filter |
| `IIRDistortion` | Distortion | Infinite impulse response filter |
| `CustomTransferDistortion` | Distortion | User-supplied $H(\omega)$ on a frequency grid |
| `CascadeDistortion` | Distortion | Series cascade of distortion stages (Rol 2020 §IV) |
| `ControlLine` | Control line | One physical control line (xy/z/readout), optionally carrying a distortion model |
| `TransferMatrix` | Transfer matrix | Frequency-domain relation across multiple Z lines, $\Phi_j(\omega)=\sum_i H_{ji}(\omega)V_i(\omega)$ |
| `ReadoutModel` | Abstract base | Common contract for readout models; the readout extension point |
| `IdealProjectiveReadout` | Readout | Projective measurement onto $|1\rangle$, returns $p_e$ |
| `IQReadoutModel` | Readout | Two Ramsey sequences with a $\pi/2$ phase offset for I/Q demodulation |

## Physical model and units

**Unit conventions**: time in ns, sample spacing `dt` from
`sqc.config.CONFIG.awg.dt`; angular frequency $\omega$ in rad/ns; flux in units
of $\Phi_0$; distortion amplitude `amplitude` is a dimensionless fraction
(typically $0.001$–$0.1$); tail time constant `tau` in ns (typically
$10$–$1000$ ns).

Physical origin of control-line distortion (following Gao 2021 §III.D on
control-line transfer functions and §V.E on calibrating distortion tails with
Cryoscope): an ideal step emitted by the AWG does not settle instantly after the
cabling and low-pass network. It carries a slowly relaxing exponential tail. The
single-exponential model is the most common description:

$$s(t) = 1 - A\,e^{-t/\tau}, \qquad H(s) = (1-A) + \frac{A}{1 + s\tau}$$

where $A$=`amplitude` and $\tau$=`tau`. The DC gain is $H(0)=1$. Every subclass
provides four self-consistent views:

- `apply(waveform, dt)`: time domain; applies the distortion to a uniformly
  sampled waveform (discretized via the bilinear transform).
- `step_response(t)`: step response $s(t)$.
- `impulse_response(t)`: impulse response $h(t)$.
- `frequency_response(omega)`: complex frequency response $H(\omega)$.

**The `smooth` flag (common to all distortion subclasses)**:

- `smooth=False` (default): keeps the direct (delta-function) feedthrough path;
  input discontinuities are preserved and the spectrum is broad; a step input
  jumps at $t=0$. This is the backward-compatible behavior.
- `smooth=True`: pure low-pass, no feedthrough path; an ideal step produces a
  smooth rise (e.g. $1-e^{-t/\tau}$) with no jump. Here `amplitude` is ignored
  (the tail weight is fixed at 1 to keep the DC gain = 1).

## DistortionModel: the distortion abstract base

The common contract for all control-line distortions, and the layer's main
extension point. The low-level interface uses numpy arrays; convenience methods
work directly on `control`-layer objects. Subclasses must implement four
abstract methods (the four views above):

- `apply(waveform: np.ndarray, dt: float) -> np.ndarray`: time-domain distortion.
- `step_response(t) -> np.ndarray` / `impulse_response(t) -> np.ndarray`.
- `frequency_response(omega) -> np.ndarray`: complex frequency response.

The base class provides two convenience wrappers (subclasses need not override
them):

- `apply_to_waveform(wf) -> Waveform`: takes a `Waveform`, returns the distorted
  new `Waveform`.
- `apply_to_signal(signal) -> FluxSignal`: takes a `FluxSignal`, returns a new
  `type=8` (user-defined raw samples) `FluxSignal`. This is the canonical path
  for reconstructed/distorted signals across the stack.

To add a custom distortion type, subclass `DistortionModel` and implement the
four abstract methods; see {doc}`../extending`.

## Built-in distortion subclasses

### SingleExponentialDistortion: single exponential tail

The most common Z-line distortion model, i.e. $s(t)=1-A e^{-t/\tau}$ above.

**Construction**

`SingleExponentialDistortion(amplitude=0.01, tau=100.0, smooth=False)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `amplitude` | float | Tail amplitude (dimensionless fraction) | `0.01` |
| `tau` | float | Time constant (ns) | `100.0` |
| `smooth` | bool | Pure low-pass mode (no feedthrough path) | `False` |

**Methods** (beyond the four inherited LTI views)

- `design_inverse(dt, formula="bilinear") -> IIRDistortion`: designs an IIR
  inverse filter that, cascaded with the forward distortion, cancels the
  exponential tail. `formula` is `"bilinear"` (bilinear transform) or
  `"rol2020"` (direct z-domain pole-zero design, Rol 2020 Eq. S22). A warning is
  raised when $|A|>0.5$ (Rol 2020 is validated for $|A|\le 0.1$).

**Output**

`apply()` returns `np.ndarray` (distorted samples); `apply_to_waveform()` returns
`Waveform`; `apply_to_signal()` returns `FluxSignal`; `design_inverse()` returns
`IIRDistortion`.

### MultiExponentialDistortion: multiple exponential tails

A sum of K single-exponential tails, describing distortion with several
coexisting timescales.

**Construction**

`MultiExponentialDistortion(amplitudes, taus, smooth=False)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `amplitudes` | `np.ndarray` | Per-component tail amplitudes (length K) | — |
| `taus` | `np.ndarray` | Per-component time constants (ns, length K) | — |
| `smooth` | bool | Pure low-pass mode | `False` |

The two arrays must be equal length. Read-only property `n_components` returns K.

**Methods** (beyond the four inherited LTI views)

- `design_inverse(dt, formula="bilinear", n_iir_stages=3, fir_taps=72,
  fir_threshold_ns=30.0) -> CascadeDistortion`: designs a cascade inverse filter
  by the iterative decomposition of Rol 2020 §IV: each slow component with $\tau$
  greater than `fir_threshold_ns` gets one IIR stage (up to `n_iir_stages`), then
  a single FIR stage corrects the fast residual.

**Output**

`apply()` returns `np.ndarray`; `design_inverse()` returns `CascadeDistortion`.

### FIRDistortion: finite impulse response

$y[n]=\sum_k b_k\,x[n-k]$.

**Construction**

`FIRDistortion(taps, smooth=False, smooth_tau=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `taps` | `np.ndarray` | FIR coefficients $b_k$ | — |
| `smooth` | bool | Pure low-pass mode | `False` |
| `smooth_tau` | float | Smoothing low-pass time constant (used only when `smooth=True`) | `None` |

**Class method**

- `FIRDistortion.design_from_residual(s_residual, dt, n_taps=72,
  ridge=1e-6) -> FIRDistortion`: solves (least squares with ridge regularization)
  for a set of taps that correct a residual step response back to a unit step;
  used by the FIR residual stage of `MultiExponentialDistortion.design_inverse`.

**Output**

`apply()` returns `np.ndarray`; `design_from_residual()` returns a `FIRDistortion`
instance.

### IIRDistortion: infinite impulse response

$y[n]=\sum_k b_k x[n-k]-\sum_{k>0} a_k y[n-k]$, $a_0=1$.

**Construction**

`IIRDistortion(b_coeffs, a_coeffs, smooth=False, smooth_tau=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `b_coeffs` | `np.ndarray` | Feedforward coefficients | — |
| `a_coeffs` | `np.ndarray` | Feedback coefficients (auto-normalized if $a_0\neq 1$) | — |
| `smooth` | bool | Pure low-pass mode | `False` |
| `smooth_tau` | float | Smoothing low-pass time constant | `None` |

This is the return type of `SingleExponentialDistortion.design_inverse`.

**Output**

`apply()` returns `np.ndarray`.

### CustomTransferDistortion: custom frequency-domain transfer function

The user supplies a complex $H(\omega)$ on a frequency grid; internally it runs
FFT → multiply by $H(\omega)$ → IFFT.

**Construction**

`CustomTransferDistortion(omega_grid, H_grid, smooth=False, smooth_tau=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `omega_grid` | `np.ndarray` | Frequency grid (rad/ns) | — |
| `H_grid` | `np.ndarray` | Corresponding complex frequency response | — |
| `smooth` | bool | Pure low-pass mode | `False` |
| `smooth_tau` | float | Smoothing low-pass time constant | `None` |

Frequencies outside the grid are extrapolated with the endpoint values.

**Output**

`apply()` returns `np.ndarray`.

### CascadeDistortion: series cascade

Applies several distortion stages in order:
$\text{out}=\text{stage}_N(\cdots\text{stage}_1(\text{in}))$; the composite
frequency response is the product of the per-stage responses; an empty `stages`
is the identity system.

**Construction**

`CascadeDistortion(stages)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `stages` | `list[DistortionModel]` | Distortion models (applied in list order) | — |

The `design_inverse` methods above return an instance of this class.

**Output**

`apply()` returns `np.ndarray`; `frequency_response()` returns the per-stage
product.

## ControlLine: a single physical control line

Models one physical control line (xy / z / readout), optionally carrying a
distortion model. It wraps a distortion model into a "source → target" link
object and additionally holds the link's physical parameters.

**Construction**

`ControlLine(name, kind, source, target, transfer_function=None, impedance=50.0, attenuation_db=20.0, delay=0.0, ...)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `name` | str | Control line identifier (e.g. `"Z0"`) | — |
| `kind` | str | Line type | `"xy"` / `"z"` / `"readout"` |
| `source` | str | Signal source identifier (e.g. `"AWG0:CH1"`) | — |
| `target` | str | Target device (e.g. `"Q0"`) | — |
| `transfer_function` | `DistortionModel` or None | Attached distortion model | `None` |
| `impedance` | float | Characteristic impedance (Ω) | `50.0` |
| `attenuation_db` | float | Total room-temperature-to-chip attenuation (dB) | `20.0` |
| `delay` | float | Propagation delay (ns) | `0.0` |
| `metadata` | dict | Additional metadata | `{}` |

**Methods**

- `apply(awg_waveform: Waveform) -> Waveform`: **forward**: AWG waveform →
  on-chip waveform. Returns the input unchanged when there is no
  `transfer_function` (with a time shift if `delay>0`).
- `predistort(target_waveform, designer) -> Waveform`: **inverse**: desired
  on-chip waveform → the AWG waveform to send, delegated to the `calibration`
  layer's `PredistortionDesigner` (see {doc}`calibration`).

**Output**

`apply()` returns `Waveform` (on-chip waveform); `predistort()` returns
`Waveform` (predistorted AWG waveform).

## TransferMatrix: multi-line transfer matrix

Describes, in the frequency domain, the relation across multiple Z lines. The
flux at a target qubit $j$ is determined jointly by all source-line voltages:

$$\Phi_j(\omega) = \sum_i H_{ji}(\omega)\,V_i(\omega)$$

The diagonal elements $H_{ii}$ are each line's response on its own target; the
off-diagonal elements $H_{ji}\,(j\neq i)$ describe the part of one line's voltage
that crosstalks into another qubit's flux.

```{note}
This page presents `TransferMatrix` only as a **forward** device (mapping AWG
voltages to on-chip flux). Reconstructing the crosstalk matrix from
measurements is not part of the public v1 scope.
```

**Construction**

`TransferMatrix(elements, frequency_axis, time_axis=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `elements` | dict | `(target_name, source_name) -> H_ji(ω)` complex-array dict | — |
| `frequency_axis` | `np.ndarray` | Frequency grid (rad/ns) | — |
| `time_axis` | `np.ndarray` | Time-domain axis (optional) | `None` |

**Read-only properties**: `sources` / `targets` (sorted source/target name lists).

**Methods**

- `H_ji(target, source) -> np.ndarray`: access a single transfer element.
- `diagonal() -> dict`: each target's self-response $H_{ii}$.
- `off_diagonal() -> dict`: the crosstalk elements $H_{ji}\,(j\neq i)$.
- `apply(source_voltages: dict[str, Waveform]) -> dict[str, FluxSignal]` —
  **forward**: for each target $j$, computes
  $\Phi_j=\text{IFFT}(\sum_i H_{ji}\cdot\text{FFT}(V_i))$; waveforms of unequal
  length are zero-padded to the longest.

**Class method**

- `TransferMatrix.from_dc_matrix(dc_matrix, source_names, target_names,
  n_freq=256) -> TransferMatrix`: builds a **frequency-flat** transfer matrix
  from a DC crosstalk matrix ($H_{ji}(\omega)$ is the constant `dc_matrix[j, i]`
  at all $\omega$).

**Output**

`apply()` returns `dict[str, FluxSignal]` (target name → on-chip flux);
`H_ji()` returns `np.ndarray`.

## Readout models

The readout extension point and two built-in implementations. Measurement
results are returned uniformly as `dict[str, float]`.

### ReadoutModel: the readout abstract base

The common contract for readout models and the **readout extension point**.
Subclasses implement the single abstract method
`measure(*args, **kwargs) -> dict[str, float]`.

### IdealProjectiveReadout: ideal projective measurement

Projective measurement onto $|1\rangle$.

**Construction**

`IdealProjectiveReadout()`

**Methods**

- `measure(state, qubit=None) -> dict[str, float]`: returns `{"p_e": float}`,
  the excited-state population $p_e=\langle 1|\rho|1\rangle$.

**Output**

`dict` with single key `"p_e"`.

### IQReadoutModel: I/Q readout

Runs two Ramsey sequences with a $\pi/2$ phase offset for I/Q demodulation
(replaces the legacy `src/protocal.py:IQ_readout`).

**Construction**

`IQReadoutModel(tau=20.0, t_rabi=None, omega_d=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `tau` | float | Free precession time (ns) | `20.0` |
| `t_rabi` | `np.ndarray` | Rabi time axis (ns) | `CONFIG.pulse.t_rabi` |
| `omega_d` | float | Drive frequency (rad·GHz) | `None` (uses qubit frequency) |

**Methods**

- `measure(qubit, **extra) -> dict[str, float]`: returns
  `{"p_e_I": float, "p_e_Q": float}`; the `qubit` must have called
  `qubit_in_mag(...)` first to populate `H_list`.

**Output**

`dict` with keys `"p_e_I"` and `"p_e_Q"`.

## Minimal example

```python
import numpy as np
from sqc.config import CONFIG
from sqc.control.waveform import Waveform
from sqc.hardware import SingleExponentialDistortion, ControlLine

# Sample spacing derives from CONFIG; build the time axis with arange (not linspace)
dt = CONFIG.awg.dt
t = np.arange(0.0, 200.0, dt)

# An ideal step waveform (what the AWG intends to send)
awg = Waveform(t_list=t, samples=np.where(t >= 20.0, 1.0, 0.0))

# Single-exponential tail distortion: A=2%, tau=50 ns
dist = SingleExponentialDistortion(amplitude=0.02, tau=50.0)

# Attach to a Z line; forward gives the real on-chip waveform
line = ControlLine(name="Z0", kind="z", source="AWG0:CH1",
                   target="Q0", transfer_function=dist)
on_chip = line.apply(awg)          # the step edge is broadened by the exponential tail

# Design an inverse filter; cascaded, it cancels the tail (predistortion)
inv = dist.design_inverse(dt, formula="rol2020")
corrected = inv.apply_to_waveform(on_chip)   # tail compensated, restored to an ideal step
```

## Physical role / extension

- The `DistortionModel` family corresponds to the real control-line transfer
  function: the LTI distortion caused jointly by cabling, attenuators, low-pass
  filters, and bias-tees. `SingleExponentialDistortion` corresponds to the
  slow-relaxation tail most commonly seen in Cryoscope calibration
  (Gao 2021 §V.E).
- `ControlLine` corresponds to one physical wire on chip (xy drive line / z flux
  line / readout line); `TransferMatrix` corresponds to the frequency-domain
  response network when multiple Z lines coexist.
- The `ReadoutModel` family corresponds to the readout chain:
  `IdealProjectiveReadout` is the noiseless upper bound of projective
  measurement, `IQReadoutModel` corresponds to real I/Q demodulated readout.
- To add a custom distortion type (e.g. a transfer function with ripple, a
  non-minimum-phase filter), subclass `DistortionModel` and implement `apply` /
  `step_response` / `impulse_response` / `frequency_response`. To add a custom
  readout, subclass `ReadoutModel` and implement `measure`. See the full
  extension guide in {doc}`../extending`.
