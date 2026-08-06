# Experiments

## What this layer provides

The `experiments` layer is the stack's protocol-orchestration layer. It
assembles the parts from the layers below (devices, flux signals, control
pulses, solvers, readout models) into a complete, runnable experiment following
some sensing protocol. Each experiment class answers a version of the same
question: given this pulse sequence and this parameter sweep, how does the
qubit's final-state population or phase vary with the parameter? It does no
numerical integration itself (that goes to `mesolve` in the {doc}`simulation`
layer) and no field reconstruction (that goes to the {doc}`reconstruction`
layer); it only orchestrates and produces raw measurement data.

One convention runs through this layer: every experiment is a `@dataclass`,
constructed with a device and optional parameters, and `run()` executes the
whole experiment in one call and returns the result. Almost all experiments
return the uniform {py:class}`~sqc.simulation.ExperimentResult` (the four dicts
`data`/`axes`/`metadata`/`config`), ready to feed the reconstruction layer. All
default parameters match the corresponding case of the legacy
`src/protocal.py`, keeping the physics regression-testable.

```{note}
The one exception: {py:class}`~sqc.experiments.RabiExperiment`'s `run()` returns
a raw QuTiP `Result` (not an `ExperimentResult`), to keep the return type
identical to the legacy case 0. The other 7 experiments all return an
`ExperimentResult`.
```

## Class overview

| Class | Role | Key outputs |
|---|---|---|
| `Experiment` | Abstract base (this layer's extension point) | — |
| `RabiExperiment` | Rabi oscillation: constant-amplitude drive, sweep duration | `p_e(t)` (raw QuTiP `Result`) |
| `RamseyExperiment` | Ramsey: $\pi/2-\tau-\pi/2$, measures free-evolution accumulated phase | `p_e(\tau)` |
| `DiffEchoExperiment` | Differential echo: $k$ repetitions accumulate phase, enhancing weak signals | `p_e(\tau)` |
| `TransientSensingExperiment` | Transient-field sensing: sliding measurement + control kernel | `delta_p`, `kernel` |
| `CryoscopeExperiment` | Cryoscope: scan truncation delay, IQ readout of $\varphi(t_d)$ | `varphi(t_d)` |
| `DelayRamseyExperiment` | Delay Ramsey: slide a short Ramsey across the falling edge to measure the tail | `varphi(t_d)` |
| `PiPulseCompensationExperiment` | Pi-pulse compensation: 2D scan of $(\tau, z)$ to measure the tail | `z_star(\tau)` |

These experiments map onto three product mainlines: waveform reconstruction
(Ramsey / DiffEcho / Transient / Cryoscope recovering $\Phi(t)$), frequency
calibration (Ramsey `f(\Phi)`), and predistortion (DelayRamsey /
PiPulseComp measuring the AWG→chip tail). Their end-to-end pipelines are in
{doc}`../examples/waveform_reconstruction`,
{doc}`../examples/frequency_calibration`, and
{doc}`../examples/predistortion`.

## Experiment: experiment abstract base class

The common contract for all sensing experiments, and this layer's extension
point. It mandates two abstract methods:

- `build_sequence()`: build the control pulse sequence for this experiment
  (some experiments build it per-parameter, in which case this returns `None`).
- `run(*args, **kwargs)`: execute the whole experiment and return the result.

It also provides one shared helper, `_route_flux(signal)`: pass a flux signal
through an optional {py:class}`~sqc.hardware.ControlLine` distortion (the
`control_line` field) before coupling it to the qubit; with `control_line=None`
(the default) it passes through. This is the predistortion mainline's uniform
entry point for injecting transfer-function distortion.

To add a custom sensing protocol, subclass `Experiment` and implement
`build_sequence()` / `run()`; `run()` should return an `ExperimentResult`. See
{doc}`../extending`.

## RabiExperiment: Rabi oscillation

Applies a constant-amplitude drive pulse and measures the excited-state
population versus time, used to calibrate $\pi$/$\pi/2$ pulse durations.

**Construction**

`RabiExperiment(qubit, t_rabi=None, omega_d=None, omega_rabi=1.0)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test | — |
| `t_rabi` | `np.ndarray` | Rabi time axis (ns) | `make_time(0, 40)` |
| `omega_d` | float | Drive frequency (rad·GHz) | Qubit sweet-spot frequency |
| `omega_rabi` | float | Constant Rabi frequency $\Omega$ (rad·GHz) | `1.0` |

Because the envelope is a fixed amplitude $\Omega$, the excited-state
population traces full Rabi oscillations $p_e = \sin^2(\Omega t / 2)$ over
`t_rabi`, roughly $\Omega\,t_\mathrm{rabi}[-1] / (2\pi)$ periods (about 6 full
swings between 0 and 1 with the defaults). Increase `omega_rabi` for a faster
drive with more periods; the $\pi$ pulse sits at the first $p_e=1$ peak,
$t_\pi = \pi/\Omega$.

**Methods**

- `build_sequence() -> Pulse`: returns a constant-amplitude
  {py:class}`~sqc.control.Pulse` (rotating frame, RWA) with
  `amplitude = omega_rabi`.
- `run() -> qutip.Result`: runs `mesolve` on the global time axis
  `CONFIG.pulse.t_global` and returns a raw QuTiP `Result`.

**Output**

Returns a QuTiP `Result` (⚠️ **not** an `ExperimentResult`; the one exception in
this layer). Take `result.expect[0]` for $p_e(t)$.

## RamseyExperiment: Ramsey interferometry

A $\pi/2-\tau-\pi/2$ sequence that reflects qubit-frequency shifts caused by an
external flux $\Phi(t)$ through the phase accumulated during free evolution; it
is the core protocol for both the waveform-reconstruction and
frequency-calibration mainlines.

**Construction**

`RamseyExperiment(qubit, flux_signal=None, omega_d=None, tau_list=None, t_global=None, phase1=π/2, phase2=0.0)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test | — |
| `flux_signal` | `FluxSignal` | Applied flux signal | Constant signal |
| `omega_d` | float | Drive frequency (rad·GHz) | Qubit sweet-spot frequency |
| `tau_list` | `np.ndarray` | Free-evolution time sweep (ns) | `CONFIG.pulse.tau_list` |
| `t_global` | `np.ndarray` | Global time axis (ns) | `CONFIG.pulse.t_global` |
| `phase1` | float | First $\pi/2$ pulse phase | `π/2` |
| `phase2` | float | Second $\pi/2$ pulse phase | `0.0` |

**Methods**

- `run() -> ExperimentResult`: first projects the flux signal onto the global
  time axis and couples it into the qubit via `qubit_in_mag` (rotating frame),
  then builds a Ramsey sequence per `tau`, runs the evolution, and collects the
  final $p_e$.

**Output**

Returns `ExperimentResult` with:

| Key | Meaning |
|---|---|
| `data["p_e"]` | $p_e(\tau)$ array |
| `data["flux_samples"]` | Flux sample array |
| `axes["tau"]` | Free-evolution time axis |
| `axes["t_flux"]` | Flux time axis |

## DiffEchoExperiment: differential echo

The differential echo protocol: sequence
$\pi/2-[\tau-\pi-\tau'-(\tau+t_\mathrm{int})-\pi-\tau'']^k-\pi/2$, which enhances
weak-signal sensitivity by accumulating phase over $k$ echo repetitions.

**Construction**

`DiffEchoExperiment(qubit, flux_signal=None, k=5, t_rabi=None, t_int=None, t_rep=None, tau_list=None, t_global=None, omega_d=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test | — |
| `flux_signal` | `FluxSignal` | Applied flux signal | Gaussian signal |
| `k` | int | Number of echo repetitions | `5` |
| `t_rabi` | `np.ndarray` | Rabi time axis (ns) | `CONFIG.pulse.t_rabi` |
| `t_int` | float | Interaction time (ns) | — |
| `t_rep` | float | Repetition period (ns) | — |
| `tau_list` | `np.ndarray` | Free-evolution sweep axis (ns) | `CONFIG.pulse.tau_list` |
| `t_global` | `np.ndarray` | Global time axis (ns) | `CONFIG.pulse.t_global` |
| `omega_d` | float | Drive frequency (rad·GHz) | Qubit sweet-spot frequency |

**Methods**

- `run() -> ExperimentResult`: replicates the base flux signal $2k$ times,
  concatenates and projects it onto the global axis and couples it into the
  qubit, then builds a differential echo pulse per `tau` and runs the evolution.

**Output**

Returns `ExperimentResult` with:

| Key | Meaning |
|---|---|
| `data["p_e"]` | $p_e(\tau)$ array |
| `axes["tau"]` | Free-evolution time axis |
| `metadata["k"]` | Number of echo repetitions |
| `metadata["t_int"]` | Interaction time |
| `metadata["omega_d"]` | Drive frequency |

## TransientSensingExperiment: transient-field sensing

Transient magnetic-field sensing based on sliding measurement: a Ramsey
control pulse (`tau=0`) is slid across the flux signal delay by delay, measuring
$p_e$ at each delay; subtracting a zero-flux reference gives $\Delta p$, and the
control kernel for deconvolution is computed alongside.

**Construction**

`TransientSensingExperiment(qubit, flux_signal=None, flux_signal_zero=None, t_rabi=None, omega_d=None, rotation_angle=π/2, rabi_rate=None, envelope="square", envelope_sigma=None, phase1=π/2, phase2=0.0, scan_list=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test | — |
| `flux_signal` | `FluxSignal` | Applied flux signal | Asymmetric impulse |
| `flux_signal_zero` | `FluxSignal` | Zero-flux reference signal | Constant zero |
| `t_rabi` | `np.ndarray` | Rabi time axis (ns) | `CONFIG.pulse.t_rabi` |
| `omega_d` | float | Drive frequency (rad·GHz) | Qubit sweet-spot frequency |
| `rotation_angle` | float or None | Target angle of each control pulse (rad) | $\pi/2$ |
| `rabi_rate` | float or None | Fixed peak Rabi rate; the envelope integral sets the angle | `None` |
| `envelope` | str or array | Square, Gaussian, or custom envelope | `"square"` |
| `envelope_sigma` | float or None | Gaussian standard deviation (ns) | One quarter of pulse duration |
| `phase1`, `phase2` | float | Rotation-axis phases of the two pulses (rad) | $\pi/2,0$ |
| `scan_list` | `np.ndarray` | Delay scan axis (ns) | Auto-generated |

**Methods**

- `run() -> ExperimentResult`: calls
  {py:class}`~sqc.simulation.SlidingMeasurementRunner` once for the signal and
  once for the zero reference, then computes the control kernel via the
  {doc}`reconstruction` layer's `KernelEstimator`.

**Output**

Returns `ExperimentResult` with:

| Key | Meaning |
|---|---|
| `data["p_e"]` | Excited-state population per delay |
| `data["delta_p"]` | $\Delta p$ (difference from zero reference) |
| `data["kernel"]` | Control kernel array |
| `axes["scan"]` | Delay scan axis |
| `axes["t_samples"]` | Kernel sample time axis |

$\Delta p$ and the kernel are the inputs to the deconvolution in waveform
reconstruction.
`metadata["rotation_angle"]` records the angle obtained by integrating the
discrete envelope. The experiment re-estimates the control kernel for the
constructed pulse, so kernels must not be reused across angle, duration, phase,
or envelope changes.

## CryoscopeExperiment: Cryoscope

The Cryoscope protocol: scan the truncation delay $t_d$, truncating and
coupling the flux signal into the qubit at each truncation, and read the
accumulated phase $\varphi(t_d)$ via IQ Ramsey; differentiating $\varphi$ with
respect to $t_d$ gives the instantaneous frequency, which is then inverted
through a calibration curve to recover $\Phi(t)$.

**Construction**

`CryoscopeExperiment(qubit, flux_signal=None, t_rabi=None, tau=20.0, trunc_list=None, omega_d=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test | — |
| `flux_signal` | `FluxSignal` | Applied flux signal | Sinusoidal signal |
| `t_rabi` | `np.ndarray` | Rabi time axis (ns) | `CONFIG.pulse.t_rabi` |
| `tau` | float | Free-precession time for IQ readout (ns) | `20.0` |
| `trunc_list` | `np.ndarray` | Truncation times (reverse order) | — |
| `omega_d` | float | Drive frequency (rad·GHz) | Qubit sweet-spot frequency |

```{note}
Any value in `trunc_list` beyond `flux_signal.t_list[-1]` raises a `ValueError`
at construction (to prevent silent truncation failures). When lengthening
`trunc_list`, remember to lengthen the `flux_signal` time window too.
```

**Methods**

- `run() -> ExperimentResult`: uses {py:class}`~sqc.hardware.IQReadoutModel` to
  measure the I/Q projections at each truncation, computes the raw phase
  $\varphi=\arctan2(0.5-p_{e,I},\,p_{e,Q}-0.5)$, then applies a model-guided
  unwrap (anchored to the theoretical accumulated phase
  $\int_0^{t_d}(\omega_q(\Phi(t))-\omega_d)\,\mathrm{d}t$) to get a continuous
  $\varphi(t_d)$.

**Output**

Returns `ExperimentResult` with:

| Key | Meaning |
|---|---|
| `data["varphi"]` | Unwrapped accumulated phase $\varphi(t_d)$ |
| `data["p_e_I"]` | I-channel excited-state population |
| `data["p_e_Q"]` | Q-channel excited-state population |
| `axes["trunc"]` | Truncation time axis (time-ascending) |

## DelayRamseyExperiment: delay Ramsey

Delay Ramsey (Ramsey tomography), used to measure flux-pulse tails (the
predistortion mainline): at a delay $t_d$ after the falling edge, slide a short
Ramsey, read the phase $\varphi(t_d)$ via IQ demodulation, and reconstruct the
tail waveform.

**Construction**

`DelayRamseyExperiment(qubit, flux_signal=None, t_d_list=None, tau_R=None, t_rabi=None, omega_d=None, t_fall=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test (biased at a flux-sensitive point) | — |
| `flux_signal` | `FluxSignal` | Flux signal with tail | Square + exponential tail |
| `t_d_list` | `np.ndarray` | Delays after falling edge (ns) | — |
| `tau_R` | float | Ramsey free-evolution time (ns) | — |
| `t_rabi` | `np.ndarray` | Rabi time axis (ns) | `CONFIG.pulse.t_rabi` |
| `omega_d` | float | Drive frequency (rad·GHz) | Qubit sweet-spot frequency |
| `t_fall` | float | Falling-edge time (ns), $t_d$ is relative to it | — |

```{note}
The `run_baseline` field is deprecated: the absolute-phase reference is now set
analytically by the model-guided unwrap (see {doc}`reconstruction`), so no
baseline measurement is needed; the field is kept only for API compatibility.
```

**Methods**

- `run() -> ExperimentResult`: windows the tail segment over the free-evolution
  interval per $t_d$ (zeroed during the two $\pi/2$ pulses to avoid detuning),
  measures the phase via IQ readout, then applies the same model-guided unwrap
  shared with calibration to get $\varphi(t_d)$.

**Output**

Returns `ExperimentResult` with:

| Key | Meaning |
|---|---|
| `data["varphi"]` | Unwrapped phase $\varphi(t_d)$ |
| `data["varphi_raw"]` | Raw (pre-unwrap) phase |
| `data["p_e_I"]` | I-channel excited-state population |
| `data["p_e_Q"]` | Q-channel excited-state population |
| `axes["t_d"]` | Delay time axis |

## PiPulseCompensationExperiment: pi-pulse compensation

The pi-pulse compensation protocol, also used to measure flux-pulse tails:
at a delay $\tau$ after the falling edge, apply a compensation flux pulse of
height $z$ together with a $\pi$ pulse; when $z$ exactly cancels the tail, the
qubit returns to resonance and flips to $|1\rangle$ ($p_e=1$). A 2D scan over
$(\tau, z)$, taking $z^*(\tau)=\arg\max_z p_e$, gives the tail waveform.

**Construction**

`PiPulseCompensationExperiment(qubit, flux_signal=None, tau_list=None, z_list=None, T_pi=None, omega_bias=None, t_rabi=None, t_fall=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test | — |
| `flux_signal` | `FluxSignal` | Flux signal (with falling-edge tail) | — |
| `tau_list` | `np.ndarray` | Delays after falling edge (ns) | — |
| `z_list` | `np.ndarray` | Compensation flux heights ($\Phi_0$) | `linspace(-0.01, 0.01, 21)` |
| `T_pi` | float | $\pi$ pulse duration (ns) | — |
| `omega_bias` | float | Bias frequency (rad·GHz) | — |
| `t_rabi` | `np.ndarray` | Rabi time axis (ns) | `CONFIG.pulse.t_rabi` |
| `t_fall` | float | Falling-edge time (ns) | — |

**Methods**

- `run() -> ExperimentResult`: runs `mesolve` at each point of the
  $(\tau, z)$ grid to get $p_e(\tau, z)$, then takes the peak per $\tau$ via
  parabolic interpolation (avoiding the staircase effect of a coarse `z_list`)
  to get a sub-resolution $z^*(\tau)$.

**Output**

Returns `ExperimentResult` with:

| Key | Meaning |
|---|---|
| `data["p_e"]` | 2D array `(n_tau, n_z)`, $p_e$ at each grid point |
| `data["z_star"]` | 1D array, optimal compensation height (tail waveform) |
| `axes["tau"]` | Delay time axis |
| `axes["z"]` | Compensation height axis |

## Minimal example

```python
import numpy as np
from sqc.devices import TransmonQubit
from sqc.experiments import RabiExperiment, RamseyExperiment

qubit = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10_000, T2=8_000)

# 1) Rabi: calibrate the π pulse — returns a raw QuTiP Result
rabi = RabiExperiment(qubit=qubit).run()
p_e_t = rabi.expect[0]                        # p_e(t)

# 2) Ramsey: measure free-evolution accumulated phase — returns an ExperimentResult
ramsey = RamseyExperiment(
    qubit=qubit,
    tau_list=np.array([0.0, 50.0, 100.0]),    # short sweep to speed up the demo
).run()
print(ramsey.data["p_e"])                     # p_e(τ)
print(ramsey.axes["tau"])                     # free-evolution time axis
```

## Physical role / extension

- Each experiment class corresponds to a complete lab pulse sequence and
  parameter sweep: Rabi calibrates gate durations, Ramsey/echo do phase
  interferometry, and Cryoscope/delay-Ramsey/pi-pulse-compensation measure
  waveforms and tails.
- The experiments layer handles orchestration, not physics or numerics: it
  calls {doc}`control` for pulses, {doc}`simulation` for evolution, and the
  {doc}`hardware` readout models for signals, and only strings these together
  per protocol and collects the raw data.
- The `ExperimentResult` it produces is the uniform input to the
  {doc}`reconstruction` / {doc}`calibration` layers.
- To add a custom sensing protocol, subclass the `Experiment` abstract base,
  implement `build_sequence()` and `run()` (returning an `ExperimentResult`),
  and, when injecting transfer-function distortion, declare a `control_line`
  field and use `_route_flux()`. Full extension guide in {doc}`../extending`.
