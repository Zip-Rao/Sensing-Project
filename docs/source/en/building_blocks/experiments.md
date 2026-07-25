# Experiments

## What this layer provides

The `experiments` layer is the stack's **protocol-orchestration layer** — it
assembles the parts from the layers below (devices, flux signals, control
pulses, solvers, readout models) into a complete, runnable experiment following
some **sensing protocol**. Each experiment class answers a version of the same
question: "given this pulse sequence and this parameter sweep, how does the
qubit's final-state population/phase vary with the parameter?" It does no
numerical integration itself (that goes to `mesolve` in the {doc}`simulation`
layer) and no field reconstruction (that goes to the {doc}`reconstruction`
layer) — it only **orchestrates and produces raw measurement data**.

One convention runs through this layer: **every experiment is a `@dataclass`,
constructed with a device and optional parameters, and `run()` executes the
whole thing in one call and returns the result**. Almost all experiments return
the uniform {py:class}`~sqc.simulation.ExperimentResult` (the four dicts
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

| Class | Sensing protocol | Key outputs |
|---|---|---|
| `Experiment` | Abstract base (this layer's extension point) | — |
| `RabiExperiment` | Rabi oscillation: constant-amplitude drive, sweep duration | `p_e(t)` (raw QuTiP `Result`) |
| `RamseyExperiment` | Ramsey: $\pi/2-\tau-\pi/2$, measures free-evolution accumulated phase | `p_e(\tau)` |
| `DiffEchoExperiment` | Differential echo: $k$ repetitions accumulate phase, enhancing weak signals | `p_e(\tau)` |
| `TransientSensingExperiment` | Transient-field sensing: sliding measurement + control kernel | `delta_p`, `kernel` |
| `CryoscopeExperiment` | Cryoscope: scan truncation delay, IQ readout of $\varphi(t_d)$ | `varphi(t_d)` |
| `DelayRamseyExperiment` | Delay Ramsey: slide a short Ramsey across the falling edge to measure the tail | `varphi(t_d)` |
| `PiPulseCompensationExperiment` | Pi-pulse compensation: 2D scan of $(\tau, z)$ to measure the tail | `z_star(\tau)` |

These experiments map onto three product mainlines: **waveform reconstruction**
(Ramsey / DiffEcho / Transient / Cryoscope recovering $\Phi(t)$), **frequency
calibration** (Ramsey `f(\Phi)`), and **predistortion** (DelayRamsey /
PiPulseComp measuring the AWG→chip tail). Their end-to-end pipelines are in
{doc}`../examples/waveform_reconstruction`,
{doc}`../examples/frequency_calibration`, and
{doc}`../examples/predistortion`.

## Experiment — experiment abstract base class

The common contract for all sensing experiments, and this layer's **extension
point**. It mandates two abstract methods:

- `build_sequence()` — build the control pulse sequence for this experiment
  (some experiments build it per-parameter, in which case this returns `None`).
- `run(*args, **kwargs)` — execute the whole experiment and return the result.

It also provides one shared helper, `_route_flux(signal)`: pass a flux signal
through an optional {py:class}`~sqc.hardware.ControlLine` distortion (the
`control_line` field) before coupling it to the qubit; with `control_line=None`
(the default) it passes through. This is the predistortion mainline's uniform
entry point for injecting transfer-function distortion.

To add a custom sensing protocol, subclass `Experiment` and implement
`build_sequence()` / `run()`; `run()` should return an `ExperimentResult`. See
{doc}`../extending`.

## RabiExperiment — Rabi oscillation

Applies a **constant-amplitude** drive pulse and measures the excited-state
population versus time, used to calibrate $\pi$/$\pi/2$ pulse durations.
`@dataclass`, main fields: `qubit`, `t_rabi` (Rabi time axis, default
`make_time(0, 40)`), `omega_d` (drive frequency, default the qubit sweet-spot
frequency).

`build_sequence()` returns a constant-amplitude {py:class}`~sqc.control.Pulse`
(rotating frame, RWA). `run()` runs `mesolve` on the global time axis
`CONFIG.pulse.t_global` and **returns a raw QuTiP `Result`** (see the opening
note); take `result.expect[0]` for $p_e(t)$.

## RamseyExperiment — Ramsey interferometry

A $\pi/2-\tau-\pi/2$ sequence that reflects qubit-frequency shifts caused by an
external flux $\Phi(t)$ through the phase accumulated during free evolution; it
is the core protocol for both the waveform-reconstruction and
frequency-calibration mainlines. `@dataclass`, main fields: `qubit`,
`flux_signal` (default constant signal), `omega_d`, `tau_list` (free-evolution
sweep, default `CONFIG.pulse.tau_list`), `t_global`, `phase1`/`phase2` (the two
$\pi/2$ pulse phases).

`run()` first projects the flux signal onto the global time axis and couples it
into the qubit via `qubit_in_mag` (rotating frame), then builds a Ramsey
sequence per `tau`, runs the evolution, and collects the final $p_e$. In the
returned `ExperimentResult`: `data["p_e"]` is $p_e(\tau)$,
`data["flux_samples"]` is the flux samples, `axes["tau"]` is the free-evolution
axis, and `axes["t_flux"]` is the flux time axis.

## DiffEchoExperiment — differential echo

The differential echo protocol: sequence
$\pi/2-[\tau-\pi-\tau'-(\tau+t_\mathrm{int})-\pi-\tau'']^k-\pi/2$, which enhances
weak-signal sensitivity by accumulating phase over $k$ echo repetitions.
`@dataclass`, main fields: `qubit`, `flux_signal` (default Gaussian), `k` (number
of echo repetitions, default 5), `t_rabi`, `t_int` (interaction time), `t_rep`
(repetition period), `tau_list`, `t_global`, `omega_d`.

`run()` replicates the base flux signal $2k$ times, concatenates and projects it
onto the global axis and couples it into the qubit, then builds a differential
echo pulse per `tau` and runs the evolution. In the returned `ExperimentResult`:
`data["p_e"]` is $p_e(\tau)$, plus `axes["tau"]` and `metadata` recording `k` /
`t_int` / `omega_d`.

## TransientSensingExperiment — transient-field sensing

Transient magnetic-field sensing based on **sliding measurement**: a Ramsey
control pulse (`tau=0`) is slid across the flux signal delay by delay, measuring
$p_e$ at each delay; subtracting a zero-flux reference gives $\Delta p$, and the
control kernel for deconvolution is computed alongside. `@dataclass`, main
fields: `qubit`, `flux_signal` (default asymmetric impulse), `flux_signal_zero`
(zero-flux reference, default constant zero), `t_rabi`, `omega_d`, `scan_list`.

`run()` calls {py:class}`~sqc.simulation.SlidingMeasurementRunner` once for the
signal and once for the zero reference, then computes the control kernel via the
{doc}`reconstruction` layer's `KernelEstimator`. In the returned
`ExperimentResult`: `data["p_e"]`, `data["delta_p"]`, `data["kernel"]`,
`axes["scan"]` (delay axis), `axes["t_samples"]` (kernel sample axis). $\Delta p$
and the kernel are exactly the inputs to the deconvolution in waveform
reconstruction.

## CryoscopeExperiment — Cryoscope

The Cryoscope protocol: scan the **truncation delay** $t_d$, truncating and
coupling the flux signal into the qubit at each truncation, and read the
accumulated phase $\varphi(t_d)$ via IQ Ramsey; differentiating $\varphi$ with
respect to $t_d$ gives the instantaneous frequency, which is then inverted
through a calibration curve to recover $\Phi(t)$. `@dataclass`, main fields:
`qubit`, `flux_signal` (default sinusoidal), `t_rabi`, `tau` (free-precession
time for IQ readout), `trunc_list` (list of truncation times, in reverse order),
`omega_d`.

`run()` uses {py:class}`~sqc.hardware.IQReadoutModel` to measure the I/Q
projections at each truncation, computes the raw phase
$\varphi=\arctan2(0.5-p_{e,I},\,p_{e,Q}-0.5)$, then applies a model-guided
unwrap (anchored to the theoretical accumulated phase
$\int_0^{t_d}(\omega_q(\Phi(t))-\omega_d)\,\mathrm{d}t$) to get a continuous
$\varphi(t_d)$. In the returned `ExperimentResult`: `data["varphi"]`,
`data["p_e_I"]`, `data["p_e_Q"]`, `axes["trunc"]` (time-ascending).

```{note}
Any value in `trunc_list` beyond `flux_signal.t_list[-1]` raises a `ValueError`
at construction (to prevent silent truncation failures). When lengthening
`trunc_list`, remember to lengthen the `flux_signal` time window too.
```

## DelayRamseyExperiment — delay Ramsey

Delay Ramsey (Ramsey tomography), used to **measure flux-pulse tails** (the
predistortion mainline): at a delay $t_d$ after the falling edge, slide a short
Ramsey, read the phase $\varphi(t_d)$ via IQ demodulation, and reconstruct the
tail waveform. `@dataclass`, main fields: `qubit` (biased at a flux-sensitive
point), `flux_signal` (default square wave + exponential tail), `t_d_list`
(delays after the falling edge), `tau_R` (Ramsey free-evolution time), `t_rabi`,
`omega_d`, `t_fall` (falling-edge time, $t_d$ is relative to it).

`run()` windows the tail segment over the free-evolution interval per $t_d$
(zeroed during the two $\pi/2$ pulses to avoid detuning), measures the phase via
IQ readout, then applies the same model-guided unwrap shared with calibration to
get $\varphi(t_d)$. In the returned `ExperimentResult`: `data["varphi"]`,
`data["varphi_raw"]`, `data["p_e_I"]`, `data["p_e_Q"]`, `axes["t_d"]`.

```{note}
The `run_baseline` field is deprecated: the absolute-phase reference is now set
analytically by the model-guided unwrap (see {doc}`reconstruction`), so no
baseline measurement is needed; the field is kept only for API compatibility.
```

## PiPulseCompensationExperiment — pi-pulse compensation

The pi-pulse compensation protocol, also used to **measure flux-pulse tails**:
at a delay $\tau$ after the falling edge, apply a compensation flux pulse of
height $z$ together with a $\pi$ pulse; when $z$ exactly cancels the tail, the
qubit returns to resonance and flips to $|1\rangle$ ($p_e=1$). A 2D scan over
$(\tau, z)$, taking $z^*(\tau)=\arg\max_z p_e$, gives the tail waveform.
`@dataclass`, main fields: `qubit`, `flux_signal`, `tau_list`, `z_list`
(compensation heights, default `linspace(-0.01, 0.01, 21)`), `T_pi`,
`omega_bias`, `t_rabi`, `t_fall`.

`run()` runs `mesolve` at each point of the $(\tau, z)$ grid to get
$p_e(\tau, z)$, then takes the peak per $\tau$ via parabolic interpolation
(avoiding the staircase effect of a coarse `z_list`) to get a sub-resolution
$z^*(\tau)$. In the returned `ExperimentResult`: `data["p_e"]` (2D `(n_tau,
n_z)`), `data["z_star"]` (the tail waveform), `axes["tau"]`, `axes["z"]`.

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

- Each experiment class corresponds to a complete lab **pulse sequence +
  parameter sweep**: Rabi calibrates gate durations, Ramsey/echo do phase
  interferometry, and Cryoscope/delay-Ramsey/pi-pulse-compensation measure
  waveforms and tails.
- The experiments layer is "orchestration", not "physics" or "numerics": it
  calls {doc}`control` for pulses, {doc}`simulation` for evolution, and the
  {doc}`hardware` readout models for signals, and only strings these together
  per protocol and collects the raw data.
- The `ExperimentResult` it produces is the uniform input to the
  {doc}`reconstruction` / {doc}`calibration` layers.
- **To add a custom sensing protocol**, subclass the `Experiment` abstract base,
  implement `build_sequence()` and `run()` (returning an `ExperimentResult`),
  and — when injecting transfer-function distortion — declare a `control_line`
  field and use `_route_flux()`. Full extension guide in {doc}`../extending`.

