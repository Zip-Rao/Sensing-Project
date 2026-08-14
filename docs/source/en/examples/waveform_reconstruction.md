# Waveform Reconstruction

## Overview

Waveform reconstruction is the platform's first product mainline: using a
superconducting transmon qubit as a time-domain magnetic flux sensor to recover
an unknown external flux waveform $\Phi(t)$ from qubit-state measurement data
alone.

Physically, external flux $\Phi(t)$ modulates the Josephson energy through the
SQUID loop and thereby shifts the qubit frequency $\omega_T(\Phi)$. A
Ramsey-type control pulse is **slid** along the time axis; at each delay
$t_d$ the qubit's excited-state population $p_e$ is measured. The coherent phase
accumulated near that delay is proportional to the instantaneous frequency shift
induced by the flux, so $p_e(t_d)$ approximates $\Phi(t)$ smoothed by the pulse
response. This smoothing is described by the **control kernel** $k(t)$—the
pulse's response function to a point flux stimulus.

The measurement and the unknown waveform are related by a convolution:

$$\Delta p(t) \approx (k * \Phi)(t),$$

where $\Delta p$ is the difference between two sliding measurements—one with
signal, one with a zero-flux reference—after subtracting the flux-independent
baseline. **Reconstruction is the inverse problem of this convolution**: given
$\Delta p$ and the kernel $k$, solve for $\Phi(t)$.

## Pipeline Architecture

The waveform-reconstruction pipeline executes in the following order: sensor
modelling → signal and probe construction → sliding measurement → experiment
orchestration → kernel estimation → inverse problem → global orchestration.
The sections below follow this logical sequence.

### Step 1: Model the Sensor — Device Layer (`TransmonQubit`)

The physical foundation of all reconstruction is the qubit itself.
{py:class}`~sqc.devices.transmon.TransmonQubit` holds the charging energy
$E_C$, Josephson energy $E_J$, decoherence times $T_1$/$T_2$, and flux bias
`flux`. The three methods most critical to waveform reconstruction are:

- `frequency()` — returns the $|0\rangle \to |1\rangle$ transition frequency
  $\omega_T(\Phi)$ at the current bias point.
- `frequency_sensitivity(flux)` — numerically computes $d\omega/d\Phi$ via
  central-difference quotient, i.e. the frequency-to-flux response sensitivity
  $\kappa$.
- `qubit_in_mag(FluxSignal)` — pre-computes `freq_coeffs` and `H_list` for use
  by `mesolve()`. **Must be called** whenever the flux signal changes to
  rebuild the time-dependent Hamiltonian.

The flux bias point directly affects sensitivity: at $\Phi=0$ (the sweet spot)
$\kappa=0$, giving zero first-order response; near $\Phi \approx 0.1\,\Phi_0$,
$\kappa$ is maximal, which is the optimal operating point for transient
sensing.

### Step 2: Construct the Flux Signal — Control Layer (`FluxSignal`)

The unknown flux waveform to be recovered is represented by
{py:class}`~sqc.control.flux_signal.FluxSignal`, a physically meaningful
subclass of `Waveform` supporting 8 signal types: zero (`type=0`), constant
(`1`), sinusoidal (`2`), Gaussian (`3`), asymmetric impulse (`4`),
double-peak (`5`), basis-expanded (`6`), complex wave-packet (`7`), and custom
raw samples (`8`). Waveform shape is controlled via `amplitude`, `width`,
`center`, `rise`, and `fall`; the `signal` property returns the `(N,)` sample
array.

In simulation, `TransientSensingExperiment.run()` automatically constructs a
zero-flux reference signal (`type=1`, `amplitude=0`) for baseline subtraction.

### Step 3: Construct the Control Pulse — Control Layer (`Pulse` + `CompositePulse`)

The probe is a Ramsey-type control pulse: two $\pi/2$ pulses separated by zero
free-evolution delay ($\tau=0$—for transient sensing the pulse slides directly
across the signal with no extra wait). Three control-layer abstractions build
this probe:

- {py:class}`~sqc.control.pulse.Pulse`: a single pulse holding I/Q envelopes
  (`Omega`, `Omega_Q`). Supports lab frame (`frame=0`) and rotating frame
  (`frame=1`), with or without the rotating-wave approximation (RWA). Its
  `hamiltonian` property returns a QuTiP list-format $[H_0, c_0(t)]$ pair
  suitable for direct passage to `mesolve`.
- {py:class}`~sqc.control.pulse.CompositePulse`: a series of concatenated
  `Pulse` objects, forming a global representation by stitching together each
  sub-pulse's time axis and Hamiltonian coefficients.
- {py:class}`~sqc.control.sequence.create_ramsey_pulse`: a factory function
  for convenient construction of a standard Ramsey pulse (two $\pi/2$ pulses +
  zero free-evolution delay).

### Step 4: Execute the Sliding Measurement — Simulation Layer (`SlidingMeasurementRunner`)

With the control pulse built,
{py:class}`~sqc.simulation.runner.SlidingMeasurementRunner` slides it point
by point along the signal time axis, performing a full quantum evolution at
each delay $t_d$:

1. Determine the time window for the control pulse at the current delay and
   construct the time-dependent Hamiltonian on that window.
2. Call QuTiP's `mesolve()` to integrate the Lindblad master equation,
   tracking the qubit density-matrix evolution.
3. Extract the excited-state population $p_e(t_d)$ from the final state.
4. Collect $p_e$ across all delays into the measurement curve.

This layer accounts for the bulk of the computational cost—each delay point
corresponds to a full `mesolve` integration.

### Step 5: Orchestrate the Experiment — Experiment Layer (`TransientSensingExperiment`)

{py:class}`~sqc.experiments.TransientSensingExperiment` chains the three steps
above into a complete measurement protocol:

1. Call `create_ramsey_pulse` to build the control pulse.
2. Instantiate `SlidingMeasurementRunner`; run the sliding measurement on the
   flux signal to obtain $p_e^{\text{sig}}$.
3. Repeat on the zero-flux reference signal to obtain $p_e^{\text{ref}}$;
   compute $\Delta p = p_e^{\text{sig}} - p_e^{\text{ref}}$.
4. Call {py:class}`~sqc.reconstruction.kernel.KernelEstimator` to estimate the
   control kernel $k(t)$ (see next step).

The returned {py:class}`~sqc.simulation.result.ExperimentResult` contains
`data["p_e"]`, `data["delta_p"]`, `data["kernel"]`, `data["flux_samples"]`
(the simulation ground truth), and axis arrays `axes["scan"]`,
`axes["t_flux"]`.

### Step 6: Kernel Estimation and Inverse Problem — Reconstruction Layer

The reconstruction layer is the algorithmic core of the inverse problem,
comprising two tightly coupled components: the kernel estimator provides the
forward model's smoothing kernel; the reconstructor performs the inversion.

#### Kernel Estimation (`KernelEstimator`)

{py:class}`~sqc.reconstruction.kernel.KernelEstimator` estimates the pulse's
point-flux response function $k(t)$: at each time point it injects a narrow
Gaussian flux stimulus (`stim_amplitude`, `stim_width`) and measures the
change in qubit excited-state population, building the kernel point by point.
`order=1` estimates only the linear kernel $k_1$; `order >= 2` additionally
estimates higher-order Volterra kernels $k_2, k_3, \ldots$ for nonlinear
inversion. The returned `KernelResult` object can be passed directly to the
reconstructor's `reconstruct()` method.

#### Inverse-Problem Solvers (`TransientReconstruction`)

{py:class}`~sqc.reconstruction.transient.TransientReconstruction` inverts
$\Delta p$ and $k$ to recover $\Phi(t)$, offering four methods ordered by
increasing accuracy and computational cost:

| `method` | Algorithm | Regime |
|---|---|---|
| `"wiener"` | Frequency-domain Wiener deconvolution $X_f = \bar{K}_f Y_f / (|K_f|^2 + \lambda^2)$ | Small signals, linear regime; fastest |
| `"hammerstein"` | Wiener deconvolution + $\omega(\Phi)$ dispersion inverse mapping | Moderate amplitude, off sweet spot |
| `"hammerstein_volterra"` | Iterative Volterra-series inversion using higher-order kernels to subtract nonlinear contributions | Large signals, significant nonlinearity |
| `"lm"` | Levenberg-Marquardt full density-matrix iterative inversion with basis-function parameterisation and adjoint Jacobian | Highest accuracy, highest cost |

The first three methods take $\Delta p$ and $k$ as input and return a
`FluxSignal`. The LM method additionally requires `qubit` and `control_pulse`
injected at construction time, as each iteration must run a forward
density-matrix simulation to compute $p_e^{\text{sim}}$ and the residual
Jacobian.

Complementing this component are basis-function utilities
({py:class}`~sqc.reconstruction.basis.generate_basis_functions`, etc.) that
provide low-dimensional parameterisation on B-spline / Fourier / Legendre
bases and smoothness regularisation matrices for LM inversion.

```{note}
Non-diagonal higher-order kernels (produced by
`KernelEstimator(extract_off_diagonal=True)`) can only be consumed by
`method="lm"`. The Wiener / Hammerstein paths accept diagonal (1-D) kernels
only; passing the wrong kind raises `ValueError`.
```

### Step 7: Global Time Axes — Configuration Layer (`CONFIG`)

{py:class}`sqc.config.CONFIG` is the single source of truth for all global
parameters. Every time axis is derived from `CONFIG.awg.dt` (the sampling time
step, $=1/\text{sample\_rate}$) or `CONFIG.pulse.make_time(start, end)`,
guaranteeing a consistent discretisation from pulse construction through
sliding measurement to reconstruction. Parameters modified by
`SensingWorkflow.configure()` are synchronised to `CONFIG` on the fly via
`_sync_config()`.

### Step 8: End-to-End Orchestration — Workflow Layer (`SensingWorkflow`)

{py:class}`~sqc.workflows.SensingWorkflow` wraps the seven steps above into a
`configure()` $\to$ `run()` $\to$ `plot()` interface—the single user-facing
entry point. It implements no physics itself but handles:

- **Configuration management**: stores parameters grouped by functional area
  (qubit, signal, reconstruction, pulse, hardware); `configure()` only mutates
  explicitly passed fields.
- **Pipeline execution**: `run()` sequentially calls `_build_qubit()` →
  `_build_flux_signal()` → `_make_experiment().run()` →
  `_run_reconstruction()`.
- **Batch operations**: `sweep("signal.amplitude", [...])` scans a single
  parameter with automatic CONFIG synchronisation;
  `compare(methods=[...])` runs multiple reconstruction algorithms on the
  same $\Delta p$, selects the best by RMSE, and does **not** re-run
  `mesolve`.
- **Visualisation**: `plot()` auto-dispatches to the appropriate figure type
  ($p_e$/$\Delta p$ measurement, reconstructed waveform vs truth, sweep
  curves, multi-method overlay + residual) based on the most recently executed
  operation.

## Usage

### End-to-End Pipeline (Recommended)

```python
from sqc.workflows import SensingWorkflow

# Configure
wf = SensingWorkflow()
wf.configure(
    protocol="transient",       # transient field sensing (sliding Ramsey)
    flux_bias=0.1,              # bias to finite κ for sensitivity
    signal_type=4,              # waveform to recover: asymmetric impulse
    signal_amplitude=0.02,      # peak flux (Φ₀)
    signal_center=100,          # impulse centre (ns)
    reconstruction="wiener",    # reconstruction algorithm
    lambda_reg=5.0,             # Wiener regularisation parameter
)

# Execute: measure + reconstruct
result = wf.run()
print("recovered points:", len(result.reconstructed_signal.signal))

# Compare three algorithms on the same measurement
cmp = wf.compare(methods=["wiener", "hammerstein", "lm"])
print("best method:", cmp.best)
for m in cmp.methods:
    print(f"  {m}: rmse={cmp.metrics[m]['rmse']:.3e}")

# Plot
wf.plot()
```

```{note}
`compare()` reuses the measurement data cached by `run()`, swapping only the
reconstruction algorithm; it does not re-run `mesolve`, so all methods are
compared on the **same** $\Delta p$. `best` is selected by RMSE; because this
is a simulation, the ground-truth `flux_samples` is known, which is what makes
RMSE a meaningful metric.
```

### Step-by-Step Pipeline

For fine-grained control, operate directly on each layer in logical pipeline
order:

```python
import numpy as np
from sqc.devices import TransmonQubit
from sqc.control.flux_signal import FluxSignal
from sqc.experiments.transient import TransientSensingExperiment
from sqc.reconstruction.transient import TransientReconstruction

# 1. Construct qubit (sensor)
qubit = TransmonQubit(
    EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
    T1=10_000, T2=8_000, flux=0.1,
)

# 2. Construct flux signal (the unknown)
t_list = np.arange(0, 400, 0.5)  # dt = 0.5 ns
flux_signal = FluxSignal(
    type=4, t_list=t_list,
    amplitude=0.02, center=100, rise=10, fall=10,
)

# 3–5. Run experiment (pulse + sliding measurement + Δp + kernel)
exp = TransientSensingExperiment(qubit=qubit, flux_signal=flux_signal)
meas = exp.run()

# 6. Solve the inverse problem
rec = TransientReconstruction(method="wiener", lambda_reg=5.0)
recovered = rec.reconstruct(meas, kernel=meas.data["kernel"])
print(recovered.signal[:5])
```

## Reading the Results

- `result.measurement.data` holds four key arrays: `p_e` (excited-state
  population), `delta_p` (baseline-subtracted signal difference), `kernel`
  (the control kernel $k$), and `flux_samples` (the ground-truth waveform
  injected in simulation; available only in simulation). Axis arrays are in
  `result.measurement.axes`: `scan` (sliding-delay axis), `t_flux` (waveform
  time axis).
- `result.reconstructed_signal` is the recovered
  {py:class}`~sqc.control.flux_signal.FluxSignal`. `wf.plot()` overlays it on
  the ground truth for direct visual quality assessment.
- The `metrics` from `compare()` give each method `rmse`, `snr`, and `peak`.
  Typical takeaway: Wiener is fastest and sufficient for small signals;
  Hammerstein accounts for the weak $\omega(\Phi)$ nonlinearity and lowers
  RMSE at large amplitudes; LM is the most accurate but most expensive (each
  iteration requires a full forward density-matrix simulation).
- The physical differences between the three algorithms are detailed in
  {doc}`../building_blocks/reconstruction` § TransientReconstruction.
- For parameter sweeping, use `wf.sweep("signal.amplitude", [...])` to obtain
  the SNR-vs-amplitude curve. See {doc}`../building_blocks/workflows` for the
  full API.
