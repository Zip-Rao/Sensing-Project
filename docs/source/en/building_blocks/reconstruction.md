# Reconstruction

## What this layer provides

The `reconstruction` layer is the stack's **inverse-problem layer** — it
**inverts** the raw measurement data produced by the {doc}`experiments` layer
(population $p_e$, phase $\varphi$, control kernels, etc.) back into the physical
flux waveform $\Phi(t)$ or magnetic field $B(t)$. It is the algorithmic core of
the waveform-reconstruction mainline: inferring "what was applied to the outside
world" from "what the qubit saw".

A **pure-function constraint** runs through this layer: `reconstruct()` **only
reads measurement data and never runs a new simulation**. If an algorithm needs
forward simulation (e.g. LM full-density-matrix inversion), the required
`qubit` / `control_pulse` are injected **at construction**, not built on the fly
inside `reconstruct()`. The output is uniformly a physical-unit
{py:class}`~sqc.control.FluxSignal` (or a $B$ as an `np.ndarray`).

```{note}
This layer contains two calibration classes, `CryoscopeCalibration` /
`DelayRamseyCalibration`. They are exported under the `sqc.reconstruction`
namespace (paired with their reconstruction classes) but conceptually belong to
**calibration** — they subclass {py:class}`~sqc.calibration.Calibration` and
produce a `CalibrationTable` for the `inversion="calibration"` lookup path. See
{doc}`calibration`.
```

## Class overview

Grouped by function into four sets:

**Extension point**

| Class | Role |
|---|---|
| `Reconstruction` | Abstract base, the common contract for all reconstructors; this layer's extension point |

**Basis-function tools** (for parameterisation, e.g. LM inversion)

| Function | Role |
|---|---|
| `generate_basis_functions` | Generate a list of B-spline / Fourier / Legendre basis functions |
| `basis_function_decomposition` | Least-squares decomposition of a signal onto basis functions |
| `regularization_matrix` / `R` | Smoothness-encouraging regularization matrix (`R` is an alias) |

**Control-kernel estimation** (for transient deconvolution)

| Class | Role |
|---|---|
| `KernelEstimator` | Probe point by point with a narrow Gaussian stimulus, estimating linear/higher-order Volterra kernels |
| `KernelResult` | Kernel-estimation result container (`k1` linear kernel, higher-order kernels, `save`/`load`) |

**Per-protocol reconstructors** (each consumes the `ExperimentResult` of its experiment)

| Class | Experiment | Inversion formula/method |
|---|---|---|
| `RamseyReconstruction` | Ramsey | IQ / unwrap phase → $B=\dot\varphi/\kappa$ |
| `EchoReconstruction` | Differential echo | $B=-\varphi/(2k\kappa t_\mathrm{int})$ |
| `TransientReconstruction` | Transient-field sensing | Wiener / Hammerstein / LM deconvolution |
| `CryoscopeReconstruction` | Cryoscope | $\dot\varphi\to h(t)$, lookup table or dispersion inversion |
| `DelayRamseyReconstruction` | Delay Ramsey | $\varphi(t_d)\to\Phi_\mathrm{tail}$ |
| `PiPulseCompReconstruction` | Pi-pulse compensation | $\Phi_\mathrm{tail}=-z^*(\tau)$ |
| `CryoscopeCalibration` | (calibration) | Square-pulse scan of $h$ to build $\varphi(h)$ lookup |
| `DelayRamseyCalibration` | (calibration) | Scan $z$, fit $\varphi_\mathrm{cal}(z)=\tau_R\kappa z$ |

## Reconstruction — reconstruction abstract base class

The common contract for all reconstructors, this layer's **extension point**. It
mandates a single abstract method:

- `reconstruct(*args, **kwargs)` — invert the physical signal from measurement
  data.

Two invariants (which any subclass must uphold): **(1) `reconstruct()` must be a
pure function** — it reads measurements and runs no new simulation; algorithms
needing forward simulation (like LM) inject `qubit`/`control_pulse` at
construction. **(2) The output is always a physical-unit `FluxSignal`** (or a
field as an `np.ndarray`).

To add a custom reconstruction algorithm, subclass `Reconstruction` and
implement `reconstruct()`. See {doc}`../extending`.

## Basis-function tools

A set of pure functions that parameterise a signal onto basis functions, used
for low-dimensional parameterisation and regularization in LM inversion and the
like.

- `generate_basis_functions(basis_type, n_basis, t_min, t_max) -> list[callable]`
  — generate `n_basis` basis functions on $[t_\min, t_\max]$. `basis_type` is
  `"bspline"` (cubic B-spline), `"fourier"` (orthonormal sine/cosine), or
  `"legendre"` (Legendre polynomials). Returns a list of callables, each taking a
  time array and returning basis values.
- `basis_function_decomposition(sig, t_array, basis_functions) -> np.ndarray`
  — least-squares decompose the signal `sig` (an array, or a `FluxSignal` with a
  `.signal` attribute) onto the basis, returning `(n_basis,)` coefficients.
- `regularization_matrix(n, basis_type) -> np.ndarray` — a smoothness-encouraging
  $n\times n$ regularization matrix: the Fourier basis penalises high
  frequencies by squared frequency order, other bases use a second-difference
  approximation of the second derivative. `R` is its backward-compatible alias.

## Control-kernel estimation

### KernelEstimator

Estimates the **control kernel** needed for transient deconvolution: at each time
point it perturbs the control pulse with a narrow Gaussian flux/frequency
stimulus and measures the qubit response, thereby estimating the linear kernel
$k_1$ and up to higher-order Volterra kernels. `@dataclass`, key fields:

- `mode` — stimulus channel: `"flux"` or `"omega"` (frequency).
- `method` — `"exp"` (experimental: run a simulation to measure the response) or
  `"sim"`.
- `order` — highest kernel order (default 1, linear kernel only; $\ge 2$
  estimates higher-order Volterra kernels).
- `extract_off_diagonal` — whether to extract off-diagonal (n-D) higher-order
  kernels (consumable only by LM inversion).
- Others: `n_levels`, `anharmonicity`, `kappa`, `stim_amplitude`, `stim_width`
  and other probe parameters.

**Methods**: `estimate(pulse, qubit=None, t_samples=None) -> (t_samples, k1)` —
returns the time axis and the linear kernel; `estimate_full(...) -> KernelResult`
— returns the full result object (including higher-order kernels).

### KernelResult

The kernel-estimation result container. Fields: `t_samples` (kernel sample time
axis), `kernels` (list of kernels per order), `mode`, `method`, `order`,
`stim_amplitude`, `units`, `off_diagonal`. Property `k1` (linear kernel =
`kernels[0]`); methods `save(path)`/`load(path)`. Can be passed directly to
`TransientReconstruction.reconstruct(measurement, kernel=...)`.

## RamseyReconstruction — Ramsey reconstruction

Inverts the magnetic field $B(\tau)$ from the Ramsey phase. `@dataclass`, fields:
`qubit`, `method` (`"iq"` or `"unwrap"`, default `"unwrap"`), `k_span` (unwrap
branch-search window, default 3).

- `method="iq"` — dual-channel $\arctan2$ phase extraction (needs
  `data["p_e_I"]`/`p_e_Q`).
- `method="unwrap"` — single-channel $\arccos$ + $k$-span branch unwrapping
  (needs `data["p_e"]`).

Both paths take the gradient of the phase and divide by the frequency
sensitivity $\kappa=\mathrm{d}\omega/\mathrm{d}\Phi$ to get $B=\dot\varphi/\kappa$.
`reconstruct(measurement) -> np.ndarray`.

## EchoReconstruction — differential echo reconstruction

Inverts the field analytically from the differential echo $p_e$:
$\varphi=\arcsin(2p_e-1)$, $B=-\varphi/(2k\kappa t_\mathrm{int})$. `@dataclass`,
fields: `qubit`, `t_int` (interaction time per echo block), `k` (number of $\pi$
pulse pairs). `reconstruct(measurement) -> np.ndarray` (reads `data["p_e"]`).

```{note}
A uniform `clip` to $[-1,1]$ precedes `arcsin`/`arccos`: `mesolve` integration
noise can push $p_e$ marginally outside $[0,1]$, which without clipping produces
silent NaNs. The Ramsey/Echo/Transient paths all have the same guard.
```

## TransientReconstruction — transient-field reconstruction

The **core deconvolver** of the waveform-reconstruction mainline: it
deconvolves the $\Delta p$ and control kernel of the transient experiment back
into $\Phi(t)$. `@dataclass`, with the key field `method` selecting one of four
algorithms:

- `"wiener"` — linear Wiener deconvolution (FFT domain), field `lambda_reg`
  (regularization).
- `"hammerstein"` — Hammerstein-Wiener nonlinear block model (needs `qubit`).
- `"hammerstein_volterra"` — iterative Hammerstein-Volterra inversion using
  higher-order kernels (needs an `order>=2` kernel; fields `max_volterra_iter`,
  `volterra_tol`).
- `"lm"` — Levenberg-Marquardt full-density-matrix inversion (needs
  `control_pulse`; fields `basis_type`, `n_basis`, `max_iter`, `tol`, `mu_init`,
  `use_adjoint`).

`reconstruct(measurement, kernel=None, **kwargs)`: `kernel` can be an
`np.ndarray` (linear kernel) or a `KernelResult` (with higher-order kernels).
Off-diagonal n-D kernels are consumable only by `method="lm"`; the
Wiener/Hammerstein paths accept diagonal (1-D) kernels only, else they raise a
`ValueError`.

## CryoscopeReconstruction — Cryoscope reconstruction

Reconstructs the waveform from phase-vs-truncation data; the core physics is
$\mathrm{d}\varphi/\mathrm{d}t=2\pi\Delta f(h(t))$. `@dataclass`, fields: `tau`
(calibration square-pulse length), `inversion`, `calibration`, `qubit`,
`use_sg_filter`/`sg_window`/`sg_poly` (Savitzky-Golay pre-smoothing of the phase
derivative). Two inversions:

- `inversion="calibration"` — invert via the $\varphi(h)$ lookup table produced
  by `CryoscopeCalibration` (needs `calibration`).
- `inversion="response"` — invert via the analytical Transmon frequency-flux
  dispersion (needs `qubit`).

`reconstruct(measurement, ...) -> FluxSignal` (reads `data["varphi"]`,
`axes["trunc"]`).

## DelayRamseyReconstruction — delay Ramsey reconstruction

Inverts the delay-Ramsey phase $\varphi(t_d)$ into the tail flux waveform (the
predistortion mainline). `@dataclass`, fields: `inversion` (`"response"`
analytical dispersion / `"calibration"` lookup), `qubit`, `calibration`,
`tau_R`. `reconstruct(measurement, ...) -> FluxSignal` (reads `data["varphi"]`,
`axes["t_d"]`).

## PiPulseCompReconstruction — pi-pulse compensation reconstruction

The most direct reconstruction: once the 2D scan gives the optimal compensation
height $z^*(\tau)$, the tail is simply $\Phi_\mathrm{tail}=-z^*(\tau)$.
`reconstruct(measurement) -> FluxSignal` (reads `data["z_star"]`, `axes["tau"]`).

```{note}
Because the $\pi$ pulse has a finite width $T$, $z^*(\tau)$ measures the moving
average of the tail over $[\tau, \tau+T]$, not the instantaneous value. For
exponential tails this only introduces a constant attenuation factor; the decay
time constant is preserved.
```

## Calibration classes (CryoscopeCalibration / DelayRamseyCalibration)

Two classes exported in this namespace but conceptually belonging to
**calibration** (subclasses of {py:class}`~sqc.calibration.Calibration`), which
build the lookup tables for their reconstructors' `"calibration"` inversion
path:

- `CryoscopeCalibration` — scan square-pulse heights $h$, IQ readout +
  model-guided unwrap to build a $\varphi(h)$ lookup (`calibrate() ->
  CalibrationTable`, `kind="phi_h"`).
- `DelayRamseyCalibration` — scan known flux heights $z$, fit the linear slope
  $\varphi_\mathrm{cal}(z)=\tau_R\kappa z$ (`calibrate() -> CalibrationTable`,
  `kind="phi_z"`).

Full calibration contract in {doc}`calibration`.

## Minimal example

```python
import numpy as np
from sqc.devices import TransmonQubit
from sqc.simulation import ExperimentResult
from sqc.reconstruction import (
    RamseyReconstruction, generate_basis_functions, basis_function_decomposition,
)

qubit = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10_000, T2=8_000, flux=0.1)

# 1) Synthesise a Ramsey measurement (normally from RamseyExperiment.run())
tau = np.linspace(0, 250, 60)
p_e = 0.5 * (1 - np.cos(0.02 * tau))
meas = ExperimentResult(data={"p_e": p_e}, axes={"tau": tau})

# 2) Pure post-processing inversion of the field B(τ) — runs no new simulation
B = RamseyReconstruction(qubit=qubit, method="unwrap").reconstruct(meas)

# 3) Decompose the signal onto a B-spline basis (parameterisation for LM, etc.)
basis = generate_basis_functions("bspline", 8, tau[0], tau[-1])
coeffs = basis_function_decomposition(p_e, tau, basis)   # (8,) coefficients
```

## Physical role / extension

- This layer corresponds to the **offline analysis** of experimental data:
  translating the measured population/phase back into the flux/field that was
  actually applied.
- Each reconstructor pairs one-to-one with an experiment class from
  {doc}`experiments`: Ramsey↔Ramsey, Echo↔differential echo,
  Transient↔transient, Cryoscope/delay-Ramsey/pi-compensation↔their tail
  measurements.
- The control kernel produced by `KernelEstimator` is a required input to
  transient deconvolution; the basis-function tools support the low-dimensional
  parameterisation in LM; the two calibration classes supply the
  `CalibrationTable` for the lookup inversion paths.
- **To add a custom reconstruction algorithm**, subclass the `Reconstruction`
  abstract base and implement `reconstruct()`, upholding the two invariants
  ("pure function, output a `FluxSignal`"), and inject `qubit` at construction if
  forward simulation is needed. Full extension guide in {doc}`../extending`.

