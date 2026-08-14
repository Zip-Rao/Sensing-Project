# Reconstruction

## What this layer provides

The `reconstruction` layer is the stack's inverse-problem layer. It
inverts the raw measurement data produced by the {doc}`experiments` layer
(population $p_e$, phase $\varphi$, control kernels, etc.) back into the physical
flux waveform $\Phi(t)$ or magnetic field $B(t)$. It is the algorithmic core of
the waveform-reconstruction mainline: inferring the applied external signal from
the qubit's measured response.

A pure-function constraint runs through this layer: `reconstruct()` only
reads measurement data and never runs a new simulation. If an algorithm needs
forward simulation (e.g. LM full-density-matrix inversion), the required
`qubit` / `control_pulse` are injected at construction, not built on the fly
inside `reconstruct()`. The output is uniformly a physical-unit
{py:class}`~sqc.control.FluxSignal` (or a $B$ as an `np.ndarray`).

```{note}
This layer contains two calibration classes, `CryoscopeCalibration` /
`DelayRamseyCalibration`. They are exported under the `sqc.reconstruction`
namespace (paired with their reconstruction classes) but conceptually belong to
calibration: they subclass {py:class}`~sqc.calibration.Calibration` and
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

| Class | Role | Inversion formula/method |
|---|---|---|
| `RamseyReconstruction` | Ramsey reconstruction | IQ / unwrap phase → $B=\dot\varphi/\kappa$ |
| `EchoReconstruction` | Differential echo reconstruction | $B=-\varphi/(2k\kappa t_\mathrm{int})$ |
| `TransientReconstruction` | Transient-field reconstruction | Wiener / Hammerstein / LM deconvolution |
| `CryoscopeReconstruction` | Cryoscope reconstruction | $\dot\varphi\to h(t)$, lookup table or dispersion inversion |
| `DelayRamseyReconstruction` | Delay Ramsey reconstruction | $\varphi(t_d)\to\Phi_\mathrm{tail}$ |
| `PiPulseCompReconstruction` | Pi-pulse compensation reconstruction | $\Phi_\mathrm{tail}=-z^*(\tau)$ |
| `CryoscopeCalibration` | (calibration) | Square-pulse scan of $h$ to build $\varphi(h)$ lookup |
| `DelayRamseyCalibration` | (calibration) | Scan $z$, fit $\varphi_\mathrm{cal}(z)=\tau_R\kappa z$ |

## Reconstruction: reconstruction abstract base class

The common contract for all reconstructors, this layer's extension point. It
mandates a single abstract method:

- `reconstruct(*args, **kwargs)`: invert the physical signal from measurement
  data.

Two invariants (which any subclass must uphold):

1. `reconstruct()` must be a pure function, reading measurements and running no
   new simulation; algorithms needing forward simulation (like LM) inject
   `qubit`/`control_pulse` at construction.
2. The output is always a physical-unit `FluxSignal` (or a field as an
   `np.ndarray`).

To add a custom reconstruction algorithm, subclass `Reconstruction` and
implement `reconstruct()`. See {doc}`../extending`.

## Basis-function tools

A set of pure functions that parameterise a signal onto basis functions, used
for low-dimensional parameterisation and regularization in LM inversion and the
like.

- `generate_basis_functions(basis_type, n_basis, t_min, t_max) -> list[callable]`:
  generate `n_basis` basis functions on $[t_\min, t_\max]$. `basis_type` is
  `"bspline"` (cubic B-spline), `"fourier"` (orthonormal sine/cosine), or
  `"legendre"` (Legendre polynomials). Returns a list of callables, each taking a
  time array and returning basis values.
- `basis_function_decomposition(sig, t_array, basis_functions) -> np.ndarray`:
  least-squares decompose the signal `sig` (an array, or a `FluxSignal` with a
  `.signal` attribute) onto the basis, returning `(n_basis,)` coefficients.
- `regularization_matrix(n, basis_type) -> np.ndarray`: a smoothness-encouraging
  $n\times n$ regularization matrix. The Fourier basis penalises high
  frequencies by squared frequency order, other bases use a second-difference
  approximation of the second derivative. `R` is its backward-compatible alias.

## Control-kernel estimation

The control kernel is the core input to transient deconvolution: it quantifies
"how much a tiny perturbation at time $t$ during the control pulse changes the
final qubit population $p_e$".  `KernelEstimator` supports a three-dimensional
design space on this perturbation→response map, covering everything from a fast
first-order linear kernel to the complete third-order off-diagonal Volterra
tensor.

### KernelEstimator

**Construction**

`KernelEstimator(mode="omega", method="exp", order=1, virtual_z_impl="math",
extract_off_diagonal=False, ...)`

**All fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `mode` | `str` | Perturbation channel: `"flux"` (magnetic) or `"omega"` (Virtual Z frequency) | `"flux"` |
| `method` | `str` | Estimation strategy: `"exp"` (experimental sim, with qubit dispersion) or `"sim"` (pure-theory Heisenberg, no qubit dispersion) | `"exp"` |
| `order` | `int` | Highest Volterra order (1 = linear only; ≥2 = higher-order) | `1` |
| `extract_off_diagonal` | `bool` | Extract full n-D off-diagonal tensors (consumable only by LM) | `False` |
| `virtual_z_impl` | `str` | VZ implementation: `"math"` (σ_z Gaussian impulse) or `"hardware"` (phase-shifted sub-pulses; CompositePulse only) | `"math"` |
| `stim_amplitude` | `float` | Probe stimulus amplitude (Φ₀ for flux mode, rad for omega mode) | `0.0215` |
| `stim_width` | `float` | Probe width (ns, flux mode only) | `3.0` |
| `auto_calibrate` | `bool` | Auto-adjust amplitude per qubit anharmonicity (flux mode) | `False` |
| `n_levels` | `int` | Hilbert-space truncation (sim mode; exp auto-detects from qubit) | `2` |
| `anharmonicity` | `float` | Anharmonicity α (GHz), used when n_levels≥3 and method="sim" | `0.0` |
| `kappa` | `float` | Linear dispersion dω/dΦ (GHz/Φ₀), only for sim-mode flux→omega conversion | `None` |
| `n_amp_samples` | `int` | Number of amplitude scan points (deprecated; fixed 5-point FD stencil now used) | `5` |
| `amp_scan_factor` | `float` | Amplitude scan range multiplier | `1.0` |
| `probe_sigma_t` | `float` or `None` | VZ Gaussian probe width σ_t (ns); `None`→`2·dt`. Smaller = less off-diagonal smear bias, but must be ≥dt | `None` |
| `richardson` | `bool` | Enable Richardson σ_t→0 extrapolation to remove probe-width bias (omega/math only) | `False` |
| `richardson_sigmas` | `tuple` or `None` | σ_t sample points (×dt) for Richardson; `None`→`(2.0, 1.5, 1.0)` | `None` |

**Methods**

- `estimate(pulse, qubit=None, t_samples=None) -> (t_samples, k1)`:
  returns the time axis and first-order linear kernel $k_1$ (backward-compat).
- `estimate_full(pulse, qubit=None, t_samples=None) -> KernelResult`:
  returns the full result object (including higher-order kernels).

### Design space: complete dispatch table

The three-dimensional design space is formed by `(mode, method, order)` plus the
`extract_off_diagonal` and `virtual_z_impl` flags.  Legal combinations and their
internal paths:

| mode | method | order | off-diag | vz_impl | internal method | cost | notes |
|---|---|---|---|---|---|---|---|
| flux | exp | 1 | — | — | `_estimate_flux` | M×1 mesolve | narrow Gaussian flux, single-sided FD |
| omega | exp | 1 | — | math | `_omega_vz_math` | M×2 mesolve | σ_z Gaussian impulse, bilateral FD |
| omega | exp | 1 | — | hardware | `_omega_vz_hardware` | M×2 mesolve | phase-shifted sub-pulses |
| omega | sim | 1 | — | — | `_heisenberg_kernels` | 1 sesolve | exact commutator formula |
| flux | exp | ≥2 | no | — | `_extract_kn_flux` | M×5 mesolve | 5-pt FD stencil |
| omega | exp | ≥2 | no | math | `_extract_kn_omega` | M×5 (or ×15 Richardson) mesolve | 5-pt FD + optional Richardson |
| omega | exp | ≥2 | no | hardware | `_extract_kn_omega` (hw branch) | M×5 mesolve | phase-shifted sub-pulses, 5-pt FD |
| omega | exp | ≥2 | yes | math | `_extract_kn_offdiag_exp` | O(M²)~O(M³) mesolve | mixed-partial FD, simultaneous σ_z kicks |
| omega | exp | ≥2 | yes | hardware | `_extract_kn_offdiag_exp` (hw branch) | O(M²)~O(M³) mesolve | simultaneous multi-time phase kicks (`with_phase_kicks`) |
| omega | sim | ≥2 | no | — | `_heisenberg_kernels` | 1 sesolve | nested commutators, any order |
| omega | sim | ≥2 | yes | — | `_heisenberg_kernels_offdiag` | 1 sesolve + O(Mⁿ) matrix ops | full n-D tensors, ≤order 3 |
| flux | sim | any | any | — | **illegal** | — | raises `ValueError` |

### Perturbation modes explained

**flux (magnetic)**: a narrow Gaussian flux signal (`FluxSignal(type=3)`) is
injected through `qubit.qubit_under_mag()` into the qubit dispersion model.
Kernel units: `1/(Φ₀·ns)`.  Probe width controlled by `stim_width`.

**omega / math VZ (σ_z Gaussian impulse)**: adds a
`(φ_z/2)·σ_z·gaussian(t−tⱼ)` term to the Hamiltonian.  Physically
equivalent to an instantaneous frequency shift, but finite σ_t convolves the
off-diagonal structure (k₃ ~15% low at σ_t=2·dt).  Richardson extrapolation
(`richardson=True`) samples at multiple σ_t and extrapolates to σ_t→0 to
remove this bias.  Kernel units: `rad⁻¹` (dimensionless).

**omega / hardware VZ (phase-shifted sub-pulses)**: does NOT add anything to
H — instead rebuilds the entire pulse sequence with all sub-pulses after the
probe time tⱼ phase-rotated by ±φ_z.  This is the closest model to actual
AWG phase-register operations.  Works on `CompositePulse`; single `Pulse`
falls back to math.  Phase kicks are instantaneous (no σ_t), so Richardson
extrapolation is not applicable.

### First-order kernel (order=1)

Fast path: only 1–2 `mesolve` calls per probe time (sim needs just 1 total
`sesolve`).  Suitable for most linear Wiener deconvolution use cases.

```python
# First-order flux (backward-compatible)
est = KernelEstimator(mode='flux', method='exp')
t, k1 = est.estimate(pulse, qubit)

# First-order hardware VZ (experiment-realistic)
est = KernelEstimator(mode='omega', method='exp', virtual_z_impl='hardware')
t, k1 = est.estimate(pulse, qubit)

# First-order sim (exact theoretical kernel, extremely fast)
est = KernelEstimator(mode='omega', method='sim')
t, k1 = est.estimate(pulse)  # qubit not needed
```

### Higher-order diagonal kernels (order≥2, off_diagonal=False)

5-point centred-difference stencil at five amplitudes
$\{-2h, -h, 0, +h, +2h\}$:

- $k_1 = (f_{-2} - 8f_{-1} + 8f_{+1} - f_{+2}) / 12h$
- $k_2 = (-f_{-2} + 16f_{-1} - 30f_0 + 16f_{+1} - f_{+2}) / 12h^2$
- $k_3 = (-f_{-2} + 2f_{-1} - 2f_{+1} + f_{+2}) / 2h^3$

5 mesolve calls per probe time; 15 with Richardson (3 σ_t values).  Output
is 1-D diagonal kernels $k_n(t, t, …, t)$.  Wiener / Hammerstein deconvolution
only consumes diagonal kernels.

```python
# Third-order exp/math diagonal + Richardson
est = KernelEstimator(mode='omega', method='exp', order=3,
                      richardson=True)
res = est.estimate_full(pulse, qubit)
k1, k2, k3 = res.kernels  # all 1-D

# Third-order exp/hardware diagonal
est = KernelEstimator(mode='omega', method='exp', order=3,
                      virtual_z_impl='hardware')
res = est.estimate_full(pulse, qubit)
```

### Higher-order off-diagonal kernels (extract_off_diagonal=True)

Returns full n-D Volterra tensors: $k_2$ as `(M,M)`, $k_3$ as `(M,M,M)`.
**Only LM full-density-matrix inversion can consume off-diagonal kernels**;
Wiener / Hammerstein paths accept diagonal kernels only.

**Sim path** (Heisenberg): exact solution, one `sesolve`, zero perturbation.
For $k_3$ this evaluates $M^3$ commutators; slower for large M.

**Exp path** (finite-difference): mixed-partial FD generalises the 5-point
stencil — $k_2(i,j)$ uses a 4-corner mixed difference
$[p(+,+)-p(+,-)-p(-,+)+p(-,-)]/4h^2$; $k_3(i,j,l)$ handles three cases
per index degeneracy (all-equal, two-equal, all-distinct).  A $p_e$ cache
avoids redundant mesolve calls.

**Hardware VZ** (`virtual_z_impl='hardware'`): applies multiple simultaneous
phase kicks via `pulse.with_phase_kicks([(t_i,+h), (t_j,-h)])` — one pulse
reconstruction with multiple phase jumps.  Single `Pulse` falls back to math.

```python
# sim off-diagonal (exact; recommended for G₃ calibration)
est = KernelEstimator(mode='omega', method='sim', order=3,
                      extract_off_diagonal=True)
res = est.estimate_full(pulse, qubit)
k1 = res.kernels[0]   # (M,)
k3 = res.kernels[2]   # (M, M, M) — full third-order tensor

# exp/hardware off-diagonal (experiment-realistic)
est = KernelEstimator(mode='omega', method='exp', order=2,
                      extract_off_diagonal=True,
                      virtual_z_impl='hardware')
res = est.estimate_full(pulse, qubit)
k1 = res.kernels[0]   # (M,)
k2 = res.kernels[1]   # (M, M)
```

### Downstream consumption constraints

| Kernel type | Wiener | Hammerstein | Hammerstein-Volterra | LM |
|---|---|---|---|---|
| 1-D diagonal ($k_1$) | ✓ | ✓ | ✓ | ✓ |
| 1-D diagonal ($k_1,k_2,k_3$) | ✗ | ✗ | ✓ | ✓ |
| n-D off-diagonal | ✗ | ✗ | ✗ | ✓ |

### KernelResult

Kernel-estimation result container.  `kernels[n-1]` is the n-th order kernel:
`kernels[0]`=$k_1$ (1-D); if extracted, `kernels[1]`=$k_2$ (diagonal 1-D or
off-diagonal 2-D), `kernels[2]`=$k_3$.

**Fields**

| Field | Type | Meaning |
|---|---|---|
| `t_samples` | `np.ndarray` | Kernel sample time axis (ns) |
| `kernels` | `list[np.ndarray]` | Per-order kernel list; `kernels[0]`=k₁ |
| `mode` | `str` | Perturbation channel |
| `method` | `str` | Estimation method |
| `order` | `int` | Highest kernel order |
| `stim_amplitude` | `float` | Stimulus amplitude |
| `units` | `str` | Kernel units (flux: `1/(Φ₀·ns)`, omega: `rad⁻¹`) |
| `off_diagonal` | `bool` | `True` if higher-order kernels are n-D tensors |

**Property**: `k1` (=`kernels[0]`, linear-kernel shortcut).

**Methods**: `save(path)` serialise to `.npz` / `load(path)` (classmethod)
deserialise.

## RamseyReconstruction: Ramsey reconstruction

Inverts the magnetic field $B(\tau)$ from the Ramsey phase.

**Construction**

`RamseyReconstruction(qubit, method="unwrap", k_span=3)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit (provides frequency sensitivity $\kappa$) | — |
| `method` | str | Phase extraction path | `"unwrap"` (or `"iq"`) |
| `k_span` | int | Unwrap branch-search half-window | `3` |

**Methods**

- `reconstruct(measurement) -> np.ndarray`:

  - `method="iq"`: dual-channel $\arctan2$ phase extraction (needs
    `data["p_e_I"]`/`p_e_Q`).
  - `method="unwrap"`: single-channel $\arccos$ + $k$-span branch unwrapping
    (needs `data["p_e"]`).

  Both paths take the gradient of the phase and divide by the frequency
  sensitivity $\kappa=\mathrm{d}\omega/\mathrm{d}\Phi$ to get
  $B=\dot\varphi/\kappa$.

**Output**

Returns `np.ndarray`, the magnetic field $B(\tau)$ array.

## EchoReconstruction: differential echo reconstruction

Inverts the field analytically from the differential echo $p_e$.

**Construction**

`EchoReconstruction(qubit, t_int, k)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit | — |
| `t_int` | float | Interaction time per echo block (ns) | — |
| `k` | int | Number of $\pi$ pulse pairs | — |

**Methods**

- `reconstruct(measurement) -> np.ndarray`: takes `data["p_e"]`,
  $\varphi=\arcsin(2p_e-1)$, $B=-\varphi/(2k\kappa t_\mathrm{int})$.

**Output**

Returns `np.ndarray`, the magnetic field $B$ array.

```{note}
A uniform `clip` to $[-1,1]$ precedes `arcsin`/`arccos`: `mesolve` integration
noise can push $p_e$ marginally outside $[0,1]$, which without clipping produces
silent NaNs. The Ramsey/Echo/Transient paths all have the same guard.
```

## TransientReconstruction: transient-field reconstruction

The core deconvolver of the waveform-reconstruction mainline: it
deconvolves the $\Delta p$ and control kernel of the transient experiment back
into $\Phi(t)$.

**Construction**

`TransientReconstruction(method="wiener", lambda_reg=1e-3, qubit=None, control_pulse=None, basis_type="bspline", n_basis=8, max_iter=50, tol=1e-6, mu_init=0.01, use_adjoint=False, max_volterra_iter=5, volterra_tol=1e-4)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `method` | str | Inversion algorithm | `"wiener"` (or `"hammerstein"`/`"hammerstein_volterra"`/`"lm"`) |
| `lambda_reg` | float | Wiener regularization parameter | `1e-3` |
| `qubit` | `TransmonQubit` | Qubit (needed by Hammerstein/LM) | `None` |
| `control_pulse` | `CompositePulse` | Control pulse (needed by LM) | `None` |
| `basis_type` | str | LM basis function type | `"bspline"` |
| `n_basis` | int | LM number of basis functions | `8` |
| `max_iter` | int | LM max iterations | `50` |
| `tol` | float | LM convergence tolerance | `1e-6` |
| `mu_init` | float | LM initial damping | `0.01` |
| `use_adjoint` | bool | LM use adjoint-method Jacobian | `False` |
| `max_volterra_iter` | int | Hammerstein-Volterra max iterations | `5` |
| `volterra_tol` | float | Volterra convergence tolerance | `1e-4` |

**Methods**

- `reconstruct(measurement, kernel=None, **kwargs) -> FluxSignal`:

  - `"wiener"`: linear Wiener deconvolution (FFT domain).
  - `"hammerstein"`: Hammerstein-Wiener nonlinear block model (needs `qubit`).
  - `"hammerstein_volterra"`: iterative inversion using higher-order kernels
    (needs `order>=2` kernel).
  - `"lm"`: Levenberg-Marquardt full-density-matrix inversion (needs
    `control_pulse`).

  `kernel` can be an `np.ndarray` (linear kernel) or a `KernelResult` (with
  higher-order kernels). Off-diagonal n-D kernels are consumable only by
  `method="lm"`; the Wiener/Hammerstein paths accept diagonal (1-D) kernels only,
  else they raise `ValueError`.

**Output**

Returns `FluxSignal`, the reconstructed flux waveform $\Phi(t)$.

## CryoscopeReconstruction: Cryoscope reconstruction

Reconstructs the waveform from phase-vs-truncation data; the core physics is
$\mathrm{d}\varphi/\mathrm{d}t=2\pi\Delta f(h(t))$.

**Construction**

`CryoscopeReconstruction(tau=None, inversion="calibration", calibration=None, qubit=None, use_sg_filter=True, sg_window=11, sg_poly=3)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `tau` | float | Calibration square-pulse length (ns) | — |
| `inversion` | str | Inversion path | `"calibration"` (or `"response"`) |
| `calibration` | `CalibrationTable` | $\varphi(h)$ table (required for `inversion="calibration"`) | `None` |
| `qubit` | `TransmonQubit` | Qubit (required for `inversion="response"`) | `None` |
| `use_sg_filter` | bool | Use Savitzky-Golay pre-smoothing | `True` |
| `sg_window` | int | SG window width | `11` |
| `sg_poly` | int | SG polynomial order | `3` |

**Methods**

- `reconstruct(measurement, ...) -> FluxSignal`: reads `data["varphi"]`,
  `axes["trunc"]`.
  - `inversion="calibration"`: invert via the $\varphi(h)$ lookup table produced
    by `CryoscopeCalibration`.
  - `inversion="response"`: invert via the analytical Transmon frequency-flux
    dispersion.

**Output**

Returns `FluxSignal`, the reconstructed flux waveform.

## DelayRamseyReconstruction: delay Ramsey reconstruction

Inverts the delay-Ramsey phase $\varphi(t_d)$ into the tail flux waveform (the
predistortion mainline).

**Construction**

`DelayRamseyReconstruction(inversion="response", qubit=None, calibration=None, tau_R=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `inversion` | str | Inversion path | `"response"` (or `"calibration"`) |
| `qubit` | `TransmonQubit` | Qubit (required for `inversion="response"`) | `None` |
| `calibration` | `CalibrationTable` | $\varphi(z)$ table (required for `inversion="calibration"`) | `None` |
| `tau_R` | float | Ramsey free-evolution time (ns) | — |

**Methods**

- `reconstruct(measurement, ...) -> FluxSignal`: reads `data["varphi"]`,
  `axes["t_d"]`.

**Output**

Returns `FluxSignal`, the tail flux waveform $\Phi_\mathrm{tail}$.

## PiPulseCompReconstruction: pi-pulse compensation reconstruction

A direct reconstruction: once the 2D scan gives the optimal compensation
height $z^*(\tau)$, the tail is $\Phi_\mathrm{tail}=-z^*(\tau)$.

**Construction**

`PiPulseCompReconstruction()`

(No extra fields needed — only `measurement` is required for inversion.)

**Methods**

- `reconstruct(measurement) -> FluxSignal`: reads `data["z_star"]`, `axes["tau"]`.

**Output**

Returns `FluxSignal`, the tail flux waveform.

```{note}
Because the $\pi$ pulse has a finite width $T$, $z^*(\tau)$ measures the moving
average of the tail over $[\tau, \tau+T]$, not the instantaneous value. For
exponential tails this only introduces a constant attenuation factor; the decay
time constant is preserved.
```

## Calibration classes (CryoscopeCalibration / DelayRamseyCalibration)

Two classes exported in this namespace but conceptually belonging to
calibration (subclasses of {py:class}`~sqc.calibration.Calibration`), which
build the lookup tables for their reconstructors' `"calibration"` inversion
path:

### CryoscopeCalibration

**Construction**

`CryoscopeCalibration(qubit, t_rabi=None, tau=20.0)`

**Methods**

- `calibrate() -> CalibrationTable`: scan square-pulse heights $h$, IQ readout +
  model-guided unwrap to build a $\varphi(h)$ lookup, `kind="phi_h"`, for
  `CryoscopeReconstruction`(inversion="calibration").

### DelayRamseyCalibration

**Construction**

`DelayRamseyCalibration(qubit, t_rabi=None, tau_R=None)`

**Methods**

- `calibrate() -> CalibrationTable`: scan known flux heights $z$, fit the linear
  slope $\varphi_\mathrm{cal}(z)=\tau_R\kappa z$, `kind="phi_z"`, for
  `DelayRamseyReconstruction`(inversion="calibration").

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

- This layer corresponds to the offline analysis of experimental data:
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
- To add a custom reconstruction algorithm, subclass the `Reconstruction`
  abstract base and implement `reconstruct()`, upholding the two invariants
  (pure function, output a `FluxSignal`), and inject `qubit` at construction if
  forward simulation is needed. Full extension guide in {doc}`../extending`.
