# Control

## What this layer provides

The `control` layer defines the **time-domain signal types** shared across the
stack and provides all the components needed to build experimental control
sequences: flux trajectories (`FluxSignal`), microwave pulse envelopes
(`Pulse`), multi-pulse experimental sequences (`CompositePulse`), and built-in
qubit-gate functions.

This layer is where signal types are defined. `Waveform` is the lowest-level
time-domain container; `FluxSignal` is its subclass (samples in units of
$\Phi_0$). The `ControlLine` / `DistortionModel` / `TransferMatrix` of the
{doc}`hardware` layer all take these two types as input and output. They are
declared there via `import` rather than redefined here, because `hardware`
depends on `control`, not the other way around.

This layer's responsibility ends at pulse design: assembling parameterized
signals into a QuTiP list-format Hamiltonian, ready to hand to the
{doc}`simulation` layer's `HamiltonianBuilder` / `MesolveRunner` for time
evolution. Transfer-function distortion between the AWG and the chip is not
handled here; that is the job of the {doc}`hardware` layer.

One convention runs through this layer: all time axes are in ns and derive from
`sqc.config.CONFIG.awg.dt` (built with `np.arange`, not `np.linspace`). Internal
time axes such as free-evolution gaps are generated from the global `dt`, so
pulse sequences project onto the unified time grid used by the simulation.

## Class overview

| Class / function | Role | Notes |
|---|---|---|
| `Waveform` | Signal primitive | Semantically neutral time-domain container $(t\_list, samples)$; the lowest-level type in this layer |
| `CompositeWaveform` | Composite signal | Multiple `Waveform`s concatenated in time order |
| `FluxSignal` | Flux signal | `Waveform` subclass, samples in $\Phi_0$; supports 9 typed constructions |
| `CompositeSignal` | Composite flux | Multiple `FluxSignal`s concatenated in order |
| `PulseBase` | Abstract base | Common contract for any control pulse; the layer's extension point |
| `Pulse` | Single pulse | One control pulse (lab/rotating frame), taking a `FluxSignal` as its Rabi envelope |
| `CompositePulse` | Composite pulse | Multiple `Pulse`s concatenated into a sequence |
| `PulseSequence` | Sequence container | A dataclass holding an ordered set of pulses |
| `create_*` functions | Sequence builders | Generate standard pulse sequences: Ramsey / echo / Cryoscope, etc. |
| `make_drag_envelope` | Envelope generator | Generates the I/Q envelope pair of a DRAG pulse |
| `ideal_*` / `simulate_*` | Two-qubit gates | Ideal matrices and simulation of iSWAP / CZ |

## Waveform: the time-domain signal primitive

The lowest-level signal type in this layer, semantically neutral: just the
container $(t\_list, samples)$, carrying no physical-unit meaning. It is the base
class of `FluxSignal`, the carrier of a `Pulse` envelope, and the unified
input/output type of the {doc}`hardware` layer's `ControlLine` /
`DistortionModel`.

**Construction**

`Waveform(t_list, samples, metadata=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `t_list` | `np.ndarray` | Time points (ns) | — |
| `samples` | `np.ndarray` | Signal value at each time | — |
| `metadata` | dict | Additional metadata | `None` |

Construction raises if the two arrays have mismatched shapes.

**Read-only properties**: `duration`, `n_points`.

**Methods**

- `value_at(t) -> float`: sample-and-hold, 0 out of range.
- `samples_on(t_global) -> np.ndarray`: linear interpolation onto a global time
  axis, 0 outside the window.
- `truncate(t_start, t_end) -> Waveform`: returns a **new** waveform, zeroed
  outside the window.
- `copy() -> Waveform`: deep copy.
- `plot(ax=None)`: quick visualisation.

`CompositeWaveform` is its subclass, holding a `components` list; the class method
`CompositeWaveform.from_components(components)` concatenates several waveforms in
time order.

## FluxSignal: flux signal

A subclass of `Waveform` with samples in units of $\Phi_0$, the physical
representation of a flux trajectory. It is compatible with the legacy
`src.signal.Signal` and supports parameterized construction by `type` (the alias
`Signal = FluxSignal` is retained).

**Construction**

`FluxSignal(type=0, t_list=None, trigger=0.0, **kwargs)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `type` | int | Signal type (0–8) | `0` |
| `t_list` | `np.ndarray` | Time axis (ns) | — |
| `trigger` | float | Trigger offset (ns) | `0.0` |
| `signal` | `np.ndarray` | Signal samples (alias for `samples`, read/write) | — |
| `params` | dict | Construction parameter record | — |

Values of `type`:

| type | Waveform | type | Waveform |
|---|---|---|---|
| 0 | zero | 5 | double peak |
| 1 | constant | 6 | basis expansion (B-spline/Fourier/Legendre) |
| 2 | sinusoid | 7 | complex wave packet |
| 3 | Gaussian | 8 | user-defined raw samples |
| 4 | asymmetric impulse (double exponential) | | |

**Parameters** (used per `type`): `amplitude`, `frequency`, `phase`, `center`,
`width`, `offset`, `rise`/`fall` (type 4), `noise_level`/`seed`; type 6 also has
`n_basis`, `b` (basis coefficients), `basis_type`. type 8 must pass
`signal=<array>`.

**Methods**

- `value_at(t) -> float` / `samples_on(t_global) -> np.ndarray`: with `trigger`
  offset.
- `truncate(t_start, t_end)`: zeroes **in place** (unlike the base class which
  returns a new object).
- `update_signal(**kwargs)`: change parameters and regenerate in place.
- `copy() -> FluxSignal`: deep copy preserving current samples.
- `plot()`: quick visualisation.

`CompositeSignal(signals)` concatenates several `FluxSignal`s in order (a
`CompositeWaveform` subclass).

**Output**

Itself is the signal object; `samples_on()` returns `np.ndarray` projected onto
the global axis.

### make_drag_envelope: DRAG envelope generator

`make_drag_envelope(t_list, amplitude, sigma, beta, phi=0.0) -> (FluxSignal,
FluxSignal)`: generates the I/Q envelope pair $(\Omega_I, \Omega_Q)$ of a DRAG
pulse. `beta` is $-1/\alpha$ (the inverse anharmonicity); `phi=0` gives an X
rotation and `phi=\pi/2` a Y rotation. It returns two `type=8` `FluxSignal`s,
which can be fed directly to a `Pulse`'s `Omega` / `Omega_Q`.

## PulseBase: the pulse abstract base

The common contract for any control pulse, and the layer's extension point.
Subclasses must provide three instance attributes: `hamiltonian` (a QuTiP
list-format time-dependent Hamiltonian), `t_list` (time axis), and `frame`
(reference frame, 0=lab, 1=rotating). To add a custom pulse type, subclass
`PulseBase` and populate these three in `__init__`; see {doc}`../extending`.

## Pulse: a single control pulse

One control pulse. It takes a `FluxSignal` (or scalar) as the Rabi envelope
`Omega` and builds a list-format Hamiltonian, in the lab frame (`frame=0`) or the
rotating frame (`frame=1`, with optional RWA).

**Construction**

`Pulse(frame=0, omega_d=0.0, phase=0.0, Omega=None, Omega_Q=None, is_rwa=True, qubit=None, trigger=0.0)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `frame` | int | Reference frame (0=lab, 1=rotating) | `0` |
| `omega_d` | float | Drive frequency (rad·GHz) | `0.0` |
| `phase` | float | Rotation-axis phase | `0.0` |
| `Omega` | `FluxSignal` or float | Rabi envelope (I component) | `None` |
| `Omega_Q` | `FluxSignal` or float | DRAG Q component (`None` = single-quadrature) | `None` |
| `is_rwa` | bool | Apply rotating-wave approximation | `True` |
| `qubit` | `TransmonQubit` | Qubit (for angle calculation) | `None` |
| `trigger` | float | Absolute start time within a sequence (ns) | `0.0` |

After construction, `hamiltonian` / `t_list` are ready.

**Methods**

- `get_Rabi_frequency(t) -> float` / `get_hamiltonian_at(t) -> Qobj`: the Rabi
  value / Hamiltonian at a single time.
- `get_hamiltonian() -> (H_list, t_list)`: the list-format Hamiltonian on the
  local time axis.
- `hamiltonian_on(t_global) -> list`: projects the local Hamiltonian onto a
  global time axis (including the `trigger` offset, 0 outside the window); the
  canonical path for splicing onto the unified simulation grid.
- `get_angle(qubit=None) -> (theta, phi)`: numerically evolves to find the
  pulse's rotation angle and phase on the qubit; `get_angle_simple() -> float`:
  the simple Rabi angle $\int\Omega\,\mathrm{d}t$.
- `with_phase_shift(phi_z) -> Pulse`: returns a **new** `Pulse` with the phase
  incremented by `phi_z` (a virtual Z gate).

**Output**

`get_hamiltonian()` returns `(H_list, t_list)`; `hamiltonian_on()` returns the
global-axis list-format Hamiltonian (ready for `mesolve`); `with_phase_shift()`
returns a new `Pulse`.

```{note}
`Pulse.get_kernel()` is **deprecated** and forwards to the {doc}`reconstruction`
layer's `KernelEstimator`. New code should use `KernelEstimator` directly.
```

## CompositePulse: a composite pulse sequence

Multiple `Pulse`s concatenated into one sequence.

**Construction**

`CompositePulse(pulses)`

**Fields**

| Field | Type | Meaning |
|---|---|---|
| `pulses` | `list[Pulse]` | Sub-pulses (in time order) |
| `frame` | int | Reference frame (taken from first sub-pulse) |
| `omega_d` | float | Drive frequency (taken from first sub-pulse) |

**Methods**

- `get_Omega(t) -> float` / `get_hamiltonian_at(t) -> Qobj`: the whole
  sequence's Rabi value / Hamiltonian at time t.
- `hamiltonian_on(t_global) -> list`: collects each sub-pulse's contribution
  (each carrying its own `trigger`) on the global axis.
- `with_phase_shift(phi_z, from_time=None) -> CompositePulse`: increments the
  phase of sub-pulses with `trigger >= from_time`, returning a **new** sequence.
- `plot()`: visualise sub-pulse envelopes.

`get_kernel()` is likewise deprecated (forwards to `KernelEstimator`).

**Output**

`hamiltonian_on()` returns the stitched global list-format Hamiltonian;
`with_phase_shift()` returns a new `CompositePulse`.

## PulseSequence: sequence container

A dataclass holding an ordered set of pulses.

**Fields**

| Field | Type | Meaning |
|---|---|---|
| `pulses` | `list` | Pulse list |
| `metadata` | dict | Additional metadata |

## Sequence builders

A set of functions for building standard experimental sequences, returning a
`CompositePulse` (`create_pulse` returns a `QobjEvo`). All free-evolution gaps
are generated from the global `dt`; each sub-pulse carries an absolute
`trigger`, so it projects onto the simulation time axis.

- `create_pulse(qubit, frame, type, t_list, omega_d, phase, angle=None,
  trigger=0.0, Omega_Q=None, **kwargs) -> QobjEvo`: a single pulse; when a target
  `angle` is given, the amplitude is auto-calibrated (for DRAG, the Q envelope is
  scaled by the same factor).
- `create_ramsey_pulse(t_rabi, tau, omega_d=0.0, phase1=π/2, phase2=0.0,
  qubit=None, trigger=0.0, rotation_angle=None, rabi_rate=None,
  envelope="square", envelope_sigma=None) -> CompositePulse`: configurable
  Ramsey sequence. The default remains a square-pulse
  $\pi/2-\tau-\pi/2$ sequence. `rotation_angle` normalizes the envelope area
  to a target angle; `rabi_rate` fixes its peak and lets duration determine the
  achieved angle. The modes are mutually exclusive. `envelope` accepts
  `"square"`, `"gaussian"`, or an array matching `t_rabi`.

```{note}
Sweeping `rotation_angle` at fixed duration changes the Rabi amplitude. To test
the short-duration, small-angle regime at fixed drive strength, keep
`rabi_rate` fixed and vary `t_rabi`. Re-estimate the response kernel for every
control configuration.
```
- `create_echo_pulse(t_rabi, tau, omega_d=0.0, phase1=0, phase2=0, phase3=0,
  trigger=0.0) -> CompositePulse`: spin echo:
  $\pi/2 - \tau - \pi - \tau - \pi/2$.
- `create_diff_echo_pulse(t_rabi, tau, t_int, t_rep, k, omega_d, ...)
  -> CompositePulse`: differential echo (k repetitions).
- `create_cpmg_pulse(t_rabi, tau, n, omega_d, phase1=0.0, phase2=π/2,
  phase3=π/2, trigger=0.0) -> CompositePulse`: CPMG sequence:
  $\pi/2 - [\tau/2 - \pi - \tau - \pi - \tau/2] - \pi/2$, with `n` π pulses.
- `create_cryoscope_pulse(t_rabi, tau, omega_d, phase1=π/2, phase2=0.0,
  trigger=0.0) -> CompositePulse`: Cryoscope: $\pi/2 - \tau - \pi/2$ (used in
  the waveform-reconstruction workflow).
- `create_pi_pulse_compensation_pulse(t_rabi, T_pi, omega_d, phase=0.0,
  qubit=None, trigger=0.0) -> CompositePulse`: the microwave drive for the
  π-pulse compensation protocol (the compensation flux is applied separately via
  `qubit_in_mag`, not part of this sequence).

## Two-qubit gates

Ideal matrices and simulation functions for two-qubit gates (in the
$|00\rangle,|01\rangle,|10\rangle,|11\rangle$ basis):

- `ideal_iSWAP() -> Qobj` / `ideal_CZ() -> Qobj`: the ideal $4\times4$ unitary
  matrices.
- `simulate_iSWAP(qubit1, qubit2, g) -> (U_eff, leakage, F)`: iSWAP simulation
  under exchange coupling `g` (GHz), returning the effective unitary, leakage,
  and average gate fidelity.
- `simulate_CZ(qubit1, qubit2, g) -> (U_eff_corr, leakage, F)`: CZ via flux
  pulsing, with phase correction.
- `simulate_cz_from_flux(t, phi_chip, qubit1, qubit2, g, ...) -> CZResult`
  (v2.6): waveform-driven CZ --- accepts an **explicit on-chip flux waveform**
  ``phi_chip(t)`` (in Φ₀), converts to detuning via the transmon dispersion
  $\omega(\Phi) = \sqrt{8 E_C E_J|\cos(\pi\Phi)|} - E_C$, builds the full
  time-dependent two-qubit Hamiltonian, and computes all gate metrics from the
  evolution operator. Returns a ``CZResult`` with ``conditional_phase``,
  ``phase_error``, ``leakage``, ``infidelity``, and key populations
  (``pop_20``, ``pop_02``). This enables **end-to-end predistortion → CZ gate**
  validation by feeding control-line-distorted (and corrected) waveforms
  directly into the gate simulation.

```{note}
The full device model of two-qubit gates (tunable coupler, multimode cavity) is
in `CoupledSystem` on the {doc}`devices` page.
```

## Minimal example

```python
import numpy as np
from sqc.config import CONFIG
from sqc.devices import TransmonQubit
from sqc.control import FluxSignal, create_ramsey_pulse

# 1) Build a flux signal directly: a Gaussian flux pulse
t = CONFIG.pulse.t_signal.copy()
flux = FluxSignal(type=3, t_list=t, amplitude=0.02, center=t[-1] / 2, width=20.0)
print(flux.samples.max())              # peak flux (Φ_0)

# 2) Build a Ramsey sequence (π/2 - τ - π/2), returning a CompositePulse
qubit = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10_000, T2=8_000)
t_rabi = CONFIG.pulse.t_rabi.copy()
ramsey = create_ramsey_pulse(t_rabi, tau=50.0, omega_d=qubit.frequency, qubit=qubit)

# 3) Project onto the global time axis → hand to the simulation layer to build H
t_global = CONFIG.pulse.t_global.copy()
H_ctrl = ramsey.hamiltonian_on(t_global)   # QuTiP list format, fed to mesolve
```

## Physical role / extension

- `FluxSignal` corresponds to the real **flux trajectory** applied to a Z line
  (the ideal, undistorted one); only after passing through the {doc}`hardware`
  layer's `ControlLine` / `DistortionModel` does it become the real on-chip
  waveform.
- `Pulse` / `CompositePulse` correspond to microwave XY drive: Rabi envelopes,
  DRAG correction, multi-pulse experimental sequences.
- The sequence builders correspond to the standard pulse sequences commonly used
  in experiments (Ramsey, echo, Cryoscope, etc.).
- To add a custom pulse type, subclass `PulseBase` and provide `hamiltonian`
  / `t_list` / `frame`. To add a custom waveform type, subclass `Waveform`
  (with physical units, as `FluxSignal` does). See the full extension guide in
  {doc}`../extending`.
