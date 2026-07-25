# Simulation

## What this layer provides

The `simulation` layer is the stack's **solver layer** — sitting between
{doc}`control` (which produces pulses and Hamiltonian material) and
{doc}`experiments` (which orchestrates full protocols), it does exactly three
things: **assemble devices and signals into a time-dependent Hamiltonian, drive
the QuTiP solver through the time evolution, and pack the output into a uniform
result container**. It holds no protocol-orchestration logic (that lives in the
experiments layer) and no reconstruction/inversion logic (that lives in the
{doc}`reconstruction` layer).

One design constraint runs through this layer: **Hamiltonian construction is a
pure function** — `HamiltonianBuilder` only reads `qubit` / `flux_signal` /
`pulse` and returns a QuTiP list-format `H_list` plus a global time axis; it
**never mutates any input object**. The same devices and signals can therefore
be reused across frames and solvers with no side effects.

All time axes are in ns and derive uniformly from `sqc.config.CONFIG.awg.dt`
(see {doc}`control`). This layer follows the same grid, so pulses, flux, and
free evolution splice together seamlessly.

## Class overview

| Class / function | Role | Notes |
|---|---|---|
| `HamiltonianBuilder` | Hamiltonian builder | Pure function: device + flux + pulse → QuTiP list-format `H_list` and global time axis |
| `RunnerBase` | Abstract base | Common contract for all solvers; this layer's extension point |
| `MesolveRunner` | Main solver | Wraps `qutip.mesolve`, takes `H_list` and runs master-equation evolution |
| `SlidingMeasurementRunner` | Sliding-measurement solver | Slides a control pulse across a flux signal, measuring excited-state population $p_e$ per delay |
| `ExperimentResult` | Result container | Serialisable `data/axes/metadata/config` with `save`/`load` |
| `MeasurementTrace` | Single trace | frozen dataclass, one 1-D sweep ($axis, p_e$) |
| `extract_expectation` | Extractor | Pull the expectation-value array from a QuTiP result |
| `extract_population` | Extractor | Pull a Fock-level population from `result.states` |
| `generate_1f_noise` | Noise generator | Generate a $1/f$ noise time series via FFT |

## HamiltonianBuilder — time-dependent Hamiltonian builder

The layer's entry point. It assembles a device, flux signal, and control pulse
into a QuTiP list-format Hamiltonian that `mesolve` can consume. **Pure function,
no side effects**: every method is a `@staticmethod` and mutates none of its
inputs.

**Interface**: `HamiltonianBuilder.build(qubit, flux_signal=None, pulse=None,
frame="rotating", omega_d=None) -> (H_list, t_global)`.

- `qubit` — a `QubitSpec` or `TransmonQubit` (the latter auto-calls `.spec()` to
  extract its inner `QubitSpec`).
- `flux_signal` — flux signal (needs `.t_list` and `.samples`/`.signal`).
  Passing `None` falls back to a **static** Hamiltonian with time axis
  `CONFIG.pulse.make_time(0, 100)`.
- `pulse` — control pulse; when non-`None`, each term of `pulse.hamiltonian` is
  appended to `H_list`.
- `frame` — `"lab"` (lab frame) or `"rotating"` (rotating frame). In the
  rotating frame the frequency coefficients are shifted by `omega_d` (defaulting
  to the sweet-spot frequency `spec.frequency()`).

**Returned `H_list` structure** (aligned term-by-term with
`src/qubit.py:qubit_in_mag`):

$$H_0 = \tfrac{\alpha}{2}(n^2 - n), \qquad
H_1 = \big(n + \tfrac12 I\big)\cdot \omega_{01}(\Phi(t))$$

Here $H_0$ is a constant-coefficient term (anharmonicity), and $H_1$'s
coefficient array comes from substituting the flux at each instant:
$\omega_{01}(\Phi(t))=$ `spec.frequency(spec.flux_bias + s)`.

```{note}
Both frames are written as `n + 0.5*I`. That `0.5*I` is a **constant energy
shift** that does not affect the dynamics; it is kept solely to stay
term-by-term identical to the output of `src/qubit.py:qubit_in_mag` (the
regression baselines depend on this).
```

## RunnerBase — solver abstract base class

The common contract for all simulation solvers, and this layer's **extension
point**. It mandates a single abstract method:

- `run(*args, **kwargs)` — execute the time evolution and return an
  `ExperimentResult`.

To plug in a new solver (e.g. `sesolve`, Monte-Carlo trajectories, an external
numerical backend), subclass `RunnerBase`, implement `run()`, and have it return
an `ExperimentResult`. See {doc}`../extending`.

## MesolveRunner — master-equation solver

The layer's workhorse solver, a direct wrapper around `qutip.mesolve`.
`@dataclass` with a single field `options` (the options dict passed to
`mesolve`; `None` uses defaults).

**Interface**: `run(H_list, psi0, t_list, c_ops, e_ops, store_states=False)
-> ExperimentResult`.

- `H_list` — the list-format Hamiltonian from `HamiltonianBuilder.build`
  (internally wrapped with `QobjEvo(H_list, tlist=t_list, order=1)`).
- `psi0` / `c_ops` / `e_ops` — initial state, collapse operators, expectation
  operators (the latter two often taken directly from
  `qubit.get_collapse_operators()` and a custom projector).
- `store_states` — whether to keep the state at each step (for use by
  `extract_population`).

In the returned `ExperimentResult`: `data["expect"]` is the expectation-value
array, `data["states"]` is the state list (`None` unless `store_states` is set),
and `axes["t"]` is the time axis.

## SlidingMeasurementRunner — sliding-measurement solver

Slides a control pulse across a flux signal **delay by delay**, running one
evolution at each delay, measuring the final excited-state population $p_e$, and
producing a 1-D trace of $p_e$ versus delay. It internalises the logic of the
legacy `src/protocal.py:Protocal.sliding_measrement` / `single_measurement`
(the legacy spelling `sliding_measrement` is retained only in the `src_mirror`
compatibility layer). `@dataclass`, with the same `options` field.

**Interface**: `run(qubit, phi_signal, control_pulse, scan_list=None)
-> ExperimentResult`.

- `qubit` — qubit object (must support `qubit_under_mag`, i.e. the per-instant
  sequence of qubits under flux).
- `phi_signal` — flux signal (needs `.t_list` and `.value_at()`).
- `control_pulse` — control pulse (`CompositePulse`).
- `scan_list` — delay scan axis; when `None`, auto-generated from
  `CONFIG.awg.dt` (covering the full overlap of pulse and signal).

In the returned `ExperimentResult`: `data["p_e"]` is the excited-state
population per delay, `axes["scan"]` is the delay axis, and
`metadata["elapsed"]` records the wall time.

```{note}
`SlidingMeasurementRunner` runs `mesolve` with tightened tolerances
(`atol=1e-10, rtol=1e-8`) at each delay; with many samples this can be slow, and
`run()` prints the total elapsed time.
```

## Result containers and extractors

### ExperimentResult

A serialisable result container, `@dataclass` with four dict fields:

| Field | Meaning |
|---|---|
| `data` | Named data arrays (`expect`, `p_e`, `states`, etc.) |
| `axes` | Named coordinate axes (`t`, `scan`, `tau`, etc.) |
| `metadata` | Free-form metadata (solver name, elapsed time, protocol type, etc.) |
| `config` | Simulation config snapshot (EC, EJ, n_levels, etc.) |

**Methods**: `save(path)` — pickle to disk (protocol 5); `ExperimentResult.load(path)`
(classmethod) — read back, raising `TypeError` on a type mismatch.

### MeasurementTrace

`@dataclass(frozen=True)` representing a single 1-D sweep trace: fields `axis`
(sweep axis), `p_e` (excited-state population), `p_e_iq` (I/Q projections, may be
`None`), and `metadata`.

### Extractors

- `extract_expectation(result, e_ops_index=0) -> np.ndarray` — pull the
  `e_ops_index`-th expectation-value array from a QuTiP `mesolve`/`sesolve`
  result; raises `ValueError` if there is no `.expect`.
- `extract_population(result, level) -> np.ndarray` — project `result.states`
  step by step to get the population of Fock `level` (requires solving with
  `store_states=True`); raises `ValueError` if there is no `.states`.

## generate_1f_noise — 1/f noise generation

`generate_1f_noise(t_list, amplitude, f_min, f_max, seed=None) -> np.ndarray`
— generate a real-valued $1/f$ noise time series via FFT: fill a complex
Gaussian spectrum with a $1/\sqrt{f}$ spectral density over $[f_{\min},
f_{\max}]$, enforce Hermitian symmetry, then inverse-transform and take the real
part. Fixing `seed` makes it reproducible. Used to inject $1/f$-type dephasing
noise into flux/frequency.

## Minimal example

```python
import numpy as np
from qutip import basis, num
from sqc.config import CONFIG
from sqc.devices import TransmonQubit
from sqc.control import FluxSignal
from sqc.simulation import HamiltonianBuilder, MesolveRunner, extract_expectation

# 1) Device + flux signal (Gaussian flux pulse, peak 0.02 Φ_0)
qubit = TransmonQubit(EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10_000, T2=8_000)
t = CONFIG.pulse.t_signal.copy()
flux = FluxSignal(type=3, t_list=t, amplitude=0.02, center=t[-1]/2, width=20.0)

# 2) Build the time-dependent Hamiltonian (rotating frame) — pure, does not touch qubit / flux
H_list, t_global = HamiltonianBuilder.build(qubit, flux_signal=flux, frame="rotating")

# 3) Hand it to the main solver, read out ⟨n⟩
n_lev = qubit.n_levels
res = MesolveRunner().run(
    H_list, psi0=basis(n_lev, 1), t_list=t_global,
    c_ops=qubit.get_collapse_operators(), e_ops=[num(n_lev)],
)
n_expect = extract_expectation(res)      # ⟨n⟩(t), shape (len(t_global),)

res.save("run.pkl")                       # persist result (read back with ExperimentResult.load)
```

## Physical role / extension

- `HamiltonianBuilder` corresponds to the step of "writing a physical device
  plus the applied flux/microwave fields as a time-dependent Hamiltonian" — the
  modelling step of a simulation, involving no numerical integration.
- `MesolveRunner` corresponds to open-system **master-equation time evolution**
  (with $T_1$/$T_2$ dissipation); `SlidingMeasurementRunner` corresponds to
  scan-style experiments like the sliding measurement window in transient-field
  sensing.
- `ExperimentResult` corresponds to the raw output record of one
  experiment/simulation, the uniform data interface for the {doc}`experiments`
  and {doc}`reconstruction` layers.
- **To plug in a custom solver**, subclass `RunnerBase`, implement `run()`, and
  return an `ExperimentResult`. Full extension guide in {doc}`../extending`.

