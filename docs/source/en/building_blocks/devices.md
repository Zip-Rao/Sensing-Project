# Devices

## What this layer provides

The `devices` layer is the **lowest layer** of the full cQED stack. It models
only the physical quantum devices themselves — superconducting Transmon qubits,
readout/coupling resonators, and the chip topology and coupled systems that
assemble them. This layer describes **physics only** (energies, frequencies,
Hilbert space, Hamiltonians, dissipation operators); it does not touch any
control pulses, time evolution, or reconstruction logic — those responsibilities
live in the `control`, `simulation`, and other layers above.

One design constraint runs through this layer: **devices are "pure parameter
objects"** — they describe physics but store no experimental state. To represent
a perturbed device (e.g. a qubit at nonzero flux), use a method that returns a
**new object** (`QubitSpec.with_flux`) rather than mutating in place. `QubitSpec`
is therefore a frozen dataclass.

## Class overview

| Class | Role | Notes |
|---|---|---|
| `Device` | Abstract base | Common contract for all physical devices |
| `QubitSpec` | Parameter kernel | Immutable Transmon parameter description; input to `HamiltonianBuilder` |
| `TransmonQubit` | Runtime device | Tunable Transmon holding QuTiP operators and Hamiltonian; the device object the upper layers accept by default |
| `Resonator` | Resonator | Multimode readout/coupling resonator |
| `CoupledSystem` | Coupled system | Two qubits coupled to a multimode cavity via a tunable coupler, for two-qubit gates |
| `ChipTopology` | Chip topology | Container for multiple qubits + resonators, lifting single-device operators into the whole-chip Hilbert space |

## Physical model and units

Transmon physics (following Gao 2021 §II.B–C, Eq. 13–20):

$$H = 4 E_C n^2 - E_J(\Phi)\cos\varphi, \qquad
E_J(\Phi) = E_{J,0}\,\lvert\cos(\pi \Phi/\Phi_0)\rvert$$

$$\omega_{01} = \sqrt{8 E_J E_C} - E_C, \qquad \alpha = -E_C$$

**Unit conventions (important)**: `EC` and `EJ` are passed in **rad·GHz**
(already including the $2\pi$ factor), flux is in units of $\Phi_0$, time is in
ns, and $\hbar = 1$. So $E_C/h = 200\,\mathrm{MHz}$ corresponds to the argument
`EC = 2*np.pi*0.2`; `frequency()` likewise returns rad·GHz (i.e. the **angular
frequency** $\omega_{01}$), and dividing by $2\pi$ gives the physical frequency
$f_{01}$ in GHz. The stack-wide default device parameters
($E_C/h=200\,\mathrm{MHz}$, $E_J/h=15\,\mathrm{GHz}$, $T_1=10\,\mu\mathrm{s}$,
$T_2=8\,\mu\mathrm{s}$) are centralized in `sqc.config.CONFIG`.

## Device — the device abstract base

The common contract for all physical devices and the layer's **extension
point**. Subclasses must be pure parameter objects (describing physics only,
storing no experimental state) and implement the following abstract members:

- `name` (property) — human-readable device identifier.
- `hilbert_dim() -> int` — the full Hilbert-space dimension.
- `hamiltonian_static() -> Qobj` — the **time-independent** Hamiltonian with no
  drive and no external field.
- `collapse_operators() -> list[Qobj]` — Lindblad collapse operators for
  dissipation/dephasing.

To add a custom device type (e.g. fluxonium), subclass `Device` and implement
these four members; see {doc}`../extending`.

## QubitSpec — the immutable parameter record

A **pure parameter description** of a tunable Transmon, and the layer's
**parameter kernel**: a `@dataclass(frozen=True)` that holds no experimental
state, with all methods being pure functions (take parameters, return values, do
not mutate `self`). The `simulation` layer's `HamiltonianBuilder` takes a
`QubitSpec` directly as input to construct the Hamiltonian needed for simulation
(when passed a `TransmonQubit`, it calls `.spec()` to extract the internal
`QubitSpec`).

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `name` | str | Device identifier | — |
| `EC` | float | Charging energy (rad·GHz) | — |
| `EJ` | float | Josephson energy at zero flux (rad·GHz) | — |
| `T1` | float | Relaxation time (ns) | — |
| `T2` | float | Dephasing time (ns) | — |
| `flux_bias` | float | Static flux bias ($\Phi_0$) | `0.0` |
| `n_levels` | int | Hilbert-space truncation dimension | `3` |

**Methods**

- `EJ_at(flux=None) -> float` — the effective Josephson energy at a given flux,
  $E_J(\Phi)=E_{J,0}\lvert\cos(\pi\Phi/\Phi_0)\rvert$; uses `flux_bias` when
  `flux=None`.
- `frequency(flux=None) -> float` — the transition angular frequency
  $\omega_{01}(\Phi)=\sqrt{8E_J(\Phi)E_C}-E_C$ (rad·GHz).
- `anharmonicity() -> float` — the anharmonicity $\alpha=-E_C$ (rad·GHz).
- `sensitivity(flux=None, delta=1e-6) -> float` — central-difference numerical
  $\mathrm{d}\omega_{01}/\mathrm{d}\Phi$ (rad·GHz/$\Phi_0$).
- `with_flux(new_flux) -> QubitSpec` — returns a **new** `QubitSpec` with
  `flux_bias=new_flux`, leaving `self` unchanged.
- `optimal_work_point() -> float` (staticmethod) — the flux where
  $\lvert\mathrm{d}f/\mathrm{d}\Phi\rvert$ is maximal, in $\Phi_0$, equal to
  $\arctan(\sqrt2)/\pi$.

## TransmonQubit — the tunable Transmon device

The **runtime device object** that enters the experimental workflow. Internally
it holds a `QubitSpec` and pre-builds the QuTiP operators, Hamiltonian, and
collapse operators, implementing the `Device` interface (`name`/`hilbert_dim`/
`hamiltonian_static`/`collapse_operators`). The `experiments`, `control`,
`calibration`, and `simulation` layers **accept a `TransmonQubit` by default**,
and the default qubit produced by `sqc.config.CONFIG` is of this type. By
contrast, `QubitSpec` describes physical parameters only; to enter the
experimental workflow, use `TransmonQubit`.

**Construction**: `TransmonQubit(EC, EJ, T1, T2, flux=0.0, state=0, n_levels=3,
name="Q")`. Note the argument here is `flux` (not `flux_bias`), and an initial
`state` may be passed (a Fock index or a `Qobj`).

**Key attributes** (readable after construction): `EC`, `EJ` (the effective
value at the current flux), `EJ_0` (the zero-flux value), `flux`, `n_levels`,
`T1`, `T2`, `frequency`, `anharmonicity`, the operators `a`/`a_dag`/`n`/`I`,
`hamiltonian`, `c_ops`, `state`.

**Accessors and the Device interface**

- `spec() -> QubitSpec` — extract the internal `QubitSpec`.
- `name` (property), `hilbert_dim()`, `hamiltonian_static()`,
  `collapse_operators()` — the `Device` contract.
- `calculate_frequency() -> float` / `calculate_anharmonicity() -> float` —
  compute $\omega_{01}$ / $\alpha$ from the current `EJ, EC` (rad·GHz).
- `frequency_sensitivity(flux, delta_flux=1e-6) -> float` — central-difference
  $\mathrm{d}\omega_{01}/\mathrm{d}\Phi$ (rad·GHz/$\Phi_0$).

**Hamiltonian and collapse operators**

- `get_hamiltonian() -> Qobj` — the lab-frame time-independent $H$.
- `get_hamiltonian_rwa(omega_d) -> Qobj` — the rotating-frame (RWA)
  $H=\Delta n+\tfrac{\alpha}{2}(n^2-n)$, $\Delta=f_{01}-\omega_d$.
- `get_collapse_operators() -> list[Qobj]` — generates
  $[\sqrt{\gamma_1}\,a,\ \sqrt{\gamma_\phi}\,n]$ from `T1, T2`.

```{note}
Placing a qubit in a flux signal and building the time-dependent Hamiltonian is
the job of the `simulation` layer's `HamiltonianBuilder`; see {doc}`simulation`.
```

**Noise, state projection, and gate simulation**

- `generate_1f_noise(t_lists, amplitude, f_min, f_max) -> np.ndarray` —
  generates a 1/f noise time series (delegates to `sqc.simulation.noise`).
- `calculate_state_projection(target_state) -> float` — the probability that the
  current `state` lands on `target_state` (a Fock index or a `Qobj`).
- `ideal_gate(theta, phi) -> Qobj` — the ideal single-qubit rotation on the
  $\{|0\rangle,|1\rangle\}$ subspace (a $2\times2$ unitary).
- `simulate_gate(theta, phi, T, sigma) -> (Qobj, float, float)` — DRAG
  single-qubit gate simulation, returning `(U_eff, leakage, fidelity)`.

## Resonator — the multimode resonator

A multimode readout/coupling resonator, subclassing `Device`. Each mode is a
Fock space, and the operators are built from the overall tensor product.

**Construction**: `Resonator(frequencies, n_levels)` — `frequencies` is a list
of per-mode frequencies (rad·GHz), `n_levels` a list of per-mode Fock truncation
dimensions (the two are equal length).

**Attributes**: `frequencies`, `M` (number of modes), `n_levels`, the mode
operator lists `a`/`adag`/`n`, the overall identity operator `I`, `hamiltonian`.

**Methods**

- `name` (property), `hilbert_dim()`, `hamiltonian_static()`,
  `collapse_operators()` (currently returns an empty list) — the `Device`
  contract.
- `get_hamiltonian() -> Qobj` — the linear Hamiltonian $\sum_i f_i n_i$.
- `get_hamiltonian_rot(omega_d) -> Qobj` — the rotating frame
  $\sum_i (f_i-\omega_d) n_i$.

## CoupledSystem — the two-qubit coupled system

Two Transmons coupled into a multimode resonator via a tunable coupler,
subclassing `Device`, used for two-qubit gate simulation.

**Construction**: `CoupledSystem(qubit1, qubit2, cavity)` — two `TransmonQubit`s
plus one `Resonator`.

**Attributes**: the tensor-product operators `a1`/`a2`/`acav`, coupling strength
`g` (2×M), sub-Hamiltonians `H_q1`/`H_q2`/`H_cav`/`H_0`, and the total `H`.

**Methods**

- `name` (property), `hilbert_dim()`, `hamiltonian_static()`,
  `collapse_operators()` — the `Device` contract.
- `initialize_g() -> list[list[float]]` — initialize the default qubit-mode
  coupling strengths.
- `control_g(g, is_time_dependent=False, f=None) -> None` — set the (optionally
  time-dependent) qubit-cavity coupling and rebuild `H`.
- `build_system(is_rwa=False, omega_d=None) -> None` — build the full system
  Hamiltonian (optionally RWA).
- `prepare_ket11(T, sigma) -> (Qobj, float)` — prepare $|11\rangle$ with a DRAG
  π pulse, returning `(final_state, fidelity)`.
- `simulate_iSWAP(g) -> (Qobj, float, float)` — iSWAP gate simulation under
  tunable coupling, returning `(U_cal, leakage, F_loc)`.

## ChipTopology — the multi-qubit chip topology

A `@dataclass` serving as a container for all devices on the whole chip and
responsible for **lifting** single-device operators into the full Hilbert space.
Subclasses `Device`.

**Fields**: `qubits` (list of qubits), `resonators` (list of resonators),
`couplings` (a `(qubit_name, resonator_name) -> g` dict), `control_lines` (a name
→ `ControlLine` dict), `transfer_matrix` (a frequency-dependent `TransferMatrix`,
may be `None`).

**Methods**

- `name` (property), `hilbert_dim()` — the whole-chip total dimension (the
  product of subsystem dimensions).
- `lift_qubit_op(op, qubit_name) -> Qobj` — embed a single-qubit operator into
  the whole-chip space (`op` in the target slot, identity elsewhere); raises
  `ValueError` if the qubit is not found.
- `hamiltonian_static() -> Qobj` — the sum of each device's static Hamiltonian
  (note: the base implementation **excludes** coupling terms; use `CoupledSystem`
  for coupled systems).
- `collapse_operators() -> list[Qobj]` — collects and embeds all devices'
  collapse operators.
- `get_qubit(name) -> TransmonQubit` — get a qubit by name; raises `ValueError`
  if not found.

## Minimal example

```python
import numpy as np
from sqc.devices import QubitSpec, TransmonQubit

# Preferred: the immutable parameter record (EC/EJ in rad·GHz, incl. the 2π factor)
spec = QubitSpec(name="Q0", EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
                 T1=10_000, T2=8_000)

print(spec.frequency() / (2 * np.pi))          # sweet-spot frequency (GHz)
print(spec.anharmonicity() / (2 * np.pi))       # anharmonicity ≈ -0.2 GHz

# Move to the most sensitive work point (returns a new object, original unchanged)
sweet = QubitSpec.optimal_work_point()           # Φ_0
biased = spec.with_flux(sweet)
print(biased.sensitivity() / (2 * np.pi))        # df/dΦ (GHz/Φ_0)

# The runtime device object: the type the experiment layers accept by default,
# holding ready-made operators and Hamiltonian
qubit = TransmonQubit(EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
                      T1=10_000, T2=8_000)
H = qubit.get_hamiltonian()                       # time-independent H (Qobj)
c_ops = qubit.get_collapse_operators()            # [√γ1·a, √γφ·n]
```

## Physical role / extension

- `TransmonQubit` / `QubitSpec` correspond to one physical qubit on a real chip;
  their parameters come directly from device characterization or fabrication
  parameters (recommended ranges in `sqc.config`: $E_J/E_C\in[40,80]$,
  $f_{01}\in[4,8]\,\mathrm{GHz}$).
- `Resonator` corresponds to a readout/coupling resonator; `CoupledSystem`
  corresponds to the two-qubit + coupler + cavity experimental unit for two-qubit
  gates; `ChipTopology` corresponds to the layout and connectivity of a whole
  multi-qubit chip.
- **To add a custom device type** (e.g. fluxonium, a multi-junction qubit),
  subclass the `Device` abstract base, implement the four members `name`,
  `hilbert_dim()`, `hamiltonian_static()`, `collapse_operators()`, and keep the
  "pure parameter object" convention (store no experimental state). See the full
  extension guide in {doc}`../extending`.
