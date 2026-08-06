# Architecture

`sqc` is organised as an **eight-layer full-stack cQED architecture**, following
the six-layer stack of Gao, Rol, Touzard & Wang (2021, *PRX Quantum* 2, 040202),
the de-facto standard for circuit-QED engineering. Each layer maps to a real
physical component or control responsibility in a superconducting qubit system.

## The eight layers

From the bottom (physical device) to the top (research workflow):

```
┌───────────────────────────────────────────────────────────────┐
│ workflows/      Top-level research pipelines                    │  composes ↓
├───────────────────────────────────────────────────────────────┤
│ calibration/    Calibration workflows → CalibrationTable        │
├───────────────────────────────────────────────────────────────┤
│ reconstruction/ Waveform-reconstruction algorithms              │
├───────────────────────────────────────────────────────────────┤
│ experiments/    Experiment protocols (device+flux+seq+readout)  │
├───────────────────────────────────────────────────────────────┤
│ simulation/     QuTiP invocation: Hamiltonian, mesolve, noise   │
├───────────────────────────────────────────────────────────────┤
│ control/        Control pulses: Waveform, FluxSignal, Pulse     │
├───────────────────────────────────────────────────────────────┤
│ hardware/       Control electronics / lines: distortion, xfer   │
├───────────────────────────────────────────────────────────────┤
│ devices/        Physical devices: Transmon, Resonator, Chip     │
└───────────────────────────────────────────────────────────────┘
```

**Dependency rule.** A layer may depend only on layers below it. Reverse
dependencies are forbidden: `reconstruction` must never import `workflows`.
This keeps each layer independently usable and testable.

Each layer exposes an **abstract base class** as its extension point
({py:class}`~sqc.devices.Device`, {py:class}`~sqc.experiments.Experiment`,
{py:class}`~sqc.reconstruction.Reconstruction`, and so on). Subclass it to add
a new device, protocol, or algorithm; see {doc}`extending`.

## Mapping to Gao 2021

| Gao 2021 layer | Covers | `sqc/` module |
|---|---|---|
| Quantum algorithms | Compilation | *(interface only)* |
| Control software | Pulse calibration, sequences | `experiments/`, `calibration/`, `workflows/` |
| Control electronics | AWG, ADC, FPGA | `hardware/electronics.py` |
| Microwave signal proc. | IQ mixer, LO, HEMT, filters | `hardware/distortion.py`, `hardware/readout.py` |
| Cryogenics + interconnect | Control lines, attenuators | `hardware/control_line.py`, `hardware/transfer_matrix.py` |
| Device | Transmon, resonator, SQUID | `devices/` |

The platform focuses on the software-simulation side of the engineering cycle
(Hamiltonian design → simulation → characterisation → feedback). It does **not**
cover chip design, fabrication, or cryogenics.

## Design invariants

- **Natural units** (ħ = 1): frequencies/energies in rad·GHz, time in ns, flux
  in Φ₀.
- **Immutable device parameters**: {py:class}`~sqc.devices.QubitSpec` is a frozen
  dataclass, so physics parameters cannot be polluted by experiment state.
- **Pure Hamiltonian construction**:
  {py:class}`~sqc.simulation.HamiltonianBuilder` is side-effect free.
- **Centralised configuration**: all time grids, AWG parameters, and qubit
  defaults derive from a single source, {py:data}`sqc.config.CONFIG`.
