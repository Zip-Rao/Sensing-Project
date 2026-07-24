# Overview

**sqc** is a QuTiP-based, composable full-stack simulation framework for
time-dependent magnetic-field sensing with superconducting transmon qubits. It
models the flux → qubit-frequency transduction chain across an eight-layer cQED
stack and reconstructs the sensed waveform from simulated measurements.

## A framework, not three apps

The platform ships three end-to-end research pipelines:

1. **Waveform reconstruction** — sense an unknown time-dependent flux and
   reconstruct it (Wiener deconvolution, Hammerstein–Wiener, or
   Levenberg–Marquardt), plus Ramsey / differential-echo / cryoscope protocols.
2. **Frequency calibration** — Ramsey-based $f(\Phi)$ / $f_{01}$ calibration.
3. **Predistortion** — model the AWG → chip transfer function and design a
   compensation filter.

These are **worked examples of composing the layers** — not the boundary of
what the platform can do. Every layer (devices, hardware, control, simulation,
experiments, reconstruction, calibration, workflows) exposes an extension
interface. You can subclass any of them to build your own sensing application:
a new protocol, a new inversion algorithm, a new device type. See
{doc}`extending`.

## Key characteristics

- **Pure Python + QuTiP** — full density-matrix time evolution via `mesolve`.
- **Natural units** (ħ = 1): rad·GHz, ns, Φ₀.
- **Immutable device parameters** and **side-effect-free** Hamiltonian
  construction.
- **Centralised configuration** — one source of truth,
  {py:data}`sqc.config.CONFIG`.
- **Eight-layer architecture** following Gao et al. (2021), the cQED
  engineering standard.

## Next steps

- {doc}`install` — get it running.
- {doc}`quickstart` — a three-line sensing pipeline.
- {doc}`architecture` — the eight-layer stack.
- {doc}`building_blocks/index` — each layer in detail.
