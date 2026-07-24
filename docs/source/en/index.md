# sqc — Superconducting Qubit Quantum Sensing Simulation Platform

```{note}
This documentation site is under construction. Pages are being filled in
incrementally; the structure below reflects the planned framework-first layout.
```

**sqc** is a QuTiP-based, composable full-stack simulation framework for
time-dependent magnetic-field sensing with superconducting Transmon qubits.
It models the flux → qubit-frequency transduction chain across an eight-layer
cQED stack and reconstructs the sensed waveform from simulated measurements.

The three built-in research pipelines (waveform reconstruction, frequency
calibration, predistortion) are **worked examples** of composing the layers —
not the boundary of what the platform can do. Every layer exposes an extension
interface so you can build your own sensing applications.

```{toctree}
:maxdepth: 2
:caption: Getting Started

overview
install
quickstart
```

```{toctree}
:maxdepth: 2
:caption: The Framework

architecture
building_blocks/index
```

```{toctree}
:maxdepth: 2
:caption: Using the Platform

examples/index
extending
tutorial
```

```{toctree}
:maxdepth: 2
:caption: Reference

api/index
theory
roadmap
```
