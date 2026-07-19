# sqc — Superconducting Qubit Quantum Sensing Simulation Platform

**English** | [中文](README.zh-CN.md)

[![tests](https://github.com/Serendipity-Zip/Sensing-Project/actions/workflows/tests.yml/badge.svg)](https://github.com/Serendipity-Zip/Sensing-Project/actions/workflows/tests.yml)
![python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)
![license](https://img.shields.io/badge/license-MIT-green)

A QuTiP-based, full-stack simulation platform for **time-dependent magnetic-field
sensing** with superconducting Transmon qubits. It models the flux → qubit-frequency
transduction chain and reconstructs the sensed waveform from simulated measurements,
following the six-layer cQED architecture of Gao, Rol, Touzard & Wang (2021).

- **Version:** 1.0.0 &nbsp;·&nbsp; **License:** MIT &nbsp;·&nbsp; **Python:** ≥ 3.10
- Natural units throughout (ħ = 1): frequencies/energies in GHz, time in ns, flux in Φ₀.

## What v1 provides

Three end-to-end research pipelines plus an interactive web demo:

1. **Waveform reconstruction** — sense an unknown time-dependent flux via a sliding
   transient protocol and reconstruct it (Wiener deconvolution, Hammerstein–Wiener,
   or Levenberg–Marquardt), plus Ramsey / differential-echo / cryoscope protocols.
2. **Qubit frequency calibration** — Ramsey-based `f(Φ)` / `f₀₁` calibration.
3. **Waveform predistortion** — measure control-line distortion and design an IIR/FIR
   inverse filter, with an end-to-end validation workflow.

A layered API (`sqc.devices → control → hardware → simulation → experiments →
reconstruction → calibration → workflows`) lets you work at whatever level you need.

## Install

Requires Python ≥ 3.10. Editable install makes `import sqc` work from anywhere:

```bash
pip install -e .            # core (numpy, scipy, qutip, matplotlib)
pip install -e ".[demo]"    # + Gradio web demo
pip install -e ".[test]"    # + pytest tooling
```

Or, for a plain dependency install: `pip install -r requirements.txt`.

## Quick start

The unified entry point is `SensingWorkflow` — configure, run, and plot:

```python
from sqc.workflows import SensingWorkflow

wf = SensingWorkflow()
wf.configure(
    protocol="transient",       # ramsey | echo | transient | cryoscope | delay_ramsey
    signal_amplitude=0.01,      # Φ₀
    reconstruction="wiener",    # wiener | hammerstein | lm
    t_rabi_duration=20,         # ns
)
result = wf.run()               # measure + reconstruct
wf.plot()                       # matplotlib figure
```

Sweep a parameter or compare reconstruction methods:

```python
sweep = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02, 0.05])

wf.run(measure=True, reconstruct=False)
cmp = wf.compare(methods=["wiener", "hammerstein", "lm"])
print("best:", cmp.best)
```

### Lower-level API

Every layer is usable directly, e.g.:

```python
from sqc.devices.transmon import TransmonQubit          # EC, EJ, T1, T2, flux
from sqc.experiments import TransientSensingExperiment
from sqc.reconstruction import KernelEstimator, TransientReconstruction
```

See [`docs/architecture.md`](docs/architecture.md) for the full module reference.

## Web demo

An interactive Gradio UI (device config, protocols, reconstruction, predistortion):

```bash
pip install -e ".[demo]"
python web_demo_v2.py
```

## Documentation

- **Tutorial notebook:** [`Simulation_sqc.ipynb`](Simulation_sqc.ipynb) — Rabi → Ramsey →
  echo → transient → cryoscope, worked end to end.
- **Architecture / module reference & roadmap:** [`docs/architecture.md`](docs/architecture.md)
  (中文, six-layer stack, per-module API, extension guide, §A2 post-v1 roadmap).

## Testing

```bash
pip install -e ".[test]"
pytest tests/ -v
pytest tests/regression -m regression      # physics regression baselines
```

## Project layout

```
sqc/                  the platform (devices, control, hardware, simulation,
                      experiments, reconstruction, calibration, workflows)
tests/                unit / integration / equivalence / regression suites
web_demo_v2.py        Gradio web demo (sqc-based)
Simulation_sqc.ipynb  end-to-end tutorial notebook
docs/                 architecture / technical documentation
src/, src_mirror/     frozen legacy reference implementation + facade
```

## Not in v1 (planned)

The following are present in the codebase but **hidden from the v1 public API** and
scheduled for a later release: two-qubit Z-crosstalk extraction, transient-based
frequency calibration, the CPMG protocol, the tunable coupler and electronics
(AWG/ADC) hardware layers. See [`docs/architecture.md`](docs/architecture.md) §A2
for the roadmap.

## Citation

If you use this platform in academic work, please cite it — see
[`CITATION.cff`](CITATION.cff).

## License

MIT © 2026 Zip. See [`LICENSE`](LICENSE).
