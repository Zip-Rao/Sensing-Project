# Roadmap

v1 deliberately ships a focused, stable public API: the three research pipelines
(waveform reconstruction, frequency calibration, predistortion) plus the web
demo. The capabilities below exist in the codebase at various maturity levels
but are **not part of the v1 public interface**. They are planned for a future
release and are hidden from the public API, the frontend, and this reference
until then.

## Planned capabilities

- **Z-crosstalk reconstruction** — multi-qubit flux crosstalk characterisation
  and compensation. A numerical scaling issue in the transient reconstruction
  kernel is being resolved before promotion.
- **Transient frequency calibration** — single-point $f_{01}$ calibration from
  an unknown transient signal (the `method="transient"` option on the frequency
  calibration classes). The Ramsey-based path is the supported v1 mainline.
- **CPMG protocol** — Carr–Purcell–Meiboom–Gill dynamical-decoupling sensing.
- **Tunable coupler** — `TunableCoupler` device for two-qubit gate schemes.
- **Electronics layer** — explicit AWG/ADC/LO abstractions
  (`hardware/electronics.py`).
- **SensingWorkflow planned methods** — `save`, `load`, `diff`, `benchmark`,
  `pipeline`, `multi_qubit`, `crosstalk`, `find_optimal_work_point`,
  `detectability_limit`, `noise_characterize`, `cross_validate`. These raise
  `NotImplementedError` today and are omitted from the API reference.

## Known limitations (v1)

- `collapse_operators()` returns an empty list for chip/resonator devices,
  silently dropping dissipation in those configurations.
- Transient waveform reconstruction on short grids recovers only a fraction of
  the amplitude; use the documented grids for quantitative results.

## Infrastructure roadmap

Planned maturity work toward a fully community-facing release: a public test
suite + CI, PyPI distribution, tagged GitHub releases, README badges, a
contributing guide, and a citable Zenodo DOI.

```{note}
Waveform *reconstruction* via the transient protocol
({py:class}`~sqc.experiments.TransientSensingExperiment`) is a **core v1
feature** and is fully supported. Only transient *frequency calibration* is
deferred — do not confuse the two.
```
