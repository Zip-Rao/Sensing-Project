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
- **Transient flux-response calibration** — fitting the $\Delta\omega(\Phi)$
  dispersion curve from an unknown transient signal, i.e.
  `FluxResponseCalibration(method="transient")`. This one method raises
  `NotImplementedError`. Transient single-point frequency *measurement*
  (`FrequencyMeasurement`, `SinglePointFrequencyCalibration`, and
  `FrequencyCalibrationWorkflow`) is implemented and part of the v1 public API.
  `FrequencyCalibrationWorkflow` accepts an arbitrary multi-stage pipeline via
  `stages=[CalibrationStage(...), ...]` (each stage picks its own measurement
  method, stepper, tolerance, and measurement axes); the transient→Ramsey
  hybrid is the default 2-stage preset. Gradient stages additionally support an
  `advance_when` predicate for adaptive stage switching on runtime signals a
  residual threshold cannot express.
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
Two transient capabilities should not be confused. Transient waveform
*reconstruction* ({py:class}`~sqc.experiments.TransientSensingExperiment`) and
transient single-point frequency *measurement* are both core v1 features and are
fully supported. Only the transient *flux-response calibration* above — fitting
the $\Delta\omega(\Phi)$ curve from an unknown signal — is deferred.
```
