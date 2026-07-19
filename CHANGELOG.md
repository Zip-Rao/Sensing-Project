# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and this project adheres to
[Semantic Versioning](https://semver.org/).

## [1.0.0] — 2026-07-17

First official release. Public API surface = the `sqc` package (three sensing/
predistortion pipelines + Gradio web demo).

### Added
- `pyproject.toml` — `pip install -e .` support; `import sqc` no longer depends on
  the working directory. Optional extras: `[demo]` (Gradio), `[test]` (pytest).
- `LICENSE` (MIT), `README.md` (English) + `README.zh-CN.md` (中文, kept in
  sync via cross-links), `CITATION.cff`.
- Public package APIs: `sqc.control` and `sqc.simulation` now re-export their
  classes; `sqc.hardware` exports `CascadeDistortion`.

### Changed
- Version unified to **1.0.0** (`sqc.__version__`, single source; docs updated).
- `requirements.txt` trimmed to the core set (numpy/scipy/qutip/matplotlib);
  removed unused `pandas` / `scikit-learn`; demo deps live in
  `requirements_demo.txt` / the `[demo]` extra.
- TODO tracker system tidied: `RELEASE_TODO.md` is the release checklist,
  `idea/refactor/_handoff_state.md` the engineering-progress source of truth;
  stale trackers annotated.

### Fixed
- `EchoReconstruction`: clip the `arcsin` argument to avoid silent NaNs when
  `p_e` drifts marginally outside [0, 1].
- `create_pulse`: rescale the DRAG `Omega_Q` envelope together with the I
  envelope; guard against a zero base-angle division.
- `create_ramsey_pulse`: advance by pulse duration `t_rabi[-1] - t_rabi[0]`,
  consistent with the echo factories.
- `SingleExponentialDistortion` bilinear inverse: identity fallback as
  amplitude → 1 (matches the `rol2020` path), avoiding an unstable inverse.
- `PredistortionDesigner._auto_method`: removed a duplicated `@staticmethod`.

### Not included (planned for a later release)
- Two-qubit Z-crosstalk extraction, transient-based frequency calibration, the
  CPMG protocol, and the tunable-coupler / electronics (AWG/ADC) hardware layers
  are present in the codebase but hidden from the v1 public API. See
  `RELEASE_TODO.md`.
