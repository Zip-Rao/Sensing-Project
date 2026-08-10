# Workflows

## What this layer provides

The `workflows` layer is the stack's top-level orchestration layer. It
composes the {doc}`experiments`, {doc}`reconstruction`, and {doc}`calibration`
layers into complete research pipelines and exposes a high-level API to the user.
It introduces no new physics; it only chains and batches the
configure, measure, reconstruct, compare, and visualise steps.

`SensingWorkflow` is the platform's single user-facing entry point: one class
where `configure()` sets parameters, `run()` runs the pipeline, `sweep()` scans a
parameter, `compare()` compares reconstruction algorithms, and `plot()`
auto-visualises. `PredistortionValidationWorkflow` is the end-to-end validation
workflow for the predistortion mainline.

### Two ways to do frequency calibration

Frequency calibration has three related classes, organised into **old and new APIs**:

| | Old API (one-shot) | New API (event-driven) |
|---|---|---|
| **Entry point** | `FrequencyCalibrationWorkflow` | `FrequencyCalibrationRuntime` |
| **Core** | Chains multiple `SinglePointFrequencyCalibration` internally | Drives `FrequencyStateMachine` internally |
| **How to run** | `wf.run()` — runs to completion | `runtime.run()` — stepwise event loop |
| **Pause/resume?** | No | Yes (`save_run` / `load_run` persistence) |
| **Recovery?** | No (stages only go forward) | Yes (Track→Reacquire→Track, Lock→Verify→Lock) |
| **Shared component** | `DampedSecantTracker` (step formula) | Same `DampedSecantTracker` |

```text
User
 │
 ├── (old) FrequencyCalibrationWorkflow.run()
 │       internal: stage1 → stage2 → ...   one-shot pipeline
 │       each stage delegates to DampedSecantTracker
 │
 └── (new) FrequencyCalibrationRuntime.run()
            │
            ├── FrequencyStateMachine  ← pure state transitions (no QuTiP)
            │     six states: Acquire → Track → Verify → Lock
            │     recovery branches: Reacquire, SafeStop
            │     DampedSecantTracker ← shared step formula
            │
            └── SQCExecutor ← command → QuTiP measurement
```

**Which one to use**:

- Want one-shot, no mid-run intervention → use `FrequencyCalibrationWorkflow`
- Want stepwise control, recovery, persistence, long-run simulation → use
  `FrequencyCalibrationRuntime` + `FrequencyStateMachine`

Both share the same `DampedSecantTracker` control law — numerical behaviour is
identical.

## Class overview

**Extension point**

| Class | Role |
|---|---|
| `Workflow` | Abstract base, the common contract for all workflows (`run() -> dict`); this layer's extension point |

**Workflows**

| Class | Role |
|---|---|
| `SensingWorkflow` | Unified research entry: configure/execute/sweep/compare/visualise |
| `PredistortionValidationWorkflow` | End-to-end predistortion validation: inject distortion → calibrate → design inverse filter → verify improvement |
| `FrequencyCalibrationWorkflow` | Multi-stage closed-loop frequency calibration (legacy staged workflow) |
| `FrequencyStateMachine` | Event-driven frequency calibration state machine V2: Acquire→Track→Verify→Lock |
| `FrequencyCalibrationRuntime` | State machine runtime: event loop + tracker + persistence |

**Result containers (active)**

| Class | Role |
|---|---|
| `WorkflowResult` | Single `run()` result (measurement, reconstructed signal, calibration table, config snapshot) |
| `SweepResult` | `sweep()` result (one `WorkflowResult` per value + metrics) |
| `CompareResult` | `compare()` result (per-method signals + metrics + best) |

**Result containers (placeholders for future features)**

| Class | Corresponding (unimplemented) method |
|---|---|
| `DiffReport` / `BenchmarkResult` / `NoiseReport` / `CVResult` | Diff report / benchmark / noise characterization / cross-validation (post-v1) |

## Workflow — workflow abstract base class

The common contract for all workflows, this layer's extension point. It
mandates a single abstract method:

- `run() -> dict`: execute the workflow and return a dict of named results.

To add a custom workflow, subclass `Workflow` and implement `run()`. See
{doc}`../extending`.

## SensingWorkflow — unified research entry point

The platform's single user-facing API: one class chaining configuration,
protocol execution, sweeping, A/B comparison, and visualisation, internally
delegating to the experiment / reconstruction / calibration layers.

**Construction**

`SensingWorkflow()`

**Methods**

- `configure(**params) -> self`: set parameters by functional group (only
  explicitly-passed values change; the rest keep their current value), chainable.
  Groups: Qubit (`EC`/`EJ`/`T1`/`T2`/`flux_bias`/`n_levels`; `EC`/`EJ` in
  GHz), Protocol & Signal (`protocol`, `signal_*`), Reconstruction
  (`reconstruction` method name plus `lambda_reg`/`lm_*`), Pulse
  (`t_rabi_duration`/`t_global_*`, `rotation_angle`/`rabi_rate`,
  `pulse_envelope`/`envelope_sigma`, and both pulse phases), Hardware
  (`sample_rate`, GSa/s).
- `run(measure=True, reconstruct=True, calibrate=False) -> WorkflowResult`: run
  a full sensing pipeline. `calibrate=True` is supported only for the
  `cryoscope`/`delay_ramsey` protocols (calibrate before measurement).
- `sweep(param_path, values) -> SweepResult`: sweep a dotted path `"group.field"`
  (e.g. `"signal.amplitude"`) over a list of values, running once per point and
  computing SNR/RMSE/peak metrics.
- `compare(methods, ...) -> CompareResult`: compare multiple reconstruction
  algorithms on the same measurement, picking the best by RMSE.
- `plot()`: auto-dispatch visualisation based on the most recent
  `run`/`sweep`/`compare`.

Protocol → (Experiment, Reconstruction) mapping: `"transient"`, `"ramsey"`,
`"echo"`, `"cryoscope"`, `"delay_ramsey"` each map to one experiment/reconstruction
pair.

**Output**

- `run()` → `WorkflowResult` (with `measurement`, `reconstructed_signal`, `config_snapshot`)
- `sweep()` → `SweepResult` (with `results`, `metrics`)
- `compare()` → `CompareResult` (with `signals`, `metrics`, `best`)

```{note}
Several advanced methods (batch benchmarking, noise characterization,
cross-validation, multi-qubit, persistence, etc.) are post-v1 placeholders that
currently raise `NotImplementedError` and are not part of the v1 public surface.
v1 covers single-pipeline runs, sweeps, and algorithm comparison via
`run`/`sweep`/`compare`.
```

## FrequencyCalibrationWorkflow — multi-stage frequency calibration (legacy)

Orchestrates `SinglePointFrequencyCalibration` as a multi-stage pipeline
(e.g. transient→Ramsey hybrid), with each stage seeded from the previous
stage's optimum. Internally delegates to `DampedSecantTracker`.

**Construction**

`FrequencyCalibrationWorkflow(qubit, f_target, stages=None, V_seed=None, ...)`

When `stages` is not given, a default two-stage hybrid is built:
coarse (transient+gradient) → fine (ramsey). See {doc}`calibration` for details.

## FrequencyStateMachine — event-driven frequency calibration V2

v2.19 introduces an **event-driven** six-state protocol:
Acquire → Track → Verify → Lock, with a Reacquire recovery branch and
SafeStop safety hold.

**Design principles**:
- Pure state machine, no QuTiP dependency — transition table is testable offline
- Only **Track** may propose a new working flux bias
- **Verify** freezes the candidate bias from entry to exit; Lock also never
  changes bias directly
- Lock drift detection does not adjust bias — small drift → Verify, large
  jump → Reacquire
- Analytic $f(\Phi)$ is a simulation oracle only — never feeds the transition
  reducer

**Usage** (via `FrequencyCalibrationRuntime`):

```python
from sqc.workflows.frequency_runtime import FrequencyCalibrationRuntime
from sqc.workflows.frequency_state_machine import FrequencyCalibrationConfig

config = FrequencyCalibrationConfig(
    epsilon_enter=2*np.pi*20e-3, epsilon_final=2*np.pi*2e-3,
    N_verify=2, max_commands=30,
)
runtime = FrequencyCalibrationRuntime(qubit=q, f_target=f_target, config=config)
result = runtime.run()
# result["state"] → "lock", result["run_status"] → "calibrated"
```

See {doc}`calibration` "Frequency calibration state machine V2" for details.

## PredistortionValidationWorkflow — end-to-end predistortion validation

The full validation workflow for the predistortion mainline. Given a target
waveform and a known distortion, it runs: inject distortion into the
control line → measure the transfer function → fit a model → design an inverse
filter → apply predistortion → re-measure → compute improvement metrics.

**Construction**

`PredistortionValidationWorkflow(target_waveform, true_distortion, designer=None, control_line_params=None, qubit=None, measurement_protocol=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `target_waveform` | `Waveform` | Desired on-chip waveform | — |
| `true_distortion` | `DistortionModel` | Ground-truth distortion model | — |
| `designer` | `PredistortionDesigner` | Inverse-filter designer | `PredistortionDesigner(method="auto")` |
| `control_line_params` | dict | Control-line construction parameters | `None` |
| `qubit` | `TransmonQubit` | Qubit (for protocol-driven measurement) | `None` |
| `measurement_protocol` | str | Measurement protocol (`None` = analytical path) | `None` |

**Methods**

- `run() -> dict`: execute the full validation pipeline.

**Output**

Returns `dict` with:

| Key | Type | Meaning |
|---|---|---|
| `target` | `Waveform` | Target waveform |
| `on_chip_uncorrected` | `Waveform` | Uncorrected on-chip waveform |
| `on_chip_corrected` | `Waveform` | Corrected on-chip waveform |
| `awg_predistorted` | `Waveform` | Predistorted AWG waveform |
| `inverse_model` | `DistortionModel` | Inverse filter model |
| `measured_model` | `DistortionModel` | Measured & fitted transfer-function model |
| `metrics` | dict | `rmse_uncorrected`, `rmse_corrected`, `improvement_factor`, settling times |

## Result containers

### WorkflowResult

Single `run()` result.

**Fields**

| Field | Type | Meaning |
|---|---|---|
| `config_snapshot` | dict | Configuration snapshot |
| `measurement` | `ExperimentResult` | Raw measurement data |
| `reconstructed_signal` | `FluxSignal` | Reconstructed flux signal |
| `reconstruction_details` | dict | Reconstruction algorithm details |
| `calibration` | `CalibrationTable` | Calibration table (optional) |

### SweepResult

`sweep()` result.

**Fields**

| Field | Type | Meaning |
|---|---|---|
| `param_path` | str | Swept parameter path |
| `values` | list | Swept value list |
| `results` | `list[WorkflowResult]` | One result per value |
| `metrics` | dict | Key → per-point metric list (e.g. `snr`, `rmse`) |

### CompareResult

`compare()` result.

**Fields**

| Field | Type | Meaning |
|---|---|---|
| `methods` | `list[str]` | Method names compared |
| `signals` | dict | Method name → `FluxSignal` |
| `metrics` | dict | Method name → metric dict |
| `best` | str | Lowest-RMSE method name |

### Placeholder containers

`DiffReport` / `BenchmarkResult` / `NoiseReport` / `CVResult`: the return
types of the diff-report, batch-benchmark, noise-characterization, and
cross-validation methods respectively; those methods are post-v1 placeholders,
so these are reserved interfaces for now.

## Minimal example

```python
import numpy as np
from sqc.workflows import SensingWorkflow
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.workflows import PredistortionValidationWorkflow

# ── 1. SensingWorkflow: configure + single run ───────────────────────
wf = SensingWorkflow()
wf.configure(
    protocol="transient",
    signal_amplitude=0.01,
    reconstruction="wiener",
    flux_bias=0.1,
)
result = wf.run()          # WorkflowResult
print(result.reconstructed_signal)   # FluxSignal

# ── 2. Sweep (signal amplitude) ──────────────────────────────────────
sweep = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02, 0.05])
print(sweep.metrics["snr"])          # SNR per amplitude

# ── 3. A/B compare reconstruction methods ────────────────────────────
wf.run(measure=True, reconstruct=False)   # measure only, cache data
cmp = wf.compare(methods=["wiener", "hammerstein", "lm"])
print(f"Best method: {cmp.best}")

# ── 4. PredistortionValidationWorkflow ────────────────────────────────
t = np.linspace(0, 500, 1000)
target = Waveform(t_list=t, samples=np.ones_like(t) * 0.05)
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)
val_wf = PredistortionValidationWorkflow(
    target_waveform=target,
    true_distortion=distortion,
)
val_result = val_wf.run()
print(f"Improvement factor: {val_result['metrics']['improvement_factor']:.1f}×")
```

## Physical role / extension

- This layer corresponds to the upper-level experiment scripts in a lab:
  pick a protocol, tune parameters, run a batch of data, produce results. It
  chains all seven layers
  {doc}`devices`→{doc}`hardware`→{doc}`control`→{doc}`simulation`→{doc}`experiments`
  →{doc}`reconstruction`→{doc}`calibration`, the convenient entry point that lets
  users avoid touching the lower layers directly.
- `FrequencyCalibrationWorkflow` is the legacy multi-stage frequency calibration
  interface; internally delegates to `DampedSecantTracker`.
- `FrequencyStateMachine` + `FrequencyCalibrationRuntime` form the V2
  event-driven interface: a six-state protocol with full guards, budgets, and
  monitor hysteresis, supporting `save_run`/`load_run` persistence and
  checkpoint recovery.
- `SensingWorkflow.sweep()` manages CONFIG synchronisation automatically (calling
  `_sync_config()` when pulse/hardware-group parameters change), so each `run()`
  uses the time axis matching the current parameters.
- For the transient protocol, `pulse.rotation_angle` and `pulse.rabi_rate` are
  mutually exclusive modes; sweeping either automatically clears the other.
  For a fixed-drive small-angle scan, set `rabi_rate` and sweep
  `pulse.t_rabi_duration`. Each point rebuilds the pulse and re-estimates its
  response kernel.
- `compare()` reuses the cached measurement across multiple reconstruction
  algorithms, avoiding re-running `mesolve`; the `best` field picks by RMSE, which
  is meaningful when a ground truth (the simulation's `flux_samples`) is present,
  and falls back to the most NaN-robust method otherwise.
- To add a custom workflow, subclass the `Workflow` abstract base and
  implement `run() -> dict`, orchestrating the experiment/reconstruction/calibration
  layers internally as needed. Full extension guide in {doc}`../extending`.
