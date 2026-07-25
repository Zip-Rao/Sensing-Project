# Workflows

## What this layer provides

The `workflows` layer is the stack's **top-level orchestration layer** — it
composes the {doc}`experiments`, {doc}`reconstruction`, and {doc}`calibration`
layers into complete research pipelines and exposes a high-level API to the user.
It introduces no new physics; it only chains and batches
"configure → measure → reconstruct → compare → visualise".

`SensingWorkflow` is the platform's **single user-facing entry point**: one class
where `configure()` sets parameters, `run()` runs the pipeline, `sweep()` scans a
parameter, `compare()` compares reconstruction algorithms, and `plot()`
auto-visualises. `PredistortionValidationWorkflow` is the end-to-end validation
workflow for the predistortion mainline.

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

The common contract for all workflows, this layer's **extension point**. It
mandates a single abstract method:

- `run() -> dict` — execute the workflow and return a dict of named results.

To add a custom workflow, subclass `Workflow` and implement `run()`. See
{doc}`../extending`.

## SensingWorkflow — unified research entry point

The platform's single user-facing API — one class chaining configuration,
protocol execution, sweeping, A/B comparison, and visualisation, internally
delegating to the experiment / reconstruction / calibration layers. Core methods:

- `configure(**params) -> self` — set parameters by functional group (only
  explicitly-passed values change; the rest keep their current value), chainable.
  Groups: **Qubit** (`EC`/`EJ`/`T1`/`T2`/`flux_bias`/`n_levels`; `EC`/`EJ` in
  GHz), **Protocol & Signal** (`protocol`, `signal_*`), **Reconstruction**
  (`reconstruction` method name plus `lambda_reg`/`lm_*`), **Pulse**
  (`t_rabi_duration`/`t_global_*`, in ns), **Hardware** (`sample_rate`, GSa/s).
- `run(measure=True, reconstruct=True, calibrate=False) -> WorkflowResult` — run
  a full sensing pipeline. `calibrate=True` is supported only for the
  `cryoscope`/`delay_ramsey` protocols (calibrate before measurement).
- `sweep(param_path, values) -> SweepResult` — sweep a dotted path `"group.field"`
  (e.g. `"signal.amplitude"`) over a list of values, running once per point and
  computing SNR/RMSE/peak metrics.
- `compare(methods, ...) -> CompareResult` — compare multiple reconstruction
  algorithms on the same measurement, picking the best by RMSE.
- `plot()` — auto-dispatch visualisation based on the most recent
  `run`/`sweep`/`compare`.

Protocol → (Experiment, Reconstruction) mapping: `"transient"`, `"ramsey"`,
`"echo"`, `"cryoscope"`, `"delay_ramsey"` each map to one experiment/reconstruction
pair.

```{note}
Several advanced methods (batch benchmarking, noise characterization,
cross-validation, multi-qubit, persistence, etc.) are post-v1 placeholders that
currently raise `NotImplementedError` and are not part of the v1 public surface.
v1 covers single-pipeline runs, sweeps, and algorithm comparison via
`run`/`sweep`/`compare`.
```

## PredistortionValidationWorkflow — end-to-end predistortion validation

The full validation workflow for the predistortion mainline (`@dataclass`). Given
a target waveform and a known distortion, it runs: inject distortion into the
control line → measure the transfer function → fit a model → design an inverse
filter → apply predistortion → re-measure → compute improvement metrics. Fields:
`target_waveform` (desired on-chip waveform), `true_distortion` (ground-truth
distortion model), `designer` (inverse-filter designer, default
`PredistortionDesigner(method="auto")`), `control_line_params`; plus `qubit` and
`measurement_protocol` for the protocol-driven measurement (`None` selects the
analytical path).

`run() -> dict` returns `target`, `on_chip_uncorrected`, `on_chip_corrected`,
`awg_predistorted`, `inverse_model`, `measured_model`, `metrics`. `metrics`
includes `rmse_uncorrected`, `rmse_corrected`, `improvement_factor`, and two
settling times.

## Result containers

- `WorkflowResult` — single `run()` result: `config_snapshot`, `measurement`
  (`ExperimentResult`), `reconstructed_signal` (`FluxSignal`),
  `reconstruction_details`, `calibration` (`CalibrationTable`).
- `SweepResult` — `sweep()` result: `param_path`, `values`, `results` (one
  `WorkflowResult` per point), `metrics` (key → per-point metric list).
- `CompareResult` — `compare()` result: `methods`, `signals`
  (method → `FluxSignal`), `metrics` (method → metric dict), `best` (lowest-RMSE
  method).
- `DiffReport` / `BenchmarkResult` / `NoiseReport` / `CVResult` — the return
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

- This layer corresponds to the **upper-level experiment scripts** in a lab:
  pick a protocol, tune parameters, run a batch of data, produce results. It
  chains all seven layers
  {doc}`devices`→{doc}`hardware`→{doc}`control`→{doc}`simulation`→{doc}`experiments`
  →{doc}`reconstruction`→{doc}`calibration`, the convenient entry point that lets
  users avoid touching the lower layers directly.
- `SensingWorkflow.sweep()` manages CONFIG synchronisation automatically (calling
  `_sync_config()` when pulse/hardware-group parameters change), so each `run()`
  uses the time axis matching the current parameters.
- `compare()` reuses the cached measurement across multiple reconstruction
  algorithms, avoiding re-running `mesolve`; the `best` field picks by RMSE, which
  is meaningful when a ground truth (the simulation's `flux_samples`) is present,
  and falls back to the most NaN-robust method otherwise.
- **To add a custom workflow**, subclass the `Workflow` abstract base and
  implement `run() -> dict`, orchestrating the experiment/reconstruction/calibration
  layers internally as needed. Full extension guide in {doc}`../extending`.
