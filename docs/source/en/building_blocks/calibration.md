# Calibration

## What this layer provides

The `calibration` layer is the stack's calibration-and-tuning layer. It
measures the current state of the qubit and control line, and when needed
drives it to a target state. It serves two mainlines: frequency calibration
(measure/tune $f_{01}(\Phi)$) and predistortion (measure control-line transfer
functions, design compensation filters). Results are uniformly packaged into a
{py:class}`~sqc.calibration.CalibrationTable` for downstream lookup inversion or
tuning.

Unlike the pure-function constraint of the {doc}`reconstruction` layer, the
calibration classes here actively run simulations to measure: each
`Calibration` subclass implements `calibrate()`, which internally drives the
qubit, runs `mesolve`, fits, and finally returns a `CalibrationTable`.

```{note}
The {doc}`reconstruction` layer also exports two classes that conceptually belong
to calibration, `CryoscopeCalibration` / `DelayRamseyCalibration` (exported
alongside their reconstructors). They likewise subclass this layer's
{py:class}`~sqc.calibration.Calibration` and produce a `CalibrationTable` for the
`inversion="calibration"` lookup path. See {doc}`reconstruction`.
```

## Class overview

Grouped by function into five sets:

**Extension point + result container**

| Class | Role |
|---|---|
| `Calibration` | Abstract base, the common contract for all calibration workflows; this layer's extension point |
| `CalibrationTable` | Calibration result container, with `evaluate` (interpolation) / `inverse` (inverse interpolation) |

**Frequency calibration / measurement / tuning**

| Class | Role |
|---|---|
| `FluxResponseCalibration` | Scan flux to build $f(\Phi)$ lookup (Ramsey point-by-point frequency measurement) |
| `FrequencyMeasurement` | Single-point $f_{01}$ measurement (Ramsey, or the advanced transient method) |
| `SinglePointFrequencyCalibration` | Single-point frequency **tuning**: closed-loop feedback drives $f_q(V)$ to target |
| `FrequencyCalibrationWorkflow` | Orchestration implementing the four-stage state machine for frequency calibration |
| `CalibrationStage` | Per-stage parameter container (used by `FrequencyCalibrationWorkflow`) |

**Waveform / control-line calibration**

| Class | Role |
|---|---|
| `WaveformCalibration` | Unified entry for transfer-function measurement + predistortion design |
| `PredistortionDesigner` | Standalone inverse-filter designer (FIR/IIR/frequency-domain inversion) |

**Scheduling**

| Class | Role |
|---|---|
| `CalibrationScheduler` | Calibration task scheduler: registry + dependency graph, executed in order |

## Calibration: calibration abstract base class

The common contract for all calibration workflows, this layer's extension
point. It mandates a single abstract method:

- `calibrate() -> CalibrationTable`: run the calibration workflow and return
  the result.

To add a custom calibration, subclass `Calibration` and implement `calibrate()`.
See {doc}`../extending`.

## CalibrationTable: calibration result container

The uniform calibration result carrier.

**Construction**

`CalibrationTable(name, qubit_name=None, kind=None, inputs=None, outputs=None, fit_params=None, metadata=None)`

**Fields**

| Field | Type | Meaning |
|---|---|---|
| `name` | str | Calibration item name |
| `qubit_name` | str | Associated qubit name |
| `kind` | str | Calibration type (`"f_phi"`, `"f01"`, `"transfer_function"`, `"predistortion"`, etc.) |
| `inputs` | `np.ndarray` | Independent variable (e.g. flux $\Phi$) |
| `outputs` | `np.ndarray` | Dependent variable (e.g. frequency $\omega$) |
| `fit_params` | dict | Fit/iteration details |
| `metadata` | dict | Additional metadata |

**Methods**

- `evaluate(x) -> np.ndarray`: interpolate `outputs` at query points `x` (e.g.
  get frequency at a given flux). Cubic spline via SciPy; auto-degrades to
  quadratic/linear for too few points, with extrapolation.
- `inverse(y) -> np.ndarray`: inverse interpolation, finding the `inputs` such
  that `outputs ≈ y` (requires `outputs` to be monotonic; internally restricts to
  the monotonic segment to build the inverse).

## FluxResponseCalibration: f(Φ) curve calibration

Scan a sequence of DC flux points, measure the frequency at each with a Ramsey
sequence, and build an $f(\Phi)$ lookup table.

**Construction**

`FluxResponseCalibration(qubit, method="ramsey", h_list=None, tau=None, t_rabi=None)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test | — |
| `method` | str | Measurement method | `"ramsey"` |
| `h_list` | `np.ndarray` | Flux scan points ($\Phi_0$) | `linspace(-0.03, 0.03, 51)` |
| `tau` | float | Free-precession time per point (ns) | — |
| `t_rabi` | `np.ndarray` | Rabi pulse time axis (ns) | `CONFIG.pulse.t_rabi` |

**Methods**

- `calibrate() -> CalibrationTable`: runs a Ramsey frequency measurement at each
  flux point, returns a `kind="f_phi"` table, `inputs=flux`, `outputs=angular
  frequency`.

**Output**

`CalibrationTable`, `kind="f_phi"`, `inputs` are the flux scan points, `outputs`
the corresponding angular frequencies.

```{note}
`method="transient"` (unknown transient signal → polynomial fit of
$\Delta\omega(\Phi)$) is a future feature that currently raises
`NotImplementedError` and is not part of the v1 public surface.
```

## FrequencyMeasurement: single-point f01 measurement

Measures $f_{01}$ at a single flux working point (read-only, no tuning).

**Construction**

`FrequencyMeasurement(qubit, flux=0.0, method="ramsey", tau_list=None, t_rabi=None, t_global=None, f_artificial=0.1)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test | — |
| `flux` | float | Measurement flux point ($\Phi_0$) | `0.0` (sweet spot) |
| `method` | str | Measurement method | `"ramsey"` (or `"transient"`) |
| `tau_list` | `np.ndarray` | Free-evolution time sweep (ns) | `CONFIG.pulse.tau_list` |
| `t_rabi` | `np.ndarray` | Rabi time axis (ns) | `CONFIG.pulse.t_rabi` |
| `t_global` | `np.ndarray` | Global time axis (ns) | `CONFIG.pulse.t_global` |
| `f_artificial` | float or None | Artificial detuning (GHz); `None`=dual-sweep mode | `0.1` |

**Methods**

- `measure(flux=None) -> float`: returns a signed angular frequency (rad·GHz).
  - `method="ramsey"` (default): Ramsey τ-sweep + FFT peak. Single-sweep mode
    (`f_artificial=0.1`, assumes $|\Delta|<0.1$ GHz) is fast; `f_artificial=None`
    uses double-sweep mode, robust for arbitrary $|\Delta|$ and signed, at 2× the
    cost.
  - `method="transient"` (advanced): τ=0 orthogonal Ramsey (R_y–R_x and
    R_y–R_{-x}) differential readout + control-kernel sensitivity
    $G=\int k_1\,\mathrm{d}t$ to invert $\Delta\omega$ directly. Cheaper than the
    τ-sweep but relies on the weak-signal linear approximation; suited to
    $|\Delta\omega|$ near zero. With `order>=3` a cubic Newton correction is added,
    the cubic coefficient $G_3$ chosen by `g3_source`: `"fit"` (odd-polynomial fit
    of $p_\mathrm{diff}(\Delta)$ with adaptive scan range) or `"kernel_full"`
    (full off-diagonal kernel triple integral $\iiint k_3$).
- `calibrate() -> CalibrationTable`: packages the single-point measurement into
  a `kind="f01"` table.

**Output**

`measure()` returns `float` (signed angular frequency); `calibrate()` returns
`CalibrationTable` (`kind="f01"`).

```{note}
The transient method is single-point frequency **measurement**, not the same as
the transient frequency **calibration** excluded by the §8 boundary
(`FluxResponseCalibration` with `method="transient"`, building the $f(\Phi)$
curve): the former is implemented and public, the latter is not. Transient
**waveform reconstruction**, meanwhile, is a v1 core feature; see
{doc}`reconstruction`.
```

## SinglePointFrequencyCalibration: single-point frequency tuning

Closed-loop feedback tuning of $f_q(V)$ to a target frequency $f_\mathrm{target}$
(Vepsalainen 2022).

**Construction**

`SinglePointFrequencyCalibration(qubit, method="closed_loop", f_target=None, V_a=None, V_b=None, epsilon_f=1e-4, max_iter=20, measure_method="ramsey", step_method="secant", bracket_tightening=True, drive_policy="sweet", ...)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `qubit` | `TransmonQubit` | Qubit under test | — |
| `method` | str | Tuning method | `"closed_loop"` |
| `f_target` | float | Target frequency (rad·GHz) | — |
| `V_a` | float | Flux-voltage bracket lower bound | — |
| `V_b` | float | Flux-voltage bracket upper bound | — |
| `epsilon_f` | float | Convergence tolerance (GHz·2π) | `1e-4` |
| `max_iter` | int | Maximum iterations | `20` |
| `measure_method` | str | Frequency measurement method | `"ramsey"` (or `"transient"`) |
| `step_method` | str | Root-search stepper | `"secant"` (or `"bisection"`/`"gradient"`) |
| `bracket_tightening` | bool | Regula-falsi bracket narrowing | `True` |
| `drive_policy` | str or callable | Drive-frequency policy | `"sweet"` (or `"target"`/`"track"`/callable) |
| `sensitivity_source` | str | Sensitivity source | `"secant"` (or `"model"`) |
| `linear_range` | float | Linear window $\Delta_\mathrm{lin}$ | — |
| `rho` | float | Out-of-range threshold ratio | `0.6` |
| `converge_streak` | int | Consecutive convergence count | `1` |
| `omega_d_seed` | float | First track drive frequency (rad·GHz) | `None` |
| `damping` | float | Gradient damping factor (`step_method="gradient"` only) | `0.8` |
| `V_seed` | float | Gradient initial flux (`step_method="gradient"` only) | — |

**Methods**

- `calibrate() -> CalibrationTable`: executes the closed-loop tuning, returns a
  `kind="f01"` table; `fit_params["history"]` holds the full iteration trace.

**Output**

`CalibrationTable`, `kind="f01"`; `fit_params["history"]` contains per-iteration
$V_k$, $f_{q,k}$, $e_k$, $f_{d,k}$, $\delta_k$, and the `out_of_range` flag.

### Two updates per iteration: flux voltage and drive frequency

Each iteration of the loop runs **two independent updates**, acting on two
different actuators:

- **Flux-voltage update** (the root search, what actually moves $f_q$ toward the
  target): the gradient step is a normalised Newton step
  $V_{k+1} = V_k - \alpha\,e_k/\hat s_k$ (with $\alpha=$ `damping`,
  $\hat s_k=\partial f_q/\partial V$); the secant / bisection steps use a
  model-free secant / bisection instead. Under the local linear model the error
  contracts as $e_{k+1} \approx (1-\alpha)e_k$, converging to the **only** fixed
  point $f_q = f_\mathrm{target}$. This is the **only** update that can move the
  convergence target.
- **Drive-frequency update** (the observer, keeping the measurement trustworthy):
  a `drive_policy` sets the drive frequency $f_d$ each iteration so the
  *current measurement* stays in its linear window. It **never enters the error
  definition**, so it cannot move the convergence target; it only keeps $e_k$
  trustworthy.

**The invariant that makes this safe:** `measure()` returns the *absolute*
frequency $f_d + \hat\delta$, and the loop computes
$e_k = \mathrm{measure}(V_k) - f_\mathrm{target}$. So the control error is always
relative to the fixed absolute target, regardless of where $f_d$ sits. Driving
$f_d$ to follow the qubit does **not** hide residual error.

`drive_policy` (gradient step only; other steppers raise) takes a string or a
`callable(state) -> omega_d`:

| policy | $f_d$ | phase | use |
|---|---|---|---|
| `"sweet"` (default) | `qubit.frequency` | **COARSE_ACQUIRE** | back-compat; fine when the target is near the sweet spot, or for wide-range coarse acquisition |
| `"target"` | `f_target` | **LOCKED** | measured detuning *is* the control error; simplest terminal form |
| `"track"` | $\hat f_{q,k} + \hat s_k\,(V_{k+1}-V_k)$ | **TRACKING** | predicts $f_q$ at the next point so $|\delta|$ stays inside the linear window as the search moves |
| `callable` | user-defined | any | arbitrary feedback / filtering / prediction laws |

### The four-stage state machine and FrequencyCalibrationWorkflow

Combining the two updates above with different measurement protocols and
`drive_policy` values, the whole closed-loop calibration forms a four-stage
state machine. Orchestration is handled by `FrequencyCalibrationWorkflow`; each
stage is one `SinglePointFrequencyCalibration`, chained in order by the
workflow, seeding each stage from the previous stage's optimum flux and
frequency estimate:

```text
COARSE_ACQUIRE
    wide-range protocol (Ramsey) obtains the initial frequency
    drive_policy="sweet"; flux-voltage update only
        ↓  (residual detuning enters the linear window)
TRACKING
    drive frequency follows the predicted qubit frequency
    local protocol (transient) measures the frequency at high precision
    gradient descent updates the flux voltage
    drive_policy="track"
        ↓  (error falls near the target tolerance)
LOCKED
    drive frequency fixed at the target frequency
    local protocol measures the target error directly
    flux-voltage update only (de-bounced convergence)
    drive_policy="target"
        ↓  (detuning out of range |δ|>ρ·Δ_lin, or measurement fails)
REACQUIRE
    re-acquire: return to the wide-range stage to re-estimate the frequency
```

- **COARSE_ACQUIRE → TRACKING → LOCKED** is the normal forward path, driven in
  turn by the `"sweet"` / `"track"` / `"target"` drive policies; each stage hands
  off once its own `epsilon_f` is met.
- **REACQUIRE** is not an automatic transition but an **orchestration-layer
  response**: when `linear_range` ($\Delta_\mathrm{lin}$) with `rho`
  (default 0.6) flags $|\delta| > \rho\,\Delta_\mathrm{lin}$ in the
  `out_of_range` signal, or a measurement fails, the caller watches that signal
  (via `stop_predicate`) and inserts a fresh wide-range stage to re-estimate the
  frequency.

#### FrequencyCalibrationWorkflow

**Construction**

`FrequencyCalibrationWorkflow(stages, seed_drive_from_prev=True)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `stages` | `list[CalibrationStage]` | Stages to execute in order | — |
| `seed_drive_from_prev` | bool | Hand each stage the previous stage's frequency estimate as `omega_d_seed` | `True` |

**Methods**

- `run() -> CalibrationTable`: executes all stages in order, returns the final
  stage's `CalibrationTable`; `fit_params` contains the full state-machine trace.

#### CalibrationStage

**Construction**

`CalibrationStage(name, calibration, stop_predicate=None)`

**Fields**

| Field | Type | Meaning |
|---|---|---|
| `name` | str | Stage name (`"COARSE_ACQUIRE"`/`"TRACKING"`/`"LOCKED"`/`"REACQUIRE"`) |
| `calibration` | `SinglePointFrequencyCalibration` | The calibration instance for this stage (with `drive_policy`, `measure_method`, etc.) |
| `stop_predicate` | callable | Stage termination condition (default: stop when own `epsilon_f` is met) |

## WaveformCalibration: waveform / control-line calibration

The unified entry for the predistortion mainline.

**Construction**

`WaveformCalibration(method="transfer_function", distortion=None, measurement_protocol=None, qubit=None, control_line=None, fit_type="single_exp", transfer_model=None, predistortion_method="auto", n_taps=72, regularization=1e-6)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `method` | str | Calibration path | `"transfer_function"` (or `"predistortion"`) |
| `distortion` | `DistortionModel` | Distortion model (analytical path: take step response directly) | `None` |
| `measurement_protocol` | str | Quantum-simulation measurement protocol | `None` (or `"cryoscope"`/`"delay_ramsey"`/`"transient"`/`"pi_pulse"`) |
| `qubit` | `TransmonQubit` | Qubit (measurement path only) | `None` |
| `control_line` | `ControlLine` | Control line (measurement path only) | `None` |
| `fit_type` | str | Transfer-function fit type | `"single_exp"` (or `"multi_exp"`/`"fir"`/`"iir"`) |
| `transfer_model` | `DistortionModel` | Measured transfer function (predistortion path only) | `None` |
| `predistortion_method` | str | Predistortion design method | `"auto"` |
| `n_taps` | int | FIR filter order | `72` |
| `regularization` | float | Ridge regularization parameter | `1e-6` |

**Methods**

- `calibrate() -> CalibrationTable`: returns a `kind="transfer_function"` or
  `"predistortion"` table depending on `method`.
- `to_distortion_model() -> DistortionModel`: converts the calibration result
  straight into a `DistortionModel`.

**Output**

`CalibrationTable`, `kind` depends on `method`; `to_distortion_model()` returns
a `DistortionModel` ready to inject into a control line.

## PredistortionDesigner: inverse-filter designer

A standalone predistortion-filter designer, usable on its own or called by
`WaveformCalibration`.

**Construction**

`PredistortionDesigner(method="auto", n_taps=72, regularization=1e-6)`

**Fields**

| Field | Type | Meaning | Default |
|---|---|---|---|
| `method` | str | Design method | `"auto"` (or `"fir_inverse"`/`"iir_inverse"`/`"frequency_inverse"`) |
| `n_taps` | int | FIR filter order | `72` |
| `regularization` | float | Ridge regularization parameter | `1e-6` |

**Methods**

- `design(transfer_model, dt) -> DistortionModel`: given a transfer-function
  model, return its inverse model. `"auto"` uses an analytical IIR inverse for
  exponential models, frequency-domain inversion otherwise.
- `predistort(target, transfer=..., or inverse_model=...) -> Waveform`: apply
  predistortion directly to a target waveform.
- `check_pole_stability(b, a) -> bool`: check whether all IIR poles lie inside
  the unit circle (stable).

**Output**

`design()` returns `DistortionModel` (inverse model); `predistort()` returns the
predistorted `Waveform`.

## CalibrationScheduler: calibration task scheduler

The calibration task scheduler: maintains a registry of named calibration tasks
with dependencies, executed in order.

**Construction**

`CalibrationScheduler()`

**Methods**

- `register(name, cal_class, depends_on=...)`: register a task.
- `register_defaults()`: load the standard registry (e.g.
  `frequency_closed_loop` depends on `flux_response_ramsey`,
  `waveform_predistortion` depends on `waveform_transfer_function`).
- `run(name, **kw)`: run by name (auto-injecting dependency results).
- `run_next()`: run the next ready task.
- `get_result(name)` / `status()`: query.

**Output**

`run()` returns `CalibrationTable`; `status()` returns
`{"ready": [...], "running": ..., "completed": [...], "failed": [...]}`.

```{note}
The Kelly 2018 DAG automation (`check_state → maintain → auto_calibrate`) is
currently an interface stub, reserved for future implementation. Scheduling today
is manual.
```

## Minimal example

```python
import numpy as np
from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.calibration import (
    CalibrationTable,
    WaveformCalibration,
    PredistortionDesigner,
    CalibrationScheduler,
)

# 1) CalibrationTable interpolation lookup
table = CalibrationTable(
    name="flux_response",
    kind="f_phi",
    inputs=np.linspace(-0.03, 0.03, 7),
    outputs=np.array([5.5, 5.7, 5.9, 6.0, 5.9, 5.7, 5.5]) * 2 * np.pi,
)
f_at_zero = table.evaluate(np.array([0.0]))           # interpolate
phi_for_target = table.inverse(np.array([5.8 * 2 * np.pi]))  # inverse interpolation

# 2) Measure transfer function → design predistortion
distortion = SingleExponentialDistortion(amplitude=0.04, tau=200.0)
wf_cal = WaveformCalibration(distortion=distortion, fit_type="single_exp")
tf_result = wf_cal.calibrate()                         # kind="transfer_function"
inv_model = wf_cal.to_distortion_model()               # convert to DistortionModel

# 3) Design an inverse filter standalone
designer = PredistortionDesigner(method="auto")
inverse = designer.design(distortion, dt=1.0)          # dt=1 ns

# 4) Run the default registry via the scheduler
scheduler = CalibrationScheduler()
scheduler.register_defaults()
print(scheduler.status())   # {"ready": ["flux_response_ramsey", ...], ...}
```

## Physical role / extension

- This layer corresponds to the calibration routines in a real control room:
  after fabrication the first step measures $f(\Phi)$, then the operating
  frequency is tuned, and finally the control-line transfer function is measured
  and a compensation filter burned in. The three steps map respectively to
  `FluxResponseCalibration` → `SinglePointFrequencyCalibration` →
  `WaveformCalibration`.
- `FrequencyCalibrationWorkflow` orchestrates `SinglePointFrequencyCalibration`
  as a four-stage state machine (COARSE_ACQUIRE → TRACKING → LOCKED →
  REACQUIRE), automatically switching between wide-range Ramsey and
  high-precision transient measurement and converging to kHz level.
- `CalibrationTable`'s `evaluate`/`inverse` underpin the
  `inversion="calibration"` path of the {doc}`reconstruction` layer: the
  reconstructor is passed the calibration table and, during inversion, calls
  `inverse(φ)` to translate the measured phase back into flux amplitude.
- `CalibrationScheduler`'s dependency-injection mechanism automatically threads
  the output of `FluxResponseCalibration` (the flux range) into
  `SinglePointFrequencyCalibration`, sparing manual parameter passing.
- To add a custom calibration workflow, subclass the `Calibration` abstract
  base and implement `calibrate() -> CalibrationTable`, then register it with
  `CalibrationScheduler` and declare its dependencies. Full extension guide in
  {doc}`../extending`.
