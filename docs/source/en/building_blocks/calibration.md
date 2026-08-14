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
| `DampedSecantTracker` | Shared stepwise secant controller: `initialize`→`propose`→`accept` |
| `FrequencyEstimate` / `TrackSnapshot` / `TrackProposal` / `TrackStepResult` | Controller data structures |
| `FrequencyCalibrationWorkflow` | Multi-stage closed-loop calibration orchestration (legacy staged workflow) |
| `CalibrationStage` | Per-stage parameter container (used by `FrequencyCalibrationWorkflow`) |

**Event-driven frequency calibration state machine (V2)**

| Class | Role |
|---|---|
| `FrequencyStateMachine` | Six-state event-driven protocol: Acquire→Track→Verify→Lock + Reacquire |
| `FrequencyCalibrationConfig` | Protocol thresholds / budgets / policy parameters |
| `FrequencyCalibrationRuntime` | Event-loop orchestration + tracker management + `save_run()`/`load_run()` persistence |
| `SQCExecutor` | Command→QuTiP measurement: Ramsey (Acquire/Reacquire/Verify), Transient (Track/Monitor) |
| `FaultInjectionExecutor` | Configurable fault injection for testing recovery paths |

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

### Stepwise secant controller: `DampedSecantTracker`

`sqc/calibration/frequency_control.py` extracts a **pure-math** controller from
`_closed_loop_gradient()` — no QuTiP dependency, no I/O; only bias step
computation, convergence checks, and best-point tracking. It is shared by
the legacy batch API and the new `FrequencyStateMachine` Track state.

**Three-step interface**:

```python
from sqc.calibration.frequency_control import DampedSecantTracker, FrequencyEstimate

tracker = DampedSecantTracker(f_target=..., damping=0.8, max_bias_step=0.02)
snapshot = tracker.initialize(seed_estimate, V_seed)    # 1. initialize
proposal = tracker.propose(snapshot)                     # 2. propose next bias
result = tracker.accept(snapshot, estimate, proposal)    # 3. accept, update state
```

**Step formula** (inside `propose`): first step uses a fixed probe
`first_bias_step·sign(e)`; subsequent steps use the damped secant
`step = damping·e·(dV/de)`, clamped to `[-max_bias_step, max_bias_step]`.
When converged, `step=0` (re-measure in place for confirmation with
`converge_streak>1`).

### Frequency calibration state machine V2: `FrequencyStateMachine`

v2.19 introduces an **event-driven** six-state protocol, replacing the old
four-stage pipeline design:

```text
Acquire → Track → Verify → Lock
             ↑        │        │
             │        │        │
           Reacquire <─────────+
```

**Six-state responsibilities (hard constraints)**:

| State | Responsibility | Can change flux? | Measurement role |
|---|---|---|---|
| Acquire | Wide-range Ramsey acquisition of absolute frequency seed | No | global acquisition |
| Track | Short-pulse local frequency measurement + secant sensitivity + flux update + drive tracking | **Yes** | local loop estimator |
| Verify | Freeze candidate bias; independent Ramsey assessment of final tolerance | **No** | independent verifier |
| Lock | Long-term frequency stabilization; low-cost drift monitoring; scheduled independent Ramsey audits | **No** | monitor only |
| Reacquire | Wide-range re-acquisition after reference / branch / local validity loss | No | global recovery |
| SafeStop | Safe hold after cancel / budget / interlock / unrecoverable fault | safe bias only | no science measurement |

**Threshold hierarchy**:

```
epsilon_hold < epsilon_final < epsilon_enter < Delta_val
   10 kHz         100 kHz         5 MHz        20 MHz
  (physics goal) (Verify pass)  (candidate)  (local window)
```

The first three thresholds act on target residual
$r=\omega_{01}-\omega_{\mathrm{tar}}$, whereas `Delta_val` acts on probe
detuning $\Delta=\omega_{01}-\omega_d$. Track therefore checks
$|\widehat\Delta|+z\sigma_\Delta\leq\Delta_\mathrm{val}$ and must not
substitute target residual.

**Command–event architecture**:

```python
from sqc.workflows.frequency_runtime import FrequencyCalibrationRuntime
from sqc.workflows.frequency_state_machine import FrequencyCalibrationConfig

config = FrequencyCalibrationConfig(
    epsilon_enter=2 * np.pi * 20e-3,    # 20 MHz
    epsilon_final=2 * np.pi * 2e-3,     # 2 MHz
    confidence_multiplier=1.0,           # preselect z_(1-beta) for experiments
    N_verify=2,
    max_commands=30,
    stop_after_lock_cycles=3,
)
runtime = FrequencyCalibrationRuntime(qubit=q, f_target=f_target, config=config)
result = runtime.run()
# bounded run: result["state"] → "safe_stop"
# result["run_status"] → "completed", result["safe_hold_confirmed"] → True
```

`CALIBRATED` marks the first successful Verify → Lock transition; it is not a
terminal runtime status. Lock monitoring and scheduled Ramsey audits continue
until a command, shot, solver-call, wall-time, or `stop_after_lock_cycles`
budget ends the run. Science-command costs are reserved before execution.
SafeStop then issues `SafeHold` with bounded retries, and
`safe_hold_confirmed` reports whether `SafeHoldApplied` was received.

**Cooperative interruption**: UI/API threads call
`runtime.request_cancel(reason, source)`, which only sets a thread-safe
`CancellationToken`. The runtime converts it to
`CancelRequested → SafeStop → SafeHold` at command boundaries and during the
interruptible monitor wait. `KeyboardInterrupt` follows the same path by
default; set `handle_keyboard_interrupt=False` to propagate it. With
`checkpoint_directory` configured, the interrupt and final SafeHold state are
persisted best-effort. A synchronous `executor.execute()` cannot be forcibly
pre-empted: cancellation received during a measurement takes effect when that
call returns, before another science command starts.

**Track validity guards** run in a fixed order: backend `out_of_range`, the
optional detuning window `linear_range - guard_margin`, optional experimental
`min_confidence`, secant sensitivity bounds `[S_min, S_max]`, and finally
`U_loop <= Delta_val`. The deterministic simulation backend reports
`uncertainty_source="deterministic_zero"`; zero there is a declared modelling
assumption, not experimental confidence. Verify is separately bounded by
`max_verify_attempts_per_episode` and `max_verify_shots`; Lock cadence is set by
`monitor_interval` and `audit_interval`.

**Lock monitor hysteresis** (prevents noise-induced state chatter):

| Condition | Transition |
|---|---|
| `U_mon ≤ epsilon_mon_clear` | Lock → Lock, clear suspect streak |
| `epsilon_mon_clear < U_mon ≤ epsilon_mon_suspect` | accumulate suspect streak; at `N_mon_suspect` → Verify |
| `epsilon_mon_suspect < U_mon < Delta_mon_reacquire` | immediate Lock → Verify |
| `U_mon ≥ Delta_mon_reacquire` or reference lost | Lock → Reacquire |

**Persistence**: `runtime.save_run(dir)` writes `config.json` / `commands.jsonl` /
`transitions.jsonl` / `checkpoint.json` / `result.json`;
`FrequencyCalibrationRuntime.load_run(dir, qubit, executor=...)` restores the
configuration, budget, tracker, Verify/Lock counters, retry state, and pending
command, then `run()` resumes directly. A pending command is replayed with its
original `command_id`, so executors must be idempotent by command ID. Recovery
therefore has explicit **at-least-once** semantics; it does not claim strict
hardware exactly-once execution. Each journal row records state before/after,
bias, drive frequency, measured frequency, residual, uncertainty,
valid/ambiguity flags, shots, elapsed time, diagnostics, and transition reason.

### FrequencyCalibrationWorkflow (legacy, kept for compatibility)

The old multi-stage staged workflow interface remains available and internally
delegates to `DampedSecantTracker`:

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
- `DampedSecantTracker` is the single implementation of the frequency calibration
  control law — shared by the legacy `_closed_loop_gradient` batch API and the new
  `FrequencyStateMachine` stepwise interface, guaranteeing identical control
  behavior.
- `FrequencyStateMachine` is the V2 event-driven protocol: six states
  (Acquire→Track→Verify→Lock + Reacquire + SafeStop) with full guards, budgets,
  and monitor hysteresis. `FrequencyCalibrationRuntime` orchestrates the
  command→event loop, manages the `DampedSecantTracker` lifecycle, and supports
  `save_run`/`load_run` persistence with checkpoint recovery.
- `FrequencyCalibrationWorkflow` (legacy) orchestrates
  `SinglePointFrequencyCalibration` as a multi-stage pipeline (e.g.
  transient→Ramsey hybrid); internally each gradient stage delegates to
  `DampedSecantTracker`.
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
