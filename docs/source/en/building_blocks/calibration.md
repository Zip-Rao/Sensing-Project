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

Grouped by function into four sets:

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

The uniform calibration result carrier (`@dataclass`). Core fields: `name`,
`qubit_name`, `kind` (e.g. `"f_phi"`, `"f01"`, `"transfer_function"`,
`"predistortion"`), `inputs`/`outputs` (independent/dependent variable, e.g.
flux↔frequency), `fit_params` (fit/iteration details), `metadata`.

Two lookup methods (cubic spline via SciPy; auto-degrades to quadratic/linear
for too few points, with extrapolation):

- `evaluate(x) -> np.ndarray`: interpolate `outputs` at query points `x` (e.g.
  get frequency at a given flux).
- `inverse(y) -> np.ndarray`: inverse interpolation, finding the `inputs` such
  that `outputs ≈ y` (requires `outputs` to be monotonic; internally restricts to
  the monotonic segment to build the inverse).

## FluxResponseCalibration: f(Φ) curve calibration

Scan a sequence of DC flux points, measure the frequency at each with a Ramsey
sequence, and build an $f(\Phi)$ lookup table. `@dataclass`, fields: `qubit`,
`method` (`"ramsey"`), `h_list` (flux scan points, default `linspace(-0.03, 0.03, 51)`),
`tau` (free-precession time per point), `t_rabi` (Rabi pulse time axis).
`calibrate()` returns a table with `kind="f_phi"`, `inputs=flux`,
`outputs=angular frequency`.

```{note}
`method="transient"` (unknown transient signal → polynomial fit of
$\Delta\omega(\Phi)$) is a future feature that currently raises
`NotImplementedError` and is not part of the v1 public surface.
```

## FrequencyMeasurement: single-point f01 measurement

Measures $f_{01}$ at a single flux working point (read-only, no tuning).
`@dataclass`, fields: `qubit`, `flux` (measurement flux point, default 0 i.e. the
sweet spot), `method`, `tau_list`, `t_rabi`, `t_global`, `f_artificial`. The core
method `measure(flux=None) -> float` returns a signed angular frequency;
`calibrate()` packages the single-point measurement into a `kind="f01"` table.

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
(Vepsalainen 2022). `@dataclass`, fields: `qubit`, `method` (`"closed_loop"`),
`f_target`, `V_a`/`V_b` (bracketing bounds, from a preceding
`FluxResponseCalibration`), `epsilon_f` (convergence tolerance, default
$10^{-4}$ GHz·2π), `max_iter`, `measure_method` (`"ramsey"` or `"transient"`),
`bracket_tightening` (regula-falsi bracket narrowing, default on), and
`step_method` (one of three):

- `"secant"` (default): secant method, superlinear convergence, typically 1–3
  iterations.
- `"bisection"`: bisection, $O(\log_2)$ convergence, bracket width halves each
  step, suited to visualisation. Auto-handles an even $f(\Phi)$ crossing the
  sweet spot (auto-splits the bracket).
- `"gradient"`: damped secant (numerical-gradient Newton step), does not
  require pre-bracketing $V_a$/$V_b$, only `V_seed`; a damping factor `damping`
  (default 0.8) suppresses overshoot; `best_V` tracks the historical best point.

`calibrate()` returns a `kind="f01"` table; `fit_params["history"]` holds the
full iteration trace.

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

Why the second update is needed: each frequency measurement is a *local*
discriminator, reporting a trustworthy detuning
$\hat\delta = \hat f_q - f_d$ only while $|\delta| < \Delta_\mathrm{lin}$
(especially the cheap transient method). If $f_d$ is pinned at the sweet spot
while the search probes flux far away, the measured error becomes biased and the
loop can diverge; the drive-frequency update removes that bias.

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

Supporting fields: `sensitivity_source` (`"secant"` — model-free slope of the
last two points, default; or `"model"` — analytic `qubit.sensitivity(V)`) feeds
the `"track"` prediction; `linear_range` ($\Delta_\mathrm{lin}$) with `rho`
(default 0.6) flags `out_of_range` in each history row and in the
`stop_predicate` state when $|\delta| > \rho\,\Delta_\mathrm{lin}$ — a caller
can watch that signal to trigger a wider **re-acquire** stage. Each `history`
row now also carries `omega_d` and `delta` for diagnostics.

`converge_streak` (default 1) requires $|e|\le\epsilon$ on that many
*consecutive* iterations before declaring convergence — it re-measures in place
to de-bounce a noisy hit, useful in the LOCKED phase. `omega_d_seed` supplies the
first `"track"` drive $f_{d,0}=\hat f_{q,0}$; without it a fresh track stage
falls back to the sweet spot on iteration 0, which is fatal when the target is
far (the transient then sees the full offset).

### The four-stage state machine

Combining the two updates above with different measurement protocols and
`drive_policy` values, the whole closed-loop calibration forms a four-stage
state machine. Each stage is one `SinglePointFrequencyCalibration` pass, chained
in order by `FrequencyCalibrationWorkflow` (below), seeding each stage from the
previous stage's optimum flux and frequency estimate:

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
  frequency. Each `history` row carries `omega_d` and `delta` for diagnostics.

Fields supporting this state machine: `sensitivity_source` (`"secant"` —
model-free slope of the last two points, default; or `"model"` — analytic
`qubit.sensitivity(V)`) feeds the TRACKING prediction; `converge_streak`
(default 1) requires $|e|\le\epsilon$ on that many *consecutive* iterations
before declaring convergence, for noise rejection in the LOCKED stage;
`omega_d_seed` supplies the first `"track"` drive $f_{d,0}=\hat f_{q,0}$ —
without it a fresh TRACKING stage falls back to the sweet spot on iteration 0,
fatal when the target is far (the transient then sees the full offset).

`FrequencyCalibrationWorkflow` (below) is the orchestration that runs this state
machine: its `CalibrationStage` exposes `drive_policy`, `sensitivity_source`,
`linear_range`/`rho`, `converge_streak`, and `seed_drive_from_prev` (default on:
hands each TRACKING stage the previous stage's frequency estimate as
`omega_d_seed`). See notebook §B2.2c for a −20 MHz far-target run where a naive
single-stage transient stalls at the wrong flux while the state machine locks to
kHz.

## WaveformCalibration: waveform / control-line calibration

The unified entry for the predistortion mainline (`@dataclass`); `method` selects
one of two paths:

- `"transfer_function"`: measure the control-line step response, fit it to a
  `DistortionModel`. Fields: `distortion` (analytical path: take the model's step
  response directly), or `measurement_protocol`
  (`"cryoscope"`/`"delay_ramsey"`/`"transient"`/`"pi_pulse"`, with
  `qubit`+`control_line`, via a true quantum-simulation measurement); `fit_type`
  (`"multi_exp"`/`"single_exp"`/`"fir"`/`"iir"`). `to_distortion_model()`
  converts the result straight into a `DistortionModel`.
- `"predistortion"`: given a measured transfer function `transfer_model`, design
  a compensation filter (internally delegating to `PredistortionDesigner`);
  fields `predistortion_method`, `n_taps`, `regularization`.

`calibrate()` returns a `kind="transfer_function"` or `"predistortion"` table
depending on `method`.

## PredistortionDesigner: inverse-filter designer

A standalone predistortion-filter designer (`@dataclass`), usable on its own or
called by `WaveformCalibration`. Fields: `method`
(`"auto"`/`"fir_inverse"`/`"iir_inverse"`/`"frequency_inverse"`), `n_taps`,
`regularization`. Core methods:

- `design(transfer_model, dt) -> DistortionModel`: given a transfer-function
  model, return its inverse model. `"auto"` uses an analytical IIR inverse for
  exponential models, frequency-domain inversion otherwise.
- `predistort(target, transfer=..., or inverse_model=...) -> Waveform`: apply
  predistortion directly to a target waveform.
- `check_pole_stability(b, a) -> bool`: check whether all IIR poles lie inside
  the unit circle (stable).

## CalibrationScheduler: calibration task scheduler

The calibration task scheduler (`@dataclass`): maintains a registry of named
calibration tasks with dependencies, executed in order. `register(name,
cal_class, depends_on=...)` registers a task; `register_defaults()` loads the
standard registry (e.g. `frequency_closed_loop` depends on
`flux_response_ramsey`, `waveform_predistortion` depends on
`waveform_transfer_function`). Execution: `run(name, **kw)` runs by name (auto-
injecting dependency results), `run_next()` runs the next ready task;
`get_result(name)`, `status()` query.

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
