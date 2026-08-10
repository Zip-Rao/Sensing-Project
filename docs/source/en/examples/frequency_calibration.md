# Frequency Calibration

## Overview

Frequency calibration is the platform's second product mainline: characterising
a superconducting transmon qubit's frequency-flux response $f_{01}(\Phi)$ to
answer two practical questions—**what is the qubit frequency at a given flux
bias**, and **how much flux bias is needed to tune the frequency to a target
value**. Unlike waveform reconstruction, frequency calibration does not invert
an external signal; it **characterises the device itself**: finding the working
point after fabrication, building a lookup table, and closed-loop tuning to a
target frequency when needed.

Physically, the transmon frequency is set by external flux modulating the
Josephson energy through the SQUID loop:

$$f_{01}(\Phi) \approx \frac{1}{2\pi}\left(\sqrt{8 E_J(\Phi)\, E_C} - E_C\right),
\qquad E_J(\Phi) = E_{J0}\,|\cos(\pi\Phi/\Phi_0)|.$$

Since $E_J \propto |\cos(\pi\Phi/\Phi_0)|$, $f_{01}(\Phi)$ is an **even function**
about $\Phi=0$, peaking at integer flux quanta (the **sweet spot**,
$\mathrm{d}f/\mathrm{d}\Phi=0$, first-order insensitive to flux noise). Biasing
to one side of the sweet spot yields sensing sensitivity
($\kappa = \mathrm{d}\omega/\mathrm{d}\Phi \neq 0$), which is precisely where
the `flux_bias` in {doc}`waveform_reconstruction` comes from.

The frequency calibration pipeline executes in the following order: single-point
measurement → sweep flux to build $f(\Phi)$ curve → lookup → closed-loop tuning
→ orchestration & run. The platform provides **two orchestration styles**: the
original one-shot pipeline (`FrequencyCalibrationWorkflow`) and the new
event-driven state machine (`FrequencyCalibrationRuntime` +
`FrequencyStateMachine`). The sections below follow this logical sequence
and conclude with a comparison to guide your choice.

## Pipeline Architecture

### Step 1: Single-Point Measurement — Calibration Layer (`FrequencyMeasurement`)

The basic operation underlying all frequency calibration is precisely measuring
$f_{01}$ at **a single flux working point**.
{py:class}`~sqc.calibration.frequency.FrequencyMeasurement` is a read-only,
no-tuning single-point frequency meter. Internally it drives a Ramsey sequence
(or transient orthogonal Ramsey), runs `mesolve`, and extracts the frequency
from the data.

The core method is `measure(flux=None, omega_d=None) -> float`: runs the
measurement at the specified flux bias `flux` and returns a signed angular
frequency (rad·GHz). The `omega_d` parameter specifies the reference drive
frequency—the measured detuning $\hat\delta = \hat f_q - f_d$ is added back to
give the absolute frequency. This is critical in subsequent closed-loop search:
feeding the previous estimate as `omega_d` keeps each measurement's detuning
inside the discriminator's linear window.

`method` selects one of two measurement protocols:

| `method` | Principle | Cost | Use when |
|---|---|---|---|
| `"ramsey"` | $\tau$-sweep + FFT peak | tens of `mesolve` calls | Robust, general-purpose |
| `"transient"` | $\tau=0$ orthogonal Ramsey + kernel sensitivity $G=\int k_1\,dt$ | 2 `mesolve` calls | Fast; suited to $\Delta\omega \approx 0$ |

Ramsey mode defaults to single-sweep (`f_artificial=0.1` GHz), assuming
$|\Delta| < 0.1$ GHz; set `f_artificial=None` for double-sweep, which is
robust for arbitrary detuning and returns a sign at 2× the cost.

### Step 2: Sweep Flux to Build $f(\Phi)$ — Calibration Layer (`FluxResponseCalibration`)

With single-point measurement capability in hand, the next step is to scan
along the flux axis, calling `FrequencyMeasurement` at each point to build a
frequency-flux lookup table.
{py:class}`~sqc.calibration.frequency.FluxResponseCalibration` handles this
process.

`@dataclass` fields: `qubit`, `h_list` (flux scan points, $\Phi_0$, default 51
points), `method` (currently only `"ramsey"`). `calibrate()` returns a
`kind="f_phi"` {py:class}`~sqc.calibration.CalibrationTable`: `inputs` are the
flux points, `outputs` the angular frequencies. This table has two downstream
uses: providing flux bounds `[V_a, V_b]` for closed-loop tuning, and direct
forward/reverse lookup.

### Step 3: Lookup and Closed-Loop Tuning — Calibration Layer (`CalibrationTable` + `SinglePointFrequencyCalibration`)

#### Forward / Reverse Lookup

{py:class}`~sqc.calibration.CalibrationTable` provides two cubic-spline-based
lookup methods:

- `evaluate(x)` — **forward**: flux → frequency (e.g. "what is $f_{01}$ at
  $\Phi=0.015$").
- `inverse(y)` — **reverse**: frequency → flux (e.g. "how much bias is needed
  for this target frequency").

Because $f(\Phi)$ is even and globally non-monotonic, `inverse` automatically
restricts to the monotonic branch; target-frequency queries must stay on one
side of the sweet spot, or the solution is not unique.

#### Closed-Loop Tuning

A lookup gives an **open-loop estimate**—interpolation infers the required flux
but does not account for real measurement noise and model mismatch.
{py:class}`~sqc.calibration.frequency.SinglePointFrequencyCalibration` drives
$f_q(V)$ to a target $f_\text{target}$ via closed-loop feedback
(Vepsalainen 2022), iteratively measuring and adjusting to converge on the
true value.

Core fields: `f_target` (target angular frequency), `V_a`/`V_b` (flux bounds
from the monotonic branch identified by `FluxResponseCalibration`),
`step_method` (root-finding method), `measure_method` (per-iteration
measurement protocol, delegated internally to `FrequencyMeasurement`).

Three stepping methods:

| `step_method` | Convergence | Needs bracket | Notes |
|---|---|---|---|
| `"secant"` | Superlinear, 1–3 iters | Yes | `bracket_tightening` (regula falsi) auto-shrinks |
| `"bisection"` | $O(\log_2)$, 10–15 iters | Yes | Bracket width halves each step; auto-handles even $f(\Phi)$ |
| `"gradient"` | Damped Newton | No | Only needs `V_seed`; `damping` (default 0.8) suppresses overshoot; `best_V` tracks optimum |

`calibrate()` returns a `kind="f01"` table; `fit_params["history"]` holds the
full iteration trace (per-iteration $V$, $f$, residual), suitable for plotting
convergence.

```{note}
Each closed-loop iteration runs two independent updates: a **flux-voltage
update** (the root search, moving $f_q$ toward the target) and a **drive-
frequency update** (the observer, `drive_policy` setting $f_d$ to keep the
measurement in its linear window). The drive frequency never enters the error
definition—$e_k = \text{measure}(V_k) - f_\text{target}$ is always relative to
the fixed target—so it only affects measurement trustworthiness, not the
convergence target. Three drive policies (`"sweet"`/`"target"`/`"track"`),
combined with different measurement protocols, form a six-state event-driven
protocol (Acquire → Track → Verify → Lock + Reacquire) orchestrated by
`FrequencyCalibrationRuntime`. See {doc}`../building_blocks/calibration` for
the full description. The old four-stage pipeline interface
(`FrequencyCalibrationWorkflow`) remains available; internally it already
delegates to the shared `DampedSecantTracker` controller.
```

### Step 4: Orchestration & Run — Workflow Layer

The previous three steps all operate at a **single flux point** or with a
**fixed strategy**. Real calibration needs to switch measurement protocols and
drive policies depending on the search phase—the coarse phase needs wide-range
Ramsey, the fine phase benefits from fast transient measurement; normal tracking
needs the drive frequency to follow the qubit, while verification requires a
frozen bias and an independent pass/fail judgement. The platform provides two
ways to orchestrate this.

#### Old API: one-shot pipeline `FrequencyCalibrationWorkflow`

{py:class}`~sqc.workflows.frequency_calibration.FrequencyCalibrationWorkflow`
chains multiple `SinglePointFrequencyCalibration` instances into an ordered
pipeline, seeding each stage from the previous stage's optimum. Because each
stage's `measure_method` is fixed at construction time, the pipeline is
**one-way, irreversible**—stages can only move forward, with no backtracking
or retry logic.

Two construction modes:

- **Default hybrid preset** (omit `stages`): auto-builds a two-stage
  transient→Ramsey pipeline. `switch_residual` (default $2\pi \cdot 5$ MHz)
  controls the coarse-to-fine handoff threshold.
- **Explicit pipeline** (pass `stages=[...]`): each stage is a
  `CalibrationStage` independently specifying `measure_method`, `step_method`,
  `epsilon_f`, `drive_policy`, etc.

`run()` returns a merged iteration history (each row tagged with
`phase`/`global_iter`/`cost`), `V_final`, `residual`, `converged`, and more.

#### New API: event-driven state machine `FrequencyStateMachine` + `FrequencyCalibrationRuntime`

The old pipeline's fundamental limitation is that it only has one path—**"keep
going forward"**. If the coarse search overshoots, the transient discriminator
loses lock, or verification fails, the pipeline cannot go back to re-acquire;
it continues to the next stage and produces meaningless results. For real
experiments that may run for hours—where the locked state must continuously
monitor drift, schedule periodic audits, and automatically recover from
loss-of-lock—a one-shot pipeline is entirely insufficient.

The V2 API solves these problems with an **event-driven six-state protocol**.
The calibration process is decomposed into discrete states, each with a clear
responsibility and explicit entry/exit guard conditions:

```text
Acquire ──→ Track ──→ Verify ──→ Lock
              ↑          │          │
              │          │          │
            Reacquire ◄─────────────┘
```

**Six states, each with one job**:

| State | Responsibility | Can change flux? |
|---|---|---|
| **Acquire** | Wide-range Ramsey acquisition of absolute frequency—gets the initial "seed" estimate | No |
| **Track** | Local transient frequency measurement + secant sensitivity estimate + flux stepping—the **only** state allowed to change bias | **Yes** |
| **Verify** | Freeze the candidate bias; run independent double-sweep Ramsey to judge whether final tolerance is met | No (frozen on entry) |
| **Lock** | Long-term stabilisation: low-cost drift monitoring + scheduled Ramsey audits | No |
| **Reacquire** | After loss-of-lock / out-of-range / low confidence, re-run wide-range acquisition | No |
| **SafeStop** | Budget exhausted / interlock / unrecoverable fault → safe hold | safe bias only |

**Key design rules**:

1. **Only Track may change the flux bias**—other states cannot even "nudge" it.
   This prevents the self-deception of adjusting bias while supposedly verifying.
2. **Verify freezes the candidate bias at entry**, holding it until exit. The
   actual applied bias must match the frozen value within `bias_freeze_tolerance`.
3. **Lock drift detection never adjusts bias directly**—small drift goes
   Lock→Verify for independent confirmation; large jumps go Lock→Reacquire for
   a fresh acquisition. This guarantees that "locked" means the bias was never
   secretly moved.
4. **Analytic $f(\Phi)$ is a simulation oracle only**—it never feeds the
   transition reducer. All decisions are based on actual measurement results.

**Threshold hierarchy**—four levels control transition tightness:

```
epsilon_hold  <  epsilon_final  <  epsilon_enter  <  Delta_val
(10 kHz)         (100 kHz)         (5 MHz)           (20 MHz)
physics bandwidth Verify pass       candidate entry   local validity window
```

- **Candidate condition** (Track → Verify): $|\text{error}| + \text{uncertainty} \le \epsilon_\text{enter}$
- **Verification condition** (Verify → Lock): $|\text{error}| + \text{uncertainty} \le \epsilon_\text{final}$, for $N_\text{verify}$ consecutive passes
- **Local validity** (must we leave Track?): $|\text{error}| + \text{uncertainty} \le \Delta_\text{val}$

**Lock monitor hysteresis**—a critical design point that prevents noise-induced
state chatter. The Lock state uses a low-cost transient monitor whose noise
characteristics differ from the independent Verify Ramsey, so it cannot simply
reuse the same threshold:

| Monitor result | Action |
|---|---|
| $U_\text{mon} \le \epsilon_\text{mon\_clear}$ | All clear, reset suspect counter |
| $\epsilon_\text{mon\_clear} < U_\text{mon} \le \epsilon_\text{mon\_suspect}$ | Grey zone: accumulate suspect count; after $N_\text{mon\_suspect}$ → Verify |
| $\epsilon_\text{mon\_suspect} < U_\text{mon} < \Delta_\text{mon\_reacquire}$ | Clear drift, immediate Lock → Verify |
| $U_\text{mon} \ge \Delta_\text{mon\_reacquire}$ or reference lost | Large jump, Lock → Reacquire directly |

The benefit: a point that just passed Verify won't exit Lock on a single noisy
monitor reading; genuine slow drift is caught after accumulating a grey-zone
streak; catastrophic jumps trigger immediate re-acquisition without wasting
verification attempts.

**Command–event architecture**—the state machine does no I/O directly; it is
decoupled from the measurement backend through commands and events:

```text
┌──────────────┐     command      ┌────────────┐     execute     ┌─────────────┐
│ StateMachine │ ───────────────→ │  Runtime   │ ──────────────→ │ SQCExecutor │
│  (pure logic) │                  │ (orchestr.)│                 │ (QuTiP backend)│
│              │ ←─────────────── │            │ ←────────────── │             │
└──────────────┘     event        └────────────┘    result       └─────────────┘
```

- **Commands** (machine → outside): `AcquireFrequency`, `TrackFrequency(bias, drive)`,
  `VerifyFrequency(frozen_bias, drive)`, `MonitorFrequency(locked_bias)`, `SafeHold(bias)`
- **Events** (outside → machine): `MeasurementSucceeded` (carries frequency,
  uncertainty, validity flags), `MeasurementTechnicalFailure`, `MeasurementRejected`,
  `TimerElapsed`, `InterlockTriggered`, `BudgetExhausted`, `CancelRequested`

Every event carries a `command_id` matching the pending command, making the
system idempotent—calling `next_command()` repeatedly returns the same command
until a matching event arrives and the state advances.

**Persistence & checkpoint recovery**: `FrequencyCalibrationRuntime` provides
`save_run(dir)` and `load_run(dir, qubit)` methods that write run state to five
files:

| File | Content |
|---|---|
| `config.json` | Full protocol configuration |
| `commands.jsonl` | Every command and its resulting event |
| `transitions.jsonl` | Every state transition (from / to / reason) |
| `checkpoint.json` | Full state machine snapshot (resumable from here) |
| `result.json` | Final result summary |

This is essential for long-running experiments—after a power loss or abnormal
exit, the run can resume from the last checkpoint without losing calibration
progress.

**Relationship to the old API**: both share the same `DampedSecantTracker`
control law (damped-secant step formula); numerical behaviour is identical.
The old `FrequencyCalibrationWorkflow` already delegates each gradient stage to
the tracker internally. The difference is at the orchestration layer—fixed
pipeline vs. event loop + transition table.

### Underpinning: Device Layer (`TransmonQubit`)

The physical foundation of all calibration operations is
{py:class}`~sqc.devices.transmon.TransmonQubit`—holding $E_C$, $E_J$,
decoherence times, and the current flux bias. Calibration classes internally
set the DC flux bias via `qubit.qubit_in_mag(FluxSignal)`, update the
Hamiltonian, and run `mesolve`. `qubit.frequency` provides the sweet-spot
frequency as the default drive reference for Ramsey measurements.

## Usage

### End-to-End Pipeline

```python
import numpy as np
from sqc.devices.transmon import TransmonQubit
from sqc.calibration import FluxResponseCalibration, FrequencyMeasurement
from sqc.calibration import SinglePointFrequencyCalibration

# ── Device: EC/EJ passed as angular frequencies (rad·GHz) ────────────
qubit = TransmonQubit(
    EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
    T1=10_000, T2=8_000, flux=0.0, n_levels=3,
)

# ── 1. Single-point measurement: f01 at the sweet spot ───────────────
fm = FrequencyMeasurement(qubit=qubit, flux=0.0, method="ramsey")
f01 = fm.measure()                       # signed angular frequency (rad·GHz)
print(f"sweet-spot f01 = {f01 / (2*np.pi):.4f} GHz")

# ── 2. Sweep flux to build f(Φ) lookup table ─────────────────────────
cal = FluxResponseCalibration(
    qubit=qubit,
    method="ramsey",
    h_list=np.linspace(-0.03, 0.03, 5),  # coarse grid for the demo
)
table = cal.calibrate()                  # CalibrationTable, kind="f_phi"

# ── 3a. Forward lookup: Φ → f ────────────────────────────────────────
f_at_bias = table.evaluate(np.array([0.015]))

# ── 3b. Reverse lookup + closed-loop tuning ──────────────────────────
f_target = table.outputs.max() * 0.999   # just below the peak (monotonic branch)
tuner = SinglePointFrequencyCalibration(
    qubit=qubit,
    f_target=f_target,
    V_a=0.0, V_b=0.03,                   # bracket bounds from the monotonic branch
    step_method="secant",                 # secant method, typically 1–3 iters
)
result = tuner.calibrate()               # CalibrationTable, kind="f01"
print("tuned bias:", result.fit_params["V_opt"],
      "converged:", result.fit_params["converged"])
```

### Multi-Stage Hybrid Pipeline

```python
from sqc.workflows.frequency_calibration import FrequencyCalibrationWorkflow

# Default two-stage hybrid: transient coarse (cubic correction) → Ramsey fine
wf = FrequencyCalibrationWorkflow(
    qubit=qubit,
    f_target=f_target,
    V_a=0.0, V_b=0.03,
    switch_residual=2 * np.pi * 5e-3,    # 5 MHz coarse-to-fine handoff
    epsilon_f=1e-4,                       # final convergence tolerance
)
hybrid_result = wf.run()
print(f"V_final={hybrid_result['V_final']:.6f}, "
      f"residual={hybrid_result['residual']/(2*np.pi)*1e3:.2f} MHz, "
      f"converged={hybrid_result['converged']}")
```

### Event-Driven State Machine (V2 API)

The `FrequencyCalibrationWorkflow` above is good for "set parameters, run, read
the result." In a real experiment, calibration may run for hours and encounter
bias drift, measurement failures, or interlock triggers. The V2 API handles
these with a **pausable, resumable, recoverable** six-state protocol.

Using it is straightforward: construct a `FrequencyCalibrationConfig` and a
`FrequencyCalibrationRuntime`, then call `run()`. The runtime internally:

1. Creates and starts the `FrequencyStateMachine`
2. Loops: `next_command()` → `SQCExecutor` performs the measurement →
   `handle(event)` advances the state
3. In Track state, automatically manages the `DampedSecantTracker` lifecycle
   (`initialize` → `propose` → `accept`)
4. Checkpoints every 10 commands

```python
from sqc.workflows.frequency_runtime import FrequencyCalibrationRuntime
from sqc.workflows.frequency_state_machine import FrequencyCalibrationConfig

# ── Protocol config: all thresholds and budgets in one place ────────────
config = FrequencyCalibrationConfig(
    # -- thresholds (all angular frequency, rad·GHz) --
    epsilon_enter=2 * np.pi * 20e-3,      # 20 MHz  — condition to enter Verify
    epsilon_final=2 * np.pi * 2e-3,       # 2 MHz   — condition to pass Verify
    Delta_val=2 * np.pi * 50e-3,          # 50 MHz  — local validity window

    # -- verification policy --
    N_verify=2,                            # 2 consecutive passes to confirm

    # -- Lock monitor hysteresis --
    epsilon_mon_clear=2 * np.pi * 5e-3,   # 5 MHz   — monitor says "clean"
    epsilon_mon_suspect=2 * np.pi * 10e-3,# 10 MHz  — monitor says "suspect"
    Delta_mon_reacquire=2 * np.pi * 30e-3,# 30 MHz  — large jump → reacquire
    N_mon_suspect=3,                       # 3 consecutive suspects → Verify

    # -- budgets --
    max_commands=200,                      # total command limit
    max_reacquire_attempts=5,             # reacquisition limit
    max_wall_time=3600.0,                 # 1 hour wall-clock

    # -- Track controller parameters --
    damping=0.8,                           # damping factor (< 1 suppresses overshoot)
    first_bias_step=0.01,                  # first probe step (Φ₀)
    max_bias_step=0.02,                    # max bias change per step (Φ₀)
)

# ── Run ────────────────────────────────────────────────────────────────
runtime = FrequencyCalibrationRuntime(qubit=qubit, f_target=f_target, config=config)
result = runtime.run()

# ── Results ─────────────────────────────────────────────────────────────
print(f"Final state: {result['state']}")            # "lock" / "safe_stop"
print(f"Run status:  {result['run_status']}")       # "calibrated" / "completed" / "failed"
print(f"Final freq:  {result['f_final']/(2*np.pi):.6f} GHz")
print(f"Candidate Φ: {result['candidate_bias']:.6f} Φ₀")
print(f"Commands: {result['n_commands']}, wall time: {result['elapsed']:.1f}s")

# ── Transition trace ────────────────────────────────────────────────────
for t in result['transition_log']:
    print(f"  v{t['version']}: {t['from']} → {t['to']}  [{t['reason']}]")

# ── Persistence ─────────────────────────────────────────────────────────
runtime.save_run("calibration_run_001")
# Writes config.json / commands.jsonl / transitions.jsonl / checkpoint.json / result.json

# ... hours later, resume from checkpoint ...
# runtime2 = FrequencyCalibrationRuntime.load_run("calibration_run_001", qubit=qubit)
# runtime2.run()  # machine resumes from the last checkpoint
```

This shows the core V2 workflow: configure thresholds → run → inspect the
trace → persist. Each entry in `transition_log` records the source state,
destination state, and the reason code that triggered the transition (e.g.
`TARGET_CANDIDATE`, `DRIFT_SUSPECTED`, `REACQUIRE_LIMIT`), enabling post-hoc
audit—why did we reacquire? How many consecutive suspect readings triggered
the Verify? Every decision is traceable.

```{note}
**Old API vs new API**: `FrequencyCalibrationWorkflow` is a one-shot pipeline
(stage1→stage2→...) suited to rapid prototyping and batch simulation.
`FrequencyCalibrationRuntime` + `FrequencyStateMachine` is the event-driven
architecture: supports state recovery (Track→Reacquire→Track), long-term Lock
monitoring with scheduled audits, and persistence with checkpoint recovery—
ideal for long-running calibration experiments and multi-scenario reliability
testing. Both share the same `DampedSecantTracker` control law (damped-secant
step formula); numerical behaviour is identical.
See {doc}`../building_blocks/workflows` for details.
```

## Reading the Results

- `fm.measure()` returns the signed angular frequency (rad·GHz) at a single
  working point; divide by $2\pi$ for GHz. Default single-sweep assumes
  $|\Delta| < 0.1$ GHz; if the point may be far from the sweet spot, set
  `f_artificial=None` for double-sweep mode.
- `table` is a `kind="f_phi"` `CalibrationTable`: `inputs` are the flux points,
  `outputs` the angular frequencies. `evaluate` does cubic-spline interpolation
  (flux→frequency); `inverse` does the reverse (frequency→flux), automatically
  restricting to the monotonic branch. A typical curve peaks at $\Phi=0$ and
  falls symmetrically on both sides—the sweet spot gives first-order flux-noise
  immunity; biasing to the side yields sensing sensitivity.
- The closed-loop `result.fit_params` holds `V_opt` (optimal flux), `converged`
  (whether within `epsilon_f` tolerance), and `history` (per-iteration
  $V$/$f$/residual). `history` can be used directly to plot convergence.
- The multi-stage hybrid `run()` returns a merged `history` (each row tagged
  with `phase`), plus `stage_boundaries` (global iteration index at each
  stage's end) and `metrics` (cumulative `mesolve` call cost).
- `FrequencyMeasurement` supports `order=3` cubic Newton correction, with
  `g3_source` selecting the cubic-coefficient source: `"fit"` (odd-polynomial
  fit of $p_\text{diff}(\Delta)$, adaptive range) or `"kernel_full"`
  (off-diagonal kernel $\iiint k_3\,dt^3$).
- **V2 state machine results**: `runtime.run()` returns `state` (final protocol
  state, e.g. `"lock"`), `run_status` (run lifecycle, e.g. `"calibrated"`
  means Lock was successfully entered), `f_final` (final frequency estimate),
  `candidate_bias` (locked bias value), `transition_log` (complete transition
  trace—each entry has `from`/`to`/`reason`/`version`), `n_commands` (commands
  consumed), and `elapsed` (wall-clock time). `transition_log` is the key to
  post-hoc audit—every state change and its trigger reason is recorded.
- See {doc}`../building_blocks/calibration` for the full field definitions and
  method options of each calibration class; see
  {doc}`../building_blocks/workflows` for V2 state machine thresholds and
  configuration.
