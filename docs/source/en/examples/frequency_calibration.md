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
→ multi-stage orchestration. The sections below follow this logical sequence.

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
combined with different measurement protocols, form a four-stage state machine
(COARSE_ACQUIRE → TRACKING → LOCKED → REACQUIRE) orchestrated by
`FrequencyCalibrationWorkflow` in the next step. See
{doc}`../building_blocks/calibration` for the full description of each policy.
```

### Step 4: Multi-Stage Orchestration — Workflow Layer (`FrequencyCalibrationWorkflow`)

A single `SinglePointFrequencyCalibration` has its `measure_method` fixed at
construction time—it cannot switch protocols mid-search, yet the coarse phase
needs the wide-range Ramsey while the fine phase benefits from the fast
transient method.
{py:class}`~sqc.workflows.frequency_calibration.FrequencyCalibrationWorkflow`
chains multiple calibration stages into an ordered pipeline, seeding each stage
from the previous stage's optimum flux and frequency estimate.

Two construction modes:

- **Default hybrid preset** (omit `stages`): auto-builds a two-stage
  transient→Ramsey pipeline. `switch_residual` (default $2\pi \cdot 5$ MHz)
  controls the coarse-to-fine handoff threshold.
- **Explicit pipeline** (pass `stages=[...]`): each stage is a
  `CalibrationStage` independently specifying `measure_method`, `step_method`,
  `epsilon_f`, `drive_policy`, etc.

`run()` returns a merged iteration history (each row tagged with
`phase`/`global_iter`/`cost`), `V_final`, `residual`, `converged`, and more.

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

```{note}
`FluxResponseCalibration` runs a full Ramsey $\tau$-sweep (tens of `mesolve`
calls) at every flux point, and the default `h_list` has 51 points, so a full
calibration is a minutes-scale task. The 5-point coarse grid above is only to
show the flow; for real runs, densify according to your accuracy needs or scan
finely only near the band of interest.
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
- See {doc}`../building_blocks/calibration` for the full field definitions and
  method options of each calibration class.
