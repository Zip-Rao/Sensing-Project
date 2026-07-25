# Frequency Calibration

## Goal

Measure the qubit's frequency response to flux, $f_{01}(\Phi)$, and use it to
answer two practical questions: **what is the qubit frequency at a given flux
bias**, and **how much flux bias is needed to tune the frequency to a target
value**. This is the platform's second product mainline: rather than
reconstructing an external signal, it **characterises the device itself** —
finding the working point after fabrication, building a lookup table, and
closed-loop tuning to a target frequency when needed. This example uses the
{doc}`../building_blocks/calibration` layer's classes directly.

## Physics

The Transmon frequency is set by external flux modulating the Josephson energy
through the SQUID loop (see {doc}`../theory`):

$$f_{01}(\Phi) \approx \frac{1}{2\pi}\left(\sqrt{8 E_J(\Phi)\, E_C} - E_C\right),
\qquad E_J(\Phi) = E_{J0}\,\lvert\cos(\pi\Phi/\Phi_0)\rvert.$$

Because $E_J\propto\lvert\cos(\pi\Phi/\Phi_0)\rvert$, $f_{01}(\Phi)$ is an **even
function** about $\Phi=0$ and peaks at integer flux quanta — the **sweet spot**,
where $\mathrm{d}f/\mathrm{d}\Phi=0$ and the qubit is first-order insensitive to
flux noise.

Single-point frequency measurement uses a Ramsey sequence: a known **artificial
detuning** $f_a$ is applied so the qubit oscillates during free precession at
$\lvert f_a - \Delta\rvert$ (where $\Delta$ is the true detuning); an FFT peak of
$p_e(\tau)$ then recovers $f_{01}$. Scanning flux point by point and measuring
the frequency at each traces out the $f_{01}(\Phi)$ curve.

With the curve in hand, the **forward** direction is interpolation
($\Phi\to f$) and the **reverse** is inverse interpolation ($f\to\Phi$). Note
that the curve is even and globally non-monotonic, so the inverse is defined only
on a **single monotonic branch** — when looking up the flux for a target
frequency, stay on one side of the sweet spot.

## End-to-end code

```python
import numpy as np
from sqc.devices.transmon import TransmonQubit
from sqc.calibration import FluxResponseCalibration, FrequencyMeasurement

# ── device: EC/EJ passed as angular frequencies (rad·GHz) ────────────
qubit = TransmonQubit(
    EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
    T1=10_000, T2=8_000, flux=0.0, n_levels=3,
)

# ── 1. single-point measurement: f01 at the sweet spot ───────────────
fm = FrequencyMeasurement(qubit=qubit, flux=0.0, method="ramsey")
f01 = fm.measure()                       # signed angular frequency (rad·GHz)
print(f"sweet-spot f01 = {f01 / (2*np.pi):.4f} GHz")

# ── 2. sweep flux to build the f(Φ) lookup table ─────────────────────
#     A denser h_list is more accurate but costs one Ramsey sweep per
#     point, so runtime grows linearly (see the note below).
cal = FluxResponseCalibration(
    qubit=qubit,
    method="ramsey",
    h_list=np.linspace(-0.03, 0.03, 5),  # coarse grid for the demo
)
table = cal.calibrate()                  # CalibrationTable, kind="f_phi"

# ── 3. forward / reverse lookup ──────────────────────────────────────
f_at_bias = table.evaluate(np.array([0.015]))          # Φ → f
f_target = table.outputs.max() * 0.999                 # just below the peak (monotonic branch)
phi_needed = table.inverse(np.array([f_target]))       # f → Φ
print(f"target f={f_target/(2*np.pi):.4f} GHz needs bias Φ={phi_needed[0]:.4f}")
```

The inverse lookup is an **open-loop estimate**. For higher accuracy, take the
flux range from $f(\Phi)$ as a bracket and **closed-loop tune** to the target
frequency with `SinglePointFrequencyCalibration`:

```python
from sqc.calibration import SinglePointFrequencyCalibration

tuner = SinglePointFrequencyCalibration(
    qubit=qubit,
    f_target=f_target,
    V_a=0.0, V_b=0.03,          # bracket bounds from the monotonic branch of f(Φ)
    step_method="secant",       # secant method, typically converges in 1–3 iterations
)
result = tuner.calibrate()      # CalibrationTable, kind="f01"
print("tuned bias:", result.fit_params["V_opt"],
      "converged:", result.fit_params["converged"])
```

```{note}
`FluxResponseCalibration` runs a full Ramsey τ-sweep (tens of `mesolve` calls) at
every flux point, and the default `h_list` has 51 points — a full calibration is
a **minutes-scale** task. The 5-point coarse grid above is only to show the flow;
for real runs, densify according to your accuracy needs or scan finely only near
the band of interest.
```

## Reading the results

- `fm.measure()` returns the signed angular frequency ($\mathrm{rad\cdot GHz}$)
  at a single working point; divide by $2\pi$ for GHz. The default single-sweep
  mode assumes $\lvert\Delta\rvert<0.1$ GHz; if the point may be far from the
  sweet spot, set `f_artificial=None` for the double-sweep mode, which is robust
  for arbitrary detuning and returns a sign.
- `table` is a `kind="f_phi"`
  {py:class}`~sqc.calibration.CalibrationTable`: `inputs` are the flux points and
  `outputs` the corresponding angular frequencies. `evaluate` does cubic-spline
  interpolation flux→frequency, and `inverse` does the reverse. Because
  $f(\Phi)$ is even, `inverse` automatically uses the monotonic branch — a target
  frequency lookup must land on one side of the sweet spot, or the solution is
  not unique.
- A typical curve peaks at $\Phi=0$ (the sweet spot) and falls off symmetrically
  on both sides; the peak is the device's maximum operating frequency. Work at
  the sweet spot for first-order flux-noise immunity, or bias to the side for
  sensing sensitivity (finite $\kappa=\mathrm{d}\omega/\mathrm{d}\Phi$) — which is
  exactly where the `flux_bias` in {doc}`waveform_reconstruction` comes from.
- The closed-loop `result.fit_params` holds a full `history` (per-iteration $V$,
  $f$, and residual) you can use to plot convergence; `converged` flags whether
  it reached the `epsilon_f` tolerance.
- For each calibration class's fields, method options, and the scheduling
  mechanism, see {doc}`../building_blocks/calibration`; for how calibration tables
  feed reconstruction's table-lookup inversion, see
  {doc}`../building_blocks/reconstruction`.
