# Waveform Reconstruction

## Goal

Given a transient magnetic flux $\Phi(t)$ acting on the qubit whose **shape is
unknown**, recover that waveform using nothing but measurements of the qubit's
excited-state population. This is the platform's first product mainline: the
qubit acts as a sensor, the external flux is the signal under test, and a single
{py:class}`~sqc.workflows.SensingWorkflow` entry point runs the whole
"configure → measure → reconstruct" pipeline.

## Physics

External flux $\Phi(t)$ modulates the Josephson energy through the SQUID loop
and thereby shifts the qubit frequency $\omega_T(\Phi)$. A Ramsey-type control
pulse is **slid** across the flux signal: at each delay the qubit accumulates a
coherent phase proportional to how the nearby flux pushes the frequency, so the
excited-state population $p_e$ traces a *blurred* image of the flux signal as a
function of the sliding delay. The blur is described by the **control kernel**
$k(t)$ — the response of the pulse to instantaneous flux (see the flux-sensing
section of {doc}`../theory`).

The measurement and the waveform are therefore related by a convolution:

$$\Delta p(t) \;\approx\; (k * \Phi)(t),$$

where $\Delta p$ is the difference between two sliding measurements — one "with
signal" and one "zero-flux reference" — with the flux-independent baseline
subtracted off. **Reconstruction is the inverse problem** of this convolution:
given $\Delta p$ and the kernel $k$, solve for $\Phi(t)$. The
{doc}`../building_blocks/reconstruction` layer offers several solvers — linear
Wiener deconvolution, Hammerstein-Wiener (which accounts for the $\omega(\Phi)$
nonlinearity), and full density-matrix Levenberg-Marquardt inversion — and this
example wires them together and runs an A/B comparison via `SensingWorkflow`.

## End-to-end code

```python
from sqc.workflows import SensingWorkflow

# ── 1. Configure a single sensing experiment ─────────────────────────
wf = SensingWorkflow()
wf.configure(
    protocol="transient",     # transient field sensing (sliding Ramsey)
    flux_bias=0.1,            # bias away from the sweet spot, finite κ
    signal_type=4,            # waveform to recover: asymmetric impulse
    signal_amplitude=0.02,    # peak flux (Φ₀)
    signal_center=100,        # impulse centre (ns)
    reconstruction="wiener",  # start with linear Wiener deconvolution
    lambda_reg=5.0,           # Wiener regularisation strength
)

# ── 2. Run the full pipeline: measure + reconstruct ──────────────────
result = wf.run()                         # WorkflowResult
rec = result.reconstructed_signal         # FluxSignal: the recovered Φ(t)
print("recovered points:", len(rec.signal))

# ── 3. Compare three reconstruction algorithms on the same data ──────
cmp = wf.compare(methods=["wiener", "hammerstein", "lm"])
print("best method:", cmp.best)           # auto-selected by RMSE
for m in cmp.methods:
    print(f"  {m}: rmse={cmp.metrics[m]['rmse']:.3e}  snr={cmp.metrics[m]['snr']:.2f}")

# ── 4. Auto-plot: top = Δp, bottom = recovered waveform vs truth ─────
wf.plot()
```

```{note}
`compare()` reuses the measurement cached by `run()` and only swaps the
reconstruction algorithm — it does not re-run `mesolve`, so all three methods
are compared on the **same** $\Delta p$. `best` is chosen by RMSE; because this
is a simulation the ground-truth waveform `flux_samples` is known, which is what
makes RMSE meaningful.
```

## Reading the results

- `result.measurement` holds four things: `p_e` (excited-state population),
  `delta_p` (baseline-subtracted signal), `kernel` (the control kernel $k$), and
  `flux_samples` (the ground-truth waveform injected in simulation, available
  only in simulation). Axis information lives in `result.measurement.axes`:
  `scan` (sliding delay) and `t_flux` (waveform time axis).
- `result.reconstructed_signal` is the recovered
  {py:class}`~sqc.control.FluxSignal` — i.e. $\Phi(t)$ inferred from population
  data alone. `wf.plot()` overlays it on the ground truth so the reconstruction
  quality is immediate.
- The `metrics` from `compare()` give each method an `rmse` / `snr` / `peak`. A
  typical takeaway: Wiener is fastest and adequate for small signals;
  Hammerstein adds the weak $\omega(\Phi)$ nonlinearity and lowers RMSE at large
  amplitudes; LM is the most accurate but the slowest (full density-matrix
  forward simulation plus iteration). The `best` field picks out the
  lowest-RMSE method for you.
- To probe the sensitivity limit, feed this same pipeline into
  `wf.sweep("signal.amplitude", [...])` to sweep amplitude and obtain the
  SNR-vs-amplitude curve — see {doc}`../building_blocks/workflows` for the full
  API, and {doc}`../building_blocks/reconstruction` for each solver's formula and
  regime of validity.
