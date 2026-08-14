# Claim-to-evidence record

This record distinguishes method definitions, deterministic numerical evidence,
and experimental evidence. Analytic simulation ground truth is not an
independent measurement.

| Claim | Evidence | Status |
|---|---|---|
| Short-pulse sequences encode locally invertible frequency information | Orthogonal zero-delay branches and F2 response curve | Deterministic numerical case study complete |
| The cubic local inverse is consistent across two coefficient routes | F2 odd-polynomial fit and full-kernel calculation at fixed drive | Deterministic numerical cross-check complete |
| The sampled cubic branch turns near 17.6 MHz | F2 frozen data and response-fold calculation | Deterministic numerical evidence complete |
| A drive-tracked short-pulse estimator can participate in local feedback | F4 seeded Track trajectory | Deterministic numerical case study complete |
| Drive-tracked full-kernel Track is cheaper than Ramsey in counted solver calls | F4 instrumented solver calls at a common analytic residual threshold | Deterministic numerical evidence complete for this benchmark |
| The reported final residual is not estimator self-consistency | Analytic transmon frequency--flux calculation at the returned bias | Deterministic ground-truth check complete; not a Verify measurement |
| The complete Acquire--Track--Verify--Lock state machine is executable | State-machine, runtime, executor, interruption, and checkpoint tests | Software capability implemented; hardware protocol not validated |
| Verify is statistically independent of Track | Separate Verify command and estimator contract | Protocol requirement implemented; no finite-shot experimental data |
| The simulation FSM reaches independently measured Verify-to-Lock completion under finite shots | Fig. 3 paired ensemble: 12/12 short-pulse and 12/12 Ramsey runs, with role-separated RNG streams | Independent simulation evidence complete; not experimental verification |
| Drive tracking is required for this short-pulse closed-loop case | S5 tracked-drive run locks; fixed-drive ablation triggers the command/actuation interlock | Deterministic simulation ablation complete for this case |
| The simulated loop can recover under the frozen drift-plus-jump stress test | Fig. 3/S8: 7/8 short-pulse and 6/8 Ramsey runs reach Verify-to-Lock | Simulation stress evidence only; not a general robustness or hardware-stability claim |
| Runtime checkpoint/resume preserves the calibration result | Fig. 3 package: zero frequency/bias difference and equal 3620-call totals | Deterministic software-recovery evidence complete |
| Lock maintains long-term frequency stability | Hardware monitor and periodic Ramsey time series | Not established |
| The protocol suppresses frequency-noise PSD or improves Allan stability | Long-duration hardware data and analysis | Not established |
| Experimental calibration is faster | Circuit repetitions, acquisition time, readout, drift, and uncertainty data | Not established |
| The protocol is robust to model, pulse, and readout error | Parameter/noise ensembles and hardware tests | Not established |
| Short-pulse sequences are generally useful across calibration tasks | Additional observables, sequences, devices, and hardware tests | Not established |

## Data roles

- **Response-calibration data** fit response coefficients and nuisance terms.
- **Response-validation data** are held out to determine the validated local interval.
- **Closed-loop independent-verification data** test a frozen candidate bias and cannot be Track convergence hits.
- **Long-term monitoring data** are collected only after verified entry to Lock.
- **Analytic simulation ground truth** scores deterministic trajectories and must not be called independent experimental verification.

## Wording boundary

Use "solver-cost reduction in deterministic numerical simulations" for the
current result. Describe the broader result as one demonstrated local
frequency-calibration case. Do not write "experimental speedup", "verified
long-term lock", "bandwidth improvement", "PSD suppression", or "Allan
stability improvement" until the corresponding hardware evidence exists.
