# Fig. 2 simulation evidence

This directory is the publication-specific evidence store for the current
simulation-only Fig. 2. It is intentionally separate from
`paper/figures/fig2/`, which contains plotting/generation code and rendered
figures only.

## Evidence products

- `config.json`: frozen protocol, device, grid, threshold, and convention data.
- `response.npz`: dissipative master-equation response and calibration fit.
- `estimator-errors.npz`: estimator outputs derived from `response.npz`, plus
  the ideal full-kernel cross-check.
- `validation-summary.npz`: held-out statistics derived from estimator errors.
- `solver-cost.npz`: Ramsey benchmark results and short-pulse call accounting.
- `solver-cost.csv`: human-readable projection of the cost comparison.
- `manifest.json`: file roles, dependency direction, and provenance.

The response calculations use dissipative `mesolve` with the configured T1 and
T2. The full-kernel coefficients use ideal `sesolve` internally through
`KernelEstimator`; this is an ideal-kernel versus dissipative-response
cross-check, not an identical-conditions comparison.

All data are deterministic (`sigma=0`), so a random seed is not applicable.
The green interval is a **simulation-defined candidate interval**, not an
experimentally validated operating range. Solver calls are an auditable
computational cost unit; they are not acquisition time or hardware speed.

Regenerate all evidence and figures from the repository root with:

```powershell
C:\Users\21034\anaconda3\envs\qutip-env\python.exe paper/figures/fig2/fig2-simulation-build-data.py
```

Re-render from frozen evidence without rerunning the simulations with each
panel script's `--plot-only` option.

## Panel and supplemental mapping

- Current simulation panel (a) / S2 reads `response.npz` and
  `estimator-errors.npz`. S2 shows the dissipative master-equation response,
  disjoint response-calibration and held-out grids, and the separate fitted
  and ideal-full-kernel cubic coefficient routes.
- Current simulation panel (b) / S3 reads `estimator-errors.npz`. The short-
  pulse estimators receive only the response. The analytic transmon dispersion
  is inverted afterward on the local `[0.5, 1]` branch to define a flux
  coordinate and scoring truth.
- Current simulation panel (c) / S7 reads `validation-summary.npz` and reports
  deterministic statistics for only the eight held-out points inside
  `|Delta|/(2*pi) <= 14 MHz`. It contains no bootstrap or finite-shot interval.
- Current simulation panel (d) / S6 reads `solver-cost.npz`; the accompanying
  `solver-cost.csv` exposes the cold-start, eight-estimate amortized, and
  marginal call counts.

The full response grid extends from -28 to +28 MHz as a valid rotating-frame
constant-detuning calculation. At the frozen reference flux
`0.9553166181245093 Phi0`, the selected analytic transmon branch does not have
a physical flux point for every positive detuning in that wider grid. Those
unmappable coordinates remain `NaN` in `estimator-errors.npz`; S3 therefore
shows the physically mappable local range and highlights the frozen +/-14 MHz
candidate interval. S2 retains the complete detuning-domain response.

The stable internal composite is `fig2-simulation-evidence.{pdf,png,svg}`.
Independent supplemental exports are:

- `figS2-response-simulation.{pdf,png,svg}`
- `figS3-frequency-recovery-simulation.{pdf,png,svg}`
- `figS6-solver-cost-simulation.{pdf,png,svg}`
- `figS7-validation-simulation.{pdf,png,svg}`

These exports are simulation evidence. Future main-text Fig. 2 requires real
response calibration and held-out data, an independent Ramsey/spectroscopy
reference, repeated finite-shot statistics, and hardware acquisition timing.
