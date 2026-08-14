# Fig. 3 simulation evidence

This directory is a reproducible, simulation-only evidence package for Fig. 3
and Supplemental Figs. S5, S6, and S8. It is independent of the historical F
figure archives and does not use their figures or NPZ data as numerical input.
It is not experimental or hardware evidence.

## Frozen protocol

- two-level transmon with `EC/h=0.2 GHz` and `EJ/h=10 GHz`;
- `T1=100 us`, `T2=50 us`, and `dt=0.5 ns`;
- 10 ns sin-squared pi/2 pulse;
- simulation candidate interval `Delta_val/(2*pi)=14 MHz`;
- analytic dispersion withheld from Acquire, Track, and Verify and used only
  for simulation scoring;
- successful completion requires an actual `Verify -> Lock` transition.

`protocol.json` embeds the upstream Fig. 2 protocol and records SHA-256 hashes
for its source config and manifest. The repository revision and dirty-worktree
flag are also recorded.

## Simulation matrix

`run-index.csv` indexes 50 freshly generated runs:

- six deterministic runs: two strategies at three target detunings;
- 24 paired finite-shot runs: 12 seeds for each strategy;
- 16 paired drift-plus-jump runs: eight seeds for each strategy;
- two deliberate solver-budget exhaustion runs;
- one drive-tracked and one fixed-drive ablation run.

Finite-shot roles are configured before execution. Acquire uses 1024 shots per
circuit. Short-pulse Track uses 1024; Ramsey Track uses 6. Verify and Monitor
use 6. A double-sweep Verify therefore uses 4800 physical shots, and the two
required consecutive Verify passes use at most 9600 shots in one episode,
below the frozen `10000`-shot limit. The assignment-error channel uses
`p(0->1)=p(1->0)=0.01`.

The process ensemble adds `0.02 MHz` drift per command and a controlled
`0.8 MHz` jump at command index 5. These parameters define a numerical stress
test; they are not inferred hardware-noise parameters.

## Files

- `runs/<run-id>/result.json`: terminal state, completion event, provenance,
  protocol config, and cost totals.
- `runs/<run-id>/measurements.jsonl`: command/event journal.
- `runs/<run-id>/state-history.jsonl`: FSM transition history.
- `runs/<run-id>/cost-ledger.jsonl`: per-command circuits, physical shots,
  `mesolve`, `sesolve`, and total solver calls.
- `aggregate.{csv,json,npz}`: state-wise descriptive summaries.
- `checkpoint-resume.json`: continuous-versus-resumed equivalence check.
- `derived/`: figure-specific CSV and NPZ projections generated only from the
  files above.
- `manifest.json`: SHA-256 for every package file except the manifest itself,
  plus hashes for the rendered figure exports.

The verified checkpoint/resume comparison has zero final-frequency difference,
zero bias difference, and identical solver-call totals. The complete package
contains no JSON `NaN` or infinity values.

## Reproduction

From the repository root, regenerate all raw runs and the manifest with:

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
C:\Users\21034\anaconda3\envs\qutip-env\python.exe paper\simulations\fig3\build_evidence.py
```

Then regenerate the derived data and all figure formats with:

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
C:\Users\21034\anaconda3\envs\qutip-env\python.exe paper\figures\fig3\generate_fig3_evidence.py
```

Rendered figures are stored under `paper/figures/fig3/` as PDF, SVG, and
600-dpi PNG files. Solver calls are a computational accounting unit, not
acquisition time or experimental speed.
