# Figure data and provenance

Place frozen, publication-specific data under this directory.  Each figure
should have:

- an immutable raw or baseline data file;
- a processed data file when transformation is nontrivial;
- the exact plotting script and command;
- parameter/configuration metadata;
- a short record of the source commit and generation date.

Do not treat images copied from `../report/` as final evidence.  Regenerate
paper figures from frozen data with English labels and column-aware dimensions.

Current publication-specific stores:

- `fig1-simulation/`: numerical placeholder evidence used by Fig. 1 sources.
- `fig2-simulation/`: the single source of truth for the simulation-only Fig. 2.
- `fig3-simulation/`: independent FSM, finite-shot, process-stress, cost, and
  checkpoint evidence for simulation-only Fig. 3 and Figs. S5/S6/S8.

Rendered PDF/PNG/SVG files never belong under `paper/data/`; evidence NPZ/CSV
files never belong under numbered directories in `paper/figures/`.
