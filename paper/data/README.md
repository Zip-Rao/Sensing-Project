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
