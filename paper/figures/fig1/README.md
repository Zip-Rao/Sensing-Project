# Fig. 1 figure assets

This directory contains the active Fig. 1 editable sources, assembly metadata,
generation scripts, and PDF/PNG/SVG exports. Numerical evidence is kept in
`paper/data/fig1-simulation/`; superseded designs remain in
`paper/figures/archive/fig1/superseded/`.

The manuscript-facing stable export is `fig1-system-protocol.pdf`. Rebuild it
with:

```powershell
C:\Users\21034\anaconda3\envs\qutip-env\python.exe paper/figures/fig1/fig1-system-protocol-assemble.py
```

`fig1c-v1-*` is retained as a numerical placeholder/evidence visualization but
is not part of the active assembled figure.

The active probe schematic uses a sine-squared 10-ns pulse window, matching the
frozen Fig. 2 simulation protocol. The retained `fig1c-v1-*` placeholder uses
its older square-pulse protocol and must not be used as evidence for Fig. 2.
