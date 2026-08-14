# Fig. 1 numerical placeholder evidence

This directory holds the numerical evidence used by Fig. 1 placeholders. The
editable SVG sources, assembly metadata, scripts, and rendered outputs live in
`paper/figures/fig1/`.

`fig1c-v1-data.npz` is retained as simulation evidence even though that panel
is no longer assembled into the active Fig. 1. Its blue curve is a direct
master-equation calculation and its interleaved markers are simulated held-out
samples. It is not experimental evidence.

The Fig. 1(a) generator currently reads the frozen F1/F4 result archives in
`result_sqc/`. Those upstream archives remain the authoritative inputs; this
directory does not duplicate them.
