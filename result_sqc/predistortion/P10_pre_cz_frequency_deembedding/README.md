# P10 pre-CZ frequency-deembedding benchmark

P10 corrects the measurement-chain confounding identified in P9. It estimates
the protocol-plus-reconstruction response `R(f)` from an independent ideal-line
calibration split, freezes the transient Wiener parameter and Cryoscope SG
window on a second D0 split, and removes `R(f)` before estimating the unknown
control-line transfer function `H(f)`.

The inverse FIR has no fallback to failed frequency bins. If fewer than eight
continuous non-zero bins pass the qualified-band test, filter construction
fails explicitly. Filter length and ridge are selected on a held-out validation
waveform subject to AWG peak, slew-rate, and noise-gain constraints. The final
short flat-top is used once as the test waveform.

The 262,144-shots-per-arm resource is selected using D0 only. A separate D4
main-line resource diagnostic repeats the frozen protocol parameters at
262,144, 524,288, and 1,048,576 shots per arm. It reports empirical continuous
design-band availability only and is excluded from resource selection,
protocol-parameter selection, and FIR selection.

This directory stops at G0--G3. It does not run or plot CZ simulations.
