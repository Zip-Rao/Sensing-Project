# Paper figures

This directory contains manuscript-ready exports and editable sources. Figure
panels remain provisional until they are assembled and checked at final print
size.

## Composite Fig. 1

`fig1-system-protocol-editable.svg` is the preferred assembled figure. Its page
is 180 mm by 106.3 mm for a REVTeX two-column `figure*`. Panels (a) and (c)
occupy the 90 mm left column; panels (b) and (d) occupy the right column at
their natural aspect ratios. The SVG keeps the background, four panels, and
uniform panel labels as separate Inkscape layers.

Use `fig1-system-protocol.pdf` for manuscript inclusion and
`fig1-system-protocol.png` as the reference preview. Rebuild the composite with
`fig1-system-protocol-assemble.py` after changing an individual panel.

## Fig. 1(a): numerical flux-calibration placeholder

The preferred global-plus-inset draft is `fig1a-v2-editable.svg`. Its page is
88 mm by 56 mm. The main axes show one complete lobe of the analytic transmon
dispersion, while the inset magnifies the actual calibration branch and shows
the visited points from the frozen drive-tracked loop. The highlighted box and
light connectors encode only geometric magnification, not a validated probe
range.

The earlier local-only v1 draft is retained under
`archive/fig1/superseded/` for comparison. The SVG files can be edited in
Inkscape 1.4 or later.

This is a numerical layout placeholder, not experimental data. The source
arrays remain in the existing F1 and F4 result archives. Rebuild both panels
with `fig1ac-v1-generate.py` from the repository root.

## Fig. 1(b): probe timing

Open `fig1b-v6-editable.svg` in Inkscape 1.4 or later. The page is 53 mm by
26.5 mm and is intended as a narrow panel in the composite Fig. 1 layout.

The SVG layers separate the background, circuit lines, pulse waveforms,
measurement symbols, and labels. It shows scanned Ramsey, fixed-delay I/Q
Ramsey, and the zero-delay phase-cycled short-pulse probe. Each pulse window
contains one waveform; phase alternatives are encoded by `R_phi` and `R_+-x`
labels and explained in the caption rather than drawn as overlapping carriers.
All measurement symbols share one aligned column. The displayed pulse window
uses a sin-squared envelope. Timing distances are schematic and not to scale.

`fig1b-v6-reference.pdf` and `fig1b-v6-reference.png` are reference exports.
Version v5 is retained under `archive/fig1/superseded/` for comparison.

## Fig. 1(c): numerical differential-response placeholder

Open `fig1c-v1-editable.svg` in Inkscape 1.4 or later. The page is 88 mm by
44 mm. The blue line is a direct master-equation calculation of
`p_d=(p_+X-p_-X)/2`; orange open circles are a separate interleaved detuning
grid retained as placeholder validation samples. The green interval is a
candidate local interval and must not be described as experimentally validated.

`fig1c-v1-data.npz` stores the plotted response arrays and protocol metadata.
The current placeholder follows the frozen numerical setup: 10 ns square
pi/2 pulses. If the sin-squared envelope drawn schematically in Fig. 1(b)
becomes the physical protocol, regenerate this response and all dependent
coefficients/results before finalizing the manuscript.

## Fig. 1(d): calibration workflow

Open `fig1d-v8-editable.svg` in Inkscape 1.4 or later. The page is 88 mm by
30 mm.

Layer order:

1. `01 Background`: white page and panel label.
2. `02 Transitions`: arrows, recovery paths, and endpoints.
3. `03 State Boxes`: the four equal-size node rectangles.
4. `04 State Text`: state names and physical-role subtitles.
5. `05 Edge Labels`: transition conditions and calibrated output label.

The main nodes are 18 mm by 7.2 mm, and the three top-row node centers are
separated by 28 mm. `fig1d-v8-reference.pdf` and
`fig1d-v8-reference.png` are reference exports.

## Editing notes

Fonts use Times New Roman with Latin Modern Roman as fallback. Keep text live
while editing; convert text to paths only in a final submission copy if the
publisher requires it.

## Archive

Superseded figure drafts are stored under `archive/` rather than deleted.
Files in the top-level figures directory are either active manuscript figures,
current editable Fig. 1 sources, or the scripts and metadata required to
rebuild them.
