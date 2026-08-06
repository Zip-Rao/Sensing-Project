# PRA paper workspace

> **开始任何论文讨论、写作或润色前，必须先阅读
> [`../WRITING_BOUNDARY.md`](../WRITING_BOUNDARY.md)。** 正式文字由用户亲自
> 撰写，AI 既有草稿仅作参考；“探究短脉冲序列参与标定的可能”是贯穿整个
> 写作过程的基本叙事，而不只是提纲层面的约束。

This directory is a working manuscript investigating how short-pulse sequences
can participate in superconducting-qubit calibration.  Frequency calibration
of a flux-tunable qubit is the concrete case study; response-kernel inversion
and drive tracking provide the local estimator and its operating-regime
control.  The manuscript uses REVTeX 4.2 with the `aps,pra,reprint` options and
the APS BibTeX style selected by REVTeX.

## Build

From this directory:

```powershell
latexmk -pdf main.tex
latexmk -pdf supplemental.tex
```

The local `.latexmkrc` places generated files in `build/`.  The PowerShell
helper runs both builds and returns a nonzero exit code on failure:

```powershell
.\build.ps1
```

## Layout

- `main.tex`: PRA manuscript entry point and author metadata.
- `sections/`: section-level source files.
- `supplemental.tex`: Supplemental Material skeleton.
- `references.bib`: curated references used by the manuscript.
- `figures/`: frozen frequency-calibration figures used by the manuscript.
- `notes/writing-plan.md`: scope and section plan.
- `notes/claim-evidence.md`: evidence gate for quantitative claims.
- `notes/pra-checklist.md`: pre-submission formatting and policy checklist.
- `data/`: publication-specific frozen data and provenance records.
- `archive/v0-framework-draft/`: superseded broad sensing-framework draft.

## Important status

This is a substantive numerical first draft, not a submission-ready manuscript.
Confirm how broadly the short-pulse framing should extend beyond frequency
calibration, together with authorship, experimental scope, bibliography
metadata, and funding.  The current evidence establishes one numerical case
study and a solver-cost reduction in simulation, not a universal calibration
advantage or an experimental speedup.
