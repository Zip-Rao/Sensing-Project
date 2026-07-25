# Worked Examples

These three pipelines ship with the platform as **worked examples** of composing
the building blocks into a complete sensing application. Each follows the same
four-part structure — goal, physics, end-to-end code, and how to read the
results — and every code block is runnable as shown.

- {doc}`waveform_reconstruction` — the qubit as a sensor: recover an unknown
  transient flux $\Phi(t)$ from excited-state population, and compare
  reconstruction algorithms, via a single `SensingWorkflow` entry point.
- {doc}`frequency_calibration` — characterise the device itself: measure the
  $f_{01}(\Phi)$ response, look up the bias for a target frequency, and
  closed-loop tune to it.
- {doc}`predistortion` — correct the control line: measure its transfer function,
  design an inverse filter, and verify that the on-chip waveform lands on target.

```{toctree}
:maxdepth: 1

waveform_reconstruction
frequency_calibration
predistortion
```
