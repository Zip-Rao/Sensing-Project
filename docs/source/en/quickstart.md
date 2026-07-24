# Quickstart

The fastest way to run an end-to-end sensing pipeline is the high-level
{py:class}`~sqc.workflows.SensingWorkflow`. It composes a qubit device, a flux
signal, a control sequence, a QuTiP simulation, and a reconstruction algorithm
behind one fluent interface.

## Sense and reconstruct a flux waveform

```python
from sqc.workflows import SensingWorkflow

result = (
    SensingWorkflow()
    .configure(
        protocol="ramsey",         # sensing protocol
        signal_type=3,             # flux signal shape (Gaussian pulse)
        signal_amplitude=0.01,     # in units of Phi_0
        reconstruction="wiener",   # reconstruction algorithm
    )
    .run(measure=True, reconstruct=True)
)

print(result)          # summary of the run
# result carries the raw measurement and the reconstructed waveform
```

`configure()` only changes the parameters you pass; everything else keeps its
default (derived from {py:data}`sqc.config.CONFIG`). It returns `self`, so calls
chain. `run()` builds the experiment, calls QuTiP `mesolve`, and reconstructs.

## Compare reconstruction algorithms

```python
wf = SensingWorkflow().configure(protocol="ramsey", signal_type=3)
comparison = wf.compare(methods=["wiener", "lm"])
```

## Sweep a parameter

```python
wf = SensingWorkflow().configure(protocol="ramsey")
sweep = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02])
```

## Where to go next

- {doc}`architecture` — the eight-layer stack the workflow is built on.
- {doc}`building_blocks/index` — use each layer directly for full control.
- {doc}`examples/index` — the three built-in pipelines in depth.
- {doc}`extending` — build your own sensing application.
