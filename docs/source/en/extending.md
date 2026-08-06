# Extending the Platform

The three built-in pipelines are examples, not limits. Every layer of the stack
exposes an abstract base class; subclass it and the new object composes with the
rest of the platform. This page gives the two most common patterns.

## Add a new experiment protocol

Use case: a new standard experiment (ALLXY, T1, CPMG) or a custom flux-sensing
protocol.

1. Create `sqc/experiments/my_protocol.py`.
2. Subclass {py:class}`~sqc.experiments.Experiment` and implement
   `build_sequence()` and `run()`:

```python
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

from sqc.config import CONFIG
from sqc.experiments.base import Experiment
from sqc.control.flux_signal import FluxSignal
from sqc.simulation.result import ExperimentResult


@dataclass
class MyProtocolExperiment(Experiment):
    """My new sensing protocol.

    Parameters
    ----------
    qubit : TransmonQubit
    flux_signal : FluxSignal or None
    param1 : float
    """
    qubit: object
    flux_signal: FluxSignal | None = None
    param1: float = 1.0
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())

    def build_sequence(self):
        ...  # return a PulseSequence

    def run(self) -> ExperimentResult:
        ...  # run mesolve, return ExperimentResult(data=..., axes=..., metadata=...)
```

An {py:class}`~sqc.simulation.ExperimentResult` carries `data`, `axes`,
`metadata`, and `config` dicts, the uniform contract every experiment returns.

## Add a new reconstruction algorithm

Use case: a new inversion algorithm (deep-learning inversion, compressed
sensing, Volterra-kernel expansion).

1. Create `sqc/reconstruction/my_algo.py`.
2. Subclass {py:class}`~sqc.reconstruction.Reconstruction`:

```python
from __future__ import annotations
import numpy as np
from sqc.reconstruction.base import Reconstruction


class MyAlgoReconstruction(Reconstruction):
    """Reconstruct Phi(t) from a measurement using my algorithm."""

    def reconstruct(self, measurement, **kwargs) -> np.ndarray:
        p_e = measurement.data["p_e"]
        # ... invert p_e -> phi -> Phi(t) ...
        return reconstructed_flux
```

Once registered, it can be selected by name through
{py:meth}`SensingWorkflow.configure(reconstruction=...) <sqc.workflows.SensingWorkflow.configure>`
and compared against the built-in methods with
{py:meth}`~sqc.workflows.SensingWorkflow.compare`.

## Other extension points

| Layer | Base class | Add… |
|---|---|---|
| devices | {py:class}`~sqc.devices.Device` | a new qubit type (e.g. fluxonium) |
| hardware | {py:class}`~sqc.hardware.DistortionModel` | a new line-distortion model |
| control | {py:class}`~sqc.control.PulseBase` | a new pulse shape |
| simulation | {py:class}`~sqc.simulation.RunnerBase` | a new solver strategy |
| calibration | {py:class}`~sqc.calibration.Calibration` | a new calibration routine |
| workflows | {py:class}`~sqc.workflows.Workflow` | a new end-to-end pipeline |

Follow the dependency rule ({doc}`architecture`): a subclass may import from
layers below it, never above.
