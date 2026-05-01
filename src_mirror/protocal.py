"""src_mirror.protocal — Facade over sqc.experiments.

Backward-compatible Protocal class. Same API as src/protocal.py,
but internal implementation delegates to sqc/experiments/ classes.

Legacy API preserved verbatim:
  Protocal(type=N).evolve(qubit) -> tuple-of-things-matching-old-behavior

R1 compliance: NEVER modifies src/. This is a NEW file in src_mirror/.

See idea/refactor/_refactor_plan.md §8.
"""
from __future__ import annotations

import numpy as np
from qutip import Qobj, basis

from sqc.devices.transmon import TransmonQubit
from sqc.control.flux_signal import FluxSignal as Signal, CompositeSignal
from sqc.control.pulse import Pulse, CompositePulse
from sqc.control.sequence import (
    create_pulse,
    create_ramsey_pulse,
    create_diff_echo_pulse,
    create_cpmg_pulse,
    create_cryoscope_pulse,
)
from sqc.hardware.readout import IQ_readout_legacy as IQ_readout
from sqc.experiments.rabi import RabiExperiment
from sqc.experiments.ramsey import RamseyExperiment
from sqc.experiments.echo import DiffEchoExperiment
from sqc.experiments.transient import TransientSensingExperiment
from sqc.reconstruction.kernel import KernelEstimator


class Protocal:
    """Legacy protocol class, now a thin facade over sqc/experiments/.

    API and return types preserved verbatim from src/protocal.py.
    """

    def __init__(self, type=0, **kwargs):
        """Initialize protocol.

        Parameters
        ----------
        type : int
            Protocol type: 0=Rabi, 1=Ramsey, 2=DiffEcho, 3=CPMG (stub),
            4=Transient, 5=Cryoscope (stub).
        **kwargs
            Protocol parameters.
        """
        self.type = type
        self.params = kwargs

    def initialize(self, qubit: TransmonQubit, state=0):
        """Initialize qubit state and set default parameters.

        Matches src/protocal.py:Protocal.initialize exactly.

        Parameters
        ----------
        qubit : TransmonQubit
            Qubit to initialize.
        state : int or Qobj
            Initial state (0=ground, 1=excited, or state vector).
        """
        n_levels = qubit.n_levels
        if isinstance(state, int) and 0 <= state < n_levels:
            qubit.state = basis(n_levels, state)
        elif (
            isinstance(state, Qobj)
            and state.dims == [[n_levels], [1]]
        ):
            qubit.state = state.unit()

        default_params = {
            "t_global": np.linspace(-50, 300, 700),
            "t_list": np.linspace(0, 100, 1000),
            "tau_list": np.linspace(0, 100, 100),
            "t_rabi": np.linspace(0, 40, 100),
        }
        for key, value in default_params.items():
            if key not in self.params:
                self.params[key] = value

    def evolve(self, qubit):
        """Evolve qubit according to protocol type.

        Returns the SAME tuple structures as src/protocal.py:Protocal.evolve.

        Parameters
        ----------
        qubit : TransmonQubit
            Qubit object (src or sqc version).

        Returns
        -------
        result : tuple or qutip.Result
            Case-dependent return value.
        """
        match self.type:
            case 0:  # Rabi oscillation
                # Legacy hardcodes t_rabi = linspace(0, 40, 1000)
                # Parameters set in initialize() are NOT used by legacy case 0.
                exp = RabiExperiment(
                    qubit=qubit,
                    t_rabi=np.linspace(0, 40, 1000),
                )
                return exp.run()

            case 1:  # Ramsey
                # Legacy hardcodes all parameters; ignores initialize() defaults.
                # Matching src/protocal.py lines 62-98 exactly.
                exp = RamseyExperiment(
                    qubit=qubit,
                    t_rabi=np.linspace(0, 20, 40),
                    tau_list=np.linspace(0, 250, 500),
                    t_global=np.linspace(-50, 300, 700),
                )
                # NOTE: Legacy calls Phi.plot() here. Optional: exp.flux_signal.plot()

                result = exp.run()
                return (
                    exp.flux_signal,
                    result.axes["tau"],
                    list(result.data["p_e"]),
                )

            case 2:  # Differential echo
                # Legacy hardcodes all parameters (lines 117-141).
                exp = DiffEchoExperiment(qubit=qubit)
                # NOTE: Legacy calls composite_phi.plot() here.
                result = exp.run()
                return (
                    exp.flux_signal,
                    result.axes["tau"],
                    list(result.data["p_e"]),
                    exp.k,
                    exp.t_int,
                )

            case 3:  # CPMG
                pass  # stub, P3+

            case 4:  # Transient sensing
                exp = TransientSensingExperiment(qubit=qubit)
                # NOTE: Legacy calls Phi.plot() here.
                result = exp.run()
                # Legacy returns (t_samples, kernel, scan_list,
                #                delta_p, p_e, Phi, control_pulse)
                return (
                    result.axes["t_samples"],
                    result.data["kernel"],
                    result.axes["scan"],
                    result.data["delta_p"],
                    result.data["p_e"],
                    exp.flux_signal,
                    exp.control_pulse,
                )

            case 5:  # Cryoscope
                # P3 will internalize; for now, raise NotImplementedError
                raise NotImplementedError(
                    "Cryoscope protocol will be internalized in Phase 3"
                )

            case _:
                raise ValueError(f"Unknown protocol type {self.type}")

    # -- Legacy helper methods -------------------------------------------------

    def single_measurement(
        self,
        qubit,
        Phi_signal,
        control_pulse,
        t_delay,
        qubit_t=None,
        H=None,
        t_evole=None,
        index=None,
    ):
        """Backward-compat single measurement.

        Legacy spelling and signature preserved. Delegates to
        SlidingMeasurementRunner internal logic.

        Parameters
        ----------
        qubit : TransmonQubit
        Phi_signal : Signal-like
        control_pulse : CompositePulse
        t_delay : float
        qubit_t, H, t_evole, index : optional
            Pre-computed values. If not provided, computed from scratch.

        Returns
        -------
        float
            Excited-state probability p_e.
        """
        from sqc.simulation.runner import SlidingMeasurementRunner

        runner = SlidingMeasurementRunner()
        return runner._single_measurement(
            qubit, Phi_signal, control_pulse, t_delay,
            qubit_t, H, t_evole,
        )

    def sliding_measrement(
        self, qubit, Phi_signal, control_pulse
    ):
        """Backward-compat sliding measurement (legacy spelling preserved).

        Parameters
        ----------
        qubit : TransmonQubit
        Phi_signal : Signal-like
        control_pulse : CompositePulse

        Returns
        -------
        tuple[np.ndarray, list[float]]
            (scan_list, p_e)
        """
        from sqc.simulation.runner import SlidingMeasurementRunner

        runner = SlidingMeasurementRunner()
        result = runner.run(qubit, Phi_signal, control_pulse)
        return result.axes["scan"], list(result.data["p_e"])


# ---------------------------------------------------------------------------
# Calibration (verbatim stub, P3 will implement)
# ---------------------------------------------------------------------------

class Calibration:
    """Calibration class (stub for P3).

    Currently mirrors src/protocal.py:Calibration interface.
    Full implementation deferred to Phase 3.
    """

    def __init__(self, qubit, type=0, **kwargs):
        self.qubit = qubit
        self.type = type
        self.params = kwargs

    def calibrate(self):
        """Execute calibration."""
        match self.type:
            case 0:  # Ramsey frequency calib
                pass
            case 1:  # f(Phi) via Ramsey
                pass
            case 2:  # Transient calib
                pass
            case 3:  # Cryoscope calib
                raise NotImplementedError(
                    "Cryoscope calibration will be implemented in Phase 3"
                )
            case _:
                raise ValueError(f"Unknown calibration type {self.type}")
