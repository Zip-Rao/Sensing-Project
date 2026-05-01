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
from sqc.experiments.cryoscope import CryoscopeExperiment
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
                # Ported from src/protocal.py case 5.
                # Creates CryoscopeExperiment with legacy defaults.
                # NOTE: Legacy calls Phi.plot() on the flux signal.
                t_rabi = np.linspace(0, 10, 20)
                tau = 100.0
                exp = CryoscopeExperiment(
                    qubit=qubit,
                    t_rabi=t_rabi,
                    tau=tau,
                )
                # Pre-plot flux signal (legacy behaviour)
                exp.flux_signal.plot()

                result = exp.run()
                # Legacy returns: trunc_list, varphi_list, Phi, p_e_list
                return (
                    result.axes["trunc"],
                    result.data["varphi"],
                    exp.flux_signal,
                    [list(result.data["p_e_I"]), list(result.data["p_e_Q"])],
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
# Calibration facade (P3c) — delegates to sqc.calibration.*
# ---------------------------------------------------------------------------

class Calibration:
    """Legacy calibration class — facade over sqc.calibration.*.

    API preserved verbatim from src/protocal.py:Calibration.

    Delegation:
      - type 0 → QubitFrequencyCalibration (Ramsey f_01)
      - type 1 → FluxResponseCalibration(method="ramsey") (f(Phi) via Ramsey)
      - type 2 → NOT IMPLEMENTED (requires Track B 1.2)
      - type 3 → NOT IMPLEMENTED (requires Track B 1.1)
    """

    def __init__(self, qubit, type=0, **kwargs):
        self.qubit = qubit
        self.type = type
        self.params = kwargs

    def calibrate(self):
        """Execute calibration via sqc/calibration/ classes.

        Returns
        -------
        - type 0: CalibrationTable (f01)
        - type 1: CalibrationTable (f_phi)
        - type 2/3: raises NotImplementedError
        """
        from sqc.calibration.qubit_frequency import QubitFrequencyCalibration
        from sqc.calibration.flux_response import FluxResponseCalibration

        match self.type:
            case 0:  # Ramsey frequency f_01 calibration
                cal = QubitFrequencyCalibration(
                    qubit=self.qubit,
                )
                return cal.calibrate()

            case 1:  # f(Phi) via Ramsey
                cal = FluxResponseCalibration(
                    qubit=self.qubit,
                    method="ramsey",
                )
                return cal.calibrate()

            case 2:  # Transient calib — requires Track B 1.2
                raise NotImplementedError(
                    "Calibration type=2 (transient frequency calibration): "
                    "requires Track B 1.2 (case 8). See _TODO_master.md 1.2."
                )

            case 3:  # Cryoscope φ(h) calib — requires Track B 1.1
                raise NotImplementedError(
                    "Calibration type=3 (cryoscope φ(h) calibration): "
                    "requires Track B 1.1 (Cryoscope case 6/7). "
                    "See _TODO_master.md 1.1. "
                    "The legacy src/protocal.py:Calibration(type=3).calibrate() "
                    "has a working implementation; use that for now."
                )

            case _:
                raise ValueError(f"Unknown calibration type {self.type}")
