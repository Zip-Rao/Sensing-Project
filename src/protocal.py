"""Compatibility facade for legacy protocol imports.

The historical class name ``Protocal`` is preserved. Implemented protocol cases
delegate to ``sqc.experiments`` while keeping legacy return shapes.
"""
from __future__ import annotations

import numpy as np
from qutip import Qobj, basis

from sqc.control.flux_signal import CompositeSignal, FluxSignal as Signal
from sqc.control.pulse import CompositePulse, Pulse
from sqc.control.sequence import (
    create_cpmg_pulse,
    create_cryoscope_pulse,
    create_diff_echo_pulse,
    create_echo_pulse,
    create_pulse,
    create_ramsey_pulse,
)
from sqc.devices.transmon import TransmonQubit


class Protocal:
    """Legacy protocol class, now a facade over sqc experiment objects."""

    def __init__(self, type=0, **kwargs):
        self.type = type
        params = dict(kwargs)
        nested = params.pop("params", None)
        if isinstance(nested, dict):
            nested.update(params)
            params = nested
        self.params = params

    def initialize(self, qubit: TransmonQubit, state=0):
        """Initialize qubit state and fill historical default parameters."""
        if isinstance(state, int) and 0 <= state < qubit.n_levels:
            qubit.state = basis(qubit.n_levels, state)
        elif isinstance(state, Qobj) and state.dims == [[qubit.n_levels], [1]]:
            qubit.state = state.unit()

        default_params = {
            "t_global": np.linspace(-50, 300, 700),
            "t_list": np.linspace(0, 100, 1000),
            "tau_list": np.linspace(0, 100, 100),
            "t_rabi": np.linspace(0, 40, 100),
        }
        for key, value in default_params.items():
            self.params.setdefault(key, value)

    def evolve(self, qubit):
        """Run the selected legacy protocol."""
        match self.type:
            case 0:
                from sqc.experiments.rabi import RabiExperiment

                result = RabiExperiment(qubit=qubit).run()
                return result.metadata["raw_result"]

            case 1:
                from sqc.experiments.ramsey import RamseyExperiment

                exp = RamseyExperiment(qubit=qubit)
                result = exp.run()
                return exp.flux_signal, result.axes["tau"], result.data["p_e"].tolist()

            case 2:
                from sqc.experiments.echo import DiffEchoExperiment

                exp = DiffEchoExperiment(qubit=qubit)
                result = exp.run()
                return (
                    exp.flux_signal,
                    result.axes["tau"],
                    result.data["p_e"].tolist(),
                    result.config["k"],
                    result.config["t_int"],
                )

            case 3:
                return None

            case 4:
                from sqc.experiments.transient import TransientSensingExperiment

                exp = TransientSensingExperiment(qubit=qubit)
                result = exp.run()
                return (
                    result.axes["t_samples"],
                    result.data["kernel"],
                    result.axes["scan"],
                    result.data["delta_p"],
                    result.data["p_e"],
                    exp.flux_signal,
                    exp.control_pulse,
                )

            case 5:
                from sqc.experiments.cryoscope import CryoscopeExperiment

                exp = CryoscopeExperiment(qubit=qubit)
                result = exp.run()
                phi = result.metadata["flux_signal"]
                return (
                    result.axes["trunc"],
                    result.data["varphi"],
                    phi,
                    [result.data["p_e_I"].tolist(), result.data["p_e_Q"].tolist()],
                )

            case 6 | 7 | 8:
                raise NotImplementedError(
                    f"Protocol case {self.type} awaits Track B implementation."
                )

            case _:
                raise ValueError(f"Unknown protocol type {self.type}")

    def single_measurement(
        self,
        qubit: TransmonQubit,
        Phi_signal: Signal,
        control_pulse: CompositePulse,
        t_delay,
        qubit_t=None,
        H=None,
        t_evole=None,
        index=None,
    ):
        """Backward-compatible single delayed measurement helper."""
        from sqc.simulation.runner import SlidingMeasurementRunner

        return SlidingMeasurementRunner().single_measurement(
            qubit,
            Phi_signal,
            control_pulse,
            t_delay,
            qubit_t=qubit_t,
            hamiltonian=H,
            t_evolve=t_evole,
        )

    def sliding_measrement(
        self, qubit: TransmonQubit, Phi_signal: Signal, control_pulse: CompositePulse
    ):
        """Backward-compatible sliding measurement helper.

        The misspelling is preserved for old notebooks.
        """
        from sqc.simulation.runner import SlidingMeasurementRunner

        result = SlidingMeasurementRunner().run(qubit, Phi_signal, control_pulse)
        return result.axes["scan"], result.data["p_e"].tolist()


class Calibration:
    """Legacy facade over sqc.calibration classes."""

    def __init__(self, qubit, type=0, **kwargs):
        self.qubit = qubit
        self.type = type
        self.params = kwargs

    def calibrate(self):
        from sqc.calibration.flux_response import FluxResponseCalibration
        from sqc.calibration.qubit_frequency import (
            QubitFrequencyCalibration,
            TransientFrequencyCalibration,
        )

        match self.type:
            case 0:
                return QubitFrequencyCalibration(qubit=self.qubit).calibrate()
            case 1:
                return FluxResponseCalibration(
                    qubit=self.qubit, method="ramsey", **self.params
                ).calibrate()
            case 2:
                return TransientFrequencyCalibration(
                    qubit=self.qubit, **self.params
                ).calibrate()
            case 3:
                calibration = FluxResponseCalibration(
                    qubit=self.qubit, method="cryoscope", **self.params
                )
                table = calibration.calibrate()
                return table.inputs, table.outputs, calibration.tau
            case _:
                raise ValueError(f"Unknown calibration type {self.type}")


def IQ_readout(qubit, type, **kwargs):
    """Legacy IQ readout function."""
    from sqc.hardware.readout import IQ_readout_legacy

    return IQ_readout_legacy(qubit, type, **kwargs)


__all__ = [
    "Protocal",
    "Calibration",
    "IQ_readout",
    "TransmonQubit",
    "Signal",
    "CompositeSignal",
    "Pulse",
    "CompositePulse",
    "create_pulse",
    "create_ramsey_pulse",
    "create_diff_echo_pulse",
    "create_echo_pulse",
    "create_cpmg_pulse",
    "create_cryoscope_pulse",
]
