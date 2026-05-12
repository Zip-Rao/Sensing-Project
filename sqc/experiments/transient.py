"""sqc.experiments.transient — TransientSensingExperiment.

Transient magnetic field sensing via sliding measurement protocol.
Replaces Protocal.evolve case 4.

Slides a Ramsey control pulse across a flux signal, measuring p_e
at each delay. Computes the control kernel for deconvolution.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.config import CONFIG
from sqc.experiments.base import Experiment
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.simulation.runner import SlidingMeasurementRunner
from sqc.simulation.result import ExperimentResult


@dataclass
class TransientSensingExperiment(Experiment):
    """Transient magnetic field sensing experiment.

    Default parameters match src/protocal.py case 4 exactly.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit object (src or sqc version).
    flux_signal : FluxSignal or None
        Flux signal. If None, creates default asymmetric impulse.
    flux_signal_zero : FluxSignal or None
        Zero-flux reference signal. If None, creates constant-zero.
    t_rabi : np.ndarray
        Rabi pulse time axis (ns). Default linspace(0, 10, 20).
    omega_d : float or None
        Drive frequency. Default qubit.frequency.
    scan_list : np.ndarray or None
        Pre-computed delay list for sliding measurement.
    """

    qubit: object  # TransmonQubit (duck typed)
    flux_signal: FluxSignal | None = None
    flux_signal_zero: FluxSignal | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    omega_d: float | None = None
    scan_list: np.ndarray | None = None

    # Populated during run()
    control_pulse: object | None = None

    def __post_init__(self):
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency

        if self.flux_signal is None:
            # Default test signal matches src/protocal.py case 4
            t_list = CONFIG.pulse.make_time(0, 200)
            self.flux_signal = FluxSignal(
                type=4,
                t_list=t_list,
                amplitude=0.01,
                rise=10,
                fall=10,
                center=100,
                noise_level=0.0001,
            )

        if self.flux_signal_zero is None:
            t_list = self.flux_signal.t_list
            self.flux_signal_zero = FluxSignal(
                type=1,
                t_list=t_list,
                amplitude=0.0,
            )

    def build_sequence(self):
        """Build the Ramsey control pulse (tau=0)."""
        return create_ramsey_pulse(
            self.t_rabi, tau=0.0, omega_d=self.omega_d,
            qubit=self.qubit,
        )

    def run(self) -> ExperimentResult:
        """Execute transient sensing experiment.

        Returns
        -------
        ExperimentResult
            With data["p_e"], data["delta_p"], data["kernel"],
            axes["scan"], axes["t_samples"].
        """
        # Build control pulse
        self.control_pulse = self.build_sequence()

        # Sliding measurement with flux signal
        runner = SlidingMeasurementRunner()
        result_sig = runner.run(
            self.qubit, self.flux_signal, self.control_pulse,
            scan_list=self.scan_list,
        )
        scan_list = result_sig.axes["scan"]
        p_e = result_sig.data["p_e"]

        # Sliding measurement with zero-flux reference
        result_base = runner.run(
            self.qubit, self.flux_signal_zero, self.control_pulse,
            scan_list=scan_list,
        )
        p_e_base = result_base.data["p_e"]

        # Compute kernel via legacy pulse.get_kernel()
        self.control_pulse.get_kernel(self.qubit)
        t_samples = self.control_pulse.t_samples
        kernel = self.control_pulse.kernel

        # Compute delta_p
        delta_p = np.array(p_e) - np.array(p_e_base)

        return ExperimentResult(
            data={
                "p_e": np.array(p_e),
                "delta_p": delta_p,
                "kernel": np.array(kernel),
                "flux_samples": self.flux_signal.signal.copy(),
            },
            axes={
                "scan": scan_list,
                "t_samples": np.array(t_samples),
                "t_flux": self.flux_signal.t_list.copy(),
            },
            metadata={
                "experiment": "TransientSensingExperiment",
                "omega_d": self.omega_d,
            },
            config={
                "t_rabi": self.t_rabi.copy(),
            },
        )
