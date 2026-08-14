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
        Rabi pulse time axis (ns). Default ``CONFIG.pulse.t_rabi``.
    omega_d : float or None
        Drive frequency. Default qubit.frequency.
    rotation_angle : float or None
        Rotation angle of each sensing pulse (rad). Default pi/2. Set to None
        when using ``rabi_rate``.
    rabi_rate : float or None
        Fixed peak Rabi rate. When set, the achieved angle follows from the
        pulse duration and envelope.
    envelope : {'square', 'gaussian'} or array-like
        Microwave pulse envelope. Default square.
    envelope_sigma : float or None
        Gaussian envelope standard deviation (ns).
    phase1, phase2 : float
        Rotation-axis phases of the two adjacent pulses (rad).
    scan_list : np.ndarray or None
        Pre-computed delay list for sliding measurement.
    """

    qubit: object  # TransmonQubit (duck typed)
    flux_signal: FluxSignal | None = None
    flux_signal_zero: FluxSignal | None = None
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())
    omega_d: float | None = None
    rotation_angle: float | None = np.pi / 2
    rabi_rate: float | None = None
    envelope: object = "square"
    envelope_sigma: float | None = None
    phase1: float = np.pi / 2
    phase2: float = 0.0

    # -- P9.B --
    control_line: object | None = None
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
            self.t_rabi,
            tau=0.0,
            omega_d=self.omega_d,
            phase1=self.phase1,
            phase2=self.phase2,
            qubit=self.qubit,
            rotation_angle=self.rotation_angle,
            rabi_rate=self.rabi_rate,
            envelope=self.envelope,
            envelope_sigma=self.envelope_sigma,
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
            self.qubit,
            self._route_flux(self.flux_signal),
            self.control_pulse,
            scan_list=self.scan_list,
        )
        scan_list = result_sig.axes["scan"]
        p_e = result_sig.data["p_e"]

        # Sliding measurement with zero-flux reference
        result_base = runner.run(
            self.qubit,
            self.flux_signal_zero,
            self.control_pulse,
            scan_list=scan_list,
        )
        p_e_base = result_base.data["p_e"]

        # Compute the control kernel directly via KernelEstimator (flux / exp,
        # order 1) — the same computation the deprecated CompositePulse.get_kernel
        # shim delegates to, but on the unified route (matching the frontend and
        # tests) and without the deprecation warning. (Local import avoids any
        # experiments<->reconstruction import cycle.)
        from sqc.reconstruction.kernel import KernelEstimator

        kernel_result = KernelEstimator(
            mode="flux", method="exp", order=1
        ).estimate_full(self.control_pulse, self.qubit)
        t_samples = kernel_result.t_samples
        kernel = list(kernel_result.k1)
        effective_angle = float(self.control_pulse.pulses[0].get_angle_simple())

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
                "rotation_angle": effective_angle,
                "rabi_rate": self.rabi_rate,
                "envelope": (
                    self.envelope if isinstance(self.envelope, str) else "custom"
                ),
                "phase1": self.phase1,
                "phase2": self.phase2,
            },
            config={
                "t_rabi": self.t_rabi.copy(),
                "rotation_angle": effective_angle,
                "rabi_rate": self.rabi_rate,
                "envelope_sigma": self.envelope_sigma,
            },
        )
