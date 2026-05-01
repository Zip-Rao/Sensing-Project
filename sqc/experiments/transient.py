"""Transient flux sensing experiment implementation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.experiments.base import Experiment
from sqc.reconstruction.kernel import KernelEstimator
from sqc.simulation.result import ExperimentResult
from sqc.simulation.runner import SlidingMeasurementRunner


@dataclass
class TransientSensingExperiment(Experiment):
    """Sliding transient-sensing protocol matching legacy protocol case 4."""

    qubit: object
    t_list: np.ndarray = field(default_factory=lambda: np.linspace(0, 200, 400))
    t_rabi: np.ndarray = field(default_factory=lambda: np.linspace(0, 10, 20))
    flux_signal: FluxSignal | None = None
    baseline_flux: FluxSignal | None = None
    runner: SlidingMeasurementRunner = field(default_factory=SlidingMeasurementRunner)
    kernel_estimator: KernelEstimator = field(default_factory=KernelEstimator)

    def __post_init__(self) -> None:
        if self.flux_signal is None:
            self.flux_signal = FluxSignal(
                type=4,
                t_list=self.t_list,
                amplitude=0.06,
                rise=10,
                fall=10,
                center=100,
                noise_level=0.0001,
            )
        if self.baseline_flux is None:
            self.baseline_flux = FluxSignal(type=1, t_list=self.t_list, amplitude=0.0)
        self.control_pulse = create_ramsey_pulse(
            self.t_rabi,
            tau=0.0,
            omega_d=self.qubit.frequency,
        )

    def build_sequence(self):
        """Return the Ramsey control pulse used for sliding measurement."""
        return self.control_pulse

    def run(self) -> ExperimentResult:
        """Run transient sensing, baseline subtraction, and kernel estimation."""
        measurement = self.runner.run(
            self.qubit, self.flux_signal, self.control_pulse
        )
        baseline = self.runner.run(
            self.qubit, self.baseline_flux, self.control_pulse
        )
        t_samples, kernel = self.kernel_estimator.estimate(
            self.control_pulse, self.qubit
        )
        self.control_pulse.t_samples = t_samples
        self.control_pulse.kernel = kernel

        p_e = measurement.data["p_e"]
        p_e_base = baseline.data["p_e"]
        delta_p = p_e - p_e_base
        return ExperimentResult(
            data={
                "kernel": np.asarray(kernel),
                "delta_p": np.asarray(delta_p),
                "p_e": np.asarray(p_e),
                "p_e_base": np.asarray(p_e_base),
                "flux_samples": np.asarray(self.flux_signal.signal),
            },
            axes={
                "t_samples": np.asarray(t_samples),
                "scan": np.asarray(measurement.axes["scan"]),
                "t_flux": np.asarray(self.flux_signal.t_list),
            },
            metadata={
                "experiment": "TransientSensingExperiment",
                "qubit_spec": self.qubit.spec(),
            },
            config={
                "t_rabi": np.asarray(self.t_rabi),
            },
        )
