"""Ramsey experiment implementation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.experiments.base import Experiment
from sqc.simulation.result import ExperimentResult


@dataclass
class RamseyExperiment(Experiment):
    """Ramsey protocol with the legacy default flux test signal."""

    qubit: object
    flux_signal: FluxSignal | None = None
    omega_d: float | None = None
    t_rabi: np.ndarray = field(default_factory=lambda: np.linspace(0, 20, 40))
    tau_list: np.ndarray = field(default_factory=lambda: np.linspace(0, 250, 500))
    t_global: np.ndarray = field(default_factory=lambda: np.linspace(-50, 300, 700))
    phase1: float = 0.0
    phase2: float = 0.0

    def __post_init__(self) -> None:
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency
        if self.flux_signal is None:
            self.flux_signal = FluxSignal(
                type=2,
                t_list=np.linspace(0, 250, 500),
                amplitude=0.001,
                frequency=0.01,
                rise=10,
                fall=10,
                center=100,
                noise_level=0.0,
            )

    def build_sequence(self):
        """Ramsey control pulses are built per delay in ``run``."""
        return None

    def run(self) -> ExperimentResult:
        """Run the Ramsey delay scan."""
        self.qubit.qubit_in_mag(self.flux_signal, frame=1, omega_d=self.omega_d)
        psi_e = basis(self.qubit.n_levels, 1)
        p_e = np.zeros(len(self.tau_list))

        for i, tau in enumerate(self.tau_list):
            control_pulse = create_ramsey_pulse(
                self.t_rabi,
                tau,
                omega_d=self.omega_d,
                phase1=self.phase1,
                phase2=self.phase2,
            )
            control_pulse.t_list = control_pulse.t_list - self.t_rabi[-1]
            h_total = QobjEvo(
                self.qubit.H_list,
                tlist=self.qubit.mag_signal.t_list,
                order=1,
            ) + QobjEvo(control_pulse.hamiltonian, tlist=control_pulse.t_list, order=1)
            result = mesolve(
                h_total,
                self.qubit.state,
                self.t_global,
                [],
                e_ops=[psi_e * psi_e.dag()],
            )
            p_e[i] = result.expect[0][-1]

        return ExperimentResult(
            data={
                "p_e": p_e,
                "flux_samples": np.asarray(self.flux_signal.signal),
            },
            axes={
                "tau": np.asarray(self.tau_list),
                "t_flux": np.asarray(self.flux_signal.t_list),
            },
            metadata={
                "experiment": "RamseyExperiment",
                "qubit_spec": self.qubit.spec(),
            },
            config={
                "omega_d": self.omega_d,
                "phase1": self.phase1,
                "phase2": self.phase2,
                "t_rabi": np.asarray(self.t_rabi),
                "t_global": np.asarray(self.t_global),
            },
        )
