"""Echo-family experiment implementations."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.control.flux_signal import CompositeSignal, FluxSignal
from sqc.control.sequence import create_diff_echo_pulse
from sqc.experiments.base import Experiment
from sqc.simulation.result import ExperimentResult


@dataclass
class DiffEchoExperiment(Experiment):
    """Differential echo protocol matching legacy protocol case 2."""

    qubit: object
    k: int = 5
    t_list: np.ndarray = field(default_factory=lambda: np.linspace(0, 100, 200))
    t_rabi: np.ndarray = field(default_factory=lambda: np.linspace(0, 10, 20))
    t_global: np.ndarray = field(default_factory=lambda: np.linspace(-10, 1010, 2020))
    flux_signal: FluxSignal | None = None
    omega_d: float | None = None

    def __post_init__(self) -> None:
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency
        if self.flux_signal is None:
            self.flux_signal = FluxSignal(
                type=3,
                t_list=self.t_list,
                amplitude=0.01,
                frequency=0.01,
                rise=10,
                fall=10,
                center=50,
                noise_level=0.0001,
            )
        self.t_int = (self.t_rabi[-1] - self.t_rabi[0]) * 0.5
        self.t_rep = self.flux_signal.t_list[-1] - self.flux_signal.t_list[0]

    def build_sequence(self):
        """Diff-echo control pulses are built per delay in ``run``."""
        return None

    def run(self) -> ExperimentResult:
        """Run the differential echo scan."""
        phi_list = [self.flux_signal.copy() for _ in range(2 * self.k)]
        composite_phi = CompositeSignal(phi_list)
        self.qubit.qubit_in_mag(composite_phi, frame=1, omega_d=self.omega_d)

        tau_list = np.asarray(self.flux_signal.t_list).copy()
        psi_e = basis(self.qubit.n_levels, 1)
        p_e = np.zeros(len(tau_list))
        for i, tau in enumerate(tau_list):
            control_pulse = create_diff_echo_pulse(
                self.t_rabi,
                tau,
                self.t_int,
                self.t_rep,
                k=self.k,
                omega_d=self.omega_d,
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
                "tau": tau_list,
                "t_flux": np.asarray(self.flux_signal.t_list),
            },
            metadata={
                "experiment": "DiffEchoExperiment",
                "qubit_spec": self.qubit.spec(),
            },
            config={
                "omega_d": self.omega_d,
                "k": self.k,
                "t_int": self.t_int,
                "t_rep": self.t_rep,
                "t_rabi": np.asarray(self.t_rabi),
                "t_global": np.asarray(self.t_global),
            },
        )


class EchoExperiment(DiffEchoExperiment):
    """Placeholder alias for an echo-family experiment."""
