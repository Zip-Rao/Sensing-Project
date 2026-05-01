"""Rabi experiment implementation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.control.sequence import create_pulse
from sqc.experiments.base import Experiment
from sqc.simulation.result import ExperimentResult


@dataclass
class RabiExperiment(Experiment):
    """Single-qubit Rabi oscillation experiment."""

    qubit: object
    t_rabi: np.ndarray = field(default_factory=lambda: np.linspace(0, 40, 1000))
    omega_d: float | None = None

    def __post_init__(self) -> None:
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency

    def build_sequence(self):
        """Build the resonant drive used by the experiment."""
        return create_pulse(
            self.qubit,
            1,
            1,
            self.t_rabi,
            self.omega_d,
            phase=0.0,
        )

    def run(self) -> ExperimentResult:
        """Run the Rabi experiment and return structured output."""
        psi_e = basis(self.qubit.n_levels, 1)
        h_0 = self.qubit.get_hamiltonian_rwa(self.qubit.frequency)
        h_pulse = self.build_sequence()
        h_rabi = QobjEvo(h_0) + QobjEvo(h_pulse, tlist=self.t_rabi)
        raw = mesolve(
            h_rabi,
            self.qubit.state,
            self.t_rabi,
            [],
            e_ops=[psi_e * psi_e.dag()],
        )
        return ExperimentResult(
            data={"p_e": np.asarray(raw.expect[0])},
            axes={"t": np.asarray(self.t_rabi)},
            metadata={
                "experiment": "RabiExperiment",
                "qubit_spec": self.qubit.spec(),
                "raw_result": raw,
            },
            config={"omega_d": self.omega_d},
        )
