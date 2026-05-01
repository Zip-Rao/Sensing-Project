"""Cryoscope experiment implementation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.experiments.base import Experiment
from sqc.hardware.readout import IQReadoutModel
from sqc.simulation.result import ExperimentResult


@dataclass
class CryoscopeExperiment(Experiment):
    """Cryoscope truncation scan matching legacy protocol case 5."""

    qubit: object
    flux_signal: FluxSignal | None = None
    t_rabi: np.ndarray = field(default_factory=lambda: np.linspace(0, 10, 20))
    tau: float = 100.0
    trunc_list: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.flux_signal is None:
            self.flux_signal = FluxSignal(
                type=2,
                t_list=np.linspace(0, 80, 160),
                amplitude=0.01,
            )
        if self.trunc_list is None:
            self.trunc_list = self.flux_signal.t_list[140:20:-1]

    def build_sequence(self):
        """Cryoscope uses IQ Ramsey pulses inside the readout model."""
        return None

    def run(self) -> ExperimentResult:
        """Run the truncation scan and return phase/IQ traces."""
        readout = IQReadoutModel(tau=self.tau, t_rabi=self.t_rabi)
        p_e_i = []
        p_e_q = []
        final_phi = None

        for trunc in self.trunc_list:
            phi = self.flux_signal.copy()
            phi.truncate(0, trunc)
            final_phi = phi
            self.qubit.qubit_in_mag(phi, frame=1, omega_d=self.qubit.frequency)
            measured = readout.measure(self.qubit)
            p_e_i.append(measured["p_e_I"])
            p_e_q.append(measured["p_e_Q"])

        p_e_i = np.asarray(p_e_i)
        p_e_q = np.asarray(p_e_q)
        varphi = np.arctan2(p_e_q - 0.5, p_e_i - 0.5)[::-1]

        return ExperimentResult(
            data={
                "varphi": varphi,
                "p_e_I": p_e_i,
                "p_e_Q": p_e_q,
                "flux_samples": np.asarray(final_phi.signal if final_phi else self.flux_signal.signal),
            },
            axes={
                "trunc": np.asarray(self.trunc_list),
                "t_flux": np.asarray(self.flux_signal.t_list),
            },
            metadata={
                "experiment": "CryoscopeExperiment",
                "qubit_spec": self.qubit.spec(),
                "flux_signal": final_phi or self.flux_signal,
            },
            config={"tau": self.tau, "t_rabi": np.asarray(self.t_rabi)},
        )
