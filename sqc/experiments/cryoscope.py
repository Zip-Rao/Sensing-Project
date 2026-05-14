"""sqc.experiments.cryoscope — CryoscopeExperiment.

Cryoscope: scan truncation delay, measure φ(t_d) via IQ Ramsey.
Replaces Protocal.evolve case 5.

Default parameters match src/protocal.py case 5 exactly:
  - Phi: type=2, t_list=linspace(0, 80, 160), amplitude=0.01
  - trunc_list: Phi.t_list[140:20:-1] (reverse order)
  - tau: 100.0
  - t_rabi: linspace(0, 10, 20)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.config import CONFIG
from sqc.experiments.base import Experiment
from sqc.hardware.readout import IQReadoutModel
from sqc.simulation.result import ExperimentResult


@dataclass
class CryoscopeExperiment(Experiment):
    """Cryoscope: scan truncation delay, measure φ(t_d) via IQ Ramsey.

    Direct port of src/protocal.py:Protocal.evolve case 5.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit object (src or sqc version).
    flux_signal : FluxSignal or None
        Flux signal to scan. If None, creates default sinusoidal signal
        matching src/protocal.py case 5.
    t_rabi : np.ndarray
        Rabi pulse time axis (ns). Default linspace(0, 10, 20).
    tau : float
        Free precession time for IQ readout (ns). Default 100.0.
    trunc_list : np.ndarray or None
        List of truncation times. If None, defaults to
        flux_signal.t_list[140:20:-1] (reverse order).
    omega_d : float or None
        Drive frequency. If None, uses qubit.frequency.
    """

    qubit: object  # TransmonQubit (duck-typed)
    flux_signal: FluxSignal | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    tau: float = field(
        default_factory=lambda: CONFIG.reconstruction.cryoscope_tau
    )
    trunc_list: np.ndarray | None = None
    omega_d: float | None = None

    def __post_init__(self):
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency
        if self.flux_signal is None:
            self.flux_signal = FluxSignal(
                type=2,
                t_list=CONFIG.pulse.make_time(0, 80),
                amplitude=0.01,
            )
        if self.trunc_list is None:
            self.trunc_list = self.flux_signal.t_list[140:20:-1]

    def build_sequence(self):
        """Cryoscope sequence is constructed per-truncation in run()."""
        return None

    def run(self) -> ExperimentResult:
        """Execute Cryoscope experiment.

        For each truncation delay, truncates the flux signal, couples
        it to the qubit, and performs IQ readout to measure the
        accumulated phase.

        Returns
        -------
        ExperimentResult
            With data["varphi"], data["p_e_I"], data["p_e_Q"],
            axes["trunc"].
        """
        readout = IQReadoutModel(
            tau=self.tau,
            t_rabi=self.t_rabi,
            omega_d=self.omega_d,
        )
        p_e_I_list: list[float] = []
        p_e_Q_list: list[float] = []

        # Iterate truncation delays in reverse order (matching legacy)
        for trunc in self.trunc_list:
            # Copy flux signal and truncate in-place (legacy semantics)
            phi_truncated = self.flux_signal.copy()
            phi_truncated.truncate(0, float(trunc))

            # Couple flux to qubit
            self.qubit.qubit_in_mag(
                phi_truncated, frame=1, omega_d=self.omega_d,
            )

            # IQ readout
            measured = readout.measure(self.qubit)
            p_e_I_list.append(measured["p_e_I"])
            p_e_Q_list.append(measured["p_e_Q"])

        p_e_I = np.asarray(p_e_I_list, dtype=float)
        p_e_Q = np.asarray(p_e_Q_list, dtype=float)

        # Compute phase via arctan2, reverse to align time axis (legacy)
        varphi = np.arctan2(0.5 - p_e_I, p_e_Q - 0.5)[::-1]
        varphi = np.unwrap(varphi)
        return ExperimentResult(
            data={
                "varphi": varphi,
                "p_e_I": p_e_I,
                "p_e_Q": p_e_Q,
            },
            axes={
                "trunc": np.asarray(self.trunc_list),
            },
            metadata={
                "experiment": "CryoscopeExperiment",
                "tau": self.tau,
            },
            config={
                "tau": self.tau,
                "t_rabi": self.t_rabi.copy(),
                "omega_d": self.omega_d,
            },
        )
