"""sqc.experiments.ramsey — RamseyExperiment.

Ramsey interferometry with optional flux signal.
Replaces Protocal.evolve case 1.

Physical model: pi/2 - tau - pi/2 sequence measures accumulated
phase during free evolution, which reflects qubit frequency shifts
caused by external flux Phi(t).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.experiments.base import Experiment
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.simulation.result import ExperimentResult


@dataclass
class RamseyExperiment(Experiment):
    """Ramsey protocol with optional flux signal.

    Default parameters match src/protocal.py case 1 exactly.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit object (src or sqc version).
    flux_signal : FluxSignal or None
        Flux signal. If None, creates default sinusoidal signal.
    omega_d : float or None
        Drive frequency. Default qubit.frequency.
    t_rabi : np.ndarray
        Rabi pulse time axis (ns). Default linspace(0, 20, 40).
    tau_list : np.ndarray
        Free evolution times (ns). Default linspace(0, 250, 500).
    t_global : np.ndarray
        Global evolution time (ns). Default linspace(-50, 300, 700).
    phase1 : float
        First pi/2 pulse phase (rad). Default 0.
    phase2 : float
        Second pi/2 pulse phase (rad). Default 0.
    """

    qubit: object  # TransmonQubit (duck typed)
    flux_signal: FluxSignal | None = None
    omega_d: float | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    tau_list: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.tau_list.copy()
    )
    t_global: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_global.copy()
    )
    phase1: float = 0.0
    phase2: float = 0.0

    def __post_init__(self):
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency
        if self.flux_signal is None:
            # Default test signal matches src/protocal.py case 1
            self.flux_signal = FluxSignal(
                type=1,
                t_list=CONFIG.pulse.t_signal.copy(),
                amplitude=0.001,
                frequency=0.01,
                rise=10,
                fall=10,
                center=100,
                noise_level=0.000,
            )

    def build_sequence(self):
        """Ramsey sequence is constructed per-tau in run()."""
        return None

    def run(self) -> ExperimentResult:
        """Execute Ramsey experiment on unified global time axis.

        Returns
        -------
        ExperimentResult
            With data["p_e"], axes["tau"], data["flux_samples"].
        """
        t_global = self.t_global

        # 1. Project flux signal onto global time axis and couple to qubit
        flux_samples_global = self.flux_signal.samples_on(t_global)
        flux_global = FluxSignal(
            type=8, t_list=t_global, signal=flux_samples_global,
            trigger=0.0,
        )
        self.qubit.qubit_in_mag(
            flux_global, frame=1, omega_d=self.omega_d
        )

        psi_e = basis(self.qubit.n_levels, 1)
        p_e_list = np.zeros(len(self.tau_list))

        for i, tau in enumerate(self.tau_list):
            # Build Ramsey pulse sequence: pi/2 - tau - pi/2.
            # Each sub-pulse carries its own global trigger, so we no
            # longer need the ``ctrl.t_list -= ...`` offset hack.
            ctrl = create_ramsey_pulse(
                self.t_rabi, tau,
                omega_d=self.omega_d,
                phase1=self.phase1,
                phase2=self.phase2,
                qubit=self.qubit,
            )

            H = (
                QobjEvo(
                    self.qubit.H_list,
                    tlist=t_global,
                    order=1,
                )
                + QobjEvo(
                    ctrl.hamiltonian_on(t_global),
                    tlist=t_global,
                    order=1,
                )
            )
            # max_step prevents adaptive stepper from skipping over
            # narrow pi/2 pulse windows on the long t_global axis.
            result = mesolve(
                H, self.qubit.state, t_global, [],
                e_ops=[psi_e * psi_e.dag()],
                options={"max_step": float(CONFIG.awg.dt)},
            )
            p_e_list[i] = result.expect[0][-1]

        return ExperimentResult(
            data={
                "p_e": p_e_list,
                "flux_samples": self.flux_signal.signal.copy(),
            },
            axes={
                "tau": self.tau_list.copy(),
                "t_flux": self.flux_signal.t_list.copy(),
            },
            metadata={
                "experiment": "RamseyExperiment",
                "omega_d": self.omega_d,
            },
            config={
                "phase1": self.phase1,
                "phase2": self.phase2,
                "t_rabi": self.t_rabi.copy(),
                "t_global": self.t_global.copy(),
            },
        )
