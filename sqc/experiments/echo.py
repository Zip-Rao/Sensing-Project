"""sqc.experiments.echo — DiffEchoExperiment.

Differential echo sensing protocol.
Replaces Protocal.evolve case 2.

Sequence: pi/2 - [tau - pi - tau' - (tau+t_int) - pi - tau'']^k - pi/2
Accumulates phase over k repetitions to enhance weak signal sensitivity.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.experiments.base import Experiment
from sqc.control.flux_signal import FluxSignal, CompositeSignal
from sqc.control.sequence import create_diff_echo_pulse
from sqc.simulation.result import ExperimentResult


@dataclass
class DiffEchoExperiment(Experiment):
    """Differential echo sensing experiment.

    Default parameters match src/protocal.py case 2 exactly.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit object (src or sqc version).
    flux_signal : FluxSignal or None
        Base flux signal. If None, creates default Gaussian.
    k : int
        Number of echo repetitions. Default 5.
    t_rabi : np.ndarray or None
        Pulse time axis (ns). Default linspace(0, 10, 20).
    t_int : float or None
        Interaction time (ns). Default computed from t_rabi.
    t_rep : float or None
        Repetition period (ns). Default computed from flux_signal.
    tau_list : np.ndarray or None
        Free evolution times (ns). Default from flux_signal.t_list.
    t_global : np.ndarray or None
        Global evolution time (ns).
    omega_d : float or None
        Drive frequency.
    """

    qubit: object  # TransmonQubit (duck typed)
    flux_signal: FluxSignal | None = None
    k: int = 5
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    t_int: float | None = None
    t_rep: float | None = None
    tau_list: np.ndarray | None = None
    t_global: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.make_time(-10, 1010)
    )
    omega_d: float | None = None

    # -- P9.B --
    control_line: object | None = None

    def __post_init__(self):
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency

        if self.flux_signal is None:
            # Default test signal matches src/protocal.py case 2
            t_list = CONFIG.pulse.make_time(0, 100)
            self.flux_signal = FluxSignal(
                type=3,
                t_list=t_list,
                amplitude=0.01,
                frequency=0.01,
                rise=10,
                fall=10,
                center=50,
                noise_level=0.0001,
            )

        if self.t_int is None:
            self.t_int = (
                self.t_rabi[-1] - self.t_rabi[0]
            ) * 0.5

        if self.t_rep is None:
            self.t_rep = (
                self.flux_signal.t_list[-1]
                - self.flux_signal.t_list[0]
            )

        if self.tau_list is None:
            self.tau_list = self.flux_signal.t_list.copy()

    def build_sequence(self):
        """Sequence is constructed per-tau in run()."""
        return None

    def run(self) -> ExperimentResult:
        """Execute differential echo experiment on unified global time axis.

        Returns
        -------
        ExperimentResult
        """
        t_global = self.t_global

        # Build composite flux signal: k*2 copies, then project onto t_global.
        flux_routed = self._route_flux(self.flux_signal)
        phi_list = [flux_routed.copy() for _ in range(2 * self.k)]
        composite_phi = CompositeSignal(phi_list)
        flux_global = FluxSignal(
            type=8, t_list=t_global,
            signal=composite_phi.samples_on(t_global),
            trigger=0.0,
        )
        self.qubit.qubit_in_mag(
            flux_global, frame=1, omega_d=self.omega_d
        )

        psi_e = basis(self.qubit.n_levels, 1)
        p_e_list = []

        for tau in self.tau_list:
            ctrl = create_diff_echo_pulse(
                self.t_rabi, tau, self.t_int, self.t_rep,
                k=self.k, omega_d=self.omega_d,
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
            # narrow pi/pi/2 pulse windows on the long t_global axis.
            result = mesolve(
                H, self.qubit.state, t_global, [],
                e_ops=[psi_e * psi_e.dag()],
                options={"max_step": float(CONFIG.awg.dt)},
            )
            p_e = result.expect[0][-1]
            p_e_list.append(p_e)

        return ExperimentResult(
            data={
                "p_e": np.array(p_e_list),
                "flux_samples": self.flux_signal.signal.copy(),
            },
            axes={
                "tau": self.tau_list.copy(),
                "t_flux": self.flux_signal.t_list.copy(),
            },
            metadata={
                "experiment": "DiffEchoExperiment",
                "k": self.k,
                "t_int": self.t_int,
                "omega_d": self.omega_d,
            },
            config={
                "k": self.k,
                "t_int": self.t_int,
                "t_rep": self.t_rep,
            },
        )
