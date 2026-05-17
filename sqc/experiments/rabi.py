"""sqc.experiments.rabi — RabiExperiment.

Rabi oscillation measurement: single constant-amplitude pulse, sweep duration.
Replaces Protocal.evolve case 0.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.experiments.base import Experiment
from sqc.control.pulse import Pulse
from sqc.control.flux_signal import FluxSignal


@dataclass
class RabiExperiment(Experiment):
    """Rabi oscillation experiment.

    Applies a constant-amplitude drive pulse and measures p_e vs time.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit (src or sqc version).
    t_rabi : np.ndarray or None
        Rabi time axis (ns). Default arange(0, 40, dt).
    omega_d : float or None
        Drive frequency (rad*GHz). Default qubit.frequency.
    """

    qubit: object  # TransmonQubit (duck typed)
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.make_time(0, 40)
    )
    omega_d: float | None = None

    def __post_init__(self):
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency

    def build_sequence(self):
        """Build Rabi pulse as a Pulse object (with trigger=0)."""
        Omega = FluxSignal(
            type=1, t_list=self.t_rabi,
            amplitude=(np.pi / 2.0) / (self.t_rabi[-1] - self.t_rabi[0]),
        )
        return Pulse(
            frame=1, omega_d=self.omega_d, phase=0.0,
            Omega=Omega, is_rwa=True, qubit=self.qubit, trigger=0.0,
        )

    def run(self):
        """Execute Rabi experiment on unified global time axis.

        Returns
        -------
        qutip.Result
            Raw mesolve result (matches legacy case 0 return type).
        """
        t_global = CONFIG.pulse.t_global
        psi_e = basis(self.qubit.n_levels, 1)
        H_0 = self.qubit.get_hamiltonian_rwa(self.omega_d)
        rabi_pulse = self.build_sequence()
        H_rabi = (
            QobjEvo(H_0, tlist=t_global)
            + QobjEvo(rabi_pulse.hamiltonian_on(t_global),
                       tlist=t_global, order=1)
        )
        # max_step prevents adaptive stepper from skipping over narrow
        # pulse windows on the long t_global axis.
        result = mesolve(
            H_rabi, self.qubit.state, t_global, [],
            e_ops=[psi_e * psi_e.dag()],
            options={"max_step": float(CONFIG.awg.dt)},
        )
        return result
