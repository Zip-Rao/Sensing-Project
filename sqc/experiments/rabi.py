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
from sqc.control.sequence import create_pulse


@dataclass
class RabiExperiment(Experiment):
    """Rabi oscillation experiment.

    Applies a constant-amplitude drive pulse and measures p_e vs time.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit (src or sqc version).
    t_rabi : np.ndarray or None
        Rabi time axis (ns). Default linspace(0, 40, 1000).
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
        """Build Rabi pulse (QobjEvo)."""
        return create_pulse(
            self.qubit, frame=1, type=1,
            t_list=self.t_rabi, omega_d=self.omega_d,
            phase=0.0,
        )

    def run(self):
        """Execute Rabi experiment.

        Returns
        -------
        qutip.Result
            Raw mesolve result (matches legacy case 0 return type).
        """
        psi_e = basis(self.qubit.n_levels, 1)
        H_0 = self.qubit.get_hamiltonian_rwa(self.omega_d)
        H_pulse = self.build_sequence()
        H_rabi = QobjEvo(H_0) + QobjEvo(H_pulse, tlist=self.t_rabi)
        result = mesolve(
            H_rabi, self.qubit.state, self.t_rabi, [],
            e_ops=[psi_e * psi_e.dag()],
        )
        return result
