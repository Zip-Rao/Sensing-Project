"""Readout model abstractions and ideal software readout models."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

import numpy as np
from qutip import QobjEvo, basis, expect, mesolve

from sqc.control.sequence import create_ramsey_pulse


class ReadoutModel(ABC):
    """Base class for qubit readout models."""

    @abstractmethod
    def measure(self, *args, **kwargs):
        """Measure and return readout data."""


@dataclass
class IdealProjectiveReadout(ReadoutModel):
    """Projective readout onto the first excited state."""

    def measure(self, state, qubit) -> dict[str, float]:
        """Return excited-state probability for ``state``."""
        psi_e = basis(qubit.n_levels, 1)
        projector = psi_e * psi_e.dag()
        return {"p_e": float(expect(projector, state).real)}


@dataclass
class IQReadoutModel(ReadoutModel):
    """IQ readout via two Ramsey sequences with phase-offset analyzers."""

    tau: float = 20.0
    t_rabi: np.ndarray = field(default_factory=lambda: np.linspace(0, 10, 20))

    def measure(self, qubit, **extra) -> dict[str, float]:
        """Run I/Q Ramsey readouts for a qubit already coupled to flux."""
        tau = extra.get("tau", self.tau)
        t_rabi = extra.get("t_rabi", self.t_rabi)
        control_pulse_i = create_ramsey_pulse(
            t_rabi,
            tau,
            omega_d=qubit.frequency,
            phase1=np.pi / 2,
            phase2=0.0,
        )
        control_pulse_q = create_ramsey_pulse(
            t_rabi,
            tau,
            omega_d=qubit.frequency,
            phase1=np.pi / 2,
            phase2=np.pi / 2,
        )
        h_i = QobjEvo(
            control_pulse_i.hamiltonian, tlist=control_pulse_i.t_list, order=1
        ) + QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1)
        h_q = QobjEvo(
            control_pulse_q.hamiltonian, tlist=control_pulse_q.t_list, order=1
        ) + QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1)

        psi_e = basis(qubit.n_levels, 1)
        e_ops = [psi_e * psi_e.dag()]
        result_i = mesolve(
            h_i,
            qubit.state,
            control_pulse_i.t_list,
            [],
            e_ops=e_ops,
            options={"store_states": True},
        )
        result_q = mesolve(
            h_q,
            qubit.state,
            control_pulse_q.t_list,
            [],
            e_ops=e_ops,
            options={"store_states": True},
        )
        return {
            "p_e_I": float(result_i.expect[0][-1]),
            "p_e_Q": float(result_q.expect[0][-1]),
        }


def IQ_readout_legacy(qubit, type, **kwargs):
    """Drop-in helper for legacy ``src.protocal.IQ_readout`` type 2/3."""
    if type not in (2, 3):
        raise NotImplementedError(f"IQ_readout type={type} is not implemented")
    model = IQReadoutModel(
        tau=kwargs.get("tau", 20.0),
        t_rabi=kwargs.get("t_rabi", np.linspace(0, 10, 20)),
    )
    result = model.measure(qubit, **kwargs)
    return result["p_e_I"], result["p_e_Q"]
