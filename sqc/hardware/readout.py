"""sqc.hardware.readout — Readout models.

IdealProjectiveReadout: projective measurement onto |1>.
IQReadoutModel: IQ readout via two Ramsey sequences with pi/2 phase offset.

Replaces src/protocal.py:IQ_readout function.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from sqc.config import CONFIG
from qutip import Qobj, QobjEvo, basis, mesolve


# ---------------------------------------------------------------------------
# ReadoutModel ABC
# ---------------------------------------------------------------------------

class ReadoutModel(ABC):
    """Abstract readout model.

    Subclasses implement measurement of qubit state.
    """

    @abstractmethod
    def measure(self, *args, **kwargs) -> dict[str, float]:
        """Perform readout and return measurement results."""
        ...


# ---------------------------------------------------------------------------
# IdealProjectiveReadout
# ---------------------------------------------------------------------------

@dataclass
class IdealProjectiveReadout(ReadoutModel):
    """Projective measurement onto |1>.

    Returns p_e = <1|rho|1>.
    """

    def measure(self, state, qubit=None) -> dict[str, float]:
        """Measure excited-state probability.

        Parameters
        ----------
        state : Qobj
            Quantum state (ket or density matrix).
        qubit : TransmonQubit or None
            Qubit (for n_levels if state dimension ambiguous).

        Returns
        -------
        dict[str, float]
            {"p_e": float}
        """
        n_levels = state.dims[0][0]
        psi_e = basis(n_levels, 1)
        proj = psi_e * psi_e.dag()
        from qutip import expect

        return {"p_e": float(expect(proj, state).real)}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _resample_hamiltonian(hamiltonian, src_tlist, dst_tlist):
    """Resample a list-format Hamiltonian onto a new time grid.

    Parameters
    ----------
    hamiltonian : list
        QuTiP list-format Hamiltonian: ``[H0, [H1, c1], [H2, c2], ...]``.
    src_tlist : array-like
        Original time points matching the coefficient arrays.
    dst_tlist : array-like
        Target time grid to interpolate onto.

    Returns
    -------
    list
        New list-format Hamiltonian with coefficients resampled onto
        ``dst_tlist``.
    """
    import numpy as np
    src = np.asarray(src_tlist, dtype=float)
    dst = np.asarray(dst_tlist, dtype=float)
    out = []
    for term in hamiltonian:
        if isinstance(term, list) and len(term) == 2:
            op, coeffs_src = term
            c_src = np.asarray(coeffs_src, dtype=complex)
            if c_src.ndim == 0 or len(c_src) == 1:
                c_dst = np.full(len(dst), complex(c_src))
            else:
                c_dst = np.interp(dst, src, c_src.real).astype(complex)
            out.append([op, c_dst])
        else:
            out.append(term)
    return out


# ---------------------------------------------------------------------------
# IQReadoutModel
# ---------------------------------------------------------------------------

@dataclass
class IQReadoutModel(ReadoutModel):
    """IQ readout via two Ramsey sequences with pi/2 phase offset.

    Replaces src/protocal.py:IQ_readout(type=2 or 3).

    Parameters
    ----------
    tau : float
        Free precession time (ns). Default 20.
    t_rabi : np.ndarray or None
        Rabi pulse time axis (ns). Default linspace(0, 10, 20).
    omega_d : float or None
        Drive frequency. If None, uses qubit frequency.
    """

    tau: float = 20.0
    t_rabi: np.ndarray | None = None
    omega_d: float | None = None

    def __post_init__(self):
        if self.t_rabi is None:
            self.t_rabi = CONFIG.pulse.t_rabi.copy()

    def measure(self, qubit, **extra) -> dict[str, float]:
        """Run two Ramsey sequences and return I, Q components.

        Requires qubit to have H_list and mag_signal populated
        (call qubit.qubit_in_mag(...) first).

        Parameters
        ----------
        qubit : TransmonQubit
            Qubit with H_list pre-computed via qubit_in_mag().
        **extra
            Override tau, t_rabi, or omega_d.

        Returns
        -------
        dict[str, float]
            {"p_e_I": float, "p_e_Q": float}
        """
        tau = extra.get("tau", self.tau)
        t_rabi = extra.get("t_rabi", self.t_rabi)
        omega_d = extra.get("omega_d", self.omega_d)
        if omega_d is None:
            omega_d = qubit.frequency

        from sqc.control.sequence import create_ramsey_pulse

        t_global = CONFIG.pulse.t_global

        ctrl_I = create_ramsey_pulse(
            t_rabi, tau, omega_d=omega_d,
            phase1=np.pi / 2, phase2=0.0,
            qubit=qubit,
        )
        ctrl_Q = create_ramsey_pulse(
            t_rabi, tau, omega_d=omega_d,
            phase1=np.pi / 2, phase2=np.pi / 2,
            qubit=qubit,
        )

        # Project control Hamiltonians onto the unified global time axis.
        H_I_on = ctrl_I.hamiltonian_on(t_global)
        H_Q_on = ctrl_Q.hamiltonian_on(t_global)

        # Resample qubit H_list onto the global time grid.  (After the
        # experiment layer is migrated in P7.3, qubit.H_list will already
        # be defined on t_global and this becomes a trivial identity.)
        t_sig = np.asarray(qubit.mag_signal.t_list, dtype=float)
        H_q_resampled = _resample_hamiltonian(qubit.H_list, t_sig, t_global)

        H_I = (
            QobjEvo(H_q_resampled, tlist=t_global, order=1)
            + QobjEvo(H_I_on, tlist=t_global, order=1)
        )
        H_Q = (
            QobjEvo(H_q_resampled, tlist=t_global, order=1)
            + QobjEvo(H_Q_on, tlist=t_global, order=1)
        )

        psi_e = basis(qubit.n_levels, 1)
        e_ops = [psi_e * psi_e.dag()]

        # max_step prevents the adaptive stepper from skipping over pulse
        # edges (where Hamiltonian coefficients change sharply).
        _dt = float(CONFIG.awg.dt)
        result_I = mesolve(
            H_I, qubit.state, t_global, [], e_ops=e_ops,
            options={"store_states": True, "max_step": _dt},
        )
        result_Q = mesolve(
            H_Q, qubit.state, t_global, [], e_ops=e_ops,
            options={"store_states": True, "max_step": _dt},
        )

        return {
            "p_e_I": float(result_I.expect[0][-1]),
            "p_e_Q": float(result_Q.expect[0][-1]),
        }


# ---------------------------------------------------------------------------
# Legacy backward-compat function
# ---------------------------------------------------------------------------

def IQ_readout_legacy(qubit, type, **kwargs):
    """Drop-in replacement for src/protocal.py:IQ_readout.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit with H_list pre-computed via qubit_in_mag().
    type : int
        Readout type. 2=cryoscope calib, 3=cryoscope measurement.
    **kwargs
        tau, t_rabi, h (for type 2).

    Returns
    -------
    tuple[float, float]
        (p_e_I, p_e_Q)
    """
    if type in (2, 3):
        readout = IQReadoutModel(
            tau=kwargs.get("tau", 20.0),
            t_rabi=kwargs.get("t_rabi", CONFIG.pulse.t_rabi.copy()),
        )
        result = readout.measure(qubit)
        return result["p_e_I"], result["p_e_Q"]
    else:
        raise NotImplementedError(f"IQ_readout type={type} not supported")
