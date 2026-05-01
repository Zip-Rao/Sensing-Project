"""sqc.simulation.runner — Simulation runner implementations.

RunnerBase: abstract interface for time evolution.
MesolveRunner: QuTiP mesolve-based solver.
SlidingMeasurementRunner: sliding measurement protocol.

Replaces Protocal.single_measurement and Protocal.sliding_measrement
(legacy spelling preserved in facade for backward compat).

See _refactor_plan.md §6.5.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from qutip import Qobj, QobjEvo, basis, expect, mesolve, qeye

from .result import ExperimentResult


# ---------------------------------------------------------------------------
# RunnerBase ABC
# ---------------------------------------------------------------------------

class RunnerBase(ABC):
    """Abstract base for all simulation runners.

    Subclasses implement run() to perform time evolution
    and return an ExperimentResult.
    """

    @abstractmethod
    def run(self, *args, **kwargs):
        """Execute the simulation and return results."""
        ...


# ---------------------------------------------------------------------------
# MesolveRunner
# ---------------------------------------------------------------------------

@dataclass
class MesolveRunner(RunnerBase):
    """QuTiP mesolve-based simulation runner.

    Parameters
    ----------
    options : dict or None
        Options dict passed to qutip.mesolve.
    """

    options: dict | None = None

    def run(
        self,
        H_list: list,
        psi0: Qobj,
        t_list: np.ndarray,
        c_ops: list[Qobj],
        e_ops: list[Qobj],
        store_states: bool = False,
    ) -> ExperimentResult:
        """Run mesolve and return ExperimentResult.

        Parameters
        ----------
        H_list : list
            QuTiP list-format Hamiltonian.
        psi0 : Qobj
            Initial state.
        t_list : np.ndarray
            Time axis.
        c_ops : list[Qobj]
            Collapse operators.
        e_ops : list[Qobj]
            Expectation operators.
        store_states : bool
            Whether to store intermediate states.

        Returns
        -------
        ExperimentResult
        """
        opts = dict(self.options or {})
        opts.setdefault("store_states", store_states)
        H = QobjEvo(H_list, tlist=t_list, order=1)
        result = mesolve(H, psi0, t_list, c_ops, e_ops, options=opts)
        return ExperimentResult(
            data={
                "expect": np.array(result.expect),
                "states": result.states if store_states else None,
            },
            axes={"t": np.asarray(t_list)},
            metadata={"runner": "MesolveRunner"},
        )


# ---------------------------------------------------------------------------
# SlidingMeasurementRunner
# ---------------------------------------------------------------------------

@dataclass
class SlidingMeasurementRunner(RunnerBase):
    """Slides a control pulse across a flux signal, measures p_e at each delay.

    Port of src/protocal.py:Protocal.sliding_measrement (lines 252-321).
    Note: legacy spelling "sliding_measrement" preserved in facade only;
    this class uses correct spelling "SlidingMeasurementRunner".

    Parameters
    ----------
    options : dict or None
        Options dict passed to qutip.mesolve.
    """

    options: dict | None = None

    def run(
        self,
        qubit,
        phi_signal,
        control_pulse,
        scan_list: np.ndarray | None = None,
    ) -> ExperimentResult:
        """Execute sliding measurement.

        Parameters
        ----------
        qubit : TransmonQubit
            Qubit object (src or sqc).
        phi_signal : Signal-like
            Flux signal with .t_list and .value_at().
        control_pulse : CompositePulse
            Control pulse object.
        scan_list : np.ndarray or None
            Pre-computed delay list. If None, auto-computed.

        Returns
        -------
        ExperimentResult
            With data["p_e"] and axes["scan"].
        """
        import time

        t_start = time.time()

        phi_tlist = np.asarray(phi_signal.t_list)
        pulse_tlist = np.asarray(control_pulse.t_list)

        # Compute scan list (delay axis)
        if scan_list is None:
            delay_start = phi_tlist[0] - 0.5 * pulse_tlist[-1]
            delay_end = phi_tlist[-1] + 0.5 * pulse_tlist[-1]
            n_samples = len(phi_tlist) + len(pulse_tlist) - 1
            scan_list = np.linspace(delay_start, delay_end, n_samples)

        # Pre-compute qubit_t array (one qubit per time point under flux)
        qubit_t = qubit.qubit_under_mag(phi_signal)

        p_e = []
        for t_delay in scan_list:
            H_total, t_evole = self._build_H_for_delay(
                qubit, phi_signal, control_pulse, qubit_t, t_delay
            )
            p_e.append(
                self._single_measurement(
                    qubit, phi_signal, control_pulse, t_delay,
                    qubit_t, H_total, t_evole,
                )
            )

        elapsed = time.time() - t_start
        print(
            f"Sliding measurement with {len(scan_list)} samples "
            f"took {elapsed:.2f} seconds"
        )

        return ExperimentResult(
            data={"p_e": np.array(p_e)},
            axes={"scan": np.array(scan_list)},
            metadata={
                "runner": "SlidingMeasurementRunner",
                "elapsed": elapsed,
            },
        )

    # -- single measurement ----------------------------------------------------

    def _single_measurement(
        self,
        qubit,
        phi_signal,
        control_pulse,
        t_delay: float,
        qubit_t=None,
        H=None,
        t_evole=None,
    ) -> float:
        """Perform a single measurement at a given delay.

        Port of src/protocal.py:Protocal.single_measurement (lines 187-250).

        Parameters
        ----------
        qubit : TransmonQubit
        phi_signal : Signal-like
        control_pulse : CompositePulse
        t_delay : float
            Time delay between signal and pulse.
        qubit_t : list or None
            Pre-computed qubit array.
        H : QobjEvo or None
            Pre-computed Hamiltonian.
        t_evole : np.ndarray or None
            Pre-computed evolution time list.

        Returns
        -------
        float
            Excited-state probability p_e.
        """
        import time

        n_levels = qubit.n_levels
        pulse_tlist = np.asarray(control_pulse.t_list)

        # Determine evolution time axis
        delta = t_delay - 0.5 * pulse_tlist[-1]
        if t_evole is None:
            t_start_ev = min(delta, 0.0)
            t_end_ev = max(
                phi_signal.t_list[-1],
                t_delay + 0.5 * pulse_tlist[-1],
            )
            N_ev = len(phi_signal.t_list) + len(pulse_tlist) - 1
            t_evole = np.linspace(t_start_ev, t_end_ev, N_ev)

        if H is None:
            H_pulse = lambda t: control_pulse.get_hamiltonian_at(t)

            def H_total(t, args):
                if t + delta < 0:
                    H_0 = (
                        qubit.get_hamiltonian()
                        if control_pulse.frame == 0
                        else qubit.get_hamiltonian_rwa(control_pulse.omega_d)
                    )
                elif 0 <= t + delta < phi_signal.t_list[-1]:
                    idx = np.searchsorted(phi_signal.t_list, t + delta)
                    qubit_current = qubit_t[idx]
                    H_0 = (
                        qubit_current.get_hamiltonian()
                        if control_pulse.frame == 0
                        else qubit_current.get_hamiltonian_rwa(
                            control_pulse.omega_d
                        )
                    )
                else:
                    H_0 = (
                        qubit.get_hamiltonian()
                        if control_pulse.frame == 0
                        else qubit.get_hamiltonian_rwa(control_pulse.omega_d)
                    )
                return H_0 + H_pulse(t)
        else:
            H_total = QobjEvo(H, tlist=t_evole, order=1)

        # mesolve options
        dt_coeff = t_evole[1] - t_evole[0]
        options = {
            "atol": 1e-10,
            "rtol": 1e-8,
            "max_step": dt_coeff / 2.0,
            "nsteps": 10000,
        }
        result = mesolve(
            H_total, qubit.state, t_evole, [], e_ops=[], options=options
        )
        final_state = result.states[-1]
        p_e_val = float(
            expect(
                basis(n_levels, 1) * basis(n_levels, 1).dag(),
                final_state,
            )
        )
        return p_e_val

    # -- build Hamiltonian for a given delay -----------------------------------

    def _build_H_for_delay(
        self,
        qubit,
        phi_signal,
        control_pulse,
        qubit_t,
        t_delay: float,
    ):
        """Build H_list and t_evole for a single delay.

        Port of the inner H(t_delay) function from
        src/protocal.py:Protocal.sliding_measrement (lines 267-309).
        """
        n_levels = qubit.n_levels
        n_op = qubit.n
        pulse_tlist = np.asarray(control_pulse.t_list)

        delta = t_delay - 0.5 * pulse_tlist[-1]
        t_start = min(delta, 0.0)
        t_end = max(
            phi_signal.t_list[-1],
            t_delay + 0.5 * pulse_tlist[-1],
        )
        N = len(phi_signal.t_list) + len(pulse_tlist) - 1
        t_evole = np.linspace(t_start, t_end, N)

        freq_coeffs = np.zeros(N)
        EJ_coeffs = np.zeros(N)

        for i, t in enumerate(t_evole):
            if 0 <= t <= phi_signal.t_list[-1]:
                idx = np.searchsorted(phi_signal.t_list, t)
                idx = min(idx, len(qubit_t) - 1)
                qubit_current = qubit_t[idx]
            else:
                qubit_current = qubit

            if control_pulse.frame == 0:
                freq_coeffs[i] = qubit_current.frequency
            else:
                freq_coeffs[i] = (
                    qubit_current.frequency - control_pulse.omega_d
                )
            EJ_coeffs[i] = -qubit_current.EJ + 0.25 * qubit_current.EC

        H_list = []
        id_arr = np.ones(N, dtype=complex)

        if control_pulse.frame == 0:
            H_list.append([
                qubit.anharmonicity * 0.5 * (n_op * n_op - n_op),
                id_arr,
            ])
            H_list.append([n_op + 0.5 * qeye(n_levels), freq_coeffs])
        else:
            H_list.append([
                qubit.anharmonicity * 0.5 * (n_op * n_op - n_op),
                id_arr,
            ])
            H_list.append([n_op, freq_coeffs])

        # Append pulse terms
        t_pulse = pulse_tlist
        for op, coeffs in control_pulse.hamiltonian:
            coeff_global = np.zeros(N, dtype=complex)
            for i, t in enumerate(t_evole):
                t_loc = t - delta
                if 0 <= t_loc <= t_pulse[-1]:
                    idx = int(np.clip(
                        np.searchsorted(t_pulse, t_loc), 0, len(t_pulse) - 1
                    ))
                    idx = min(idx, len(coeffs) - 1)
                    coeff_global[i] = coeffs[idx]
            H_list.append([op, coeff_global])

        return H_list, t_evole
