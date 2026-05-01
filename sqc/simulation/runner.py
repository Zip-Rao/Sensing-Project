"""Simulation runners for QuTiP-based experiments."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from qutip import QobjEvo, basis, expect, mesolve, qeye

from sqc.simulation.result import ExperimentResult


class RunnerBase(ABC):
    """Base class for simulation runners."""

    @abstractmethod
    def run(self, *args, **kwargs) -> ExperimentResult:
        """Run a simulation and return structured results."""


@dataclass
class MesolveRunner(RunnerBase):
    """Thin ExperimentResult wrapper around :func:`qutip.mesolve`."""

    options: dict | None = None

    def run(
        self,
        hamiltonian,
        psi0,
        t_list,
        c_ops=None,
        e_ops=None,
        store_states: bool = False,
        metadata: dict | None = None,
    ) -> ExperimentResult:
        """Run mesolve and pack expectation values into ExperimentResult."""
        t_array = np.asarray(t_list, dtype=float)
        opts = dict(self.options or {})
        opts.setdefault("store_states", store_states)
        c_ops = [] if c_ops is None else c_ops
        e_ops = [] if e_ops is None else e_ops

        if isinstance(hamiltonian, QobjEvo):
            h = hamiltonian
        else:
            h = QobjEvo(hamiltonian, tlist=t_array, order=1)

        result = mesolve(h, psi0, t_array, c_ops, e_ops, options=opts)
        md = {"runner": "MesolveRunner", "raw_result": result}
        if metadata:
            md.update(metadata)
        return ExperimentResult(
            data={"expect": np.asarray(result.expect)},
            axes={"t": t_array},
            metadata=md,
        )


@dataclass
class SlidingMeasurementRunner(RunnerBase):
    """Slide a control pulse over a flux signal and measure final p_e."""

    options: dict | None = None

    def run(
        self,
        qubit,
        flux_signal,
        control_pulse,
        scan_list: np.ndarray | None = None,
    ) -> ExperimentResult:
        """Return p_e for each pulse/flux delay."""
        if scan_list is None:
            delay_start = flux_signal.t_list[0] - 0.5 * control_pulse.t_list[-1]
            delay_end = flux_signal.t_list[-1] + 0.5 * control_pulse.t_list[-1]
            n_samples = len(flux_signal.t_list) + len(control_pulse.t_list) - 1
            scan_list = np.linspace(delay_start, delay_end, n_samples)
        else:
            scan_list = np.asarray(scan_list, dtype=float)

        qubit_t = qubit.qubit_under_mag(flux_signal)
        p_e = []
        for t_delay in scan_list:
            h_total, t_evolve = self._build_delayed_hamiltonian(
                qubit, flux_signal, control_pulse, qubit_t, t_delay
            )
            p_e.append(
                self.single_measurement(
                    qubit,
                    flux_signal,
                    control_pulse,
                    t_delay,
                    qubit_t=qubit_t,
                    hamiltonian=h_total,
                    t_evolve=t_evolve,
                )
            )

        return ExperimentResult(
            data={"p_e": np.asarray(p_e)},
            axes={"scan": scan_list},
            metadata={"runner": "SlidingMeasurementRunner"},
        )

    def single_measurement(
        self,
        qubit,
        flux_signal,
        control_pulse,
        t_delay,
        qubit_t=None,
        hamiltonian=None,
        t_evolve=None,
    ) -> float:
        """Run one delayed measurement and return final excited population."""
        if qubit_t is None:
            qubit_t = qubit.qubit_under_mag(flux_signal)

        if hamiltonian is None:
            hamiltonian, t_evolve = self._build_delayed_hamiltonian(
                qubit, flux_signal, control_pulse, qubit_t, t_delay
            )

        h_total = QobjEvo(hamiltonian, tlist=t_evolve, order=1)
        dt_coeff = t_evolve[1] - t_evolve[0]
        options = {
            "atol": 1e-10,
            "rtol": 1e-8,
            "max_step": dt_coeff / 2,
            "nsteps": 10000,
        }
        options.update(self.options or {})
        result = mesolve(h_total, qubit.state, t_evolve, [], e_ops=[], options=options)
        final_state = result.states[-1]
        projector = basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()
        return float(expect(projector, final_state))

    def _build_delayed_hamiltonian(
        self,
        qubit,
        flux_signal,
        control_pulse,
        qubit_t,
        t_delay: float,
    ) -> tuple[list, np.ndarray]:
        """Build the delayed Hamiltonian used by legacy sliding measurement."""
        delta = t_delay - 0.5 * control_pulse.t_list[-1]
        t_start = min(delta, 0)
        t_end = max(flux_signal.t_list[-1], t_delay + 0.5 * control_pulse.t_list[-1])
        n_points = len(flux_signal.t_list) + len(control_pulse.t_list) - 1
        t_evolve = np.linspace(t_start, t_end, n_points)

        freq_coeffs = np.zeros(n_points)
        for i, t in enumerate(t_evolve):
            if 0 <= t <= flux_signal.t_list[-1]:
                index = np.searchsorted(flux_signal.t_list, t)
                index = min(index, len(qubit_t) - 1)
                qubit_current = qubit_t[index]
            else:
                qubit_current = qubit

            if control_pulse.frame == 0:
                freq_coeffs[i] = qubit_current.frequency
            else:
                freq_coeffs[i] = qubit_current.frequency - control_pulse.omega_d

        h_list = []
        identity_coeff = np.ones(n_points, dtype=complex)
        anharm = qubit.anharmonicity * 0.5 * (qubit.n * qubit.n - qubit.n)
        h_list.append([anharm, identity_coeff])
        if control_pulse.frame == 0:
            h_list.append([qubit.n + 0.5 * qeye(qubit.n_levels), freq_coeffs])
        else:
            h_list.append([qubit.n, freq_coeffs])

        t_pulse = control_pulse.t_list
        for op, coeffs in control_pulse.hamiltonian:
            coeff_global = np.zeros(n_points, dtype=complex)
            for i, t in enumerate(t_evolve):
                t_loc = t - delta
                if 0 <= t_loc <= t_pulse[-1]:
                    index = np.clip(
                        np.searchsorted(t_pulse, t_loc), 0, len(t_pulse) - 1
                    )
                    index = min(index, len(coeffs) - 1)
                    coeff_global[i] = coeffs[index]
            h_list.append([op, coeff_global])

        return h_list, t_evolve
