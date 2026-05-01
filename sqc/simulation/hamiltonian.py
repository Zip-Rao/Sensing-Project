"""sqc.simulation.hamiltonian — HamiltonianBuilder.

Pure functional builder: takes device + flux + pulse specs, returns
QuTiP list-format Hamiltonian. Never mutates input objects.

See _refactor_plan.md §4.2 (design principle #3).
"""
from __future__ import annotations

from typing import Literal

import numpy as np
from qutip import qeye, num


class HamiltonianBuilder:
    """Build time-dependent Hamiltonian in QuTiP list format.

    Pure functional: never mutates qubit, flux_signal, or pulse.

    Usage
    -----
    H_list, t_global = HamiltonianBuilder.build(
        qubit=spec, flux_signal=sig, pulse=p, frame="rotating", omega_d=omega
    )
    """

    @staticmethod
    def build(
        qubit,
        flux_signal=None,
        pulse=None,
        frame: Literal["lab", "rotating"] = "rotating",
        omega_d: float | None = None,
    ) -> tuple[list, np.ndarray]:
        """Build H_list and global time axis for mesolve.

        Parameters
        ----------
        qubit : QubitSpec or TransmonQubit
            Qubit spec. If TransmonQubit, .spec() is called to get QubitSpec.
        flux_signal : Signal-like or None
            Flux signal with .t_list and .signal/.samples attributes.
            If None, returns static Hamiltonian.
        pulse : Pulse-like or None
            Control pulse (appended to H_list if given).
        frame : {"lab", "rotating"}
            Reference frame.
        omega_d : float or None
            Drive/reference frequency for rotating frame. Defaults to
            qubit frequency if None.

        Returns
        -------
        H_list : list
            QuTiP list-format Hamiltonian: [op_scalar, [op, coeff_array], ...]
        t_global : np.ndarray
            Global time axis (from flux_signal.t_list if given, else linspace).

        Notes
        -----
        The H_list format matches src/qubit.py:qubit_in_mag output:
            H_list[0] = (alpha/2) * (n^2 - n)   [scalar → constant coeff = 1]
            H_list[1] = [n + 0.5*I, freq_coeffs] or [n, freq_coeffs]
        Additional pulse terms from pulse.hamiltonian are appended.
        """
        # Resolve QubitSpec
        if hasattr(qubit, "spec"):
            spec = qubit.spec()
        else:
            spec = qubit

        n_lev = spec.n_levels
        n_op = num(n_lev)

        # Determine time axis and frequency coefficients
        if flux_signal is None:
            t_global = np.linspace(0, 100, 100)
            freq_coeffs = np.full(len(t_global), spec.frequency())
        else:
            t_global = np.asarray(flux_signal.t_list, dtype=float)
            samples = np.asarray(
                flux_signal.samples
                if hasattr(flux_signal, "samples")
                else flux_signal.signal
            ,
                dtype=float,
            )
            freq_coeffs = np.array(
                [spec.frequency(spec.flux_bias + s) for s in samples],
                dtype=float,
            )

        # Apply rotating frame shift
        if frame == "rotating":
            if omega_d is None:
                omega_d = spec.frequency()
            freq_coeffs = freq_coeffs - omega_d

        # Build H_list matching src/qubit.py:qubit_in_mag format
        anharmonicity = spec.anharmonicity()
        H_list = [
            anharmonicity * 0.5 * (n_op * n_op - n_op),
        ]

        # WHY: use n+0.5*I in both frames to match src/qubit.py:qubit_in_mag exactly.
        # The 0.5*I term is a constant energy shift — does not affect dynamics.
        H_list.append([n_op + 0.5 * qeye(n_lev), freq_coeffs])

        # Append pulse terms if given
        if pulse is not None:
            for op, coeff in pulse.hamiltonian:
                H_list.append([op, coeff])

        return H_list, t_global
