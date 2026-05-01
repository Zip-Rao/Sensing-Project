"""Hamiltonian builders for device + control objects."""
from __future__ import annotations

from typing import Literal

import numpy as np
from qutip import num, qeye


class HamiltonianBuilder:
    """Pure builder for QuTiP list-format Hamiltonians."""

    @staticmethod
    def build(
        qubit,
        flux_signal=None,
        pulse=None,
        frame: Literal["lab", "rotating"] = "rotating",
        omega_d: float | None = None,
    ) -> tuple[list, np.ndarray]:
        """Return (H_list, t_list) without mutating inputs."""
        spec = qubit.spec() if hasattr(qubit, "spec") else qubit
        n_levels = spec.n_levels
        n_op = num(n_levels)

        if flux_signal is None:
            t_list = np.linspace(0, 100, 100)
            samples = np.zeros_like(t_list)
        else:
            t_list = np.asarray(flux_signal.t_list)
            samples = np.asarray(
                flux_signal.samples if hasattr(flux_signal, "samples") else flux_signal.signal
            )

        freq_coeffs = np.zeros(len(t_list))
        for i, sample in enumerate(samples):
            flux = spec.flux_bias + sample
            frequency = spec.frequency(flux)
            if frame == "lab":
                freq_coeffs[i] = frequency
            else:
                if omega_d is None:
                    omega_d = spec.frequency()
                freq_coeffs[i] = frequency - omega_d

        h_list = [
            spec.anharmonicity() * 0.5 * (n_op * n_op - n_op),
            [n_op + 0.5 * qeye(n_levels), freq_coeffs],
        ]

        if pulse is not None:
            h_list.extend(pulse.hamiltonian)

        return h_list, t_list
