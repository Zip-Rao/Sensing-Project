"""Control-kernel estimation utilities."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from qutip import QobjEvo, basis, mesolve, sigmaz

from sqc.control.flux_signal import FluxSignal


@dataclass
class KernelEstimator:
    """Estimate pulse response kernels by narrow Gaussian perturbations.

    The default path matches the legacy ``CompositePulse.get_kernel`` behavior:
    each sample is perturbed by a Gaussian flux signal, evolved, and normalized
    by the stimulus area.
    """

    stim_amplitude: float = 0.0215
    stim_width: float = 3.0
    auto_calibrate: bool = False

    def estimate(self, pulse: Any, qubit) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(t_samples, kernel)`` for a pulse or composite pulse."""
        t_list = np.asarray(pulse.t_list, dtype=float)
        psi_e = basis(qubit.n_levels, 1)

        if pulse.frame == 0:
            h_0 = QobjEvo(qubit.get_hamiltonian())
        else:
            h_0 = QobjEvo(qubit.get_hamiltonian_rwa(qubit.frequency))

        h_pulse = QobjEvo(pulse.hamiltonian, tlist=t_list, order=1)
        h_base = h_0 + h_pulse
        result_base = mesolve(
            h_base, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()]
        )
        p_e_base = result_base.expect[0][-1]

        kernel = np.zeros(len(t_list), dtype=float)
        for i, t_i in enumerate(t_list):
            stim = FluxSignal(
                type=3,
                t_list=t_list,
                amplitude=self._amplitude_for(qubit),
                center=t_i,
                width=self.stim_width,
            )
            stim_area = np.trapezoid(stim.signal, stim.t_list)
            qubit_t = qubit.qubit_under_mag(stim)
            h_stim = QobjEvo(
                qubit.qubit_under_mag_hamiltonian(
                    qubit_t, stim.t_list, pulse.frame, pulse.omega_d
                ),
                tlist=stim.t_list,
                order=1,
            )
            result = mesolve(
                h_base + h_stim,
                qubit.state,
                t_list,
                [],
                e_ops=[psi_e * psi_e.dag()],
            )
            kernel[i] = (result.expect[0][-1] - p_e_base) / stim_area

        return t_list.copy(), kernel

    def estimate_legacy_single_pulse(self, pulse: Any, qubit) -> tuple[np.ndarray, np.ndarray]:
        """Match the older ``Pulse.get_kernel`` sigma-z stimulus path."""
        t_list = np.asarray(pulse.t_list, dtype=float)
        psi_e = basis(qubit.n_levels, 1)

        if pulse.frame == 0:
            h_0 = QobjEvo(qubit.get_hamiltonian())
        else:
            h_0 = QobjEvo(qubit.get_hamiltonian_rwa(qubit.frequency))

        h_pulse = QobjEvo(pulse.hamiltonian, tlist=t_list, order=1)
        h_base = h_0 + h_pulse
        result_base = mesolve(
            h_base, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()]
        )
        p_e_base = result_base.expect[0][-1]

        kernel = np.zeros(len(t_list), dtype=float)
        for i, t_i in enumerate(t_list):
            stim = FluxSignal(
                type=3,
                t_list=t_list,
                amplitude=self.stim_amplitude,
                center=t_i,
                width=self.stim_width,
            )
            h_stim = QobjEvo([sigmaz(), stim.signal], tlist=stim.t_list, order=1)
            stim_area = (
                self.stim_amplitude * self.stim_width / 4.0 * np.sqrt(2 * np.pi)
            )
            result = mesolve(
                h_base + h_stim,
                qubit.state,
                t_list,
                [],
                e_ops=[psi_e * psi_e.dag()],
            )
            kernel[i] = (result.expect[0][-1] - p_e_base) / stim_area

        return t_list.copy(), kernel

    def _amplitude_for(self, qubit) -> float:
        """Return the perturbation amplitude for this qubit."""
        if not self.auto_calibrate:
            return self.stim_amplitude
        kappa = qubit.frequency_sensitivity(qubit.flux)
        if kappa == 0:
            return self.stim_amplitude
        target_shift = 0.01 * abs(qubit.anharmonicity)
        return float(target_shift / abs(kappa))
