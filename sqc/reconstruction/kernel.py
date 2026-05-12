"""sqc.reconstruction.kernel — KernelEstimator.

Unified kernel estimation that replaces three duplicate implementations:
  - src/pulse.py:Pulse.get_kernel
  - src/pulse.py:CompositePulse.get_kernel
  - src/analysis.py:Analysis.get_kernel

Estimates the control kernel by perturbing with a narrow Gaussian stimulus
at each time point.

See _refactor_plan.md §7.3 and §14.1 (debt D2).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.control.pulse import PulseBase, CompositePulse


@dataclass
class KernelEstimator:
    """Estimate the control kernel of a pulse (or composite pulse).

    Perturbs the qubit with a narrow Gaussian stimulus at each time point
    of the pulse, measuring the change in final p_e.

    Parameters
    ----------
    stim_amplitude : float
        Amplitude of the probe stimulus (Phi_0).
        Default: ``CONFIG.reconstruction.stim_amplitude``.
    stim_width : float
        Width of the Gaussian probe (ns).
        Default: ``CONFIG.reconstruction.stim_width``.
    auto_calibrate : bool
        If True, adjust amplitude per qubit. Default False.
    """

    stim_amplitude: float = field(
        default_factory=lambda: CONFIG.reconstruction.stim_amplitude
    )
    stim_width: float = field(
        default_factory=lambda: CONFIG.reconstruction.stim_width
    )
    auto_calibrate: bool = False

    def estimate(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Estimate kernel for a pulse + qubit system.

        Parameters
        ----------
        pulse : Pulse, CompositePulse, or PulseBase
            Control pulse with .hamiltonian, .t_list, .frame, .omega_d.
        qubit : TransmonQubit
            Qubit for simulation (src or sqc version).
        t_samples : np.ndarray or None
            Time points to evaluate kernel at. If None, uses pulse.t_list.

        Returns
        -------
        t_samples : np.ndarray
            Time points where kernel was evaluated.
        kernel : np.ndarray
            Kernel values (dimensionless, d(p_e) / d(flux_integral)).
        """
        import time

        t_start = time.time()

        pulse_tlist = np.asarray(pulse.t_list)
        n_levels = qubit.n_levels
        psi_e = basis(n_levels, 1)

        # Build static H_0
        if hasattr(pulse, "frame") and pulse.frame == 0:
            H_0 = QobjEvo(qubit.get_hamiltonian())
        else:
            H_0 = QobjEvo(qubit.get_hamiltonian_rwa(
                qubit.frequency if hasattr(qubit, "frequency") else 0.0
            ))

        # Build H_pulse
        H_pulse = QobjEvo(pulse.hamiltonian, tlist=pulse_tlist, order=1)
        H_base = H_0 + H_pulse

        # Baseline result
        result_base = mesolve(
            H_base, qubit.state, pulse_tlist, [],
            e_ops=[psi_e * psi_e.dag()],
        )
        p_e_base = result_base.expect[0][-1]

        # Kernel evaluation
        if t_samples is None:
            t_samples = pulse_tlist.copy()

        frame_val = getattr(pulse, "frame", 1)
        omega_d_val = getattr(pulse, "omega_d", qubit.frequency)

        kernel = np.zeros(len(t_samples))
        for i, t_i in enumerate(t_samples):
            # Build narrow Gaussian stimulus at t_i
            stim = FluxSignal(
                type=3,
                t_list=pulse_tlist,
                amplitude=self._amplitude_for(qubit),
                center=t_i,
                width=self.stim_width,
            )
            stim_area = np.trapezoid(stim.signal, stim.t_list)

            # Per-timestep qubit array
            qubit_t = qubit.qubit_under_mag(stim)

            # Build stimulus Hamiltonian (CompositePulse-style, matches legacy)
            H_stim = QobjEvo(
                qubit.qubit_under_mag_hamiltonian(
                    qubit_t,
                    stim.t_list,
                    frame_val,
                    omega_d_val,
                ),
                tlist=stim.t_list,
                order=1,
            )

            H_total = H_0 + H_pulse + H_stim
            result_stim = mesolve(
                H_total, qubit.state, pulse_tlist, [],
                e_ops=[psi_e * psi_e.dag()],
            )
            p_e_stim = result_stim.expect[0][-1]
            kernel[i] = (p_e_stim - p_e_base) / stim_area

        elapsed = time.time() - t_start
        print(f"Kernel computation took {elapsed:.2f} seconds")

        return t_samples, kernel

    def _amplitude_for(self, qubit) -> float:
        """Determine stimulus amplitude for the given qubit.

        Parameters
        ----------
        qubit : TransmonQubit

        Returns
        -------
        float
            Stimulus amplitude.
        """
        if not self.auto_calibrate:
            return self.stim_amplitude
        # Auto-calibration: pick amplitude such that the stimulus produces
        # a frequency shift ~ 1% of qubit anharmonicity.
        try:
            kappa = qubit.frequency_sensitivity(
                qubit.flux if hasattr(qubit, "flux") else 0.0
            )
        except Exception:
            kappa = 0.0
        if kappa == 0.0:
            return self.stim_amplitude
        target_shift = 0.01 * abs(qubit.anharmonicity)
        return target_shift / abs(kappa)
