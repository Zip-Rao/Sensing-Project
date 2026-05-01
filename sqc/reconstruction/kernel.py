"""sqc.reconstruction.kernel — KernelEstimator.

Estimates the control kernel for deconvolution-based reconstruction.
Phase 1: verbatim port from Pulse.get_kernel + CompositePulse.get_kernel.
Phase 2: refactored into a standalone class.

See _refactor_plan.md §7.3.
"""
from __future__ import annotations

import numpy as np
from qutip import sigmaz, basis, QobjEvo, mesolve


class KernelEstimator:
    """Estimate control kernel by perturbing the qubit with narrow stimuli.

    This is a Phase 1 verbatim port of Pulse.get_kernel() and
    CompositePulse.get_kernel(). Phase 2 will refactor into a clean
    standalone estimator.

    Parameters
    ----------
    stimulus_amplitude : float
        Amplitude of the probe stimulus.
    stimulus_width : float
        Width of the Gaussian probe (ns).
    """

    def __init__(
        self,
        stimulus_amplitude: float = 0.0215,
        stimulus_width: float = 3.0,
    ) -> None:
        self.stimulus_amplitude = stimulus_amplitude
        self.stimulus_width = stimulus_width

    def estimate(
        self, pulse, qubit
    ) -> tuple[np.ndarray, np.ndarray]:
        """Estimate kernel for a single pulse.

        Parameters
        ----------
        pulse : Pulse-like
            Control pulse with .hamiltonian, .t_list, .frame attributes.
        qubit : TransmonQubit
            Qubit used for simulation.

        Returns
        -------
        t_samples : np.ndarray
            Time points where kernel was evaluated.
        kernel : np.ndarray
            Kernel values.
        """
        raise NotImplementedError("Phase 2 will refactor from Pulse.get_kernel")
