"""Levenberg-Marquardt reconstruction interface.

The current repository marks the LM Jacobian convergence fix as unfinished in
``idea/_TODO_master.md``. This module provides the Phase 3 public interface and
keeps the implementation boundary explicit until the fixed Track B algorithm is
available for porting.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.reconstruction.base import Reconstruction
from sqc.reconstruction.basis import basis_function_decomposition


@dataclass
class LMReconstruction(Reconstruction):
    """Levenberg-Marquardt waveform reconstruction interface."""

    qubit: object
    control_pulse: object
    basis_type: Literal["bspline", "fourier", "legendre"] = "fourier"
    n_basis: int = 100
    lambda_reg: float = 100.0
    max_iter: int = 10
    tol: float = 1e-6
    mu_init: float = 1e-3
    use_adjoint: bool = True

    def reconstruct(
        self, measurement, kernel=None, calibration=None, initial_guess=None
    ) -> tuple[FluxSignal, dict]:
        """Run LM reconstruction once Track B's fixed implementation is present."""
        t_list = np.asarray(measurement.axes["t_signal"], dtype=float)
        if initial_guess is None:
            initial_guess = np.zeros_like(t_list)
        b_init_signal = FluxSignal(
            type=6,
            t_list=t_list,
            n_basis=self.n_basis,
            basis_type=self.basis_type,
        )
        b_init = basis_function_decomposition(
            initial_guess, t_list, b_init_signal.basis_functions
        )
        b_init_signal.update_signal(b=b_init)
        raise NotImplementedError(
            "LMReconstruction is blocked until Track B _TODO_master.md 0.3 "
            "provides the fixed Jacobian implementation and baseline."
        )


def forward_simulation(*args, **kwargs):
    """Compatibility placeholder for the old module-level LM helper."""
    raise NotImplementedError("forward_simulation awaits Track B LM fix")


def compute_jacobian(*args, **kwargs):
    """Compatibility placeholder for the old module-level LM helper."""
    raise NotImplementedError("compute_jacobian awaits Track B LM fix")


def compute_jacobian_finite_difference(*args, **kwargs):
    """Compatibility placeholder for the old module-level LM helper."""
    raise NotImplementedError("compute_jacobian_finite_difference awaits Track B LM fix")


def levenberg_marquardt(*args, **kwargs):
    """Compatibility placeholder for the old module-level LM helper."""
    raise NotImplementedError("levenberg_marquardt awaits Track B LM fix")
