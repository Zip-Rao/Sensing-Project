"""src_mirror.analysis — Facade over sqc.reconstruction.

Backward-compatible Analysis class. Same API as src/analysis.py,
but internal implementation delegates to sqc/ classes.

R1 compliance: NEVER modifies src/. This is a NEW file in src_mirror/.

Methods implemented in P3a:
  - get_expectation_values → sqc.simulation.result.extract_expectation
  - get_population          → sqc.simulation.result.extract_population
  - get_signal_from_ramsey_by_iq     → RamseyIQReconstruction
  - get_signal_from_ramsey_by_unwrap → RamseyUnwrapReconstruction
  - get_signal_from_diff_echo        → DiffEchoReconstruction
  - get_kernel              → KernelEstimator
  - wiener_deconvolution    → WienerReconstruction
  - hammerstein_wiener_deconvolution → HammersteinWienerReconstruction

Methods deferred to P3b/P3c (raise NotImplementedError):
  - numerical_inverse       (needs Track B 0.3)
  - get_h_from_phi          (needs Track B 1.1)
  - get_signal_from_cryoscope (needs Track B 1.1)
  - get_volterra_kernel     (stub, not yet implemented)

Module-level re-exports (for backward compat with src/analysis.py):
  - generate_basis_functions
  - basis_function_decomposition
  - R (alias for regularization_matrix)

See idea/refactor/_refactor_plan.md §8, phase_3_handbook.md §3.9.
"""
from __future__ import annotations

import numpy as np
from qutip import basis as _qbasis, expect as _qexpect

from sqc.reconstruction.basis import (
    generate_basis_functions,
    basis_function_decomposition,
    regularization_matrix,
    R,
)
from sqc.reconstruction.kernel import KernelEstimator
from sqc.reconstruction.wiener import (
    WienerReconstruction,
    RamseyIQReconstruction,
    RamseyUnwrapReconstruction,
    DiffEchoReconstruction,
)
from sqc.reconstruction.hammerstein import HammersteinWienerReconstruction
from sqc.simulation.result import (
    ExperimentResult,
    extract_expectation,
    extract_population,
)


# ---------------------------------------------------------------------------
# Module-level helpers (re-exported for backward compat)
# ---------------------------------------------------------------------------

__all__ = [
    "Analysis",
    "generate_basis_functions",
    "basis_function_decomposition",
    "R",
    "forward_simulation",
    "compute_jacobian",
    "compute_jacobian_finite_difference",
    "levenberg_marquardt",
]


def forward_simulation(*args, **kwargs):
    """Stub — requires Track B 0.3 (LM convergence fix)."""
    raise NotImplementedError(
        "forward_simulation: P3b — requires Track B 0.3 (LM convergence fix)"
    )


def compute_jacobian(*args, **kwargs):
    """Stub — requires Track B 0.3 (LM convergence fix)."""
    raise NotImplementedError(
        "compute_jacobian: P3b — requires Track B 0.3 (LM convergence fix)"
    )


def compute_jacobian_finite_difference(*args, **kwargs):
    """Stub — requires Track B 0.3 (LM convergence fix)."""
    raise NotImplementedError(
        "compute_jacobian_finite_difference: P3b — "
        "requires Track B 0.3 (LM convergence fix)"
    )


def levenberg_marquardt(*args, **kwargs):
    """Stub — requires Track B 0.3 (LM convergence fix)."""
    raise NotImplementedError(
        "levenberg_marquardt: P3b — requires Track B 0.3 (LM convergence fix)"
    )


# ---------------------------------------------------------------------------
# Analysis facade class
# ---------------------------------------------------------------------------

class Analysis:
    """Legacy facade over sqc.reconstruction.* classes.

    Same API as src/analysis.py:Analysis. Methods that have been
    internalised delegate to sqc/ classes; methods that depend on
    Track B raise NotImplementedError with a clear message.

    Usage (backward compatible)::

        ana = Analysis()
        B = ana.get_signal_from_ramsey_by_iq(qubit, tau_list, I, Q)
    """

    def __init__(self):
        pass

    # ------------------------------------------------------------------
    # Data extraction (no change from legacy)
    # ------------------------------------------------------------------

    def get_expectation_values(self, result, e_ops_index):
        """Extract expectation values from mesolve result.

        Delegates to sqc.simulation.result.extract_expectation.

        Parameters
        ----------
        result : qutip.Result
            QuTiP solver result.
        e_ops_index : int
            Index into result.expect.

        Returns
        -------
        np.ndarray
        """
        return extract_expectation(result, e_ops_index)

    def get_population(self, result, level):
        """Extract Fock-level population from mesolve result.

        Delegates to sqc.simulation.result.extract_population.

        Parameters
        ----------
        result : qutip.Result
            QuTiP solver result (store_states=True required).
        level : int
            Fock level to project onto.

        Returns
        -------
        np.ndarray
        """
        return extract_population(result, level)

    # ------------------------------------------------------------------
    # Ramsey protocol analysis
    # ------------------------------------------------------------------

    def get_signal_from_ramsey_by_iq(
        self, qubit, tau_list, p_e_list_I, p_e_list_Q
    ):
        """Reconstruct B(tau) from IQ-demodulated Ramsey data.

        Delegates to RamseyIQReconstruction.

        Parameters
        ----------
        qubit : TransmonQubit or QubitSpec
        tau_list : array-like
            Delay times.
        p_e_list_I : array-like
            I-channel excited-state probability.
        p_e_list_Q : array-like
            Q-channel excited-state probability.

        Returns
        -------
        np.ndarray
            Reconstructed B(tau).
        """
        recon = RamseyIQReconstruction(qubit=qubit)
        meas = ExperimentResult(
            data={
                "p_e_I": np.asarray(p_e_list_I, dtype=float),
                "p_e_Q": np.asarray(p_e_list_Q, dtype=float),
            },
            axes={"tau": np.asarray(tau_list, dtype=float)},
        )
        return recon.reconstruct(meas)

    def get_signal_from_ramsey_by_unwrap(
        self, qubit, tau_list, p_e_list, k_span=3
    ):
        """Reconstruct B(tau) via phase unwrapping.

        Delegates to RamseyUnwrapReconstruction.

        Parameters
        ----------
        qubit : TransmonQubit or QubitSpec
        tau_list : array-like
            Delay times.
        p_e_list : array-like
            Excited-state probability at each tau.
        k_span : int
            Branch-search window. Default 3.

        Returns
        -------
        np.ndarray
            Reconstructed B(tau).
        """
        recon = RamseyUnwrapReconstruction(qubit=qubit, k_span=k_span)
        meas = ExperimentResult(
            data={"p_e": np.asarray(p_e_list, dtype=float)},
            axes={"tau": np.asarray(tau_list, dtype=float)},
        )
        return recon.reconstruct(meas)

    # ------------------------------------------------------------------
    # Differential echo analysis
    # ------------------------------------------------------------------

    def get_signal_from_diff_echo(self, qubit, p_e_list, t_int, k):
        """Reconstruct B from differential echo p_e.

        Delegates to DiffEchoReconstruction.

        Parameters
        ----------
        qubit : TransmonQubit or QubitSpec
        p_e_list : array-like
            Excited-state probability at each scan point.
        t_int : float
            Integration time (ns).
        k : int
            Number of pi-pulse pairs.

        Returns
        -------
        np.ndarray
            Reconstructed B.
        """
        recon = DiffEchoReconstruction(qubit=qubit, t_int=t_int, k=k)
        meas = ExperimentResult(
            data={"p_e": np.asarray(p_e_list, dtype=float)},
        )
        return recon.reconstruct(meas)

    # ------------------------------------------------------------------
    # Kernel estimation
    # ------------------------------------------------------------------

    def get_kernel(self, control_pulse, qubit):
        """Estimate control kernel for a pulse.

        Delegates to KernelEstimator.

        Parameters
        ----------
        control_pulse : CompositePulse
            Control pulse.
        qubit : TransmonQubit

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (t_samples, kernel)
        """
        estimator = KernelEstimator()
        return estimator.estimate(control_pulse, qubit)

    def get_volterra_kernel(self, control_pulse, qubit):
        """Stub for second-order Volterra kernel.

        Not yet implemented in src/analysis.py either.
        """
        raise NotImplementedError(
            "get_volterra_kernel: not implemented (stub in src/analysis.py)"
        )

    # ------------------------------------------------------------------
    # Wiener / Hammerstein-Wiener deconvolution
    # ------------------------------------------------------------------

    def wiener_deconvolution(self, delta_p, kernel, dt, lambdas):
        """Linear Wiener deconvolution.

        Delegates to WienerReconstruction.

        Parameters
        ----------
        delta_p : np.ndarray
            Measured probability difference.
        kernel : np.ndarray
            Control kernel.
        dt : float
            Time step.
        lambdas : float
            Regularisation parameter.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (t_list, reconstructed_signal)
        """
        recon = WienerReconstruction(lambda_reg=lambdas)
        result = recon.reconstruct(delta_p, kernel, dt=dt)
        return result.t_list, result.signal

    def hammerstein_wiener_deconvolution(
        self, qubit, delta_p, kernel, dt, lambdas
    ):
        """Hammerstein-Wiener nonlinear deconvolution.

        Delegates to HammersteinWienerReconstruction.

        Parameters
        ----------
        qubit : TransmonQubit or QubitSpec
        delta_p : np.ndarray
            Measured probability difference.
        kernel : np.ndarray
            Control kernel.
        dt : float
            Time step.
        lambdas : float
            Regularisation parameter.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (B_list, B) — time axis and reconstructed B-field.
        """
        recon = HammersteinWienerReconstruction(
            qubit=qubit, lambda_reg=lambdas,
        )
        # Wiener step gives us the time axis; Hammerstein returns B
        wiener = WienerReconstruction(lambda_reg=lambdas)
        omega_signal = wiener.reconstruct(delta_p, kernel, dt=dt)
        omega_lists = omega_signal.t_list
        B = recon.reconstruct(delta_p, kernel, dt=dt)
        return omega_lists, B

    # ------------------------------------------------------------------
    # LM numerical inversion (P3b)
    # ------------------------------------------------------------------

    def numerical_inverse(
        self, qubit, control_pulse, p_meas, t_list, B_guess,
        basis_type="fourier", n_basis=100, lambdas=100.0,
        max_iter=10, tol=1e-6,
    ):
        """Levenberg-Marquardt numerical inversion.

        Not yet internalised — requires Track B 0.3 (LM convergence fix).

        Raises
        ------
        NotImplementedError
            Always in P3a.
        """
        raise NotImplementedError(
            "numerical_inverse: P3b — requires Track B 0.3 (LM convergence fix). "
            "Use src.analysis.Analysis.numerical_inverse for now."
        )

    # ------------------------------------------------------------------
    # Cryoscope (P3c)
    # ------------------------------------------------------------------

    def get_h_from_phi(self, h_list, phi_list):
        """Build phi(h) and h(phi) interpolation functions.

        Not yet internalised — requires Track B 1.1 (Cryoscope).

        Raises
        ------
        NotImplementedError
            Always in P3a.
        """
        raise NotImplementedError(
            "get_h_from_phi: P3c — requires Track B 1.1 (Cryoscope). "
            "Use src.analysis.Analysis.get_h_from_phi for now."
        )

    def get_signal_from_cryoscope(
        self, qubit, trunc_list, varphi_meas, h_of_phi, tau, dt
    ):
        """Cryoscope waveform reconstruction.

        Not yet internalised — requires Track B 1.1 (Cryoscope).

        Raises
        ------
        NotImplementedError
            Always in P3a.
        """
        raise NotImplementedError(
            "get_signal_from_cryoscope: P3c — requires Track B 1.1 (Cryoscope). "
            "Use src.analysis.Analysis.get_signal_from_cryoscope for now."
        )
