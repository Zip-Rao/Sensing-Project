"""Compatibility facade for legacy analysis imports."""
from __future__ import annotations

import numpy as np

from sqc.reconstruction.basis import (
    R,
    basis_function_decomposition,
    generate_basis_functions,
)
from sqc.reconstruction.cryoscope import build_calibration_curve
from sqc.reconstruction.numerical_inverse import (
    compute_jacobian,
    compute_jacobian_finite_difference,
    forward_simulation,
    levenberg_marquardt,
)
from sqc.simulation.result import ExperimentResult, extract_expectation, extract_population


class Analysis:
    """Legacy facade over sqc.reconstruction and sqc.simulation helpers."""

    def get_expectation_values(self, result, e_ops_index):
        """Extract expectation values from a QuTiP result."""
        return extract_expectation(result, e_ops_index)

    def get_population(self, result, level):
        """Extract a level population trace from a QuTiP result."""
        return extract_population(result, level)

    def get_signal_from_ramsey_by_iq(self, qubit, tau_list, p_e_list_I, p_e_list_Q):
        """Delegate Ramsey IQ reconstruction."""
        from sqc.reconstruction.wiener import RamseyIQReconstruction

        spec = qubit.spec() if hasattr(qubit, "spec") else qubit
        measurement = ExperimentResult(
            data={
                "p_e_I": np.asarray(p_e_list_I),
                "p_e_Q": np.asarray(p_e_list_Q),
            },
            axes={"tau": np.asarray(tau_list)},
        )
        return RamseyIQReconstruction(qubit=spec).reconstruct(measurement)

    def get_signal_from_ramsey_by_unwrap(
        self, qubit, tau_list, p_e_list, k_span=3
    ):
        """Delegate single-channel Ramsey unwrap reconstruction."""
        from sqc.reconstruction.wiener import RamseyUnwrapReconstruction

        spec = qubit.spec() if hasattr(qubit, "spec") else qubit
        measurement = ExperimentResult(
            data={"p_e": np.asarray(p_e_list)},
            axes={"tau": np.asarray(tau_list)},
        )
        return RamseyUnwrapReconstruction(qubit=spec, k_span=k_span).reconstruct(
            measurement
        )

    def get_signal_from_diff_echo(self, qubit, p_e_list, t_int, k):
        """Delegate differential echo reconstruction."""
        from sqc.reconstruction.wiener import DiffEchoReconstruction

        spec = qubit.spec() if hasattr(qubit, "spec") else qubit
        measurement = ExperimentResult(data={"p_e": np.asarray(p_e_list)})
        return DiffEchoReconstruction(qubit=spec, t_int=t_int, k=k).reconstruct(
            measurement
        )

    def get_kernel(self, control_pulse, qubit):
        """Delegate kernel estimation."""
        from sqc.reconstruction.kernel import KernelEstimator

        return KernelEstimator().estimate(control_pulse, qubit)

    def wiener_deconvolution(self, delta_p, kernel, dt, lambdas):
        """Delegate Wiener deconvolution and return legacy tuple."""
        from sqc.reconstruction.wiener import WienerReconstruction

        measurement = ExperimentResult(
            data={"delta_p": np.asarray(delta_p)},
            axes={"scan": np.arange(len(delta_p)) * dt},
        )
        signal = WienerReconstruction(lambda_reg=lambdas).reconstruct(
            measurement, np.asarray(kernel), dt=dt
        )
        return signal.t_list, signal.signal

    def hammerstein_wiener_deconvolution(self, qubit, delta_p, kernel, dt, lambdas):
        """Delegate Hammerstein-Wiener deconvolution."""
        from sqc.reconstruction.hammerstein import HammersteinWienerReconstruction

        measurement = ExperimentResult(
            data={"delta_p": np.asarray(delta_p)},
            axes={"scan": np.arange(len(delta_p)) * dt},
        )
        return HammersteinWienerReconstruction(
            qubit=qubit, lambda_reg=lambdas
        ).reconstruct(measurement, np.asarray(kernel), dt=dt)

    def numerical_inverse(
        self,
        qubit,
        control_pulse,
        p_meas,
        t_list,
        B_guess,
        basis_type="fourier",
        n_basis=100,
        lambdas=100.0,
        max_iter=10,
        tol=1e-6,
    ):
        """Delegate LM reconstruction once Track B's fixed implementation lands."""
        from sqc.reconstruction.numerical_inverse import LMReconstruction

        measurement = ExperimentResult(
            data={"p_meas": np.asarray(p_meas)},
            axes={"t_signal": np.asarray(t_list)},
        )
        return LMReconstruction(
            qubit=qubit,
            control_pulse=control_pulse,
            basis_type=basis_type,
            n_basis=n_basis,
            lambda_reg=lambdas,
            max_iter=max_iter,
            tol=tol,
        ).reconstruct(measurement, initial_guess=np.asarray(B_guess))

    def get_h_from_phi(self, h_list, phi_list):
        """Build cryoscope calibration interpolation callables."""
        return build_calibration_curve(h_list, phi_list)

    def build_calibration_curve(self, h_list, phi_cal):
        """Alias for new cryoscope calibration-curve helper."""
        return build_calibration_curve(h_list, phi_cal)

    def get_signal_from_cryoscope(self, qubit, trunc_list, varphi_meas, h_of_phi, tau, dt):
        """Legacy cryoscope reconstruction by calibration inverse."""
        delta_phi = np.diff(varphi_meas, prepend=0)
        delta_phi_norm = delta_phi * tau / dt
        return trunc_list, h_of_phi(delta_phi_norm)

    def reconstruct_waveform_cryoscope(
        self,
        qubit,
        tau_list,
        phi_meas,
        h_of_phi=None,
        T_bias=None,
        method="sg",
        sg_window=5,
        sg_poly=2,
    ):
        """Reconstruct a waveform from cryoscope phase data."""
        from scipy.signal import savgol_filter

        tau = np.asarray(tau_list, dtype=float)
        phi = np.asarray(phi_meas, dtype=float)
        dt = tau[1] - tau[0]

        if method == "calibration" and h_of_phi is not None and T_bias is not None:
            delta_phi = np.gradient(phi, tau) * T_bias
            return tau, h_of_phi(delta_phi)

        if method == "sg":
            if sg_window >= len(tau):
                sg_window = max(3, len(tau) // 2 * 2 - 1)
            dphi_dt = savgol_filter(phi, sg_window, sg_poly, deriv=1, delta=dt)
        elif method == "diff":
            dphi_dt = np.gradient(phi, tau)
        else:
            raise ValueError(f"Unknown method: {method}")

        delta_f = dphi_dt / (2 * np.pi)
        kappa = qubit.frequency_sensitivity(qubit.flux)
        if abs(kappa) > 1e-6:
            h_recon = delta_f / kappa
        else:
            f_q = qubit.frequency - delta_f
            ratio = (f_q + qubit.EC) ** 2 / (8 * qubit.EC * qubit.EJ_0)
            h_recon = np.arccos(np.clip(ratio, 0, 1)) / np.pi
        return tau, h_recon


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
