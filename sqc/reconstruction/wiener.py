"""sqc.reconstruction.wiener — linear reconstruction methods.

Port from src/analysis.py:
  - Analysis.wiener_deconvolution           → WienerReconstruction
  - Analysis.get_signal_from_ramsey_by_iq   → RamseyIQReconstruction
  - Analysis.get_signal_from_ramsey_by_unwrap → RamseyUnwrapReconstruction
  - Analysis.get_signal_from_diff_echo      → DiffEchoReconstruction
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.reconstruction.base import Reconstruction


def _get_sensitivity(qubit) -> float:
    """Extract frequency sensitivity kappa = dω/dΦ.

    Works with both sqc.devices.transmon.QubitSpec (P1+) and
    legacy src.qubit.TransmonQubit.
    """
    # QubitSpec (frozen, has flux_bias attribute)
    if hasattr(qubit, "flux_bias"):
        return qubit.sensitivity()

    # Legacy TransmonQubit (mutable, has flux attribute)
    if hasattr(qubit, "frequency_sensitivity") and hasattr(qubit, "flux"):
        return qubit.frequency_sensitivity(qubit.flux)

    raise TypeError(
        f"qubit {type(qubit).__name__} has neither flux_bias (QubitSpec) nor "
        f"frequency_sensitivity + flux (legacy TransmonQubit)"
    )


# ---------------------------------------------------------------------------
# Wiener deconvolution
# ---------------------------------------------------------------------------

@dataclass
class WienerReconstruction(Reconstruction):
    """Linear Wiener deconvolution.

    Replaces Analysis.wiener_deconvolution.

    Parameters
    ----------
    lambda_reg : float
        Regularisation parameter for the Wiener filter. Default 1.0.
    """

    lambda_reg: float = 1.0

    def reconstruct(
        self,
        measurement,
        kernel: np.ndarray,
        calibration=None,
        dt: float | None = None,
    ) -> FluxSignal:
        """Linear Wiener deconvolution.

        Parameters
        ----------
        measurement : ExperimentResult or ndarray
            If ExperimentResult: extracts "delta_p" from measurement.data.
            If ndarray: treated directly as delta_p.
        kernel : np.ndarray
            Control kernel (kernel_len).
        calibration : ignored
            Reserved for interface compatibility.
        dt : float or None
            Time step. If None, inferred from measurement axes.

        Returns
        -------
        FluxSignal
            Reconstructed signal with type=8 (custom), t_list and signal set.
        """
        # --- extract delta_p ---
        if hasattr(measurement, "data") and isinstance(measurement.data, dict):
            delta_p = np.asarray(measurement.data["delta_p"], dtype=float)
        else:
            delta_p = np.asarray(measurement, dtype=float)

        if dt is None:
            if hasattr(measurement, "axes") and "scan" in measurement.axes:
                scan = np.asarray(measurement.axes["scan"])
                dt = scan[1] - scan[0]
            else:
                dt = 1.0

        # --- Wiener deconvolution (verbatim from src/analysis.py:221-251) ---
        N_del = len(delta_p)
        N_ker = len(kernel)
        N = N_del - N_ker + 1

        N_fft = N_del

        p_pad = np.zeros(N_fft)
        p_pad[:N_del] = delta_p

        k_pad = np.zeros(N_fft)
        k_pad[:N_ker] = kernel

        Y = np.fft.fft(p_pad)
        H = np.fft.fft(k_pad)

        H_conj = np.conj(H)
        G = H_conj / (np.abs(H) ** 2 + self.lambda_reg**2)

        X_w = Y * G / dt
        x_rec = np.real(np.fft.ifft(X_w))

        x_final = x_rec[:N]
        x_lists = np.linspace(0, (N - 1) * dt, N)

        return FluxSignal(type=8, t_list=x_lists, signal=x_final)


# ---------------------------------------------------------------------------
# Ramsey IQ reconstruction
# ---------------------------------------------------------------------------

@dataclass
class RamseyIQReconstruction(Reconstruction):
    """Ramsey waveform reconstruction from I/Q demodulated data.

    Replaces Analysis.get_signal_from_ramsey_by_iq.

    p_e(τ) = 0.5*(1 - cos(φ(τ))); two measurements with π/2 phase offset
    directly give sin(φ) and cos(φ).

    Parameters
    ----------
    qubit
        QubitSpec or legacy TransmonQubit providing sensitivity().
    """

    qubit: object  # QubitSpec or legacy TransmonQubit

    def reconstruct(
        self,
        measurement,
        kernel=None,
        calibration=None,
        **kwargs,
    ) -> np.ndarray:
        """Reconstruct B(τ) from IQ Ramsey data.

        Parameters
        ----------
        measurement : ExperimentResult
            Must have data["p_e_I"], data["p_e_Q"] and axes["tau"].
        kernel, calibration : ignored
            Reserved for interface compatibility.

        Returns
        -------
        np.ndarray
            Reconstructed magnetic field B(tau) at each tau point.
        """
        p_e_list_I = np.asarray(measurement.data["p_e_I"], dtype=float)
        p_e_list_Q = np.asarray(measurement.data["p_e_Q"], dtype=float)
        tau_list = np.asarray(measurement.axes["tau"], dtype=float)

        C_I = np.max(p_e_list_I) - np.min(p_e_list_I)
        C_Q = np.max(p_e_list_Q) - np.min(p_e_list_Q)

        cosphi = (1 - 2 * p_e_list_I) / C_I
        sinphi = (1 - 2 * p_e_list_Q) / C_Q

        phi = np.unwrap(np.arctan2(sinphi, cosphi))
        B = np.gradient(phi, tau_list)
        kappa = _get_sensitivity(self.qubit)
        B /= kappa

        return B


# ---------------------------------------------------------------------------
# Ramsey phase-unwrapping reconstruction
# ---------------------------------------------------------------------------

@dataclass
class RamseyUnwrapReconstruction(Reconstruction):
    """Ramsey waveform reconstruction via phase unwrapping.

    Replaces Analysis.get_signal_from_ramsey_by_unwrap.

    From p_e(τ) = 0.5*(1 - cos(φ(τ))), uses φ(0)=0 prior and greedy
    tracking of the smoothest path through candidate branches
    {2kπ ± arccos(cosφ)}.

    Parameters
    ----------
    qubit
        QubitSpec or legacy TransmonQubit providing sensitivity().
    k_span : int
        Search window for integer branch k. Default 3.
    """

    qubit: object
    k_span: int = 3

    def reconstruct(
        self,
        measurement,
        kernel=None,
        calibration=None,
        **kwargs,
    ) -> np.ndarray:
        """Reconstruct B(τ) via phase unwrapping.

        Parameters
        ----------
        measurement : ExperimentResult
            Must have data["p_e"] and axes["tau"].
        kernel, calibration : ignored

        Returns
        -------
        np.ndarray
            Reconstructed magnetic field B(tau) at each tau point.
        """
        p_e = np.asarray(measurement.data["p_e"], dtype=float)
        tau = np.asarray(measurement.axes["tau"], dtype=float)

        N = len(p_e)

        # cos(φ) = 1 - 2*p_e
        cos_phi = 1 - 2 * p_e
        theta = np.arccos(np.clip(cos_phi, -1, 1))

        ks = np.arange(-self.k_span, self.k_span + 1)

        def _candidates(i):
            return np.array(
                [2 * k * np.pi + s * theta[i] for k in ks for s in [1, -1]]
            )

        phi = np.zeros(N)

        # Step 0: φ(0) = 0 (choose candidate closest to 0)
        cands = _candidates(0)
        phi[0] = cands[np.argmin(np.abs(cands))]

        # Step 1: closest to phi[0]
        cands = _candidates(1)
        phi[1] = cands[np.argmin(np.abs(cands - phi[0]))]

        # Steps 2..N-1: linear extrapolation
        for i in range(2, N):
            cands = _candidates(i)
            expected = 2 * phi[i - 1] - phi[i - 2]
            phi[i] = cands[np.argmin(np.abs(cands - expected))]

        B = np.gradient(phi, tau) / _get_sensitivity(self.qubit)
        return B


# ---------------------------------------------------------------------------
# Differential echo reconstruction
# ---------------------------------------------------------------------------

@dataclass
class DiffEchoReconstruction(Reconstruction):
    """Differential echo waveform reconstruction.

    Replaces Analysis.get_signal_from_diff_echo.

    varphi(τ) = arcsin(2*p_e - 1)
    B(τ) = -varphi / (2*k*kappa*t_int)

    Parameters
    ----------
    qubit
        QubitSpec or legacy TransmonQubit providing sensitivity().
    t_int : float
        Integration time (ns).
    k : int
        Number of π-pulse pairs.
    """

    qubit: object
    t_int: float
    k: int

    def reconstruct(
        self,
        measurement,
        kernel=None,
        calibration=None,
        **kwargs,
    ) -> np.ndarray:
        """Reconstruct B from differential echo p_e data.

        Parameters
        ----------
        measurement : ExperimentResult
            Must have data["p_e"].
        kernel, calibration : ignored

        Returns
        -------
        np.ndarray
            Reconstructed magnetic field B.
        """
        p_e_list = np.asarray(measurement.data["p_e"], dtype=float)
        varphi = np.arcsin(2 * p_e_list - 1)
        kappa = _get_sensitivity(self.qubit)
        B = -varphi / (2 * self.k * kappa * self.t_int)
        return B
