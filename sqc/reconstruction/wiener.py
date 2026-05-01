"""Linear and Ramsey-family waveform reconstruction methods."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.reconstruction.base import Reconstruction


@dataclass
class WienerReconstruction(Reconstruction):
    """Linear Wiener deconvolution for sliding Ramsey measurements."""

    lambda_reg: float = 1.0

    def reconstruct(self, measurement, kernel, calibration=None, dt=None) -> FluxSignal:
        """Return a FluxSignal reconstructed from ``delta_p`` and ``kernel``."""
        delta_p = np.asarray(measurement.data["delta_p"], dtype=float)
        kernel = np.asarray(kernel, dtype=float)
        if dt is None:
            dt = measurement.axes["scan"][1] - measurement.axes["scan"][0]

        n_del = len(delta_p)
        n_ker = len(kernel)
        n_signal = n_del - n_ker + 1
        n_fft = n_del

        p_pad = np.zeros(n_fft)
        p_pad[:n_del] = delta_p
        k_pad = np.zeros(n_fft)
        k_pad[:n_ker] = kernel

        measured_fft = np.fft.fft(p_pad)
        kernel_fft = np.fft.fft(k_pad)
        gain = np.conj(kernel_fft) / (np.abs(kernel_fft) ** 2 + self.lambda_reg**2)
        recovered = np.real(np.fft.ifft(measured_fft * gain / dt))

        t_list = np.linspace(0, (n_signal - 1) * dt, n_signal)
        return FluxSignal(type=8, t_list=t_list, signal=recovered[:n_signal])


@dataclass
class RamseyIQReconstruction(Reconstruction):
    """Magnetic-flux reconstruction from Ramsey IQ probabilities."""

    qubit: object

    def reconstruct(self, measurement, kernel=None, calibration=None) -> np.ndarray:
        """Return B(tau) from I/Q Ramsey traces."""
        tau_list = np.asarray(measurement.axes["tau"], dtype=float)
        p_e_i = np.asarray(measurement.data["p_e_I"], dtype=float)
        p_e_q = np.asarray(measurement.data["p_e_Q"], dtype=float)
        c_i = np.max(p_e_i) - np.min(p_e_i)
        c_q = np.max(p_e_q) - np.min(p_e_q)
        cosphi = (1 - 2 * p_e_i) / c_i
        sinphi = (1 - 2 * p_e_q) / c_q
        phi = np.unwrap(np.arctan2(sinphi, cosphi))
        return np.gradient(phi, tau_list) / self._sensitivity()

    def _sensitivity(self) -> float:
        if hasattr(self.qubit, "sensitivity"):
            return self.qubit.sensitivity()
        return self.qubit.frequency_sensitivity(self.qubit.flux)


@dataclass
class RamseyUnwrapReconstruction(Reconstruction):
    """Magnetic-flux reconstruction from single-channel Ramsey data."""

    qubit: object
    k_span: int = 3

    def reconstruct(self, measurement, kernel=None, calibration=None) -> np.ndarray:
        """Return B(tau) by branch-tracking the Ramsey phase."""
        tau = np.asarray(measurement.axes["tau"], dtype=float)
        p_e = np.asarray(measurement.data["p_e"], dtype=float)
        n_points = len(p_e)

        cos_phi = 1 - 2 * p_e
        theta = np.arccos(np.clip(cos_phi, -1, 1))
        ks = np.arange(-self.k_span, self.k_span + 1)

        def candidates(i):
            return np.array([2 * k * np.pi + sign * theta[i] for k in ks for sign in [1, -1]])

        phi = np.zeros(n_points)
        cands = candidates(0)
        phi[0] = cands[np.argmin(np.abs(cands))]
        if n_points > 1:
            cands = candidates(1)
            phi[1] = cands[np.argmin(np.abs(cands - phi[0]))]
        for i in range(2, n_points):
            cands = candidates(i)
            expected = 2 * phi[i - 1] - phi[i - 2]
            phi[i] = cands[np.argmin(np.abs(cands - expected))]

        return np.gradient(phi, tau) / self._sensitivity()

    def _sensitivity(self) -> float:
        if hasattr(self.qubit, "sensitivity"):
            return self.qubit.sensitivity()
        return self.qubit.frequency_sensitivity(self.qubit.flux)


@dataclass
class DiffEchoReconstruction(Reconstruction):
    """Direct differential-echo reconstruction."""

    qubit: object
    t_int: float
    k: int

    def reconstruct(self, measurement, kernel=None, calibration=None) -> np.ndarray:
        """Return B(tau) from differential echo p_e."""
        p_e = np.asarray(measurement.data["p_e"], dtype=float)
        varphi = np.arcsin(2 * p_e - 1)
        if hasattr(self.qubit, "sensitivity"):
            kappa = self.qubit.sensitivity()
        else:
            kappa = self.qubit.frequency_sensitivity(self.qubit.flux)
        return -varphi / (2 * self.k * kappa * self.t_int)
