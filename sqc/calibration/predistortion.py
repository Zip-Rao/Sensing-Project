"""sqc.calibration.predistortion — PredistortionDesigner.

Designs FIR/IIR/frequency-domain predistortion filters given a measured
transfer function. Creates an inverse DistortionModel that, when applied
to the target waveform, produces the AWG waveform needed to achieve that
target on-chip.

Per phase_4_handbook.md §3.4.

For Single/MultiExponentialDistortion, uses an analytical IIR inverse
derived from the continuous-time transfer function, discretised via
bilinear transform. This guarantees that the forward + inverse cascade
is numerically exact (to within floating-point precision) when both
use the same dt.

For general DistortionModel subclasses, uses a frequency-domain inverse
(CustomTransferDistortion) via FFT convolution.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional

import numpy as np
from scipy.signal import bilinear, lfilter


@dataclass
class PredistortionDesigner:
    """Designs predistortion filter given measured transfer function.

    For exponential-distortion models (single/multi), uses analytical IIR
    inverse.  For general models, uses frequency-domain inversion.

    Parameters
    ----------
    method : Literal["auto", "fir_inverse", "iir_inverse", "frequency_inverse"]
        Inverse design method. "auto" detects the distortion type and
        chooses the best method. Default "auto".
    n_taps : int
        Number of FIR taps for fir_inverse fallback. Default 64.
    regularization : float
        Regularization parameter for frequency-domain inverse. Default 1e-4.
    """

    method: Literal["auto", "fir_inverse", "iir_inverse", "frequency_inverse"] = "auto"
    n_taps: int = 64
    regularization: float = 1e-4

    def design(
        self,
        transfer_model: object,
        dt: float,
    ) -> object:
        """Design an inverse DistortionModel for the given transfer function.

        Parameters
        ----------
        transfer_model : DistortionModel
            The measured transfer function.
        dt : float
            Sample spacing (ns) — determines the discrete-time filter.

        Returns
        -------
        DistortionModel
            Inverse model.  Type depends on method and transfer_model class.
        """
        resolved_method = self.method
        if resolved_method == "auto":
            resolved_method = self._auto_method(transfer_model)

        if resolved_method == "iir_inverse":
            return self._design_iir_inverse(transfer_model, dt)
        elif resolved_method == "fir_inverse":
            return self._design_fir_inverse(transfer_model, dt)
        elif resolved_method == "frequency_inverse":
            return self._design_frequency_inverse(transfer_model, dt)
        else:
            raise ValueError(f"Unknown method: {resolved_method}")

    @staticmethod
    def _auto_method(transfer_model: object) -> str:
        """Heuristic: use IIR for exponential models, frequency for others."""
        cls_name = type(transfer_model).__name__
        if "Exponential" in cls_name:
            return "iir_inverse"
        return "frequency_inverse"

    # ------------------------------------------------------------------
    # IIR inverse (analytical, for exponential distortion)
    # ------------------------------------------------------------------

    def _design_iir_inverse(
        self, transfer_model: object, dt: float,
    ) -> object:
        """Design IIR inverse analytically from the transfer function.

        For SingleExponentialDistortion:
            H(s) = (1 - amp) + amp / (1 + s*tau)
                 = (1 + s*tau*(1-amp)) / (1 + s*tau)
            H_inv(s) = (1 + s*tau) / (1 + s*tau*(1-amp))
            Disc retised via bilinear transform -> (b_inv, a_inv).

        For MultiExponentialDistortion:
            TODO: analytical multi-exp inverse (P4 uses numerical fallback).
        """
        from sqc.hardware.distortion import (
            IIRDistortion,
            SingleExponentialDistortion,
            MultiExponentialDistortion,
        )

        if isinstance(transfer_model, SingleExponentialDistortion):
            return self._single_exp_to_iir_inverse(transfer_model, dt)

        if isinstance(transfer_model, MultiExponentialDistortion):
            return self._multi_exp_to_iir_inverse(transfer_model, dt)

        # Generic fallback: numerical frequency-domain inverse
        return self._design_frequency_inverse(transfer_model, dt)

    @staticmethod
    def _single_exp_to_iir_inverse(
        model: "SingleExponentialDistortion", dt: float,
    ) -> "IIRDistortion":
        """Analytical IIR inverse for SingleExponentialDistortion.

        H(s) = (1 + s*tau*(1-amp)) / (1 + s*tau)
        H_inv(s) = (1 + s*tau) / (1 + s*tau*(1-amp))

        Bilinear transform to discrete time.
        """
        from sqc.hardware.distortion import IIRDistortion

        amp = model.amplitude
        tau = model.tau

        # Numerator: 1 + s*tau  →  [tau, 1] in s-domain
        # Denominator: 1 + s*tau*(1-amp)  →  [tau*(1-amp), 1]
        b_cont = [tau, 1.0]
        a_cont = [tau * (1.0 - amp), 1.0]

        b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
        return IIRDistortion(b_coeffs=b, a_coeffs=a)

    @staticmethod
    def _multi_exp_to_iir_inverse(
        model: "MultiExponentialDistortion", dt: float,
    ) -> object:
        """Approximate IIR inverse for MultiExponentialDistortion.

        Strategy: cascade (in parallel) the inverse of each single-exp
        component. Since the components are parallel in the forward model,
        their inverses are also parallel. But the exact inverse of a sum
        of transfer functions does not decompose into parallel inverses.

        For P4: use frequency-domain fallback for multi-exp.
        """
        from sqc.hardware.distortion import CustomTransferDistortion

        # Build frequency-domain inverse
        n_fft = 4096
        omega_grid = 2.0 * np.pi * np.fft.fftfreq(n_fft, d=dt)
        H = model.frequency_response(omega_grid)
        H_inv = np.conj(H) / (np.abs(H) ** 2 + 1e-4)
        return CustomTransferDistortion(
            omega_grid=omega_grid.copy(),
            H_grid=H_inv.copy(),
        )

    # ------------------------------------------------------------------
    # FIR inverse (numerical, empirical)
    # ------------------------------------------------------------------

    def _design_fir_inverse(
        self, transfer_model: object, dt: float,
    ) -> object:
        """Design FIR inverse via empirical impulse response measurement.

        1. Feed impulse through forward model to get h_fwd[n].
        2. FFT to get H_emp(omega).
        3. Regularized inverse in frequency domain.
        4. IFFT to get h_inv[n], truncated to n_taps.
        """
        from sqc.hardware.distortion import FIRDistortion

        n_fft = 4096
        # Impulse input
        imp = np.zeros(n_fft)
        imp[0] = 1.0 / dt
        h_fwd = np.asarray(transfer_model.apply(imp, dt), dtype=float)

        # Frequency-domain inverse
        H_emp = np.fft.fft(h_fwd)
        H_inv = np.conj(H_emp) / (np.abs(H_emp) ** 2 + self.regularization ** 2)

        # IFFT to time-domain kernel
        h_inv_full = np.fft.ifft(H_inv).real

        # Truncate to n_taps
        n_keep = min(self.n_taps, len(h_inv_full))
        h_trunc = h_inv_full[:n_keep]
        if len(h_trunc) < self.n_taps:
            padded = np.zeros(self.n_taps)
            padded[:len(h_trunc)] = h_trunc
            h_trunc = padded

        return FIRDistortion(taps=h_trunc)

    # ------------------------------------------------------------------
    # Frequency-domain inverse (general, FFT convolution)
    # ------------------------------------------------------------------

    def _design_frequency_inverse(
        self, transfer_model: object, dt: float,
    ) -> object:
        """Design frequency-domain inverse as CustomTransferDistortion.

        Uses the continuous-time frequency response H(omega) evaluated
        on the discrete FFT grid. The inverse is applied via
        FFT -> multiply -> IFFT.
        """
        from sqc.hardware.distortion import CustomTransferDistortion

        n_fft = 4096
        omega_grid = 2.0 * np.pi * np.fft.fftfreq(n_fft, d=dt)
        H = np.asarray(
            transfer_model.frequency_response(omega_grid), dtype=complex,
        )
        H_inv = np.conj(H) / (np.abs(H) ** 2 + self.regularization ** 2)
        return CustomTransferDistortion(
            omega_grid=omega_grid.copy(),
            H_grid=H_inv.copy(),
        )

    # ------------------------------------------------------------------
    # Convenience: apply predistortion to a Waveform
    # ------------------------------------------------------------------

    def predistort(
        self,
        target: "Waveform",
        transfer: Optional[object] = None,
        inverse_model: Optional[object] = None,
    ) -> "Waveform":
        """Apply predistortion to a target waveform.

        Either provide a transfer function (which is inverted internally)
        or a pre-computed inverse model.

        Parameters
        ----------
        target : Waveform
            Desired on-chip waveform.
        transfer : DistortionModel or None
            The transfer function to invert.
        inverse_model : DistortionModel or None
            Pre-computed inverse model.

        Returns
        -------
        Waveform
            Predistorted waveform to send to AWG.
        """
        if inverse_model is None:
            if transfer is None:
                raise ValueError(
                    "Either transfer or inverse_model must be provided."
                )
            dt = float(target.t_list[1] - target.t_list[0])
            inverse_model = self.design(transfer, dt=dt)

        return inverse_model.apply_to_waveform(target)

    @staticmethod
    def check_pole_stability(b_coeffs: np.ndarray, a_coeffs: np.ndarray) -> bool:
        """Check if IIR filter is stable (all poles inside unit circle)."""
        roots = np.roots(a_coeffs)
        return np.all(np.abs(roots) < 1.0 - 1e-10)
