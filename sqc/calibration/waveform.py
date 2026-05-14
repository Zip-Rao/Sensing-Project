"""sqc.calibration.waveform — Waveform / control-line calibration.

Two components:
  - WaveformCalibration: unified entry for transfer-function measurement
    and predistortion filter design.
  - PredistortionDesigner: standalone inverse-filter designer (can be
    used independently or via WaveformCalibration).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import bilinear

from sqc.calibration.base import Calibration, CalibrationTable
from sqc.control.waveform import Waveform


# ===================================================================
# WaveformCalibration — unified entry
# ===================================================================

@dataclass
class WaveformCalibration(Calibration):
    """Waveform / control-line calibration.

    Methods
    -------
    - ``"transfer_function"``: measure step response, fit DistortionModel.
    - ``"predistortion"``: design predistortion filter. Requires a prior
      transfer_function result (passed via ``transfer_model``).

    Transfer-function parameters
    ----------------------------
    distortion : DistortionModel
        The distortion model to characterise (in simulation mode, its
        step response is measured directly).
    method : str
        "transfer_function" or "predistortion".
    fit_type : str
        Model to fit: "multi_exp", "single_exp", "fir", "iir".
    n_exp_components : int
        Number of exponential components for multi_exp fit.
    t_max : float
        Maximum time for step response (ns).
    n_points : int
        Number of time points.

    Predistortion parameters
    ------------------------
    transfer_model : DistortionModel or None
        Measured transfer function (from a prior transfer_function run).
    predistortion_method : str
        Inverse design method: "auto", "fir_inverse", "iir_inverse",
        "frequency_inverse".
    n_taps : int
        FIR tap count for fir_inverse fallback.
    regularization : float
        Regularisation for frequency-domain inverse.
    dt : float or None
        Sample spacing (ns). If None, inferred from target waveform.
    """

    # Shared
    method: Literal["transfer_function", "predistortion"] = "transfer_function"

    # -- transfer_function params --
    distortion: object | None = None
    fit_type: Literal["multi_exp", "single_exp", "fir", "iir"] = "multi_exp"
    n_exp_components: int = 3
    t_max: float = 500.0
    n_points: int = 2000

    # -- predistortion params --
    transfer_model: object | None = None
    predistortion_method: Literal[
        "auto", "fir_inverse", "iir_inverse", "frequency_inverse"
    ] = "auto"
    n_taps: int = 64
    regularization: float = 1e-4
    dt: float | None = None

    # ------------------------------------------------------------------
    def calibrate(self) -> CalibrationTable:
        match self.method:
            case "transfer_function" | "simulation":
                return self._calibrate_transfer_function()
            case "predistortion":
                return self._calibrate_predistortion()
            case _:
                raise ValueError(
                    f"Unknown method '{self.method}'. "
                    f"Valid: 'transfer_function', 'predistortion'."
                )

    # ------------------------------------------------------------------
    # Transfer function measurement
    # ------------------------------------------------------------------

    def _calibrate_transfer_function(self) -> CalibrationTable:
        t, step = self._measure_step_response()
        fit_params = self._fit_step_response(t, step)

        return CalibrationTable(
            name="transfer_function",
            qubit_name="",
            kind="transfer_function",
            inputs=t,
            outputs=step,
            fit_params=fit_params,
            metadata={
                "method": "simulation",
                "fit_type": self.fit_type,
                "class": type(self.distortion).__name__ if self.distortion else "",
            },
        )

    def _measure_step_response(self) -> tuple[np.ndarray, np.ndarray]:
        t = np.linspace(0, self.t_max, self.n_points)
        step = self.distortion.step_response(t)
        return t, step

    def _fit_step_response(self, t: np.ndarray, step: np.ndarray) -> dict:
        if self.fit_type == "single_exp":
            return self._fit_single_exp(t, step)
        elif self.fit_type == "multi_exp":
            return self._fit_multi_exp(t, step)
        elif self.fit_type == "fir":
            return self._fit_fir(t, step)
        elif self.fit_type == "iir":
            return self._fit_iir(t, step)
        raise ValueError(f"Unknown fit_type: {self.fit_type}")

    def _fit_single_exp(self, t: np.ndarray, step: np.ndarray) -> dict:
        mask = t >= 0

        def model(t_data, amp, tau):
            return 1.0 - amp * np.exp(-np.maximum(t_data, 0) / max(tau, 1e-9))

        p0 = [0.05, 50.0]
        try:
            popt, _ = curve_fit(
                model, t[mask], step[mask],
                p0=p0, bounds=(0, [1.0, 10000.0]),
                maxfev=5000,
            )
            return {"amplitude": float(popt[0]), "tau": float(popt[1])}
        except Exception:
            mask2 = (t > 0) & (step < 1.0)
            if np.any(mask2):
                diff = np.log(np.maximum(1.0 - step[mask2], 1e-15))
                A = np.column_stack([np.ones_like(t[mask2]), -t[mask2]])
                coeffs, _, _, _ = np.linalg.lstsq(A, diff, rcond=None)
                return {
                    "amplitude": float(np.exp(coeffs[0])),
                    "tau": float(1.0 / max(coeffs[1], 1e-9)),
                }
            return {"amplitude": 0.0, "tau": 1.0}

    def _fit_multi_exp(self, t: np.ndarray, step: np.ndarray) -> dict:
        K = self.n_exp_components
        amps, taus = [], []
        residual = step.copy()

        for _ in range(K):
            fit = self._fit_single_exp(t, residual)
            amp = fit["amplitude"]
            tau = fit["tau"]
            if amp < 1e-6:
                break
            amps.append(amp)
            taus.append(tau)
            residual = residual - amp * (1.0 - np.exp(-np.maximum(t, 0) / tau))

        if len(amps) == 0:
            return {"amplitudes": np.array([0.0]), "taus": np.array([1.0])}

        mask = t >= 0
        try:
            def model(t_data, *params):
                K2 = len(params) // 2
                result = np.ones_like(t_data)
                for k in range(K2):
                    result -= params[2 * k] * np.exp(
                        -np.maximum(t_data, 0) / max(params[2 * k + 1], 1e-9)
                    )
                return result

            p0 = []
            bounds_low, bounds_high = [], []
            for a, tau_val in zip(amps, taus):
                p0.extend([a, tau_val])
                bounds_low.extend([0, 0.1])
                bounds_high.extend([1.0, 10000.0])

            popt, _ = curve_fit(
                model, t[mask], step[mask],
                p0=p0, bounds=(bounds_low, bounds_high),
                maxfev=10000,
            )
            K_fit = len(popt) // 2
            amps_fit = np.array([popt[2 * k] for k in range(K_fit)])
            taus_fit = np.array([popt[2 * k + 1] for k in range(K_fit)])
        except Exception:
            amps_fit = np.array(amps)
            taus_fit = np.array(taus)

        return {
            "amplitudes": amps_fit, "taus": taus_fit,
            "n_components": len(amps_fit),
        }

    def _fit_fir(self, t: np.ndarray, step: np.ndarray) -> dict:
        if len(t) < 2:
            return {"taps": np.array([1.0])}
        dt_val = t[1] - t[0]
        imp = np.gradient(step, dt_val)
        n_taps = min(64, len(t))
        taps = imp[:n_taps]
        taps = taps / np.sum(np.abs(taps)) if np.sum(np.abs(taps)) > 0 else taps
        return {"taps": taps, "dt": float(dt_val)}

    def _fit_iir(self, t: np.ndarray, step: np.ndarray) -> dict:
        exp_fit = self._fit_multi_exp(t, step)
        if "amplitudes" in exp_fit and len(exp_fit["amplitudes"]) >= 1:
            amp = exp_fit["amplitudes"][0]
            tau_val = exp_fit["taus"][0]
            if len(t) >= 2:
                dt_val = t[1] - t[0]
                b_cont = [tau_val * (1.0 - amp), 1.0]
                a_cont = [tau_val, 1.0]
                b, a = bilinear(b_cont, a_cont, fs=1.0 / dt_val)
                return {"b": b, "a": a, "dt": float(dt_val)}
        return {"b": np.array([1.0]), "a": np.array([1.0]), "dt": 1.0}

    def to_distortion_model(self) -> object:
        """Convert calibration result into a DistortionModel."""
        from sqc.hardware.distortion import (
            SingleExponentialDistortion,
            MultiExponentialDistortion,
            FIRDistortion,
            IIRDistortion,
        )

        table = self.calibrate()

        if self.fit_type == "single_exp":
            return SingleExponentialDistortion(
                amplitude=table.fit_params.get("amplitude", 0.01),
                tau=table.fit_params.get("tau", 100.0),
            )
        elif self.fit_type == "multi_exp":
            return MultiExponentialDistortion(
                amplitudes=table.fit_params.get("amplitudes", np.array([0.01])),
                taus=table.fit_params.get("taus", np.array([100.0])),
            )
        elif self.fit_type == "fir":
            return FIRDistortion(
                taps=table.fit_params.get("taps", np.array([1.0])),
            )
        elif self.fit_type == "iir":
            return IIRDistortion(
                b_coeffs=table.fit_params.get("b", np.array([1.0])),
                a_coeffs=table.fit_params.get("a", np.array([1.0])),
            )
        raise ValueError(f"Unknown fit_type: {self.fit_type}")

    # ------------------------------------------------------------------
    # Predistortion (wraps PredistortionDesigner)
    # ------------------------------------------------------------------

    def _calibrate_predistortion(self) -> CalibrationTable:
        if self.transfer_model is None:
            raise ValueError(
                "transfer_model is required for method='predistortion'. "
                "Run transfer_function calibration first."
            )
        dt_val = self.dt or CONFIG.awg.dt
        designer = PredistortionDesigner(
            method=self.predistortion_method,
            n_taps=self.n_taps,
            regularization=self.regularization,
        )
        inverse = designer.design(self.transfer_model, dt=dt_val)
        return CalibrationTable(
            name="predistortion",
            qubit_name="",
            kind="predistortion",
            inputs=np.array([]),
            outputs=np.array([]),
            fit_params={"inverse_model": inverse, "dt": float(dt_val)},
            metadata={"method": self.predistortion_method},
        )


# ===================================================================
# PredistortionDesigner — standalone inverse-filter designer
# ===================================================================

@dataclass
class PredistortionDesigner:
    """Design a predistortion filter given a measured transfer function.

    For exponential-distortion models (single/multi), uses analytical IIR
    inverse.  For general models, uses frequency-domain inversion.

    Parameters
    ----------
    method : str
        "auto", "fir_inverse", "iir_inverse", or "frequency_inverse".
    n_taps : int
        FIR tap count for fir_inverse fallback.
    regularization : float
        Regularisation for frequency-domain inverse.
    """

    method: Literal[
        "auto", "fir_inverse", "iir_inverse", "frequency_inverse"
    ] = "auto"
    n_taps: int = 64
    regularization: float = 1e-4

    def design(self, transfer_model: object, dt: float) -> object:
        """Design an inverse DistortionModel.

        Parameters
        ----------
        transfer_model : DistortionModel
            The measured transfer function.
        dt : float
            Sample spacing (ns).

        Returns
        -------
        DistortionModel
            Inverse model.
        """
        resolved = self.method
        if resolved == "auto":
            resolved = self._auto_method(transfer_model)

        if resolved == "iir_inverse":
            return self._design_iir_inverse(transfer_model, dt)
        elif resolved == "fir_inverse":
            return self._design_fir_inverse(transfer_model, dt)
        elif resolved == "frequency_inverse":
            return self._design_frequency_inverse(transfer_model, dt)
        raise ValueError(f"Unknown method: {resolved}")

    @staticmethod
    def _auto_method(transfer_model: object) -> str:
        cls_name = type(transfer_model).__name__
        return "iir_inverse" if "Exponential" in cls_name else "frequency_inverse"

    # -- IIR inverse -------------------------------------------------------

    def _design_iir_inverse(self, transfer_model: object, dt: float) -> object:
        from sqc.hardware.distortion import (
            IIRDistortion,
            SingleExponentialDistortion,
            MultiExponentialDistortion,
        )
        if isinstance(transfer_model, SingleExponentialDistortion):
            return self._single_exp_to_iir_inverse(transfer_model, dt)
        if isinstance(transfer_model, MultiExponentialDistortion):
            return self._multi_exp_to_iir_inverse(transfer_model, dt)
        return self._design_frequency_inverse(transfer_model, dt)

    @staticmethod
    def _single_exp_to_iir_inverse(
        model: "SingleExponentialDistortion", dt: float,
    ) -> "IIRDistortion":
        from sqc.hardware.distortion import IIRDistortion

        amp = model.amplitude
        tau = model.tau
        b_cont = [tau, 1.0]
        a_cont = [tau * (1.0 - amp), 1.0]
        b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
        return IIRDistortion(b_coeffs=b, a_coeffs=a)

    @staticmethod
    def _multi_exp_to_iir_inverse(
        model: "MultiExponentialDistortion", dt: float,
    ) -> object:
        from sqc.hardware.distortion import CustomTransferDistortion

        n_fft = 4096
        omega_grid = 2.0 * np.pi * np.fft.fftfreq(n_fft, d=dt)
        H = model.frequency_response(omega_grid)
        H_inv = np.conj(H) / (np.abs(H) ** 2 + 1e-4)
        return CustomTransferDistortion(
            omega_grid=omega_grid.copy(), H_grid=H_inv.copy(),
        )

    # -- FIR inverse -------------------------------------------------------

    def _design_fir_inverse(self, transfer_model: object, dt: float) -> object:
        from sqc.hardware.distortion import FIRDistortion

        n_fft = 4096
        imp = np.zeros(n_fft)
        imp[0] = 1.0 / dt
        h_fwd = np.asarray(transfer_model.apply(imp, dt), dtype=float)
        H_emp = np.fft.fft(h_fwd)
        H_inv = np.conj(H_emp) / (np.abs(H_emp) ** 2 + self.regularization ** 2)
        h_inv_full = np.fft.ifft(H_inv).real

        n_keep = min(self.n_taps, len(h_inv_full))
        h_trunc = h_inv_full[:n_keep]
        if len(h_trunc) < self.n_taps:
            padded = np.zeros(self.n_taps)
            padded[:len(h_trunc)] = h_trunc
            h_trunc = padded
        return FIRDistortion(taps=h_trunc)

    # -- Frequency-domain inverse ------------------------------------------

    def _design_frequency_inverse(
        self, transfer_model: object, dt: float,
    ) -> object:
        from sqc.hardware.distortion import CustomTransferDistortion

        n_fft = 4096
        omega_grid = 2.0 * np.pi * np.fft.fftfreq(n_fft, d=dt)
        H = np.asarray(
            transfer_model.frequency_response(omega_grid), dtype=complex,
        )
        H_inv = np.conj(H) / (np.abs(H) ** 2 + self.regularization ** 2)
        return CustomTransferDistortion(
            omega_grid=omega_grid.copy(), H_grid=H_inv.copy(),
        )

    # -- Convenience: apply predistortion to a Waveform --------------------

    def predistort(
        self,
        target: Waveform,
        transfer: object | None = None,
        inverse_model: object | None = None,
    ) -> Waveform:
        """Apply predistortion to a target waveform.

        Provide either a transfer function (inverted internally) or a
        pre-computed inverse model.
        """
        if inverse_model is None:
            if transfer is None:
                raise ValueError(
                    "Either transfer or inverse_model must be provided."
                )
            dt_val = float(target.t_list[1] - target.t_list[0])
            inverse_model = self.design(transfer, dt=dt_val)
        return inverse_model.apply_to_waveform(target)

    @staticmethod
    def check_pole_stability(
        b_coeffs: np.ndarray, a_coeffs: np.ndarray,
    ) -> bool:
        """Check if IIR filter is stable (all poles inside unit circle)."""
        roots = np.roots(a_coeffs)
        return bool(np.all(np.abs(roots) < 1.0 - 1e-10))


# ---------------------------------------------------------------------------
# Module-level convenience import
# ---------------------------------------------------------------------------
from sqc.config import CONFIG
