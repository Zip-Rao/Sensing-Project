"""sqc.calibration.transfer_function — TransferFunctionCalibration.

Measures step response of a control line via simulation (or Cryoscope
when Track B 1.1 is complete) and fits a DistortionModel.

Per phase_4_handbook.md §3.3.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional

import numpy as np
from scipy.optimize import curve_fit

from sqc.calibration.base import Calibration, CalibrationTable
from sqc.control.waveform import Waveform


@dataclass
class TransferFunctionCalibration(Calibration):
    """Measure step response of a control line, fit transfer function H(omega).

    In simulation mode: directly measures the step response from a known
    DistortionModel. When Track B 1.1 Cryoscope is complete, this will
    measure via CryoscopeExperiment for hardware-based calibration.

    Strategy:
    1. Apply known step via ControlLine distortion model.
    2. Fit IIR / multi-exp / FIR model to the step response.
    3. Return CalibrationTable with fitted parameters.

    Parameters
    ----------
    distortion : DistortionModel
        The distortion model to calibrate (in simulation, we measure
        its step response directly).
    method : Literal["simulation", "cryoscope", "transient", "hybrid"]
        Measurement method. "simulation" uses the distortion model directly.
        Hardware methods (cryoscope/transient/hybrid) require Track B 1.1/1.2.
    fit_type : Literal["multi_exp", "single_exp", "fir", "iir"]
        Type of model to fit.
    n_exp_components : int
        Number of exponential components for multi_exp fit. Default 3.
    t_max : float
        Maximum time for step response (ns). Default 500.0.
    n_points : int
        Number of time points. Default 2000.
    """

    distortion: object  # DistortionModel (duck-typed)
    method: Literal["simulation", "cryoscope", "transient", "hybrid"] = "simulation"
    fit_type: Literal["multi_exp", "single_exp", "fir", "iir"] = "multi_exp"
    n_exp_components: int = 3
    t_max: float = 500.0
    n_points: int = 2000

    def calibrate(self) -> CalibrationTable:
        """Run the calibration workflow and return results.

        Returns
        -------
        CalibrationTable
            Calibration results with fitted transfer function parameters.
        """
        # 1. Measure step response
        t, step = self._measure_step_response()

        # 2. Fit chosen model
        fit_params = self._fit_step_response(t, step)

        return CalibrationTable(
            name="transfer_function",
            qubit_name="",  # not qubit-specific
            kind="transfer_function",
            inputs=t,
            outputs=step,
            fit_params=fit_params,
            metadata={
                "method": self.method,
                "fit_type": self.fit_type,
                "class": type(self.distortion).__name__,
            },
        )

    def _measure_step_response(self) -> tuple[np.ndarray, np.ndarray]:
        """Measure the step response.

        In simulation mode, uses the distortion model directly.
        """
        if self.method == "simulation":
            return self._measure_simulation()
        elif self.method == "cryoscope":
            return self._measure_cryoscope()
        elif self.method == "transient":
            return self._measure_transient()
        elif self.method == "hybrid":
            return self._measure_hybrid()
        else:
            raise ValueError(f"Unknown method: {self.method}")

    def _measure_simulation(self) -> tuple[np.ndarray, np.ndarray]:
        """Direct measurement from the distortion model."""
        t = np.linspace(0, self.t_max, self.n_points)
        step = self.distortion.step_response(t)
        return t, step

    def _measure_cryoscope(self) -> tuple[np.ndarray, np.ndarray]:
        """Measure step response via CryoscopeExperiment.

        Requires Track B 1.1 (Cryoscope case 6/7).
        """
        raise NotImplementedError(
            "Cryoscope-based transfer function measurement requires "
            "Track B 1.1 (Cryoscope case 6/7). Use method='simulation'."
        )

    def _measure_transient(self) -> tuple[np.ndarray, np.ndarray]:
        """Measure step response via sliding measurement.

        Requires Track B 1.2 (transient calibration case 8).
        """
        raise NotImplementedError(
            "Transient-based transfer function measurement requires "
            "Track B 1.2 (transient calibration case 8). Use method='simulation'."
        )

    def _measure_hybrid(self) -> tuple[np.ndarray, np.ndarray]:
        """Hybrid: Cryoscope for slow tails, transient for fast edge."""
        raise NotImplementedError(
            "Hybrid measurement requires Track B 1.1 and 1.2. "
            "Use method='simulation'."
        )

    def _fit_step_response(
        self, t: np.ndarray, step: np.ndarray,
    ) -> dict:
        """Fit the step response to the chosen model type."""
        if self.fit_type == "single_exp":
            return self._fit_single_exp(t, step)
        elif self.fit_type == "multi_exp":
            return self._fit_multi_exp(t, step)
        elif self.fit_type == "fir":
            return self._fit_fir(t, step)
        elif self.fit_type == "iir":
            return self._fit_iir(t, step)
        else:
            raise ValueError(f"Unknown fit_type: {self.fit_type}")

    def _fit_single_exp(self, t: np.ndarray, step: np.ndarray) -> dict:
        """Fit s(t) = 1 - amp * exp(-t/tau)."""
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
            # Fallback: simple exponential fit via log
            mask2 = (t > 0) & (step < 1.0)
            if np.any(mask2):
                diff = np.log(np.maximum(1.0 - step[mask2], 1e-15))
                # Fit line: log(1-s(t)) = log(amp) - t/tau
                A = np.column_stack([np.ones_like(t[mask2]), -t[mask2]])
                coeffs, _, _, _ = np.linalg.lstsq(A, diff, rcond=None)
                return {
                    "amplitude": float(np.exp(coeffs[0])),
                    "tau": float(1.0 / max(coeffs[1], 1e-9)),
                }
            return {"amplitude": 0.0, "tau": 1.0}

    def _fit_multi_exp(self, t: np.ndarray, step: np.ndarray) -> dict:
        """Fit s(t) = 1 - sum_k amp_k * exp(-t/tau_k).

        Uses iterative fitting: fit single exponential, subtract, repeat.
        """
        K = self.n_exp_components
        amps = []
        taus = []
        residual = step.copy()

        for _ in range(K):
            fit = self._fit_single_exp(t, residual)
            amp = fit["amplitude"]
            tau = fit["tau"]
            if amp < 1e-6:
                break
            amps.append(amp)
            taus.append(tau)
            # Subtract fitted component and refit
            residual = residual - amp * (1.0 - np.exp(-np.maximum(t, 0) / tau))

        if len(amps) == 0:
            return {"amplitudes": np.array([0.0]), "taus": np.array([1.0])}

        # Refine with joint curve_fit
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
            bounds_low = []
            bounds_high = []
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
            "amplitudes": amps_fit,
            "taus": taus_fit,
            "n_components": len(amps_fit),
        }

    def _fit_fir(self, t: np.ndarray, step: np.ndarray) -> dict:
        """Fit FIR model to step response."""
        # Step response to impulse response via gradient
        if len(t) < 2:
            return {"taps": np.array([1.0])}
        dt = t[1] - t[0]
        imp = np.gradient(step, dt)
        # Take first N significant taps
        n_taps = min(64, len(t))
        taps = imp[:n_taps]
        # Normalize
        taps = taps / np.sum(np.abs(taps)) if np.sum(np.abs(taps)) > 0 else taps
        return {"taps": taps, "dt": float(dt)}

    def _fit_iir(self, t: np.ndarray, step: np.ndarray) -> dict:
        """Fit IIR model to step response.

        For now, fit a multi-exp model and convert to IIR.
        """
        exp_fit = self._fit_multi_exp(t, step)
        # Convert multi-exp to IIR: each exp component is a first-order filter
        # The parallel connection of K first-order filters -> K-th order
        # For simplicity, return coefficients for a first-order approximation
        if "amplitudes" in exp_fit and len(exp_fit["amplitudes"]) >= 1:
            amp = exp_fit["amplitudes"][0]
            tau = exp_fit["taus"][0]
            if len(t) >= 2:
                dt = t[1] - t[0]
                from scipy.signal import bilinear

                b_cont = [tau * (1.0 - amp), 1.0]
                a_cont = [tau, 1.0]
                b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
                return {"b": b, "a": a, "dt": float(dt)}
        return {"b": np.array([1.0]), "a": np.array([1.0]), "dt": 1.0}

    def to_distortion_model(self) -> object:
        """Convert the calibration result into a DistortionModel.

        Returns
        -------
        DistortionModel
            The fitted distortion model.

        Raises
        ------
        RuntimeError
            If calibrate() has not been called.
        """
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
        else:
            raise ValueError(f"Unknown fit_type: {self.fit_type}")
