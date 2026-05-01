"""Transfer-function calibration for control-line distortion models."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from sqc.calibration.base import Calibration, CalibrationTable
from sqc.control.waveform import Waveform
from sqc.hardware.distortion import (
    CustomTransferDistortion,
    DistortionModel,
    FIRDistortion,
    IIRDistortion,
    MultiExponentialDistortion,
)


def _uniform_dt(t_list: np.ndarray) -> float:
    t = np.asarray(t_list, dtype=float)
    if t.ndim != 1 or len(t) < 2:
        raise ValueError("transfer calibration requires at least two time samples")
    diffs = np.diff(t)
    dt = float(np.median(diffs))
    if dt <= 0 or not np.allclose(diffs, dt, rtol=1e-5, atol=1e-12):
        raise ValueError("transfer calibration requires a uniform time grid")
    return dt


def _qubit_name(qubit) -> str:
    if qubit is None:
        return "unknown"
    if hasattr(qubit, "spec"):
        return qubit.spec().name
    if hasattr(qubit, "name"):
        return str(qubit.name)
    return type(qubit).__name__


@dataclass
class TransferFunctionCalibration(Calibration):
    """Calibrate a line transfer function from a measured or simulated step.

    Real cryoscope/transient acquisition is still outside this phase's migrated
    code. This class therefore accepts either a measured step response or a
    known simulation model. If neither is provided it raises a clear blocker
    instead of inventing measurement data.
    """

    qubit: object | None = None
    method: Literal["cryoscope", "transient", "hybrid", "simulated"] = "simulated"
    fit_type: Literal["multi_exp", "iir", "fir", "custom"] = "custom"
    n_exp_components: int = 3
    measured_step_response: Waveform | None = None
    transfer_model: DistortionModel | None = None
    t_list: np.ndarray | None = None

    def calibrate(self) -> CalibrationTable:
        step = self._measure_step_response()
        fit_params = self._fit(step)
        return CalibrationTable(
            qubit_name=_qubit_name(self.qubit),
            kind="transfer_function",
            inputs=step.t_list,
            outputs=step.samples,
            fit_params=fit_params,
            metadata={
                "method": self.method,
                "fit_type": self.fit_type,
                "source": "measured"
                if self.measured_step_response is not None
                else "simulated_model",
            },
        )

    def _measure_step_response(self) -> Waveform:
        if self.measured_step_response is not None:
            return self.measured_step_response.copy()

        if self.transfer_model is not None:
            t = (
                np.asarray(self.t_list, dtype=float)
                if self.t_list is not None
                else np.linspace(0.0, 100.0, 1001)
            )
            return Waveform(
                t_list=t.copy(),
                samples=np.asarray(self.transfer_model.step_response(t), dtype=float),
                metadata={"source_model": type(self.transfer_model).__name__},
            )

        raise NotImplementedError(
            "TransferFunctionCalibration needs measured_step_response or "
            "transfer_model. Live cryoscope/transient acquisition is a later "
            "Track B dependency."
        )

    def _fit(self, step: Waveform) -> dict:
        if self.fit_type == "custom":
            return self._fit_custom(step)
        if self.fit_type == "fir":
            return self._fit_fir(step)
        if self.fit_type == "iir":
            return self._fit_iir(step)
        if self.fit_type == "multi_exp":
            return self._fit_multi_exp(step)
        raise ValueError(f"Unknown fit_type: {self.fit_type}")

    def _fit_custom(self, step: Waveform) -> dict:
        dt = _uniform_dt(step.t_list)
        omega = 2.0 * np.pi * np.fft.fftfreq(len(step.samples), d=dt)
        if self.transfer_model is not None and hasattr(
            self.transfer_model, "discrete_frequency_response"
        ):
            h = self.transfer_model.discrete_frequency_response(omega, dt)  # type: ignore[attr-defined]
        elif self.transfer_model is not None:
            h = self.transfer_model.frequency_response(omega)
        else:
            impulse = np.diff(step.samples, prepend=0.0)
            h = np.fft.fft(impulse)
        order = np.argsort(omega)
        return {
            "model_type": "custom",
            "dt": dt,
            "omega": omega[order],
            "H_real": h.real[order],
            "H_imag": h.imag[order],
        }

    def _fit_fir(self, step: Waveform) -> dict:
        dt = _uniform_dt(step.t_list)
        taps = np.diff(step.samples, prepend=0.0)
        return {"model_type": "fir", "dt": dt, "taps": taps}

    def _fit_iir(self, step: Waveform) -> dict:
        dt = _uniform_dt(step.t_list)
        if self.transfer_model is not None and hasattr(
            self.transfer_model, "iir_coefficients"
        ):
            b, a = self.transfer_model.iir_coefficients(dt)  # type: ignore[attr-defined]
        elif isinstance(self.transfer_model, IIRDistortion):
            b, a = self.transfer_model.b, self.transfer_model.a
        else:
            params = self._fit_multi_exp(step)
            model = MultiExponentialDistortion(
                amplitudes=params["amplitudes"], taus=params["taus"]
            )
            b, a = model.iir_coefficients(dt)
        return {"model_type": "iir", "dt": dt, "b": np.asarray(b), "a": np.asarray(a)}

    def _fit_multi_exp(self, step: Waveform) -> dict:
        from scipy.optimize import curve_fit

        t = np.asarray(step.t_list, dtype=float)
        y = np.asarray(step.samples, dtype=float)
        t0 = t[0]
        x = t - t0
        dt = _uniform_dt(t)
        n_components = max(1, int(self.n_exp_components))

        def model(x_data, *params):
            amps = np.asarray(params[:n_components])
            log_taus = np.asarray(params[n_components:])
            taus = np.exp(log_taus)
            out = np.ones_like(x_data, dtype=float)
            for amp, tau in zip(amps, taus):
                out -= amp * np.exp(-x_data / tau)
            return out

        residual0 = float(np.clip(1.0 - y[0], -2.0, 2.0))
        amp0 = np.full(n_components, residual0 / n_components)
        duration = max(float(x[-1] - x[0]), dt)
        tau0 = np.geomspace(max(dt, duration / 200.0), max(2 * dt, duration / 4.0), n_components)
        p0 = np.concatenate([amp0, np.log(tau0)])
        lower = np.concatenate([np.full(n_components, -2.0), np.log(np.full(n_components, dt / 10.0))])
        upper = np.concatenate([np.full(n_components, 2.0), np.log(np.full(n_components, duration * 20.0))])

        try:
            popt, _ = curve_fit(
                model,
                x,
                y,
                p0=p0,
                bounds=(lower, upper),
                maxfev=20000,
            )
        except Exception:
            popt = p0

        amps = np.asarray(popt[:n_components], dtype=float)
        taus = np.exp(np.asarray(popt[n_components:], dtype=float))
        return {
            "model_type": "multi_exp",
            "dt": dt,
            "t0": t0,
            "amplitudes": amps,
            "taus": taus,
        }

    @staticmethod
    def model_from_table(table: CalibrationTable) -> DistortionModel:
        """Rebuild a distortion model from a transfer-function table."""
        params = dict(table.fit_params)
        model_type = params.get("model_type")
        if model_type == "custom":
            return CustomTransferDistortion(
                omega=np.asarray(params["omega"], dtype=float),
                H=np.asarray(params["H_real"], dtype=float)
                + 1j * np.asarray(params["H_imag"], dtype=float),
            )
        if model_type == "fir":
            return FIRDistortion(
                taps=np.asarray(params["taps"], dtype=float),
                dt=float(params["dt"]),
            )
        if model_type == "iir":
            return IIRDistortion(
                b=np.asarray(params["b"], dtype=float),
                a=np.asarray(params["a"], dtype=float),
                dt=float(params["dt"]),
            )
        if model_type == "multi_exp":
            return MultiExponentialDistortion(
                amplitudes=np.asarray(params["amplitudes"], dtype=float),
                taus=np.asarray(params["taus"], dtype=float),
            )

        dt = _uniform_dt(table.inputs)
        taps = np.diff(np.asarray(table.outputs, dtype=float), prepend=0.0)
        return FIRDistortion(taps=taps, dt=dt)


__all__ = ["TransferFunctionCalibration"]
