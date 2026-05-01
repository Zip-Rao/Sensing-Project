"""Predistortion filter design for calibrated control lines."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from sqc.control.waveform import Waveform
from sqc.hardware.distortion import (
    CustomTransferDistortion,
    DistortionModel,
    FIRDistortion,
    IIRDistortion,
    discrete_iir_frequency_response,
)


def _uniform_dt(t_list: np.ndarray) -> float:
    t = np.asarray(t_list, dtype=float)
    if t.ndim != 1 or len(t) < 2:
        raise ValueError("predistortion requires at least two time samples")
    diffs = np.diff(t)
    dt = float(np.median(diffs))
    if dt <= 0 or not np.allclose(diffs, dt, rtol=1e-5, atol=1e-12):
        raise ValueError("predistortion requires a uniform time grid")
    return dt


def _next_pow2(n: int) -> int:
    return 1 << max(1, int(n - 1).bit_length())


@dataclass
class PredistortionDesigner:
    """Design an inverse filter for a measured control-line transfer function."""

    method: Literal["iir_inverse", "fir_inverse", "frequency_inverse"] = "iir_inverse"
    n_taps: int = 64
    regularization: float = 1e-4
    n_freq: int = 4096

    def __post_init__(self) -> None:
        if self.n_taps <= 0:
            raise ValueError("n_taps must be positive")
        if self.regularization < 0:
            raise ValueError("regularization must be non-negative")
        if self.n_freq < 2:
            raise ValueError("n_freq must be at least two")

    def design(
        self,
        transfer_model: DistortionModel,
        omega_grid: np.ndarray | None = None,
        dt: float | None = None,
    ) -> DistortionModel:
        """Return a model that approximates the inverse transfer function."""
        if self.method == "iir_inverse":
            inverse = self._design_iir_inverse(transfer_model, dt)
            if inverse is not None:
                return inverse
            return self._design_frequency_inverse(transfer_model, omega_grid, dt)
        if self.method == "fir_inverse":
            return self._design_fir_inverse(transfer_model, dt)
        if self.method == "frequency_inverse":
            return self._design_frequency_inverse(transfer_model, omega_grid, dt)
        raise ValueError(f"Unknown predistortion method: {self.method}")

    def predistort(
        self,
        target: Waveform,
        transfer: DistortionModel | None,
    ) -> Waveform:
        """Apply the designed inverse filter to a target on-chip waveform."""
        if transfer is None:
            return target.copy()
        dt = _uniform_dt(target.t_list)
        inverse_model = self.design(transfer, dt=dt)
        awg = inverse_model.apply(target)
        awg.metadata.update(
            {
                "predistortion": self.method,
                "target_distortion": type(transfer).__name__,
            }
        )
        return awg

    def _response(
        self,
        transfer_model: DistortionModel,
        omega: np.ndarray,
        dt: float | None,
    ) -> np.ndarray:
        if dt is not None and hasattr(transfer_model, "discrete_frequency_response"):
            return transfer_model.discrete_frequency_response(omega, dt)  # type: ignore[attr-defined]
        return transfer_model.frequency_response(omega)

    def _regularized_inverse(self, h: np.ndarray) -> np.ndarray:
        return np.conj(h) / (np.abs(h) ** 2 + self.regularization**2)

    def _omega_grid(self, dt: float | None) -> np.ndarray:
        effective_dt = 1.0 if dt is None else float(dt)
        n = _next_pow2(max(self.n_freq, 8 * self.n_taps))
        return 2.0 * np.pi * np.fft.fftfreq(n, d=effective_dt)

    def _design_frequency_inverse(
        self,
        transfer_model: DistortionModel,
        omega_grid: np.ndarray | None,
        dt: float | None,
    ) -> CustomTransferDistortion:
        if omega_grid is None:
            omega_grid = self._omega_grid(dt)
        omega_grid = np.asarray(omega_grid, dtype=float)
        h = self._response(transfer_model, omega_grid, dt)
        h_inv = self._regularized_inverse(h)
        order = np.argsort(omega_grid)
        return CustomTransferDistortion(omega=omega_grid[order], H=h_inv[order])

    def _design_fir_inverse(
        self,
        transfer_model: DistortionModel,
        dt: float | None,
    ) -> FIRDistortion:
        if dt is None:
            dt = getattr(transfer_model, "dt", None)
        if dt is None:
            dt = 1.0
        n_fft = _next_pow2(max(self.n_freq, 16 * self.n_taps))
        omega = 2.0 * np.pi * np.fft.fftfreq(n_fft, d=float(dt))
        h = self._response(transfer_model, omega, float(dt))
        h_inv = self._regularized_inverse(h)
        taps = np.fft.ifft(h_inv).real[: self.n_taps]
        return FIRDistortion(taps=taps, dt=float(dt))

    def _design_iir_inverse(
        self,
        transfer_model: DistortionModel,
        dt: float | None,
    ) -> IIRDistortion | None:
        if isinstance(transfer_model, IIRDistortion):
            return IIRDistortion(
                b=transfer_model.a.copy(),
                a=transfer_model.b.copy(),
                dt=transfer_model.dt,
            )
        if isinstance(transfer_model, FIRDistortion):
            roots = np.roots(transfer_model.taps)
            if len(roots) == 0 or np.all(np.abs(roots) < 1.0):
                return IIRDistortion(
                    b=np.array([1.0]),
                    a=transfer_model.taps.copy(),
                    dt=transfer_model.dt,
                )
            return None
        if dt is not None and hasattr(transfer_model, "iir_coefficients"):
            b, a = transfer_model.iir_coefficients(float(dt))  # type: ignore[attr-defined]
            if not np.all(np.isfinite(b)) or not np.all(np.isfinite(a)):
                return None
            inverse_roots = np.roots(b)
            if len(inverse_roots) and np.any(np.abs(inverse_roots) >= 1.0):
                return None
            return IIRDistortion(b=a, a=b, dt=float(dt))
        return None


__all__ = ["PredistortionDesigner"]
