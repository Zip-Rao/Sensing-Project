"""Control-chain distortion models.

The models in this module describe the linear response between an AWG output
and the flux waveform that reaches the chip. Time is measured in ns throughout
the project, so angular frequencies are rad/ns.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from scipy.signal import lfilter

from sqc.control.waveform import Waveform


def _uniform_dt(t_list: np.ndarray) -> float:
    """Return the grid spacing, rejecting non-uniform sample grids."""
    t = np.asarray(t_list, dtype=float)
    if t.ndim != 1 or len(t) < 2:
        raise ValueError("DistortionModel requires at least two time samples")
    diffs = np.diff(t)
    dt = float(np.median(diffs))
    if dt <= 0 or not np.allclose(diffs, dt, rtol=1e-5, atol=1e-12):
        raise ValueError("DistortionModel requires a uniformly increasing t_list")
    return dt


def _check_dt(actual: float, expected: float) -> None:
    if not np.isclose(actual, expected, rtol=1e-5, atol=1e-12):
        raise ValueError(
            f"waveform dt={actual:.6g} does not match model dt={expected:.6g}"
        )


def _polyval_zm1(coeffs: np.ndarray, omega: np.ndarray, dt: float) -> np.ndarray:
    """Evaluate sum_k coeffs[k] z^-k with z=exp(i omega dt)."""
    coeffs = np.asarray(coeffs, dtype=float)
    powers = np.exp(-1j * np.outer(np.asarray(omega, dtype=float) * dt, np.arange(len(coeffs))))
    return powers @ coeffs


def discrete_iir_frequency_response(
    b: np.ndarray,
    a: np.ndarray,
    omega: np.ndarray,
    dt: float,
) -> np.ndarray:
    """Evaluate a discrete-time IIR response at continuous angular frequency."""
    numerator = _polyval_zm1(np.asarray(b, dtype=float), omega, dt)
    denominator = _polyval_zm1(np.asarray(a, dtype=float), omega, dt)
    return numerator / denominator


class DistortionModel(ABC):
    """Base class for AWG-to-chip distortion models."""

    @abstractmethod
    def apply(self, waveform: Waveform) -> Waveform:
        """Apply the distortion model to a waveform."""

    @abstractmethod
    def step_response(self, t: np.ndarray) -> np.ndarray:
        """Return the model step response."""

    @abstractmethod
    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        """Return the model impulse response."""

    @abstractmethod
    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """Return complex frequency response."""


@dataclass
class SingleExponentialDistortion(DistortionModel):
    """Single slow-tail control-line response.

    The continuous model is

        H(s) = (1 - amplitude) + amplitude / (1 + s tau)

    so the DC gain is one and a positive amplitude creates a slow settling tail.
    """

    amplitude: float
    tau: float

    def __post_init__(self) -> None:
        self.amplitude = float(self.amplitude)
        self.tau = float(self.tau)
        if not np.isfinite(self.amplitude):
            raise ValueError("amplitude must be finite")
        if self.tau <= 0 or not np.isfinite(self.tau):
            raise ValueError("tau must be positive and finite")

    def iir_coefficients(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        """Return discrete IIR coefficients for the sampled model."""
        alpha = float(np.exp(-dt / self.tau))
        b = np.array(
            [1.0 - self.amplitude * alpha, -(1.0 - self.amplitude) * alpha],
            dtype=float,
        )
        a = np.array([1.0, -alpha], dtype=float)
        return b, a

    def apply(self, waveform: Waveform) -> Waveform:
        dt = _uniform_dt(waveform.t_list)
        b, a = self.iir_coefficients(dt)
        out = lfilter(b, a, waveform.samples)
        return Waveform(
            t_list=waveform.t_list.copy(),
            samples=np.asarray(out, dtype=float),
            metadata={**waveform.metadata, "distortion": "single_exponential"},
        )

    def step_response(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        if len(t) < 2:
            return 1.0 - self.amplitude * np.exp(-np.maximum(t, 0.0) / self.tau)
        dt = _uniform_dt(t)
        b, a = self.iir_coefficients(dt)
        return lfilter(b, a, np.ones_like(t, dtype=float))

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        impulse = np.zeros_like(t, dtype=float)
        if len(impulse) == 0:
            return impulse
        impulse[0] = 1.0
        dt = _uniform_dt(t)
        b, a = self.iir_coefficients(dt)
        return lfilter(b, a, impulse)

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        omega = np.asarray(omega, dtype=float)
        return (1.0 - self.amplitude) + self.amplitude / (
            1.0 + 1j * omega * self.tau
        )

    def discrete_frequency_response(
        self, omega: np.ndarray, dt: float
    ) -> np.ndarray:
        b, a = self.iir_coefficients(dt)
        return discrete_iir_frequency_response(b, a, omega, dt)


@dataclass
class MultiExponentialDistortion(DistortionModel):
    """Parallel sum of exponential control-line tails with unit DC gain."""

    amplitudes: np.ndarray
    taus: np.ndarray

    def __post_init__(self) -> None:
        self.amplitudes = np.asarray(self.amplitudes, dtype=float)
        self.taus = np.asarray(self.taus, dtype=float)
        if self.amplitudes.ndim != 1 or self.taus.ndim != 1:
            raise ValueError("amplitudes and taus must be one-dimensional")
        if self.amplitudes.shape != self.taus.shape:
            raise ValueError("amplitudes and taus must have the same shape")
        if len(self.amplitudes) == 0:
            raise ValueError("at least one exponential component is required")
        if not np.all(np.isfinite(self.amplitudes)):
            raise ValueError("amplitudes must be finite")
        if np.any(self.taus <= 0) or not np.all(np.isfinite(self.taus)):
            raise ValueError("taus must be positive and finite")

    def iir_coefficients(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        """Return one equivalent IIR filter for the parallel branches."""
        alphas = np.exp(-dt / self.taus)
        betas = 1.0 - alphas
        denominator = np.array([1.0], dtype=float)
        for alpha in alphas:
            denominator = np.convolve(denominator, np.array([1.0, -alpha]))

        direct = 1.0 - float(np.sum(self.amplitudes))
        numerator = direct * denominator
        for k, (amp, beta) in enumerate(zip(self.amplitudes, betas)):
            branch = np.array([1.0], dtype=float)
            for j, alpha in enumerate(alphas):
                if j != k:
                    branch = np.convolve(branch, np.array([1.0, -alpha]))
            padded = np.pad(branch, (0, len(denominator) - len(branch)))
            numerator = numerator + amp * beta * padded
        return numerator, denominator

    def apply(self, waveform: Waveform) -> Waveform:
        dt = _uniform_dt(waveform.t_list)
        b, a = self.iir_coefficients(dt)
        out = lfilter(b, a, waveform.samples)
        return Waveform(
            t_list=waveform.t_list.copy(),
            samples=np.asarray(out, dtype=float),
            metadata={**waveform.metadata, "distortion": "multi_exponential"},
        )

    def step_response(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        dt = _uniform_dt(t)
        b, a = self.iir_coefficients(dt)
        return lfilter(b, a, np.ones_like(t, dtype=float))

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        impulse = np.zeros_like(t, dtype=float)
        if len(impulse) == 0:
            return impulse
        impulse[0] = 1.0
        dt = _uniform_dt(t)
        b, a = self.iir_coefficients(dt)
        return lfilter(b, a, impulse)

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        omega = np.asarray(omega, dtype=float)
        h = np.ones_like(omega, dtype=complex) * (1.0 - np.sum(self.amplitudes))
        for amp, tau in zip(self.amplitudes, self.taus):
            h += amp / (1.0 + 1j * omega * tau)
        return h

    def discrete_frequency_response(
        self, omega: np.ndarray, dt: float
    ) -> np.ndarray:
        b, a = self.iir_coefficients(dt)
        return discrete_iir_frequency_response(b, a, omega, dt)


@dataclass
class FIRDistortion(DistortionModel):
    """Generic finite-impulse-response distortion."""

    taps: np.ndarray
    dt: float

    def __post_init__(self) -> None:
        self.taps = np.asarray(self.taps, dtype=float)
        self.dt = float(self.dt)
        if self.taps.ndim != 1 or len(self.taps) == 0:
            raise ValueError("taps must be a non-empty one-dimensional array")
        if self.dt <= 0 or not np.isfinite(self.dt):
            raise ValueError("dt must be positive and finite")

    def apply(self, waveform: Waveform) -> Waveform:
        _check_dt(_uniform_dt(waveform.t_list), self.dt)
        out = lfilter(self.taps, np.array([1.0]), waveform.samples)
        return Waveform(
            t_list=waveform.t_list.copy(),
            samples=np.asarray(out, dtype=float),
            metadata={**waveform.metadata, "distortion": "fir"},
        )

    def step_response(self, t: np.ndarray) -> np.ndarray:
        _check_dt(_uniform_dt(t), self.dt)
        return lfilter(self.taps, np.array([1.0]), np.ones_like(t, dtype=float))

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        _check_dt(_uniform_dt(t), self.dt)
        impulse = np.zeros_like(t, dtype=float)
        if len(impulse):
            impulse[0] = 1.0
        return lfilter(self.taps, np.array([1.0]), impulse)

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        return discrete_iir_frequency_response(self.taps, np.array([1.0]), omega, self.dt)


@dataclass
class IIRDistortion(DistortionModel):
    """Generic infinite-impulse-response distortion."""

    b: np.ndarray
    a: np.ndarray
    dt: float

    def __post_init__(self) -> None:
        self.b = np.asarray(self.b, dtype=float)
        self.a = np.asarray(self.a, dtype=float)
        self.dt = float(self.dt)
        if self.b.ndim != 1 or len(self.b) == 0:
            raise ValueError("b must be a non-empty one-dimensional array")
        if self.a.ndim != 1 or len(self.a) == 0:
            raise ValueError("a must be a non-empty one-dimensional array")
        if self.a[0] == 0:
            raise ValueError("a[0] must be non-zero")
        if self.dt <= 0 or not np.isfinite(self.dt):
            raise ValueError("dt must be positive and finite")
        if not np.isclose(self.a[0], 1.0):
            scale = self.a[0]
            self.b = self.b / scale
            self.a = self.a / scale

    def apply(self, waveform: Waveform) -> Waveform:
        _check_dt(_uniform_dt(waveform.t_list), self.dt)
        out = lfilter(self.b, self.a, waveform.samples)
        return Waveform(
            t_list=waveform.t_list.copy(),
            samples=np.asarray(out, dtype=float),
            metadata={**waveform.metadata, "distortion": "iir"},
        )

    def step_response(self, t: np.ndarray) -> np.ndarray:
        _check_dt(_uniform_dt(t), self.dt)
        return lfilter(self.b, self.a, np.ones_like(t, dtype=float))

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        _check_dt(_uniform_dt(t), self.dt)
        impulse = np.zeros_like(t, dtype=float)
        if len(impulse):
            impulse[0] = 1.0
        return lfilter(self.b, self.a, impulse)

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        return discrete_iir_frequency_response(self.b, self.a, omega, self.dt)

    def discrete_frequency_response(
        self, omega: np.ndarray, dt: float
    ) -> np.ndarray:
        _check_dt(float(dt), self.dt)
        return self.frequency_response(omega)


@dataclass
class CustomTransferDistortion(DistortionModel):
    """Frequency-domain distortion defined on a fixed angular-frequency grid."""

    omega: np.ndarray
    H: np.ndarray

    def __post_init__(self) -> None:
        omega = np.asarray(self.omega, dtype=float)
        h = np.asarray(self.H, dtype=complex)
        if omega.ndim != 1 or h.ndim != 1 or omega.shape != h.shape:
            raise ValueError("omega and H must be one-dimensional arrays of equal shape")
        if len(omega) < 2:
            raise ValueError("at least two frequency samples are required")
        if not np.all(np.isfinite(omega)) or not np.all(np.isfinite(h.real + h.imag)):
            raise ValueError("omega and H must be finite")
        order = np.argsort(omega)
        omega = omega[order]
        h = h[order]
        unique = np.diff(omega, prepend=-np.inf) > 0
        self.omega = omega[unique]
        self.H = h[unique]
        if len(self.omega) < 2:
            raise ValueError("omega grid must contain at least two unique points")

    def apply(self, waveform: Waveform) -> Waveform:
        dt = _uniform_dt(waveform.t_list)
        n = len(waveform.samples)
        n_fft = 1 << max(1, (2 * n - 1).bit_length())
        omega_fft = 2.0 * np.pi * np.fft.fftfreq(n_fft, d=dt)
        filtered = np.fft.ifft(
            np.fft.fft(waveform.samples, n_fft) * self.frequency_response(omega_fft)
        ).real[:n]
        return Waveform(
            t_list=waveform.t_list.copy(),
            samples=np.asarray(filtered, dtype=float),
            metadata={**waveform.metadata, "distortion": "custom_transfer"},
        )

    def step_response(self, t: np.ndarray) -> np.ndarray:
        return self.apply(Waveform(t_list=np.asarray(t, dtype=float), samples=np.ones_like(t))).samples

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        impulse = np.zeros_like(t, dtype=float)
        if len(impulse):
            impulse[0] = 1.0
        return self.apply(Waveform(t_list=np.asarray(t, dtype=float), samples=impulse)).samples

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        omega = np.asarray(omega, dtype=float)
        real = np.interp(omega, self.omega, self.H.real)
        imag = np.interp(omega, self.omega, self.H.imag)
        return real + 1j * imag

    def discrete_frequency_response(
        self, omega: np.ndarray, dt: float
    ) -> np.ndarray:
        return self.frequency_response(omega)


__all__ = [
    "DistortionModel",
    "SingleExponentialDistortion",
    "MultiExponentialDistortion",
    "FIRDistortion",
    "IIRDistortion",
    "CustomTransferDistortion",
    "discrete_iir_frequency_response",
]
