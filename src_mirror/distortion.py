"""src_mirror.distortion — DistortionModel implementations (Track B 1.3).

Models the transfer function from AWG voltage to on-chip flux/signal.
Each subclass implements apply(), step_response(), impulse_response(),
and frequency_response().

These are developed in src_mirror/ per R1 (never modify src/).
P4 will internalize them into sqc/hardware/distortion.py.

Physical reference: Gao 2021 §III.D (control-line transfer functions),
§V.E (Cryoscope calibration of distortion tails).
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

import numpy as np
from scipy.signal import lfilter, freqz, bilinear


# ---------------------------------------------------------------------------
# Base class (compatible with sqc.hardware.distortion.DistortionModel ABC)
# ---------------------------------------------------------------------------

class DistortionModel(ABC):
    """Abstract model of AWG → on-chip signal distortion.

    Subclasses implement linear time-invariant (LTI) distortion
    via convolution kernels, IIR/FIR filters, or frequency-domain
    transfer functions.
    """

    @abstractmethod
    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        """Apply distortion to uniformly-sampled waveform data.

        Parameters
        ----------
        waveform : np.ndarray, shape (N,)
            Input samples.
        dt : float
            Sample spacing (ns).

        Returns
        -------
        np.ndarray, shape (N,)
            Distorted output samples.
        """
        ...

    @abstractmethod
    def step_response(self, t: np.ndarray) -> np.ndarray:
        """Step response s(t) = output when input is Θ(t)."""
        ...

    @abstractmethod
    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        """Impulse response h(t) (continuous-time kernel)."""
        ...

    @abstractmethod
    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """Complex frequency response H(ω).

        Parameters
        ----------
        omega : np.ndarray
            Angular frequencies (rad/ns).

        Returns
        -------
        np.ndarray, complex
            H(ω) at each frequency.
        """
        ...

    # Convenience wrappers for Waveform objects
    def apply_to_waveform(self, wf):
        """Apply distortion to a sqc.control.waveform.Waveform.

        Returns a new Waveform with distorted samples.
        """
        from sqc.control.waveform import Waveform

        dt = float(wf.t_list[1] - wf.t_list[0])
        out = self.apply(wf.samples, dt)
        return Waveform(
            t_list=wf.t_list.copy(),
            samples=out,
            metadata={**wf.metadata, "distortion": type(self).__name__},
        )


# ---------------------------------------------------------------------------
# Single-exponential distortion (most common flux-line tail model)
# ---------------------------------------------------------------------------

@dataclass
class SingleExponentialDistortion(DistortionModel):
    """Single-exponential charge-redistribution tail.

    Continuous-time model:
        h(t) = (1 - amp)·δ(t) + (amp / τ)·exp(-t / τ)·Θ(t)

    where Θ(t) is the Heaviside step.  The DC gain is always 1 (the
    delta-function term compensates the tail so the step response settles
    to unity).

    This is the standard model for slow flux-line distortion observed
    in Cryoscope calibration (Gao 2021 §V.E).

    Parameters
    ----------
    amplitude : float
        Tail amplitude (dimensionless fraction, typically 0.001–0.1).
    tau : float
        Tail time constant (ns, typically 10–1000 ns).
    """

    amplitude: float = 0.01
    tau: float = 100.0

    # ---- discrete-time IIR coefficients (bilinear transform) ----
    def _iir_coeffs(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        """Return (b, a) for the discrete-time filter.

        Derivation: the continuous-time transfer function is
            H(s) = (1 - amp) + amp / (1 + s·τ) = 1 - amp·s·τ / (1 + s·τ)

        Bilinear transform s → 2/dt · (1 - z⁻¹)/(1 + z⁻¹) gives
        a first-order IIR filter.
        """
        if dt <= 0:
            raise ValueError(f"dt must be positive, got {dt}")
        # Continuous-time coefficients (numerator, denominator in s)
        # H(s) = (1-amp)*(1+sτ) + amp  /  (1+sτ)
        #      = (1 + sτ - amp·sτ) / (1 + sτ)
        #      = (1 + sτ(1-amp)) / (1 + sτ)
        b_cont = [self.tau * (1.0 - self.amplitude), 1.0]  # s¹, s⁰
        a_cont = [self.tau, 1.0]                             # s¹, s⁰
        b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
        return b, a

    # ---- DistortionModel interface ----

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        b, a = self._iir_coeffs(dt)
        return lfilter(b, a, waveform)

    def step_response(self, t: np.ndarray) -> np.ndarray:
        """s(t) = 1 - amp·exp(-t/τ)  for t ≥ 0."""
        out = np.ones_like(t, dtype=float)
        mask = t >= 0
        out[mask] = 1.0 - self.amplitude * np.exp(-t[mask] / self.tau)
        return out

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        """h(t) = (1-amp)·δ(t) + (amp/τ)·exp(-t/τ)·Θ(t).

        The delta is approximated by placing its weight in the first bin.
        """
        h = np.zeros_like(t, dtype=float)
        mask = t >= 0
        h[mask] = (self.amplitude / self.tau) * np.exp(-t[mask] / self.tau)
        # place delta mass in first sample
        if len(t) > 0:
            dt_approx = t[1] - t[0] if len(t) > 1 else 1.0
            h[0] += (1.0 - self.amplitude) / dt_approx
        return h

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """H(ω) = (1 - amp) + amp / (1 + j·ω·τ)."""
        return (1.0 - self.amplitude) + self.amplitude / (1.0 + 1j * omega * self.tau)


# ---------------------------------------------------------------------------
# Multi-exponential distortion
# ---------------------------------------------------------------------------

@dataclass
class MultiExponentialDistortion(DistortionModel):
    """Sum of K single-exponential tails.

    H(s) = 1 + Σₖ ampₖ · (1/(1 + s·τₖ) - 1)
         = 1 - Σₖ ampₖ + Σₖ ampₖ/(1 + s·τₖ)

    The DC gain is 1 by construction.

    Parameters
    ----------
    amplitudes : np.ndarray, shape (K,)
        Tail amplitudes.
    taus : np.ndarray, shape (K,)
        Tail time constants (ns).
    """

    amplitudes: np.ndarray = field(default_factory=lambda: np.array([0.01]))
    taus: np.ndarray = field(default_factory=lambda: np.array([100.0]))

    def __post_init__(self) -> None:
        self.amplitudes = np.asarray(self.amplitudes, dtype=float)
        self.taus = np.asarray(self.taus, dtype=float)
        if self.amplitudes.shape != self.taus.shape:
            raise ValueError("amplitudes and taus must have the same length")
        # pre-compute global scale factor
        self._dc_gain = 1.0 - np.sum(self.amplitudes)

    @property
    def n_components(self) -> int:
        return len(self.amplitudes)

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        """Parallel sum: direct path + K first-order IIR tails.

        H(s) = (1 - Σamp_k) + Σ_k amp_k/(1 + s·τ_k)

        Each tail is discretised via bilinear transform independently,
        then all contributions are summed.
        """
        if self.n_components == 0:
            return waveform.copy()
        # Direct path (DC gain correction so total DC gain = 1)
        out = self._dc_gain * waveform.astype(float)
        # Add each exponential tail contribution
        for amp, tau in zip(self.amplitudes, self.taus):
            # Continuous-time: H_k(s) = amp / (1 + s·τ)
            b_cont = [0.0, amp]   # numerator in s: amp
            a_cont = [tau, 1.0]   # denominator in s: τ·s + 1
            b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
            out += lfilter(b, a, waveform)
        return out

    def step_response(self, t: np.ndarray) -> np.ndarray:
        out = np.ones_like(t, dtype=float)
        mask = t >= 0
        for amp, tau in zip(self.amplitudes, self.taus):
            out[mask] -= amp * np.exp(-t[mask] / tau)
        return out

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        h = np.zeros_like(t, dtype=float)
        mask = t >= 0
        for amp, tau in zip(self.amplitudes, self.taus):
            h[mask] += (amp / tau) * np.exp(-t[mask] / tau)
        if len(t) > 0:
            dt_approx = t[1] - t[0] if len(t) > 1 else 1.0
            h[0] += self._dc_gain / dt_approx
        return h

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        H = np.full_like(omega, self._dc_gain, dtype=complex)
        for amp, tau in zip(self.amplitudes, self.taus):
            H += amp / (1.0 + 1j * omega * tau)
        return H


# ---------------------------------------------------------------------------
# FIR filter distortion
# ---------------------------------------------------------------------------

@dataclass
class FIRDistortion(DistortionModel):
    """Finite impulse response (FIR) distortion model.

    y[n] = Σₖ b[k] · x[n - k]

    Parameters
    ----------
    taps : np.ndarray, shape (M,)
        FIR filter coefficients b[0..M-1].
    """

    taps: np.ndarray = field(default_factory=lambda: np.array([1.0]))

    def __post_init__(self) -> None:
        self.taps = np.asarray(self.taps, dtype=float)

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        return lfilter(self.taps, [1.0], waveform)

    def step_response(self, t: np.ndarray) -> np.ndarray:
        """Cumulative sum of impulse response."""
        h = self.impulse_response(t)
        return np.cumsum(h) * (t[1] - t[0]) if len(t) > 1 else h.copy()

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        """Place FIR taps at the first M sample positions."""
        h = np.zeros_like(t, dtype=float)
        if len(t) > 0:
            dt = t[1] - t[0] if len(t) > 1 else 1.0
            n_place = min(len(self.taps), len(h))
            h[:n_place] = self.taps[:n_place] / dt  # continuous-time scaling
        return h

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """H(ω) = Σₖ b[k] · exp(-j·ω·k·dt).

        Uses nominal dt=1 for normalised frequency.
        """
        return np.polyval(self.taps[::-1], np.exp(-1j * omega))


# ---------------------------------------------------------------------------
# IIR filter distortion
# ---------------------------------------------------------------------------

@dataclass
class IIRDistortion(DistortionModel):
    """Infinite impulse response (IIR) distortion model.

    y[n] = Σₖ b[k]·x[n-k] - Σₖ a[k]·y[n-k],  a[0] = 1

    Parameters
    ----------
    b : np.ndarray, shape (M,)
        Feed-forward coefficients.
    a : np.ndarray, shape (M,)
        Feed-back coefficients (a[0] must be 1).
    """

    b: np.ndarray = field(default_factory=lambda: np.array([1.0]))
    a: np.ndarray = field(default_factory=lambda: np.array([1.0]))

    def __post_init__(self) -> None:
        self.b = np.asarray(self.b, dtype=float)
        self.a = np.asarray(self.a, dtype=float)
        if self.a[0] != 1.0:
            self.b = self.b / self.a[0]
            self.a = self.a / self.a[0]

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        return lfilter(self.b, self.a, waveform)

    def step_response(self, t: np.ndarray) -> np.ndarray:
        if len(t) < 2:
            return np.ones_like(t)
        dt = t[1] - t[0]
        n = len(t)
        step_in = np.ones(n)
        return lfilter(self.b, self.a, step_in)

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        if len(t) < 2:
            return np.zeros_like(t)
        dt = t[1] - t[0]
        n = len(t)
        imp_in = np.zeros(n)
        imp_in[0] = 1.0 / dt
        return lfilter(self.b, self.a, imp_in)

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """H(ω) evaluated via polynomial evaluation on unit circle."""
        H = np.zeros(len(omega), dtype=complex)
        z = np.exp(1j * omega)
        for i, bk in enumerate(self.b):
            H += bk * z ** (-i)
        denom = np.zeros(len(omega), dtype=complex)
        for i, ak in enumerate(self.a):
            denom += ak * z ** (-i)
        return H / denom


# ---------------------------------------------------------------------------
# Custom frequency-domain transfer function
# ---------------------------------------------------------------------------

@dataclass
class CustomTransferDistortion(DistortionModel):
    """User-supplied complex frequency-domain transfer function H(ω).

    Parameters
    ----------
    omega_grid : np.ndarray, shape (F,)
        Frequency grid (rad/ns).
    H_grid : np.ndarray, shape (F,), complex
        H(ω) evaluated on omega_grid.
    """

    omega_grid: np.ndarray = field(default_factory=lambda: np.linspace(-np.pi, np.pi, 1024))
    H_grid: np.ndarray = field(default_factory=lambda: np.ones(1024, dtype=complex))

    def __post_init__(self) -> None:
        self.omega_grid = np.asarray(self.omega_grid, dtype=float)
        self.H_grid = np.asarray(self.H_grid, dtype=complex)

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        """FFT → multiply by H(ω) → IFFT."""
        n = len(waveform)
        omega_wf = 2.0 * np.pi * np.fft.fftfreq(n, d=dt)
        # Interpolate H onto FFT grid
        H_interp = self._interpolate_H(omega_wf)
        X = np.fft.fft(waveform)
        Y = H_interp * X
        return np.fft.ifft(Y).real

    def _interpolate_H(self, omega_target: np.ndarray) -> np.ndarray:
        """Linear interpolation of H(ω) onto target frequency grid."""
        re = np.interp(omega_target, self.omega_grid, self.H_grid.real,
                       left=self.H_grid.real[0], right=self.H_grid.real[-1])
        im = np.interp(omega_target, self.omega_grid, self.H_grid.imag,
                       left=self.H_grid.imag[0], right=self.H_grid.imag[-1])
        return re + 1j * im

    def step_response(self, t: np.ndarray) -> np.ndarray:
        if len(t) < 2:
            return np.ones_like(t)
        dt = t[1] - t[0]
        step_in = np.where(t >= 0, 1.0, 0.0)
        return self.apply(step_in, dt)

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        if len(t) < 2:
            return np.zeros_like(t)
        dt = t[1] - t[0]
        imp_in = np.zeros_like(t)
        imp_in[0] = 1.0 / dt
        return self.apply(imp_in, dt)

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        return self._interpolate_H(omega)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = [
    "DistortionModel",
    "SingleExponentialDistortion",
    "MultiExponentialDistortion",
    "FIRDistortion",
    "IIRDistortion",
    "CustomTransferDistortion",
]
