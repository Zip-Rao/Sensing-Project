"""sqc.hardware.distortion — DistortionModel ABC and subclasses.

Models the transfer function from AWG voltage to on-chip flux/signal.
Each subclass implements linear time-invariant (LTI) distortion
via convolution kernels, IIR/FIR filters, or frequency-domain
transfer functions.

Physical reference: Gao 2021 §III.D (control-line transfer functions),
§V.E (Cryoscope calibration of distortion tails).

All subclasses support a ``smooth`` flag:
  - smooth=False (default): includes the direct (delta-function) path,
    preserving input discontinuities -> broad spectrum.
  - smooth=True: pure-lowpass behaviour, no direct path. An ideal step
    input produces a smooth ramp (e.g. 1 - exp(-t/tau)) with no jump.

ABC and 6 subclasses:
  - SingleExponentialDistortion  (most common flux-line tail model)
  - MultiExponentialDistortion   (sum of K single-exponential tails)
  - FIRDistortion                (finite impulse response)
  - IIRDistortion                (infinite impulse response)
  - CustomTransferDistortion     (user-supplied H(omega) on frequency grid)
  - CascadeDistortion            (series cascade, Rol 2020 Sec.IV)
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

import numpy as np
from scipy.signal import lfilter, bilinear


# ---------------------------------------------------------------------------
# DistortionModel ABC
# ---------------------------------------------------------------------------

class DistortionModel(ABC):
    """Abstract model of control-line distortion.

    Subclasses implement linear time-invariant (LTI) distortion.
    The low-level interface uses numpy arrays; the convenience method
    `apply_to_waveform()` works with sqc.control.waveform.Waveform.
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
        """Step response s(t) = output when input is Theta(t)."""
        ...

    @abstractmethod
    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        """Impulse response h(t) (continuous-time kernel)."""
        ...

    @abstractmethod
    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """Complex frequency response H(omega).

        Parameters
        ----------
        omega : np.ndarray
            Angular frequencies (rad/ns).

        Returns
        -------
        np.ndarray, complex
            H(omega) at each frequency.
        """
        ...

    # Convenience wrapper for Waveform objects
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

    def apply_to_signal(self, signal):
        """Apply distortion to a FluxSignal, preserving its type.

        Uses ``type=8`` (user-defined raw samples), which is the
        canonical pathway for reconstructed / distorted signals
        throughout sqc (see delay_ramsey.py:138, cryoscope.py).

        Parameters
        ----------
        signal : FluxSignal

        Returns
        -------
        FluxSignal (type=8) with distorted samples.
        """
        from sqc.control.flux_signal import FluxSignal

        dt = float(signal.t_list[1] - signal.t_list[0])
        out = self.apply(np.asarray(signal.samples), dt)
        return FluxSignal(type=8, t_list=signal.t_list.copy(), signal=out)


# ---------------------------------------------------------------------------
# Single-exponential distortion (most common flux-line tail model)
# ---------------------------------------------------------------------------

@dataclass
class SingleExponentialDistortion(DistortionModel):
    """Single-exponential charge-redistribution tail.

    Continuous-time model (smooth=False, default):
        h(t) = (1 - amp)·delta(t) + (amp / tau)·exp(-t / tau)·Theta(t)
        s(t) = 1 - amp·exp(-t/tau)            (has jump at t=0)

    Continuous-time model (smooth=True):
        h(t) = (1 / tau)·exp(-t / tau)·Theta(t)
        s(t) = 1 - exp(-t/tau)                 (smooth ramp from 0)
        H(s) = 1 / (1 + s·tau)

    where Theta(t) is the Heaviside step.  The DC gain is always 1.

    This is the standard model for slow flux-line distortion observed
    in Cryoscope calibration (Gao 2021 §V.E).

    Parameters
    ----------
    amplitude : float
        Tail amplitude (dimensionless fraction, typically 0.001–0.1).
        In smooth mode, amplitude is ignored (the tail weight is fixed
        at unity to give DC gain = 1).
    tau : float
        Tail time constant (ns, typically 10–1000 ns).
    smooth : bool
        If True, use pure-lowpass model with no delta-function direct path.
        Default False (backward-compatible).
    """

    amplitude: float = 0.01
    tau: float = 100.0
    smooth: bool = False

    # ---- discrete-time IIR coefficients (bilinear transform) ----
    def _iir_coeffs(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        """Return (b, a) for the discrete-time filter.

        smooth=False: H(s) = (1 - amp) + amp / (1 + s·tau)
                      = (1 + s·tau·(1-amp)) / (1 + s·tau)

        smooth=True:  H(s) = 1 / (1 + s·tau)

        Bilinear transform s -> 2/dt · (1 - z⁻¹)/(1 + z⁻¹) gives
        a first-order IIR filter.
        """
        if dt <= 0:
            raise ValueError(f"dt must be positive, got {dt}")
        if self.smooth:
            b_cont = [0.0, 1.0]
        else:
            b_cont = [self.tau * (1.0 - self.amplitude), 1.0]
        a_cont = [self.tau, 1.0]
        b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
        return b, a

    # ---- DistortionModel interface ----

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        b, a = self._iir_coeffs(dt)
        return lfilter(b, a, waveform)

    def step_response(self, t: np.ndarray) -> np.ndarray:
        """s(t) for t >= 0. With smooth=True: 1 - exp(-t/tau)."""
        out = np.ones_like(t, dtype=float)
        mask = t >= 0
        if self.smooth:
            out[mask] = 1.0 - np.exp(-t[mask] / self.tau)
        else:
            out[mask] = 1.0 - self.amplitude * np.exp(-t[mask] / self.tau)
        return out

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        """h(t). With smooth=True: (1/tau)·exp(-t/tau) (no delta)."""
        h = np.zeros_like(t, dtype=float)
        mask = t >= 0
        if self.smooth:
            h[mask] = (1.0 / self.tau) * np.exp(-t[mask] / self.tau)
        else:
            h[mask] = (self.amplitude / self.tau) * np.exp(-t[mask] / self.tau)
            # place delta mass in first sample
            if len(t) > 0:
                dt_approx = t[1] - t[0] if len(t) > 1 else 1.0
                h[0] += (1.0 - self.amplitude) / dt_approx
        return h

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """H(omega). With smooth=True: 1 / (1 + j·omega·tau)."""
        if self.smooth:
            return 1.0 / (1.0 + 1j * omega * self.tau)
        return (1.0 - self.amplitude) + self.amplitude / (1.0 + 1j * omega * self.tau)


    def design_inverse(
        self, dt: float,
        formula: str = "bilinear",
    ) -> "IIRDistortion":
        """Design an IIR inverse filter that compensates this distortion.

        The forward model step response is s(t) = 1 - A*exp(-t/tau).
        The inverse filter, when cascaded with the forward model,
        cancels the exponential tail.

        Parameters
        ----------
        dt : float
            Sample spacing (ns).
        formula : str
            ``"bilinear"`` -- bilinear transform (continuous-time inverse).
            ``"rol2020"`` -- direct z-domain pole-zero design.

        Returns
        -------
        IIRDistortion
            Inverse filter.
        """
        import warnings

        amp = self.amplitude
        if abs(amp) > 0.5:
            warnings.warn(
                f"|amplitude|={abs(amp):.3f} > 0.5; "
                f"Rol 2020 verified range is |A| <= 0.1. "
                f"Results may degrade.",
                stacklevel=2,
            )

        if formula == "bilinear":
            return self._design_inverse_bilinear(dt)
        elif formula == "rol2020":
            return self._design_inverse_rol2020(dt)
        raise ValueError(
            f"Unknown formula '{formula}'. Valid: 'bilinear', 'rol2020'."
        )

    def _design_inverse_bilinear(self, dt: float) -> "IIRDistortion":
        """Bilinear-transform inverse: H_inv(s) = (1+s*tau)/(1+s*tau*(1-A))."""
        tau = self.tau
        amp = self.amplitude
        # As A -> 1 the forward DC gain -> 0, so a stable inverse does not
        # exist: a_cont's leading coefficient tau*(1-A) -> 0 collapses the
        # denominator order and produces an unbounded high-frequency inverse.
        # Degrade to identity, matching _design_inverse_rol2020's guard.
        if abs(1.0 - amp) < 1e-15:
            return IIRDistortion(
                b_coeffs=np.array([1.0]), a_coeffs=np.array([1.0]),
            )
        b_cont = [tau, 1.0]
        a_cont = [tau * (1.0 - amp), 1.0]
        b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
        return IIRDistortion(b_coeffs=b, a_coeffs=a)

    def _design_inverse_rol2020(self, dt: float) -> "IIRDistortion":
        """Direct z-domain inverse (Rol 2020 Eq. S22).

        For step response s(t) = 1 - A*exp(-t/tau), the discretised
        system has pole alpha = exp(-dt/tau).  The inverse filter places
        a zero at alpha and a pole at (alpha-A)/(1-A).
        """
        amp = self.amplitude
        tau = self.tau

        if abs(1.0 - amp) < 1e-15:
            return IIRDistortion(
                b_coeffs=np.array([1.0]), a_coeffs=np.array([1.0]),
            )

        alpha = np.exp(-dt / max(tau, 1e-30))
        b0 = 1.0 / (1.0 - amp)
        b1 = -alpha / (1.0 - amp)
        a1 = (amp - alpha) / (1.0 - amp)

        return IIRDistortion(
            b_coeffs=np.array([b0, b1]),
            a_coeffs=np.array([1.0, a1]),
        )


# ---------------------------------------------------------------------------
# Multi-exponential distortion
# ---------------------------------------------------------------------------

@dataclass
class MultiExponentialDistortion(DistortionModel):
    """Sum of K single-exponential tails.

    smooth=False (default):
        H(s) = 1 - sum_k amp_k + sum_k amp_k/(1 + s·tau_k)
        DC gain = 1 (direct path compensates tails).

    smooth=True:
        H(s) = sum_k w_k / (1 + s·tau_k)
        where w_k = amp_k / sum(amp_k).
        DC gain = 1 (via normalisation), no direct path.
        Pure multi-time-constant lowpass — ideal step → smooth ramp.

    Parameters
    ----------
    amplitudes : np.ndarray, shape (K,)
        Tail amplitudes.
    taus : np.ndarray, shape (K,)
        Tail time constants (ns).
    smooth : bool
        If True, use pure-lowpass model with no direct path.
        Default False (backward-compatible).
    """

    amplitudes: np.ndarray = field(default_factory=lambda: np.array([0.01]))
    taus: np.ndarray = field(default_factory=lambda: np.array([100.0]))
    smooth: bool = False

    def __post_init__(self) -> None:
        self.amplitudes = np.asarray(self.amplitudes, dtype=float)
        self.taus = np.asarray(self.taus, dtype=float)
        if self.amplitudes.shape != self.taus.shape:
            raise ValueError("amplitudes and taus must have the same length")
        # pre-compute global scale factor (used in non-smooth mode)
        self._dc_gain = 1.0 - np.sum(self.amplitudes)
        # pre-compute normalised weights (used in smooth mode)
        total = np.sum(self.amplitudes)
        self._weights = self.amplitudes / total if total != 0 else self.amplitudes.copy()

    @property
    def n_components(self) -> int:
        return len(self.amplitudes)

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        """Parallel sum of K first-order IIR tails.

        smooth=False: direct path (1 - sum amp) + sum of tails.
        smooth=True:  weighted sum of pure-lowpass tails, no direct path.

        Each tail is discretised via bilinear transform independently.
        """
        if self.n_components == 0:
            return waveform.copy()
        if self.smooth:
            out = np.zeros_like(waveform, dtype=float)
            for weight, tau in zip(self._weights, self.taus):
                # H_k(s) = weight / (1 + s·tau)
                b_cont = [0.0, weight]
                a_cont = [tau, 1.0]
                b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
                out += lfilter(b, a, waveform)
            return out
        else:
            # Direct path (DC gain correction so total DC gain = 1)
            out = self._dc_gain * waveform.astype(float)
            # Add each exponential tail contribution
            for amp, tau in zip(self.amplitudes, self.taus):
                # Continuous-time: H_k(s) = amp / (1 + s·tau)
                b_cont = [0.0, amp]
                a_cont = [tau, 1.0]
                b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
                out += lfilter(b, a, waveform)
            return out

    def step_response(self, t: np.ndarray) -> np.ndarray:
        out = np.ones_like(t, dtype=float)
        mask = t >= 0
        if self.smooth:
            out[mask] = 0.0
            for w, tau in zip(self._weights, self.taus):
                out[mask] += w * (1.0 - np.exp(-t[mask] / tau))
        else:
            for amp, tau in zip(self.amplitudes, self.taus):
                out[mask] -= amp * np.exp(-t[mask] / tau)
        return out

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        h = np.zeros_like(t, dtype=float)
        mask = t >= 0
        if self.smooth:
            for w, tau in zip(self._weights, self.taus):
                h[mask] += (w / tau) * np.exp(-t[mask] / tau)
        else:
            for amp, tau in zip(self.amplitudes, self.taus):
                h[mask] += (amp / tau) * np.exp(-t[mask] / tau)
            if len(t) > 0:
                dt_approx = t[1] - t[0] if len(t) > 1 else 1.0
                h[0] += self._dc_gain / dt_approx
        return h

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        if self.smooth:
            H = np.zeros_like(omega, dtype=complex)
            for w, tau in zip(self._weights, self.taus):
                H += w / (1.0 + 1j * omega * tau)
            return H
        else:
            H = np.full_like(omega, self._dc_gain, dtype=complex)
            for amp, tau in zip(self.amplitudes, self.taus):
                H += amp / (1.0 + 1j * omega * tau)
            return H


    def design_inverse(
        self,
        dt: float,
        formula: str = "bilinear",
        n_iir_stages: int = 3,
        fir_taps: int = 72,
        fir_threshold_ns: float = 30.0,
    ) -> "CascadeDistortion":
        """Design a cascade inverse for this multi-exponential distortion.

        Each exponential component with tau > *fir_threshold_ns* gets an
        IIR inverse stage (up to *n_iir_stages*).  A residual FIR stage
        corrects fast (< *fir_threshold_ns*) errors.

        This implements the Rol 2020 Sec.IV iterative decomposition:
        IIR_i -> ... -> IIR_N -> FIR.

        Parameters
        ----------
        dt : float
            Sample spacing (ns).
        formula : str
            ``"bilinear"`` or ``"rol2020"``, passed to each component.
        n_iir_stages : int
            Maximum number of IIR stages.
        fir_taps : int
            Number of FIR taps for residual correction (0 = skip).
        fir_threshold_ns : float
            Only components with tau > this value get an IIR stage.

        Returns
        -------
        CascadeDistortion
            Cascade of IIR inverse stages + optional FIR.
        """
        if self.n_components == 0:
            return CascadeDistortion(stages=[])

        # 1. Sort by tau descending
        order = np.argsort(self.taus)[::-1]
        amps_sorted = self.amplitudes[order]
        taus_sorted = self.taus[order]

        # 2. Select components with tau > threshold, capped at n_iir_stages
        sel = taus_sorted > fir_threshold_ns
        if not np.any(sel):
            t_step = np.arange(0, fir_threshold_ns * 4, dt)
            step_in = np.where(t_step >= 0, 1.0, 0.0)
            s_fwd = self.apply(step_in, dt)
            if fir_taps > 0:
                fir = FIRDistortion.design_from_residual(
                    s_fwd, dt, n_taps=fir_taps,
                )
                return CascadeDistortion(stages=[fir])
            return CascadeDistortion(stages=[])

        sel_indices = np.where(sel)[0][:n_iir_stages]
        amps_sel = amps_sorted[sel_indices]
        taus_sel = taus_sorted[sel_indices]

        # 3. Design IIR inverse for each selected component
        iir_stages = []
        for amp_k, tau_k in zip(amps_sel, taus_sel):
            comp = SingleExponentialDistortion(
                amplitude=float(amp_k), tau=float(tau_k),
            )
            iir_stages.append(comp.design_inverse(dt, formula=formula))

        # 4. Compute residual after IIR cascade on unit step
        t_step = np.arange(
            0, float(max(np.max(taus_sorted) * 4, fir_threshold_ns * 2)), dt,
        )
        step_in = np.where(t_step >= 0, 1.0, 0.0)
        s_fwd = self.apply(step_in, dt)
        partial_cascade = CascadeDistortion(stages=iir_stages)
        s_corrected = partial_cascade.apply(s_fwd, dt)

        # 5. FIR residual correction (fast components)
        if fir_taps > 0:
            fir = FIRDistortion.design_from_residual(
                s_corrected, dt, n_taps=fir_taps,
            )
            return CascadeDistortion(stages=iir_stages + [fir])

        return CascadeDistortion(stages=iir_stages)


# ---------------------------------------------------------------------------
# FIR filter distortion
# ---------------------------------------------------------------------------

@dataclass
class FIRDistortion(DistortionModel):
    """Finite impulse response (FIR) distortion model.

    y[n] = sum_k b[k] · x[n - k]

    smooth=False (default): taps applied directly via lfilter.  If b[0] != 0
        there is direct feedthrough and a step input produces a jump.

    smooth=True: the FIR output is cascaded through a first-order lowpass
        with time constant ``smooth_tau`` to remove any sharp edges.

    Parameters
    ----------
    taps : np.ndarray, shape (M,)
        FIR filter coefficients b[0..M-1].
    smooth : bool
        If True, cascade through a smoothing lowpass after the FIR.
        Default False.
    smooth_tau : float
        Time constant (ns) for the smoothing lowpass, used only when
        smooth=True.  Default 1.0.
    """

    taps: np.ndarray = field(default_factory=lambda: np.array([1.0]))
    smooth: bool = False
    smooth_tau: float = 1.0

    def __post_init__(self) -> None:
        self.taps = np.asarray(self.taps, dtype=float)

    @property
    def dt(self) -> float:
        """Nominal dt for FIR taps (for frequency_response)."""
        return 1.0

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        out = lfilter(self.taps, [1.0], waveform)
        if self.smooth:
            b_cont = [0.0, 1.0]
            a_cont = [self.smooth_tau, 1.0]
            b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
            out = lfilter(b, a, out)
        return out

    def step_response(self, t: np.ndarray) -> np.ndarray:
        """Cumulative sum of impulse response."""
        if self.smooth:
            if len(t) < 2:
                return np.ones_like(t)
            dt_val = t[1] - t[0]
            step_in = np.where(t >= 0, 1.0, 0.0)
            return self.apply(step_in, dt_val)
        h = self.impulse_response(t)
        return np.cumsum(h) * (t[1] - t[0]) if len(t) > 1 else h.copy()

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        """Place FIR taps at the first M sample positions."""
        if self.smooth:
            if len(t) < 2:
                return np.zeros_like(t)
            dt_val = t[1] - t[0]
            imp_in = np.zeros_like(t)
            imp_in[0] = 1.0 / dt_val
            return self.apply(imp_in, dt_val)
        h = np.zeros_like(t, dtype=float)
        if len(t) > 0:
            dt_val = t[1] - t[0] if len(t) > 1 else 1.0
            n_place = min(len(self.taps), len(h))
            h[:n_place] = self.taps[:n_place] / dt_val  # continuous-time scaling
        return h

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """H(omega) = sum_k b[k] · exp(-j·omega·k·dt).

        Uses nominal dt=1 for normalised frequency.
        """
        H = np.polyval(self.taps[::-1], np.exp(-1j * omega))
        if self.smooth:
            H *= 1.0 / (1.0 + 1j * omega * self.smooth_tau)
        return H


    @classmethod
    def design_from_residual(
        cls,
        s_residual: np.ndarray,
        dt: float,
        n_taps: int = 72,
        ridge: float = 1e-6,
    ) -> "FIRDistortion":
        """Design FIR taps that correct a residual step response.

        Finds h such that h * s_residual approx u(t) in the least-squares
        sense over the first *n_taps* samples of the step window.
        Small ridge regularisation stabilises the solution for
        near-singular Toeplitz matrices.

        Parameters
        ----------
        s_residual : np.ndarray, shape (N,)
            Residual step response after partial correction.
        dt : float
            Sample spacing (ns).  Used only for consistency.
        n_taps : int
            Number of FIR taps.
        ridge : float
            Tikhonov regularisation.

        Returns
        -------
        FIRDistortion
            Optimal FIR correction filter.
        """
        s = np.asarray(s_residual, dtype=float)
        n = min(len(s), n_taps * 4)

        # Build Toeplitz convolution matrix S of shape (n, n_taps)
        S_mat = np.zeros((n, n_taps), dtype=float)
        for i in range(n):
            for k in range(n_taps):
                idx = i - k
                S_mat[i, k] = s[idx] if idx >= 0 else 0.0

        # Target: unit step
        target = np.ones(n, dtype=float)

        # Regularised least squares: (S^T S + lambda*I) h = S^T target
        A = S_mat.T @ S_mat + ridge * np.eye(n_taps)
        b = S_mat.T @ target
        h = np.linalg.solve(A, b)

        # Normalise: keep DC gain at 1
        dc = np.sum(h)
        if abs(dc) > 1e-15:
            h = h / dc

        return cls(taps=h)


# ---------------------------------------------------------------------------
# IIR filter distortion
# ---------------------------------------------------------------------------

@dataclass
class IIRDistortion(DistortionModel):
    """Infinite impulse response (IIR) distortion model.

    y[n] = sum_k b[k]·x[n-k] - sum_{k>0} a[k]·y[n-k],  a[0] = 1

    smooth=False (default): coefficients applied directly via lfilter.
    smooth=True: the IIR output is cascaded through a first-order lowpass
        with time constant ``smooth_tau`` to remove any sharp edges.

    Parameters
    ----------
    b : np.ndarray, shape (M,)
        Feed-forward coefficients.
    a : np.ndarray, shape (M,)
        Feed-back coefficients (a[0] must be 1).
    smooth : bool
        If True, cascade through a smoothing lowpass after the IIR.
        Default False.
    smooth_tau : float
        Time constant (ns) for the smoothing lowpass, used only when
        smooth=True.  Default 1.0.
    """

    b_coeffs: np.ndarray = field(default_factory=lambda: np.array([1.0]))
    a_coeffs: np.ndarray = field(default_factory=lambda: np.array([1.0]))
    smooth: bool = False
    smooth_tau: float = 1.0

    def __post_init__(self) -> None:
        self.b_coeffs = np.asarray(self.b_coeffs, dtype=float)
        self.a_coeffs = np.asarray(self.a_coeffs, dtype=float)
        if self.a_coeffs[0] != 1.0:
            self.b_coeffs = self.b_coeffs / self.a_coeffs[0]
            self.a_coeffs = self.a_coeffs / self.a_coeffs[0]

    @property
    def dt(self) -> float:
        """Nominal dt for IIR coefficients (for frequency_response)."""
        return 1.0

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        out = lfilter(self.b_coeffs, self.a_coeffs, waveform)
        if self.smooth:
            b_cont = [0.0, 1.0]
            a_cont = [self.smooth_tau, 1.0]
            b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
            out = lfilter(b, a, out)
        return out

    def step_response(self, t: np.ndarray) -> np.ndarray:
        if self.smooth:
            if len(t) < 2:
                return np.ones_like(t)
            dt_val = t[1] - t[0]
            return self.apply(np.ones(len(t)), dt_val)
        if len(t) < 2:
            return np.ones_like(t)
        n = len(t)
        step_in = np.ones(n)
        return lfilter(self.b_coeffs, self.a_coeffs, step_in)

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        if self.smooth:
            if len(t) < 2:
                return np.zeros_like(t)
            dt_val = t[1] - t[0]
            imp_in = np.zeros(len(t))
            imp_in[0] = 1.0 / dt_val
            return self.apply(imp_in, dt_val)
        if len(t) < 2:
            return np.zeros_like(t)
        dt_val = t[1] - t[0]
        n = len(t)
        imp_in = np.zeros(n)
        imp_in[0] = 1.0 / dt_val
        return lfilter(self.b_coeffs, self.a_coeffs, imp_in)

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """H(omega) evaluated via polynomial evaluation on unit circle."""
        H = np.zeros(len(omega), dtype=complex)
        z = np.exp(1j * omega)
        for i, bk in enumerate(self.b_coeffs):
            H += bk * z ** (-i)
        denom = np.zeros(len(omega), dtype=complex)
        for i, ak in enumerate(self.a_coeffs):
            denom += ak * z ** (-i)
        H = H / denom
        if self.smooth:
            H *= 1.0 / (1.0 + 1j * omega * self.smooth_tau)
        return H


# ---------------------------------------------------------------------------
# Custom frequency-domain transfer function
# ---------------------------------------------------------------------------

@dataclass
class CustomTransferDistortion(DistortionModel):
    """User-supplied complex frequency-domain transfer function H(omega).

    smooth=False (default): FFT -> multiply by H(omega) -> IFFT applied
        directly.  If H(omega) does not decay at high frequencies,
        input discontinuities are preserved.

    smooth=True: the frequency-domain output is cascaded through a
        first-order lowpass with time constant ``smooth_tau``.

    Parameters
    ----------
    omega_grid : np.ndarray, shape (F,)
        Frequency grid (rad/ns).
    H_grid : np.ndarray, shape (F,), complex
        H(omega) evaluated on omega_grid.
    smooth : bool
        If True, cascade through a smoothing lowpass after the transfer.
        Default False.
    smooth_tau : float
        Time constant (ns) for the smoothing lowpass, used only when
        smooth=True.  Default 1.0.
    """

    omega_grid: np.ndarray = field(
        default_factory=lambda: np.linspace(-np.pi, np.pi, 1024)
    )
    H_grid: np.ndarray = field(
        default_factory=lambda: np.ones(1024, dtype=complex)
    )
    smooth: bool = False
    smooth_tau: float = 1.0

    def __post_init__(self) -> None:
        self.omega_grid = np.asarray(self.omega_grid, dtype=float)
        self.H_grid = np.asarray(self.H_grid, dtype=complex)

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        """FFT -> multiply by H(omega) -> IFFT."""
        n = len(waveform)
        omega_wf = 2.0 * np.pi * np.fft.fftfreq(n, d=dt)
        H_interp = self._interpolate_H(omega_wf)
        X = np.fft.fft(waveform)
        Y = H_interp * X
        out = np.fft.ifft(Y).real
        if self.smooth:
            b_cont = [0.0, 1.0]
            a_cont = [self.smooth_tau, 1.0]
            b, a = bilinear(b_cont, a_cont, fs=1.0 / dt)
            out = lfilter(b, a, out)
        return out

    def _interpolate_H(self, omega_target: np.ndarray) -> np.ndarray:
        """Linear interpolation of H(omega) onto target frequency grid."""
        re = np.interp(
            omega_target, self.omega_grid, self.H_grid.real,
            left=self.H_grid.real[0], right=self.H_grid.real[-1],
        )
        im = np.interp(
            omega_target, self.omega_grid, self.H_grid.imag,
            left=self.H_grid.imag[0], right=self.H_grid.imag[-1],
        )
        return re + 1j * im

    def step_response(self, t: np.ndarray) -> np.ndarray:
        if len(t) < 2:
            return np.ones_like(t)
        dt_val = t[1] - t[0]
        step_in = np.where(t >= 0, 1.0, 0.0)
        return self.apply(step_in, dt_val)

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        if len(t) < 2:
            return np.zeros_like(t)
        dt_val = t[1] - t[0]
        imp_in = np.zeros_like(t)
        imp_in[0] = 1.0 / dt_val
        return self.apply(imp_in, dt_val)

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        H = self._interpolate_H(omega)
        if self.smooth:
            H *= 1.0 / (1.0 + 1j * omega * self.smooth_tau)
        return H


# ---------------------------------------------------------------------------
# Cascade distortion (Rol 2020 Sec.IV: series-cascaded IIR + FIR stages)
# ---------------------------------------------------------------------------

@dataclass
class CascadeDistortion(DistortionModel):
    """Series-cascaded distortion stages.

    Each stage is applied sequentially::

        output = stage_N( ... stage_2(stage_1(input)) ... )

    The composite frequency response is the product of individual
    stage responses.  An empty ``stages`` list is the identity system.

    Parameters
    ----------
    stages : list of DistortionModel
        Distortion stages applied in order.
    """

    stages: list = field(default_factory=list)

    def apply(self, waveform: np.ndarray, dt: float) -> np.ndarray:
        out = np.asarray(waveform, dtype=float)
        for s in self.stages:
            out = s.apply(out, dt)
        return out

    def step_response(self, t: np.ndarray) -> np.ndarray:
        if len(t) < 2:
            return np.ones_like(t)
        dt = float(t[1] - t[0])
        step_in = np.where(t >= 0, 1.0, 0.0)
        return self.apply(step_in, dt)

    def impulse_response(self, t: np.ndarray) -> np.ndarray:
        if len(t) < 2:
            return np.zeros_like(t)
        dt = float(t[1] - t[0])
        imp_in = np.zeros_like(t)
        imp_in[0] = 1.0 / dt
        return self.apply(imp_in, dt)

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        if not self.stages:
            return np.ones_like(omega, dtype=complex)
        H = self.stages[0].frequency_response(omega)
        for s in self.stages[1:]:
            H = H * s.frequency_response(omega)
        return H


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
    "CascadeDistortion",
]
