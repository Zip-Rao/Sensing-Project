"""sqc.control.flux_signal — FluxSignal(Waveform) with type-based construction.

Backward-compatible with src.signal.Signal: supports type 0-8,
value_at, truncate (in-place), update_signal, params, copy, plot.

Physical interpretation: samples are in units of Phi_0.

Alias: Signal = FluxSignal (for backward compat in mirror layer).
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from sqc.config import CONFIG
from .waveform import Waveform, CompositeWaveform


# ---------------------------------------------------------------------------
# FluxSignal
# ---------------------------------------------------------------------------

class FluxSignal(Waveform):
    """Flux signal in units of Phi_0.

    Backward-compatible with src.signal.Signal: supports type-based
    construction (type 0-8) with keyword parameters.

    Parameters
    ----------
    type : int
        Signal type:
        0 = zero, 1 = constant, 2 = sinusoidal, 3 = Gaussian,
        4 = asymmetric impulse, 5 = double-peak, 6 = basis-expanded,
        7 = complex wave-packet, 8 = user-defined.
    t_list : array-like or None
        Time axis (ns).
    **kwargs
        Signal parameters (amplitude, frequency, phase, center, width,
        rise, fall, offset, noise_level, seed, etc.).
    """

    def __init__(
        self, type: int = 0, t_list=None, trigger: float = 0.0, **kwargs
    ) -> None:
        # Store type and parameters
        self._type: int = type
        self._params: dict = self._fill_default_params(type, kwargs)
        self.trigger: float = trigger

        # Build basis functions for type 6 BEFORE generating samples
        self._basis_functions: list = []
        if type == 6 and t_list is not None:
            self._build_basis_functions(t_list)

        # Generate samples
        t_arr = np.asarray(t_list, dtype=float) if t_list is not None else CONFIG.pulse.t_signal.copy()

        if type == 6 and self._basis_functions:
            # Generate with basis functions applied
            noise_level = self._params.get("noise_level", 0.0)
            seed = self._params.get("seed", 42)
            np.random.seed(seed)
            signal_arr = np.random.normal(0, noise_level, size=len(t_arr))
            offset = self._params.get("offset", 0.0)
            b = self._params.get("b", np.ones(self._params.get("n_basis", 10)))
            for i, phi_i in enumerate(self._basis_functions):
                if i < len(b):
                    signal_arr += b[i] * phi_i(t_arr)
            samples = signal_arr + offset
        else:
            samples = self._generate_signal(type, t_arr, self._params)

        # Initialize parent Waveform (dataclass __init__ + __post_init__)
        Waveform.__init__(
            self,  # type: ignore[arg-type]
            t_list=t_arr,
            samples=samples,
            metadata={"type": type, "params": dict(kwargs)},
        )

    # -- Legacy property accessors ------------------------------------------

    @property
    def signal(self) -> np.ndarray:
        """Backward-compat alias for samples."""
        return self.samples

    @signal.setter
    def signal(self, value: np.ndarray) -> None:
        self.samples = np.asarray(value, dtype=float)

    @property
    def type(self) -> int:
        return self._type

    @property
    def params(self) -> dict:
        return self._params

    @property
    def basis_functions(self) -> list:
        return self._basis_functions

    # -- Parameter defaults -------------------------------------------------

    @staticmethod
    def _fill_default_params(type: int, params: dict) -> dict:
        """Fill missing params with defaults matching src/signal.py:init()."""
        defaults = {
            "amplitude": 1.0,
            "frequency": 0.01,
            "phase": 0.0,
            "center": 0.0,
            "width": 30.0,
            "offset": 0.0,
            "t_rabi": CONFIG.pulse.make_time(0, 40),
            "tau": 20.0,
            "rise": 5.0,
            "fall": 2.0,
            "noise_level": 0.0,
            "seed": 42,
        }
        for key, value in defaults.items():
            if key not in params:
                params[key] = value
        if type == 6:
            if "n_basis" not in params:
                params["n_basis"] = 10
            if "b" not in params:
                params["b"] = np.ones(params["n_basis"])
            if "basis_type" not in params:
                params["basis_type"] = "bspline"
        return params

    # -- Signal generation --------------------------------------------------

    @staticmethod
    def _generate_signal(
        type: int, t_array: np.ndarray, params: dict
    ) -> np.ndarray:
        """Generate signal samples (without basis functions).

        Verbatim port of src/signal.py:Signal.generate().
        """
        noise_level = params.get("noise_level", 0.0)
        seed = params.get("seed", 42)
        np.random.seed(seed)
        signal = np.random.normal(0, noise_level, size=len(t_array))

        offset = params.get("offset", 0.0)
        amplitude = params.get("amplitude", 1.0)

        if type == 0:  # zero
            return signal + offset
        elif type == 1:  # constant
            signal += amplitude * np.ones_like(t_array)
            return signal + offset
        elif type == 2:  # sinusoidal
            omega = 2 * np.pi * params.get("frequency", 0.01)
            signal += amplitude * np.sin(
                omega * t_array + params.get("phase", 0.0)
            )
            return signal + offset
        elif type == 3:  # Gaussian
            sigma = params.get("width", 30.0) / 4.0
            signal += amplitude * np.exp(
                -((t_array - params.get("center", 0.0)) ** 2)
                / (2 * sigma**2)
            )
            return signal + offset
        elif type == 4:  # asymmetric impulse (double-exponential)
            rise = params.get("rise", 5.0)
            fall = params.get("fall", 2.0)
            t0 = params.get("center", 0.0)
            # Matches src/signal.py line 89:
            # amplitude * exp(-(t-t0)/fall - exp(-(t-t0)/rise))
            signal += amplitude * np.exp(
                -(t_array - t0) / fall
                - np.exp(-(t_array - t0) / rise)
            )
            return signal + offset
        elif type == 5:  # double-peak
            sigma = params.get("width", 30.0) / 4.0
            center = params.get("center", 0.0)
            w = params.get("width", 30.0)
            center1 = center - w / 2.0
            center2 = center + w / 2.0
            signal += amplitude * (
                np.exp(-((t_array - center1) ** 2) / (2 * sigma**2))
                + np.exp(-((t_array - center2) ** 2) / (2 * sigma**2))
            )
            return signal + offset
        elif type == 7:  # complex wave-packet
            def wave_packet(t, A, t0, w, f, phi):
                return (
                    0.0001
                    * A
                    * np.exp(-((t - t0) ** 2) / (2 * w**2))
                    * np.cos(2 * np.pi * f * (t - t0) + phi)
                )
            p1 = wave_packet(t_array, A=8, t0=50, w=12, f=0.06, phi=0)
            p2 = wave_packet(
                t_array, A=25, t0=90, w=18, f=0.03, phi=np.pi / 4
            )
            p3 = wave_packet(
                t_array, A=15, t0=120, w=20, f=0.02, phi=np.pi / 2
            )
            signal += p1 + p2 + p3
            return signal + offset
        elif type == 8:  # user-defined
            if "signal" not in params:
                raise ValueError("Custom signal not provided for type=8.")
            return np.asarray(params["signal"]) + offset
        else:
            raise ValueError(f"Unknown signal type: {type}")

    # -- basis functions ----------------------------------------------------

    def _build_basis_functions(self, t_list) -> None:
        """Build basis functions for type=6 (verbatim from src/signal.py)."""
        t_array = np.asarray(t_list, dtype=float)
        basis_type = self._params.get("basis_type", "bspline")
        n_basis = self._params.get("n_basis", 10)

        if basis_type == "bspline":
            from scipy.interpolate import BSpline

            degree = 3
            knots = np.linspace(
                t_array[0], t_array[-1], n_basis - degree + 1
            )
            knots = np.r_[[t_array[0]] * degree, knots, [t_array[-1]] * degree]
            self._basis_functions = []
            for i in range(n_basis):
                coeffs = np.zeros(n_basis)
                coeffs[i] = 1.0
                spline = BSpline(knots, coeffs, degree)
                t_lo, t_hi = t_array[0], t_array[-1]

                def make_safe_spline(spl, lo, hi):
                    def f(t):
                        t = np.asarray(t, dtype=float)
                        result = spl(t)
                        result[(t < lo) | (t > hi)] = 0.0
                        return result
                    return f

                self._basis_functions.append(
                    make_safe_spline(spline, t_lo, t_hi)
                )
        elif basis_type == "fourier":
            T_val = t_array[-1] - t_array[0]
            self._basis_functions = []
            self._basis_functions.append(
                lambda t, T_val=T_val: np.ones_like(t) / np.sqrt(max(T_val, 1e-12))
            )
            n_max = (n_basis - 1) // 2
            for n in range(1, n_max + 1):
                self._basis_functions.append(
                    lambda t, n=n, T_val=T_val, t0=t_array[0]: np.sqrt(2 / max(T_val, 1e-12))
                    * np.sin(2 * np.pi * n * (t - t0) / T_val)
                )
                self._basis_functions.append(
                    lambda t, n=n, T_val=T_val, t0=t_array[0]: np.sqrt(2 / max(T_val, 1e-12))
                    * np.cos(2 * np.pi * n * (t - t0) / T_val)
                )
        elif basis_type == "legendre":
            from scipy.special import legendre

            self._basis_functions = []
            for n_val in range(n_basis):
                Pn = legendre(n_val)
                t0 = t_array[0]
                T_val = t_array[-1] - t_array[0]

                def make_legendre(P, t0=t0, T_val=T_val):
                    def f(t):
                        return P(
                            2 * (t - t0) / max(T_val, 1e-12) - 1
                        )
                    return f

                self._basis_functions.append(make_legendre(Pn))
        else:
            raise ValueError(f"Unsupported basis type: {basis_type}")

    # -- Legacy methods -----------------------------------------------------

    def update_signal(self, **kwargs) -> None:
        """Update signal parameters and regenerate samples (in-place).

        Parameters
        ----------
        **kwargs
            Parameters to update (e.g., amplitude=2.0).
        """
        for key, value in kwargs.items():
            self._params[key] = value
        if self._type == 6 and self._basis_functions:
            t_array = self.t_list
            noise_level = self._params.get("noise_level", 0.0)
            seed = self._params.get("seed", 42)
            np.random.seed(seed)
            signal_arr = np.random.normal(0, noise_level, size=len(t_array))
            offset = self._params.get("offset", 0.0)
            b = self._params.get("b", np.ones(self._params.get("n_basis", 10)))
            for i, phi_i in enumerate(self._basis_functions):
                if i < len(b):
                    signal_arr += b[i] * phi_i(t_array)
            self.samples = signal_arr + offset
        else:
            self.samples = self._generate_signal(
                self._type, self.t_list, self._params
            )

    def truncate(self, t_start: float, t_end: float) -> None:
        """In-place truncation (legacy behavior).

        Overrides Waveform.truncate which returns a new object.
        Sets samples outside [t_start, t_end] to zero.

        Parameters
        ----------
        t_start : float
            Start time (ns).
        t_end : float
            End time (ns).
        """
        t_array = np.asarray(self.t_list)
        mask = (t_array < t_start) | (t_array > t_end)
        self.samples[mask] = 0.0

    def value_at(self, t: float) -> float:
        """Get signal value at time t. Returns 0 if out of range.

        Parameters
        ----------
        t : float
            Query time (ns).

        Returns
        -------
        float
        """
        signal_arr = self.samples
        t_array = np.asarray(self.t_list)
        idx = int(np.abs(t_array - t).argmin())
        if t_array[0] <= t <= t_array[-1]:
            return float(signal_arr[idx])
        return 0.0

    def samples_on(self, t_global: np.ndarray) -> np.ndarray:
        """Project signal samples onto a global time axis.

        Computes local time t_loc = t_global - self.trigger, then
        linearly interpolates self.samples onto t_global within the
        signal's local time window [0, t_list[-1]].  Outside that
        window the contribution is zero.

        Parameters
        ----------
        t_global : np.ndarray
            Global time axis (ns).

        Returns
        -------
        np.ndarray
            Signal values at each global time point, shape (len(t_global),).
        """
        out = np.zeros(len(t_global), dtype=float)
        t_loc = t_global - self.trigger
        mask = (t_loc >= self.t_list[0]) & (t_loc <= self.t_list[-1])
        out[mask] = np.interp(t_loc[mask], self.t_list, self.samples)
        return out

    def copy(self) -> "FluxSignal":
        """Create a deep copy of this signal.

        Returns
        -------
        FluxSignal
        """
        return FluxSignal(
            type=self._type, t_list=self.t_list.copy(), **self._params
        )

    def plot(self, ax=None, **kwargs):
        """Plot the signal waveform.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
        **kwargs
            Passed to ax.plot().

        Returns
        -------
        matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(10, 4))
        ax.plot(self.t_list, self.samples, **kwargs)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Signal Amplitude")
        ax.set_title("Signal Waveform")
        ax.grid(True)
        return ax


# ---------------------------------------------------------------------------
# Backward-compat alias
# ---------------------------------------------------------------------------

Signal = FluxSignal


# ---------------------------------------------------------------------------
# CompositeSignal
# ---------------------------------------------------------------------------

class CompositeSignal(CompositeWaveform):
    """Backward-compat for src.signal.CompositeSignal.

    Parameters
    ----------
    signals : list[FluxSignal]
        List of FluxSignal objects to concatenate in time order.
    """

    def __init__(self, signals: list[FluxSignal]) -> None:
        self.signals = list(signals)
        t_list = self._compute_t_list()
        samples = self._compute_samples()
        CompositeWaveform.__init__(
            self,  # type: ignore[arg-type]
            t_list=np.array(t_list, dtype=float),
            samples=np.array(samples, dtype=float),
            components=list(signals),
        )

    @property
    def signal(self) -> np.ndarray:
        """Backward-compat alias for samples."""
        return self.samples

    def _compute_t_list(self) -> np.ndarray:
        """Compute merged time list from sub-signals."""
        t_list = []
        curr = 0.0
        for sig in self.signals:
            pulse_list = [t + curr for t in sig.t_list]
            t_list.extend(pulse_list)
            if pulse_list:
                curr = pulse_list[-1]  # P7.5: removed 1e-9 separator
        return np.array(t_list, dtype=float)

    def _compute_samples(self) -> np.ndarray:
        """Concatenate sub-signal samples."""
        return np.concatenate([sig.signal for sig in self.signals])

    def plot(self, ax=None, **kwargs):
        """Plot the composite signal."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(10, 4))
        ax.plot(self.t_list, self.samples, **kwargs)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Composite Signal Amplitude")
        ax.set_title("Composite Signal Waveform")
        ax.grid(True)
        return ax
