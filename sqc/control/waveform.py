"""sqc.control.waveform — Generic time-domain waveform data structures.

Waveform: atomic time-domain signal (mutable).
CompositeWaveform: concatenation of multiple Waveforms.

See _refactor_plan.md §5.2.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np


@dataclass
class Waveform:
    """Generic time-domain waveform, semantically neutral.

    For physical interpretation (e.g., flux signal), use FluxSignal.

    Attributes
    ----------
    t_list : np.ndarray
        Time points (ns).
    samples : np.ndarray
        Signal values at each time point.
    metadata : dict
        Arbitrary key-value metadata.
    """

    t_list: np.ndarray
    samples: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.t_list = np.asarray(self.t_list, dtype=float)
        self.samples = np.asarray(self.samples, dtype=float)
        if self.t_list.shape != self.samples.shape:
            raise ValueError(
                f"shape mismatch: t_list {self.t_list.shape} vs "
                f"samples {self.samples.shape}"
            )

    @property
    def duration(self) -> float:
        """Total duration (ns)."""
        return float(self.t_list[-1] - self.t_list[0])

    @property
    def n_points(self) -> int:
        """Number of time points."""
        return len(self.t_list)

    def value_at(self, t: float) -> float:
        """Sample-and-hold at time t. Returns 0 if t out of range.

        Parameters
        ----------
        t : float
            Query time (ns).

        Returns
        -------
        float
            Signal value at nearest time point, or 0 if out of range.
        """
        if t < self.t_list[0] or t > self.t_list[-1]:
            return 0.0
        idx = int(np.argmin(np.abs(self.t_list - t)))
        return float(self.samples[idx])

    def samples_on(self, t_global: np.ndarray) -> np.ndarray:
        """Project samples onto a global time axis.

        Linearly interpolates self.samples onto *t_global* within
        ``[self.t_list[0], self.t_list[-1]]``; zero outside.

        Parameters
        ----------
        t_global : np.ndarray
            Global time points (ns).

        Returns
        -------
        np.ndarray
            Interpolated values, shape ``(len(t_global),)``.
        """
        out = np.zeros(len(t_global), dtype=float)
        mask = (t_global >= self.t_list[0]) & (t_global <= self.t_list[-1])
        out[mask] = np.interp(t_global[mask], self.t_list, self.samples)
        return out

    def truncate(self, t_start: float, t_end: float) -> "Waveform":
        """Return a NEW waveform with samples zeroed outside [t_start, t_end].

        Parameters
        ----------
        t_start : float
            Start time (ns).
        t_end : float
            End time (ns).

        Returns
        -------
        Waveform
            New Waveform with zeroed edges.
        """
        mask = (self.t_list >= t_start) & (self.t_list <= t_end)
        new_samples = np.where(mask, self.samples, 0.0)
        return type(self)(
            t_list=self.t_list.copy(),
            samples=new_samples,
            metadata=dict(self.metadata),
        )

    def copy(self) -> "Waveform":
        """Return a deep copy of this waveform."""
        return type(self)(
            t_list=self.t_list.copy(),
            samples=self.samples.copy(),
            metadata=dict(self.metadata),
        )

    def plot(self, ax=None, **kwargs):
        """Plot the waveform.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. Creates new figure if None.
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
        ax.set_ylabel("Amplitude")
        ax.grid(True)
        return ax


@dataclass
class CompositeWaveform(Waveform):
    """Concatenation of multiple Waveforms in time.

    Attributes
    ----------
    components : list[Waveform]
        Constituent waveforms concatenated in time order.
    """

    components: list[Waveform] = field(default_factory=list)

    @classmethod
    def from_components(cls, components: list[Waveform]) -> "CompositeWaveform":
        """Build a CompositeWaveform from a list of Waveform objects.

        Parameters
        ----------
        components : list[Waveform]
            Waveforms to concatenate in order.

        Returns
        -------
        CompositeWaveform
            New waveform with concatenated time and sample arrays.
        """
        t_list = []
        samples = []
        offset = 0.0
        for w in components:
            if len(w.t_list) == 0:
                continue
            t_list.extend(float(t + offset) for t in w.t_list)
            samples.extend(w.samples)
            offset = t_list[-1] + 1e-9  # tiny gap to avoid duplicate time points
        return cls(
            t_list=np.array(t_list, dtype=float),
            samples=np.array(samples, dtype=float),
            components=list(components),
        )
