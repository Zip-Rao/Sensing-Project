"""Generic waveform data structures."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class Waveform:
    """Generic time-domain waveform without physical semantics."""

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
        """Waveform duration."""
        return float(self.t_list[-1] - self.t_list[0])

    @property
    def n_points(self) -> int:
        """Number of samples."""
        return int(len(self.t_list))

    def value_at(self, t: float) -> float:
        """Return the nearest sample, or zero outside the support."""
        if t < self.t_list[0] or t > self.t_list[-1]:
            return 0.0
        idx = int(np.argmin(np.abs(self.t_list - t)))
        return float(self.samples[idx])

    def truncate(self, t_start: float, t_end: float) -> "Waveform":
        """Return a copy with samples outside [t_start, t_end] zeroed."""
        mask = (self.t_list >= t_start) & (self.t_list <= t_end)
        return type(self)(
            t_list=self.t_list.copy(),
            samples=np.where(mask, self.samples, 0.0),
            metadata=dict(self.metadata),
        )

    def copy(self) -> "Waveform":
        """Return a deep copy of arrays and metadata."""
        return type(self)(
            t_list=self.t_list.copy(),
            samples=self.samples.copy(),
            metadata=dict(self.metadata),
        )


@dataclass
class CompositeWaveform(Waveform):
    """Concatenation of multiple waveforms."""

    components: list[Waveform] = field(default_factory=list)

    @classmethod
    def from_components(cls, components: list[Waveform]) -> "CompositeWaveform":
        """Build a composite waveform by concatenating components."""
        t_list = []
        samples = []
        offset = 0.0
        for waveform in components:
            t_list.extend(t + offset for t in waveform.t_list)
            samples.extend(waveform.samples)
            if len(waveform.t_list) > 0:
                offset = t_list[-1] + 1e-9
        return cls(
            t_list=np.asarray(t_list),
            samples=np.asarray(samples),
            components=list(components),
        )
