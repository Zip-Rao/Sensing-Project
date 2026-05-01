"""Multi-line transfer matrices for Z-control crosstalk."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.control.waveform import Waveform


def _uniform_dt(t_list: np.ndarray) -> float:
    t = np.asarray(t_list, dtype=float)
    if t.ndim != 1 or len(t) < 2:
        raise ValueError("TransferMatrix requires at least two time samples")
    diffs = np.diff(t)
    dt = float(np.median(diffs))
    if dt <= 0 or not np.allclose(diffs, dt, rtol=1e-5, atol=1e-12):
        raise ValueError("TransferMatrix requires a uniform time grid")
    return dt


@dataclass
class TransferMatrix:
    """Frequency-domain transfer matrix H_ji(omega).

    The convention is ``Phi_j(omega) = sum_i H_ji(omega) V_i(omega)``.
    Target names label on-chip flux channels, source names label AWG/control
    lines. Angular frequencies are in rad/ns.
    """

    elements: dict[tuple[str, str], np.ndarray] = field(default_factory=dict)
    frequency_axis: np.ndarray = field(default_factory=lambda: np.zeros(0))
    time_axis: np.ndarray | None = None

    def __post_init__(self) -> None:
        self.frequency_axis = np.asarray(self.frequency_axis, dtype=float)
        if self.frequency_axis.ndim != 1:
            raise ValueError("frequency_axis must be one-dimensional")
        if self.elements and len(self.frequency_axis) == 0:
            raise ValueError("frequency_axis is required when elements are set")

        order = np.argsort(self.frequency_axis)
        self.frequency_axis = self.frequency_axis[order]

        normalized: dict[tuple[str, str], np.ndarray] = {}
        for key, value in self.elements.items():
            arr = np.asarray(value, dtype=complex)
            if arr.ndim == 0:
                arr = np.full_like(self.frequency_axis, complex(arr), dtype=complex)
            if arr.ndim != 1:
                raise ValueError(f"transfer element {key} must be one-dimensional")
            if arr.shape != self.frequency_axis.shape:
                raise ValueError(
                    f"transfer element {key} has shape {arr.shape}, expected "
                    f"{self.frequency_axis.shape}"
                )
            normalized[key] = arr[order]
        self.elements = normalized

    @property
    def sources(self) -> list[str]:
        """Sorted source/control-line names."""
        return sorted({source for _, source in self.elements})

    @property
    def targets(self) -> list[str]:
        """Sorted on-chip target names."""
        return sorted({target for target, _ in self.elements})

    def H_ji(self, target: str, source: str) -> np.ndarray:
        """Return the transfer element from source to target."""
        return self.elements[(target, source)]

    def diagonal(self) -> dict[str, np.ndarray]:
        """Return self-response elements H_ii."""
        return {
            target: self.elements[(target, target)]
            for target in self.targets
            if (target, target) in self.elements
        }

    def off_diagonal(self) -> dict[tuple[str, str], np.ndarray]:
        """Return crosstalk elements H_ji where target != source."""
        return {
            (target, source): h
            for (target, source), h in self.elements.items()
            if target != source
        }

    def apply(self, source_voltages: dict[str, Waveform]) -> dict[str, FluxSignal]:
        """Apply the matrix to AWG source waveforms."""
        if not source_voltages:
            return {}

        t_ref: np.ndarray | None = None
        dt_ref: float | None = None
        n_ref: int | None = None
        spectra: dict[str, np.ndarray] = {}

        for source, waveform in source_voltages.items():
            t = np.asarray(waveform.t_list, dtype=float)
            dt = _uniform_dt(t)
            samples = np.asarray(waveform.samples, dtype=float)
            if t_ref is None:
                t_ref = t
                dt_ref = dt
                n_ref = len(samples)
            else:
                if len(samples) != n_ref or not np.allclose(t, t_ref):
                    raise ValueError("all source waveforms must share the same t_list")
                if not np.isclose(dt, dt_ref):
                    raise ValueError("all source waveforms must share the same dt")
            spectra[source] = np.fft.fft(samples)

        assert t_ref is not None and dt_ref is not None and n_ref is not None
        omega_fft = 2.0 * np.pi * np.fft.fftfreq(n_ref, d=dt_ref)
        result: dict[str, FluxSignal] = {}

        for target in self.targets:
            phi_omega = np.zeros(n_ref, dtype=complex)
            has_contribution = False
            for source, source_fft in spectra.items():
                if (target, source) not in self.elements:
                    continue
                h_ji = self.interpolate(target, source, omega_fft)
                phi_omega += h_ji * source_fft
                has_contribution = True
            if not has_contribution:
                continue
            phi = np.fft.ifft(phi_omega).real
            result[target] = FluxSignal(type=8, t_list=t_ref.copy(), signal=phi)

        return result

    def interpolate(
        self,
        target: str,
        source: str,
        omega_target: np.ndarray,
    ) -> np.ndarray:
        """Interpolate H_ji onto another angular-frequency grid."""
        h_arr = self.elements[(target, source)]
        omega_target = np.asarray(omega_target, dtype=float)
        real = np.interp(omega_target, self.frequency_axis, h_arr.real)
        imag = np.interp(omega_target, self.frequency_axis, h_arr.imag)
        return real + 1j * imag

    @classmethod
    def from_dc_matrix(
        cls,
        dc_matrix: np.ndarray,
        source_names: list[str],
        target_names: list[str],
        frequency_axis: np.ndarray | None = None,
    ) -> "TransferMatrix":
        """Build a frequency-flat transfer matrix from DC crosstalk values."""
        dc = np.asarray(dc_matrix, dtype=complex)
        if dc.shape != (len(target_names), len(source_names)):
            raise ValueError(
                "dc_matrix shape must be (len(target_names), len(source_names))"
            )
        if frequency_axis is None:
            frequency_axis = np.linspace(-np.pi, np.pi, 1024)
        omega = np.asarray(frequency_axis, dtype=float)
        elements: dict[tuple[str, str], np.ndarray] = {}
        for j, target in enumerate(target_names):
            for i, source in enumerate(source_names):
                elements[(target, source)] = np.full(
                    omega.shape,
                    dc[j, i],
                    dtype=complex,
                )
        return cls(elements=elements, frequency_axis=omega)


__all__ = ["TransferMatrix"]
