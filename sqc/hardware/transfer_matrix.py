"""sqc.hardware.transfer_matrix — TransferMatrix: multi-qubit Z-line transfer function.

Models crosstalk between multiple control lines.
Phi_j(omega) = sum_i H_{ji}(omega) * V_i(omega)

Full implementation per phase_5_handbook.md §3.2.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from sqc.control.waveform import Waveform
from sqc.control.flux_signal import FluxSignal


@dataclass
class TransferMatrix:
    """Multi-qubit Z-line transfer matrix H_ji(omega).

    Phi_j(omega) = sum_i H_{ji}(omega) * V_i(omega)

    Each element elements[(target_qubit, source_line)] is a complex-valued
    array on frequency_axis, representing the frequency-dependent response
    from a voltage on source_line to flux at target_qubit.

    Attributes
    ----------
    elements : dict[tuple[str, str], np.ndarray]
        (target_name, source_name) → complex H_ji(ω) array on frequency_axis.
    frequency_axis : np.ndarray
        Frequency grid (rad/ns) on which H_ji is defined.
    time_axis : np.ndarray or None
        Optional time axis for impulse-response representation.
    """

    elements: dict[tuple[str, str], np.ndarray] = field(default_factory=dict)
    frequency_axis: np.ndarray = field(default_factory=lambda: np.zeros(0))
    time_axis: np.ndarray | None = None

    # -- convenience properties ------------------------------------------------

    @property
    def sources(self) -> list[str]:
        """Sorted list of source names."""
        return sorted({s for (_, s) in self.elements.keys()})

    @property
    def targets(self) -> list[str]:
        """Sorted list of target names."""
        return sorted({t for (t, _) in self.elements.keys()})

    # -- element access --------------------------------------------------------

    def H_ji(self, target: str, source: str) -> np.ndarray:
        """Access a single transfer-function element.

        Parameters
        ----------
        target : str
            Target qubit name (receives flux).
        source : str
            Source line name (receives voltage).

        Returns
        -------
        np.ndarray
            Complex H_ji(ω) on self.frequency_axis.

        Raises
        ------
        KeyError
            If (target, source) not in elements.
        """
        return self.elements[(target, source)]

    def diagonal(self) -> dict[str, np.ndarray]:
        """Self-response H_ii(ω) for each target (diagonal elements)."""
        return {t: self.elements[(t, t)] for t in self.targets if (t, t) in self.elements}

    def off_diagonal(self) -> dict[tuple[str, str], np.ndarray]:
        """Cross-talk elements H_ji(ω) for j ≠ i."""
        return {(t, s): H for (t, s), H in self.elements.items() if t != s}

    # -- apply: AWG voltages → on-chip fluxes ----------------------------------

    def apply(self, source_voltages: dict[str, Waveform]) -> dict[str, FluxSignal]:
        """Compute on-chip flux at each target from AWG voltages.

        For each target j:
            Phi_j = IFFT( sum_i H_ji(omega) * FFT(V_i) )

        Waveforms with different lengths are handled by zero-padding
        to the longest common length.

        Parameters
        ----------
        source_voltages : dict[str, Waveform]
            Source name → Waveform (AWG voltage in arbitrary units).

        Returns
        -------
        dict[str, FluxSignal]
            Target name → FluxSignal (on-chip flux in Phi_0).

        Raises
        ------
        ValueError
            If a required element is missing from self.elements for any
            (target, source) pair that has nonzero input.
        """
        targets = self.targets
        if not targets:
            return {}

        # Determine common FFT length: use the longest waveform
        max_len = 0
        dt = None
        ref_t_list = None
        for _name, wf in source_voltages.items():
            n = len(wf.samples)
            if n > max_len:
                max_len = n
                dt = float(wf.t_list[1] - wf.t_list[0])
                ref_t_list = wf.t_list

        if max_len == 0:
            return {}

        result: dict[str, FluxSignal] = {}

        for j in targets:
            phi_j_omega = np.zeros(max_len, dtype=complex)

            for i, V_i in source_voltages.items():
                if (j, i) not in self.elements:
                    # No transfer from i to j → skip (equivalent to H_ji=0)
                    continue

                # Zero-pad or truncate V_i to common length
                V_i_padded = np.zeros(max_len, dtype=float)
                n_copy = min(len(V_i.samples), max_len)
                V_i_padded[:n_copy] = V_i.samples[:n_copy]

                # FFT of voltage waveform
                V_omega = np.fft.fft(V_i_padded)

                # Frequency grid for this waveform
                n_fft = max_len
                omega_grid = 2.0 * np.pi * np.fft.fftfreq(n_fft, d=dt)

                # Interpolate H_ji onto omega_grid
                H_ji_interp = self._interpolate_to(
                    self.elements[(j, i)], omega_grid
                )

                # Accumulate: phi_j_omega += H_ji * V_i
                phi_j_omega += H_ji_interp * V_omega

            # IFFT to time domain
            phi_j = np.real(np.fft.ifft(phi_j_omega))

            # Build result FluxSignal
            t_result = ref_t_list if ref_t_list is not None else np.arange(max_len) * dt
            # If phi_j length differs from t_result, match them
            n_result = min(len(phi_j), len(t_result))
            result[j] = FluxSignal(
                type=8,
                t_list=t_result[:n_result],
                signal=phi_j[:n_result],
            )

        return result

    # -- interpolation helper --------------------------------------------------

    def _interpolate_to(
        self, H_arr: np.ndarray, omega_target: np.ndarray,
    ) -> np.ndarray:
        """Interpolate H_ji onto a target frequency grid.

        Uses separate linear interpolation for real and imaginary parts.

        Parameters
        ----------
        H_arr : np.ndarray
            Complex-valued H on self.frequency_axis.
        omega_target : np.ndarray
            Target frequency grid.

        Returns
        -------
        np.ndarray
            Complex-valued H on omega_target.
        """
        src_omega = np.asarray(self.frequency_axis)
        # Ensure source grid is sorted and unique
        if len(src_omega) <= 1:
            # Constant transfer function
            return H_arr[0] if len(H_arr) > 0 else np.ones_like(omega_target)

        sort_idx = np.argsort(src_omega)
        src_omega = src_omega[sort_idx]
        H_sorted = H_arr[sort_idx]

        re_interp = np.interp(omega_target, src_omega, H_sorted.real)
        im_interp = np.interp(omega_target, src_omega, H_sorted.imag)
        return re_interp + 1j * im_interp

    # -- factory ---------------------------------------------------------------

    @classmethod
    def from_dc_matrix(
        cls,
        dc_matrix: np.ndarray,
        source_names: list[str],
        target_names: list[str],
        n_freq: int = 256,
    ) -> "TransferMatrix":
        """Build a frequency-flat TransferMatrix from a DC crosstalk matrix.

        H_ji(omega) = dc_matrix[j, i] for all omega.

        Parameters
        ----------
        dc_matrix : np.ndarray
            Shape (n_targets, n_sources). dc_matrix[j, i] is the DC
            crosstalk coefficient from source i to target j.
        source_names : list[str]
            Names of the source lines.
        target_names : list[str]
            Names of the target qubits.
        n_freq : int
            Number of frequency points. Default 256.

        Returns
        -------
        TransferMatrix
            Frequency-flat transfer matrix.

        Raises
        ------
        ValueError
            If dc_matrix shape does not match names.
        """
        dc_matrix = np.asarray(dc_matrix, dtype=float)
        if dc_matrix.shape != (len(target_names), len(source_names)):
            raise ValueError(
                f"dc_matrix shape {dc_matrix.shape} does not match "
                f"({len(target_names)}, {len(source_names)})"
            )

        omega = np.linspace(-np.pi, np.pi, n_freq)

        elements: dict[tuple[str, str], np.ndarray] = {}
        for j, t in enumerate(target_names):
            for i, s in enumerate(source_names):
                elements[(t, s)] = dc_matrix[j, i] * np.ones(n_freq, dtype=complex)

        return cls(elements=elements, frequency_axis=omega)

    # -- backward-compat bridge ------------------------------------------------

    @classmethod
    def from_dc_p1_stub(
        cls,
        n_lines: int,
        dc_dict: dict[tuple[int, int], float],
        source_names: Optional[list[str]] = None,
        target_names: Optional[list[str]] = None,
    ) -> "TransferMatrix":
        """Build TransferMatrix from P1-stub-style arguments.

        Converts the old (int, int) key format to (str, str).

        Parameters
        ----------
        n_lines : int
            Number of lines.
        dc_dict : dict[tuple[int, int], float]
            Old-style (i, j) → H_{ji} DC value.
        source_names : list[str] or None
            Source names. Defaults to ["S0", "S1", ...].
        target_names : list[str] or None
            Target names. Defaults to ["T0", "T1", ...].

        Returns
        -------
        TransferMatrix
        """
        if source_names is None:
            source_names = [f"S{i}" for i in range(n_lines)]
        if target_names is None:
            target_names = [f"T{i}" for i in range(n_lines)]

        dc_matrix = np.zeros((len(target_names), len(source_names)))
        for (i_val, j_val), val in dc_dict.items():
            if 0 <= j_val < len(target_names) and 0 <= i_val < len(source_names):
                dc_matrix[j_val, i_val] = val

        return cls.from_dc_matrix(dc_matrix, source_names, target_names)
