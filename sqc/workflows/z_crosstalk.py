"""Two-qubit Z-line crosstalk extraction and compensation workflow."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.control.waveform import Waveform
from sqc.devices.chip import ChipTopology
from sqc.hardware.transfer_matrix import TransferMatrix
from sqc.workflows.base import Workflow


@dataclass
class ZCrosstalkWorkflow(Workflow):
    """End-to-end Z-crosstalk demo for a two-qubit chip.

    A pulse is applied on qubit A's Z source. The workflow reconstructs the
    parasitic flux on qubit B, estimates H_BA, then applies a cancellation pulse
    on qubit B's own Z source to reduce the parasitic phase.
    """

    chip: ChipTopology
    flux_pulse_on_A: Waveform
    true_transfer_matrix: TransferMatrix
    qubit_A_name: str = "QA"
    qubit_B_name: str = "QB"
    regularization: float = 1e-9
    band_threshold: float = 1e-3

    def run(self) -> dict:
        source_voltages = {self.qubit_A_name: self.flux_pulse_on_A}
        on_chip_fluxes = self.true_transfer_matrix.apply(source_voltages)
        phi_A = on_chip_fluxes[self.qubit_A_name]
        phi_B_true = on_chip_fluxes[self.qubit_B_name]

        phi_B_reconstructed = self._reconstruct_phi_B(phi_B_true)
        omega = self._omega_grid(self.flux_pulse_on_A)
        v_a_omega = np.fft.fft(self.flux_pulse_on_A.samples)
        phi_b_omega = np.fft.fft(phi_B_reconstructed.samples)
        h_ba_estimated = self._estimate_transfer(phi_b_omega, v_a_omega)
        h_ba_true = self.true_transfer_matrix.interpolate(
            self.qubit_B_name,
            self.qubit_A_name,
            omega,
        )

        compensation_pulse = self._design_compensation_pulse(phi_b_omega, omega)
        compensated_fluxes = self.true_transfer_matrix.apply(
            {
                self.qubit_A_name: self.flux_pulse_on_A,
                self.qubit_B_name: compensation_pulse,
            }
        )
        phi_B_after = compensated_fluxes[self.qubit_B_name]

        phase_uncomp = self._parasitic_phase(phi_B_true)
        phase_comp = self._parasitic_phase(phi_B_after)
        compensation_factor = phase_uncomp / max(phase_comp, 1e-30)
        fit_error_db = self._fit_error_db(
            h_ba_estimated,
            h_ba_true,
            v_a_omega,
        )

        return {
            "omega": omega,
            "phi_A": phi_A,
            "phi_B_true": phi_B_true,
            "phi_B_reconstructed": phi_B_reconstructed,
            "H_BA_estimated": h_ba_estimated,
            "H_BA_true": h_ba_true,
            "fit_error_dB": fit_error_db,
            "compensation_pulse": compensation_pulse,
            "phi_B_after_compensation": phi_B_after,
            "parasitic_phase_uncompensated": phase_uncomp,
            "parasitic_phase_compensated": phase_comp,
            "compensation_factor": compensation_factor,
            "metrics": {
                "fit_error_dB": fit_error_db,
                "compensation_factor": compensation_factor,
                "phase_uncompensated": phase_uncomp,
                "phase_compensated": phase_comp,
            },
        }

    def _reconstruct_phi_B(self, phi_B_true: FluxSignal) -> FluxSignal:
        """Placeholder measurement reconstruction for the deterministic demo."""
        return FluxSignal(
            type=8,
            t_list=phi_B_true.t_list.copy(),
            signal=np.asarray(phi_B_true.samples, dtype=float).copy(),
        )

    def _estimate_transfer(
        self,
        output_omega: np.ndarray,
        input_omega: np.ndarray,
    ) -> np.ndarray:
        denom = np.abs(input_omega) ** 2 + self.regularization**2
        return output_omega * np.conj(input_omega) / denom

    def _design_compensation_pulse(
        self,
        phi_b_omega: np.ndarray,
        omega: np.ndarray,
    ) -> Waveform:
        if (self.qubit_B_name, self.qubit_B_name) not in self.true_transfer_matrix.elements:
            raise KeyError(
                "B-line self response H_BB is required for compensation"
            )
        h_bb = self.true_transfer_matrix.interpolate(
            self.qubit_B_name,
            self.qubit_B_name,
            omega,
        )
        v_b_comp_omega = -phi_b_omega * np.conj(h_bb) / (
            np.abs(h_bb) ** 2 + self.regularization**2
        )
        samples = np.fft.ifft(v_b_comp_omega).real
        return Waveform(
            t_list=self.flux_pulse_on_A.t_list.copy(),
            samples=samples,
            metadata={
                "purpose": "z_crosstalk_compensation",
                "target": self.qubit_B_name,
                "source": self.qubit_B_name,
            },
        )

    def _parasitic_phase(self, phi: Waveform) -> float:
        sensitivity = self._qubit_sensitivity(self.qubit_B_name)
        t = np.asarray(phi.t_list, dtype=float)
        samples = np.asarray(phi.samples, dtype=float)
        dt = float(np.median(np.diff(t)))
        phase_trace = np.cumsum(samples * sensitivity * dt)
        return float(np.max(np.abs(phase_trace)))

    def _qubit_sensitivity(self, qubit_name: str) -> float:
        qubit = self.chip.qubit(qubit_name)
        if hasattr(qubit, "frequency_sensitivity"):
            sensitivity = float(qubit.frequency_sensitivity(qubit.flux))
        elif hasattr(qubit, "sensitivity"):
            sensitivity = float(qubit.sensitivity())
        else:
            sensitivity = 1.0
        if abs(sensitivity) < 1e-12:
            return 1.0
        return sensitivity

    def _fit_error_db(
        self,
        h_estimated: np.ndarray,
        h_true: np.ndarray,
        input_omega: np.ndarray,
    ) -> float:
        weights = np.abs(input_omega)
        if np.max(weights) == 0:
            return 0.0
        mask = weights >= self.band_threshold * np.max(weights)
        if not np.any(mask):
            mask = np.ones_like(weights, dtype=bool)
        numerator = np.sqrt(np.mean(np.abs(h_estimated[mask] - h_true[mask]) ** 2))
        denominator = np.sqrt(np.mean(np.abs(h_true[mask]) ** 2))
        return float(20.0 * np.log10(numerator / max(denominator, 1e-30) + 1e-30))

    @staticmethod
    def _omega_grid(waveform: Waveform) -> np.ndarray:
        t = np.asarray(waveform.t_list, dtype=float)
        dt = float(np.median(np.diff(t)))
        return 2.0 * np.pi * np.fft.fftfreq(len(waveform.samples), d=dt)


__all__ = ["ZCrosstalkWorkflow"]
