"""End-to-end predistortion validation workflow."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.calibration.predistortion import PredistortionDesigner
from sqc.calibration.transfer_function import TransferFunctionCalibration
from sqc.control.waveform import Waveform
from sqc.hardware.control_line import ControlLine
from sqc.hardware.distortion import DistortionModel
from sqc.workflows.base import Workflow


@dataclass
class PredistortionValidationWorkflow(Workflow):
    """Simulate a measured line, design predistortion, and score improvement."""

    qubit: object | None
    target_waveform: Waveform
    true_distortion: DistortionModel
    designer: PredistortionDesigner = field(
        default_factory=lambda: PredistortionDesigner(method="iir_inverse")
    )
    fit_type: str = "iir"
    tolerance: float = 1e-3

    def run(self) -> dict:
        line = ControlLine(
            name="Z0",
            kind="z",
            source="AWG0:CH1",
            target=self._target_name(),
            transfer_function=self.true_distortion,
        )

        on_chip_uncorrected = line.apply(self.target_waveform)
        calibration = TransferFunctionCalibration(
            qubit=self.qubit,
            method="simulated",
            fit_type=self.fit_type,  # type: ignore[arg-type]
            measured_step_response=None,
            transfer_model=self.true_distortion,
            t_list=self.target_waveform.t_list,
        )
        table = calibration.calibrate()
        measured_model = TransferFunctionCalibration.model_from_table(table)

        awg_predistorted = self.designer.predistort(
            self.target_waveform, measured_model
        )
        on_chip_corrected = line.apply(awg_predistorted)

        rmse_uncorrected = self._rmse(on_chip_uncorrected, self.target_waveform)
        rmse_corrected = self._rmse(on_chip_corrected, self.target_waveform)
        metrics = {
            "rmse_uncorrected": rmse_uncorrected,
            "rmse_corrected": rmse_corrected,
            "improvement_factor": rmse_uncorrected / max(rmse_corrected, 1e-30),
            "max_abs_uncorrected": self._max_abs_error(
                on_chip_uncorrected, self.target_waveform
            ),
            "max_abs_corrected": self._max_abs_error(
                on_chip_corrected, self.target_waveform
            ),
            "settling_uncorrected_ns": self._settling_time(
                on_chip_uncorrected, self.target_waveform, self.tolerance
            ),
            "settling_corrected_ns": self._settling_time(
                on_chip_corrected, self.target_waveform, self.tolerance
            ),
        }

        return {
            "target": self.target_waveform,
            "on_chip_uncorrected": on_chip_uncorrected,
            "awg_predistorted": awg_predistorted,
            "on_chip_corrected": on_chip_corrected,
            "transfer_table": table,
            "measured_model": measured_model,
            "metrics": metrics,
        }

    def _target_name(self) -> str:
        if self.qubit is None:
            return "Q0"
        if hasattr(self.qubit, "spec"):
            return self.qubit.spec().name
        if hasattr(self.qubit, "name"):
            return str(self.qubit.name)
        return "Q0"

    @staticmethod
    def _rmse(measured: Waveform, target: Waveform) -> float:
        return float(
            np.sqrt(np.mean((np.asarray(measured.samples) - target.samples) ** 2))
        )

    @staticmethod
    def _max_abs_error(measured: Waveform, target: Waveform) -> float:
        return float(np.max(np.abs(np.asarray(measured.samples) - target.samples)))

    @staticmethod
    def _settling_time(
        measured: Waveform,
        target: Waveform,
        tolerance: float = 1e-3,
    ) -> float:
        error = np.abs(np.asarray(measured.samples) - target.samples)
        scale = max(float(np.max(np.abs(target.samples))), 1.0)
        within = error <= tolerance * scale
        if np.all(within):
            return 0.0
        for idx in range(len(within)):
            if np.all(within[idx:]):
                return float(target.t_list[idx] - target.t_list[0])
        return float(target.t_list[-1] - target.t_list[0])


__all__ = ["PredistortionValidationWorkflow"]
