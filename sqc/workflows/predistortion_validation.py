"""sqc.workflows.predistortion_validation — PredistortionValidationWorkflow.

End-to-end predistortion validation workflow:

1. Inject known distortion into a ControlLine.
2. Measure step response via TransferFunctionCalibration.
3. Fit H(omega) model.
4. Design inverse via PredistortionDesigner.
5. Apply predistortion and verify improvement.

Returns a metrics dict with rmse_uncorrected, rmse_corrected,
improvement_factor, and settling times.

Per phase_4_handbook.md §3.5.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from sqc.control.waveform import Waveform
from sqc.workflows.base import Workflow


@dataclass
class PredistortionValidationWorkflow(Workflow):
    """End-to-end predistortion validation.

    1. Build control line with true distortion.
    2. Forward-model the target through the line (uncorrected).
    3. Calibrate transfer function from the known distortion.
    4. Design inverse filter via PredistortionDesigner.
    5. Apply predistortion.
    6. Re-measure with predistorted AWG waveform.
    7. Compute residuals and improvement metrics.

    Parameters
    ----------
    target_waveform : Waveform
        Desired on-chip waveform.
    true_distortion : DistortionModel
        Ground-truth distortion model (in simulation, we know this).
    designer : PredistortionDesigner or None
        Inverse filter designer. Uses default if None.
    control_line_params : dict
        Extra parameters for ControlLine construction (name, kind, etc.).
    """

    target_waveform: Waveform
    true_distortion: object  # DistortionModel

    # Optional overrides
    designer: Optional[object] = None  # PredistortionDesigner
    control_line_params: dict = field(default_factory=dict)

    # -- P9.B: protocol-driven measurement --
    qubit: object | None = None
    measurement_protocol: str | None = None  # None => analytical path

    def run(self) -> dict:
        """Execute the predistortion validation workflow.

        Returns
        -------
        dict
            Keys: target, on_chip_uncorrected, on_chip_corrected,
            awg_predistorted, metrics.
        """
        from sqc.calibration.waveform import PredistortionDesigner
        from sqc.calibration.waveform import WaveformCalibration

        # Resolve designer
        designer = self.designer
        if designer is None:
            designer = PredistortionDesigner(method="auto")

        # 1. Build control line
        line = self._build_control_line()

        # 2. Without predistortion: AWG -> on-chip is distorted
        on_chip_uncorrected = line.apply(self.target_waveform)

        # 3. Calibrate transfer function
        #    Protocol-driven path:  measurement_protocol → quantum simulation.
        #    Analytical path:       distortion.step_response() directly.
        fit_type = self._infer_fit_type(self.true_distortion)
        if self.measurement_protocol is not None:
            cal = WaveformCalibration(
                qubit=self.qubit,
                control_line=line,
                measurement_protocol=self.measurement_protocol,
                method="transfer_function",
                fit_type=fit_type,
                n_exp_components=3,
                t_max=float(self.target_waveform.t_list[-1]) * 2,
            )
        else:
            cal = WaveformCalibration(
                distortion=self.true_distortion,
                method="simulation",
                fit_type=fit_type,
                n_exp_components=3,
            )
        measured_table = cal.calibrate()
        measured_model = cal.to_distortion_model()

        # 4. Design inverse
        dt = float(self.target_waveform.t_list[1] - self.target_waveform.t_list[0])
        inverse_model = designer.design(measured_model, dt=dt)

        # 5. Predistort
        awg_predistorted = inverse_model.apply_to_waveform(self.target_waveform)

        # 6. Apply forward distortion to predistorted AWG waveform
        on_chip_corrected = line.apply(awg_predistorted)

        # 7. Compute metrics
        metrics = self._compute_metrics(
            on_chip_uncorrected, on_chip_corrected, self.target_waveform,
        )
        metrics["calibration_fit_type"] = measured_table.fit_params
        metrics["inverse_model_type"] = type(inverse_model).__name__

        return {
            "target": self.target_waveform,
            "on_chip_uncorrected": on_chip_uncorrected,
            "on_chip_corrected": on_chip_corrected,
            "awg_predistorted": awg_predistorted,
            "inverse_model": inverse_model,
            "measured_model": measured_model,
            "metrics": metrics,
        }

    def _build_control_line(self):
        """Build a ControlLine with the true distortion."""
        from sqc.hardware.control_line import ControlLine

        params = {
            "name": "Z0",
            "kind": "z",
            "source": "AWG0",
            "target": "Q0",
            "transfer_function": self.true_distortion,
        }
        params.update(self.control_line_params)
        return ControlLine(**params)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _infer_fit_type(distortion: object) -> str:
        """Infer the best fit_type from the distortion class name."""
        cls_name = type(distortion).__name__
        if "Single" in cls_name:
            return "single_exp"
        elif "Multi" in cls_name:
            return "multi_exp"
        elif "FIR" in cls_name:
            return "fir"
        elif "IIR" in cls_name:
            return "iir"
        return "multi_exp"

    # ------------------------------------------------------------------
    # Metrics
    # ------------------------------------------------------------------

    def _compute_metrics(
        self,
        uncorrected: Waveform,
        corrected: Waveform,
        target: Waveform,
    ) -> dict:
        """Compute improvement metrics."""
        rmse_uncorrected = self._rmse(uncorrected, target)
        rmse_corrected = self._rmse(corrected, target)
        return {
            "rmse_uncorrected": rmse_uncorrected,
            "rmse_corrected": rmse_corrected,
            "improvement_factor": rmse_uncorrected / max(rmse_corrected, 1e-30),
            "settling_uncorrected_ns": self._settling_time(uncorrected, target),
            "settling_corrected_ns": self._settling_time(corrected, target),
        }

    @staticmethod
    def _rmse(measured: Waveform, target: Waveform) -> float:
        """Root-mean-square error between measured and target."""
        return float(np.sqrt(np.mean((measured.samples - target.samples) ** 2)))

    @staticmethod
    def _settling_time(
        w: Waveform, target: Waveform, tolerance: float = 0.001,
    ) -> float:
        """Time when the waveform settles within tolerance of target.

        Returns the last time index where abs(error) > tolerance,
        or 0.0 if immediately settled.
        """
        error = np.abs(w.samples - target.samples)
        # Find the last index where error exceeds tolerance
        bad = np.where(error > tolerance)[0]
        if len(bad) == 0:
            return 0.0
        last_bad_idx = bad[-1]
        if last_bad_idx + 1 < len(w.t_list):
            return float(w.t_list[last_bad_idx])
        return float(w.t_list[-1])
