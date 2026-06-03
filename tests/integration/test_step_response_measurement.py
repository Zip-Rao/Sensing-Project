"""Integration tests: protocol-driven step response measurement (P9.B)."""

import numpy as np
import pytest

from sqc.hardware.distortion import SingleExponentialDistortion
from sqc.hardware.control_line import ControlLine


def _rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


class TestStepResponseMeasurementPlumbing:
    """Fast tests: StepResponseMeasurement construction and integration."""

    @staticmethod
    def _make_line():
        return ControlLine(
            name="Z0", kind="z",
            source="AWG0", target="Q0",
            transfer_function=SingleExponentialDistortion(
                amplitude=0.05, tau=50.0,
            ),
        )

    def test_construction_cryoscope(self, qubit_default):
        """StepResponseMeasurement constructs with cryoscope protocol."""
        from sqc.calibration.waveform import StepResponseMeasurement

        meas = StepResponseMeasurement(
            qubit=qubit_default,
            control_line=self._make_line(),
            protocol="cryoscope",
            step_amplitude=0.05,
            t_max=200.0,
        )
        assert meas.protocol == "cryoscope"
        assert meas.step_amplitude == 0.05

    def test_construction_delay_ramsey(self, qubit_default):
        """StepResponseMeasurement constructs with delay_ramsey protocol."""
        from sqc.calibration.waveform import StepResponseMeasurement

        meas = StepResponseMeasurement(
            qubit=qubit_default,
            control_line=None,
            protocol="delay_ramsey",
            step_amplitude=0.03,
            t_max=300.0,
        )
        assert meas.protocol == "delay_ramsey"

    def test_unknown_protocol_raises(self, qubit_default):
        """Unknown protocol name raises ValueError."""
        from sqc.calibration.waveform import StepResponseMeasurement

        with pytest.raises(ValueError, match="Unknown protocol"):
            StepResponseMeasurement(
                qubit=qubit_default, protocol="nonexistent",
            )

    def test_flux_bias_warning_cryoscope(self, qubit_default):
        """Warns if qubit.flux is not at sweet spot for cryoscope."""
        import warnings

        from sqc.calibration.waveform import StepResponseMeasurement

        qubit_default.flux = 0.15  # off sweet spot
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            StepResponseMeasurement(
                qubit=qubit_default, protocol="cryoscope",
            )
            cryo_warnings = [
                x for x in w
                if "sweet spot" in str(x.message).lower()
            ]
            assert len(cryo_warnings) > 0

    def test_waveform_calibration_with_measurement(self, qubit_default):
        """WaveformCalibration delegates to StepResponseMeasurement when set."""
        from sqc.calibration.waveform import StepResponseMeasurement
        from sqc.calibration.waveform import WaveformCalibration

        meas = StepResponseMeasurement(
            qubit=qubit_default,
            control_line=self._make_line(),
            protocol="cryoscope",
            t_max=200.0,
        )
        cal = WaveformCalibration(
            measurement=meas,
            method="transfer_function",
            fit_type="single_exp",
            n_exp_components=1,
        )
        result = cal.calibrate()
        assert result.name == "transfer_function"
        assert len(result.inputs) > 0

    def test_predistortion_workflow_protocol_field(self):
        """PredistortionValidationWorkflow accepts measurement_protocol."""
        from sqc.workflows.predistortion_validation import (
            PredistortionValidationWorkflow,
        )
        from sqc.control.waveform import Waveform

        wf = Waveform(
            t_list=np.arange(0, 100, 0.5),
            samples=np.ones(200),
        )
        dist = SingleExponentialDistortion(amplitude=0.05, tau=50.0)

        workflow = PredistortionValidationWorkflow(
            target_waveform=wf,
            true_distortion=dist,
            measurement_protocol=None,  # default analytical path
        )
        result = workflow.run()
        assert "metrics" in result
        assert result["metrics"]["improvement_factor"] > 10
