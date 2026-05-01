"""Integration tests for transfer-function calibration."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.calibration.transfer_function import TransferFunctionCalibration
from sqc.devices.transmon import TransmonQubit
from sqc.hardware.distortion import IIRDistortion, SingleExponentialDistortion


pytestmark = pytest.mark.integration


def _make_qubit():
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        n_levels=2,
    )


def test_transfer_function_calibration_rebuilds_iir_model():
    t = np.linspace(0.0, 100.0, 1001)
    true_model = SingleExponentialDistortion(amplitude=0.08, tau=15.0)
    calibration = TransferFunctionCalibration(
        qubit=_make_qubit(),
        method="simulated",
        fit_type="iir",
        transfer_model=true_model,
        t_list=t,
    )

    table = calibration.calibrate()
    rebuilt = TransferFunctionCalibration.model_from_table(table)

    omega = np.linspace(-2.0, 2.0, 101)
    assert table.kind == "transfer_function"
    assert isinstance(rebuilt, IIRDistortion)
    assert np.max(
        np.abs(
            rebuilt.frequency_response(omega)
            - true_model.discrete_frequency_response(omega, t[1] - t[0])
        )
    ) < 1e-12


def test_transfer_function_calibration_requires_measurement_or_model():
    with pytest.raises(NotImplementedError):
        TransferFunctionCalibration(qubit=_make_qubit()).calibrate()
