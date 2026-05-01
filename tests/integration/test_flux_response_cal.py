"""Integration tests for flux-response calibration."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.calibration.flux_response import FluxResponseCalibration
from sqc.calibration.qubit_frequency import TransientFrequencyCalibration
from sqc.devices.transmon import TransmonQubit


pytestmark = pytest.mark.integration


def _make_qubit():
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000,
        T2=8000,
        n_levels=2,
    )


def test_flux_response_cryoscope_returns_table():
    q = _make_qubit()
    cal = FluxResponseCalibration(
        qubit=q,
        method="cryoscope",
        h_list=np.linspace(-0.01, 0.01, 5),
        tau=20.0,
    )
    table = cal.calibrate()

    assert table.kind == "phi_h"
    assert table.inputs.shape == (5,)
    assert table.outputs.shape == (5,)
    assert np.all(np.isfinite(table.outputs))


def test_transient_frequency_calibration_reports_track_b_blocker():
    with pytest.raises(NotImplementedError):
        TransientFrequencyCalibration(qubit=_make_qubit()).calibrate()
