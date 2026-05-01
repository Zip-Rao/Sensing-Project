"""Unit tests for CalibrationTable interpolation."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.calibration.base import CalibrationTable


pytestmark = pytest.mark.unit


def test_calibration_table_evaluate_and_inverse():
    table = CalibrationTable(
        qubit_name="Q0",
        kind="linear",
        inputs=np.linspace(0, 3, 4),
        outputs=np.linspace(0, 6, 4),
    )

    assert table.evaluate(1.5) == pytest.approx(3.0)
    assert table.inverse(3.0) == pytest.approx(1.5)


def test_calibration_table_inverse_rejects_flat_outputs():
    table = CalibrationTable(
        qubit_name="Q0",
        kind="flat",
        inputs=np.array([0.0, 1.0, 2.0]),
        outputs=np.array([1.0, 1.0, 1.0]),
    )

    with pytest.raises(ValueError):
        table.inverse(1.0)
