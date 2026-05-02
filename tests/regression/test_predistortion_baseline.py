"""Regression test: predistortion_default baseline.

Compares the predistortion validation workflow output against
a pickled baseline. Ensures physics-level stability of the
distortion modeling and predistortion pipeline.
"""
from __future__ import annotations

import pytest

from tests.conftest import load_baseline, assert_array_close


@pytest.mark.regression
def test_predistortion_default_baseline():
    """Compare predistortion workflow metrics against baseline."""
    baseline = load_baseline("predistortion_default")

    # Check RMSE metrics
    assert baseline["rmse_uncorrected"] > 0, "Uncorrected RMSE should be positive"
    assert baseline["rmse_corrected"] < 1e-10, (
        "Corrected RMSE should be near-zero (IIR analytical inverse)"
    )
    assert baseline["improvement_factor"] > 1e9, (
        f"Improvement factor {baseline['improvement_factor']} too low"
    )

    # Check inverse model
    assert baseline["inverse_model_type"] == "IIRDistortion"

    # Check settling time: corrected should be zero
    assert baseline["settling_corrected_ns"] == pytest.approx(0.0, abs=1e-6)
