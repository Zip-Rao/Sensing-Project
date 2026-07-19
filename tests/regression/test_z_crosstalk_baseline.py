"""tests.regression.test_z_crosstalk_baseline — Z-crosstalk regression baseline.

Verifies that the ZCrosstalkWorkflow produces consistent numerical
results for a fixed set of parameters.
"""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import assert_array_close, load_baseline


@pytest.fixture(scope="module")
def z_crosstalk_baseline():
    """Load z_crosstalk_default baseline."""
    return load_baseline("z_crosstalk_default")


@pytest.mark.slow
def test_z_crosstalk_default_baseline(z_crosstalk_baseline):
    """Z-crosstalk workflow matches baseline."""
    from sqc.control.waveform import Waveform
    from sqc.devices.chip import ChipTopology
    from sqc.devices.transmon import TransmonQubit
    from sqc.hardware.transfer_matrix import TransferMatrix
    from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

    np.random.seed(42)

    # Reconstruct qubits with exact baseline parameters
    qA = TransmonQubit(
        EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
        T1=10000, T2=8000, flux=0.0, n_levels=2, name="QA",
    )
    qB = TransmonQubit(
        EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
        T1=10000, T2=8000, flux=0.0, n_levels=2, name="QB",
    )
    chip = ChipTopology(qubits=[qA, qB])

    # Reconstruct pulse from baseline
    t_pulse = z_crosstalk_baseline["t_pulse"]
    pulse_samples = z_crosstalk_baseline["pulse_samples"]
    pulse = Waveform(t_list=t_pulse, samples=pulse_samples)

    # Reconstruct transfer matrix
    tm = TransferMatrix.from_dc_matrix(
        dc_matrix=np.array([[1.0, 0.0], [0.04, 1.0]]),
        source_names=["QA", "QB"],
        target_names=["QA", "QB"],
    )

    # Reconstruct t_rabi
    t_rabi = z_crosstalk_baseline["t_rabi"]

    wf = ZCrosstalkWorkflow(
        chip=chip,
        flux_pulse_on_A=pulse,
        true_transfer_matrix=tm,
        qubit_A_name="QA",
        qubit_B_name="QB",
        t_rabi=t_rabi,
        wiener_lambda_reg=z_crosstalk_baseline["wiener_lambda_reg"],
        deconv_lambda_reg=z_crosstalk_baseline["deconv_lambda_reg"],
    )

    result = wf.run()

    # Compare key outputs
    assert_array_close(
        result["phi_B_true"].samples,
        z_crosstalk_baseline["phi_B_true_samples"],
        name="phi_B_true",
    )

    assert_array_close(
        result["phi_B_reconstructed"].samples,
        z_crosstalk_baseline["phi_B_reconstructed_samples"],
        name="phi_B_reconstructed",
    )

    assert_array_close(
        result["H_BA_estimated"].real,
        z_crosstalk_baseline["H_BA_estimated_real"],
        name="H_BA_estimated_real",
    )

    assert_array_close(
        result["H_BA_estimated"].imag,
        z_crosstalk_baseline["H_BA_estimated_imag"],
        name="H_BA_estimated_imag",
    )

    assert result["compensation_factor"] == pytest.approx(
        z_crosstalk_baseline["compensation_factor"],
        rel=1e-6,
        abs=1e-9,
    )

    assert result["fit_error_dB"] == pytest.approx(
        z_crosstalk_baseline["fit_error_dB"],
        rel=1e-6,
        abs=1e-9,
    )
