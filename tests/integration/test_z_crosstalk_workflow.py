"""tests.integration.test_z_crosstalk_workflow — Integration tests for ZCrosstalkWorkflow.

Verifies:
- Algorithmic correctness (H_BA extraction, compensation) with synthetic data.
- End-to-end workflow structural correctness (runs without error).
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import Waveform
from sqc.control.flux_signal import FluxSignal
from sqc.devices.chip import ChipTopology
from sqc.devices.transmon import TransmonQubit
from sqc.hardware.transfer_matrix import TransferMatrix


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _make_minimal_chip(n_levels: int = 2) -> ChipTopology:
    """Build minimal 2-qubit chip for testing."""
    qA = TransmonQubit(
        EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
        T1=10000, T2=8000, flux=0.0, n_levels=n_levels, name="QA",
    )
    qB = TransmonQubit(
        EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
        T1=10000, T2=8000, flux=0.0, n_levels=n_levels, name="QB",
    )
    return ChipTopology(qubits=[qA, qB])


def _make_flux_pulse(n: int = 160) -> Waveform:
    """Create a step-like flux pulse."""
    t = np.linspace(0, 80, n)
    samples = np.where((t > 20) & (t < 60), 1.0, 0.0)
    return Waveform(t_list=t, samples=samples)


# ---------------------------------------------------------------------------
# fast tests (algorithmic only, no mesolve)
# ---------------------------------------------------------------------------

class TestZCrosstalkWorkflowAlgorithmic:
    """Verify algorithmic components with synthetic / perfect data."""

    def test_construction(self):
        """Workflow can be constructed."""
        from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

        chip = _make_minimal_chip()
        pulse = _make_flux_pulse()
        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[1.0, 0.0], [0.03, 1.0]]),
            source_names=["QA", "QB"],
            target_names=["QA", "QB"],
        )

        wf = ZCrosstalkWorkflow(
            chip=chip,
            flux_pulse_on_A=pulse,
            true_transfer_matrix=tm,
            qubit_A_name="QA",
            qubit_B_name="QB",
        )
        assert wf.qubit_A_name == "QA"
        assert wf.qubit_B_name == "QB"

    def test_H_BA_extraction_from_perfect_reconstruction(self):
        """H_BA extracted from perfect phi_B matches truth at DC.

        Uses TransferMatrix.apply() to generate perfect phi_B_true
        (no Wiener noise), then verifies H_BA extraction accuracy.
        """
        from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

        chip = _make_minimal_chip()
        pulse = _make_flux_pulse()
        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[1.0, 0.0], [0.04, 1.0]]),
            source_names=["QA", "QB"],
            target_names=["QA", "QB"],
        )

        wf = ZCrosstalkWorkflow(
            chip=chip,
            flux_pulse_on_A=pulse,
            true_transfer_matrix=tm,
            qubit_A_name="QA",
            qubit_B_name="QB",
            deconv_lambda_reg=1e-6,  # near-zero regularization
        )

        # Compute perfect phi_B from transfer matrix
        fluxes = tm.apply({"QA": pulse})
        phi_B_perfect = fluxes["QB"]

        # Extract H_BA
        H_est, omega = wf._extract_H_BA(phi_B_perfect)
        H_true = wf._interpolate_H_true(omega)

        # Check DC component
        dc_idx = np.argmin(np.abs(omega))
        dc_est = np.abs(H_est[dc_idx])
        dc_true = np.abs(H_true[dc_idx])
        err_dc = abs(dc_est - dc_true) / (dc_true + 1e-12)
        assert err_dc < 0.02, f"DC H_BA error {err_dc:.4f} exceeds 2%"

        # Check fit error at significant frequencies
        V_A_fft = np.fft.fft(np.asarray(pulse.samples))
        V_power = np.abs(V_A_fft) ** 2
        significant = V_power > 0.01 * V_power.max()
        if np.sum(significant) > 0:
            err_sig = wf._compute_fit_error(H_est[significant], H_true[significant])
            assert err_sig < -20.0, (
                f"Fit error at significant bins: {err_sig:.2f} dB exceeds -20 dB"
            )

    def test_H_BA_from_resampled_perfect_reconstruction(self):
        """H_BA extraction handles different time grids via resampling.

        Tests _extract_H_BA when phi_B_rec has a different time grid
        from V_A (simulating Wiener output mismatch).
        """
        from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

        chip = _make_minimal_chip()
        pulse = _make_flux_pulse(n=100)  # V_A has 100 points
        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[1.0, 0.0], [0.05, 1.0]]),
            source_names=["QA", "QB"],
            target_names=["QA", "QB"],
        )

        wf = ZCrosstalkWorkflow(
            chip=chip,
            flux_pulse_on_A=pulse,
            true_transfer_matrix=tm,
            qubit_A_name="QA",
            qubit_B_name="QB",
            deconv_lambda_reg=1e-6,
        )

        # Create phi_B on a DIFFERENT time grid (simulating Wiener output)
        fluxes = tm.apply({"QA": pulse})
        phi_B_perfect = fluxes["QB"]  # has 100 points on same grid

        # Create a version with different time grid (shorter, different dt)
        t_alt = np.linspace(5, 75, 70)  # offset and different spacing
        phi_arr = np.interp(t_alt, phi_B_perfect.t_list, phi_B_perfect.samples)
        phi_B_alt = FluxSignal(type=8, t_list=t_alt, signal=phi_arr)

        # Extract H_BA — should handle resampling automatically
        H_est, omega = wf._extract_H_BA(phi_B_alt)
        H_true = wf._interpolate_H_true(omega)

        dc_idx = np.argmin(np.abs(omega))
        dc_est = np.abs(H_est[dc_idx])
        dc_true = np.abs(H_true[dc_idx])
        err_dc = abs(dc_est - dc_true) / (dc_true + 1e-12)
        assert err_dc < 0.15, (
            f"DC H_BA error after resampling: {err_dc:.4f} exceeds 15% "
            "(interpolation introduces some error)"
        )

    def test_fit_error_metric(self):
        """_compute_fit_error behaves correctly (better estimate → lower dB)."""
        from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

        n = 100
        H_true = np.ones(n, dtype=complex) * 0.03
        H_bad = np.ones(n, dtype=complex) * 0.06  # 2x off
        H_good = np.ones(n, dtype=complex) * 0.0301  # 0.33% off

        err_bad = ZCrosstalkWorkflow._compute_fit_error(H_bad, H_true)
        err_good = ZCrosstalkWorkflow._compute_fit_error(H_good, H_true)

        assert err_good < err_bad  # good estimate has lower dB error
        assert err_good < 0.0  # less than 0 dB = better than factor-of-2

    def test_compensation_factor_with_perfect_reconstruction(self):
        """Compensation with perfect reconstruction gives factor > 5."""
        from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

        chip = _make_minimal_chip()
        pulse = _make_flux_pulse()
        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[1.0, 0.0], [0.05, 1.0]]),
            source_names=["QA", "QB"],
            target_names=["QA", "QB"],
        )

        wf = ZCrosstalkWorkflow(
            chip=chip,
            flux_pulse_on_A=pulse,
            true_transfer_matrix=tm,
            qubit_A_name="QA",
            qubit_B_name="QB",
        )
        wf._qubit_A = chip.get_qubit("QA")
        wf._qubit_B = chip.get_qubit("QB")

        fluxes = tm.apply({"QA": pulse})
        wf._phi_A = fluxes["QA"]
        wf._phi_B_true = fluxes["QB"]

        result = wf._compensate(wf._phi_B_true)  # perfect rec = perfect comp
        assert result["compensation_factor"] > 100.0, (
            f"Perfect compensation should give very high factor, got {result['compensation_factor']}"
        )
        assert result["parasitic_phase_compensated"] < 1e-10


# ---------------------------------------------------------------------------
# end-to-end test (runs mesolve, marked slow)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestZCrosstalkWorkflowEndToEnd:
    """Full end-to-end ZCrosstalkWorkflow.run() test.

    Uses moderate parameters for reasonable runtime.
    The Wiener reconstruction quality depends heavily on kernel length,
    so accuracy requirements here are relaxed — algorithmic accuracy
    is verified in TestZCrosstalkWorkflowAlgorithmic above.
    """

    def test_full_workflow_runs(self):
        """Verify workflow completes and produces correctly structured output."""
        from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

        np.random.seed(42)

        chip = _make_minimal_chip(n_levels=2)

        t_pulse = np.linspace(0, 60, 60)
        samples = np.where((t_pulse > 12) & (t_pulse < 48), 1.0, 0.0)
        pulse = Waveform(t_list=t_pulse, samples=samples)

        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[1.0, 0.0], [0.04, 1.0]]),
            source_names=["QA", "QB"],
            target_names=["QA", "QB"],
        )

        wf = ZCrosstalkWorkflow(
            chip=chip,
            flux_pulse_on_A=pulse,
            true_transfer_matrix=tm,
            qubit_A_name="QA",
            qubit_B_name="QB",
            t_rabi=np.linspace(0, 10, 8),
            wiener_lambda_reg=1e-2,
            deconv_lambda_reg=1e-3,
        )

        result = wf.run()

        # ---- structural checks ----
        required_keys = [
            "phi_A", "phi_B_true", "phi_B_reconstructed",
            "H_BA_estimated", "H_BA_true", "omega_grid",
            "fit_error_dB", "parasitic_phase_uncompensated",
            "parasitic_phase_compensated", "compensation_factor",
            "phi_B_after_compensation",
        ]
        for key in required_keys:
            assert key in result, f"Missing key: {key}"

        assert isinstance(result["phi_B_true"], FluxSignal)
        assert isinstance(result["phi_B_reconstructed"], FluxSignal)
        assert isinstance(result["H_BA_estimated"], np.ndarray)
        assert isinstance(result["H_BA_true"], np.ndarray)
        assert len(result["H_BA_estimated"]) == len(result["omega_grid"])
        assert len(result["H_BA_true"]) == len(result["omega_grid"])
        assert isinstance(result["compensation_factor"], float)
        assert result["compensation_factor"] > 0  # positive

        # ---- basic sanity: phi_B_reconstructed has nonzero content ----
        assert np.max(np.abs(result["phi_B_reconstructed"].samples)) > 1e-12

        # ---- DC H_BA sanity: should not be complete garbage ----
        H_est = result["H_BA_estimated"]
        omega = result["omega_grid"]
        dc_idx = np.argmin(np.abs(omega))
        # DC value should be nonzero (at least 1% of expected)
        assert np.abs(H_est[dc_idx]) > 0.0004, (
            f"DC H_BA estimate too small: {np.abs(H_est[dc_idx]):.6f}"
        )

    def test_multiple_crosstalk_levels(self):
        """Workflow works for different crosstalk coefficients."""
        from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

        np.random.seed(42)

        for crosstalk in [0.02, 0.06, 0.10]:
            chip = _make_minimal_chip(n_levels=2)
            t_pulse = np.linspace(0, 60, 60)
            samples = np.where((t_pulse > 12) & (t_pulse < 48), 1.0, 0.0)
            pulse = Waveform(t_list=t_pulse, samples=samples)

            tm = TransferMatrix.from_dc_matrix(
                dc_matrix=np.array([[1.0, 0.0], [crosstalk, 1.0]]),
                source_names=["QA", "QB"],
                target_names=["QA", "QB"],
            )

            wf = ZCrosstalkWorkflow(
                chip=chip,
                flux_pulse_on_A=pulse,
                true_transfer_matrix=tm,
                qubit_A_name="QA",
                qubit_B_name="QB",
                t_rabi=np.linspace(0, 10, 8),
                wiener_lambda_reg=1e-2,
                deconv_lambda_reg=1e-3,
            )

            result = wf.run()
            assert result["phi_B_true"].samples.max() == pytest.approx(
                crosstalk, rel=0.01
            ), f"phi_B_true max mismatch for crosstalk={crosstalk}"
            assert result["compensation_factor"] > 0
