"""Unit tests for CryoscopeExperiment.

Smoke test: verify CryoscopeExperiment runs without error and produces
expected output shapes.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.devices.transmon import TransmonQubit
from sqc.control.flux_signal import FluxSignal
from sqc.experiments.cryoscope import CryoscopeExperiment


def _make_qubit(n_levels=2):
    """Create a fresh TransmonQubit with default test parameters."""
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=n_levels,
    )


class TestCryoscopeExperiment:
    """Smoke tests for CryoscopeExperiment."""

    def test_creates_default_signal(self):
        """CryoscopeExperiment creates default flux signal when none given."""
        q = _make_qubit()
        exp = CryoscopeExperiment(qubit=q, tau=20.0, t_rabi=np.linspace(0, 5, 10))
        assert exp.flux_signal is not None
        assert len(exp.flux_signal.t_list) == 160
        assert len(exp.trunc_list) > 0

    def test_accepts_custom_signal(self):
        """CryoscopeExperiment accepts custom flux signal and trunc_list."""
        q = _make_qubit()
        phi = FluxSignal(
            type=2,
            t_list=np.linspace(0, 40, 80),
            amplitude=0.005,
        )
        trunc_list = phi.t_list[60:10:-1]
        exp = CryoscopeExperiment(
            qubit=q,
            flux_signal=phi,
            trunc_list=trunc_list,
            tau=20.0,
            t_rabi=np.linspace(0, 5, 10),
        )
        assert exp.flux_signal is phi
        assert len(exp.trunc_list) == len(trunc_list)

    def test_run_smoke(self):
        """CryoscopeExperiment.run() completes without error.

        Uses a minimal truncation set to keep test runtime manageable.
        """
        q = _make_qubit()
        # Short signal, few truncation points
        phi = FluxSignal(
            type=2,
            t_list=np.linspace(0, 20, 40),
            amplitude=0.01,
        )
        # Only 5 truncation points for speed
        trunc_list = phi.t_list[35:20:-3][:5]
        exp = CryoscopeExperiment(
            qubit=q,
            flux_signal=phi,
            trunc_list=trunc_list,
            tau=20.0,
            t_rabi=np.linspace(0, 5, 10),
        )
        result = exp.run()

        # Check output shapes
        n_trunc = len(trunc_list)
        assert result.data["varphi"].shape == (n_trunc,)
        assert result.data["p_e_I"].shape == (n_trunc,)
        assert result.data["p_e_Q"].shape == (n_trunc,)
        assert len(result.axes["trunc"]) == n_trunc

        # varphi should be in [-pi, pi] range (arctan2)
        assert np.all(np.abs(result.data["varphi"]) <= np.pi + 1e-10)

        # p_e values should be in [0, 1]
        assert np.all(result.data["p_e_I"] >= 0) and np.all(result.data["p_e_I"] <= 1)
        assert np.all(result.data["p_e_Q"] >= 0) and np.all(result.data["p_e_Q"] <= 1)

    def test_run_with_different_omega_d(self):
        """CryoscopeExperiment with explicit omega_d."""
        q = _make_qubit()
        phi = FluxSignal(
            type=1,  # constant
            t_list=np.linspace(0, 20, 40),
            amplitude=0.01,
        )
        trunc_list = phi.t_list[30:15:-5]
        exp = CryoscopeExperiment(
            qubit=q,
            flux_signal=phi,
            trunc_list=trunc_list,
            tau=10.0,
            t_rabi=np.linspace(0, 5, 10),
            omega_d=q.frequency * 0.99,  # slight detuning
        )
        result = exp.run()
        assert len(result.data["varphi"]) == len(trunc_list)


class TestCryoscopeReconstruction:
    """Smoke tests for CryoscopeReconstruction (unified API)."""

    @staticmethod
    def _make_mock_measurement(varphi, trunc):
        """Create a minimal ExperimentResult-like object for testing."""
        from sqc.simulation.result import ExperimentResult

        return ExperimentResult(
            data={"varphi": np.asarray(varphi, dtype=float)},
            axes={"trunc": np.asarray(trunc, dtype=float)},
            metadata={},
            config={},
        )

    def test_calibration_inversion(self):
        """inversion='calibration': linear cal + linear phase → correct h."""
        from sqc.reconstruction.cryoscope import CryoscopeReconstruction
        from sqc.calibration.base import CalibrationTable

        h_vals = np.array([-0.01, 0.0, 0.01])
        phi_vals = 2.0 * h_vals
        cal = CalibrationTable(name="test", kind="phi_h", inputs=h_vals, outputs=phi_vals)

        tau = 50.0
        dt = 1.0
        t = np.arange(0, 100, dt)
        varphi = 0.02 * t

        meas = self._make_mock_measurement(varphi, t)

        recon = CryoscopeReconstruction(calibration=cal, tau=tau)
        result = recon.reconstruct(meas, dt=dt)

        expected_h = (0.02 * tau) / 2.0
        assert isinstance(result, FluxSignal)
        assert result.t_list.shape == t.shape
        assert np.allclose(result.signal[1:], expected_h, rtol=1e-10, atol=1e-10)

    def test_calibration_with_sg(self):
        """inversion='calibration' + use_sg_filter=True."""
        from sqc.reconstruction.cryoscope import CryoscopeReconstruction
        from sqc.calibration.base import CalibrationTable

        h_vals = np.linspace(-0.01, 0.01, 11)
        phi_vals = 2.0 * h_vals
        cal = CalibrationTable(name="test", kind="phi_h", inputs=h_vals, outputs=phi_vals)

        tau = 50.0
        dt = 1.0
        t = np.arange(0, 100, dt)
        varphi = 0.02 * t

        meas = self._make_mock_measurement(varphi, t)

        recon = CryoscopeReconstruction(
            calibration=cal, tau=tau,
            inversion="calibration", use_sg_filter=True, sg_window=5,
        )
        result = recon.reconstruct(meas, dt=dt)

        assert isinstance(result, FluxSignal)
        expected_h = (0.02 * tau) / 2.0
        assert np.allclose(result.signal[5:-5], expected_h, rtol=1e-10, atol=1e-10)

    def test_response_inversion(self):
        """inversion='response': uses qubit dispersion to map Δf → h."""
        from sqc.reconstruction.cryoscope import CryoscopeReconstruction
        from sqc.devices.transmon import TransmonQubit

        q = TransmonQubit(
            EC=0.2 * 2 * np.pi, EJ=10.0 * 2 * np.pi,
            T1=100e3, T2=50e3,
            flux=np.arctan(np.sqrt(2)) / np.pi,  # optimal sensitivity point
            n_levels=3,
        )

        tau = 50.0
        dt = 1.0
        t = np.arange(0, 100, dt)
        # Constant Δf = 0.01 GHz → h ≈ 0.01 / |κ| (linear regime)
        delta_f = 0.01  # GHz
        varphi = 2.0 * np.pi * delta_f * t  # φ = 2π·Δf·t

        meas = self._make_mock_measurement(varphi, t)

        recon = CryoscopeReconstruction(qubit=q, tau=tau, inversion="response")
        result = recon.reconstruct(meas, dt=dt)

        assert isinstance(result, FluxSignal)
        assert len(result.signal) == len(t)
        # dφ/dt = 2π·Δf, so Δf reconstructed ≈ 0.01 GHz
        # h should be positive and small (qubit at work point with negative sensitivity)
        assert np.all(np.abs(result.signal[1:-1]) < 0.1)


class TestFluxResponseCalibration:
    """Verify FluxResponseCalibration methods."""

    def test_cryoscope_smoke(self):
        """method='cryoscope' runs and returns CalibrationTable(kind='phi_h')."""
        q = _make_qubit()
        from sqc.reconstruction.cryoscope import CryoscopeCalibration

        # Use a small h_list for speed
        f = CryoscopeCalibration(
            qubit=q,
            h_list=np.linspace(-0.01, 0.01, 3),
            tau=20.0,
        )
        table = f.calibrate()

        assert table.kind == "phi_h"
        assert len(table.inputs) == 3
        assert len(table.outputs) == 3
        assert table.fit_params["method"] == "cryoscope"
        assert table.fit_params["tau"] == 20.0

    def test_transient_stub_raises(self):
        """method='transient' raises NotImplementedError."""
        q = _make_qubit()
        from sqc.calibration.frequency import FluxResponseCalibration

        f = FluxResponseCalibration(qubit=q, method="transient")
        with pytest.raises(NotImplementedError, match="Track B 1.2"):
            f.calibrate()


class TestTransientFrequencyCalibrationStub:
    """Verify that TransientFrequencyCalibration stub raises correctly."""

    def test_calibrate_raises_not_implemented(self):
        """TransientFrequencyCalibration.calibrate raises NotImplementedError."""
        q = _make_qubit()
        from sqc.calibration.frequency import SinglePointFrequencyCalibration

        tc = SinglePointFrequencyCalibration(qubit=q, method="transient")
        with pytest.raises(NotImplementedError, match="Track B 1.2"):
            tc.calibrate()
