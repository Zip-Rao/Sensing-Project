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
    """Smoke tests for CryoscopeReconstruction."""

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

    def test_calib_inverse_smoke(self):
        """Linear calibration + linear phase → correct h reconstruction."""
        from sqc.reconstruction.cryoscope import CryoscopeReconstruction
        from sqc.calibration.base import CalibrationTable

        # Linear calibration: φ(h) = k * h, with k = 2.0 rad/Φ₀
        h_vals = np.array([-0.01, 0.0, 0.01])
        phi_vals = 2.0 * h_vals  # φ = 2.0 * h
        cal = CalibrationTable(
            name="test", kind="phi_h",
            inputs=h_vals, outputs=phi_vals,
        )

        # Linear phase vs time: φ(t) = α * t, where α = 0.02 rad/ns
        # Then dφ/dt = 0.02, φ_cal(h) = dφ/dt * tau = 0.02 * tau
        # So h = inverse(0.02 * tau) = (0.02 * tau) / 2.0
        tau = 50.0
        dt = 1.0
        t = np.arange(0, 100, dt)
        varphi = 0.02 * t  # linear phase accumulation

        meas = self._make_mock_measurement(varphi, t)

        recon = CryoscopeReconstruction(calibration=cal, tau=tau, method="calib_inverse")
        result = recon.reconstruct(meas, dt=dt)

        expected_h = (0.02 * tau) / 2.0  # = 0.5 Φ₀
        assert isinstance(result, FluxSignal)
        assert result.t_list.shape == t.shape
        # All values should be close to expected (linear in → linear out)
        assert np.allclose(result.signal[1:], expected_h, rtol=1e-10, atol=1e-10)

    def test_sg_diff_smoke(self):
        """SG_diff method runs and returns FluxSignal."""
        from sqc.reconstruction.cryoscope import CryoscopeReconstruction
        from sqc.calibration.base import CalibrationTable

        h_vals = np.linspace(-0.01, 0.01, 11)
        phi_vals = 2.0 * h_vals
        cal = CalibrationTable(
            name="test", kind="phi_h",
            inputs=h_vals, outputs=phi_vals,
        )

        tau = 50.0
        dt = 1.0
        t = np.arange(0, 100, dt)
        varphi = 0.02 * t

        meas = self._make_mock_measurement(varphi, t)

        recon = CryoscopeReconstruction(calibration=cal, tau=tau, method="SG_diff", sg_window=5)
        result = recon.reconstruct(meas, dt=dt)

        assert isinstance(result, FluxSignal)
        assert len(result.signal) == len(t)
        # SG should give similar result to calib_inverse for clean data
        expected_h = (0.02 * tau) / 2.0
        assert np.allclose(result.signal[5:-5], expected_h, rtol=1e-10, atol=1e-10)

    def test_diff_smoke(self):
        """diff method runs and returns FluxSignal."""
        from sqc.reconstruction.cryoscope import CryoscopeReconstruction
        from sqc.calibration.base import CalibrationTable

        h_vals = np.linspace(-0.01, 0.01, 11)
        phi_vals = 2.0 * h_vals
        cal = CalibrationTable(
            name="test", kind="phi_h",
            inputs=h_vals, outputs=phi_vals,
        )

        tau = 50.0
        dt = 1.0
        t = np.arange(0, 100, dt)
        varphi = 0.02 * t

        meas = self._make_mock_measurement(varphi, t)

        recon = CryoscopeReconstruction(calibration=cal, tau=tau, method="diff")
        result = recon.reconstruct(meas, dt=dt)

        assert isinstance(result, FluxSignal)
        assert len(result.signal) == len(t)
        expected_h = (0.02 * tau) / 2.0
        assert np.allclose(result.signal[1:-1], expected_h, rtol=1e-10, atol=1e-10)


class TestFluxResponseCalibration:
    """Verify FluxResponseCalibration methods."""

    def test_cryoscope_smoke(self):
        """method='cryoscope' runs and returns CalibrationTable(kind='phi_h')."""
        q = _make_qubit()
        from sqc.calibration.flux_response import FluxResponseCalibration

        # Use a small h_list for speed
        f = FluxResponseCalibration(
            qubit=q, method="cryoscope",
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
        from sqc.calibration.flux_response import FluxResponseCalibration

        f = FluxResponseCalibration(qubit=q, method="transient")
        with pytest.raises(NotImplementedError, match="Track B 1.2"):
            f.calibrate()


class TestTransientFrequencyCalibrationStub:
    """Verify that TransientFrequencyCalibration stub raises correctly."""

    def test_calibrate_raises_not_implemented(self):
        """TransientFrequencyCalibration.calibrate raises NotImplementedError."""
        q = _make_qubit()
        from sqc.calibration.qubit_frequency import TransientFrequencyCalibration

        tc = TransientFrequencyCalibration(qubit=q)
        with pytest.raises(NotImplementedError, match="Track B 1.2"):
            tc.calibrate()
