"""Unit tests for DelayRamseyExperiment.

Smoke test: verify DelayRamseyExperiment runs without error and produces
expected output shapes. Integration test: calibration → experiment → reconstruction.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.devices.transmon import TransmonQubit
from sqc.control.flux_signal import FluxSignal
from sqc.experiments.delay_ramsey import DelayRamseyExperiment


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


class TestDelayRamseyExperiment:
    """Smoke tests for DelayRamseyExperiment."""

    def test_creates_default_signal(self):
        """Creates default flux signal when none given."""
        q = _make_qubit()
        exp = DelayRamseyExperiment(
            qubit=q, tau_R=10.0, t_rabi=np.linspace(0, 5, 10),
        )
        assert exp.flux_signal is not None
        assert len(exp.flux_signal.t_list) > 0
        assert len(exp.t_d_list) > 0

    def test_accepts_custom_signal(self):
        """Accepts custom FluxSignal and t_d_list."""
        q = _make_qubit()
        phi = FluxSignal(
            type=8,
            t_list=np.linspace(0, 40, 80),
            signal=np.zeros(80),
        )
        t_d_list = np.arange(0, 50, 2.0)
        exp = DelayRamseyExperiment(
            qubit=q,
            flux_signal=phi,
            t_d_list=t_d_list,
            tau_R=10.0,
            t_rabi=np.linspace(0, 5, 10),
        )
        assert exp.flux_signal is phi
        assert len(exp.t_d_list) == len(t_d_list)

    def test_run_smoke(self):
        """run() completes without error, produces expected output shapes."""
        q = _make_qubit()
        phi = FluxSignal(
            type=8,
            t_list=np.linspace(0, 30, 60),
            signal=np.zeros(60),
        )
        t_d_list = np.arange(0, 10, 2.0)[:5]
        exp = DelayRamseyExperiment(
            qubit=q,
            flux_signal=phi,
            t_d_list=t_d_list,
            tau_R=10.0,
            t_rabi=np.linspace(0, 5, 10),
            run_baseline=True,
        )
        result = exp.run()

        n_td = len(t_d_list)
        assert result.data["varphi"].shape == (n_td,)
        assert result.data["varphi_raw"].shape == (n_td,)
        assert result.data["p_e_I"].shape == (n_td,)
        assert result.data["p_e_Q"].shape == (n_td,)
        assert len(result.axes["t_d"]) == n_td

        # p_e values in [0, 1]
        assert np.all(result.data["p_e_I"] >= 0)
        assert np.all(result.data["p_e_I"] <= 1)
        assert np.all(result.data["p_e_Q"] >= 0)
        assert np.all(result.data["p_e_Q"] <= 1)

    def test_run_no_baseline(self):
        """run() with run_baseline=False still completes."""
        q = _make_qubit()
        phi = FluxSignal(
            type=8,
            t_list=np.linspace(0, 30, 60),
            signal=np.zeros(60),
        )
        t_d_list = np.arange(0, 10, 2.0)[:3]
        exp = DelayRamseyExperiment(
            qubit=q,
            flux_signal=phi,
            t_d_list=t_d_list,
            tau_R=10.0,
            t_rabi=np.linspace(0, 5, 10),
            run_baseline=False,
        )
        result = exp.run()
        assert result.data["varphi_base"] is None
        # varphi should equal varphi_raw
        np.testing.assert_array_almost_equal(
            result.data["varphi"], result.data["varphi_raw"],
        )


class TestDelayRamseyCalibration:
    """Smoke tests for DelayRamseyCalibration."""

    def test_calibrate_smoke(self):
        """calibrate() completes and returns CalibrationTable."""
        from sqc.reconstruction.delay_ramsey import DelayRamseyCalibration

        q = _make_qubit()
        cal = DelayRamseyCalibration(
            qubit=q,
            z_list=np.linspace(-0.01, 0.01, 5),
            tau_R=10.0,
            t_rabi=np.linspace(0, 5, 10),
        )
        table = cal.calibrate()
        assert table.kind == "phi_z"
        assert len(table.inputs) == 5
        assert len(table.outputs) == 5
        assert "k" in table.fit_params


class TestTailReconstruction:
    """Smoke tests for TailReconstruction."""

    def test_delay_ramsey_reconstruct(self):
        """TailReconstruction with delay_ramsey method via calibration."""
        from sqc.reconstruction.delay_ramsey import DelayRamseyReconstruction
        from sqc.reconstruction.pi_pulse_comp import PiPulseCompReconstruction
        from sqc.calibration.base import CalibrationTable
        from sqc.simulation.result import ExperimentResult

        # Synthetic calibration: phi = k * z with k = 2.0
        cal = CalibrationTable(
            name="test",
            kind="phi_z",
            inputs=np.array([0.0, 0.01]),
            outputs=np.array([0.0, 0.02]),
            fit_params={"k": 2.0, "tau_R": 20.0},
        )

        # Synthetic measurement: phi = [0, 0.01, 0.02]
        meas = ExperimentResult(
            data={"varphi": np.array([0.0, 0.01, 0.02])},
            axes={"t_d": np.array([0.0, 5.0, 10.0])},
            config={"tau_R": 20.0},
        )

        recon = DelayRamseyReconstruction(
            calibration=cal, inversion="calibration",
        )
        flux = recon.reconstruct(meas)

        assert flux.type == 8
        np.testing.assert_array_almost_equal(flux.signal, np.array([0.0, 0.005, 0.01]))

    def test_pi_pulse_comp_reconstruct(self):
        """TailReconstruction with pi_pulse_comp method."""
        from sqc.reconstruction.delay_ramsey import DelayRamseyReconstruction
        from sqc.reconstruction.pi_pulse_comp import PiPulseCompReconstruction
        from sqc.simulation.result import ExperimentResult

        # z_star = [0.005, 0.003, 0.001] → flux = -z_star
        meas = ExperimentResult(
            data={"z_star": np.array([0.005, 0.003, 0.001])},
            axes={"tau": np.array([0.0, 10.0, 20.0])},
        )

        recon = PiPulseCompReconstruction()
        flux = recon.reconstruct(meas)

        assert flux.type == 8
        np.testing.assert_array_almost_equal(
            flux.signal, np.array([-0.005, -0.003, -0.001]),
        )
