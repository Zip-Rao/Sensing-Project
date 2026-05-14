"""Unit tests for PiPulseCompensationExperiment.

Smoke test: verify PiPulseCompensationExperiment runs without error
and produces expected 2D output shapes.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.devices.transmon import TransmonQubit
from sqc.control.flux_signal import FluxSignal
from sqc.experiments.pi_pulse_comp import PiPulseCompensationExperiment


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


class TestPiPulseCompensationExperiment:
    """Smoke tests for PiPulseCompensationExperiment."""

    def test_creates_default_signal(self):
        """Creates default flux signal when none given."""
        q = _make_qubit()
        exp = PiPulseCompensationExperiment(
            qubit=q, T_pi=5.0, t_rabi=np.linspace(0, 5, 10),
        )
        assert exp.flux_signal is not None
        assert len(exp.flux_signal.t_list) > 0
        assert len(exp.tau_list) > 0
        assert len(exp.z_list) > 0

    def test_accepts_custom_parameters(self):
        """Accepts custom tau_list, z_list, T_pi."""
        q = _make_qubit()
        phi = FluxSignal(
            type=8,
            t_list=np.linspace(0, 40, 80),
            signal=np.zeros(80),
        )
        tau_list = np.arange(0, 20, 2.0)
        z_list = np.linspace(-0.005, 0.005, 7)
        exp = PiPulseCompensationExperiment(
            qubit=q,
            flux_signal=phi,
            tau_list=tau_list,
            z_list=z_list,
            T_pi=5.0,
            t_rabi=np.linspace(0, 5, 10),
        )
        assert exp.flux_signal is phi
        assert len(exp.tau_list) == len(tau_list)
        assert len(exp.z_list) == len(z_list)

    def test_run_smoke(self):
        """run() completes, produces expected 2D output shapes."""
        q = _make_qubit()
        phi = FluxSignal(
            type=8,
            t_list=np.linspace(0, 30, 60),
            signal=np.zeros(60),
        )
        tau_list = np.arange(0, 12, 4.0)[:3]   # 3 points
        z_list = np.linspace(-0.005, 0.005, 5)  # 5 points
        exp = PiPulseCompensationExperiment(
            qubit=q,
            flux_signal=phi,
            tau_list=tau_list,
            z_list=z_list,
            T_pi=5.0,
            t_rabi=np.linspace(0, 5, 10),
        )
        result = exp.run()

        # Check shapes
        assert result.data["p_e"].shape == (3, 5)
        assert result.data["z_star"].shape == (3,)
        assert len(result.axes["tau"]) == 3
        assert len(result.axes["z"]) == 5

        # p_e values in [0, 1]
        assert np.all(result.data["p_e"] >= 0)
        assert np.all(result.data["p_e"] <= 1)

        # z_star should be within z_list range
        assert np.all(result.data["z_star"] >= z_list[0])
        assert np.all(result.data["z_star"] <= z_list[-1])

    def test_run_with_explicit_omega_bias(self):
        """run() with explicit omega_bias."""
        q = _make_qubit()
        phi = FluxSignal(
            type=8,
            t_list=np.linspace(0, 30, 60),
            signal=np.zeros(60),
        )
        tau_list = np.arange(0, 6, 4.0)[:2]
        z_list = np.linspace(-0.005, 0.005, 3)
        exp = PiPulseCompensationExperiment(
            qubit=q,
            flux_signal=phi,
            tau_list=tau_list,
            z_list=z_list,
            T_pi=5.0,
            t_rabi=np.linspace(0, 5, 10),
            omega_bias=q.frequency * 1.01,
        )
        result = exp.run()
        assert result.data["p_e"].shape == (2, 3)
