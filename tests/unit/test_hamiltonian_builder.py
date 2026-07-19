"""Unit tests for HamiltonianBuilder.

Tests: no-flux case, lab frame, rotating frame, equivalence with
old qubit_in_mag output.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.devices.transmon import QubitSpec
from sqc.simulation.hamiltonian import HamiltonianBuilder
from sqc.control.flux_signal import FluxSignal


class TestHamiltonianBuilder:
    """HamiltonianBuilder.build correctness."""

    @pytest.fixture
    def spec(self) -> QubitSpec:
        return QubitSpec(
            name="Q0", EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
            T1=10000, T2=8000, n_levels=2,
        )

    @pytest.fixture
    def flux_signal(self) -> FluxSignal:
        return FluxSignal(
            type=2,
            t_list=np.linspace(0, 250, 500),
            amplitude=0.01,
            frequency=0.005,
        )

    def test_no_flux_lab_frame(self, spec):
        """build with flux_signal=None produces valid H_list."""
        H_list, t_global = HamiltonianBuilder.build(
            qubit=spec, flux_signal=None, pulse=None, frame="lab"
        )
        assert len(H_list) == 2
        assert len(t_global) > 0

    def test_no_flux_rotating_frame(self, spec):
        """build with flux_signal=None, rotating frame."""
        H_list, t_global = HamiltonianBuilder.build(
            qubit=spec, flux_signal=None, pulse=None,
            frame="rotating", omega_d=spec.frequency(),
        )
        # Freq_coeffs should be ~0 (omega_d == omega_01)
        assert np.allclose(H_list[1][1], 0.0, atol=1e-6)

    def test_with_flux_lab_frame(self, spec, flux_signal):
        """Lab frame H_list with flux signal."""
        H_list, t_global = HamiltonianBuilder.build(
            qubit=spec, flux_signal=flux_signal, pulse=None, frame="lab"
        )
        assert len(H_list) == 2
        # freq_coeffs should vary with time (sinusoidal flux)
        coeffs = H_list[1][1]
        assert np.std(coeffs) > 0.0  # not constant

    def test_with_flux_rotating_frame(self, spec, flux_signal):
        """Rotating frame H_list with flux signal."""
        H_list, t_global = HamiltonianBuilder.build(
            qubit=spec, flux_signal=flux_signal, pulse=None,
            frame="rotating", omega_d=spec.frequency(),
        )
        assert len(H_list) == 2
        # freq_coeffs should be centered around 0
        coeffs = H_list[1][1]
        assert abs(np.mean(coeffs)) < 0.1

    def test_equivalence_with_old_qubit_in_mag(self, spec, flux_signal):
        """HamiltonianBuilder output must match old TransmonQubit.qubit_in_mag."""
        from src.qubit import TransmonQubit as OldQ

        # Old path
        q_old = OldQ(
            EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15, T1=10000, T2=8000,
            n_levels=2,
        )
        q_old.qubit_in_mag(flux_signal, frame=0)

        # New path (using spec directly)
        H_list_new, _ = HamiltonianBuilder.build(
            qubit=spec, flux_signal=flux_signal, pulse=None, frame="lab"
        )

        # Compare H_list structure
        assert np.allclose(
            q_old.H_list[0].full(), H_list_new[0].full(), atol=1e-12
        ), "H_list[0] mismatch"
        assert np.allclose(
            q_old.H_list[1][0].full(), H_list_new[1][0].full(), atol=1e-12
        ), "H_list[1][0] operator mismatch"
        assert np.allclose(
            q_old.H_list[1][1], H_list_new[1][1], rtol=1e-12, atol=1e-12
        ), "freq_coeffs mismatch"
