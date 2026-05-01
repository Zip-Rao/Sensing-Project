"""Unit tests for QubitSpec — frozen parameter record.

Tests: immutability, frequency/anharmonicity/sensitivity, with_flux,
optimal_work_point units.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.devices.transmon import QubitSpec


class TestQubitSpec:
    """QubitSpec correctness and immutability."""

    @pytest.fixture
    def spec(self) -> QubitSpec:
        return QubitSpec(
            name="Q0", EC=2 * np.pi * 0.2, EJ=2 * np.pi * 15,
            T1=10000, T2=8000,
        )

    # -- immutability -------------------------------------------------------

    def test_frozen_dataclass(self, spec):
        """QubitSpec is frozen; setting an attribute must raise."""
        with pytest.raises(Exception):
            spec.EC = 2 * np.pi * 0.3  # type: ignore[misc]

    def test_with_flux_returns_new_object(self, spec):
        """with_flux must return a different object with updated flux_bias."""
        spec2 = spec.with_flux(0.1)
        assert spec2 is not spec
        assert spec2.flux_bias == 0.1
        assert spec.flux_bias == 0.0  # original unchanged

    # -- spectroscopic values -----------------------------------------------

    def test_frequency_value(self, spec):
        """Frequency should be ~4.7 GHz for default params."""
        f = spec.frequency() / (2 * np.pi)
        assert 4.0 < f < 8.0, f"f_01={f} GHz out of recommended range"
        expected = np.sqrt(8 * spec.EJ * spec.EC) - spec.EC
        assert f == pytest.approx(expected / (2 * np.pi), rel=1e-12)

    def test_anharmonicity_value(self, spec):
        """Anharmonicity = -EC."""
        assert spec.anharmonicity() == pytest.approx(-spec.EC, rel=1e-12)
        assert spec.anharmonicity() / (2 * np.pi) == pytest.approx(-0.2, rel=1e-6)

    def test_sensitivity_at_zero(self, spec):
        """Sensitivity at zero flux should be small (symmetric point)."""
        sens = spec.sensitivity(flux=0.0)
        # At flux=0, derivative of |cos(pi*Phi)| is 0
        assert abs(sens) < 1e-6

    def test_sensitivity_nonzero_flux(self, spec):
        """Sensitivity should be non-zero away from sweet spot."""
        sens = spec.sensitivity(flux=0.15)
        assert abs(sens) > 0.01  # significant sensitivity away from zero

    # -- optimal work point -------------------------------------------------

    def test_optimal_work_point_units(self, spec):
        """optimal_work_point must return a value in Phi_0 (not radians)."""
        owp = spec.optimal_work_point()
        # arctan(sqrt(2)) / pi ≈ 0.304 in Phi_0
        expected = np.arctan(np.sqrt(2)) / np.pi
        assert owp == pytest.approx(expected, rel=1e-12)
        assert 0.1 < owp < 0.5  # sanity: must be in Phi_0 range

    # -- EJ_at -------------------------------------------------------------

    def test_EJ_at_zero_flux(self, spec):
        """EJ at zero flux should equal nominal EJ."""
        assert spec.EJ_at(0.0) == pytest.approx(spec.EJ, rel=1e-12)

    def test_EJ_at_half_flux(self, spec):
        """EJ at Phi=0.5 should be zero."""
        assert spec.EJ_at(0.5) == pytest.approx(0.0, abs=1e-12)
