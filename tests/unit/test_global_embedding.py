"""tests.unit.test_global_embedding — Tests for Phase 7 trigger + *_on API.

Covers:
- Pulse.trigger + hamiltonian_on
- CompositePulse.hamiltonian_on
- FluxSignal.trigger + samples_on
- Factory functions trigger parameter
"""
import numpy as np
import pytest

from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.control.pulse import Pulse, CompositePulse
from sqc.control.sequence import (
    create_ramsey_pulse,
    create_echo_pulse,
    create_cryoscope_pulse,
    create_cpmg_pulse,
    create_diff_echo_pulse,
    create_pulse,
    create_pi_pulse_compensation_pulse,
)

_GT = CONFIG.awg.dt
_t_rabi = np.arange(0, 40, _GT)
_t_global = CONFIG.pulse.t_global


# ---------------------------------------------------------------------------
# FluxSignal.samples_on
# ---------------------------------------------------------------------------

class TestFluxSignalGlobalEmbedding:
    """samples_on correctly projects a local signal onto t_global."""

    def test_trigger_zero_gives_full_signal(self):
        """trigger=0: samples_on returns signal padded with zeros for t > t_list[-1]."""
        sig = FluxSignal(type=2, t_list=np.arange(0, 10, 0.5),
                         amplitude=1.0, frequency=0.1, trigger=0.0)
        samples = sig.samples_on(_t_global)
        assert len(samples) == len(_t_global)
        # Within local window, values should match
        mask = _t_global <= sig.t_list[-1]
        expected = np.interp(_t_global[mask], sig.t_list, sig.samples)
        np.testing.assert_allclose(
            samples[mask], expected,
            err_msg="values inside local window should interpolate correctly"
        )
        # Outside local window, values should be zero
        assert np.all(samples[~mask] == 0.0), \
            "values outside local window should be zero"

    def test_trigger_positive_shifts_signal(self):
        """trigger=50: samples_on should have signal shifted to t >= 50."""
        offset = 50.0
        sig = FluxSignal(type=2, t_list=np.arange(0, 20, 0.5),
                         amplitude=1.0, frequency=0.05, trigger=offset)
        samples = sig.samples_on(_t_global)
        # t < offset should be zero
        assert np.all(samples[_t_global < offset] == 0.0), \
            "samples before trigger should be zero"
        # t >= offset within signal window should have non-zero
        mask = (_t_global >= offset) & (_t_global <= offset + sig.t_list[-1])
        assert np.any(samples[mask] != 0.0), \
            "samples after trigger within window should be non-zero"

    def test_trigger_negative_clips_signal(self):
        """Negative trigger: samples_on clips the tail of the signal."""
        sig = FluxSignal(type=2, t_list=np.arange(0, 40, 0.5),
                         amplitude=1.0, frequency=0.1, trigger=-20.0)
        samples = sig.samples_on(_t_global)
        # For t < 0 in global, signal should be zero
        # t_global starts at -50, so t_loc = -50 - (-20) = -30 which is < 0
        # So first points of t_global have t_loc < 0
        early_mask = _t_global < sig.trigger
        assert np.all(samples[early_mask] == 0.0), \
            "samples where t_global < trigger should be zero"

    def test_interpolation_accuracy(self):
        """Linear interpolation for points between grid points."""
        t_local = np.array([0.0, 1.0, 2.0, 3.0])
        vals = np.array([0.0, 1.0, 0.0, -1.0])
        sig = FluxSignal(type=8, t_list=t_local, signal=vals, trigger=0.0)
        t_global = np.array([0.5, 1.5, 2.5])
        samples = sig.samples_on(t_global)
        expected = np.array([0.5, 0.5, -0.5])
        np.testing.assert_allclose(samples, expected,
                                   err_msg="linear interpolation should be exact at midpoints")

    def test_default_trigger_is_zero(self):
        """Default trigger for FluxSignal is 0.0."""
        sig = FluxSignal(type=1, t_list=np.arange(0, 10, 0.5))
        assert sig.trigger == 0.0


# ---------------------------------------------------------------------------
# Pulse.hamiltonian_on
# ---------------------------------------------------------------------------

class TestPulseGlobalEmbedding:
    """hamiltonian_on correctly projects a pulse onto t_global."""

    def test_hamiltonian_on_length_matches_t_global(self):
        """hamiltonian_on output coeffs should have length len(t_global)."""
        Omega = FluxSignal(
            type=1, t_list=_t_rabi,
            amplitude=(np.pi / 2.0) / (_t_rabi[-1] - _t_rabi[0]),
        )
        p = Pulse(frame=1, omega_d=0.0, phase=0.0, Omega=Omega,
                  is_rwa=True, trigger=0.0)
        H_on = p.hamiltonian_on(_t_global)
        assert len(H_on) == 1  # RWA: single term
        _, coeffs = H_on[0]
        assert len(coeffs) == len(_t_global), \
            f"coeffs length {len(coeffs)} != t_global length {len(_t_global)}"

    def test_coeffs_zero_outside_window(self):
        """hamiltonian_on coeffs should be zero outside pulse window."""
        Omega = FluxSignal(
            type=1, t_list=_t_rabi,
            amplitude=(np.pi / 2.0) / (_t_rabi[-1] - _t_rabi[0]),
        )
        trigger = 100.0
        p = Pulse(frame=1, omega_d=0.0, phase=0.0, Omega=Omega,
                  is_rwa=True, trigger=trigger)
        H_on = p.hamiltonian_on(_t_global)
        _, coeffs = H_on[0]
        mask_inside = (_t_global >= trigger) & (_t_global <= trigger + _t_rabi[-1])
        assert np.any(coeffs[mask_inside] != 0.0), \
            "coeffs inside pulse window should be non-zero"
        assert np.all(coeffs[~mask_inside] == 0.0), \
            "coeffs outside pulse window should be zero"

    def test_lab_frame_complex_coeffs(self):
        """Lab-frame hamiltonian_on handles potentially real coeffs."""
        Omega = FluxSignal(
            type=1, t_list=_t_rabi,
            amplitude=(np.pi / 2.0) / (_t_rabi[-1] - _t_rabi[0]),
        )
        p = Pulse(frame=0, omega_d=5.0, phase=0.0, Omega=Omega,
                  is_rwa=True, trigger=0.0)
        H_on = p.hamiltonian_on(_t_global)
        assert len(H_on) >= 1
        for _, coeffs in H_on:
            assert len(coeffs) == len(_t_global)

    def test_default_trigger_is_zero(self):
        """Default trigger for Pulse is 0.0."""
        Omega = FluxSignal(type=1, t_list=_t_rabi, amplitude=0.0)
        p = Pulse(frame=1, omega_d=0.0, phase=0.0, Omega=Omega, is_rwa=True)
        assert p.trigger == 0.0

    def test_non_rwa_complex_coeffs(self):
        """Non-RWA rotating frame has complex coeffs; interp handles them."""
        Omega = FluxSignal(
            type=1, t_list=_t_rabi,
            amplitude=np.pi / (_t_rabi[-1] - _t_rabi[0]),
        )
        p = Pulse(frame=1, omega_d=0.2, phase=0.0, Omega=Omega,
                  is_rwa=False, trigger=50.0)
        H_on = p.hamiltonian_on(_t_global)
        assert len(H_on) == 3  # H_rwa_co + H_cr1 + H_cr2
        for _, coeffs in H_on:
            assert len(coeffs) == len(_t_global)


# ---------------------------------------------------------------------------
# CompositePulse.hamiltonian_on
# ---------------------------------------------------------------------------

class TestCompositePulseGlobalEmbedding:
    """CompositePulse.hamiltonian_on collects child pulse projections."""

    def test_ramsey_composite_projection(self):
        """Ramsey CompositePulse hamiltonian_on with trigger=30."""
        ctrl = create_ramsey_pulse(
            t_rabi=_t_rabi, tau=20.0, omega_d=0.0,
            phase1=np.pi / 2, phase2=0.0,
            trigger=30.0,
        )
        H_on = ctrl.hamiltonian_on(_t_global)
        # Should collect from all 3 sub-pulses (pi/2, gap, pi/2)
        assert len(H_on) == 3, \
            f"expected 3 Hamiltonian terms from 3 sub-pulses, got {len(H_on)}"
        for _, coeffs in H_on:
            assert len(coeffs) == len(_t_global)

    def test_echo_composite_projection(self):
        """Echo CompositePulse hamiltonian_on."""
        ctrl = create_echo_pulse(
            t_rabi=_t_rabi, tau=30.0, omega_d=0.0,
            trigger=10.0,
        )
        H_on = ctrl.hamiltonian_on(_t_global)
        # 5 sub-pulses: pi/2, tau, pi, tau, pi/2
        assert len(H_on) == 5
        for _, coeffs in H_on:
            assert len(coeffs) == len(_t_global)

    def test_cryoscope_composite_projection(self):
        """Cryoscope CompositePulse hamiltonian_on."""
        ctrl = create_cryoscope_pulse(
            t_rabi=_t_rabi, tau=40.0, omega_d=0.0,
            trigger=0.0,
        )
        H_on = ctrl.hamiltonian_on(_t_global)
        assert len(H_on) == 3  # pi/2, gap, pi/2
        for _, coeffs in H_on:
            assert len(coeffs) == len(_t_global)

    def test_cpmg_composite_projection(self):
        """CPMG CompositePulse hamiltonian_on."""
        ctrl = create_cpmg_pulse(
            t_rabi=_t_rabi, tau=40.0, n=2, omega_d=0.0,
            trigger=0.0,
        )
        H_on = ctrl.hamiltonian_on(_t_global)
        # pi/2 + tau/2 + n*(pi + tau) + tau/2 + pi/2 = 2 + 2*n + 2 = 4 + 2*n
        assert len(H_on) == 4 + 2 * 2  # 8 for n=2
        for _, coeffs in H_on:
            assert len(coeffs) == len(_t_global)

    def test_diff_echo_composite_projection(self):
        """Differential echo CompositePulse hamiltonian_on."""
        ctrl = create_diff_echo_pulse(
            t_rabi=_t_rabi, tau=10.0, t_int=20.0, t_rep=60.0, k=2,
            omega_d=0.0, trigger=0.0,
        )
        H_on = ctrl.hamiltonian_on(_t_global)
        # pi/2 + k*(5 pulses) + pi/2 = 2 + 5*k
        assert len(H_on) == 2 + 5 * 2  # 12 for k=2
        for _, coeffs in H_on:
            assert len(coeffs) == len(_t_global)


# ---------------------------------------------------------------------------
# Factory function trigger defaults
# ---------------------------------------------------------------------------

class TestFactoryTriggerDefaults:
    """All factory functions accept trigger=0.0 without breaking."""

    def test_create_pulse_honours_trigger(self):
        result = create_pulse(
            qubit=None, frame=1, type=1, t_list=_t_rabi,
            omega_d=0.0, phase=0.0, angle=np.pi, trigger=50.0,
        )
        assert result is not None

    def test_create_ramsey_trigger_default(self):
        ctrl = create_ramsey_pulse(t_rabi=_t_rabi, tau=20.0)
        # Default trigger=0 should create valid CompositePulse
        assert len(ctrl.pulses) == 3  # pi/2, gap, pi/2

    def test_create_pi_pulse_comp_trigger(self):
        result = create_pi_pulse_compensation_pulse(
            t_rabi=_t_rabi, T_pi=40, omega_d=0.0, trigger=10.0,
        )
        assert result is not None
