"""Unit tests for KernelEstimator."""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import assert_array_close


@pytest.fixture
def qubit():
    from src.qubit import TransmonQubit

    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=2,
    )


def test_kernel_estimator_default_params():
    """KernelEstimator with default parameters can be created."""
    from sqc.reconstruction.kernel import KernelEstimator

    ke = KernelEstimator()
    assert ke.stim_amplitude == 0.0215
    assert ke.stim_width == 3.0
    assert not ke.auto_calibrate


def test_kernel_estimator_estimate_ramsey_pulse(qubit):
    """KernelEstimator.estimate returns correct shape for a ramsey pulse."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse

    t_rabi = np.linspace(0, 10, 10)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    ke = KernelEstimator()
    t_samples, kernel = ke.estimate(pulse, qubit)

    assert isinstance(t_samples, np.ndarray)
    assert isinstance(kernel, np.ndarray)
    assert len(t_samples) == len(kernel)
    assert len(t_samples) > 0
    # Kernel values should be finite
    assert np.all(np.isfinite(kernel))


def test_kernel_estimator_estimate_nonzero_kernel(qubit):
    """KernelEstimator produces non-trivial (non-zero) kernel values."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse

    t_rabi = np.linspace(0, 10, 20)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    ke = KernelEstimator()
    t_samples, kernel = ke.estimate(pulse, qubit)

    # Kernel should not be all zeros
    assert np.max(np.abs(kernel)) > 0.0
    # There should be variation across time points
    assert np.max(kernel) != np.min(kernel) or len(kernel) == 1


def test_kernel_estimator_matches_legacy(qubit):
    """KernelEstimator output matches legacy CompositePulse.get_kernel."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    from sqc.control.pulse import CompositePulse

    t_rabi = np.linspace(0, 10, 20)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    # Legacy get_kernel
    pulse.get_kernel(qubit)
    legacy_ts = np.array(pulse.t_samples)
    legacy_k = np.array(pulse.kernel)

    # New KernelEstimator
    ke = KernelEstimator()
    new_ts, new_k = ke.estimate(pulse, qubit)

    assert_array_close(new_k, legacy_k, name="kernel_vs_legacy")
    assert_array_close(new_ts, legacy_ts, name="t_samples_vs_legacy")


# ═════════════════════════════════════════════════════════════════════════════
# Phase 10.1 — mode / virtual_z_impl / KernelResult
# ═════════════════════════════════════════════════════════════════════════════


def test_flux_mode_requires_qubit():
    """KernelEstimator(mode='flux') with qubit=None raises ValueError."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 5)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=0.0)

    ke = KernelEstimator(mode='flux')
    with pytest.raises(ValueError, match="requires a qubit"):
        ke.estimate(pulse, None)


def test_omega_mode_requires_qubit():
    """KernelEstimator(mode='omega') with qubit=None raises ValueError."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 5)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=0.0)

    ke = KernelEstimator(mode='omega')
    with pytest.raises(ValueError, match="requires a qubit"):
        ke.estimate(pulse, None)


def test_omega_math_vz_kernel_shape(qubit):
    """Omega mode + math VZ produces correct shaped output."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 10)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    ke = KernelEstimator(mode='omega', virtual_z_impl='math')
    t_samples, kernel = ke.estimate(pulse, qubit)

    assert isinstance(t_samples, np.ndarray)
    assert isinstance(kernel, np.ndarray)
    assert len(t_samples) == len(kernel)
    assert len(t_samples) > 0
    # Kernel values should be finite
    assert np.all(np.isfinite(kernel))
    # Omega kernel should be non-trivial
    assert np.max(np.abs(kernel)) > 0.0


def test_flux_omega_kappa_conversion():
    """In linear regime, ∫k_flux dt ≈ κ · ∫k_omega dt.

    κ = dω/dΦ (frequency_sensitivity) evaluated at the qubit's
    operating point.  The relationship holds because a flux area A
    produces a frequency shift κ·A, which over the pulse duration
    produces a phase shift φ = κ·A.  Hence:

        ∫k_flux dt  ≈  κ · ∫k_omega dt

    NOTE: Uses a qubit at flux=0.15 (away from sweet spot), because
    at flux=0, κ = dω/dΦ = 0.
    """
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    from src.qubit import TransmonQubit
    import numpy as np

    # Qubit at flux=0.15 — away from sweet spot, good sensitivity
    _q = TransmonQubit(
        EC=0.2 * 2 * np.pi, EJ=15 * 2 * np.pi,
        T1=10000, T2=5000, n_levels=2, flux=0.15,
    )
    t_rabi = np.linspace(0, 10, 8)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=_q.frequency)

    # Flux kernel (single-sided difference, legacy mode)
    ke_flux = KernelEstimator(mode='flux')
    t_flux, k_flux = ke_flux.estimate(pulse, _q)

    # Omega kernel (math VZ)
    ke_omega = KernelEstimator(
        mode='omega', virtual_z_impl='math', stim_amplitude=0.01,
    )
    t_omega, k_omega = ke_omega.estimate(pulse, _q)

    # dω/dΦ at operating point (GHz / Φ₀)
    kappa = abs(_q.frequency_sensitivity(_q.flux))
    assert kappa > 1.0, f"kappa too small ({kappa}) — flux bias needed"

    # Compare integrated kernels G = ∫k(t) dt — the quantity used in
    # Wiener deconvolution (see frequency.py:294).
    G_flux = np.trapezoid(k_flux, t_flux)
    G_omega = np.trapezoid(k_omega, t_omega)
    G_predicted = kappa * G_omega
    rel_diff = abs(G_flux - G_predicted) / max(abs(G_flux), abs(G_predicted))
    assert rel_diff < 0.15, (
        f"Integrated kernel mismatch: G_flux={G_flux:.6g}, "
        f"kappa*G_omega={G_predicted:.6g}, rel_diff={rel_diff:.3e}"
    )


def test_backward_compat_default(qubit):
    """KernelEstimator() (no args) produces same result as explicit mode='flux'."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 8)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    # Default estimator (mode='flux' internally)
    ke_default = KernelEstimator()
    t_def, k_def = ke_default.estimate(pulse, qubit)

    # Explicit flux mode
    ke_explicit = KernelEstimator(mode='flux')
    t_exp, k_exp = ke_explicit.estimate(pulse, qubit)

    assert_array_close(k_def, k_exp, name="default_vs_explicit_flux")
    assert_array_close(t_def, t_exp, name="t_default_vs_explicit")


def test_kernel_result_dataclass(qubit):
    """KernelResult fields are populated correctly and k1 property works."""
    from sqc.reconstruction.kernel import KernelEstimator, KernelResult
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 6)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    ke = KernelEstimator(mode='flux')
    res = ke.estimate_full(pulse, qubit)

    assert isinstance(res, KernelResult)
    assert isinstance(res.t_samples, np.ndarray)
    assert isinstance(res.kernels, list)
    assert len(res.kernels) == 1  # order == 1
    assert res.mode == 'flux'
    assert res.method == 'exp'
    assert res.order == 1
    assert res.stim_amplitude == 0.0215
    assert res.units == '1/(Φ₀·ns)'

    # k1 property
    k1 = res.k1
    assert isinstance(k1, np.ndarray)
    assert len(k1) == len(res.t_samples)
    assert_array_close(k1, res.kernels[0], name="k1_vs_kernels[0]")

    # Omega mode units
    ke_o = KernelEstimator(mode='omega')
    res_o = ke_o.estimate_full(pulse, qubit)
    assert res_o.mode == 'omega'
    assert res_o.units == 'rad⁻¹'


def test_vz_math_vs_hardware_consistency(qubit):
    """VZ hardware implementation produces non-trivial kernel, correlated with math.

    Hardware VZ phase-shifts sub-pulses whose trigger >= t_j.  The math VZ
    inserts a σ_z impulse.  These operations are NOT numerically identical —
    the hardware path applies the phase kick at sub-pulse granularity while
    the math path applies an instantaneous kick — but they should produce
    kernels of similar magnitude and shape for simple pulses.

    This is a sanity check, not a precision equivalence test.
    """
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    from sqc.control.pulse import CompositePulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 6)
    tau = 5.0  # ns free evolution
    pulse = create_ramsey_pulse(
        t_rabi, tau=tau, omega_d=qubit.frequency,
    )

    # Use full t_list (dense sampling) to capture kernel shape fully
    t_samples = pulse.t_list.copy()

    ke_math = KernelEstimator(
        mode='omega', virtual_z_impl='math', stim_amplitude=0.01,
    )
    ke_hw = KernelEstimator(
        mode='omega', virtual_z_impl='hardware', stim_amplitude=0.01,
    )

    t_m, k_m = ke_math.estimate(pulse, qubit, t_samples=t_samples)
    t_h, k_h = ke_hw.estimate(pulse, qubit, t_samples=t_samples)

    # Both kernels should be non-trivial
    assert np.max(np.abs(k_m)) > 0.0, "Math VZ kernel is all zeros"
    assert np.max(np.abs(k_h)) > 0.0, "Hardware VZ kernel is all zeros"

    # Integrated absolute kernels should be within 50%
    G_m = np.trapezoid(np.abs(k_m), t_m)
    G_h = np.trapezoid(np.abs(k_h), t_h)
    rel_diff = abs(G_m - G_h) / max(G_m, G_h)
    assert rel_diff < 0.50, (
        f"VZ math vs hardware magnitude mismatch: "
        f"∫|k_math|={G_m:.4f}, ∫|k_hw|={G_h:.4f}, rel_diff={rel_diff:.3e}"
    )

    # Signs should be predominantly consistent (dot product > 0)
    dot = np.dot(k_m, k_h)
    assert dot > 0, (
        f"VZ math and hardware kernels have opposite sign patterns: "
        f"dot(k_math, k_hw) = {dot:.4f}"
    )
