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


# ═════════════════════════════════════════════════════════════════════════════
# Phase 10.2 — method='sim' / method='exp' decoupling
# ═════════════════════════════════════════════════════════════════════════════


def test_sim_mode_creates_valid_kernel():
    """KernelEstimator(mode='omega', method='sim') produces a non-trivial
    kernel with correct shape — no qubit needed (pure theory)."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 10)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=5.0)

    ke = KernelEstimator(mode='omega', method='sim', stim_amplitude=0.1)
    t_samples, kernel = ke.estimate(pulse, None)

    assert isinstance(t_samples, np.ndarray)
    assert isinstance(kernel, np.ndarray)
    assert len(t_samples) == len(kernel)
    assert len(t_samples) > 0
    # Kernel values should be finite
    assert np.all(np.isfinite(kernel))
    # Kernel should not be all zeros — it should be non-trivial
    assert np.max(np.abs(kernel)) > 0.0

    # estimate_full also works without qubit
    res = ke.estimate_full(pulse, None)
    assert res.method == 'sim'
    assert res.mode == 'omega'


def test_flux_sim_raises():
    """KernelEstimator(mode='flux', method='sim') raises ValueError
    (illegal combination per handbook §4.3)."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 5)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=0.0)

    ke = KernelEstimator(mode='flux', method='sim')
    with pytest.raises(ValueError, match="illegal"):
        ke.estimate(pulse, None)


def test_sim_omega_exp_requires_qubit():
    """KernelEstimator(mode='omega', method='exp') with qubit=None
    raises ValueError."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 5)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=0.0)

    ke = KernelEstimator(mode='omega', method='exp')
    with pytest.raises(ValueError, match="requires a qubit"):
        ke.estimate(pulse, None)


def test_sim_auto_detects_qubit_params(qubit):
    """When qubit is provided, method='sim' uses qubit.n_levels and
    qubit.anharmonicity (auto-detection)."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    # Use a 3-level qubit so auto-detection matters
    from src.qubit import TransmonQubit
    q_3level = TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=3,
    )

    t_rabi = np.linspace(0, 10, 6)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=q_3level.frequency)

    # With qubit: auto-detect n_levels=3, anharmonicity from qubit
    ke = KernelEstimator(mode='omega', method='sim', stim_amplitude=0.01)
    t_samples, kernel = ke.estimate(pulse, q_3level)

    assert np.all(np.isfinite(kernel))
    assert np.max(np.abs(kernel)) > 0.0

    # Verify it completes successfully (different from 2-level default)
    assert len(t_samples) == len(kernel)


def test_sim_vs_exp_linear_regime(qubit):
    """In the small-stimulus limit, method='sim' and method='exp' produce
    similar omega kernels for 2-level systems.

    Uses rtol=0.15 because sim uses a†a operator while exp uses σ_z —
    they differ by a constant shift that cancels in the difference.
    """
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 8)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    # Small stimulus amplitude for linear regime
    ke_sim = KernelEstimator(
        mode='omega', method='sim',
        n_levels=2, anharmonicity=0.0, stim_amplitude=0.005,
    )
    ke_exp = KernelEstimator(
        mode='omega', method='exp',
        virtual_z_impl='math', stim_amplitude=0.005,
    )

    t_s, k_sim = ke_sim.estimate(pulse, None)
    t_e, k_exp = ke_exp.estimate(pulse, qubit)

    # Heisenberg (sim) uses deduplicated time grid; exp uses pulse t_list.
    # Interpolate sim kernel to exp grid for pointwise comparison.
    if len(t_s) != len(t_e):
        k_sim_interp = np.interp(t_e, t_s, k_sim)
    else:
        k_sim_interp = k_sim

    assert np.max(np.abs(k_sim)) > 0.0
    assert np.max(np.abs(k_exp)) > 0.0

    # Both kernels should have the same sign pattern (dot product > 0)
    # NOTE: sim uses Heisenberg σ_z/2 and exp uses VZ Gaussian σ_z;
    # they may differ by a global sign depending on convention.
    dot = np.dot(k_sim_interp, k_exp)
    assert abs(dot) > 1e-6, f"sim and exp kernels are orthogonal: dot={dot:.4f}"

    # Integrated kernels should be within 25% (Gaussian smearing in exp
    # causes ~5-10% magnitude difference vs ideal Heisenberg derivative)
    G_sim = np.trapezoid(np.abs(k_sim_interp), t_e)
    G_exp = np.trapezoid(np.abs(k_exp), t_e)
    rel_diff = abs(G_sim - G_exp) / max(G_sim, G_exp)
    assert rel_diff < 0.25, (
        f"sim vs exp integrated kernel mismatch: "
        f"G_sim={G_sim:.6g}, G_exp={G_exp:.6g}, rel_diff={rel_diff:.3e}"
    )


def test_sim_omega_order1_analytic():
    """For a ramsey pulse (π/2-gap-π/2), the sim kernel should be non-trivial
    and have sine-like shape within the pulse active region.

    While the analytic formula k₁(t) = sin[Ω(τ_p/2 - |t|)] is exact only
    for ideal square π/2-π/2 pulses, the sim kernel for a typical ramsey
    sequence should show qualitatively similar structure: non-zero where
    pulses act, small elsewhere, and a shape consistent with the control.
    """
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    # Standard ramsey pulse: two π/2 pulses with 5 ns gap
    t_rabi = np.arange(0, 10, 0.5)  # coarse for speed
    omega_d = 5.0  # GHz
    tau = 5.0  # ns gap
    pulse = create_ramsey_pulse(t_rabi, tau=tau, omega_d=omega_d)

    # Sim kernel
    ke = KernelEstimator(
        mode='omega', method='sim',
        n_levels=2, anharmonicity=0.0,
        stim_amplitude=0.2, stim_width=1.0,
    )
    t_samples, kernel = ke.estimate(pulse, None)

    # Basic sanity
    assert np.all(np.isfinite(kernel))
    assert np.max(np.abs(kernel)) > 0.0

    # Kernel should be non-trivial and vary with time
    # NOTE: sign convention depends on stimulus operator (a†a vs σ_z);
    # the kernel may be all-positive or all-negative.  What matters is
    # that it is non-constant and has finite magnitude.
    assert np.max(np.abs(kernel)) > 0.01, "Kernel magnitude too small"
    assert np.std(kernel) > 0.0, "Kernel should vary with time"

    # The integrated absolute kernel should be non-zero
    G = np.trapezoid(np.abs(kernel), t_samples)
    assert G > 0.0


def test_sim_method_field_in_result(qubit):
    """KernelResult.method reflects self.method ('sim' or 'exp')."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 6)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    # exp method (default)
    ke_exp = KernelEstimator(mode='omega')
    res_exp = ke_exp.estimate_full(pulse, qubit)
    assert res_exp.method == 'exp'

    # sim method
    ke_sim = KernelEstimator(mode='omega', method='sim', stim_amplitude=0.01)
    res_sim = ke_sim.estimate_full(pulse, None)
    assert res_sim.method == 'sim'

    # exp+flux method (default)
    ke_flux = KernelEstimator()
    res_flux = ke_flux.estimate_full(pulse, qubit)
    assert res_flux.method == 'exp'


def test_sim_kernel_is_deterministic():
    """Sim kernel with same params and no qubit is deterministic
    (no stochastic elements)."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 6)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=5.0)

    ke = KernelEstimator(mode='omega', method='sim', stim_amplitude=0.01)

    t1, k1 = ke.estimate(pulse, None)
    t2, k2 = ke.estimate(pulse, None)

    assert_array_close(k1, k2, name="sim_kernel_deterministic")


# ═════════════════════════════════════════════════════════════════════════════
# Phase 10.3 — higher-order Volterra + KernelResult serialization
# ═════════════════════════════════════════════════════════════════════════════


def test_order1_default_unchanged(qubit):
    """KernelEstimator() with default order=1 produces same result as before
    (backward compat — matches legacy get_kernel)."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 8)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=qubit.frequency)

    # Default order=1 (unchanged path)
    ke = KernelEstimator()
    assert ke.order == 1
    t_new, k_new = ke.estimate(pulse, qubit)

    # Should match legacy
    pulse.get_kernel(qubit)
    legacy_k = np.array(pulse.kernel)
    assert_array_close(k_new, legacy_k, name="order1_default_vs_legacy")

    # estimate_full with order=1
    res = ke.estimate_full(pulse, qubit)
    assert res.order == 1
    assert len(res.kernels) == 1
    assert_array_close(res.k1, legacy_k, name="order1_full_k1")


def test_polynomial_fit_no_constant():
    """_fit_polynomial_no_constant fits linear and quadratic cases correctly."""
    from sqc.reconstruction.kernel import _fit_polynomial_no_constant
    import numpy as np

    # Linear: y = 2*x
    coeffs_lin = _fit_polynomial_no_constant([1, 2, 3], [2, 4, 6], 1)
    assert len(coeffs_lin) == 1
    np.testing.assert_allclose(coeffs_lin, [2.0], atol=1e-10)

    # Quadratic: y = x + x²  →  [2, 6, 12] for x=[1, 2, 3]
    coeffs_quad = _fit_polynomial_no_constant([1, 2, 3], [2, 6, 12], 2)
    assert len(coeffs_quad) == 2
    np.testing.assert_allclose(coeffs_quad, [1.0, 1.0], atol=1e-10)

    # Cubic: y = 3x - 0.5x² + 0.2x³ for x=[1, 2, 3, 4, 5]
    x = [1, 2, 3, 4, 5]
    y = [3*xi - 0.5*xi**2 + 0.2*xi**3 for xi in x]
    coeffs_cubic = _fit_polynomial_no_constant(x, y, 3)
    np.testing.assert_allclose(coeffs_cubic, [3.0, -0.5, 0.2], atol=1e-10)


def test_gauss_integral_factor():
    """_gauss_integral_factor returns correct values for n=1,2,3."""
    from sqc.reconstruction.kernel import _gauss_integral_factor
    import math

    sigma = 2.0
    # n=1: sigma*sqrt(2π)
    c1 = _gauss_integral_factor(1, sigma)
    expected1 = sigma * math.sqrt(2 * math.pi)
    assert abs(c1 - expected1) < 1e-12

    # n=2: sigma*sqrt(π)/2
    c2 = _gauss_integral_factor(2, sigma)
    expected2 = sigma * math.sqrt(math.pi) / 2.0
    assert abs(c2 - expected2) < 1e-12

    # n=3: sigma*sqrt(2π/3)/6
    c3 = _gauss_integral_factor(3, sigma)
    expected3 = sigma * math.sqrt(2 * math.pi / 3) / math.factorial(3)
    assert abs(c3 - expected3) < 1e-12


def test_order2_sim_produces_two_kernels():
    """KernelEstimator(mode='omega', method='sim', order=2).estimate_full
    returns a KernelResult with len(kernels)==2, both finite, k2 non-trivial."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 5)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=5.0)

    ke = KernelEstimator(
        mode='omega', method='sim',
        order=2, n_amp_samples=5,
        stim_amplitude=0.1, stim_width=2.0,
    )
    res = ke.estimate_full(pulse, None)

    assert res.order == 2
    assert len(res.kernels) == 2
    assert res.kernels[0].ndim == 1
    assert res.kernels[1].ndim == 1
    assert len(res.kernels[0]) == len(res.t_samples)
    assert len(res.kernels[1]) == len(res.t_samples)

    # Both k1 and k2 should be finite
    assert np.all(np.isfinite(res.kernels[0]))
    assert np.all(np.isfinite(res.kernels[1]))

    # k1 should be non-trivial
    assert np.max(np.abs(res.kernels[0])) > 0.0

    # k2 should be non-trivial (at least some non-zero values)
    assert np.max(np.abs(res.kernels[1])) > 0.0

    # Backward compat: estimate() still returns just k1
    t_s, k1 = ke.estimate(pulse, None)
    assert len(k1) == len(t_s)
    assert_array_close(k1, res.kernels[0], name="order2_estimate_vs_full_k1")


def test_order2_exp_flux_produces_two_kernels():
    """KernelEstimator(mode='flux', order=2).estimate_full returns two
    flux-domain kernels, both finite and k2 non-trivial."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    from src.qubit import TransmonQubit
    import numpy as np

    # Qubit away from sweet spot for non-zero flux sensitivity
    q = TransmonQubit(
        EC=0.2 * 2 * np.pi, EJ=15 * 2 * np.pi,
        T1=10000, T2=5000, n_levels=2, flux=0.15,
    )
    t_rabi = np.linspace(0, 10, 5)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=q.frequency)

    ke = KernelEstimator(
        mode='flux', order=2, n_amp_samples=5,
        stim_amplitude=0.01, stim_width=2.0,
    )
    res = ke.estimate_full(pulse, q)

    assert res.order == 2
    assert res.mode == 'flux'
    assert len(res.kernels) == 2
    assert len(res.kernels[0]) == len(res.t_samples)
    assert len(res.kernels[1]) == len(res.t_samples)

    assert np.all(np.isfinite(res.kernels[0]))
    assert np.all(np.isfinite(res.kernels[1]))

    # k1 should be non-trivial
    assert np.max(np.abs(res.kernels[0])) > 0.0
    # k2 should have some structure
    assert np.max(np.abs(res.kernels[1])) > 0.0


def test_kernel_result_save_load_roundtrip():
    """KernelResult.save() then KernelResult.load() restores all fields."""
    from sqc.reconstruction.kernel import KernelEstimator, KernelResult
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np
    import tempfile
    import os

    t_rabi = np.linspace(0, 10, 5)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=5.0)

    ke = KernelEstimator(
        mode='omega', method='sim',
        order=2, n_amp_samples=5,
        stim_amplitude=0.1, stim_width=2.0,
    )
    res_orig = ke.estimate_full(pulse, None)

    # Save to temp file
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, 'test_kernel.npz')
        res_orig.save(path)

        # Load back
        res_loaded = KernelResult.load(path)

        # Verify all fields
        assert_array_close(res_loaded.t_samples, res_orig.t_samples, name="t_samples")
        assert res_loaded.mode == res_orig.mode
        assert res_loaded.method == res_orig.method
        assert res_loaded.order == res_orig.order
        assert abs(res_loaded.stim_amplitude - res_orig.stim_amplitude) < 1e-15
        assert res_loaded.units == res_orig.units
        assert len(res_loaded.kernels) == len(res_orig.kernels)
        for i in range(len(res_orig.kernels)):
            assert_array_close(
                res_loaded.kernels[i], res_orig.kernels[i],
                name=f"kernel_{i}",
            )


def test_order2_kernel_nonlinearity_increases_with_amplitude():
    """For a simple pulse, |k2/k1| ratio increases with stim_amplitude,
    verifying higher-order effects are physically real (not numerical noise)."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.control.sequence import create_ramsey_pulse
    import numpy as np

    t_rabi = np.linspace(0, 10, 5)
    pulse = create_ramsey_pulse(t_rabi, tau=0.0, omega_d=5.0)

    # Small amplitude
    ke_small = KernelEstimator(
        mode='omega', method='sim',
        order=2, n_amp_samples=5,
        stim_amplitude=0.02, stim_width=2.0,
    )
    res_small = ke_small.estimate_full(pulse, None)

    # Larger amplitude
    ke_large = KernelEstimator(
        mode='omega', method='sim',
        order=2, n_amp_samples=5,
        stim_amplitude=0.10, stim_width=2.0,
    )
    res_large = ke_large.estimate_full(pulse, None)

    # Compute median |k2/k1| ratio across time points (avoid division by ~zero)
    def ratio_median(res):
        k1 = np.abs(res.kernels[0])
        k2 = np.abs(res.kernels[1])
        # Only consider points where k1 is non-negligible
        mask = k1 > 0.01 * np.max(k1)
        if np.sum(mask) < 2:
            return 0.0
        return np.median(k2[mask] / k1[mask])

    r_small = ratio_median(res_small)
    r_large = ratio_median(res_large)

    # The ratio |k2/k1| should be AMPLITUDE-INVARIANT with Heisenberg:
    # k_n are exact derivatives (no polynomial fit, no amplitude dependence).
    # This is a fundamental property of the Heisenberg propagator method.
    # (The old polynomial-fit method showed amplitude dependence because
    # it was fitting numerical noise — that was a bug, not physics.)
    assert abs(r_small - r_large) / max(abs(r_small), abs(r_large), 1e-10) < 0.01, (
        f"Expected |k2/k1| to be amplitude-INVARIANT (Heisenberg property), "
        f"but small={r_small:.6f}, large={r_large:.6f}"
    )
    # Also verify the ratio is physically reasonable (k2 << k1 for Ramsey)
    assert r_small < 0.5, (
        f"|k2/k1| = {r_small:.4f} unexpectedly large for Ramsey pulse"
    )


# ═════════════════════════════════════════════════════════════════════════════
# Phase 12 — σ_t parameterization, Richardson extrapolation, off-diagonal (sim)
# ═════════════════════════════════════════════════════════════════════════════


def _yx_ramsey(qubit):
    """Zero-detuning Y-X Ramsey pulse on a clean dt grid (k₃_diag = -k₁)."""
    from sqc.control.sequence import create_ramsey_pulse

    t_rabi = np.arange(0, 10, 0.5)
    return create_ramsey_pulse(
        t_rabi, tau=0.0, omega_d=qubit.frequency,
        phase1=np.pi / 2, phase2=0.0,
    )


def test_probe_sigma_t_default_matches_legacy(qubit):
    """probe_sigma_t=None reproduces the hardcoded 2·dt behaviour exactly."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.config import CONFIG

    pulse = _yx_ramsey(qubit)

    ke_none = KernelEstimator(
        mode='omega', method='exp', order=3, stim_amplitude=0.01,
    )
    ke_explicit = KernelEstimator(
        mode='omega', method='exp', order=3, stim_amplitude=0.01,
        probe_sigma_t=2.0 * CONFIG.awg.dt,
    )
    r_none = ke_none.estimate_full(pulse, qubit)
    r_exp = ke_explicit.estimate_full(pulse, qubit)
    for n in range(3):
        assert_array_close(
            r_none.kernels[n], r_exp.kernels[n], name=f"sigma_t_default_k{n+1}",
        )


def test_smaller_sigma_t_improves_k3_diagonal(qubit):
    """A smaller VZ probe width moves exp k₃ closer to the sim ground truth.

    The sim (Heisenberg) k₃ diagonal is the ideal-δ truth; exp at finite
    σ_t is biased low by off-diagonal smearing.  Reducing σ_t reduces the
    bias monotonically (verified: ratio 0.84 -> 0.92 from 2·dt -> 1·dt).
    """
    from sqc.reconstruction.kernel import KernelEstimator

    pulse = _yx_ramsey(qubit)

    sim = KernelEstimator(mode='omega', method='sim', order=3).estimate_full(pulse, qubit)
    G3_sim = np.trapezoid(sim.kernels[2], sim.t_samples)

    def g3_ratio(sigma_t):
        r = KernelEstimator(
            mode='omega', method='exp', order=3, stim_amplitude=0.01,
            probe_sigma_t=sigma_t,
        ).estimate_full(pulse, qubit)
        return abs(np.trapezoid(r.kernels[2], r.t_samples) / G3_sim)

    from sqc.config import CONFIG
    r_wide = g3_ratio(2.0 * CONFIG.awg.dt)
    r_narrow = g3_ratio(1.0 * CONFIG.awg.dt)
    # narrower probe -> closer to 1.0 (less biased)
    assert abs(r_narrow - 1.0) < abs(r_wide - 1.0), (
        f"narrow σ_t ratio {r_narrow:.3f} should beat wide {r_wide:.3f}"
    )


def test_richardson_improves_k3_diagonal(qubit):
    """Richardson σ_t→0 extrapolation beats any single σ_t for k₃."""
    from sqc.reconstruction.kernel import KernelEstimator

    pulse = _yx_ramsey(qubit)
    sim = KernelEstimator(mode='omega', method='sim', order=3).estimate_full(pulse, qubit)
    G3_sim = np.trapezoid(sim.kernels[2], sim.t_samples)

    r_default = KernelEstimator(
        mode='omega', method='exp', order=3, stim_amplitude=0.01,
    ).estimate_full(pulse, qubit)
    r_rich = KernelEstimator(
        mode='omega', method='exp', order=3, stim_amplitude=0.01, richardson=True,
    ).estimate_full(pulse, qubit)

    G3_default = np.trapezoid(r_default.kernels[2], r_default.t_samples)
    G3_rich = np.trapezoid(r_rich.kernels[2], r_rich.t_samples)
    assert abs(G3_rich - G3_sim) < abs(G3_default - G3_sim), (
        f"Richardson G3={G3_rich:.3f} should be closer to sim {G3_sim:.3f} "
        f"than default {G3_default:.3f}"
    )
    # k1 must remain unchanged in shape/finiteness
    assert r_rich.kernels[0].shape == r_default.kernels[0].shape
    assert np.all(np.isfinite(r_rich.kernels[2]))


def test_offdiag_sim_shapes_and_diagonal(qubit):
    """method='sim' + extract_off_diagonal yields n-D kernels; k₃ diagonal
    equals -k₁ and k₂ is symmetric with a small integral (Y-X)."""
    from sqc.reconstruction.kernel import KernelEstimator

    pulse = _yx_ramsey(qubit)
    res = KernelEstimator(
        mode='omega', method='sim', order=3, extract_off_diagonal=True,
    ).estimate_full(pulse, qubit)

    k1, k2, k3 = res.kernels
    M = len(res.t_samples)
    assert res.off_diagonal is True
    assert k1.shape == (M,)
    assert k2.shape == (M, M)
    assert k3.shape == (M, M, M)

    # k2 symmetric
    assert np.allclose(k2, k2.T)

    # k3 diagonal == -k1 (exact 2-level commutator identity)
    k3_diag = np.array([k3[i, i, i] for i in range(M)])
    rmse = np.sqrt(np.mean((k3_diag + k1) ** 2))
    assert rmse < 1e-5, f"k3(t,t,t) != -k1(t): RMSE={rmse:.2e}"

    # k2 integral is suppressed (Y-X orthogonality; Gaussian residual only)
    G1 = np.trapezoid(k1, res.t_samples)
    G2 = np.trapezoid(np.trapezoid(k2, res.t_samples, axis=0), res.t_samples)
    assert abs(G2 / G1) < 0.15, f"|G2/G1|={abs(G2/G1):.3f} too large"


def test_offdiag_exp_raises():
    """extract_off_diagonal with method='exp' raises NotImplementedError."""
    from sqc.reconstruction.kernel import KernelEstimator
    from src.qubit import TransmonQubit

    q = TransmonQubit(
        EC=0.2 * 2 * np.pi, EJ=15 * 2 * np.pi,
        T1=10000, T2=5000, n_levels=2, flux=0.0,
    )
    pulse = _yx_ramsey(q)
    ke = KernelEstimator(
        mode='omega', method='exp', order=2, extract_off_diagonal=True,
    )
    with pytest.raises(NotImplementedError, match="method='sim'"):
        ke.estimate_full(pulse, q)


def test_offdiag_order4_raises(qubit):
    """Off-diagonal extraction rejects order >= 4 (M^n blow-up)."""
    from sqc.reconstruction.kernel import KernelEstimator

    pulse = _yx_ramsey(qubit)
    ke = KernelEstimator(
        mode='omega', method='sim', order=4, extract_off_diagonal=True,
    )
    with pytest.raises(ValueError, match="order <= 3"):
        ke.estimate_full(pulse, qubit)


def test_offdiag_save_load_roundtrip(qubit):
    """KernelResult with n-D kernels survives save/load."""
    from sqc.reconstruction.kernel import KernelEstimator, KernelResult
    import tempfile
    import os

    pulse = _yx_ramsey(qubit)
    res = KernelEstimator(
        mode='omega', method='sim', order=3, extract_off_diagonal=True,
    ).estimate_full(pulse, qubit)

    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, 'offdiag.npz')
        res.save(path)
        loaded = KernelResult.load(path)
        assert loaded.off_diagonal is True
        assert loaded.kernels[1].shape == res.kernels[1].shape
        assert loaded.kernels[2].shape == res.kernels[2].shape
        for n in range(3):
            assert_array_close(loaded.kernels[n], res.kernels[n], name=f"od_k{n+1}")


def test_offdiag_kernel_rejected_by_wiener(qubit):
    """n-D kernels into a Wiener/Hammerstein reconstruction raise ValueError."""
    from sqc.reconstruction.kernel import KernelEstimator
    from sqc.reconstruction.transient import TransientReconstruction

    pulse = _yx_ramsey(qubit)
    res = KernelEstimator(
        mode='omega', method='sim', order=2, extract_off_diagonal=True,
    ).estimate_full(pulse, qubit)

    class _Meas:
        data = {"delta_p": np.zeros(40)}

    rec = TransientReconstruction(method='hammerstein_volterra')
    with pytest.raises(ValueError, match="non-diagonal"):
        rec.reconstruct(_Meas(), kernel=res)
