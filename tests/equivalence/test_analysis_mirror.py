"""Equivalence tests: src_mirror.analysis.Analysis == src.analysis.Analysis.

For each analysis method that has been ported to the mirror layer,
verify that the old and new Analysis classes produce identical
numeric results within the physics regression tolerance
(rtol=1e-6, atol=1e-9).

Methods tested:
  - get_signal_from_ramsey_by_unwrap
  - get_signal_from_diff_echo
  - wiener_deconvolution
  - hammerstein_wiener_deconvolution
  - get_expectation_values
  - get_population

Note: get_signal_from_ramsey_by_iq is NOT tested for equivalence
because both old and new call plt.show() which triggers Qt crash.
The internal logic (arctan2 + unwrap + gradient) is tested via
RamseyUnwrap equivalence instead.
"""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import assert_array_close


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_qubit(n_levels=2, flux=0.01):
    """Create a legacy TransmonQubit with test parameters."""
    from src.qubit import TransmonQubit

    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=flux,
        state=0,
        n_levels=n_levels,
    )


# ---------------------------------------------------------------------------
# Ramsey unwrap
# ---------------------------------------------------------------------------

def test_analysis_mirror_ramsey_unwrap():
    """RamseyUnwrap: old vs new produce identical B(tau)."""
    from src.analysis import Analysis as OldAnalysis
    from src_mirror.analysis import Analysis as NewAnalysis

    q_old = _make_qubit()
    q_new = _make_qubit()

    tau = np.linspace(0, 100, 100)
    phi_true = 0.2 * np.sin(2 * np.pi * tau / 100)
    p_e = 0.5 * (1 - np.cos(phi_true))

    old_ana = OldAnalysis()
    B_old = old_ana.get_signal_from_ramsey_by_unwrap(q_old, tau, p_e, k_span=3)

    new_ana = NewAnalysis()
    B_new = new_ana.get_signal_from_ramsey_by_unwrap(q_new, tau, p_e, k_span=3)

    assert_array_close(B_old, B_new, name="ramsey_unwrap_B")


# ---------------------------------------------------------------------------
# Differential echo
# ---------------------------------------------------------------------------

def test_analysis_mirror_diff_echo():
    """DiffEcho: old vs new produce identical B."""
    from src.analysis import Analysis as OldAnalysis
    from src_mirror.analysis import Analysis as NewAnalysis

    q_old = _make_qubit()
    q_new = _make_qubit()

    p_e = 0.5 + 0.1 * np.sin(np.linspace(0, np.pi, 50))
    t_int = 10.0
    k = 2

    old_ana = OldAnalysis()
    B_old = old_ana.get_signal_from_diff_echo(q_old, p_e, t_int, k)

    new_ana = NewAnalysis()
    B_new = new_ana.get_signal_from_diff_echo(q_new, p_e, t_int, k)

    assert_array_close(B_old, B_new, name="diff_echo_B")


# ---------------------------------------------------------------------------
# Wiener deconvolution
# ---------------------------------------------------------------------------

def test_analysis_mirror_wiener():
    """Wiener deconvolution: old vs new produce identical (t_list, signal)."""
    from src.analysis import Analysis as OldAnalysis
    from src_mirror.analysis import Analysis as NewAnalysis

    np.random.seed(42)
    delta_p = np.random.randn(100) * 0.01
    kernel = np.exp(-np.linspace(0, 5, 30) ** 2)
    kernel = kernel / np.abs(kernel).max()
    dt = 1.0
    lambdas = 10.0

    old_ana = OldAnalysis()
    t_old, B_old = old_ana.wiener_deconvolution(delta_p, kernel, dt, lambdas)

    new_ana = NewAnalysis()
    t_new, B_new = new_ana.wiener_deconvolution(delta_p, kernel, dt, lambdas)

    assert_array_close(t_old, t_new, name="wiener_t_list")
    assert_array_close(B_old, B_new, name="wiener_signal")


# ---------------------------------------------------------------------------
# Hammerstein-Wiener deconvolution
# ---------------------------------------------------------------------------

@pytest.mark.xfail(
    reason=(
        "src/analysis.py hammerstein_wiener_deconvolution has a pre-existing "
        "NaN bug (arccos domain error when argument ∉ [−1,1]). "
        "sqc/reconstruction/transient.py correctly clips with np.clip(…, 0, 1). "
        "src/ frozen by R1 — divergence is intentional."
    ),
    strict=True,
)
def test_analysis_mirror_hammerstein():
    """HammersteinWiener: old vs new produce identical (B_list, B).

    XFAIL: src/ has a NaN bug in arccos that sqc/ correctly fixes.
    """
    from src.analysis import Analysis as OldAnalysis
    from src_mirror.analysis import Analysis as NewAnalysis

    np.random.seed(43)
    q_old = _make_qubit()
    q_new = _make_qubit()

    delta_p = np.random.randn(100) * 0.01
    kernel = np.exp(-np.linspace(0, 5, 30) ** 2)
    kernel = kernel / np.abs(kernel).max()
    dt = 1.0
    lambdas = 10.0

    old_ana = OldAnalysis()
    B_list_old, B_old = old_ana.hammerstein_wiener_deconvolution(
        q_old, delta_p, kernel, dt, lambdas,
    )

    new_ana = NewAnalysis()
    B_list_new, B_new = new_ana.hammerstein_wiener_deconvolution(
        q_new, delta_p, kernel, dt, lambdas,
    )

    assert_array_close(B_list_old, B_list_new, name="hw_B_list")
    assert_array_close(B_old, B_new, name="hw_B")


# ---------------------------------------------------------------------------
# Expectation values and population extraction
# ---------------------------------------------------------------------------

def test_analysis_mirror_expectation_and_population():
    """Expectation and population extraction produce identical results.

    Note: src/analysis.py:Analysis.get_population has a pre-existing bug
    (checks result.state instead of result.states), so it always returns
    None. We verify the new facade matches the correct sqc functions
    directly, and that get_expectation_values matches the old code
    (which works correctly).
    """
    from qutip import basis, destroy, mesolve

    from src.analysis import Analysis as OldAnalysis
    from src_mirror.analysis import Analysis as NewAnalysis
    from sqc.simulation.result import (
        extract_expectation,
        extract_population,
    )

    a = destroy(2)
    n_op = a.dag() * a
    H = 0.5 * n_op
    psi0 = basis(2, 0)
    t_list = np.linspace(0, 10, 100)
    e_ops = [n_op]

    result = mesolve(H, psi0, t_list, c_ops=[], e_ops=e_ops,
                     options={"store_states": True})

    # Expectation: old and new should match (this path works in old code)
    old_ana = OldAnalysis()
    new_ana = NewAnalysis()

    expect_old = old_ana.get_expectation_values(result, 0)
    expect_new = new_ana.get_expectation_values(result, 0)
    assert_array_close(expect_old, expect_new, name="expectation")
    # Both should match the direct function
    expect_direct = extract_expectation(result, 0)
    assert_array_close(expect_new, expect_direct, name="expectation_vs_direct")

    # Population: new facade should match direct sqc function
    # (old code has bug: checks result.state instead of result.states
    #  and returns None; not testable for equivalence)
    pop_new = new_ana.get_population(result, 0)
    pop_direct = extract_population(result, 0)
    assert_array_close(pop_new, pop_direct, name="population_vs_direct")


# ---------------------------------------------------------------------------
# Not-implemented methods
# ---------------------------------------------------------------------------

def test_analysis_mirror_not_implemented():
    """Methods dependent on Track B raise NotImplementedError.

    numerical_inverse implemented in P3b; get_h_from_phi in P3c.
    """
    from src_mirror.analysis import Analysis

    ana = Analysis()

    for method_name, args in [
        ("get_signal_from_cryoscope", (None, None, None, None, None, None)),
        ("get_volterra_kernel", (None, None)),
    ]:
        with pytest.raises(NotImplementedError):
            getattr(ana, method_name)(*args)


def test_analysis_mirror_get_h_from_phi():
    """get_h_from_phi implemented in P3c — pure interpolation, no Track B req."""
    import numpy as np
    from src_mirror.analysis import Analysis
    from src.analysis import Analysis as OldAnalysis

    h_list = np.linspace(-0.03, 0.03, 21)
    phi_list = 3.0 * h_list + 0.1 * h_list**3  # nonlinear monotonic

    ana_new = Analysis()
    ana_old = OldAnalysis()

    phi_new, h_new = ana_new.get_h_from_phi(h_list, phi_list)
    phi_old, h_old = ana_old.get_h_from_phi(h_list, phi_list)

    # Test forward: phi(h) should be close to original
    h_test = np.linspace(-0.02, 0.02, 11)
    from tests.conftest import assert_array_close
    assert_array_close(phi_new(h_test), phi_old(h_test), name="get_h_from_phi_forward")
    # Test inverse: h(phi) roundtrip
    phi_test = np.linspace(-0.06, 0.06, 11)
    assert_array_close(h_new(phi_test), h_old(phi_test), name="get_h_from_phi_inverse")


# ---------------------------------------------------------------------------
# Module-level stubs
# ---------------------------------------------------------------------------

def test_analysis_mirror_module_stubs():
    """Module-level LM helpers are no longer stubs (P3b implemented).

    Verify they are callable (need actual arguments — just check
    they are defined and not raising NotImplementedError on import).
    """
    from src_mirror.analysis import (
        forward_simulation,
        compute_jacobian,
        compute_jacobian_finite_difference,
        levenberg_marquardt,
    )

    for fn in [
        forward_simulation,
        compute_jacobian,
        compute_jacobian_finite_difference,
        levenberg_marquardt,
    ]:
        assert callable(fn), f"{fn.__name__} is not callable"


# ---------------------------------------------------------------------------
# Module-level re-exports
# ---------------------------------------------------------------------------

def test_analysis_mirror_module_exports():
    """Module-level helpers are correctly re-exported."""
    from src_mirror.analysis import (
        generate_basis_functions,
        basis_function_decomposition,
        R,
    )
    from sqc.reconstruction.basis import (
        generate_basis_functions as gbf_orig,
        basis_function_decomposition as bfd_orig,
        R as R_orig,
    )

    assert generate_basis_functions is gbf_orig
    assert basis_function_decomposition is bfd_orig
    assert R is R_orig
