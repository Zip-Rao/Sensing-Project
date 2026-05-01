"""Generate physics regression baselines.

Run once, before any refactor work begins:

    python -m tests.regression.generate_baselines

Re-run only when intentionally changing physics (require justification
in commit message).
"""
from __future__ import annotations

import sys
from pathlib import Path

# allow standalone execution
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
from qutip import basis

from src.qubit import TransmonQubit
from src.signal import Signal, CompositeSignal
from src.pulse import (
    create_ramsey_pulse,
    create_diff_echo_pulse,
)
from src.protocal import Protocal
from tests.conftest import save_baseline


def _make_default_qubit() -> TransmonQubit:
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=2,
    )


def _baseline_qubit_static() -> dict:
    """Frequency, anharmonicity, sensitivity at default flux."""
    q = _make_default_qubit()
    return {
        "frequency": float(q.frequency),
        "anharmonicity": float(q.anharmonicity),
        "sensitivity_at_zero": float(q.frequency_sensitivity(0.0)),
        "sensitivity_at_optimal": float(
            q.frequency_sensitivity(q.optimal_work_point() / np.pi)
        ),
    }


def _baseline_ramsey() -> dict:
    """Run Protocal.evolve case 1 and snapshot outputs."""
    q = _make_default_qubit()
    proto = Protocal(type=1)
    proto.initialize(q, state=0)

    # NOTE: Protocal.evolve case 1 returns Phi, tau_list, p_e_list
    Phi, tau_list, p_e_list = proto.evolve(q)
    return {
        "Phi_signal": np.asarray(Phi.signal),
        "Phi_t_list": np.asarray(Phi.t_list),
        "tau_list": np.asarray(tau_list),
        "p_e_list": np.asarray(p_e_list),
    }


def _baseline_diff_echo() -> dict:
    """Run Protocal.evolve case 2."""
    q = _make_default_qubit()
    proto = Protocal(type=2)
    proto.initialize(q, state=0)
    Phi, tau_list, p_e_list, k, t_int = proto.evolve(q)
    return {
        "Phi_signal": np.asarray(Phi.signal),
        "Phi_t_list": np.asarray(Phi.t_list),
        "tau_list": np.asarray(tau_list),
        "p_e_list": np.asarray(p_e_list),
        "k": int(k),
        "t_int": float(t_int),
    }


def _baseline_transient() -> dict:
    """Run Protocal.evolve case 4 (sliding measurement)."""
    q = _make_default_qubit()
    proto = Protocal(type=4)
    proto.initialize(q, state=0)
    t_samples, kernel, scan_list, delta_p, p_e, Phi, ctrl = proto.evolve(q)
    return {
        "t_samples": np.asarray(t_samples),
        "kernel": np.asarray(kernel),
        "scan_list": np.asarray(scan_list),
        "delta_p": np.asarray(delta_p),
        "p_e": np.asarray(p_e),
        "Phi_signal": np.asarray(Phi.signal),
        "Phi_t_list": np.asarray(Phi.t_list),
    }


def _baseline_lm() -> dict:
    """Generate LM reconstruction baseline.

    Runs the LM numerical inversion with minimal parameters (n_basis=5,
    max_iter=2) on a short time-domain signal. The derived p_meas and
    the LM optimal coefficients b_opt are saved as the baseline.

    Fast by design: ~2-3 seconds total runtime.
    """
    import numpy as np
    from qutip import qeye

    from src.qubit import TransmonQubit
    from src.signal import Signal
    from src.pulse import create_ramsey_pulse
    from src.analysis import (
        forward_simulation,
        basis_function_decomposition,
        levenberg_marquardt as lm_old,
    )

    np.random.seed(42)

    # ---- qubit ---------------------------------------------------------
    q = TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=2,
    )

    # ---- signal (known truth) ------------------------------------------
    t_list = np.linspace(0, 20, 11)
    B_true = Signal(type=2, t_list=t_list, amplitude=0.01, frequency=0.1)
    q.qubit_in_mag(B_true, frame=0, omega_d=q.frequency)

    # ---- control pulse -------------------------------------------------
    cp = create_ramsey_pulse(
        t_rabi=np.linspace(0, 5, 6), tau=5.0, omega_d=q.frequency,
    )

    # ---- measurement time axis -----------------------------------------
    meas_start = t_list[0] - 0.5 * cp.t_list[-1]
    meas_end = t_list[-1] + 0.5 * cp.t_list[-1]
    t_meas = np.linspace(
        meas_start, meas_end, len(t_list) + len(cp.t_list) - 1,
    )

    # ---- build H_list (replicates inner H function in LM) --------------
    n_meas = len(t_meas)
    H_list: list = [[] for _ in range(n_meas)]
    t_evolve_list: list = [[] for _ in range(n_meas)]

    for i, t_delay in enumerate(t_meas):
        delta = t_delay - 0.5 * cp.t_list[-1]
        t_start = min(t_list[0], delta)
        t_end = max(t_list[-1], delta + cp.t_list[-1])
        N_e = len(t_list) + len(cp.t_list) - 1
        t_evolve = np.linspace(t_start, t_end, N_e)
        t_evolve_list[i] = t_evolve

        freq_coeffs = np.zeros(N_e, dtype=float)
        for j, t in enumerate(t_evolve):
            if t_list[0] <= t <= t_list[-1]:
                index = np.searchsorted(t_list, t)
                index = min(index, len(t_list) - 1)
                freq_coeffs[j] = (
                    q.freq_coeffs[index]
                    if cp.frame == 0
                    else q.freq_coeffs[index] - cp.omega_d
                )
            else:
                freq_coeffs[j] = (
                    q.frequency
                    if cp.frame == 0
                    else q.frequency - cp.omega_d
                )

        id_coeffs = np.ones(N_e, dtype=complex)
        if cp.frame == 0:
            H_list[i].append(
                [
                    q.anharmonicity * q.n * (q.n - 1) * 0.5,
                    id_coeffs,
                ]
            )
            H_list[i].append(
                [q.n + 0.5 * qeye(q.n_levels), freq_coeffs]
            )
        else:
            H_list[i].append(
                [
                    q.anharmonicity * q.n * (q.n - 1) * 0.5,
                    id_coeffs,
                ]
            )
            H_list[i].append([q.n, freq_coeffs])

        for op, coeffs in cp.hamiltonian:
            coeff_global = np.zeros(N_e, dtype=complex)
            for j, t in enumerate(t_evolve):
                t_loc = t - delta
                if 0 <= t_loc <= cp.t_list[-1]:
                    index = np.searchsorted(cp.t_list, t_loc)
                    index = min(index, len(cp.t_list) - 1)
                    coeff_global[j] = coeffs[index]
            H_list[i].append([op, coeff_global])

    # ---- forward simulation to get p_meas -------------------------------
    results = forward_simulation(
        q, cp, B_true, t_meas, H_list, t_evolve_list,
    )
    p_meas = np.array([res.expect[0][-1] for res in results])

    # ---- LM reconstruction (minimal params for speed) --------------------
    n_basis = 5
    basis_type = "fourier"
    lambdas = 100.0
    max_iter = 2
    tol = 1e-3

    B_init = Signal(
        type=6, t_list=t_list, n_basis=n_basis, basis_type=basis_type,
    )
    B_guess = np.zeros_like(t_list)
    b_init = basis_function_decomposition(
        B_guess, t_list, B_init.basis_functions,
    )
    B_init.update_signal(b=b_init)

    b_opt, history = lm_old(
        q, p_meas, t_list, cp, b_init, B_init, lambdas, max_iter, tol,
    )

    return {
        "t_list": np.asarray(t_list),
        "t_meas": np.asarray(t_meas),
        "p_meas": np.asarray(p_meas),
        "b_init": np.asarray(b_init),
        "b_opt": np.asarray(b_opt),
        "n_basis": n_basis,
        "basis_type": basis_type,
        "B_true_signal": np.asarray(B_true.signal),
        "history_res_final": np.asarray(history["res"][-1]),
        "history_mu": np.asarray(history["mu"]),
    }


BASELINES = {
    "qubit_static": _baseline_qubit_static,
    "ramsey_default": _baseline_ramsey,
    "diff_echo_default": _baseline_diff_echo,
    "transient_default": _baseline_transient,
    "lm_default": _baseline_lm,
}


def main() -> None:
    print("Generating physics regression baselines ...")
    for name, fn in BASELINES.items():
        print(f"  * {name} ... ", end="", flush=True)
        data = fn()
        path = save_baseline(name, data)
        print(f"saved to {path.relative_to(path.parent.parent.parent)}")
    print("Done.")


if __name__ == "__main__":
    main()
