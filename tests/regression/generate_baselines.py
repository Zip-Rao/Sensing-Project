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


BASELINES = {
    "qubit_static": _baseline_qubit_static,
    "ramsey_default": _baseline_ramsey,
    "diff_echo_default": _baseline_diff_echo,
    "transient_default": _baseline_transient,
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
