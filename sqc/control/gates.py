"""sqc.control.gates — Two-qubit gate functions.

Verbatim port of src/qubit.py module-level gate functions:
    ideal_iSWAP, simulate_iSWAP, ideal_CZ, simulate_CZ.
"""
from __future__ import annotations

import numpy as np
from qutip import Qobj, destroy, propagator, qeye, tensor


def ideal_iSWAP() -> Qobj:
    """Ideal iSWAP gate matrix.

    Returns
    -------
    Qobj
        4×4 unitary in |00>,|01>,|10>,|11> basis.
    """
    iSWAP = Qobj(
        [
            [1, 0, 0, 0],
            [0, 0, 1j, 0],
            [0, 1j, 0, 0],
            [0, 0, 0, 1],
        ],
        dims=[[2, 2], [2, 2]],
    )
    return iSWAP


def simulate_iSWAP(qubit1, qubit2, g: float):
    """Simulate iSWAP gate (simple model, no coupler).

    Parameters
    ----------
    qubit1, qubit2 : TransmonQubit
    g : float
        Exchange coupling (GHz).

    Returns
    -------
    U_eff : Qobj
        4×4 effective unitary.
    leakage : float
    F : float
        Average gate fidelity to ideal iSWAP.
    """
    H_0 = tensor(
        qubit1.get_hamiltonian_rwa(qubit1.frequency),
        qeye(qubit2.n_levels),
    ) + tensor(
        qeye(qubit1.n_levels),
        qubit2.get_hamiltonian_rwa(qubit2.frequency),
    )
    H_int = g * (
        tensor(
            destroy(qubit1.n_levels), destroy(qubit2.n_levels).dag()
        )
        + tensor(
            destroy(qubit1.n_levels).dag(), destroy(qubit2.n_levels)
        )
    )
    H = H_0 + H_int
    T = np.pi / (2 * g)
    U = propagator(H, T)
    idx = [
        0, 1, qubit2.n_levels, qubit2.n_levels + 1
    ]
    U_eff = Qobj(
        U.full()[np.ix_(idx, idx)], dims=[[2, 2], [2, 2]]
    )
    U_eff = U_eff.dag()  # TODO: verify convention
    U_ideal = ideal_iSWAP()
    leakage = 1.0 - np.mean(np.sum(np.abs(U_eff.full()) ** 2, axis=0))
    F = (abs((U_ideal.dag() @ U_eff).tr()) ** 2 + 4) / (4 * 5)
    return U_eff, leakage, F


def ideal_CZ() -> Qobj:
    """Ideal CZ gate matrix.

    Returns
    -------
    Qobj
        4×4 unitary.
    """
    CZ = Qobj(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, -1],
        ],
        dims=[[2, 2], [2, 2]],
    )
    return CZ


def simulate_CZ(qubit1, qubit2, g: float):
    """Simulate CZ gate via flux pulsing.

    Parameters
    ----------
    qubit1, qubit2 : TransmonQubit
    g : float
        Exchange coupling (GHz).

    Returns
    -------
    U_eff_corr : Qobj
        4×4 effective unitary with phase correction.
    leakage : float
    F : float
    """
    H_0 = tensor(
        qubit1.get_hamiltonian_rwa(qubit1.frequency),
        qeye(qubit2.n_levels),
    ) + tensor(
        qeye(qubit1.n_levels),
        qubit2.get_hamiltonian_rwa(qubit2.frequency),
    )
    H_int = g * (
        tensor(
            destroy(qubit1.n_levels), destroy(qubit2.n_levels).dag()
        )
        + tensor(
            destroy(qubit1.n_levels).dag(), destroy(qubit2.n_levels)
        )
    )
    H_1 = H_0 + H_int

    def Delta_1(t, args):
        T, T1, T2 = args["T"], args["T1"], args["T2"]
        A = args["A"]
        if t < 0 or t > T:
            return 0.0
        elif t < T1:
            return A * np.sin(0.5 * np.pi * t / T1) ** 2
        elif t < T - T2:
            return A
        else:
            return A * np.sin(0.5 * np.pi * (T - t) / T2) ** 2

    T = np.pi / (np.sqrt(2.0) * g)
    T1 = T / 4
    T2 = T / 4
    A = -qubit2.anharmonicity
    args = {"T": T, "T1": T1, "T2": T2, "A": A}
    n_1 = tensor(qubit1.n, qeye(qubit2.n_levels))
    H = [H_1, [n_1, Delta_1]]
    U = propagator(H, T, args=args)
    idx = [
        0, 1, qubit2.n_levels, qubit2.n_levels + 1
    ]
    U_eff = Qobj(
        U.full()[np.ix_(idx, idx)], dims=[[2, 2], [2, 2]]
    )
    phi00 = np.angle(U_eff[0, 0])
    phi10 = np.angle(U_eff[1, 1])
    phi01 = np.angle(U_eff[2, 2])
    phi1 = phi10 - phi00
    phi2 = phi01 - phi00
    D = Qobj(
        np.diag(
            [
                np.exp(-1j * phi00),
                np.exp(-1j * phi10),
                np.exp(-1j * phi01),
                np.exp(-1j * (phi10 + phi01 - phi00)),
            ]
        ),
        dims=[[2, 2], [2, 2]],
    )
    U_eff_corr = D @ U_eff
    U_ideal = ideal_CZ()
    leakage = 1.0 - 0.25 * (np.sum(np.abs(U_eff_corr.full()) ** 2))
    F = (abs((U_ideal.dag() @ U_eff_corr).tr()) ** 2 + 4) / (4 * 5)
    return U_eff_corr, leakage, F
