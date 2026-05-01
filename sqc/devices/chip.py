"""sqc.devices.chip — CoupledSystem: two qubits + coupler + multi-mode cavity.

Verbatim port of src/qubit.py:Coupled_System → CoupledSystem.
Renamed: spelling corrected per _refactor_plan.md §7.1.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from qutip import Qobj, basis, propagator, sesolve, tensor

from .base import Device


class CoupledSystem(Device):
    """Two Transmon qubits coupled via tunable coupler to multi-mode resonator.

    Parameters
    ----------
    qubit1 : TransmonQubit
        First qubit.
    qubit2 : TransmonQubit
        Second qubit.
    cavity : Resonator
        Multi-mode resonator.
    """

    def __init__(
        self, qubit1, qubit2, cavity
    ) -> None:
        self.qubit1 = qubit1
        self.qubit2 = qubit2
        self.cavity = cavity

        # Tensor-product operators
        self.a1: Qobj = tensor(
            self.qubit1.a, self.qubit2.I, self.cavity.I
        )
        self.a2: Qobj = tensor(
            self.qubit1.I, self.qubit2.a, self.cavity.I
        )
        self.acav: list[Qobj] = [
            tensor(self.qubit1.I, self.qubit2.I, self.cavity.a[m])
            for m in range(self.cavity.M)
        ]
        self.g: list[list[float]] = self.initialize_g()

        # Build system Hamiltonian
        self.build_system()

    # -- Device ABC ---------------------------------------------------------

    @property
    def name(self) -> str:
        return "CoupledSystem"

    def hilbert_dim(self) -> int:
        return int(np.prod([
            self.qubit1.n_levels,
            self.qubit2.n_levels,
            *self.cavity.n_levels,
        ]))

    def hamiltonian_static(self) -> Qobj:
        return self.H

    def collapse_operators(self) -> list[Qobj]:
        # Not implemented for coupled system
        return []

    # -- Coupling initialization --------------------------------------------

    def initialize_g(self) -> list[list[float]]:
        """Initialize default coupling strengths.

        Returns
        -------
        list[list[float]]
            g[0] = qubit1-mode couplings, g[1] = qubit2-mode couplings.
        """
        g = [[], []]
        for _m in range(self.cavity.M):
            g[0].append(2 * np.pi * 0.01)
            g[1].append(2 * np.pi * 0.01)
        return g

    def control_g(
        self, g, is_time_dependent: bool = False, f=None
    ) -> None:
        """Control qubit-qubit exchange coupling strength.

        Parameters
        ----------
        g : list[list[float]]
            2 × M coupling strengths.
        is_time_dependent : bool
            Whether coupling is time-dependent.
        f : callable or None
            Time-dependent coupling function f(t, args) → list[list[float]].
        """
        self.g = g
        if not is_time_dependent:
            H_qc1 = sum(
                self.g[0][m]
                * (self.a1.dag() * self.acav[m] + self.a1 * self.acav[m].dag())
                for m in range(self.cavity.M)
            )
            H_qc2 = sum(
                self.g[1][m]
                * (self.a2.dag() * self.acav[m] + self.a2 * self.acav[m].dag())
                for m in range(self.cavity.M)
            )
            self.H = self.H_0 + H_qc1 + H_qc2
        else:
            self.H = [self.H_0]
            for M_idx in range(self.cavity.M):
                H_qc1_m = self.g[0][M_idx] * (
                    self.a1.dag() * self.acav[M_idx]
                    + self.a1 * self.acav[M_idx].dag()
                )
                H_qc2_m = self.g[1][M_idx] * (
                    self.a2.dag() * self.acav[M_idx]
                    + self.a2 * self.acav[M_idx].dag()
                )

                def f1(t, args, idx=M_idx, f_ref=f):
                    return f_ref(t, args)[0][idx]

                def f2(t, args, idx=M_idx, f_ref=f):
                    return f_ref(t, args)[1][idx]

                self.H.append([H_qc1_m, f1])
                self.H.append([H_qc2_m, f2])

    def build_system(
        self, is_rwa: bool = False, omega_d: Optional[float] = None
    ) -> None:
        """Build full system Hamiltonian.

        Parameters
        ----------
        is_rwa : bool
            Whether to use rotating wave approximation.
        omega_d : float or None
            Drive reference frequency.
        """
        if is_rwa:
            self.H_q1: Qobj = tensor(
                self.qubit1.get_hamiltonian_rwa(omega_d),
                self.qubit2.I,
                self.cavity.I,
            )
            self.H_q2: Qobj = tensor(
                self.qubit1.I,
                self.qubit2.get_hamiltonian_rwa(omega_d),
                self.cavity.I,
            )
            self.H_cav: Qobj = tensor(
                self.qubit1.I,
                self.qubit2.I,
                self.cavity.get_hamiltonian_rot(omega_d),
            )
        else:
            self.H_q1 = tensor(
                self.qubit1.get_hamiltonian(),
                self.qubit2.I,
                self.cavity.I,
            )
            self.H_q2 = tensor(
                self.qubit1.I,
                self.qubit2.get_hamiltonian(),
                self.cavity.I,
            )
            self.H_cav = tensor(
                self.qubit1.I,
                self.qubit2.I,
                self.cavity.get_hamiltonian(),
            )

        self.H_0: Qobj = self.H_q1 + self.H_q2 + self.H_cav

        # Equivalent qubit-cavity coupling
        H_qc1 = sum(
            self.g[0][m]
            * (self.a1.dag() * self.acav[m] + self.a1 * self.acav[m].dag())
            for m in range(self.cavity.M)
        )
        H_qc2 = sum(
            self.g[1][m]
            * (self.a2.dag() * self.acav[m] + self.a2 * self.acav[m].dag())
            for m in range(self.cavity.M)
        )

        self.H = self.H_q1 + self.H_q2 + self.H_cav + H_qc1 + H_qc2

    # -- State preparation --------------------------------------------------

    def prepare_ket11(self, T: float, sigma: float):
        """Prepare |11> state with DRAG pi pulses.

        Parameters
        ----------
        T : float
            Gate duration (ns).
        sigma : float
            Gaussian width (ns).

        Returns
        -------
        final_state : Qobj
            Prepared state.
        fidelity : float
        """
        t_lists = np.linspace(0, T, 1000)

        def Omega_I(t, args):
            return args["omega_I"] * np.exp(
                -0.5 * ((t - args["t0"]) / args["sigma"]) ** 2
            )

        def Omega_Q(t, args):
            return (
                args["omega_Q"]
                * np.exp(
                    -0.5 * ((t - args["t0"]) / args["sigma"]) ** 2
                )
                * (-(t - args["t0"]) / (args["sigma"] ** 2))
            )

        args1 = {
            "omega_I": np.pi / (sigma * np.sqrt(2 * np.pi)),
            "omega_Q": -np.pi
            / (sigma * np.sqrt(2 * np.pi) * self.qubit1.anharmonicity),
            "t0": T / 2,
            "sigma": sigma,
        }
        args2 = {
            "omega_I": np.pi / (sigma * np.sqrt(2 * np.pi)),
            "omega_Q": -np.pi
            / (sigma * np.sqrt(2 * np.pi) * self.qubit2.anharmonicity),
            "t0": T / 2,
            "sigma": sigma,
        }
        H_d = [
            [0.5 * (self.a1 + self.a1.dag()), Omega_I],
            [0.5 * (-1j) * (self.a1 - self.a1.dag()), Omega_Q],
            [0.5 * (self.a2 + self.a2.dag()), Omega_I],
            [0.5 * (-1j) * (self.a2 - self.a2.dag()), Omega_Q],
        ]
        # Shut off coupling
        g_shut = [
            [0 for _m in range(self.cavity.M)],
            [0 for _m in range(self.cavity.M)],
        ]
        self.control_g(g_shut)

        H = [self.H] + H_d
        psi0 = tensor(
            self.qubit1.state,
            self.qubit2.state,
            *[basis(self.cavity.n_levels[m], 0) for m in range(self.cavity.M)],
        )
        result = sesolve(H, psi0, t_lists, args=args1)
        final_state = result.states[-1]
        target_state = tensor(
            basis(self.qubit1.n_levels, 1),
            basis(self.qubit2.n_levels, 1),
            *[basis(self.cavity.n_levels[m], 0) for m in range(self.cavity.M)],
        )
        fidelity = np.abs(target_state.overlap(final_state)) ** 2
        return final_state, fidelity

    # -- iSWAP gate ---------------------------------------------------------

    def simulate_iSWAP(self, g: float):
        """Simulate iSWAP gate with tunable coupling.

        Parameters
        ----------
        g : float
            Coupling strength.

        Returns
        -------
        U_cal : Qobj
            Calibrated gate unitary.
        leakage : float
        F_loc : float
            Fidelity up to local Z rotations.
        """

        def f(t, args):
            T, T1, T2 = args["T"], args["T1"], args["T2"]
            if t < 0 or t > T:
                return [
                    [0 for _m in range(self.cavity.M)],
                    [0 for _m in range(self.cavity.M)],
                ]
            elif t < T1:
                return [
                    [np.sin(0.5 * np.pi * t / T1) ** 2 for _m in range(self.cavity.M)],
                    [np.sin(0.5 * np.pi * t / T1) ** 2 for _m in range(self.cavity.M)],
                ]
            elif t < T - T2:
                return [
                    [1 for _m in range(self.cavity.M)],
                    [1 for _m in range(self.cavity.M)],
                ]
            else:
                return [
                    [np.sin(0.5 * np.pi * (T - t) / T2) ** 2 for _m in range(self.cavity.M)],
                    [np.sin(0.5 * np.pi * (T - t) / T2) ** 2 for _m in range(self.cavity.M)],
                ]

        g_lists = [
            [g for _m in range(self.cavity.M)],
            [g for _m in range(self.cavity.M)],
        ]
        self.control_g(g_lists, is_time_dependent=True, f=f)
        J = sum(
            self.g[0][m]
            * self.g[1][m]
            / np.abs(self.qubit1.frequency - self.cavity.frequencies[m])
            for m in range(self.cavity.M)
        )
        T_ideal = np.pi / (2 * J)
        args = {"T": T_ideal, "T1": T_ideal / 10, "T2": T_ideal / 10}
        U = propagator(
            self.H, T_ideal, args=args, options={"nsteps": 1000000}
        )
        # Compute subspace unitary
        kets = [
            tensor(
                basis(self.qubit1.n_levels, i),
                basis(self.qubit2.n_levels, j),
                *[basis(self.cavity.n_levels[m], 0) for m in range(self.cavity.M)],
            )
            for i in range(2)
            for j in range(2)
        ]
        U_c = np.zeros((4, 4), dtype=complex)
        for i, ket_i in enumerate(kets):
            for j, ket_j in enumerate(kets):
                U_c[i, j] = ket_i.dag() * U * ket_j

        U_eff = Qobj(U_c, dims=[[2, 2], [2, 2]])
        from sqc.control.gates import ideal_iSWAP

        U_ideal = ideal_iSWAP()
        leakage = 1.0 - np.mean(np.sum(np.abs(U_eff.full()) ** 2, axis=0))
        F = (abs((U_ideal.dag() @ U_eff).tr()) ** 2 + 4) / (4 * 5)

        import scipy.linalg as la

        Uu, _ = la.polar(U_eff.full())
        Uu = Qobj(Uu, dims=U_eff.dims)

        def Rz(phi):
            return Qobj(
                [
                    [np.exp(-1j * phi / 2), 0],
                    [0, np.exp(1j * phi / 2)],
                ]
            )

        def avg_gate_fidelity(U_mat, V_mat):
            d = U_mat.shape[0]
            tr_val = (V_mat.dag() * U_mat).tr()
            return (abs(tr_val) ** 2 + d) / (d * (d + 1))

        from scipy.optimize import minimize

        def cost(x):
            a, b, c, d = x
            L = tensor(Rz(a), Rz(b))
            R = tensor(Rz(c), Rz(d))
            Uc = L * Uu * R
            return -avg_gate_fidelity(Uc, U_ideal)

        res = minimize(cost, x0=np.zeros(4))
        F_loc = -res.fun
        a, b, c, d_ = res.x
        L_opt = tensor(Rz(a), Rz(b))
        R_opt = tensor(Rz(c), Rz(d_))
        U_cal = L_opt * Uu * R_opt
        return U_cal, leakage, F_loc
