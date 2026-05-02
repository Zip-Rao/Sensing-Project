"""sqc.devices.chip — Multi-qubit chip topology and coupled systems.

ChipTopology: generic multi-qubit + resonator chip model (P5).
CoupledSystem: two qubits + coupler + multi-mode cavity (P1 port).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from qutip import Qobj, basis, destroy, propagator, qeye, sesolve, tensor

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


# ===========================================================================
# ChipTopology — multi-qubit chip descriptor (P5)
# ===========================================================================


@dataclass
class ChipTopology(Device):
    """Multi-qubit chip topology: qubits + resonators + couplings + control lines.

    Acts as a container for all devices on a chip and provides methods
    to embed single-device operators into the full Hilbert space.

    Per phase_5_handbook.md §3.1.

    Parameters
    ----------
    qubits : list[TransmonQubit]
        All qubits on the chip.
    resonators : list[Resonator]
        All readout/coupling resonators.
    couplings : dict[tuple[str, str], float]
        (qubit_name, resonator_name) → coupling g (rad/ns).
    control_lines : dict[str, ControlLine]
        Name → ControlLine mapping.
    transfer_matrix : TransferMatrix or None
        Frequency-dependent Z-line crosstalk matrix.
    """

    qubits: list = field(default_factory=list)  # list[TransmonQubit]
    resonators: list = field(default_factory=list)  # list[Resonator]
    couplings: dict = field(default_factory=dict)  # dict[tuple[str, str], float]
    control_lines: dict = field(default_factory=dict)
    transfer_matrix: Optional[object] = None  # TransferMatrix | None

    # -- Device ABC ----------------------------------------------------------

    @property
    def name(self) -> str:
        q_names = [q.name for q in self.qubits]
        return f"Chip({','.join(q_names)})" if q_names else "Chip(empty)"

    def hilbert_dim(self) -> int:
        """Total Hilbert space dimension (product of all subsystem dims)."""
        d = 1
        for q in self.qubits:
            d *= getattr(q, "n_levels", q.hilbert_dim())
        for r in self.resonators:
            d *= r.hilbert_dim()
        return d

    def lift_qubit_op(self, op: Qobj, qubit_name: str) -> Qobj:
        """Embed a single-qubit operator into the full chip Hilbert space.

        Builds a tensor product where the target qubit slot gets `op`
        and all other slots get identity.

        Parameters
        ----------
        op : Qobj
            Operator on the single-qubit Hilbert space (n_levels × n_levels).
        qubit_name : str
            Name of the target qubit.

        Returns
        -------
        Qobj
            Operator on the full chip Hilbert space.

        Raises
        ------
        ValueError
            If qubit_name not found in self.qubits.
        """
        # Collect all subsystem ops in order: qubits then resonators
        all_ops: list[Qobj] = []
        found = False
        for q in self.qubits:
            nq = getattr(q, "n_levels", q.hilbert_dim())
            if q.name == qubit_name:
                all_ops.append(op)
                found = True
            else:
                all_ops.append(qeye(nq))

        if not found:
            raise ValueError(
                f"Qubit '{qubit_name}' not found in chip topology. "
                f"Available: {[q.name for q in self.qubits]}"
            )

        for r in self.resonators:
            for n in r.n_levels:
                all_ops.append(qeye(n))

        return tensor(*all_ops)

    def hamiltonian_static(self) -> Qobj:
        """Sum of static (bare) Hamiltonians for all devices.

        H = sum_i lift(H_i) + sum_r lift(H_r)

        Note: coupling terms (g * a_q† a_r + h.c.) are NOT included
        in this basic implementation. For coupled systems, use the
        dedicated CoupledSystem class or extend this method.

        Returns
        -------
        Qobj
            Total static Hamiltonian.
        """
        H = None
        for q in self.qubits:
            H_q = getattr(q, "hamiltonian_static", q.get_hamiltonian)()
            H_q_full = self.lift_qubit_op(H_q, q.name)
            H = H_q_full if H is None else H + H_q_full

        # Resonator terms (each mode is lifted into the full space)
        if H is None:
            H = 0 * self._identity()

        for r in self.resonators:
            H_r = r.hamiltonian_static()
            H_r_full = self._lift_resonator_op(H_r, r)
            H += H_r_full

        return H

    def collapse_operators(self) -> list[Qobj]:
        """All Lindblad collapse operators for the chip.

        Lifts each qubit's collapse operators into the full space.

        Returns
        -------
        list[Qobj]
            Collapse operators for mesolve.
        """
        c_ops_all: list[Qobj] = []
        for q in self.qubits:
            c_ops_q = getattr(q, "collapse_operators", q.get_collapse_operators)()
            for c in c_ops_q:
                c_ops_all.append(self.lift_qubit_op(c, q.name))

        for r in self.resonators:
            for c in r.collapse_operators():
                c_ops_all.append(self._lift_resonator_op(c, r))

        return c_ops_all

    # -- helpers --------------------------------------------------------------

    def _identity(self) -> Qobj:
        """Build the identity operator on the full Hilbert space."""
        all_ops: list[Qobj] = []
        for q in self.qubits:
            nq = getattr(q, "n_levels", q.hilbert_dim())
            all_ops.append(qeye(nq))
        for r in self.resonators:
            for n in r.n_levels:
                all_ops.append(qeye(n))
        return tensor(*all_ops) if all_ops else qeye(1)

    def _lift_resonator_op(self, op: Qobj, resonator) -> Qobj:
        """Embed a resonator operator into the full chip space."""
        all_ops: list[Qobj] = []
        for q in self.qubits:
            nq = getattr(q, "n_levels", q.hilbert_dim())
            all_ops.append(qeye(nq))
        # The resonator's op is already on its own full tensor-product space,
        # so we can't easily decompose it into per-mode identities.
        # Instead, we tensor with identities for qubit spaces.
        qubit_id = tensor(*[qeye(getattr(q, "n_levels", q.hilbert_dim()))
                           for q in self.qubits])
        return tensor(qubit_id, op) if qubit_id.dims != [[1], [1]] else op

    # -- factory methods ------------------------------------------------------

    @classmethod
    def from_legacy_coupled_system(cls, coupled_system) -> "ChipTopology":
        """Build a ChipTopology from a legacy Coupled_System / CoupledSystem.

        Extracts qubits and resonator from the coupled system object.

        Parameters
        ----------
        coupled_system : CoupledSystem or src.qubit.Coupled_System
            Legacy coupled system object.

        Returns
        -------
        ChipTopology
        """
        qubits = []
        for attr_name in ("qubit1", "qubit2"):
            if hasattr(coupled_system, attr_name):
                qubits.append(getattr(coupled_system, attr_name))

        resonator = None
        for attr_name in ("cavity", "resonator"):
            if hasattr(coupled_system, attr_name):
                resonator = getattr(coupled_system, attr_name)
                break

        resonators = [resonator] if resonator is not None else []

        return cls(qubits=qubits, resonators=resonators)

    def get_qubit(self, name: str):
        """Get a qubit by name.

        Parameters
        ----------
        name : str
            Qubit name.

        Returns
        -------
        TransmonQubit

        Raises
        ------
        ValueError
            If no qubit with that name is found.
        """
        for q in self.qubits:
            if q.name == name:
                return q
        raise ValueError(
            f"Qubit '{name}' not found. Available: {[q.name for q in self.qubits]}"
        )
