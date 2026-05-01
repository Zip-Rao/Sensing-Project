"""sqc.devices.transmon — Flux-tunable Transmon qubit models.

QubitSpec: frozen parameter-only description (preferred for new code).
TransmonQubit: backward-compatible mutable class wrapping QubitSpec.

Physical model (cf. Gao 2021 §II.B-C, Eq. 13-20):
    H = 4 EC n^2 - EJ(Phi) cos(varphi)
    EJ(Phi) = EJ_0 |cos(pi Phi/Phi_0)|             (SQUID, Eq. 20)
    omega_T = sqrt(8 EJ EC) - EC                   (Eq. 18)
    alpha    = -EC                                  (Eq. 18, sign convention)

Units: EC, EJ in rad*GHz (includes 2*pi factor); flux in Phi_0.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Optional

import numpy as np
from qutip import Qobj, basis, destroy, expect, num, qeye, sesolve


# ---------------------------------------------------------------------------
# QubitSpec — frozen parameter record
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class QubitSpec:
    """Pure parameter description of a flux-tunable Transmon qubit.

    Holds NO experimental state. All methods are pure: they take parameters
    and return values without mutating self.

    Parameters
    ----------
    name : str
        Human-readable identifier.
    EC : float
        Charging energy (rad*GHz).
    EJ : float
        Josephson energy at zero flux (rad*GHz).
    T1 : float
        Relaxation time (ns).
    T2 : float
        Dephasing time (ns).
    flux_bias : float
        Static flux bias (Phi_0). Default 0.
    n_levels : int
        Hilbert space truncation. Default 3.
    """

    name: str
    EC: float
    EJ: float
    T1: float
    T2: float
    flux_bias: float = 0.0
    n_levels: int = 3

    # -- flux-dependent EJ --------------------------------------------------

    def EJ_at(self, flux: float | None = None) -> float:
        """EJ at given flux, EJ(Phi) = EJ_0 * |cos(pi * Phi/Phi_0)|.

        Parameters
        ----------
        flux : float or None
            Flux in Phi_0. Uses self.flux_bias if None.

        Returns
        -------
        float
            Effective Josephson energy at the given flux.
        """
        f = self.flux_bias if flux is None else flux
        return self.EJ * abs(math.cos(math.pi * f))

    # -- spectroscopic properties -------------------------------------------

    def frequency(self, flux: float | None = None) -> float:
        """f_01(Phi) = sqrt(8 EJ(Phi) EC) - EC, in rad*GHz.

        Parameters
        ----------
        flux : float or None
            Flux in Phi_0. Uses self.flux_bias if None.

        Returns
        -------
        float
            Qubit transition frequency (rad*GHz).
        """
        return float(np.sqrt(8 * self.EJ_at(flux) * self.EC) - self.EC)

    def anharmonicity(self) -> float:
        """alpha = -EC (rad*GHz)."""
        return -self.EC

    def sensitivity(
        self, flux: float | None = None, delta: float = 1e-6
    ) -> float:
        """Numerical domega/dPhi via central difference.

        Parameters
        ----------
        flux : float or None
            Flux in Phi_0 at which to evaluate. Uses flux_bias if None.
        delta : float
            Small flux step (Phi_0).

        Returns
        -------
        float
            df_01/dPhi (rad*GHz / Phi_0).
        """
        f = self.flux_bias if flux is None else flux
        f_plus = self.frequency(f + delta)
        f_minus = self.frequency(f - delta)
        return (f_plus - f_minus) / (2 * delta)

    # -- helpers ------------------------------------------------------------

    def with_flux(self, new_flux: float) -> "QubitSpec":
        """Return a new QubitSpec with flux_bias=new_flux. Does not mutate self.

        Parameters
        ----------
        new_flux : float
            New flux bias in Phi_0.

        Returns
        -------
        QubitSpec
        """
        return QubitSpec(
            name=self.name,
            EC=self.EC,
            EJ=self.EJ,
            T1=self.T1,
            T2=self.T2,
            flux_bias=new_flux,
            n_levels=self.n_levels,
        )

    @staticmethod
    def optimal_work_point() -> float:
        """Flux at which |df/dPhi| is maximum, in units of Phi_0.

        Returns
        -------
        float
            Optimal flux bias in Phi_0.
        """
        return float(np.arctan(np.sqrt(2)) / np.pi)


# ---------------------------------------------------------------------------
# TransmonQubit — backward-compatible mutable class
# ---------------------------------------------------------------------------

class TransmonQubit:
    """Backward-compatible Transmon qubit class.

    Wraps QubitSpec and maintains mutable state for legacy callers
    (src/protocal.py, Simulation.ipynb). New code should use QubitSpec
    + HamiltonianBuilder directly.

    Parameters
    ----------
    EC : float
        Charging energy (rad*GHz).
    EJ : float
        Josephson energy at zero flux (rad*GHz).
    T1 : float
        Relaxation time (ns).
    T2 : float
        Dephasing time (ns).
    flux : float
        Static flux bias (Phi_0). Default 0.
    state : int or Qobj
        Initial state. Default 0 (ground).
    n_levels : int
        Hilbert space truncation. Default 3.
    name : str
        Qubit identifier. Default "Q".
    """

    def __init__(
        self,
        EC: float,
        EJ: float,
        T1: float,
        T2: float,
        flux: float = 0.0,
        state: int | Qobj = 0,
        n_levels: int = 3,
        name: str = "Q",
    ) -> None:
        # Internal QubitSpec
        self._spec = QubitSpec(
            name=name,
            EC=EC,
            EJ=EJ,
            T1=T1,
            T2=T2,
            flux_bias=flux,
            n_levels=n_levels,
        )

        # Legacy attributes (preserved verbatim for backward compatibility)
        self.EC = EC
        self.EJ = EJ * abs(math.cos(math.pi * flux))  # at current flux
        self.EJ_0 = EJ
        self.flux = flux
        self.n_levels = n_levels
        self.T1 = T1
        self.T2 = T2

        # Frequency and anharmonicity
        self.frequency: float = self._spec.frequency()
        self.anharmonicity: float = self._spec.anharmonicity()

        # Operators
        self.a: Qobj = destroy(n_levels)
        self.a_dag: Qobj = self.a.dag()
        self.n: Qobj = num(n_levels)
        self.I: Qobj = qeye(n_levels)

        # Hamiltonian and collapse ops
        self.hamiltonian: Qobj = self.get_hamiltonian()
        self.c_ops: list[Qobj] = self.get_collapse_operators()

        # Initial state
        if isinstance(state, int) and 0 <= state < n_levels:
            self.state: Qobj = basis(n_levels, state)
        elif isinstance(state, Qobj) and state.dims == [[n_levels], [1]]:
            self.state = state.unit()
        else:
            self.state = basis(n_levels, 0)

        # isinmag flag (set by qubit_in_mag)
        self.isinmag: bool = False

    # -- accessor ------------------------------------------------------------

    def spec(self) -> QubitSpec:
        """Return the internal QubitSpec."""
        return self._spec

    # -- Device ABC compatibility --------------------------------------------

    @property
    def name(self) -> str:
        return self._spec.name

    def hilbert_dim(self) -> int:
        return self.n_levels

    def hamiltonian_static(self) -> Qobj:
        return self.hamiltonian

    def collapse_operators(self) -> list[Qobj]:
        return self.c_ops

    # -- Static methods ------------------------------------------------------

    def calculate_frequency(self) -> float:
        """Calculate 0-1 transition frequency from current EJ, EC."""
        return self._spec.frequency()

    def calculate_anharmonicity(self) -> float:
        """Calculate anharmonicity = -EC."""
        return self._spec.anharmonicity()

    def frequency_sensitivity(self, flux: float, delta_flux: float = 1e-6) -> float:
        """Numerical df_01/dPhi via central difference.

        Parameters
        ----------
        flux : float
            Flux bias (Phi_0) at which to evaluate.
        delta_flux : float
            Small step for finite difference.

        Returns
        -------
        float
            d(omega_01)/dPhi (rad*GHz / Phi_0).
        """
        flux_plus = flux + delta_flux
        EJ_plus = self.EJ_0 * abs(math.cos(math.pi * flux_plus))
        f_plus = np.sqrt(8 * EJ_plus * self.EC) - self.EC

        flux_minus = flux - delta_flux
        EJ_minus = self.EJ_0 * abs(math.cos(math.pi * flux_minus))
        f_minus = np.sqrt(8 * EJ_minus * self.EC) - self.EC

        return (f_plus - f_minus) / (2 * delta_flux)

    # -- Hamiltonian ---------------------------------------------------------

    def get_hamiltonian(self) -> Qobj:
        """Time-independent Hamiltonian (lab frame).

        Returns
        -------
        Qobj
            H = (-EJ + 0.25*EC)*I + f_01*(n + 0.5*I) + (alpha/2)*(n^2 - n)
        """
        n_levels = self.n_levels
        H_0 = (-self.EJ + 0.25 * self.EC) * qeye(n_levels)
        H_1 = self.frequency * (self.n + 0.5 * qeye(n_levels))
        H_2 = (self.anharmonicity / 2.0) * (self.n * self.n - self.n)
        return H_0 + H_1 + H_2

    def get_hamiltonian_rwa(self, omega_d: float) -> Qobj:
        """Hamiltonian in rotating frame (RWA).

        Parameters
        ----------
        omega_d : float
            Drive frequency (rad*GHz).

        Returns
        -------
        Qobj
            H_rwa = Delta * n + (alpha/2) * (n^2 - n)
        """
        Delta = self.frequency - omega_d
        H_0 = Delta * self.n
        H_1 = (self.anharmonicity / 2.0) * (self.n * self.n - self.n)
        return H_0 + H_1

    # -- Collapse operators --------------------------------------------------

    def get_collapse_operators(self) -> list[Qobj]:
        """Lindblad collapse operators for T1, T2.

        Returns
        -------
        list[Qobj]
            [sqrt(gamma_1) * a, sqrt(gamma_phi) * n]
        """
        gamma_1 = 1.0 / self.T1
        gamma_phi = 1.0 / self.T2 - 0.5 * gamma_1
        c_ops = []
        c_ops.append(np.sqrt(gamma_1) * self.a)
        c_ops.append(np.sqrt(gamma_phi) * self.n)
        return c_ops

    # -- Noise ---------------------------------------------------------------

    def generate_1f_noise(
        self,
        t_lists: np.ndarray,
        amplitude: float,
        f_min: float,
        f_max: float,
    ) -> np.ndarray:
        """Generate 1/f noise time series.

        Parameters
        ----------
        t_lists : np.ndarray
            Time array (ns).
        amplitude : float
            Noise amplitude.
        f_min : float
            Minimum frequency (Hz).
        f_max : float
            Maximum frequency (Hz).

        Returns
        -------
        np.ndarray
            Noise time series.
        """
        from sqc.simulation.noise import generate_1f_noise

        return generate_1f_noise(t_lists, amplitude, f_min, f_max)

    # -- State projection ----------------------------------------------------

    def calculate_state_projection(self, target_state: int | Qobj) -> float:
        """Probability that the qubit is in target_state.

        Parameters
        ----------
        target_state : int or Qobj
            Target Fock state index or state vector.

        Returns
        -------
        float
            Projection probability.
        """
        if isinstance(target_state, int) and 0 <= target_state < self.n_levels:
            target = basis(self.n_levels, target_state)
        elif (
            isinstance(target_state, Qobj)
            and target_state.dims == [[self.n_levels], [1]]
        ):
            target = target_state.unit()
        else:
            raise ValueError("Invalid target_state")
        proj = target.dag() * target
        prob = expect(proj, self.state)
        return float(prob.real)

    # -- qubit_under_mag (legacy per-timestep) --------------------------------

    def qubit_under_mag(self, Phi_signal, is_noise: bool = False):
        """Per-timestep qubit array under flux signal.

        Parameters
        ----------
        Phi_signal : Signal-like
            Flux signal with .t_list and .value_at(t).
        is_noise : bool
            If True, add 1/f noise to the flux.

        Returns
        -------
        list[TransmonQubit]
            One new qubit per time point.
        """
        qubit = []
        noise_arr = np.zeros(len(Phi_signal.t_list))
        if is_noise:
            noise_arr = self.generate_1f_noise(
                Phi_signal.t_list, amplitude=0.001, f_min=1e-3, f_max=1e3
            )
        for i, t in enumerate(Phi_signal.t_list):
            q = TransmonQubit(
                EC=self.EC,
                EJ=self.EJ_0,
                T1=self.T1,
                T2=self.T2,
                flux=self.flux + Phi_signal.value_at(t) + noise_arr[i],
                state=self.state,
                n_levels=self.n_levels,
            )
            qubit.append(q)
        return qubit

    def qubit_under_mag_hamiltonian(
        self,
        qubit_t: list,
        t_list: np.ndarray,
        frame: int = 0,
        omega_d: float | None = None,
    ):
        """Build time-dependent H_list from per-timestep qubit array.

        Parameters
        ----------
        qubit_t : list[TransmonQubit]
            Pre-computed per-timestep qubits.
        t_list : np.ndarray
            Time array.
        frame : int
            0=lab, else rotating.
        omega_d : float or None
            Drive frequency for rotating frame.

        Returns
        -------
        list
            QuTiP H_list format.
        """
        H_list = []
        freq_coeffs = np.zeros(len(t_list))

        for i, _t in enumerate(t_list):
            qubit_current = qubit_t[i]
            if frame == 0:
                freq_coeffs[i] = qubit_current.frequency
            else:
                freq_coeffs[i] = qubit_current.frequency - omega_d

        H_list.append(
            qubit_current.anharmonicity * 0.5 * (self.n * self.n - self.n)
        )
        H_list.append([self.n + 0.5 * qeye(self.n_levels), freq_coeffs])
        return H_list

    # -- qubit_in_mag (optimised, delegates to HamiltonianBuilder) -----------

    def qubit_in_mag(
        self, Phi_signal, frame: int = 0, omega_d: float | None = None
    ) -> None:
        """Pre-compute H_list under flux signal (in-place, sets attributes).

        Legacy interface; delegates to HamiltonianBuilder.build().
        Sets self.isinmag, self.mag_signal, self.freq_coeffs, self.H_list.

        Parameters
        ----------
        Phi_signal : Signal-like
            Flux signal with .t_list and .signal attributes.
        frame : int
            0=lab, else rotating.
        omega_d : float or None
            Drive frequency for rotating frame.
        """
        from sqc.simulation.hamiltonian import HamiltonianBuilder

        H_list, _ = HamiltonianBuilder.build(
            qubit=self._spec,
            flux_signal=Phi_signal,
            pulse=None,
            frame="lab" if frame == 0 else "rotating",
            omega_d=omega_d,
        )
        self.isinmag = True
        self.mag_signal = Phi_signal
        # Extract freq_coeffs from H_list[1][1] for backward compat
        self.freq_coeffs = np.asarray(H_list[1][1])
        self.H_list = H_list

    # -- change_flux ---------------------------------------------------------

    def change_flux(self, flux: float) -> None:
        """Mutate flux bias in-place (legacy).

        Parameters
        ----------
        flux : float
            New flux bias (Phi_0).
        """
        self.flux = flux
        self.EJ = self.EJ_0 * abs(math.cos(math.pi * flux))
        self._spec = self._spec.with_flux(flux)
        self.frequency = self._spec.frequency()
        self.anharmonicity = self._spec.anharmonicity()
        self.hamiltonian = self.get_hamiltonian()

    # -- optimal_work_point --------------------------------------------------

    def optimal_work_point(self) -> float:
        """Return optimal flux working point in radians.

        Legacy: returns arctan(sqrt(2)) in radians (not Phi_0).
        Equivalent to QubitSpec.optimal_work_point() * pi.
        """
        return float(np.arctan(np.sqrt(2)))

    # -- Gate simulation ------------------------------------------------------

    def ideal_gate(self, theta: float, phi: float) -> Qobj:
        """Ideal single-qubit rotation in the {|0>, |1>} subspace.

        Parameters
        ----------
        theta : float
            Rotation angle (rad).
        phi : float
            Rotation axis phase (rad).

        Returns
        -------
        Qobj
            2x2 unitary operator.
        """
        a2 = destroy(2)
        U = (
            -1j
            * theta
            / 2.0
            * (
                (a2 + a2.dag()) * np.cos(phi)
                - 1j * (-a2 + a2.dag()) * np.sin(phi)
            )
        ).expm()
        return U

    def simulate_gate(
        self, theta: float, phi: float, T: float, sigma: float
    ):
        """DRAG single-qubit gate simulation.

        Parameters
        ----------
        theta : float
            Rotation angle (rad).
        phi : float
            Rotation axis phase (rad).
        T : float
            Gate duration (ns).
        sigma : float
            Gaussian width (ns).

        Returns
        -------
        tuple[Qobj, float, float]
            (U_eff, leakage, fidelity)
        """

        def Omega_I(t, args):
            A = args["A"]
            t0 = args["t0"]
            sig = args["sigma"]
            beta = args["beta"]
            ph = args["phi"]
            env = A * np.exp(-0.5 * ((t - t0) / sig) ** 2)
            drag = (
                -A
                * (t - t0)
                / (sig**2)
                * np.exp(-0.5 * ((t - t0) / sig) ** 2)
            )
            return env * np.cos(ph) - beta * drag * np.sin(ph)

        def Omega_Q(t, args):
            A = args["A"]
            t0 = args["t0"]
            sig = args["sigma"]
            beta = args["beta"]
            ph = args["phi"]
            env = A * np.exp(-0.5 * ((t - t0) / sig) ** 2)
            drag = (
                -A
                * (t - t0)
                / (sig**2)
                * np.exp(-0.5 * ((t - t0) / sig) ** 2)
            )
            return env * np.sin(ph) + beta * drag * np.cos(ph)

        t_0 = T / 2.0
        t_list = np.linspace(0, T, 1000)
        beta = -1.0 / self.anharmonicity
        args = {
            "A": theta / (sigma * np.sqrt(2 * np.pi)),
            "t0": t_0,
            "sigma": sigma,
            "beta": beta,
            "phi": phi,
        }
        ket_0 = basis(self.n_levels, 0)
        ket_1 = basis(self.n_levels, 1)
        a = destroy(self.n_levels)
        H = [
            self.get_hamiltonian_rwa(self.frequency),
            [0.5 * (a + a.dag()), Omega_I],
            [0.5 * (-1j) * (a - a.dag()), Omega_Q],
        ]
        result_0 = sesolve(H, ket_0, t_list, args=args)
        result_1 = sesolve(H, ket_1, t_list, args=args)
        final_state_0 = result_0.states[-1]
        final_state_1 = result_1.states[-1]

        amp00 = ket_0.overlap(final_state_0)
        amp01 = ket_0.overlap(final_state_1)
        amp10 = ket_1.overlap(final_state_0)
        amp11 = ket_1.overlap(final_state_1)

        U_eff = Qobj(
            np.array(
                [[amp00, amp01], [amp10, amp11]], dtype=complex
            ),
            dims=[[2], [2]],
        )
        U_ideal = self.ideal_gate(theta, phi)
        leakage = 1.0 - 0.5 * (np.sum(np.abs(U_eff.full()) ** 2))
        F = (abs((U_ideal.dag() @ U_eff).tr()) ** 2 + 2) / (2 * 3)
        return U_eff, leakage, F
