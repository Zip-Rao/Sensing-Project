"""sqc.control.pulse — PulseBase ABC + Pulse + CompositePulse.

Pulse: single control pulse (lab or rotating frame).
CompositePulse: concatenation of multiple Pulses.

Phase 1: verbatim port from src/pulse.py. get_kernel() retained
as deprecated method; Phase 2 will refactor to KernelEstimator.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
from qutip import (
    Qobj, QobjEvo, basis, create, destroy, expect, mesolve, qeye, sigmaz,
)


# ---------------------------------------------------------------------------
# PulseBase ABC
# ---------------------------------------------------------------------------

class PulseBase(ABC):
    """Abstract base for any control pulse.

    Subclasses must implement hamiltonian (QuTiP list format),
    t_list (time axis), and frame (reference frame).
    """

    @property
    @abstractmethod
    def hamiltonian(self) -> list:
        """QuTiP list-format Hamiltonian: [[op, coeff_arr], ...]."""
        ...

    @property
    @abstractmethod
    def t_list(self):
        """Time axis (np.ndarray) for this pulse."""
        ...

    @property
    @abstractmethod
    def frame(self) -> int:
        """Reference frame: 0 = lab, 1 = rotating."""
        ...


# ---------------------------------------------------------------------------
# Pulse — full implementation (verbatim port from src/pulse.py)
# ---------------------------------------------------------------------------

class Pulse(PulseBase):
    """Single control pulse.

    Parameters
    ----------
    frame : int
        Reference frame: 0 = lab, 1 = rotating.
    omega_d : float
        Drive frequency (rad*GHz).
    phase : float
        Rotation axis phase (rad).
    Omega : Signal-like or float
        Rabi frequency envelope (GHz).
    is_rwa : bool
        Whether to use rotating wave approximation.
    qubit : TransmonQubit or None
        Qubit for n_levels resolution. Defaults to 2-level if None.
    """

    def __init__(
        self,
        frame: int = 0,
        omega_d: float = 0.0,
        phase: float = 0.0,
        Omega=None,
        is_rwa: bool = True,
        qubit=None,
    ) -> None:
        self.frame: int = frame
        self.omega_d: float = omega_d
        self.phase: float = phase
        self.Omega = Omega
        self.is_rwa: bool = is_rwa
        self.qubit = qubit

        if qubit is not None:
            self.n_levels: int = qubit.n_levels
        else:
            self.n_levels: int = 2  # default: backward-compatible

        self.hamiltonian, self.t_list = self.get_hamiltonian()

    # -- Rabi frequency -----------------------------------------------------

    def get_Rabi_frequency(self, t: float) -> float:
        """Get Rabi frequency at time t.

        Parameters
        ----------
        t : float
            Query time (ns).

        Returns
        -------
        float
            Omega(t) in GHz.
        """
        Omega = self.Omega
        if hasattr(Omega, "value_at"):
            return Omega.value_at(t)
        elif isinstance(Omega, (int, float)):
            return float(Omega)
        return 0.0

    # -- Hamiltonian at a single time ----------------------------------------

    def get_hamiltonian_at(self, t: float) -> Qobj:
        """Hamiltonian at a single time point.

        Parameters
        ----------
        t : float
            Time (ns).

        Returns
        -------
        Qobj
        """
        n = self.n_levels
        a = destroy(n)
        adag = create(n)
        Omega_t = self.get_Rabi_frequency(t)
        if self.frame == 0:  # lab frame
            phase_term = self.omega_d * t + self.phase
            return Omega_t * np.cos(phase_term) * (a + adag)
        else:  # rotating frame
            H_rwa_co = Omega_t / 2.0 * (
                a * np.exp(1j * self.phase)
                + adag * np.exp(-1j * self.phase)
            )
            if self.is_rwa:
                return H_rwa_co
            else:
                phase_cr = 2 * self.omega_d * t + self.phase
                H_cr = Omega_t / 2.0 * (
                    a * np.exp(-1j * phase_cr)
                    + adag * np.exp(1j * phase_cr)
                )
                return H_rwa_co + H_cr

    # -- Full Hamiltonian (list format) -------------------------------------

    def get_hamiltonian(self):
        """Build full H_list for QuTiP consumption.

        Returns
        -------
        H_list : list
            QuTiP list-format Hamiltonian.
        t_list : np.ndarray
            Time axis.
        """
        n = self.n_levels
        a = destroy(n)
        adag = create(n)
        t_list = np.asarray(self.Omega.t_list)

        if hasattr(self.Omega, "signal"):
            Omega = np.array(self.Omega.signal, dtype=float)
        else:
            Omega = float(self.Omega)

        if self.frame == 0:  # lab frame
            coeff = Omega * np.cos(self.omega_d * t_list + self.phase)
            H = [[a + adag, coeff]]
        else:  # rotating frame
            H1 = 0.5 * (
                a * np.exp(1j * self.phase) + adag * np.exp(-1j * self.phase)
            )
            coeff_rwa = np.ones_like(t_list) if np.isscalar(Omega) else Omega
            if not np.isscalar(Omega):
                coeff_rwa = Omega.astype(complex)
            else:
                coeff_rwa = np.full(len(t_list), float(Omega), dtype=complex)
            if self.is_rwa:
                H = [[H1, coeff_rwa]]
            else:
                H_cr1 = 0.5 * a * np.exp(-1j * self.phase)
                H_cr2 = 0.5 * adag * np.exp(1j * self.phase)
                coeff_cr1 = (
                    Omega * np.exp(-2j * self.omega_d * t_list)
                    if isinstance(self.Omega, (int, float))
                    else Omega * np.exp(-2j * self.omega_d * t_list)
                )
                coeff_cr2 = (
                    Omega.astype(complex) * np.exp(2j * self.omega_d * t_list)
                    if not np.isscalar(Omega)
                    else np.full(len(t_list), float(Omega), dtype=complex)
                    * np.exp(2j * self.omega_d * t_list)
                )
                H = [[H1, coeff_rwa], [H_cr1, coeff_cr1], [H_cr2, coeff_cr2]]
        return H, t_list

    # -- Rotation angle -----------------------------------------------------

    def get_angle(self, qubit=None):
        """Get rotation angle from pulse action on qubit.

        Parameters
        ----------
        qubit : TransmonQubit or None

        Returns
        -------
        theta : float
            Rotation angle (rad).
        phi : float
            Rotation phase (rad).
        """
        if qubit is not None:
            self.qubit = qubit
            self.n_levels = qubit.n_levels

        n_levels = self.n_levels
        psi_0 = self.qubit.state if self.qubit is not None else basis(n_levels, 0)
        t_list = self.Omega.t_list
        dt = t_list[1] - t_list[0]
        U = qeye(n_levels)
        for i in range(len(t_list) - 1):
            t_i = (t_list[i] + t_list[i + 1]) / 2.0
            H_i = self.get_hamiltonian_at(t_i)
            U_i = (-1j * H_i * dt).expm()
            U = U_i * U
        psi_final = U * psi_0
        a_mat = U[0, 0]
        a_mat = np.real(a_mat)
        angle = 2 * np.arccos(a_mat)
        c0 = psi_final[0, 0]
        c1 = psi_final[1, 0]
        norm = np.sqrt(abs(c0) ** 2 + abs(c1) ** 2)
        p1 = abs(c1) ** 2 / norm**2
        theta = 2 * np.arcsin(np.sqrt(p1))
        phi = np.angle(c1) - np.angle(c0)
        return theta, phi

    def get_angle_simple(self) -> float:
        """Simple Rabi angle = integral of Omega(t).

        Returns
        -------
        float
        """
        num_points = 1000
        t_list = np.linspace(
            self.Omega.t_list[0], self.Omega.t_list[-1], num_points
        )
        Omegas = [self.get_Rabi_frequency(t) for t in t_list]
        angle = np.trapezoid(Omegas, t_list)
        return angle

    # -- get_kernel (DEPRECATED, verbatim port for backward compat) ---------

    def get_kernel(self, qubit):
        """DEPRECATED. Phase 2 will use KernelEstimator.

        Estimate control kernel by perturbing with narrow stimuli.

        Parameters
        ----------
        qubit : TransmonQubit
            Qubit for simulation.

        Notes
        -----
        Sets self.t_samples and self.kernel in-place.
        """
        kernel = []
        t_list = self.Omega.t_list
        samples = range(0, len(t_list), 1)
        t_samples = t_list[::1]
        psi_e = basis(qubit.n_levels, 1)

        if self.frame == 0:
            H_0 = QobjEvo(qubit.hamiltonian)
        else:
            H_0 = QobjEvo(qubit.get_hamiltonian_rwa(qubit.frequency))
        H_pulse = QobjEvo(self.hamiltonian, tlist=t_list, order=1)
        H_base = H_0 + H_pulse

        result_base = mesolve(
            H_base, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()]
        )
        p_e_base = result_base.expect[0][-1]

        # Import Signal from sqc for constructing stimulus
        from sqc.control.flux_signal import FluxSignal

        for i in samples:
            t_i = t_list[i]
            stim_pulse = FluxSignal(
                type=3,
                t_list=t_list,
                amplitude=self._kernel_stimulus_amplitude(),
                center=t_i,
                width=self._kernel_stimulus_width(),
            )
            H_stim = QobjEvo(
                [sigmaz(), stim_pulse.signal],
                tlist=stim_pulse.t_list,
                order=1,
            )
            stim_area = (
                stim_pulse.params["amplitude"]
                * stim_pulse.params["width"]
                / 4.0
                * np.sqrt(2 * np.pi)
            )
            H_total = H_0 + H_pulse + H_stim
            result_stim = mesolve(
                H_total, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()]
            )
            p_e_stim = result_stim.expect[0][-1]
            kernel.append((p_e_stim - p_e_base) / stim_area)

        self.t_samples = t_samples
        self.kernel = kernel

    def _kernel_stimulus_amplitude(self) -> float:
        """Stimulus amplitude for kernel estimation."""
        return 1.0  # match original Pulse.get_kernel

    def _kernel_stimulus_width(self) -> float:
        """Stimulus width for kernel estimation."""
        return 2.0  # match original Pulse.get_kernel


# ---------------------------------------------------------------------------
# CompositePulse — verbatim port from src/pulse.py
# ---------------------------------------------------------------------------

class CompositePulse(PulseBase):
    """Sequence of multiple Pulses concatenated in time.

    Parameters
    ----------
    pulses : list[Pulse]
        Pulses to concatenate in order.
    """

    def __init__(self, pulses: list[Pulse]) -> None:
        self.pulses: list[Pulse] = list(pulses)
        self.t_list: np.ndarray = self.get_t_list()
        self.frame: int = pulses[0].frame
        self.omega_d: float = pulses[0].omega_d
        self.hamiltonian, _ = self.get_hamiltonian()

    def get_t_list(self) -> np.ndarray:
        """Compute merged time axis from sub-pulses.

        Returns
        -------
        np.ndarray
        """
        t_list = []
        curr = 0.0
        for pulse in self.pulses:
            pulse_list = [t + curr for t in pulse.Omega.t_list]
            t_list.extend(pulse_list)
            if pulse_list:
                curr = pulse_list[-1] + 1e-9
        return np.array(t_list)

    def get_Omega(self, t: float) -> float:
        """Rabi frequency at time t.

        Parameters
        ----------
        t : float
            Query time (ns).

        Returns
        -------
        float
        """
        curr = 0.0
        for pulse in self.pulses:
            dur = pulse.Omega.t_list[-1]
            if curr <= t <= curr + dur:
                return pulse.get_Rabi_frequency(t - curr)
            curr += dur
        return 0.0

    def get_hamiltonian_at(self, t: float) -> Qobj:
        """Hamiltonian at time t.

        Parameters
        ----------
        t : float
            Query time.

        Returns
        -------
        Qobj
        """
        curr = 0.0
        for pulse in self.pulses:
            dur = pulse.Omega.t_list[-1]
            if curr <= t <= curr + dur:
                return pulse.get_hamiltonian_at(t - curr)
            curr += dur
        return Qobj(np.zeros((pulse.n_levels, pulse.n_levels)))

    def get_hamiltonian(self):
        """Build full H_list on merged time axis.

        Returns
        -------
        H_list : list
            QuTiP list-format.
        t_global : np.ndarray
            Merged time axis.
        """
        t_global = np.asarray(self.t_list)
        N = len(t_global)
        hamiltonian = []
        curr = 0.0
        for pulse in self.pulses:
            term_curr, t_local = pulse.hamiltonian, pulse.t_list
            duration = t_local[-1] - t_local[0]
            for op, coeff_local in term_curr:
                coeff_global = np.zeros(N, dtype=complex)
                t_loc = t_global - curr
                mask = (t_loc >= t_local[0]) & (t_loc <= t_local[-1])
                indices = np.clip(
                    np.searchsorted(t_local, t_loc[mask]) - 1,
                    0,
                    len(coeff_local) - 1,
                )
                coeff_global[mask] = coeff_local[indices]
                hamiltonian.append([op, coeff_global])
            curr += duration
        return hamiltonian, t_global

    def plot(self):
        """Plot the composite pulse envelope."""
        import matplotlib.pyplot as plt

        t_list = self.t_list
        Omega_values = [self.get_Omega(t) for t in t_list]
        plt.figure(figsize=(8, 4))
        plt.plot(t_list, Omega_values)
        plt.xlabel("Time (ns)")
        plt.ylabel("Rabi Frequency (GHz)")
        plt.title("Composite Pulse Rabi Frequency Envelope")
        plt.grid()
        plt.show()

    # -- get_kernel (DEPRECATED, verbatim port) -----------------------------

    def get_kernel(self, qubit):
        """DEPRECATED. Phase 2 will use KernelEstimator.

        Estimate control kernel for composite pulse.

        Parameters
        ----------
        qubit : TransmonQubit

        Notes
        -----
        Sets self.t_samples and self.kernel in-place.
        """
        from sqc.control.flux_signal import FluxSignal

        kernel = []
        t_list = self.t_list
        for i in range(len(t_list) - 1):
            if t_list[i] == t_list[i + 1]:
                print(
                    f"Warning: Duplicate time points at index {i} and "
                    f"{i+1} with time {t_list[i]}. "
                    "This may cause issues in kernel calculation."
                )
        samples = range(0, len(t_list), 1)
        t_samples = t_list[::1]
        psi_e = basis(qubit.n_levels, 1)

        if self.frame == 0:
            H_0 = QobjEvo(qubit.hamiltonian)
        else:
            H_0 = QobjEvo(qubit.get_hamiltonian_rwa(qubit.frequency))
        H_pulse = QobjEvo(self.hamiltonian, tlist=t_list, order=1)
        H_base = H_0 + H_pulse
        result_base = mesolve(
            H_base, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()]
        )
        p_e_base = result_base.expect[0][-1]
        for i in samples:
            t_i = t_list[i]
            stim_pulse = FluxSignal(
                type=3,
                t_list=t_list,
                amplitude=0.0215,
                center=t_i,
                width=3,
            )
            stim_area = np.trapezoid(stim_pulse.signal, stim_pulse.t_list)
            qubit_t = qubit.qubit_under_mag(stim_pulse)
            H_stim = QobjEvo(
                qubit.qubit_under_mag_hamiltonian(
                    qubit_t, stim_pulse.t_list, self.frame, self.omega_d
                ),
                tlist=stim_pulse.t_list,
                order=1,
            )
            H_total = H_0 + H_pulse + H_stim
            result_stim = mesolve(
                H_total, qubit.state, t_list, [], e_ops=[psi_e * psi_e.dag()]
            )
            p_e_stim = result_stim.expect[0][-1]
            kernel.append((p_e_stim - p_e_base) / stim_area)

        self.t_samples = t_samples
        self.kernel = kernel
