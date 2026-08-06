"""sqc.control.pulse — PulseBase ABC + Pulse + CompositePulse.

Pulse: single control pulse (lab or rotating frame).
CompositePulse: concatenation of multiple Pulses.

Phase 1: verbatim port from src/pulse.py.
Phase 10.5: get_kernel() deprecated → KernelEstimator shim.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
from qutip import (
    Qobj, basis, create, destroy, qeye,
)


# ---------------------------------------------------------------------------
# PulseBase ABC
# ---------------------------------------------------------------------------

class PulseBase(ABC):
    """Abstract base for any control pulse.

    Subclasses must provide the following *instance attributes*:
      - hamiltonian : QuTiP list-format Hamiltonian
      - t_list : time axis (np.ndarray)
      - frame : reference frame (0=lab, 1=rotating)

    This base class does NOT enforce these via @abstractmethod because
    Python's ABC machinery cannot check instance attributes set in
    __init__. Instead, subclasses set these attributes in __init__ and
    the base class verifies them in __post_init_checks__().
    """

    def __init_subclass__(cls, **kwargs):
        """Verify that subclass defines required attributes."""
        super().__init_subclass__(**kwargs)
        # Defer check to instance creation time via __post_init_checks__


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
        Rabi frequency envelope — I (in-phase) component (GHz).
    Omega_Q : Signal-like, float, or None
        Q (quadrature) component for DRAG pulses. If None, the
        pulse is single-quadrature (backward-compatible).
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
        Omega_Q=None,
        is_rwa: bool = True,
        qubit=None,
        trigger: float = 0.0,
    ) -> None:
        self.frame: int = frame
        self.omega_d: float = omega_d
        self.phase: float = phase
        self.Omega = Omega
        self.Omega_Q = Omega_Q
        self.is_rwa: bool = is_rwa
        self.qubit = qubit
        self.trigger: float = trigger

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
        Omega_I_t = self.get_Rabi_frequency(t)

        # Q component
        Omega_Q_t = 0.0
        if self.Omega_Q is not None:
            if hasattr(self.Omega_Q, "value_at"):
                Omega_Q_t = self.Omega_Q.value_at(t)
            elif isinstance(self.Omega_Q, (int, float)):
                Omega_Q_t = float(self.Omega_Q)

        if self.frame == 0:  # lab frame
            coeff = (Omega_I_t * np.cos(self.omega_d * t + self.phase)
                     + Omega_Q_t * np.sin(self.omega_d * t + self.phase))
            return coeff * (a + adag)
        else:  # rotating frame
            H_rwa = Omega_I_t / 2.0 * (
                a * np.exp(1j * self.phase)
                + adag * np.exp(-1j * self.phase)
            )
            if Omega_Q_t != 0.0:
                H_rwa += Omega_Q_t / 2.0 * (
                    -1j * a * np.exp(1j * self.phase)
                    + 1j * adag * np.exp(-1j * self.phase)
                )
            if self.is_rwa:
                return H_rwa
            else:
                phase_cr = 2 * self.omega_d * t + self.phase
                H_cr = Omega_I_t / 2.0 * (
                    a * np.exp(-1j * phase_cr)
                    + adag * np.exp(1j * phase_cr)
                )
                return H_rwa + H_cr

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
            Omega_I = np.array(self.Omega.signal, dtype=float)
        else:
            Omega_I = float(self.Omega)

        # -- Q quadrature (DRAG) ----------------------------------------------
        has_q = self.Omega_Q is not None
        if has_q:
            if hasattr(self.Omega_Q, "signal"):
                Omega_Q = np.array(self.Omega_Q.signal, dtype=float)
            else:
                Omega_Q = float(self.Omega_Q)

        if self.frame == 0:  # lab frame
            coeff = Omega_I * np.cos(self.omega_d * t_list + self.phase)
            if has_q:
                coeff = coeff + Omega_Q * np.sin(
                    self.omega_d * t_list + self.phase
                )
            H = [[a + adag, coeff]]
        else:  # rotating frame
            # I operator: σ_φ / 2
            H_I = 0.5 * (
                a * np.exp(1j * self.phase) + adag * np.exp(-1j * self.phase)
            )
            # Q operator: σ_{φ+π/2} / 2
            H_Q = 0.5 * (
                -1j * a * np.exp(1j * self.phase)
                + 1j * adag * np.exp(-1j * self.phase)
            )

            def _coeff(arr_or_scalar):
                if np.isscalar(arr_or_scalar):
                    return np.full(len(t_list), float(arr_or_scalar),
                                   dtype=complex)
                return arr_or_scalar.astype(complex)

            if self.is_rwa:
                H = [[H_I, _coeff(Omega_I)]]
                if has_q:
                    H.append([H_Q, _coeff(Omega_Q)])
            else:
                if has_q:
                    raise NotImplementedError(
                        "DRAG with is_rwa=False is not yet supported."
                    )
                H_cr1 = 0.5 * a * np.exp(-1j * self.phase)
                H_cr2 = 0.5 * adag * np.exp(1j * self.phase)
                coeff_cr1 = (
                    Omega_I * np.exp(-2j * self.omega_d * t_list)
                    if isinstance(self.Omega, (int, float))
                    else Omega_I * np.exp(-2j * self.omega_d * t_list)
                )
                coeff_cr2 = (
                    Omega_I.astype(complex) * np.exp(2j * self.omega_d * t_list)
                    if not np.isscalar(Omega_I)
                    else np.full(len(t_list), float(Omega_I), dtype=complex)
                    * np.exp(2j * self.omega_d * t_list)
                )
                H = [[H_I, _coeff(Omega_I)],
                     [H_cr1, coeff_cr1], [H_cr2, coeff_cr2]]
        return H, t_list

    # -- Global time axis projection ----------------------------------------

    def hamiltonian_on(self, t_global: np.ndarray) -> list:
        """Project local Hamiltonian onto global time axis.

        Computes t_loc = t_global - self.trigger, then linearly
        interpolates the local coefficients onto t_global within the
        pulse's local time window [0, t_list[-1]].  Outside that
        window the contribution is zero.

        Complex coefficients (rotating frame) are interpolated
        separately on real and imaginary parts.

        Parameters
        ----------
        t_global : np.ndarray
            Global time axis (ns).

        Returns
        -------
        list
            QuTiP list-format Hamiltonian [[op, coeffs], ...] where
            each coeffs array has length len(t_global).
        """
        H_local = self.hamiltonian
        result = []
        t_loc = t_global - self.trigger
        mask = (t_loc >= self.t_list[0]) & (t_loc <= self.t_list[-1])
        for op, coeff_local in H_local:
            if np.iscomplexobj(coeff_local):
                coeff_global = np.zeros(len(t_global), dtype=complex)
                coeff_global[mask] = (
                    np.interp(t_loc[mask], self.t_list, coeff_local.real)
                    + 1j * np.interp(t_loc[mask], self.t_list, coeff_local.imag)
                )
            else:
                coeff_global = np.zeros(len(t_global), dtype=float)
                coeff_global[mask] = np.interp(
                    t_loc[mask], self.t_list, coeff_local
                )
            result.append([op, coeff_global])
        return result

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
        t_list = np.asarray(self.Omega.t_list, dtype=float)
        Omegas = [self.get_Rabi_frequency(t) for t in t_list]
        angle = np.trapezoid(Omegas, t_list)
        return angle

    # -- Phase shift (Virtual Z) -------------------------------------------

    def with_phase_shift(self, phi_z: float) -> "Pulse":
        """Return a new Pulse with phase incremented by *phi_z*.

        Parameters
        ----------
        phi_z : float
            Phase increment (rad).

        Returns
        -------
        Pulse
            New pulse with ``phase = old_phase + phi_z``.
        """
        return Pulse(
            frame=self.frame,
            omega_d=self.omega_d,
            phase=self.phase + phi_z,
            Omega=self.Omega,
            Omega_Q=self.Omega_Q,
            is_rwa=self.is_rwa,
            qubit=self.qubit,
            trigger=self.trigger,
        )

    # -- Copy helper -------------------------------------------------------

    def _copy(self) -> "Pulse":
        """Return an exact copy of this Pulse (no phase change)."""
        return Pulse(
            frame=self.frame,
            omega_d=self.omega_d,
            phase=self.phase,
            Omega=self.Omega,
            Omega_Q=self.Omega_Q,
            is_rwa=self.is_rwa,
            qubit=self.qubit,
            trigger=self.trigger,
        )

    # -- Multi-kick phase shift (Virtual Z) --------------------------------

    def with_phase_kicks(self, kicks: list[tuple[float, float]]) -> "Pulse":
        """Return new Pulse with cumulative phase from multiple kicks.

        A single Pulse has no sub-pulse structure, so all kicks are applied
        as a single cumulative phase shift.

        Parameters
        ----------
        kicks : list[tuple[float, float]]
            List of ``(t_kick, delta_phi)`` tuples.

        Returns
        -------
        Pulse
        """
        total = sum(dphi for _, dphi in kicks)
        return self.with_phase_shift(total) if total != 0.0 else self._copy()

    # -- get_kernel (DEPRECATED shim → KernelEstimator, Phase 10.5) --------

    def get_kernel(self, qubit, t_samples=None, *,
                   _suppress_deprecation: bool = False):
        """[DEPRECATED] Use ``KernelEstimator(mode='flux', method='exp').estimate()``.

        This shim forwards to :class:`sqc.reconstruction.KernelEstimator` for
        backward compatibility.  New code should use ``KernelEstimator``
        directly.

        Parameters
        ----------
        qubit : TransmonQubit
            Qubit for simulation.
        t_samples : np.ndarray or None
            Time points to evaluate kernel at.  If None, uses pulse time axis.
        _suppress_deprecation : bool
            Internal flag to silence the deprecation warning when called
            from within the shim itself.

        Returns
        -------
        t_samples : np.ndarray
        kernel : np.ndarray

        Notes
        -----
        Sets ``self.t_samples`` and ``self.kernel`` in-place for backward compat.
        """
        import warnings
        if not _suppress_deprecation:
            warnings.warn(
                "Pulse.get_kernel() is deprecated. Use "
                "sqc.reconstruction.KernelEstimator(mode='flux', method='exp', "
                "order=1).estimate(pulse, qubit) instead.",
                DeprecationWarning, stacklevel=2,
            )
        from sqc.reconstruction.kernel import KernelEstimator
        est = KernelEstimator(
            mode='flux', method='exp', order=1,
            deprecation_warn_legacy=False,
        )
        result = est.estimate_full(self, qubit, t_samples=t_samples)
        # Store on self for backward compat
        self.t_samples = result.t_samples
        self.kernel = list(result.k1)  # preserve list type for compat
        return result.t_samples, result.k1


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

        *Deprecated*: prefer ``hamiltonian_on(t_global)`` instead.
        This method is retained for backward compatibility only.

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
                curr = pulse_list[-1]  # P7.5: removed 1e-9 separator
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

    def hamiltonian_on(self, t_global: np.ndarray) -> list:
        """Collect all child pulse contributions on the global time axis.

        Each child Pulse carries its own ``trigger``; this method simply
        delegates to each child's ``hamiltonian_on(t_global)`` and merges
        the resulting H-lists.

        Parameters
        ----------
        t_global : np.ndarray
            Global time axis (ns).

        Returns
        -------
        list
            QuTiP list-format Hamiltonian (merged from all child pulses).
        """
        merged = []
        for pulse in self.pulses:
            merged.extend(pulse.hamiltonian_on(t_global))
        return merged

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

    # -- Phase shift (Virtual Z) -------------------------------------------

    def with_phase_shift(
        self, phi_z: float, from_time: float | None = None
    ) -> "CompositePulse":
        """Return a new CompositePulse with sub-pulses phase-shifted.

        Parameters
        ----------
        phi_z : float
            Phase increment (rad).
        from_time : float or None
            Only sub-pulses whose ``trigger >= from_time`` are shifted.
            If None, all sub-pulses are shifted.

        Returns
        -------
        CompositePulse
        """
        new_pulses = []
        for p in self.pulses:
            if from_time is None or p.trigger >= from_time:
                new_pulses.append(p.with_phase_shift(phi_z))
            else:
                # Exact copy without phase change
                new_pulses.append(Pulse(
                    frame=p.frame,
                    omega_d=p.omega_d,
                    phase=p.phase,
                    Omega=p.Omega,
                    Omega_Q=p.Omega_Q if hasattr(p, 'Omega_Q') else None,
                    is_rwa=p.is_rwa,
                    qubit=p.qubit,
                    trigger=p.trigger,
                ))
        return CompositePulse(new_pulses)

    # -- Multi-kick phase shift (Virtual Z) --------------------------------

    def with_phase_kicks(
        self, kicks: list[tuple[float, float]]
    ) -> "CompositePulse":
        """Return new CompositePulse with sub-pulses phase-shifted by
        multiple kicks.

        Generalises :meth:`with_phase_shift` (the single-kick special
        case).  For each sub-pulse, phase contributions from all kicks
        whose time ≤ the sub-pulse's ``trigger`` are summed.

        Parameters
        ----------
        kicks : list[tuple[float, float]]
            List of ``(t_kick, delta_phi)`` tuples, sorted by ``t_kick``
            ascending.  A sub-pulse at trigger ``t`` accumulates
            ``sum(delta_phi for t_k, delta_phi in kicks if t >= t_k)``.

        Returns
        -------
        CompositePulse
        """
        new_pulses = []
        for p in self.pulses:
            cumulative = sum(
                dphi for t_k, dphi in kicks if p.trigger >= t_k
            )
            if cumulative != 0.0:
                new_pulses.append(p.with_phase_shift(cumulative))
            else:
                new_pulses.append(p._copy())
        return CompositePulse(new_pulses)

    # -- get_kernel (DEPRECATED shim → KernelEstimator, Phase 10.5) --------

    def get_kernel(self, qubit, t_samples=None, *,
                   _suppress_deprecation: bool = False):
        """[DEPRECATED] Use ``KernelEstimator(mode='flux', method='exp').estimate()``.

        This shim forwards to :class:`sqc.reconstruction.KernelEstimator` for
        backward compatibility.  New code should use ``KernelEstimator``
        directly.

        Parameters
        ----------
        qubit : TransmonQubit
            Qubit for simulation.
        t_samples : np.ndarray or None
            Time points to evaluate kernel at.  If None, uses pulse time axis.
        _suppress_deprecation : bool
            Internal flag to silence the deprecation warning when called
            from within the shim itself.

        Returns
        -------
        t_samples : np.ndarray
        kernel : np.ndarray

        Notes
        -----
        Sets ``self.t_samples`` and ``self.kernel`` in-place for backward compat.
        """
        import warnings
        if not _suppress_deprecation:
            warnings.warn(
                "CompositePulse.get_kernel() is deprecated. Use "
                "sqc.reconstruction.KernelEstimator(mode='flux', method='exp', "
                "order=1).estimate(pulse, qubit) instead.",
                DeprecationWarning, stacklevel=2,
            )
        from sqc.reconstruction.kernel import KernelEstimator
        est = KernelEstimator(
            mode='flux', method='exp', order=1,
            deprecation_warn_legacy=False,
        )
        result = est.estimate_full(self, qubit, t_samples=t_samples)
        # Store on self for backward compat
        self.t_samples = result.t_samples
        self.kernel = list(result.k1)  # preserve list type for compat
        return result.t_samples, result.k1
