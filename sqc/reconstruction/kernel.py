"""sqc.reconstruction.kernel — KernelEstimator + KernelResult.

Unified kernel estimation that replaces three duplicate implementations:
  - src/pulse.py:Pulse.get_kernel
  - src/pulse.py:CompositePulse.get_kernel
  - src/analysis.py:Analysis.get_kernel

Estimates the control kernel by perturbing the qubit at each time point
and measuring the change in final p_e.

Supports two perturbation modes:
  - flux:  narrow Gaussian flux stimulus (legacy; Φ₀ units)
  - omega: Virtual Z phase kick via σ_z impulse (math impl) or
           phase-shifted pulse sequence (hardware impl); dimensionless

See _refactor_plan.md §7.3, §14.1 (debt D2), and phase_10 handbook.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Literal, Optional

import numpy as np
import qutip
from qutip import Qobj, QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.control.pulse import PulseBase, Pulse, CompositePulse


# ═════════════════════════════════════════════════════════════════════════════
# KernelResult dataclass
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class KernelResult:
    """Container for kernel estimation results with metadata.

    Parameters
    ----------
    t_samples : np.ndarray
        Time points where the kernel was evaluated (ns).
    kernels : list of np.ndarray
        Kernel values.  ``kernels[0]`` is the first-order kernel k₁.
        Length equals ``order``.
    mode : str
        Perturbation mode: ``'flux'`` or ``'omega'``.
    method : str
        How the kernel was obtained: ``'exp'`` (experimental/simulation).
    order : int
        Highest kernel order (currently always 1).
    stim_amplitude : float
        Amplitude of the probe stimulus (units depend on mode).
    units : str
        Physical units, e.g. ``'1/(Φ₀·ns)'`` (flux) or ``'rad⁻¹'`` (omega).
    """

    t_samples: np.ndarray
    kernels: list          # list[np.ndarray]; kernels[0] = k₁
    mode: str              # 'flux' | 'omega'
    method: str            # 'exp' (simulation-based)
    order: int             # kernel order (1 for now)
    stim_amplitude: float
    units: str

    @property
    def k1(self) -> np.ndarray:
        """First-order kernel (backward-compat shortcut)."""
        return self.kernels[0]


# ═════════════════════════════════════════════════════════════════════════════
# KernelEstimator
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class KernelEstimator:
    """Estimate the control kernel of a pulse (or composite pulse).

    Perturbs the qubit at each time point of the pulse and measures the
    change in final |e⟩ population p_e.

    Two perturbation modes are supported:

    - **flux** (default, backward-compatible):
      A narrow Gaussian flux stimulus (Φ₀) is added via
      ``qubit_under_mag``.  The kernel has units 1/(Φ₀·ns) and is a
      single-sided finite difference:  k₁ = Δp_e / (flux area).

    - **omega** (Virtual Z):
      A phase kick φ_z (radians) is applied, modelled either as a σ_z
      impulse Hamiltonian (``virtual_z_impl='math'``) or by
      phase-shifting subsequent sub-pulses (``virtual_z_impl='hardware'``).
      The kernel uses a **bilateral** finite difference:
      k₁ = (p₊ − p₋) / (2·φ_z).  Units are rad⁻¹ (dimensionless).

    Parameters
    ----------
    stim_amplitude : float
        Amplitude of the probe stimulus.  Interpretation depends on mode:
        Φ₀ for ``'flux'``, radians for ``'omega'``.
        Default: ``CONFIG.reconstruction.stim_amplitude``.
    stim_width : float
        Width of the Gaussian probe (ns).  Used only in flux mode.
        Default: ``CONFIG.reconstruction.stim_width``.
    auto_calibrate : bool
        If True, adjust amplitude per qubit (flux mode only).
        Default: False.
    mode : {'flux', 'omega'}
        Perturbation mode.  Default ``'flux'`` for backward compatibility.
    virtual_z_impl : {'math', 'hardware'}
        How the Virtual Z perturbation is modelled when ``mode='omega'``.
        ``'math'`` adds a σ_z impulse Hamiltonian term (narrow Gaussian).
        ``'hardware'`` rebuilds the pulse sequence with phase-shifted
        sub-pulses after each probe time.  Default ``'math'``.
    """

    stim_amplitude: float = field(
        default_factory=lambda: CONFIG.reconstruction.stim_amplitude
    )
    stim_width: float = field(
        default_factory=lambda: CONFIG.reconstruction.stim_width
    )
    auto_calibrate: bool = False
    mode: Literal['flux', 'omega'] = 'flux'
    virtual_z_impl: Literal['math', 'hardware'] = 'math'

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def estimate(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Estimate kernel for a pulse + qubit system.

        Parameters
        ----------
        pulse : Pulse, CompositePulse, or PulseBase
            Control pulse with .hamiltonian, .t_list, .frame, .omega_d.
        qubit : TransmonQubit or None
            Qubit for simulation.  Required when ``mode='flux'`` or
            ``mode='omega'``; raises ValueError if None.
        t_samples : np.ndarray or None
            Time points to evaluate kernel at.  If None, uses
            ``pulse.t_list``.

        Returns
        -------
        t_samples : np.ndarray
            Time points where kernel was evaluated.
        kernel : np.ndarray
            Kernel values.  Units depend on ``self.mode``.
        """
        self._validate_inputs(qubit)

        if self.mode == 'flux':
            return self._estimate_flux(pulse, qubit, t_samples)
        else:  # omega
            return self._estimate_omega(pulse, qubit, t_samples)

    def estimate_full(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> KernelResult:
        """Full estimate returning ``KernelResult`` with metadata.

        Parameters
        ----------
        pulse : Pulse, CompositePulse, or PulseBase
        qubit : TransmonQubit
        t_samples : np.ndarray or None

        Returns
        -------
        KernelResult
        """
        t_samples_out, kernel = self.estimate(pulse, qubit, t_samples)

        if self.mode == 'flux':
            units = '1/(Φ₀·ns)'
        else:
            units = 'rad⁻¹'

        return KernelResult(
            t_samples=t_samples_out,
            kernels=[kernel],
            mode=self.mode,
            method='exp',
            order=1,
            stim_amplitude=self.stim_amplitude,
            units=units,
        )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate_inputs(self, qubit) -> None:
        """Validate (mode, qubit) combination.

        Validation matrix
        -----------------
        ===============  ==========
        (mode, qubit)     Behavior
        ===============  ==========
        (flux, None)      RAISE ValueError
        (flux, qubit)     OK
        (omega, None)     RAISE ValueError
        (omega, qubit)    OK
        ===============  ==========
        """
        if qubit is None:
            raise ValueError(
                f"KernelEstimator(mode='{self.mode}') requires a qubit; "
                f"got qubit=None."
            )

    # ------------------------------------------------------------------
    # Flux-mode estimation (legacy logic, preserved exactly)
    # ------------------------------------------------------------------

    def _estimate_flux(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Kernel estimation via narrow Gaussian flux stimulus.

        Single-sided finite difference:  k₁(tⱼ) = (p_e^stim − p_e^base) / area.

        This is the original flux-based algorithm, moved here from the
        monolithic ``estimate()`` method in Phase 10.1.
        """
        import time

        t_start = time.time()

        pulse_tlist = np.asarray(pulse.t_list)
        n_levels = qubit.n_levels
        psi_e = basis(n_levels, 1)

        # Build static H_0
        if hasattr(pulse, "frame") and pulse.frame == 0:
            H_0 = QobjEvo(qubit.get_hamiltonian())
        else:
            H_0 = QobjEvo(qubit.get_hamiltonian_rwa(
                qubit.frequency if hasattr(qubit, "frequency") else 0.0
            ))

        # Build H_pulse
        H_pulse = QobjEvo(pulse.hamiltonian, tlist=pulse_tlist, order=1)
        H_base = H_0 + H_pulse

        # Baseline result
        result_base = mesolve(
            H_base, qubit.state, pulse_tlist, [],
            e_ops=[psi_e * psi_e.dag()],
        )
        p_e_base = result_base.expect[0][-1]

        # Kernel evaluation
        if t_samples is None:
            t_samples = pulse_tlist.copy()

        frame_val = getattr(pulse, "frame", 1)
        omega_d_val = getattr(pulse, "omega_d", qubit.frequency)

        kernel = np.zeros(len(t_samples))
        for i, t_i in enumerate(t_samples):
            # Build narrow Gaussian stimulus at t_i
            stim = FluxSignal(
                type=3,
                t_list=pulse_tlist,
                amplitude=self._amplitude_for(qubit),
                center=t_i,
                width=self.stim_width,
            )
            stim_area = np.trapezoid(stim.signal, stim.t_list)

            # Per-timestep qubit array
            qubit_t = qubit.qubit_under_mag(stim)

            # Build stimulus Hamiltonian (CompositePulse-style, matches legacy)
            H_stim = QobjEvo(
                qubit.qubit_under_mag_hamiltonian(
                    qubit_t,
                    stim.t_list,
                    frame_val,
                    omega_d_val,
                ),
                tlist=stim.t_list,
                order=1,
            )

            H_total = H_0 + H_pulse + H_stim
            result_stim = mesolve(
                H_total, qubit.state, pulse_tlist, [],
                e_ops=[psi_e * psi_e.dag()],
            )
            p_e_stim = result_stim.expect[0][-1]
            kernel[i] = (p_e_stim - p_e_base) / stim_area

        elapsed = time.time() - t_start
        print(f"Kernel computation took {elapsed:.2f} seconds")

        return t_samples, kernel

    # ------------------------------------------------------------------
    # Omega-mode (Virtual Z) estimation
    # ------------------------------------------------------------------

    def _estimate_omega(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Kernel estimation via Virtual Z phase kick.

        Bilateral finite difference:
            k₁(tⱼ) = (p_e⁺ − p_e⁻) / (2·φ_z)

        where p_e⁺ is the |e⟩ population with +φ_z applied at tⱼ and
        p_e⁻ with −φ_z.

        The stimulus amplitude ``self.stim_amplitude`` is interpreted as
        radians (φ_z) in this mode.  The kernel is dimensionless (rad⁻¹).
        """
        import time

        t_start = time.time()

        pulse_tlist = np.asarray(pulse.t_list)
        n_levels = qubit.n_levels
        psi_e = basis(n_levels, 1)

        # Build static H_0
        if hasattr(pulse, "frame") and pulse.frame == 0:
            H_0 = QobjEvo(qubit.get_hamiltonian())
        else:
            H_0 = QobjEvo(qubit.get_hamiltonian_rwa(
                qubit.frequency if hasattr(qubit, "frequency") else 0.0
            ))

        # Build H_pulse (baseline)
        H_pulse = QobjEvo(pulse.hamiltonian, tlist=pulse_tlist, order=1)

        if t_samples is None:
            t_samples = pulse_tlist.copy()

        phi_z = self.stim_amplitude  # interpreted as radians

        if self.virtual_z_impl == 'math':
            kernel = self._omega_vz_math(
                pulse_tlist, n_levels, psi_e,
                H_0, H_pulse, qubit, t_samples, phi_z,
            )
        else:  # hardware
            kernel = self._omega_vz_hardware(
                pulse, pulse_tlist, psi_e,
                H_0, qubit, t_samples, phi_z,
            )

        elapsed = time.time() - t_start
        print(f"Kernel computation (omega, {self.virtual_z_impl}) "
              f"took {elapsed:.2f} seconds")

        return t_samples, kernel

    # -- math VZ: σ_z impulse Hamiltonian ---------------------------------

    def _omega_vz_math(
        self,
        pulse_tlist: np.ndarray,
        n_levels: int,
        psi_e,
        H_0: QobjEvo,
        H_pulse: QobjEvo,
        qubit,
        t_samples: np.ndarray,
        phi_z: float,
    ) -> np.ndarray:
        """Virtual Z via σ_z impulse (narrow Gaussian)."""
        _DT = CONFIG.awg.dt
        sigma_t = 2.0 * _DT

        # σ_z in n-level Fock basis: diag(1, -1, 0, …, 0)
        diag_vals = np.array([1.0, -1.0] + [0.0] * (n_levels - 2))
        sigma_z_op = qutip.Qobj(np.diag(diag_vals))

        # Common normalization
        norm = phi_z / (sigma_t * np.sqrt(2 * np.pi))

        kernel = np.zeros(len(t_samples))
        H_base = H_0 + H_pulse

        for i, t_j in enumerate(t_samples):
            # Gaussian impulse at t_j
            gauss = np.exp(-0.5 * ((pulse_tlist - t_j) / sigma_t) ** 2)
            coeff_plus = +norm * gauss
            coeff_minus = -norm * gauss

            H_vz_plus = QobjEvo(
                [[0.5 * sigma_z_op, coeff_plus]],
                tlist=pulse_tlist, order=1,
            )
            H_vz_minus = QobjEvo(
                [[0.5 * sigma_z_op, coeff_minus]],
                tlist=pulse_tlist, order=1,
            )

            result_plus = mesolve(
                H_base + H_vz_plus, qubit.state, pulse_tlist, [],
                e_ops=[psi_e * psi_e.dag()],
            )
            result_minus = mesolve(
                H_base + H_vz_minus, qubit.state, pulse_tlist, [],
                e_ops=[psi_e * psi_e.dag()],
            )

            p_plus = result_plus.expect[0][-1]
            p_minus = result_minus.expect[0][-1]

            kernel[i] = (p_plus - p_minus) / (2.0 * phi_z)

        return kernel

    # -- hardware VZ: phase-shifted pulse sequence ------------------------

    def _omega_vz_hardware(
        self,
        pulse,
        pulse_tlist: np.ndarray,
        psi_e,
        H_0: QobjEvo,
        qubit,
        t_samples: np.ndarray,
        phi_z: float,
    ) -> np.ndarray:
        """Virtual Z via phase-shifted sub-pulses (hardware-realistic).

        For each probe time tⱼ, rebuilds the pulse sequence with
        sub-pulses whose trigger > tⱼ phase-shifted by ±φ_z.

        Falls back to the math implementation when the pulse is a single
        ``Pulse`` (no sub-pulse splitting possible).
        """
        # If not a CompositePulse, fall back to math with a warning.
        if not isinstance(pulse, CompositePulse):
            warnings.warn(
                "virtual_z_impl='hardware' on a single Pulse falls back "
                "to math implementation (σ_z impulse)."
            )
            n_levels = qubit.n_levels
            # Rebuild H_pulse for consistency
            H_pulse_fb = QobjEvo(pulse.hamiltonian, tlist=pulse_tlist, order=1)
            return self._omega_vz_math(
                pulse_tlist, n_levels, psi_e,
                H_0, H_pulse_fb, qubit, t_samples, phi_z,
            )

        # CompositePulse: phase-shift sub-pulses after each t_j
        kernel = np.zeros(len(t_samples))
        for i, t_j in enumerate(t_samples):
            pulse_plus = pulse.with_phase_shift(+phi_z, from_time=t_j)
            pulse_minus = pulse.with_phase_shift(-phi_z, from_time=t_j)

            H_pulse_plus = QobjEvo(
                pulse_plus.hamiltonian, tlist=pulse_tlist, order=1,
            )
            H_pulse_minus = QobjEvo(
                pulse_minus.hamiltonian, tlist=pulse_tlist, order=1,
            )

            result_plus = mesolve(
                H_0 + H_pulse_plus, qubit.state, pulse_tlist, [],
                e_ops=[psi_e * psi_e.dag()],
            )
            result_minus = mesolve(
                H_0 + H_pulse_minus, qubit.state, pulse_tlist, [],
                e_ops=[psi_e * psi_e.dag()],
            )

            p_plus = result_plus.expect[0][-1]
            p_minus = result_minus.expect[0][-1]

            kernel[i] = (p_plus - p_minus) / (2.0 * phi_z)

        return kernel

    # ------------------------------------------------------------------
    # Stimulus amplitude helpers
    # ------------------------------------------------------------------

    def _amplitude_for(self, qubit) -> float:
        """Determine stimulus amplitude for the given qubit.

        Parameters
        ----------
        qubit : TransmonQubit

        Returns
        -------
        float
            Stimulus amplitude.
        """
        if not self.auto_calibrate:
            return self.stim_amplitude
        # Auto-calibration: pick amplitude such that the stimulus produces
        # a frequency shift ~ 1% of qubit anharmonicity.
        try:
            kappa = qubit.frequency_sensitivity(
                qubit.flux if hasattr(qubit, "flux") else 0.0
            )
        except Exception:
            kappa = 0.0
        if kappa == 0.0:
            return self.stim_amplitude
        target_shift = 0.01 * abs(qubit.anharmonicity)
        return target_shift / abs(kappa)
