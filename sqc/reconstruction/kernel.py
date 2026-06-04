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

import math
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

    # -- serialization ----------------------------------------------------

    def save(self, path: str) -> None:
        """Save KernelResult to a ``.npz`` file via :func:`numpy.savez`.

        Parameters
        ----------
        path : str
            File path.  ``.npz`` is appended if not present.
        """
        if not path.endswith('.npz'):
            path += '.npz'
        save_dict = dict(
            t_samples=self.t_samples,
            mode=np.array(self.mode, dtype=object),
            method=np.array(self.method, dtype=object),
            order=np.array(self.order),
            stim_amplitude=np.array(self.stim_amplitude),
            units=np.array(self.units, dtype=object),
        )
        # Store each kernel individually to avoid numpy stacking issues
        for i, k in enumerate(self.kernels):
            save_dict[f'kernels_{i}'] = k
        np.savez(path, **save_dict)

    @classmethod
    def load(cls, path: str) -> 'KernelResult':
        """Load a KernelResult from a ``.npz`` file.

        Parameters
        ----------
        path : str
            File path to a ``.npz`` file saved by :meth:`save`.

        Returns
        -------
        KernelResult
        """
        if not path.endswith('.npz'):
            path += '.npz'
        data = np.load(path, allow_pickle=True)
        order = int(data['order'].item())
        kernels_list = [data[f'kernels_{i}'] for i in range(order)]
        return cls(
            t_samples=data['t_samples'],
            kernels=kernels_list,
            mode=str(data['mode'].item()),
            method=str(data['method'].item()),
            order=order,
            stim_amplitude=float(data['stim_amplitude'].item()),
            units=str(data['units'].item()),
        )


# ═════════════════════════════════════════════════════════════════════════════
# Module-level helpers
# ═════════════════════════════════════════════════════════════════════════════


def _fit_polynomial_no_constant(x, y, N):
    """Fit ``y = a_1*x + a_2*x^2 + ... + a_N*x^N`` (no constant term).

    Uses :func:`numpy.linalg.lstsq` on a Vandermonde matrix without the
    constant column.  For ``N >= 4``, the Vandermonde matrix becomes
    ill-conditioned; the function emits a :class:`RuntimeWarning` if the
    condition number exceeds ``1e8``.

    Parameters
    ----------
    x : array_like
        Independent variable values (1-D).
    y : array_like
        Dependent variable values (1-D).
    N : int
        Highest polynomial order to fit (>= 1).

    Returns
    -------
    coeffs : np.ndarray
        Coefficients ``[a_1, a_2, ..., a_N]``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # Vandermonde: columns are x^N, x^(N-1), ..., x, 1
    # Drop constant column (last) and reverse to get x^1, x^2, ..., x^N
    V = np.vander(x, N=N + 1)[:, :-1][:, ::-1]
    if N >= 4:
        cond = np.linalg.cond(V)
        if cond > 1e8:
            warnings.warn(
                f"Vandermonde condition number {cond:.2e} > 1e8; "
                f"polynomial fit may be unreliable for N={N}.",
                RuntimeWarning,
            )
    coeffs, *_ = np.linalg.lstsq(V, y, rcond=None)
    return coeffs


def _gauss_integral_factor(n: int, sigma: float) -> float:
    """Compute the factor ``c_n = I_n / ε^n / n!`` for an n-th power Gaussian.

    For a Gaussian stimulus ``ε·exp(-t²/(2σ²))``, the integral of its
    n-th power is::

        I_n = ∫ [ε·exp(-t²/(2σ²))]^n dt = ε^n · σ · √(2π/n)

    The Volterra expansion relates the n-th diagonal kernel ``k_n^diag``
    to the polynomial coefficient ``a_n`` (from fitting ``Δp_e`` vs ``ε``)
    via::

        a_n = k_n^diag · I_n / n! = k_n^diag · ε^n · c_n

    Therefore ``k_n^diag = a_n / c_n`` where::

        c_n = σ · √(2π/n) / n!

    Parameters
    ----------
    n : int
        Kernel order (>= 1).
    sigma : float
        Gaussian width (ns).

    Returns
    -------
    float
        The factor ``c_n``.
    """
    return sigma * math.sqrt(2 * math.pi / n) / math.factorial(n)


# ═════════════════════════════════════════════════════════════════════════════
# KernelEstimator
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class KernelEstimator:
    """Estimate the control kernel of a pulse (or composite pulse).

    Perturbs the qubit at each time point of the pulse and measures the
    change in final |e⟩ population p_e.

    Three orthogonal dimensions:
    - ``mode``: physical quantity of the stimulus (``'flux'`` | ``'omega'``)
    - ``method``: simulation strategy (``'sim'`` | ``'exp'``)
    - ``order``: highest Volterra order (always 1 for now)

    **Perturbation modes**:

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

    **Simulation methods**:

    - **exp** (default, backward-compatible):
      Full qubit Hamiltonian simulation via ``qubit.qubit_under_mag()``
      (flux mode) or Virtual Z via σ_z impulse (omega mode).  Requires
      a qubit object.

    - **sim** (new):
      Direct frequency stimulus using the a†a operator — no qubit
      dispersion model, no ``qubit_under_mag()``.  This is the "pure
      theory" kernel.  Only supports ``mode='omega'`` (``mode='flux'``
      is illegal with ``method='sim'`` — use ``mode='omega'`` instead
      and convert via κ manually if needed).

    Parameters
    ----------
    stim_amplitude : float
        Amplitude of the probe stimulus.  Interpretation depends on mode:
        Φ₀ for ``'flux'``, radians for ``'omega'``.
        Default: ``CONFIG.reconstruction.stim_amplitude``.
    stim_width : float
        Width of the Gaussian probe (ns).
        Default: ``CONFIG.reconstruction.stim_width``.
    auto_calibrate : bool
        If True, adjust amplitude per qubit (flux mode only).
        Default: False.
    mode : {'flux', 'omega'}
        Perturbation mode.  Default ``'flux'`` for backward compatibility.
    method : {'sim', 'exp'}
        Simulation strategy.  Default ``'exp'`` for backward compatibility.
        ``'sim'`` uses a pure frequency stimulus (a†a operator) without
        qubit dispersion; only valid with ``mode='omega'``.
        ``'exp'`` uses the full qubit Hamiltonian (requires a qubit).
    virtual_z_impl : {'math', 'hardware'}
        How the Virtual Z perturbation is modelled when ``mode='omega'``.
        ``'math'`` adds a σ_z impulse Hamiltonian term (narrow Gaussian).
        ``'hardware'`` rebuilds the pulse sequence with phase-shifted
        sub-pulses after each probe time.  Default ``'math'``.
    n_levels : int
        Hilbert space dimension when ``method='sim'``.  Default 2.
        Auto-detected from qubit when provided.
    anharmonicity : float
        Anharmonicity α (GHz) used when n_levels >= 3 and
        ``method='sim'``.  Default 0.  Auto-detected from qubit
        when provided.
    kappa : float or None
        Linear dispersion dω/dΦ (GHz/Φ₀).  Only used in sim mode for
        flux→omega conversion; ignored otherwise.  Default None.
    """

    stim_amplitude: float = field(
        default_factory=lambda: CONFIG.reconstruction.stim_amplitude
    )
    stim_width: float = field(
        default_factory=lambda: CONFIG.reconstruction.stim_width
    )
    auto_calibrate: bool = False
    mode: Literal['flux', 'omega'] = 'flux'
    method: Literal['sim', 'exp'] = 'exp'
    virtual_z_impl: Literal['math', 'hardware'] = 'math'

    # -- sim-mode-only parameters -----------------------------------------
    n_levels: int = 2           #: Hilbert space dim when method='sim'
    anharmonicity: float = 0.0  #: Used when n_levels >= 3 and method='sim'
    kappa: float | None = None  #: dω/dΦ, only used in sim mode (for flux→omega conversion)

    # -- higher-order parameters (Phase 10.3) ------------------------------
    order: int = 1              #: Highest Volterra order (1 = linear, backward compat)
    n_amp_samples: int = 5      #: Number of amplitude scan points when order >= 2
    amp_scan_factor: float = 1.0  #: Scan range = [-factor, +factor] * stim_amplitude
    extract_off_diagonal: bool = False  #: Full k_n(t_i, t_j, ...) instead of diagonal

    # -- backward-compat shim (Phase 10.5) ---------------------------------
    deprecation_warn_legacy: bool = True  #: Reserved; suppresses recursion guard
    #  When False, internal create of KernelEstimator from the legacy
    #  Pulse.get_kernel() shim bypasses any future circular-guard logic.

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def estimate(
        self,
        pulse,
        qubit=None,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Estimate kernel for a pulse + qubit system.

        Parameters
        ----------
        pulse : Pulse, CompositePulse, or PulseBase
            Control pulse with .hamiltonian, .t_list, .frame, .omega_d.
        qubit : TransmonQubit or None
            Qubit for simulation.  Required when ``method='exp'`` or
            ``mode='flux'``.  Optional when ``method='sim'`` and
            ``mode='omega'`` (pure theory kernel).
        t_samples : np.ndarray or None
            Time points to evaluate kernel at.  If None, uses
            ``pulse.t_list``.

        Returns
        -------
        t_samples : np.ndarray
            Time points where kernel was evaluated.
        kernel : np.ndarray
            First-order kernel values (backward compat).  Units depend on
            ``self.mode``.  When ``order >= 2``, higher-order kernels are
            computed internally but only ``k1`` is returned; use
            :meth:`estimate_full` to access all orders.
        """
        self._validate_inputs(qubit)

        # -- order=1 fast path (existing behaviour, unchanged) ----------
        if self.order == 1:
            if self.method == 'sim':
                return self._estimate_sim(pulse, qubit, t_samples)
            elif self.mode == 'flux':
                return self._estimate_flux(pulse, qubit, t_samples)
            else:  # omega, exp
                return self._estimate_omega(pulse, qubit, t_samples)

        # -- order >= 2: polynomial-fit amplitude scan ------------------
        t_samples_out, kernels_list = self._extract_kn(pulse, qubit, t_samples)
        return t_samples_out, kernels_list[0]

    def estimate_full(
        self,
        pulse,
        qubit=None,
        t_samples: np.ndarray | None = None,
    ) -> KernelResult:
        """Full estimate returning ``KernelResult`` with all kernel orders.

        Parameters
        ----------
        pulse : Pulse, CompositePulse, or PulseBase
        qubit : TransmonQubit or None
            Optional when ``method='sim'`` and ``mode='omega'``.
        t_samples : np.ndarray or None

        Returns
        -------
        KernelResult
            Result with ``kernels`` list of length ``self.order``.
        """
        if self.order == 1:
            t_samples_out, kernel = self.estimate(pulse, qubit, t_samples)
            kernels_list = [kernel]
        else:
            t_samples_out, kernels_list = self._extract_kn(pulse, qubit, t_samples)

        if self.mode == 'flux':
            units = '1/(Φ₀·ns)'
        else:
            units = 'rad⁻¹'

        return KernelResult(
            t_samples=t_samples_out,
            kernels=kernels_list,
            mode=self.mode,
            method=self.method,
            order=self.order,
            stim_amplitude=self.stim_amplitude,
            units=units,
        )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate_inputs(self, qubit) -> None:
        """Validate (mode, method, qubit, kappa) combination.

        Validation matrix
        -----------------
        ================== ================== ================ ===================
        (mode, method, qubit, kappa)           Behavior
        ================== ================== ================ ===================
        (flux, sim, *, *)                      RAISE ValueError
        (flux, exp, None, *)                   RAISE ValueError
        (flux, exp, qubit, *)                  OK (existing behavior)
        (omega, sim, None, *)                  OK (pure theory, n_levels=2)
        (omega, sim, qubit, *)                 OK (auto-detect params)
        (omega, exp, None, *)                  RAISE ValueError
        (omega, exp, qubit, *)                 OK (existing behavior)
        ================== ================== ================ ===================
        """
        # flux+sim is illegal — flux must go through qubit dispersion
        if self.mode == 'flux' and self.method == 'sim':
            raise ValueError(
                "KernelEstimator(mode='flux', method='sim') is illegal. "
                "flux+sim requires a qubit dispersion model which sim mode "
                "deliberately bypasses. Use mode='omega' with method='sim' "
                "instead, and convert to flux units manually via κ=dω/dΦ "
                "if needed."
            )

        # exp method always requires a qubit (flux needs qubit_under_mag,
        # omega needs qubit Hamiltonian)
        if self.method == 'exp' and qubit is None:
            raise ValueError(
                f"KernelEstimator(mode='{self.mode}', method='exp') "
                f"requires a qubit; got qubit=None."
            )

        # sim method with omega mode and no qubit — OK, pure theory

    # ------------------------------------------------------------------
    # Sim-mode estimation (pure frequency stimulus, no qubit dispersion)
    # ------------------------------------------------------------------

    def _estimate_sim(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Kernel estimation via Heisenberg-picture propagator.

        Uses a single ``sesolve`` to compute the control propagator U(t),
        then evaluates the first-order kernel via the commutator formula:

            k₁(t) = i ⟨0|[W(t), Q]|0⟩

        where W(t) = U†(t)σ_z U(t)/2 and Q = U†(T)M U(T).

        This is the "pure theory" kernel — no perturbation, no mesolve
        loop, no polynomial fitting.  Only valid with ``mode='omega'``
        (flux+sim is rejected by ``_validate_inputs``).

        When ``qubit`` is provided, ``n_levels`` is auto-detected;
        otherwise ``self.n_levels`` is used.
        """
        t_samples_out, kernels_list = self._heisenberg_kernels(
            pulse, qubit, t_samples, order=1,
        )
        return t_samples_out, kernels_list[0]

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
    # Higher-order Volterra extraction (Phase 10.3)
    # ------------------------------------------------------------------

    def _extract_kn(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, list]:
        """Dispatch to the correct higher-order extraction method.

        Returns ``(t_samples, kernels_list)`` where ``kernels_list`` has
        length ``self.order``.
        """
        if self.method == 'sim':
            return self._extract_kn_sim(pulse, qubit, t_samples)
        elif self.mode == 'flux':
            return self._extract_kn_flux(pulse, qubit, t_samples)
        else:  # omega, exp
            return self._extract_kn_omega(pulse, qubit, t_samples)

    # -- sim: higher-order via a†a operator ---------------------------------

    def _heisenberg_kernels(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
        order: int = 1,
    ) -> tuple[np.ndarray, list]:
        """Compute k₁…k_N via Heisenberg-picture propagator (method='sim').

        **Algorithm** (see ``idea/refactor/phase_10_kernel_theory_and_fix.md``):

        1. ``sesolve(qeye(n), t_grid)`` → propagator U(t) at all times.
        2. For each t_j, form W(t_j) = U†(t_j) σ_z U(t_j) / 2.
        3. Evaluate nested commutators:

           - k₁ = i ⟨0| [W, Q] |0⟩
           - k₂ = −⟨0| [W, [W, Q]] |0⟩
           - k₃ = −i ⟨0| [W, [W, [W, Q]]] |0⟩
           - …

        No perturbation, no mesolve loop, no polynomial fitting.
        Works for any ``order`` ≥ 1 at constant cost (one sesolve total).

        Parameters
        ----------
        pulse : Pulse or CompositePulse
        qubit : TransmonQubit or None
            If None, uses ``self.n_levels`` (default 2).
        t_samples : np.ndarray or None
            Time points where kernels are evaluated.  If None, uses the
            deduplicated pulse time grid.
        order : int
            Highest kernel order (1..N).

        Returns
        -------
        (t_samples, [k1, k2, ..., kN]) : tuple
        """
        import time

        t_start = time.time()

        # -- resolve Hilbert space -----------------------------------------
        pulse_ops = pulse.hamiltonian
        if isinstance(pulse_ops, list) and len(pulse_ops) > 0:
            first_op = pulse_ops[0][0] if isinstance(pulse_ops[0], list) else pulse_ops[0]
            pulse_dim = first_op.shape[0] if hasattr(first_op, 'shape') else 2
        else:
            pulse_dim = 2

        if qubit is not None:
            n = getattr(qubit, 'n_levels', self.n_levels)
        else:
            n = self.n_levels

        n = min(n, pulse_dim)
        if n < getattr(qubit, 'n_levels', self.n_levels):
            warnings.warn(
                f"KernelEstimator(method='sim'): pulse operators are "
                f"{pulse_dim}x{pulse_dim}, but qubit has n_levels="
                f"{getattr(qubit, 'n_levels', self.n_levels)}. "
                f"Using n_levels={n} for compatibility."
            )

        # -- build Hamiltonian on clean (deduplicated) time grid -----------
        # pulse.t_list may have duplicates at sub-pulse boundaries.
        # sesolve requires unique time points.
        raw_tlist = np.asarray(pulse.t_list, dtype=float)
        t_grid, unique_idx = np.unique(raw_tlist, return_index=True)

        # RWA Hamiltonian at qubit frequency (or 0 if no qubit)
        omega_ref = (
            qubit.frequency if qubit is not None and hasattr(qubit, 'frequency')
            else 0.0
        )
        if qubit is not None and hasattr(pulse, 'frame') and pulse.frame == 0:
            H0 = QobjEvo(qubit.get_hamiltonian())
        elif n == 2:
            H0 = QobjEvo(0 * qutip.qeye(2))
        else:
            alpha = getattr(qubit, 'anharmonicity', self.anharmonicity) if qubit is not None else self.anharmonicity
            a = qutip.destroy(n)
            H0 = QobjEvo((alpha / 2.0) * (a.dag() * a.dag() * a * a))

        H_pulse = QobjEvo(pulse.hamiltonian_on(t_grid), tlist=t_grid, order=1)
        H_full = H0 + H_pulse

        # -- single sesolve: U(t) as n x n operator -----------------------
        I_op = qutip.qeye(n)
        res = qutip.sesolve(
            H_full, I_op, t_grid,
            options={'max_step': float(CONFIG.awg.dt), 'store_states': True},
        )
        U_list = res.states  # U(t_i, 0) at each t_i

        # -- Heisenberg operators at final time ----------------------------
        sigma_z = qutip.Qobj(np.diag([1.0, -1.0] + [0.0] * (n - 2)))
        M = qutip.basis(n, 1) * qutip.basis(n, 1).dag()   # |1><1|
        U_T = U_list[-1]
        Q = U_T.dag() * M * U_T                            # Heisenberg M

        state0 = qutip.basis(n, 0)

        # -- resolve output time points -----------------------------------
        if t_samples is None:
            t_out = t_grid.copy()
        else:
            t_out = np.asarray(t_samples, dtype=float)

        # -- evaluate k1…kN at each t_j via nested commutators ------------
        kernels = [np.zeros(len(t_out)) for _ in range(order)]

        for i_out, t_j in enumerate(t_out):
            # Find nearest propagator on t_grid
            idx = int(np.argmin(np.abs(t_grid - t_j)))
            U_t = U_list[idx]
            Z_t = U_t.dag() * sigma_z * U_t     # Heisenberg sigma_z
            W = Z_t / 2.0                         # effective perturbation op

            # Compute nested commutators iteratively
            # comm[n] = ad_W^n (Q) = [W, [W, ..., [W, Q]...]] (n nestings)
            comm = Q  # comm[0] = Q
            for n_order in range(1, order + 1):
                comm = W * comm - comm * W        # comm[n] = [W, comm[n-1]]
                val = state0.dag() * comm * state0
                scalar = val[0, 0] if hasattr(val, 'shape') and val.shape == (1, 1) else complex(val)
                # k_n = i^n * <0|ad_W^n(Q)|0>
                kernels[n_order - 1][i_out] = float(np.real(
                    (1j) ** n_order * scalar
                ))

        elapsed = time.time() - t_start
        print(f"Heisenberg kernel (sim, order={order}) "
              f"took {elapsed:.2f} seconds")

        return t_out, kernels

    def _extract_kn_sim(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, list]:
        """Higher-order sim kernels via Heisenberg propagator.

        Thin wrapper around :meth:`_heisenberg_kernels` for order ≥ 2.
        """
        return self._heisenberg_kernels(pulse, qubit, t_samples, order=self.order)

    # -- flux exp: higher-order via FluxSignal scan -------------------------

    def _extract_kn_flux(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, list]:
        """Higher-order flux kernel extraction via polynomial-fit amplitude scan.

        Plan A from handbook §3.5.3.1: direct polynomial scan in flux
        domain.  For each ``t_j``, scan flux amplitude ``ε`` (Φ₀) across
        a symmetric grid, build a ``FluxSignal``, run mesolve, fit
        ``Δp_e`` vs ``ε``, and convert to diagonal flux kernel values.

        Returns ``(t_samples, [k1_flux, ..., kN_flux])``.
        """
        import time

        t_start = time.time()

        pulse_tlist = np.asarray(pulse.t_list)
        n_levels = qubit.n_levels
        psi_e = basis(n_levels, 1)
        N = self.order

        # Build H_0 and H_pulse (same as _estimate_flux)
        if hasattr(pulse, "frame") and pulse.frame == 0:
            H_0 = QobjEvo(qubit.get_hamiltonian())
        else:
            H_0 = QobjEvo(qubit.get_hamiltonian_rwa(
                qubit.frequency if hasattr(qubit, "frequency") else 0.0,
            ))

        H_pulse = QobjEvo(pulse.hamiltonian, tlist=pulse_tlist, order=1)
        H_base = H_0 + H_pulse

        # Baseline p_e
        result_base = mesolve(
            H_base, qubit.state, pulse_tlist, [],
            e_ops=[psi_e * psi_e.dag()],
        )
        p_e_base = result_base.expect[0][-1]

        if t_samples is None:
            t_samples = pulse_tlist.copy()

        frame_val = getattr(pulse, "frame", 1)
        omega_d_val = getattr(pulse, "omega_d", qubit.frequency)

        _DT = CONFIG.awg.dt
        sigma = max(self.stim_width, 2.0 * _DT)

        # Symmetric amplitude grid (flux units Φ₀)
        M = self.n_amp_samples
        eps_grid = np.linspace(
            -self.amp_scan_factor, self.amp_scan_factor, M,
        ) * self._amplitude_for(qubit)

        kernels = [np.zeros(len(t_samples)) for _ in range(N)]

        for i, t_i in enumerate(t_samples):
            delta_p = np.zeros(M)

            for j, eps in enumerate(eps_grid):
                stim = FluxSignal(
                    type=3,
                    t_list=pulse_tlist,
                    amplitude=eps,
                    center=t_i,
                    width=sigma,
                )

                qubit_t = qubit.qubit_under_mag(stim)

                H_stim = QobjEvo(
                    qubit.qubit_under_mag_hamiltonian(
                        qubit_t, stim.t_list, frame_val, omega_d_val,
                    ),
                    tlist=stim.t_list, order=1,
                )

                H_total = H_0 + H_pulse + H_stim
                result_stim = mesolve(
                    H_total, qubit.state, pulse_tlist, [],
                    e_ops=[psi_e * psi_e.dag()],
                )
                delta_p[j] = result_stim.expect[0][-1] - p_e_base

            # Fit polynomial: Δp_e = a₁·ε + a₂·ε² + ... + a_N·ε^N
            coeffs = _fit_polynomial_no_constant(eps_grid, delta_p, N)

            # Convert to diagonal flux kernel values
            for n in range(1, N + 1):
                c_n = _gauss_integral_factor(n, sigma)
                kernels[n - 1][i] = coeffs[n - 1] / c_n

        elapsed = time.time() - t_start
        print(f"Higher-order kernel computation (flux, exp, order={N}) "
              f"took {elapsed:.2f} seconds")

        return t_samples, kernels

    # -- omega exp: higher-order via VZ amplitude scan ----------------------

    def _extract_kn_omega(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, list]:
        """Higher-order omega kernel extraction via Virtual-Z amplitude scan.

        Uses the math VZ approach (narrow Gaussian σ_z impulse) to scan
        ``φ_z`` across a symmetric grid.  Fits ``Δp_e`` vs ``φ_z`` and
        converts polynomial coefficients to diagonal omega kernel values.

        The conversion accounts for the finite width ``σ_t = 2·dt`` of
        the Gaussian approximation to the delta function.  For ``n=1``
        this reduces to the same bilateral-difference formula used in
        the order-1 fast path.

        Returns ``(t_samples, [k1, ..., kN])``.
        """
        import time

        t_start = time.time()

        pulse_tlist = np.asarray(pulse.t_list)
        n_levels = qubit.n_levels
        psi_e = basis(n_levels, 1)
        N = self.order

        # Build H_0 and H_pulse (same as _estimate_omega)
        if hasattr(pulse, "frame") and pulse.frame == 0:
            H_0 = QobjEvo(qubit.get_hamiltonian())
        else:
            H_0 = QobjEvo(qubit.get_hamiltonian_rwa(
                qubit.frequency if hasattr(qubit, "frequency") else 0.0,
            ))

        H_pulse = QobjEvo(pulse.hamiltonian, tlist=pulse_tlist, order=1)
        H_base = H_0 + H_pulse

        # Baseline p_e
        result_base = mesolve(
            H_base, qubit.state, pulse_tlist, [],
            e_ops=[psi_e * psi_e.dag()],
        )
        p_e_base = result_base.expect[0][-1]

        if t_samples is None:
            t_samples = pulse_tlist.copy()

        # σ_z in n-level Fock basis: diag(1, -1, 0, …, 0)
        _DT = CONFIG.awg.dt
        sigma_t = 2.0 * _DT  # narrow Gaussian width for delta approximation
        diag_vals = np.array([1.0, -1.0] + [0.0] * (n_levels - 2))
        sigma_z_op = qutip.Qobj(np.diag(diag_vals))

        # Symmetric phi_z grid (radians)
        M = self.n_amp_samples
        phi_z_grid = np.linspace(
            -self.amp_scan_factor, self.amp_scan_factor, M,
        ) * self.stim_amplitude

        kernels = [np.zeros(len(t_samples)) for _ in range(N)]

        for i, t_j in enumerate(t_samples):
            # Gaussian shape at t_j (shared across all scan amplitudes)
            gauss = np.exp(-0.5 * ((pulse_tlist - t_j) / sigma_t) ** 2)
            norm_factor = 1.0 / (sigma_t * np.sqrt(2 * np.pi))

            delta_p = np.zeros(M)
            for j, phi_z in enumerate(phi_z_grid):
                coeff = phi_z * norm_factor * gauss
                H_vz = QobjEvo(
                    [[0.5 * sigma_z_op, coeff]],
                    tlist=pulse_tlist, order=1,
                )

                result = mesolve(
                    H_base + H_vz, qubit.state, pulse_tlist, [],
                    e_ops=[psi_e * psi_e.dag()],
                )
                delta_p[j] = result.expect[0][-1] - p_e_base

            # Fit polynomial: Δp_e = a₁·φ_z + a₂·φ_z² + ... + a_N·φ_z^N
            coeffs = _fit_polynomial_no_constant(phi_z_grid, delta_p, N)

            # Convert to diagonal omega kernel values.
            # For an ideal delta: k_n^diag = a_n * n!.
            # For the narrow Gaussian approximating the delta, the
            # n-th power integral I_n = ∫[f(t)]^n dt differs from φ_z^n,
            # introducing a factor:
            #   k_n^diag = a_n * n! * σ_t^(n-1) * (2π)^((n-1)/2) * √n
            for n in range(1, N + 1):
                factor = (
                    math.factorial(n)
                    * sigma_t ** (n - 1)
                    * (2 * math.pi) ** ((n - 1) / 2)
                    * math.sqrt(n)
                )
                kernels[n - 1][i] = coeffs[n - 1] * factor

        elapsed = time.time() - t_start
        print(f"Higher-order kernel computation (omega, exp, order={N}) "
              f"took {elapsed:.2f} seconds")

        return t_samples, kernels

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
