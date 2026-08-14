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
    off_diagonal: bool = False  # True -> kernels[n-1] is an n-D array
    #  (shape (M,)*n) holding the full k_n(t_i, t_j, ...); only diagonal
    #  kernels (all 1-D) are accepted by Wiener/Hammerstein reconstruction.

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
            off_diagonal=np.array(self.off_diagonal),
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
        off_diag = bool(data['off_diagonal'].item()) if 'off_diagonal' in data else False
        return cls(
            t_samples=data['t_samples'],
            kernels=kernels_list,
            mode=str(data['mode'].item()),
            method=str(data['method'].item()),
            order=order,
            stim_amplitude=float(data['stim_amplitude'].item()),
            units=str(data['units'].item()),
            off_diagonal=off_diag,
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

    # -- probe-width / Richardson parameters (Phase 12) --------------------
    probe_sigma_t: float | None = None  #: VZ Gaussian width (ns) for omega exp;
    #  None -> 2.0 * CONFIG.awg.dt (legacy default, numerically unchanged).
    #  Smaller sigma_t reduces the off-diagonal smearing bias of the exp
    #  high-order diagonal (k₃ is ~15% low at 2·dt) but sigma_t < dt is
    #  under-resolved on the time grid and blows up — keep sigma_t >= dt.
    richardson: bool = False    #: order>=2 exp: sample several sigma_t and
    #  extrapolate sigma_t -> 0 per (order, time-point) to remove the smear bias.
    richardson_sigmas: tuple[float, ...] | None = None  #: sigma_t values as
    #  multiples of dt; None -> (2.0, 1.5, 1.0) (all >= dt to avoid grid blow-up).

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
        self._validate_inputs(qubit)
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
            off_diagonal=bool(self.extract_off_diagonal and self.order >= 2),
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

        # off-diagonal extraction — sim route is exact (Heisenberg n-D kernels);
        # exp route uses mixed-partial FD (generalizing the diagonal 5-pt stencil).
        # order > 3 is rejected regardless of method (M^n blow-up).
        if self.extract_off_diagonal and self.order > 3:
            raise ValueError(
                f"extract_off_diagonal supports order <= 3 (got {self.order}); "
                f"k_n for n>=4 is an M^n tensor (combinatorial blow-up)."
            )

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
        sigma_t = self._sigma_t()

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
        elif self.extract_off_diagonal:
            # exp mixed-partial FD → full n-D tensors
            return self._extract_kn_offdiag_exp(pulse, qubit, t_samples)
        elif self.mode == 'flux':
            return self._extract_kn_flux(pulse, qubit, t_samples)
        else:  # omega, exp, diagonal
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
        When ``self.extract_off_diagonal`` is set, dispatches to
        :meth:`_heisenberg_kernels_offdiag` which returns full n-D kernels
        instead of diagonal-only.
        """
        if self.extract_off_diagonal:
            return self._heisenberg_kernels_offdiag(
                pulse, qubit, t_samples, order=self.order,
            )
        return self._heisenberg_kernels(pulse, qubit, t_samples, order=self.order)

    # -- sim: full off-diagonal kernels via Heisenberg propagator -----------

    def _heisenberg_kernels_offdiag(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
        order: int = 2,
    ) -> tuple[np.ndarray, list]:
        """Full off-diagonal Volterra kernels k₁…k_N via Heisenberg sim.

        Unlike :meth:`_heisenberg_kernels` (diagonal only), this returns
        ``kernels[n-1]`` as an **n-D** ndarray of shape ``(M,)*n`` holding
        the fully-symmetric Volterra kernel evaluated on every combination
        of probe times (``M = len(t_samples)``).

        **Algorithm** (theory §2.3–2.4; verified in
        ``kernel/verify_full_kernel.py``).  With ``W(t) = U†(t) σ_z U(t)/2``
        and ``Q = U†(T) M U(T)``:

        - k₁(t)            = i ⟨0|[W(t), Q]|0⟩
        - k₂(t_>, t_<)     = −⟨0|[W(t_<), [W(t_>), Q]]|0⟩
        - k₃(t₁≥t₂≥t₃)     = −i ⟨0|[W(t₃), [W(t₂), [W(t₁), Q]]]|0⟩

        time-ordered, then symmetrized over all permutations.  One
        ``sesolve`` total; the n-D fill is pure numpy matrix algebra.

        Cost scales as ``Mⁿ`` commutator evaluations — only ``order ≤ 3``
        is supported (order ≥ 4 raises).  ``mode='omega'`` only; ``exp``
        is rejected upstream in :meth:`_validate_inputs`.

        Returns ``(t_samples, [k1, k2, ..., kN])`` with k1 1-D and
        k_n (n≥2) n-D.
        """
        import time

        if order < 2:
            # Off-diagonal is meaningless for order 1; fall back to diagonal.
            return self._heisenberg_kernels(pulse, qubit, t_samples, order=1)
        if order > 3:
            raise ValueError(
                f"extract_off_diagonal supports order <= 3 (got {order}); "
                f"k_n for n>=4 is an M^n tensor (combinatorial blow-up)."
            )

        t_start = time.time()

        # -- resolve Hilbert space (mirror _heisenberg_kernels) ------------
        pulse_ops = pulse.hamiltonian
        if isinstance(pulse_ops, list) and len(pulse_ops) > 0:
            first_op = pulse_ops[0][0] if isinstance(pulse_ops[0], list) else pulse_ops[0]
            pulse_dim = first_op.shape[0] if hasattr(first_op, 'shape') else 2
        else:
            pulse_dim = 2
        n = getattr(qubit, 'n_levels', self.n_levels) if qubit is not None else self.n_levels
        n = min(n, pulse_dim)

        # -- Hamiltonian on deduplicated grid + single sesolve -------------
        raw_tlist = np.asarray(pulse.t_list, dtype=float)
        t_grid, _ = np.unique(raw_tlist, return_index=True)

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

        res = qutip.sesolve(
            H_full, qutip.qeye(n), t_grid,
            options={'max_step': float(CONFIG.awg.dt), 'store_states': True},
        )
        U_list = res.states

        sigma_z = np.diag([1.0, -1.0] + [0.0] * (n - 2)).astype(complex)
        M_proj = np.zeros((n, n), dtype=complex)
        M_proj[1, 1] = 1.0
        U_T = np.asarray(U_list[-1].full())
        Q = U_T.conj().T @ M_proj @ U_T
        state0 = np.zeros(n, dtype=complex); state0[0] = 1.0

        # -- output time points + W(t) at each (numpy, for speed) ----------
        if t_samples is None:
            t_out = t_grid.copy()
        else:
            t_out = np.asarray(t_samples, dtype=float)
        M = len(t_out)

        W_out = np.empty((M, n, n), dtype=complex)
        for i, t_j in enumerate(t_out):
            idx = int(np.argmin(np.abs(t_grid - t_j)))
            U_t = np.asarray(U_list[idx].full())
            W_out[i] = (U_t.conj().T @ sigma_z @ U_t) / 2.0

        def comm(A, B):
            return A @ B - B @ A

        def expect0(op):
            return complex(state0.conj() @ op @ state0)

        if M ** order > 2_000_000:
            warnings.warn(
                f"_heisenberg_kernels_offdiag: M^order = {M}^{order} = "
                f"{M**order} commutator evals — this may be slow. "
                f"Reduce t_samples for order={order}."
            )

        kernels: list = []

        # -- k1 (1-D) ------------------------------------------------------
        k1 = np.zeros(M)
        for i in range(M):
            k1[i] = float(np.real(1j * expect0(comm(W_out[i], Q))))
        kernels.append(k1)

        # -- k2 (M, M) -----------------------------------------------------
        if order >= 2:
            k2 = np.zeros((M, M))
            for i in range(M):
                for j in range(M):
                    # t_> = later time, t_< = earlier time
                    if t_out[i] >= t_out[j]:
                        Wg, Wl = W_out[i], W_out[j]
                    else:
                        Wg, Wl = W_out[j], W_out[i]
                    c2 = comm(Wl, comm(Wg, Q))
                    k2[i, j] = float(np.real(-expect0(c2)))
            kernels.append(k2)

        # -- k3 (M, M, M) --------------------------------------------------
        if order >= 3:
            k3 = np.zeros((M, M, M))
            for i in range(M):
                for j in range(M):
                    for l in range(M):
                        # time-order descending: W(t_(1)>=t_(2)>=t_(3))
                        order_idx = sorted((i, j, l), key=lambda x: -t_out[x])
                        a, b, c = order_idx
                        cc = comm(W_out[a], Q)        # innermost: largest time
                        cc = comm(W_out[b], cc)
                        cc = comm(W_out[c], cc)        # outermost: smallest time
                        k3[i, j, l] = float(np.real(-1j * expect0(cc)))
            kernels.append(k3)

        elapsed = time.time() - t_start
        print(f"Heisenberg off-diagonal kernel (sim, order={order}, M={M}) "
              f"took {elapsed:.2f} seconds")
        return t_out, kernels

    # -- flux exp: higher-order via FluxSignal scan -------------------------

    def _extract_kn_flux(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, list]:
        """Higher-order flux kernel via 5-point finite-difference stencil.

        Replaces the old polynomial-fit method (which suffered from
        Vandermonde ill-conditioning).  Samples ``p_e`` at exactly 5
        flux amplitudes ``{-2h, -h, 0, +h, +2h}`` and applies fixed-
        coefficient stencils for k₁…k₃.

        Flux perturbation uses ``FluxSignal`` + ``qubit_under_mag()``,
        identical to the legacy order-1 path.

        Returns ``(t_samples, [k1, ..., kN])``.
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

        if t_samples is None:
            t_samples = pulse_tlist.copy()

        frame_val = getattr(pulse, "frame", 1)
        omega_d_val = getattr(pulse, "omega_d", qubit.frequency)
        h = self._amplitude_for(qubit) * self.amp_scan_factor
        # Clamp h to avoid huge flux perturbations
        h = min(h, 0.05)

        kernels = [np.zeros(len(t_samples)) for _ in range(N)]

        for i, t_i in enumerate(t_samples):
            # Evaluate p_e at 5 flux amplitudes
            p_vals = np.zeros(5)
            for j, mult in enumerate([-2, -1, 0, 1, 2]):
                eps = mult * h
                stim = FluxSignal(
                    type=3, t_list=pulse_tlist,
                    amplitude=eps, center=t_i,
                    width=self.stim_width,
                )
                qubit_t = qubit.qubit_under_mag(stim)
                H_stim = QobjEvo(
                    qubit.qubit_under_mag_hamiltonian(
                        qubit_t, stim.t_list, frame_val, omega_d_val,
                    ),
                    tlist=stim.t_list, order=1,
                )
                result = mesolve(
                    H_base + H_stim, qubit.state, pulse_tlist, [],
                    e_ops=[psi_e * psi_e.dag()],
                )
                p_vals[j] = result.expect[0][-1]

            fm2, fm1, f0, fp1, fp2 = p_vals

            # 5-point FD stencils
            if N >= 1:
                kernels[0][i] = (fm2 - 8*fm1 + 8*fp1 - fp2) / (12 * h)
            if N >= 2:
                kernels[1][i] = (-fm2 + 16*fm1 - 30*f0 + 16*fp1 - fp2) / (12 * h**2)
            if N >= 3:
                kernels[2][i] = (-fm2 + 2*fm1 - 2*fp1 + fp2) / (2 * h**3)

        elapsed = time.time() - t_start
        print(f"FD kernel (flux exp, order={N}) took {elapsed:.2f} seconds")
        return t_samples, kernels

    # -- omega exp: FD stencil via Virtual Z ----------------------------------

    def _extract_kn_omega(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, list]:
        """Higher-order omega kernel via 5-point finite-difference stencil.

        Replaces the old polynomial-fit method.  Samples ``p_e`` at
        exactly 5 VZ amplitudes ``{-2h, -h, 0, +h, +2h}`` (radians)
        and applies fixed-coefficient stencils.

        The VZ Gaussian probe width σ_t is resolved by :meth:`_sigma_t`
        (default ``2·dt``).  A finite σ_t convolves the off-diagonal
        kernel and biases the high-order diagonal (k₃ ~15% low at 2·dt).
        When ``self.richardson`` is set and ``order >= 2``, the stencil is
        run at several σ_t (``richardson_sigmas`` × dt) and extrapolated
        to σ_t → 0 per (order, time-point) to remove that bias; ``k₁`` is
        taken from the smallest σ_t (it is already accurate to O(σ_t²)).

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

        if t_samples is None:
            t_samples = pulse_tlist.copy()

        diag_vals = np.array([1.0, -1.0] + [0.0] * (n_levels - 2))
        sigma_z_op = qutip.Qobj(np.diag(diag_vals))
        h = self.stim_amplitude * self.amp_scan_factor

        # High-order stencils divide by h^n (n up to 3); the default
        # mesolve tolerance (~1e-8) then dominates k₃ as noise/h³.  Use a
        # tight tolerance so the FD truncation error — not integrator
        # noise — sets the accuracy (see phase_10_kernel_theory_and_fix §5.3).
        _fd_opts = {
            "atol": 1e-12, "rtol": 1e-10, "max_step": float(CONFIG.awg.dt),
        }

        def _fd_kernels_at(sigma_t: float) -> list:
            """5-point FD stencil for k₁…k_N at a fixed probe width σ_t."""
            norm_f = 1.0 / (sigma_t * np.sqrt(2 * np.pi))
            ks = [np.zeros(len(t_samples)) for _ in range(N)]
            for i, t_j in enumerate(t_samples):
                gauss = np.exp(-0.5 * ((pulse_tlist - t_j) / sigma_t) ** 2)
                p_vals = np.zeros(5)
                for j, mult in enumerate([-2, -1, 0, 1, 2]):
                    coeff = (mult * h) * norm_f * gauss
                    H_vz = QobjEvo(
                        [[0.5 * sigma_z_op, coeff]],
                        tlist=pulse_tlist, order=1,
                    )
                    result = mesolve(
                        H_base + H_vz, qubit.state, pulse_tlist, [],
                        e_ops=[psi_e * psi_e.dag()],
                        options=_fd_opts,
                    )
                    p_vals[j] = result.expect[0][-1]
                fm2, fm1, f0, fp1, fp2 = p_vals
                if N >= 1:
                    ks[0][i] = (fm2 - 8*fm1 + 8*fp1 - fp2) / (12 * h)
                if N >= 2:
                    ks[1][i] = (-fm2 + 16*fm1 - 30*f0 + 16*fp1 - fp2) / (12 * h**2)
                if N >= 3:
                    ks[2][i] = (-fm2 + 2*fm1 - 2*fp1 + fp2) / (2 * h**3)
            return ks

        def _fd_kernels_at_hw() -> list:
            """5-point FD stencil via hardware VZ (phase-shifted sub-pulses).

            No sigma_t dependence — phase kicks are instantaneous.
            """
            ks = [np.zeros(len(t_samples)) for _ in range(N)]
            for i, t_j in enumerate(t_samples):
                p_vals = np.zeros(5)
                for j, mult in enumerate([-2, -1, 0, 1, 2]):
                    pulse_kicked = pulse.with_phase_kicks(
                        [(t_j, mult * h)]
                    )
                    H_pulse_kicked = QobjEvo(
                        pulse_kicked.hamiltonian_on(pulse_tlist),
                        tlist=pulse_tlist, order=1,
                    )
                    result = mesolve(
                        H_0 + H_pulse_kicked, qubit.state, pulse_tlist,
                        [], e_ops=[psi_e * psi_e.dag()],
                        options=_fd_opts,
                    )
                    p_vals[j] = result.expect[0][-1]
                fm2, fm1, f0, fp1, fp2 = p_vals
                if N >= 1:
                    ks[0][i] = (fm2 - 8*fm1 + 8*fp1 - fp2) / (12 * h)
                if N >= 2:
                    ks[1][i] = (-fm2 + 16*fm1 - 30*f0 + 16*fp1 - fp2) / (12 * h**2)
                if N >= 3:
                    ks[2][i] = (-fm2 + 2*fm1 - 2*fp1 + fp2) / (2 * h**3)
            return ks

        # -- dispatch: math vs hardware VZ ---------------------------------
        if self.virtual_z_impl == 'hardware':
            if not isinstance(pulse, CompositePulse):
                warnings.warn(
                    "virtual_z_impl='hardware' on a single Pulse falls back "
                    "to math implementation (σ_z impulse) for higher-order "
                    "diagonal kernels."
                )
                if not (self.richardson and N >= 2):
                    kernels = _fd_kernels_at(self._sigma_t())
                else:
                    sigmas = self._richardson_sigma_values()
                    per_sigma = [_fd_kernels_at(st) for st in sigmas]
                    kernels = self._richardson_extrapolate(sigmas, per_sigma)
                    i_min = int(np.argmin(sigmas))
                    kernels[0] = per_sigma[i_min][0]
            else:
                if self.richardson and N >= 2:
                    warnings.warn(
                        "Richardson extrapolation is not applicable to "
                        "virtual_z_impl='hardware' (phase kicks are "
                        "instantaneous, no probe width to extrapolate); "
                        "ignoring richardson=True.",
                        RuntimeWarning,
                    )
                kernels = _fd_kernels_at_hw()
        elif not (self.richardson and N >= 2):
            # Single probe width (legacy / parameterized).
            kernels = _fd_kernels_at(self._sigma_t())
        else:
            # Multi-σ_t Richardson extrapolation toward σ_t -> 0.
            sigmas = self._richardson_sigma_values()
            per_sigma = [_fd_kernels_at(st) for st in sigmas]   # list over σ_t
            kernels = self._richardson_extrapolate(sigmas, per_sigma)
            # k₁ kept from the smallest σ_t (already O(σ_t²) accurate).
            i_min = int(np.argmin(sigmas))
            kernels[0] = per_sigma[i_min][0]

        elapsed = time.time() - t_start
        print(f"FD kernel (omega exp, order={N}, "
              f"richardson={self.richardson and N >= 2}) "
              f"took {elapsed:.2f} seconds")
        return t_samples, kernels

    # -- exp: full off-diagonal kernels via mixed-partial FD --------------

    def _extract_kn_offdiag_exp(
        self,
        pulse,
        qubit,
        t_samples: np.ndarray | None = None,
    ) -> tuple[np.ndarray, list]:
        """Full off-diagonal Volterra kernels k₁…k_N via mixed-partial FD.

        Generalises the diagonal 5-point FD stencil to n simultaneous kicks
        at distinct probe times, recovering the full n-D symmetric tensor
        from measurable ``p_e`` alone — the hardware-measurable analog of
        :meth:`_heisenberg_kernels_offdiag`.

        **Method** (order 2, omega mode):

        - k₁(i): 5-point single-axis stencil (existing diagonal).
        - k₂(i,i): 5-point 2nd derivative on the diagonal.
        - k₂(i,j), i≠j: 4-corner mixed central difference
          ``[p(+,+) − p(+,−) − p(−,+) + p(−,−)] / 4h²``.
        - Mirror to ``(j,i)`` for full symmetry.

        For order 3 the tensor is computed over the sorted-index wedge
        ``i ≤ j ≤ l`` with three cases per the degeneracy of the index
        triple, then mirrored to all 6 permutations.

        A per-session ``p_e`` cache keyed by ``((idx, mult), …)`` reuses
        shared corners (base ``p(0)``, single-kick, multi-kick) across
        many tensor entries, keeping the cost at ~O(M²···M³) mesolve calls
        rather than a naive O(Mⁿ).

        ``mode='flux'`` builds summed Gaussian ``FluxSignal`` objects and
        passes them through ``qubit_under_mag``, correctly folding the
        pointwise κ-nonlinearity into the kernel (genuine flux kernel,
        distinct from the sim omega kernel).

        Parameters
        ----------
        pulse : Pulse or CompositePulse
        qubit : TransmonQubit
        t_samples : np.ndarray or None

        Returns
        -------
        (t_samples, [k1, (k2), (k3)]) : tuple
            k1 is 1-D shape (M,); k2 is (M,M); k3 is (M,M,M).
        """
        import time
        t_start_full = time.time()

        pulse_tlist = np.asarray(pulse.t_list)
        n_levels = qubit.n_levels
        psi_e = basis(n_levels, 1)
        N = min(self.order, 3)  # order > 3 already rejected in _validate_inputs

        # -- build H_base once ---------------------------------------------
        if hasattr(pulse, "frame") and pulse.frame == 0:
            H_0 = QobjEvo(qubit.get_hamiltonian())
        else:
            H_0 = QobjEvo(qubit.get_hamiltonian_rwa(
                qubit.frequency if hasattr(qubit, "frequency") else 0.0,
            ))
        H_pulse = QobjEvo(pulse.hamiltonian, tlist=pulse_tlist, order=1)
        H_base = H_0 + H_pulse

        if t_samples is None:
            t_samples = pulse_tlist.copy()
        M = len(t_samples)

        # -- cost warning (mirror sim offdiag) -----------------------------
        _cost_thresh = 2_000_000
        if M ** N > _cost_thresh:
            warnings.warn(
                f"_extract_kn_offdiag_exp: M^order = {M}^{N} = "
                f"{M**N} combinations — this will be slow. "
                f"Reduce t_samples for order={N}."
            )

        # -- stencil engine (pure function of h, pe) ------------------------
        def _offdiag_stencils(M_val, h_val, _pe):
            """Return ``[k1, (k2), (k3)]`` using the given FD step and
            ``p_e`` evaluator.  All index loops and permutation fills live
            here so the same code serves omega (any σ_t) and flux modes."""
            p_base = _pe()

            # k₁ (M,)
            k1 = np.zeros(M_val, dtype=float)
            for i in range(M_val):
                fm2 = _pe((i, -2))
                fm1 = _pe((i, -1))
                fp1 = _pe((i, +1))
                fp2 = _pe((i, +2))
                k1[i] = (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h_val)
            ks: list = [k1]

            # k₂ (M, M)
            if N >= 2:
                k2 = np.zeros((M_val, M_val), dtype=float)
                h2 = h_val * h_val
                for i in range(M_val):
                    fm2 = _pe((i, -2))
                    fm1 = _pe((i, -1))
                    fp1 = _pe((i, +1))
                    fp2 = _pe((i, +2))
                    k2[i, i] = (
                        -fm2 + 16.0 * fm1 - 30.0 * p_base
                        + 16.0 * fp1 - fp2
                    ) / (12.0 * h2)
                    for j in range(i + 1, M_val):
                        pp = _pe((i, +1), (j, +1))
                        pm = _pe((i, +1), (j, -1))
                        mp = _pe((i, -1), (j, +1))
                        mm = _pe((i, -1), (j, -1))
                        val = (pp - pm - mp + mm) / (4.0 * h2)
                        k2[i, j] = val
                        k2[j, i] = val
                ks.append(k2)

            # k₃ (M, M, M)
            if N >= 3:
                k3 = np.zeros((M_val, M_val, M_val), dtype=float)
                denom_2h3 = 2.0 * h_val * h_val * h_val
                denom_8h3 = 8.0 * h_val * h_val * h_val

                def _k3_diag(idx):
                    fm2 = _pe((idx, -2))
                    fm1 = _pe((idx, -1))
                    fp1 = _pe((idx, +1))
                    fp2 = _pe((idx, +2))
                    return (-fm2 + 2.0 * fm1 - 2.0 * fp1 + fp2) / denom_2h3

                def _k3_two_equal(idx_double, idx_single):
                    s = 0.0
                    s -= _pe((idx_single, -1), (idx_double, -1))
                    s += _pe((idx_single, -1), (idx_double, 0)) * 2.0
                    s -= _pe((idx_single, -1), (idx_double, +1))
                    s += _pe((idx_single, +1), (idx_double, -1))
                    s -= _pe((idx_single, +1), (idx_double, 0)) * 2.0
                    s += _pe((idx_single, +1), (idx_double, +1))
                    return s / denom_2h3

                def _k3_distinct(a, b, c):
                    s = 0.0
                    for si in [-1, +1]:
                        for sj in [-1, +1]:
                            for sk in [-1, +1]:
                                s += si * sj * sk * _pe(
                                    (a, si), (b, sj), (c, sk),
                                )
                    return s / denom_8h3

                for i in range(M_val):
                    k3[i, i, i] = _k3_diag(i)
                    for j in range(i, M_val):
                        for l in range(j, M_val):
                            if i == j == l:
                                continue
                            if i == j and j < l:
                                val = _k3_two_equal(i, l)
                            elif i < j and j == l:
                                val = _k3_two_equal(j, i)
                            elif i < j and j < l:
                                val = _k3_distinct(i, j, l)
                            else:
                                continue
                            for a, b, c in [
                                (i, j, l), (i, l, j), (j, i, l),
                                (j, l, i), (l, i, j), (l, j, i),
                            ]:
                                k3[a, b, c] = val
                ks.append(k3)

            return ks

        # -- FD parameters (shared) ------------------------------------------
        _fd_opts = {
            "atol": 1e-12, "rtol": 1e-10,
            "max_step": float(CONFIG.awg.dt),
        }

        if self.mode == 'flux':
            h_flux = self._amplitude_for(qubit) * self.amp_scan_factor
            h_flux = min(h_flux, 0.05)
            frame_val = getattr(pulse, "frame", 1)
            omega_d_val = getattr(pulse, "omega_d", qubit.frequency)
        else:
            # omega: h is same across all sigma_t
            h_omega = self.stim_amplitude * self.amp_scan_factor

        # -- omega inner: build + evaluate at a given sigma_t ---------------
        def _compute_offdiag_omega(sigma_t: float):
            """Return ``[k1, (k2), (k3)]`` at a fixed VZ probe width."""
            norm_f = 1.0 / (sigma_t * np.sqrt(2 * np.pi))
            diag_vals_z = np.array([1.0, -1.0] + [0.0] * (n_levels - 2))
            sigma_z_op = qutip.Qobj(np.diag(diag_vals_z))
            h = h_omega
            _pe_cache: dict[tuple, float] = {}

            def _pe_kicks(kicks_tuple: tuple) -> float:
                if kicks_tuple in _pe_cache:
                    return _pe_cache[kicks_tuple]
                coeff_sum = np.zeros(len(pulse_tlist), dtype=float)
                for idx, mult in kicks_tuple:
                    t_p = t_samples[idx]
                    gauss = np.exp(
                        -0.5 * ((pulse_tlist - t_p) / sigma_t) ** 2
                    )
                    coeff_sum += (mult * h) * norm_f * gauss
                H_vz = QobjEvo(
                    [[0.5 * sigma_z_op, coeff_sum]],
                    tlist=pulse_tlist, order=1,
                )
                result = mesolve(
                    H_base + H_vz, qubit.state, pulse_tlist, [],
                    e_ops=[psi_e * psi_e.dag()],
                    options=_fd_opts,
                )
                _pe_cache[kicks_tuple] = float(result.expect[0][-1])
                return _pe_cache[kicks_tuple]

            def _pe(*specs):
                return _pe_kicks(tuple(sorted(specs, key=lambda x: x[0])))

            return _offdiag_stencils(M, h, _pe)

        # -- omega hardware inner: build + evaluate (no sigma_t knob) ------
        def _compute_offdiag_omega_hardware():
            """Return ``[k1, (k2), (k3)]`` via hardware VZ (phase-shifted
            sub-pulses).  No sigma_t — phase kicks are instantaneous."""
            h = h_omega
            _pe_cache: dict[tuple, float] = {}

            def _pe_kicks(kicks_tuple: tuple) -> float:
                if kicks_tuple in _pe_cache:
                    return _pe_cache[kicks_tuple]
                kicks_list = [
                    (t_samples[idx], mult * h) for idx, mult in kicks_tuple
                ]
                kicks_list.sort(key=lambda x: x[0])
                pulse_kicked = pulse.with_phase_kicks(kicks_list)
                H_pulse_kicked = QobjEvo(
                    pulse_kicked.hamiltonian_on(pulse_tlist),
                    tlist=pulse_tlist, order=1,
                )
                result = mesolve(
                    H_0 + H_pulse_kicked, qubit.state, pulse_tlist, [],
                    e_ops=[psi_e * psi_e.dag()],
                    options=_fd_opts,
                )
                _pe_cache[kicks_tuple] = float(result.expect[0][-1])
                return _pe_cache[kicks_tuple]

            def _pe(*specs):
                return _pe_kicks(tuple(sorted(specs, key=lambda x: x[0])))

            return _offdiag_stencils(M, h, _pe)

        # -- flux inner: build + evaluate (no sigma_t knob) -----------------
        def _compute_offdiag_flux():
            """Return ``[k1, (k2), (k3)]`` via summed Gaussian flux kicks."""
            h = h_flux
            _pe_cache: dict[tuple, float] = {}

            def _pe_kicks(kicks_tuple: tuple) -> float:
                if kicks_tuple in _pe_cache:
                    return _pe_cache[kicks_tuple]
                stim_signal = np.zeros(len(pulse_tlist), dtype=float)
                for idx, mult in kicks_tuple:
                    t_p = t_samples[idx]
                    gauss = (mult * h) * np.exp(
                        -0.5 * ((pulse_tlist - t_p) / self.stim_width) ** 2
                    )
                    stim_signal += gauss
                stim = FluxSignal(type=8, t_list=pulse_tlist,
                                  signal=stim_signal)
                qubit_t = qubit.qubit_under_mag(stim)
                H_stim = QobjEvo(
                    qubit.qubit_under_mag_hamiltonian(
                        qubit_t, stim.t_list, frame_val, omega_d_val,
                    ),
                    tlist=stim.t_list, order=1,
                )
                result = mesolve(
                    H_base + H_stim, qubit.state, pulse_tlist, [],
                    e_ops=[psi_e * psi_e.dag()],
                    options=_fd_opts,
                )
                _pe_cache[kicks_tuple] = float(result.expect[0][-1])
                return _pe_cache[kicks_tuple]

            def _pe(*specs):
                return _pe_kicks(tuple(sorted(specs, key=lambda x: x[0])))

            return _offdiag_stencils(M, h, _pe)

        # -- run -----------------------------------------------------------------
        if self.mode == 'flux':
            kernels = _compute_offdiag_flux()
        elif self.virtual_z_impl == 'hardware':
            if not isinstance(pulse, CompositePulse):
                warnings.warn(
                    "virtual_z_impl='hardware' on a single Pulse falls back "
                    "to math implementation (σ_z impulse) for off-diagonal "
                    "extraction."
                )
                if not (self.richardson and N >= 2):
                    kernels = _compute_offdiag_omega(self._sigma_t())
                else:
                    sigmas = self._richardson_sigma_values()
                    per_sigma = [_compute_offdiag_omega(st) for st in sigmas]
                    kernels = self._richardson_extrapolate(sigmas, per_sigma)
                    i_min = int(np.argmin(sigmas))
                    kernels[0] = per_sigma[i_min][0]
            else:
                if self.richardson and N >= 2:
                    warnings.warn(
                        "Richardson extrapolation is not applicable to "
                        "virtual_z_impl='hardware' (phase kicks are "
                        "instantaneous, no probe width to extrapolate); "
                        "ignoring richardson=True.",
                        RuntimeWarning,
                    )
                kernels = _compute_offdiag_omega_hardware()
        elif not (self.richardson and N >= 2):
            kernels = _compute_offdiag_omega(self._sigma_t())
        else:
            # Richardson: compute at multiple sigma_t → extrapolate each element
            sigmas = self._richardson_sigma_values()
            per_sigma = [_compute_offdiag_omega(st) for st in sigmas]
            kernels = self._richardson_extrapolate(sigmas, per_sigma)
            # k1 kept from smallest sigma_t (already O(σ_t²) accurate)
            i_min = int(np.argmin(sigmas))
            kernels[0] = per_sigma[i_min][0]

        elapsed = time.time() - t_start_full
        print(
            f"FD off-diagonal kernel (exp {self.mode}, order={N}, M={M}) "
            f"took {elapsed:.1f}s"
        )
        return t_samples, kernels

    # ------------------------------------------------------------------
    # Stimulus amplitude helpers
    # ------------------------------------------------------------------

    def _sigma_t(self) -> float:
        """Resolve the VZ Gaussian probe width σ_t (ns) for omega exp paths.

        Returns ``self.probe_sigma_t`` when set, else the legacy default
        ``2.0 * CONFIG.awg.dt``.  Smaller σ_t reduces the off-diagonal
        smearing bias of the high-order diagonal kernel, but σ_t < dt is
        under-resolved on the time grid; callers should keep σ_t >= dt.
        """
        if self.probe_sigma_t is not None:
            return float(self.probe_sigma_t)
        return 2.0 * CONFIG.awg.dt

    def _richardson_sigma_values(self) -> np.ndarray:
        """Resolve the σ_t sample points (ns) for Richardson extrapolation.

        ``richardson_sigmas`` is interpreted as multiples of ``dt``;
        ``None`` -> ``(2.0, 1.5, 1.0) × dt``.  All values are kept
        ``>= dt`` (a Gaussian narrower than the grid step is under-
        resolved and the stencil blows up).
        """
        dt = CONFIG.awg.dt
        mults = self.richardson_sigmas or (2.0, 1.5, 1.0)
        sig = np.array([float(m) * dt for m in mults], dtype=float)
        return np.maximum(sig, dt)

    @staticmethod
    def _richardson_extrapolate(sigmas: np.ndarray, per_sigma: list) -> list:
        """Extrapolate kernels to σ_t -> 0 per (order, element).

        Supports both 1-D diagonal kernels ``(M,)`` and n-D off-diagonal
        kernels ``(M,M)``, ``(M,M,M)`` — each tensor element is independently
        extrapolated via ``np.polyfit`` to σ_t = 0.

        Parameters
        ----------
        sigmas : np.ndarray
            σ_t sample points (ns), length ``S``.
        per_sigma : list
            ``per_sigma[s]`` is ``[k1, ..., kN]`` at σ_t ``sigmas[s]``.
            Each ``k_n`` may be shape ``(M,)``, ``(M,M)``, or ``(M,M,M)``.

        Returns
        -------
        list
            ``[k1, ..., kN]`` extrapolated to σ_t = 0, same shapes as input.
            With a single σ_t the input is returned unchanged.
        """
        sigmas = np.asarray(sigmas, dtype=float)
        S = len(sigmas)
        N = len(per_sigma[0])
        if S == 1:
            return list(per_sigma[0])
        deg = S - 1
        out = []
        for n in range(N):
            # stack across sigma: (S, *kernel_shape)
            stack = np.stack([per_sigma[s][n] for s in range(S)], axis=0)
            kernel_shape = stack.shape[1:]
            flat = stack.reshape(S, -1)               # (S, K)
            ext_flat = np.zeros(flat.shape[1], dtype=float)
            for k in range(flat.shape[1]):
                coeffs = np.polyfit(sigmas, flat[:, k], deg=deg)
                ext_flat[k] = coeffs[-1]              # σ_t = 0 (constant term)
            out.append(ext_flat.reshape(kernel_shape))
        return out

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
