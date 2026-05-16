"""sqc.reconstruction.transient — transient-field sensing reconstruction.

    method="wiener":      linear Wiener deconvolution
    method="hammerstein": Hammerstein-Wiener nonlinear block model
    method="lm":          Levenberg-Marquardt full density-matrix inversion
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from qutip import Qobj, QobjEvo, basis, expect, mesolve, qeye

from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.control.pulse import CompositePulse
from sqc.reconstruction.base import Reconstruction
from sqc.reconstruction.basis import (
    R as _R,
    BasisType,
    basis_function_decomposition,
    generate_basis_functions,
    regularization_matrix,
)


# ---------------------------------------------------------------------------
# shared helper
# ---------------------------------------------------------------------------

def _get_qubit_params(qubit) -> tuple[float, float, float]:
    """Extract (EC, EJ, frequency) from QubitSpec or legacy TransmonQubit."""
    EC = qubit.EC
    EJ = getattr(qubit, "EJ_0", qubit.EJ)
    freq = qubit.frequency() if callable(qubit.frequency) else qubit.frequency
    return EC, EJ, freq


# ===================================================================
# TransientReconstruction
# ===================================================================

@dataclass
class TransientReconstruction(Reconstruction):
    """Transient-field sensing reconstruction.

    Parameters
    ----------
    method : str
        - ``"wiener"``: linear Wiener deconvolution.
        - ``"hammerstein"``: Hammerstein-Wiener nonlinear block model.
        - ``"lm"``: Levenberg-Marquardt full-density-matrix inversion.

    wiener / hammerstein params
    ---------------------------
    lambda_reg : float
        Regularisation parameter.  Default from CONFIG.
    qubit : QubitSpec or TransmonQubit, optional
        Required for hammerstein; unused by wiener.

    lm params
    ---------
    control_pulse : CompositePulse, optional
        Control pulse for the sliding measurement.  Required for lm.
    basis_type : str
        Basis type for signal parameterisation (lm only).
    n_basis : int
        Number of basis functions (lm only).
    max_iter : int
        Maximum LM iterations.
    tol : float
        Convergence tolerance (relative ||db||).
    mu_init : float
        Initial damping parameter.
    use_adjoint : bool
        Use adjoint Jacobian (True) or finite-difference (False).
    """

    method: Literal["wiener", "hammerstein", "lm"] = "wiener"

    # -- wiener / hammerstein --
    lambda_reg: float = field(
        default_factory=lambda: CONFIG.reconstruction.lambda_reg
    )
    qubit: object | None = None

    # -- lm --
    control_pulse: CompositePulse | None = None
    basis_type: BasisType = field(
        default_factory=lambda: CONFIG.reconstruction.lm_basis_type
    )
    n_basis: int = field(
        default_factory=lambda: CONFIG.reconstruction.lm_n_basis
    )
    # lm.lambda_reg aliases self.lambda_reg above (reused)
    max_iter: int = field(
        default_factory=lambda: CONFIG.reconstruction.lm_max_iter
    )
    tol: float = field(
        default_factory=lambda: CONFIG.reconstruction.lm_tol
    )
    mu_init: float = field(
        default_factory=lambda: CONFIG.reconstruction.lm_mu_init
    )
    use_adjoint: bool = True

    # ------------------------------------------------------------------
    def reconstruct(self, measurement, kernel=None, **kwargs):
        match self.method:
            case "wiener":
                return self._reconstruct_wiener(measurement, kernel, **kwargs)
            case "hammerstein":
                return self._reconstruct_hammerstein(measurement, kernel, **kwargs)
            case "lm":
                return self._reconstruct_lm(measurement, **kwargs)
            case _:
                raise ValueError(f"Unknown method: {self.method}")

    # ==================================================================
    # wiener
    # ==================================================================

    def _reconstruct_wiener(
        self, measurement, kernel: np.ndarray, dt: float | None = None, **__
    ) -> FluxSignal:
        if hasattr(measurement, "data") and isinstance(measurement.data, dict):
            delta_p = np.asarray(measurement.data["delta_p"], dtype=float)
        else:
            delta_p = np.asarray(measurement, dtype=float)

        if dt is None:
            if hasattr(measurement, "axes") and "scan" in measurement.axes:
                scan = np.asarray(measurement.axes["scan"])
                dt = scan[1] - scan[0]
            else:
                dt = 1.0

        N_del = len(delta_p)
        N_ker = len(kernel)
        N = N_del - N_ker + 1
        N_fft = N_del

        p_pad = np.zeros(N_fft); p_pad[:N_del] = delta_p
        k_pad = np.zeros(N_fft); k_pad[:N_ker] = kernel

        Y = np.fft.fft(p_pad)
        H = np.fft.fft(k_pad)
        G = np.conj(H) / (np.abs(H) ** 2 + self.lambda_reg ** 2)
        X_w = Y * G / dt
        x_rec = np.real(np.fft.ifft(X_w))[:N]

        return FluxSignal(
            type=8,
            t_list=np.linspace(0, (N - 1) * dt, N),
            signal=x_rec,
        )

    # ==================================================================
    # hammerstein
    # ==================================================================

    def _reconstruct_hammerstein(
        self, measurement, kernel, dt: float | None = None, **__
    ) -> FluxSignal:
        # Step 1: linear wiener → ω(t)
        omega_signal = self._reconstruct_wiener(measurement, kernel, dt=dt)
        omega = np.asarray(omega_signal.signal, dtype=float)

        # Step 2: inverse Transmon dispersion
        EC, EJ, freq = _get_qubit_params(self.qubit)
        ratio = np.clip((omega + freq + EC) ** 2 / (8 * EC * EJ), 0.0, 1.0)
        B = (1.0 / np.pi) * np.arccos(ratio)

        # Subtract working-point flux bias to recover h = Φ - Φ_bias
        bias = getattr(self.qubit, "flux_bias",
                       getattr(self.qubit, "flux", 0.0))
        B = B - bias

        return FluxSignal(
            type=8, t_list=omega_signal.t_list.copy(), signal=B,
        )

    # ==================================================================
    # lm  (Levenberg-Marquardt)
    # ==================================================================

    def _reconstruct_lm(self, measurement, **kwargs) -> tuple[FluxSignal, dict]:
        p_meas = np.asarray(measurement.data["p_meas"], dtype=float)
        t_list = np.asarray(measurement.axes["t_signal"], dtype=float)

        B_init = FluxSignal(
            type=6, t_list=t_list,
            n_basis=self.n_basis, basis_type=self.basis_type,
        )
        initial_guess = np.asarray(
            kwargs.get("initial_guess", np.zeros_like(t_list)), dtype=float,
        )
        b_init = basis_function_decomposition(
            initial_guess, t_list, B_init.basis_functions,
        )
        B_init.update_signal(b=b_init)

        b_opt, history = self._levenberg_marquardt(p_meas, t_list, b_init, B_init)
        B_init.update_signal(b=b_opt)
        return B_init, history

    # ------------------------------------------------------------------
    # _build_h_for_signal
    # ------------------------------------------------------------------

    def _build_h_for_signal(self, t_list, t_meas):
        qubit = self.qubit
        cp = self.control_pulse
        n_meas = len(t_meas)
        H_list = [[] for _ in range(n_meas)]
        t_evolve_list = [[] for _ in range(n_meas)]  # type: ignore[assignment]

        for i, t_delay in enumerate(t_meas):
            delta = t_delay - 0.5 * cp.t_list[-1]
            t_start = min(t_list[0], delta)
            t_end = max(t_list[-1], delta + cp.t_list[-1])
            N_e = len(t_list) + len(cp.t_list) - 1
            t_evolve = np.linspace(t_start, t_end, N_e)
            t_evolve_list[i] = t_evolve

            freq_coeffs = np.zeros(N_e, dtype=float)
            for j, t in enumerate(t_evolve):
                if t_list[0] <= t <= t_list[-1]:
                    index = min(np.searchsorted(t_list, t), len(t_list) - 1)
                    freq_coeffs[j] = (
                        qubit.freq_coeffs[index] if cp.frame == 0
                        else qubit.freq_coeffs[index] - cp.omega_d
                    )
                else:
                    freq_coeffs[j] = (
                        qubit.frequency if cp.frame == 0
                        else qubit.frequency - cp.omega_d
                    )

            id_coeffs = np.ones(N_e, dtype=complex)
            if cp.frame == 0:
                H_list[i].append([
                    qubit.anharmonicity * qubit.n * (qubit.n - 1) * 0.5,
                    id_coeffs,
                ])
                H_list[i].append([
                    qubit.n + 0.5 * qeye(qubit.n_levels), freq_coeffs,
                ])
            else:
                H_list[i].append([
                    qubit.anharmonicity * qubit.n * (qubit.n - 1) * 0.5,
                    id_coeffs,
                ])
                H_list[i].append([qubit.n, freq_coeffs])

            for op, coeffs in cp.hamiltonian:
                coeff_global = np.zeros(N_e, dtype=complex)
                for j, t in enumerate(t_evolve):
                    t_loc = t - delta
                    if 0 <= t_loc <= cp.t_list[-1]:
                        index = min(np.searchsorted(cp.t_list, t_loc), len(cp.t_list) - 1)
                        coeff_global[j] = coeffs[index]
                H_list[i].append([op, coeff_global])

        return H_list, t_evolve_list

    # ------------------------------------------------------------------
    # _forward_simulation
    # ------------------------------------------------------------------

    def _forward_simulation(self, B_curr, t_meas, H_list, t_evolve_list):
        import time
        _ = B_curr
        start = time.time()
        result = []
        qubit = self.qubit

        for i, _t_i in enumerate(t_meas):
            H = H_list[i]
            t_evolve = t_evolve_list[i]
            dt_evolve = t_evolve[1] - t_evolve[0]
            options = {
                "store_states": True,
                "atol": 1e-10, "rtol": 1e-8,
                "max_step": dt_evolve / 2, "nsteps": 10000,
            }
            H_total = QobjEvo(H, tlist=t_evolve, order=1)
            result.append(mesolve(
                H_total, qubit.state, t_evolve, [],
                e_ops=[basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()],
                options=options,
            ))
        print(f"Forward simulation time: {time.time() - start:.2f} s")
        return result

    # ------------------------------------------------------------------
    # _compute_jacobian_adjoint
    # ------------------------------------------------------------------

    def _compute_jacobian_adjoint(self, B_curr, t_meas, results, H_list, t_evolve_list):
        import time
        from scipy.integrate import solve_ivp

        start = time.time()
        qubit = self.qubit
        N = len(t_meas)
        M = len(B_curr.params["b"])
        J = np.zeros((N, M))
        dim = qubit.state.shape[0]

        tB_list = B_curr.t_list
        sensitivity = [
            qubit.frequency_sensitivity(qubit.flux + B_curr.value_at(t))
            for t in tB_list
        ]

        states = []
        for res in results:
            states.append(np.array([(s * s.dag()).full() for s in res.states]))

        c_ops_list = [qubit.c_ops[k].full() for k in range(len(qubit.c_ops))]
        ops = [op.full() for op, _ in H_list[0]]
        mu_0 = basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()
        mu_0 = mu_0.full().flatten()
        G_mat = qubit.n.full()

        for i in range(N):
            t_evolve = t_evolve_list[i]
            t_m, t_M = t_evolve[0], t_evolve[-1]
            N_e = len(t_evolve)

            t_mapped = np.empty_like(t_evolve)
            for s in range(N_e):
                if t_evolve[s] < tB_list[0]:
                    t_mapped[s] = tB_list[0]
                elif t_evolve[s] >= tB_list[-1]:
                    t_mapped[s] = tB_list[-1]
                else:
                    idx = min(np.searchsorted(tB_list, t_evolve[s]), len(tB_list) - 1)
                    t_mapped[s] = tB_list[idx]

            H_evolve = np.zeros((N_e, dim, dim), dtype=complex)
            for j in range(N_e):
                for m_idx, (op_mat, coeffs) in enumerate(H_list[i]):
                    c = coeffs[j] if isinstance(coeffs, np.ndarray) else coeffs
                    H_evolve[j] += ops[m_idx] * c

            def adjoint(s, mu, *args):
                mu = mu.reshape((dim, dim))
                t_curr = t_M - s
                j_float = (t_curr - t_m) / (t_M - t_m) * (N_e - 1)
                j_lo = max(int(np.floor(j_float)), 0)
                j_hi = min(j_lo + 1, N_e - 1)
                alpha = j_float - j_lo
                H = (1 - alpha) * H_evolve[j_lo] + alpha * H_evolve[j_hi]
                return (1j * (H @ mu - mu @ H)).flatten()

            s_list = t_M - t_evolve[::-1]
            ds = s_list[1] - s_list[0]
            sol = solve_ivp(
                adjoint, t_span=(s_list[0], s_list[-1]), y0=mu_0,
                method="RK45", t_eval=s_list,
                rtol=1e-6, atol=1e-8, max_step=ds * 10,
            )

            N_s = len(sol.t)
            lambda_t = np.zeros((N_s, dim, dim), dtype=complex)
            for n in range(N_s):
                lambda_t[N_s - 1 - n] = sol.y[:, n].reshape((dim, dim))

            rho_stack = states[i][:N_s]
            commutator = G_mat @ rho_stack - rho_stack @ G_mat
            trace_vals = np.trace(lambda_t @ commutator, axis1=-2, axis2=-1)
            g_values = -1j * trace_vals

            for s in range(N_s):
                if t_evolve[s] < tB_list[0] or t_evolve[s] >= tB_list[-1]:
                    g_values[s] *= 0.0
                else:
                    idx = min(np.searchsorted(tB_list, t_evolve[s]), len(sensitivity) - 1)
                    g_values[s] *= sensitivity[idx]

            for k, phi_k in enumerate(B_curr.basis_functions):
                J[i, k] = float(np.real(
                    np.trapezoid(g_values * phi_k(t_mapped), t_evolve)
                ))

        print(f"Jacobian computation time: {time.time() - start:.2f} s")
        return J

    # ------------------------------------------------------------------
    # _compute_jacobian_fd  (fixed: reuses _build_h_for_signal)
    # ------------------------------------------------------------------

    def _compute_jacobian_fd(
        self, B_curr, t_meas, p_sim, H_list, t_evolve_list, epsilon=1e-6,
    ):
        qubit = self.qubit
        N = len(t_meas)
        M = len(B_curr.params["b"])
        J_fd = np.zeros((N, M))
        t_signal = np.asarray(B_curr.t_list, dtype=float)

        for k in range(M):
            b_perturb = np.asarray(B_curr.params["b"], dtype=float).copy()
            b_perturb[k] += epsilon
            B_perturb = B_curr.copy()
            B_perturb.update_signal(b=b_perturb)

            qubit.qubit_in_mag(B_perturb)
            H_perturb, t_evolve_perturb = self._build_h_for_signal(t_signal, t_meas)
            result_perturb = self._forward_simulation(
                B_perturb, t_meas, H_perturb, t_evolve_perturb,
            )
            p_sim_perturb = np.array([res.expect[0][-1] for res in result_perturb])
            J_fd[:, k] = (p_sim_perturb - p_sim) / epsilon

        qubit.qubit_in_mag(B_curr)
        return J_fd

    # ------------------------------------------------------------------
    # _levenberg_marquardt
    # ------------------------------------------------------------------

    def _levenberg_marquardt(self, p_meas, t_list, b_init, B_init):
        qubit = self.qubit
        cp = self.control_pulse
        n_basis = len(b_init)

        b = np.asarray(b_init, dtype=float).copy()
        B_curr = B_init.copy()
        basis_type_val = B_curr.params["basis_type"]
        mu = self.mu_init
        reg = self.lambda_reg

        qubit.qubit_in_mag(B_curr)

        history = {"b": [b.copy()], "res": [], "mu": [mu]}

        meas_start = t_list[0] - 0.5 * cp.t_list[-1]
        meas_end = t_list[-1] + 0.5 * cp.t_list[-1]
        t_meas = np.linspace(meas_start, meas_end, len(t_list) + len(cp.t_list) - 1)

        def build_H():
            return self._build_h_for_signal(t_list, t_meas)

        H_curr_list, t_evolve_list = build_H()

        for iteration in range(self.max_iter):
            results = self._forward_simulation(B_curr, t_meas, H_curr_list, t_evolve_list)
            p_sim = np.array([res.expect[0][-1] for res in results])
            res = p_meas - p_sim
            history["res"].append(res.copy())

            if self.use_adjoint:
                J = self._compute_jacobian_adjoint(B_curr, t_meas, results, H_curr_list, t_evolve_list)
            else:
                J = self._compute_jacobian_fd(B_curr, t_meas, p_sim, H_curr_list, t_evolve_list)

            A = J.T @ J + mu * np.eye(n_basis) + reg * regularization_matrix(n_basis, basis_type_val)
            delta_b = np.linalg.solve(A, J.T @ res)

            b_trial = b + delta_b
            B_trial = B_curr.copy()
            B_trial.update_signal(b=b_trial)

            qubit.qubit_in_mag(B_trial)
            H_trial_list, _ = build_H()
            result_trial = self._forward_simulation(B_trial, t_meas, H_trial_list, t_evolve_list)
            p_sim_trial = np.array([res.expect[0][-1] for res in result_trial])
            res_trial = p_meas - p_sim_trial

            if np.linalg.norm(res) > np.linalg.norm(res_trial):
                b = b_trial; B_curr = B_trial; H_curr_list = H_trial_list
                mu = max(mu / 2, 1e-8)
                history["b"].append(b.copy()); history["mu"].append(mu)
                if np.linalg.norm(delta_b) / (np.linalg.norm(b) + 1e-8) < self.tol:
                    print(f"Converged at iteration {iteration}")
                    break
            else:
                mu = min(mu * 2, 1e8)
                qubit.qubit_in_mag(B_curr)
                H_curr_list, _ = build_H()
                history["b"].append(b.copy()); history["mu"].append(mu)

            print(f"Iter {iteration}: res={np.linalg.norm(res):.6f}, "
                  f"res_trial={np.linalg.norm(res_trial):.6f}, mu={mu:.6e}")

        return b, history
