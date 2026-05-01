"""sqc.reconstruction.numerical_inverse — Levenberg-Marquardt waveform reconstruction.

Verbatim port from src/analysis.py:272-825 (Track B 0.3 fixed version).

LMReconstruction implements full-density-matrix numerical inversion using
basis-function parameterisation of the unknown flux signal. The optimiser
minimises ||p_meas - p_sim(b)||^2 + reg * ||R b||^2 where p_sim is
computed by QuTiP mesolve at each delay time.

Key algorithms:
  - _forward_simulation: run mesolve for each measurement delay
  - _compute_jacobian_adjoint: adjoint-state method (solve_ivp backward)
  - _compute_jacobian_fd: finite-difference fallback
  - _levenberg_marquardt: the main LM optimisation loop

See _refactor_plan.md §7.5, phase_3_handbook.md §3.4.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from qutip import Qobj, QobjEvo, basis, expect, mesolve, qeye

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


@dataclass
class LMReconstruction(Reconstruction):
    """Levenberg-Marquardt full-density-matrix waveform reconstruction.

    Reconstructs the magnetic flux signal B(t) (in units of Phi_0) from
    a sliding-measurement excited-state probability trace p_meas(tau).

    The unknown signal is parameterised as a linear combination of basis
    functions (B-spline, Fourier, or Legendre). The optimiser uses
    Levenberg-Marquardt with either adjoint or finite-difference Jacobian.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit instance (must support qubit_in_mag, frequency_sensitivity, etc.).
    control_pulse : CompositePulse
        Control pulse used for the sliding measurement.
    basis_type : {"bspline", "fourier", "legendre"}
        Basis function type for signal parameterisation.
    n_basis : int
        Number of basis functions.
    lambda_reg : float
        Regularisation strength (smoothness penalty).
    max_iter : int
        Maximum Levenberg-Marquardt iterations.
    tol : float
        Convergence tolerance (relative change in ||b||).
    mu_init : float
        Initial LM damping parameter.
    use_adjoint : bool
        If True, use adjoint method for Jacobian; else finite difference.

    References
    ----------
    Gao 2021 §V.C (readout) — not directly, but our reconstruction is
    patterned after transmon frequency-flux dispersion (Gao 2021 §II.C).
    """

    qubit: object  # TransmonQubit (duck-typed for src.qubit compat)
    control_pulse: CompositePulse
    basis_type: BasisType = "fourier"
    n_basis: int = 100
    lambda_reg: float = 100.0
    max_iter: int = 10
    tol: float = 1e-6
    mu_init: float = 1e-3
    use_adjoint: bool = True

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def reconstruct(
        self,
        measurement,
        kernel=None,
        calibration=None,
        initial_guess: np.ndarray | None = None,
    ) -> tuple[FluxSignal, dict]:
        """Reconstruct B(t) from measurement data p_meas(tau).

        Parameters
        ----------
        measurement : ExperimentResult
            Must contain measurement.data["p_meas"] (1-D p_e array)
            and measurement.axes["t_signal"] (time axis).
        kernel : None
            Not used by LM (kernel is implicit in full-density simulation).
        calibration : None
            Not used by LM (frequency response is simulated directly).
        initial_guess : np.ndarray or None
            Initial guess for B(t). If None, uses zeros.

        Returns
        -------
        tuple[FluxSignal, dict]
            (B_opt, history) where B_opt is the reconstructed flux signal
            and history contains "b", "res", "mu" arrays per iteration.
        """
        p_meas = np.asarray(measurement.data["p_meas"], dtype=float)
        t_list = np.asarray(measurement.axes["t_signal"], dtype=float)

        # Build initial basis-expanded signal
        B_init = FluxSignal(
            type=6,
            t_list=t_list,
            n_basis=self.n_basis,
            basis_type=self.basis_type,
        )
        if initial_guess is None:
            initial_guess = np.zeros_like(t_list)
        else:
            initial_guess = np.asarray(initial_guess, dtype=float)

        b_init = basis_function_decomposition(
            initial_guess, t_list, B_init.basis_functions
        )
        B_init.update_signal(b=b_init)

        # Run LM optimisation
        b_opt, history = self._levenberg_marquardt(
            p_meas, t_list, b_init, B_init,
        )

        B_init.update_signal(b=b_opt)
        return B_init, history

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _build_h_for_signal(self, t_list, t_meas):
        """Build Hamiltonian list and evolve time axes for each delay.

        This is the ``H(qubit)`` inner function from the original
        ``levenberg_marquardt`` (src/analysis.py:699-734), verbatim port.

        Uses ``self.qubit.freq_coeffs`` (set by ``qubit_in_mag`` before
        calling this method) to build frame-corrected frequency coefficients.

        Parameters
        ----------
        t_list : np.ndarray
            Signal time axis (same as B_curr.t_list).
        t_meas : np.ndarray
            Measurement (delay) time axis.

        Returns
        -------
        H_list : list[list]
            H_list[i] is the QuTiP list-format Hamiltonian for delay i.
        t_evolve_list : list[np.ndarray]
            t_evolve_list[i] is the time axis for delay i.
        """
        qubit = self.qubit
        cp = self.control_pulse
        n_meas = len(t_meas)

        H_list: list = [[] for _ in range(n_meas)]
        t_evolve_list: list = [[] for _ in range(n_meas)]  # type: ignore[assignment]

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
                    index = np.searchsorted(t_list, t)
                    index = min(index, len(t_list) - 1)
                    freq_coeffs[j] = (
                        qubit.freq_coeffs[index]
                        if cp.frame == 0
                        else qubit.freq_coeffs[index] - cp.omega_d
                    )
                else:
                    freq_coeffs[j] = (
                        qubit.frequency
                        if cp.frame == 0
                        else qubit.frequency - cp.omega_d
                    )

            id_coeffs = np.ones(N_e, dtype=complex)
            if cp.frame == 0:
                H_list[i].append(
                    [
                        qubit.anharmonicity
                        * qubit.n
                        * (qubit.n - 1)
                        * 0.5,
                        id_coeffs,
                    ]
                )
                H_list[i].append(
                    [
                        qubit.n + 0.5 * qeye(qubit.n_levels),
                        freq_coeffs,
                    ]
                )
            else:
                H_list[i].append(
                    [
                        qubit.anharmonicity
                        * qubit.n
                        * (qubit.n - 1)
                        * 0.5,
                        id_coeffs,
                    ]
                )
                H_list[i].append([qubit.n, freq_coeffs])

            for op, coeffs in cp.hamiltonian:
                coeff_global = np.zeros(N_e, dtype=complex)
                for j, t in enumerate(t_evolve):
                    t_loc = t - delta
                    if 0 <= t_loc <= cp.t_list[-1]:
                        index = np.searchsorted(cp.t_list, t_loc)
                        index = min(index, len(cp.t_list) - 1)
                        coeff_global[j] = coeffs[index]
                H_list[i].append([op, coeff_global])

        return H_list, t_evolve_list

    def _forward_simulation(
        self, B_curr, t_meas, H_list, t_evolve_list
    ):
        """Forward simulation: mesolve for each measurement delay.

        Verbatim port of src/analysis.py:402-446.

        Parameters
        ----------
        B_curr : FluxSignal
            Current flux signal guess (used only for time axis reference).
        t_meas : np.ndarray
            Measurement delay times.
        H_list : list
            Pre-built Hamiltonian lists for each delay.
        t_evolve_list : list[np.ndarray]
            Pre-built evolve time axes for each delay.

        Returns
        -------
        list[qutip.Result]
            mesolve results for each delay time.
        """
        import time

        _ = B_curr  # retained for backward-compat signature
        start = time.time()
        result: list = []

        qubit = self.qubit

        for i, _t_i in enumerate(t_meas):
            H = H_list[i]
            t_evolve = t_evolve_list[i]
            dt_evolve = t_evolve[1] - t_evolve[0]
            options = {
                "store_states": True,
                "atol": 1e-10,
                "rtol": 1e-8,
                "max_step": dt_evolve / 2,
                "nsteps": 10000,
            }

            H_total = QobjEvo(H, tlist=t_evolve, order=1)

            result.append(
                mesolve(
                    H_total,
                    qubit.state,
                    t_evolve,
                    [],
                    e_ops=[
                        basis(qubit.n_levels, 1)
                        * basis(qubit.n_levels, 1).dag()
                    ],
                    options=options,
                )
            )
        end = time.time()
        print(f"Forward simulation time: {end - start:.2f} s")

        return result

    def _compute_jacobian_adjoint(
        self, B_curr, t_meas, results, H_list, t_evolve_list
    ):
        """Adjoint-state Jacobian computation.

        Verbatim port of src/analysis.py:448-597.

        Solves the backward-in-time adjoint equation via solve_ivp for
        each measurement row, then contracts with the forward-trajectory
        commutator to obtain d(p_e_i)/d(b_k).

        Parameters
        ----------
        B_curr : FluxSignal
            Current signal (type=6, with .basis_functions, .params["b"]).
        t_meas : np.ndarray
            Measurement delay times.
        results : list[qutip.Result]
            Forward-simulation results (store_states=True required).
        H_list : list
            Pre-built Hamiltonian lists.
        t_evolve_list : list[np.ndarray]
            Pre-built evolve time axes.

        Returns
        -------
        np.ndarray (N, M)
            Jacobian matrix: J[i, k] = d(p_e_i)/d(b_k).
        """
        import time

        from scipy.integrate import solve_ivp

        start = time.time()
        qubit = self.qubit

        N = len(t_meas)
        M = len(B_curr.params["b"])
        J = np.zeros((N, M))
        dim = qubit.state.shape[0]

        B = B_curr
        tB_list = B.t_list

        # ---- sensitivity at each signal time point -----------------------
        sensitivity = [
            qubit.frequency_sensitivity(qubit.flux + B.value_at(t))
            for t in tB_list
        ]

        # ---- extract states from forward simulation ----------------------
        result = results
        states: list = [[] for _ in range(len(result))]
        for i, res in enumerate(result):
            states[i] = np.array(
                [(s * s.dag()).full() for s in res.states]
            )

        # ---- pre-compute constant matrices -------------------------------
        c_ops_list = [qubit.c_ops[k].full() for k in range(len(qubit.c_ops))]
        # c_ops_dag currently unused (dissipation commented out in adjoint rhs)
        c_ops_dag_list = [
            c_ops_list[k].conj().T for k in range(len(qubit.c_ops))
        ]
        ops = [op.full() for op, _ in H_list[0]]
        mu_0 = basis(qubit.n_levels, 1) * basis(qubit.n_levels, 1).dag()
        mu_0 = mu_0.full().flatten()

        G_mat = qubit.n.full()

        # ---- loop over measurement times ---------------------------------
        for i in range(N):
            start_i = time.time()
            _t_i = t_meas[i]  # noqa: F841 — retained for debug reference
            t_evolve = t_evolve_list[i]
            t_m = t_evolve[0]
            t_M = t_evolve[-1]
            N_e = len(t_evolve)

            # Map evolve times to signal times (t_mapped)
            t_mapped = np.empty_like(t_evolve)
            for s in range(N_e):
                if t_evolve[s] < tB_list[0]:
                    t_mapped[s] = tB_list[0]
                elif t_evolve[s] >= tB_list[-1]:
                    t_mapped[s] = tB_list[-1]
                else:
                    index = np.searchsorted(tB_list, t_evolve[s])
                    index = min(index, len(tB_list) - 1)
                    t_mapped[s] = tB_list[index]

            # Build H_evolve matrix array
            H_evolve = np.zeros((N_e, dim, dim), dtype=complex)
            for j in range(N_e):
                for m, (op, coeffs) in enumerate(H_list[i]):
                    coeff = (
                        coeffs[j]
                        if isinstance(coeffs, np.ndarray)
                        else coeffs
                    )
                    H_evolve[j] += ops[m] * coeff

            # Adjoint ODE right-hand side
            def adjoint(s, mu, *args):
                mu = mu.reshape((dim, dim))
                t_curr = t_M - s
                j_float = (t_curr - t_m) / (t_M - t_m) * (N_e - 1)
                j_lo = max(int(np.floor(j_float)), 0)
                j_hi = min(j_lo + 1, N_e - 1)
                alpha = j_float - j_lo
                H = (1 - alpha) * H_evolve[j_lo] + alpha * H_evolve[j_hi]

                rhs = 1j * (H @ mu - mu @ H)
                # Dissipation currently omitted (matches original code):
                # for k in range(len(qubit.c_ops)):
                #     rhs += c_ops_dag_list[k] @ mu @ c_ops_list[k] \
                #         - 0.5 * (c_ops_dag_list[k] @ c_ops_list[k] @ mu
                #                + mu @ c_ops_dag_list[k] @ c_ops_list[k])
                return rhs.flatten()

            s_list = t_M - t_evolve[::-1]
            ds = s_list[1] - s_list[0]

            sol = solve_ivp(
                adjoint,
                t_span=(s_list[0], s_list[-1]),
                y0=mu_0,
                method="RK45",
                t_eval=s_list,
                rtol=1e-6,
                atol=1e-8,
                max_step=ds * 10,
            )

            N_s = len(sol.t)
            lambda_t = np.zeros((N_s, dim, dim), dtype=complex)
            for n in range(N_s):
                lambda_t[N_s - 1 - n] = sol.y[:, n].reshape((dim, dim))

            # Gradient contribution
            rho_stack = states[i][:N_s]
            commutator = (
                G_mat @ rho_stack - rho_stack @ G_mat
            )
            trace_values = np.trace(
                lambda_t @ commutator, axis1=-2, axis2=-1
            )
            g_values = -1j * trace_values

            # Zero out contribution outside signal domain
            for s in range(N_s):
                if t_evolve[s] < tB_list[0]:
                    g_values[s] *= 0.0
                elif t_evolve[s] >= tB_list[-1]:
                    g_values[s] *= 0.0
                else:
                    index = np.searchsorted(tB_list, t_evolve[s])
                    index = min(index, len(sensitivity) - 1)
                    g_values[s] *= sensitivity[index]

            # Project onto basis function derivatives
            for k, phi_k in enumerate(B_curr.basis_functions):
                integrand_k = g_values * phi_k(t_mapped)
                integral = np.trapezoid(integrand_k, t_evolve)
                J[i, k] = np.real(integral)

            end_i = time.time()
            # print(f"Jacobian row {i} computation time: {end_i - start_i:.2f} s")

        end = time.time()
        print(f"Jacobian computation time: {end - start:.2f} s")

        return J

    def _compute_jacobian_fd(
        self, B_curr, t_meas, p_sim, H_list, t_evolve_list, epsilon=1e-6
    ):
        """Finite-difference Jacobian (reference / fallback).

        Verbatim port of src/analysis.py:599-662.

        Parameters
        ----------
        B_curr : FluxSignal
            Current signal.
        t_meas : np.ndarray
            Measurement delay times.
        p_sim : np.ndarray
            Simulated p_e values at b_curr.
        H_list : list
            Pre-built Hamiltonians (not used — rebuilt for each perturbation).
        t_evolve_list : list[np.ndarray]
            Pre-built evolve axes (not used — rebuilt).
        epsilon : float
            Finite-difference step size.

        Returns
        -------
        np.ndarray (N, M)
            Finite-difference Jacobian.
        """
        qubit = self.qubit
        cp = self.control_pulse

        N = len(t_meas)
        M = len(B_curr.params["b"])
        J_fd = np.zeros((N, M))
        t_signal = np.asarray(B_curr.t_list, dtype=float)

        def _build_h_for_current_qubit():
            """Inner H-list builder using current qubit freq_coeffs."""
            H_local: list = [[] for _ in range(N)]
            t_evolve_local: list = [[] for _ in range(N)]

            for i, t_delay in enumerate(t_meas):
                delta = t_delay - 0.5 * cp.t_list[-1]
                t_start = min(t_signal[0], delta)
                t_end = max(t_signal[-1], delta + cp.t_list[-1])
                N_e = len(t_signal) + len(cp.t_list) - 1
                t_evolve = np.linspace(t_start, t_end, N_e)
                t_evolve_local[i] = t_evolve

                freq_coeffs = np.zeros(N_e, dtype=float)
                for j, t in enumerate(t_evolve):
                    if t_signal[0] <= t <= t_signal[-1]:
                        index = np.searchsorted(t_signal, t)
                        index = min(index, len(t_signal) - 1)
                        freq_coeffs[j] = (
                            qubit.freq_coeffs[index]
                            if cp.frame == 0
                            else qubit.freq_coeffs[index] - cp.omega_d
                        )
                    else:
                        freq_coeffs[j] = (
                            qubit.frequency
                            if cp.frame == 0
                            else qubit.frequency - cp.omega_d
                        )

                id_coeffs = np.ones(N_e, dtype=complex)
                if cp.frame == 0:
                    H_local[i].append(
                        [
                            qubit.anharmonicity
                            * qubit.n
                            * (qubit.n - 1)
                            * 0.5,
                            id_coeffs,
                        ]
                    )
                    H_local[i].append(
                        [
                            qubit.n + 0.5 * qeye(qubit.n_levels),
                            freq_coeffs,
                        ]
                    )
                else:
                    H_local[i].append(
                        [
                            qubit.anharmonicity
                            * qubit.n
                            * (qubit.n - 1)
                            * 0.5,
                            id_coeffs,
                        ]
                    )
                    H_local[i].append([qubit.n, freq_coeffs])

                for op, coeffs in cp.hamiltonian:
                    coeff_global = np.zeros(N_e, dtype=complex)
                    for j, t in enumerate(t_evolve):
                        t_loc = t - delta
                        if 0 <= t_loc <= cp.t_list[-1]:
                            index = np.searchsorted(cp.t_list, t_loc)
                            index = min(index, len(cp.t_list) - 1)
                            coeff_global[j] = coeffs[index]
                    H_local[i].append([op, coeff_global])

            return H_local, t_evolve_local

        for k in range(M):
            b_perturb = np.asarray(B_curr.params["b"], dtype=float).copy()
            b_perturb[k] += epsilon
            B_perturb = B_curr.copy()
            B_perturb.update_signal(b=b_perturb)

            qubit.qubit_in_mag(B_perturb)
            H_perturb, t_evolve_perturb = _build_h_for_current_qubit()
            result_perturb = self._forward_simulation(
                B_perturb, t_meas, H_perturb, t_evolve_perturb
            )
            p_sim_perturb = np.array(
                [res.expect[0][-1] for res in result_perturb]
            )
            J_fd[:, k] = (p_sim_perturb - p_sim) / epsilon

        # Restore baseline state
        qubit.qubit_in_mag(B_curr)
        return J_fd

    def _levenberg_marquardt(
        self, p_meas, t_list, b_init, B_init
    ):
        """Levenberg-Marquardt optimisation loop.

        Verbatim port of src/analysis.py:666-825.

        Minimises ||p_meas - p_sim(b)||^2 + lambda_reg * ||R b||^2
        over basis coefficients b.

        Parameters
        ----------
        p_meas : np.ndarray
            Measured p_e values at each delay time.
        t_list : np.ndarray
            Signal time axis.
        b_init : np.ndarray
            Initial basis coefficients (n_basis,).
        B_init : FluxSignal
            Initial flux signal (type=6).

        Returns
        -------
        b_opt : np.ndarray
            Optimal basis coefficients.
        history : dict
            Dictionary with keys "b", "res", "mu" tracking iteration state.
        """
        qubit = self.qubit
        cp = self.control_pulse
        n_basis = len(b_init)

        b = np.asarray(b_init, dtype=float).copy()
        B_curr = B_init.copy()
        basis_type = B_curr.params["basis_type"]
        mu = self.mu_init
        reg = self.lambda_reg

        qubit.qubit_in_mag(B_curr)

        history: dict = {
            "b": [b.copy()],
            "res": [],
            "mu": [mu],
        }

        # ---- calibration for measurement time axis -----------------------
        meas_start = t_list[0] - 0.5 * cp.t_list[-1]
        meas_end = t_list[-1] + 0.5 * cp.t_list[-1]
        t_meas = np.linspace(
            meas_start,
            meas_end,
            len(t_list) + len(cp.t_list) - 1,
        )

        # ---- H-list builder (same as _build_h_for_signal) ----------------
        def build_H():
            return self._build_h_for_signal(t_list, t_meas)

        H_curr_list, t_evolve_list = build_H()

        for iteration in range(self.max_iter):
            # ---- forward simulation ----------------------------------
            results = self._forward_simulation(
                B_curr, t_meas, H_curr_list, t_evolve_list
            )
            p_sim = np.array(
                [result.expect[0][-1] for result in results]
            )
            res = p_meas - p_sim
            history["res"].append(res.copy())

            # ---- Jacobian --------------------------------------------
            if self.use_adjoint:
                J = self._compute_jacobian_adjoint(
                    B_curr, t_meas, results, H_curr_list, t_evolve_list
                )
            else:
                J = self._compute_jacobian_fd(
                    B_curr, t_meas, p_sim, H_curr_list, t_evolve_list
                )

            # ---- LM update step --------------------------------------
            A = (
                J.T @ J
                + mu * np.eye(n_basis)
                + reg * regularization_matrix(n_basis, basis_type)
            )
            delta_b = np.linalg.solve(A, J.T @ res)

            b_trial = b + delta_b
            B_trial = B_curr.copy()
            B_trial.update_signal(b=b_trial)

            qubit.qubit_in_mag(B_trial)
            H_trial_list, _ = build_H()
            result_trial = self._forward_simulation(
                B_trial, t_meas, H_trial_list, t_evolve_list
            )
            p_sim_trial = np.array(
                [result.expect[0][-1] for result in result_trial]
            )
            res_trial = p_meas - p_sim_trial

            if np.linalg.norm(res) > np.linalg.norm(res_trial):
                # Accept step
                b = b_trial
                B_curr = B_trial
                H_curr_list = H_trial_list
                mu = max(mu / 2, 1e-8)
                history["b"].append(b.copy())
                history["mu"].append(mu)

                if (
                    np.linalg.norm(delta_b) / (np.linalg.norm(b) + 1e-8)
                    < self.tol
                ):
                    print(f"Converged at iteration {iteration}")
                    break
            else:
                # Reject step, increase damping
                mu = min(mu * 2, 1e8)
                qubit.qubit_in_mag(B_curr)
                H_curr_list, _ = build_H()
                history["b"].append(b.copy())
                history["mu"].append(mu)

            # Print iteration info
            print(
                f"Iter {iteration}: "
                f"res={np.linalg.norm(res):.6f}, "
                f"res_trial={np.linalg.norm(res_trial):.6f}, "
                f"mu={mu:.6e}"
            )

        return b, history
