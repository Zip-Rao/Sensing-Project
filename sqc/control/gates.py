"""sqc.control.gates — Two-qubit gate functions.

Verbatim port of src/qubit.py module-level gate functions:
    ideal_iSWAP, simulate_iSWAP, ideal_CZ, simulate_CZ.

Extended for waveform-driven CZ simulation (P7):
    simulate_cz_from_flux — accepts explicit time and on-chip flux arrays.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.interpolate import interp1d
from qutip import Qobj, destroy, propagator, qeye, tensor, mesolve, Options, basis


def ideal_iSWAP() -> Qobj:
    """Ideal iSWAP gate matrix.

    Returns
    -------
    Qobj
        4×4 unitary in |00>,|01>,|10>,|11> basis.
    """
    iSWAP = Qobj(
        [
            [1, 0, 0, 0],
            [0, 0, 1j, 0],
            [0, 1j, 0, 0],
            [0, 0, 0, 1],
        ],
        dims=[[2, 2], [2, 2]],
    )
    return iSWAP


def simulate_iSWAP(qubit1, qubit2, g: float):
    """Simulate iSWAP gate (simple model, no coupler).

    Parameters
    ----------
    qubit1, qubit2 : TransmonQubit
    g : float
        Exchange coupling (GHz).

    Returns
    -------
    U_eff : Qobj
        4×4 effective unitary.
    leakage : float
    F : float
        Average gate fidelity to ideal iSWAP.
    """
    H_0 = tensor(
        qubit1.get_hamiltonian_rwa(qubit1.frequency),
        qeye(qubit2.n_levels),
    ) + tensor(
        qeye(qubit1.n_levels),
        qubit2.get_hamiltonian_rwa(qubit2.frequency),
    )
    H_int = g * (
        tensor(
            destroy(qubit1.n_levels), destroy(qubit2.n_levels).dag()
        )
        + tensor(
            destroy(qubit1.n_levels).dag(), destroy(qubit2.n_levels)
        )
    )
    H = H_0 + H_int
    T = np.pi / (2 * g)
    U = propagator(H, T)
    idx = [
        0, 1, qubit2.n_levels, qubit2.n_levels + 1
    ]
    U_eff = Qobj(
        U.full()[np.ix_(idx, idx)], dims=[[2, 2], [2, 2]]
    )
    U_eff = U_eff.dag()  # TODO: verify convention
    U_ideal = ideal_iSWAP()
    leakage = 1.0 - np.mean(np.sum(np.abs(U_eff.full()) ** 2, axis=0))
    F = (abs((U_ideal.dag() @ U_eff).tr()) ** 2 + 4) / (4 * 5)
    return U_eff, leakage, F


def ideal_CZ() -> Qobj:
    """Ideal CZ gate matrix.

    Returns
    -------
    Qobj
        4×4 unitary.
    """
    CZ = Qobj(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, -1],
        ],
        dims=[[2, 2], [2, 2]],
    )
    return CZ


def _projected_average_gate_fidelity(
    target: np.ndarray, projected_map: np.ndarray,
) -> float:
    """Average target overlap for a possibly trace-decreasing projection."""
    target = np.asarray(target, dtype=complex)
    projected_map = np.asarray(projected_map, dtype=complex)
    if target.shape != projected_map.shape or target.ndim != 2:
        raise ValueError("target and projected_map must be same-size matrices")
    d = target.shape[0]
    survival_trace = float(np.trace(
        projected_map.conj().T @ projected_map
    ).real)
    return float(
        (abs(np.trace(target.conj().T @ projected_map)) ** 2
         + survival_trace) / (d * (d + 1))
    )


def simulate_CZ(qubit1, qubit2, g: float):
    """Simulate CZ gate via flux pulsing.

    Parameters
    ----------
    qubit1, qubit2 : TransmonQubit
    g : float
        Exchange coupling (GHz).

    Returns
    -------
    U_eff_corr : Qobj
        4×4 effective unitary with phase correction.
    leakage : float
    F : float
    """
    H_0 = tensor(
        qubit1.get_hamiltonian_rwa(qubit1.frequency),
        qeye(qubit2.n_levels),
    ) + tensor(
        qeye(qubit1.n_levels),
        qubit2.get_hamiltonian_rwa(qubit2.frequency),
    )
    H_int = g * (
        tensor(
            destroy(qubit1.n_levels), destroy(qubit2.n_levels).dag()
        )
        + tensor(
            destroy(qubit1.n_levels).dag(), destroy(qubit2.n_levels)
        )
    )
    H_1 = H_0 + H_int

    def Delta_1(t, args):
        T, T1, T2 = args["T"], args["T1"], args["T2"]
        A = args["A"]
        if t < 0 or t > T:
            return 0.0
        elif t < T1:
            return A * np.sin(0.5 * np.pi * t / T1) ** 2
        elif t < T - T2:
            return A
        else:
            return A * np.sin(0.5 * np.pi * (T - t) / T2) ** 2

    T = np.pi / (np.sqrt(2.0) * g)
    T1 = T / 4
    T2 = T / 4
    A = -qubit2.anharmonicity
    args = {"T": T, "T1": T1, "T2": T2, "A": A}
    n_1 = tensor(qubit1.n, qeye(qubit2.n_levels))
    H = [H_1, [n_1, Delta_1]]
    U = propagator(H, T, args=args)
    idx = [
        0, 1, qubit2.n_levels, qubit2.n_levels + 1
    ]
    U_eff = Qobj(
        U.full()[np.ix_(idx, idx)], dims=[[2, 2], [2, 2]]
    )
    phi00 = np.angle(U_eff[0, 0])
    phi10 = np.angle(U_eff[1, 1])
    phi01 = np.angle(U_eff[2, 2])
    phi1 = phi10 - phi00
    phi2 = phi01 - phi00
    D = Qobj(
        np.diag(
            [
                np.exp(-1j * phi00),
                np.exp(-1j * phi10),
                np.exp(-1j * phi01),
                np.exp(-1j * (phi10 + phi01 - phi00)),
            ]
        ),
        dims=[[2, 2], [2, 2]],
    )
    U_eff_corr = D @ U_eff
    U_ideal = ideal_CZ()
    leakage = 1.0 - 0.25 * (np.sum(np.abs(U_eff_corr.full()) ** 2))
    F = (abs((U_ideal.dag() @ U_eff_corr).tr()) ** 2 + 4) / (4 * 5)
    return U_eff_corr, leakage, F


# ===================================================================
# waveform-driven CZ simulation (P7 — explicit on-chip flux)
# ===================================================================


@dataclass
class CZResult:
    """Results of a waveform-driven CZ gate simulation.

    Attributes
    ----------
    U_eff : np.ndarray
        4×4 effective unitary in the computational subspace.
    conditional_phase : float
        Conditional phase φ_CZ (rad).
    phase_error : float
        Absolute deviation |φ_CZ − π| (rad).
    leakage : float
        Total leakage out of the computational subspace,
        1 − Tr(P_comp ρ) averaged over computational basis states.
    populations : dict[str, float]
        Key non-computational populations at final time
        (e.g. ``pop_20``, ``pop_02``).
    infidelity : float
        Average gate infidelity (1 − F_avg).
    phase_trajectory : np.ndarray or None
        Conditional phase vs time, if *store_trajectories* is True.
    leakage_trajectory : np.ndarray or None
        Leakage vs time, if *store_trajectories* is True.
    """

    U_eff: np.ndarray
    conditional_phase: float
    phase_error: float
    leakage: float
    populations: dict
    infidelity: float
    phase_trajectory: np.ndarray | None = None
    leakage_trajectory: np.ndarray | None = None


def simulate_cz_from_flux(
    t: np.ndarray,
    phi_chip: np.ndarray,
    qubit1,
    qubit2,
    g: float,
    n_levels: int = 3,
    store_trajectories: bool = False,
    solver_options: Optional[dict] = None,
) -> CZResult:
    """Simulate a CZ gate driven by an explicit on-chip flux waveform.

    Converts ``phi_chip(t)`` → detuning Δ(t) via the qubit dispersion
    relation, builds the time-dependent two-qubit Hamiltonian, and
    computes the full evolution operator.  Gate metrics are extracted
    from the computational-subspace projection.

    Parameters
    ----------
    t : np.ndarray
        Time axis (ns).  Must be uniformly sampled.
    phi_chip : np.ndarray
        On-chip flux at each time point (Φ₀).  Same length as *t*.
    qubit1 : TransmonQubit
        Flux-pulsed qubit (the one whose frequency is modulated).
    qubit2 : TransmonQubit
        Static qubit (idle frequency).
    g : float
        Exchange coupling strength (GHz).
    n_levels : int
        Hilbert-space truncation per qubit.  Default 3.
    store_trajectories : bool
        If True, record phase and leakage vs time (uses mesolve).
    solver_options : dict or None
        Options forwarded to ``qutip.mesolve`` or ``qutip.propagator``.

    Returns
    -------
    CZResult
    """
    from scipy.interpolate import interp1d

    # -- 1. Convert on-chip flux → detuning via qubit dispersion ----------
    flux_bias = float(getattr(qubit1, "flux", 0.0))
    spec = getattr(qubit1, "_spec", None)
    if spec is None:
        raise ValueError("qubit1 must have a _spec (QubitSpec) attribute.")

    # Total flux = static bias + pulse
    phi_total = flux_bias + np.asarray(phi_chip, dtype=float)

    # Compute instantaneous frequency (rad·GHz — same units as Hamiltonian)
    omega_t = np.array([spec.frequency(float(p)) for p in phi_total],
                       dtype=float)
    omega_idle = spec.frequency(flux_bias)

    # Detuning on qubit1 only (rad·GHz)
    detuning = omega_t - omega_idle

    # -- 2. Build the two-qubit Hamiltonian --------------------------------
    # Interpolate detuning for QuTiP function-coefficient format
    detune_interp = interp1d(
        t, detuning, kind="linear",
        bounds_error=False, fill_value=0.0,
    )

    def _detune_func(t_val, args=None):
        return float(detune_interp(t_val))

    # Identity on the other qubit
    id2 = qeye(qubit2.n_levels)
    n1_op = tensor(qubit1.n, id2)

    # Static Hamiltonian
    H_0 = (
        tensor(qubit1.get_hamiltonian_rwa(omega_idle), id2)
        + tensor(qeye(qubit1.n_levels), qubit2.get_hamiltonian_rwa(omega_idle))
    )
    H_int = g * (
        tensor(destroy(qubit1.n_levels), destroy(qubit2.n_levels).dag())
        + tensor(destroy(qubit1.n_levels).dag(), destroy(qubit2.n_levels))
    )
    H_const = H_0 + H_int

    # Time-dependent part: detuning on qubit1
    H = [H_const, [n1_op, _detune_func]]

    T = float(t[-1] - t[0])

    # -- 3. Evolve ---------------------------------------------------------
    opts = dict(solver_options or {})
    if store_trajectories:
        # mesolve path: gives access to intermediate states
        opts.setdefault("store_states", True)
        opts.setdefault("atol", 1e-10)
        opts.setdefault("rtol", 1e-8)
        dt_ev = float(t[1] - t[0])
        opts.setdefault("max_step", dt_ev / 2.0)
        opts.setdefault("nsteps", 10000)

        # Prepare four initial states (computational basis)
        from qutip import tensor as qt_tensor
        psi00 = qt_tensor(basis(qubit1.n_levels, 0), basis(qubit2.n_levels, 0))
        psi01 = qt_tensor(basis(qubit1.n_levels, 0), basis(qubit2.n_levels, 1))
        psi10 = qt_tensor(basis(qubit1.n_levels, 1), basis(qubit2.n_levels, 0))
        psi11 = qt_tensor(basis(qubit1.n_levels, 1), basis(qubit2.n_levels, 1))

        U_full = np.zeros((qubit1.n_levels * qubit2.n_levels,
                           qubit1.n_levels * qubit2.n_levels), dtype=complex)

        from qutip import QobjEvo
        H_evo = QobjEvo(H, tlist=t)
        traj_phase = np.zeros(len(t))
        traj_leak = np.zeros(len(t))

        comp_columns = [0, 1, qubit2.n_levels, qubit2.n_levels + 1]
        for idx, psi0 in enumerate([psi00, psi01, psi10, psi11]):
            result = mesolve(H_evo, psi0, t, [], [], options=Options(**opts))
            U_full[:, comp_columns[idx]] = result.states[-1].full().flatten()

            if idx == 3:  # |11⟩ — track conditional phase and leakage
                for k, state in enumerate(result.states):
                    # Phase: arg(⟨11|U|11⟩) relative to ⟨00|U|00⟩
                    # leakage: 1 - |⟨comp|state⟩|²
                    proj_comp = 0.0
                    for a in range(2):
                        for b_val in range(2):
                            bra = qt_tensor(
                                basis(qubit1.n_levels, a),
                                basis(qubit2.n_levels, b_val),
                            )
                            proj_comp += abs(bra.overlap(state)) ** 2
                    traj_leak[k] = 1.0 - float(proj_comp)
                    # Phase tracking uses arg of ⟨11|state⟩
                    bra11 = qt_tensor(
                        basis(qubit1.n_levels, 1),
                        basis(qubit2.n_levels, 1),
                    )
                    traj_phase[k] = float(np.angle(bra11.overlap(state)))
    else:
        # propagator path: faster, no intermediate states
        U_full = propagator(H, T, options=Options(**opts)).full()
        traj_phase = traj_leak = None

    # -- 4. Extract gate metrics -------------------------------------------
    # Computational subspace indices: |00⟩, |01⟩, |10⟩, |11⟩
    n2 = qubit2.n_levels
    idx = [0, 1, n2, n2 + 1]
    U_eff_arr = U_full[np.ix_(idx, idx)]
    U_eff = Qobj(U_eff_arr, dims=[[2, 2], [2, 2]])

    # Virtual-Z correction (same as simulate_CZ)
    phi00 = np.angle(U_eff_arr[0, 0])
    phi10 = np.angle(U_eff_arr[1, 1])
    phi01 = np.angle(U_eff_arr[2, 2])
    phi11 = np.angle(U_eff_arr[3, 3])

    # Conditional phase
    cond_phase = phi11 - phi01 - phi10 + phi00
    # Wrap to [-π, π]
    cond_phase = float(np.arctan2(np.sin(cond_phase), np.cos(cond_phase)))

    # Apply virtual-Z diagonal correction
    D_arr = np.diag([
        np.exp(-1j * phi00),
        np.exp(-1j * phi10),
        np.exp(-1j * phi01),
        np.exp(-1j * (phi10 + phi01 - phi00)),
    ])
    U_eff_corr_arr = D_arr @ U_eff_arr

    # Leakage: averaged over 4 computational input states
    proj_comp_avg = np.mean(np.sum(np.abs(U_eff_corr_arr) ** 2, axis=0))
    leakage = float(1.0 - proj_comp_avg)

    # Average gate fidelity
    U_ideal = ideal_CZ().full()
    # For a projected (generally non-unitary) computational-subspace map M,
    # the Haar-averaged overlap is
    # (|Tr(U_target^dag M)|^2 + Tr(M^dag M)) / [d(d+1)].
    # Replacing the second term by d is valid only when M is unitary and
    # would otherwise hide part of the leakage error.
    F_avg = _projected_average_gate_fidelity(U_ideal, U_eff_corr_arr)
    infidelity = 1.0 - F_avg

    # Key populations (from the full U matrix, for |11⟩ input)
    col11 = U_full[:, n2 + 1]
    pop_20 = float(abs(col11[2 * n2 + 0]) ** 2)  # |20⟩
    pop_02 = float(abs(col11[0 * n2 + 2]) ** 2)  # |02⟩

    return CZResult(
        U_eff=U_eff_corr_arr,
        conditional_phase=abs(cond_phase),
        phase_error=abs(abs(cond_phase) - np.pi),
        leakage=leakage,
        populations={"pop_20": pop_20, "pop_02": pop_02},
        infidelity=infidelity,
        phase_trajectory=traj_phase,
        leakage_trajectory=traj_leak,
    )
