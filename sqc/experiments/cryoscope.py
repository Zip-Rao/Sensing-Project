"""sqc.experiments.cryoscope — CryoscopeExperiment.

Cryoscope: scan truncation delay, measure φ(t_d) via IQ Ramsey.
Replaces Protocal.evolve case 5.

Default parameters match src/protocal.py case 5 exactly:
  - Phi: type=2, t_list=linspace(0, 80, 160), amplitude=0.01
  - trunc_list: Phi.t_list[140:20:-1] (reverse order)
  - tau: 100.0
  - t_rabi: linspace(0, 10, 20)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.config import CONFIG
from sqc.experiments.base import Experiment
from sqc.hardware.readout import IQReadoutModel
from sqc.reconstruction.dispersion import (
    cumulative_phase_theory,
    unwrap_phase_with_model,
)
from sqc.simulation.result import ExperimentResult


@dataclass
class CryoscopeExperiment(Experiment):
    """Cryoscope: scan truncation delay, measure φ(t_d) via IQ Ramsey.

    Direct port of src/protocal.py:Protocal.evolve case 5.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit object (src or sqc version).
    flux_signal : FluxSignal or None
        Flux signal to scan. If None, creates default sinusoidal signal
        matching src/protocal.py case 5.
    t_rabi : np.ndarray
        Rabi pulse time axis (ns). Default linspace(0, 10, 20).
    tau : float
        Free precession time for IQ readout (ns). Default 100.0.
    trunc_list : np.ndarray or None
        List of truncation times. If None, defaults to
        flux_signal.t_list[140:20:-1] (reverse order).
    omega_d : float or None
        Drive frequency. If None, uses qubit.frequency.
    """

    qubit: object  # TransmonQubit (duck-typed)
    flux_signal: FluxSignal | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    tau: float = field(
        default_factory=lambda: CONFIG.reconstruction.cryoscope_tau
    )
    trunc_list: np.ndarray | None = None
    omega_d: float | None = None

    # -- P9.B: control-line distortion injection --
    control_line: object | None = None

    def __post_init__(self):
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency
        if self.flux_signal is None:
            # Extended to 100 ns to prevent silent truncation failures
            # when the user lengthens trunc_list beyond the signal window.
            self.flux_signal = FluxSignal(
                type=2,
                t_list=CONFIG.pulse.make_time(0, 100),
                amplitude=0.01,
            )
        if self.trunc_list is None:
            print("Warning: trunc_list not provided, defaulting to flux_signal.t_list[140:20:-")
            self.trunc_list = self.flux_signal.t_list[180:40:-1]

        # -- boundary sanity check -------------------------------------------
        _t_max = float(self.flux_signal.t_list[-1])
        _bad = [t for t in self.trunc_list if t > _t_max]
        if _bad:
            raise ValueError(
                f"trunc_list contains values ({_bad}) beyond "
                f"flux_signal.t_list[-1] ({_t_max} ns). "
                f"Extend flux_signal.t_list or shorten trunc_list."
            )

    def build_sequence(self):
        """Cryoscope sequence is constructed per-truncation in run()."""
        return None

    def run(self) -> ExperimentResult:
        """Execute Cryoscope experiment.

        For each truncation delay, truncates the flux signal, couples
        it to the qubit, and performs IQ readout to measure the
        accumulated phase.

        Returns
        -------
        ExperimentResult
            With data["varphi"], data["p_e_I"], data["p_e_Q"],
            axes["trunc"].
        """
        readout = IQReadoutModel(
            tau=self.tau,
            t_rabi=self.t_rabi,
            omega_d=self.omega_d,
        )
        p_e_I_list: list[float] = []
        p_e_Q_list: list[float] = []
        print("trunc_list:", self.trunc_list)
        # Iterate truncation delays in reverse order (matching legacy)
        for trunc in self.trunc_list:
            # Copy flux signal and truncate in-place (legacy semantics)
            phi_truncated = self.flux_signal.copy()
            phi_truncated.truncate(0, float(trunc))

            # Couple flux to qubit
            phi_truncated = self._route_flux(phi_truncated)
            self.qubit.qubit_in_mag(
                phi_truncated, frame=1, omega_d=self.omega_d,
            )

            # IQ readout
            measured = readout.measure(self.qubit)
            p_e_I_list.append(measured["p_e_I"])
            p_e_Q_list.append(measured["p_e_Q"])

        p_e_I = np.asarray(p_e_I_list, dtype=float)
        p_e_Q = np.asarray(p_e_Q_list, dtype=float)
        trunc_arr = np.asarray(self.trunc_list, dtype=float)

        # Sort so trunc / varphi are time-ascending
        order = np.argsort(trunc_arr)
        trunc_sorted = trunc_arr[order]
        varphi_raw = np.arctan2(0.5 - p_e_I, p_e_Q - 0.5)[order]

        # Model-guided unwrap anchored to ∫₀^{t_d} (ω_q(Φ(t)) − ω_d) dt
        # so the absolute-phase convention matches CryoscopeCalibration.
        varphi_theory = cumulative_phase_theory(
            self.qubit,
            np.asarray(self.flux_signal.t_list, dtype=float),
            np.asarray(self.flux_signal.signal, dtype=float),
            trunc_sorted, omega_d=self.omega_d,
        )
        varphi = unwrap_phase_with_model(varphi_raw, varphi_theory)
        return ExperimentResult(
            data={
                "varphi": varphi,
                "p_e_I": p_e_I,
                "p_e_Q": p_e_Q,
            },
            axes={
                "trunc": trunc_sorted,
            },
            metadata={
                "experiment": "CryoscopeExperiment",
                "tau": self.tau,
            },
            config={
                "tau": self.tau,
                "t_rabi": self.t_rabi.copy(),
                "omega_d": self.omega_d,
            },
        )
