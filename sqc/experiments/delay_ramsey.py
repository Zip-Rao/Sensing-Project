"""sqc.experiments.delay_ramsey — DelayRamseyExperiment.

Delay Ramsey (Ramsey tomography): measure flux pulse tails by sliding
a short Ramsey sequence across the falling edge and reading phase via
IQ demodulation.

Protocol:
  1. Apply a square wave flux pulse with tail (distorted).
  2. At delay t_d after falling edge, run short Ramsey (tau_R).
  3. IQ readout extracts phase φ(t_d).
  4. Subtract baseline phase (no-flux reference).
  5. Reconstruct tail flux via φ → Φ mapping (calibration slope k).
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


def _default_tail_signal():
    """Create default square-wave + exponential tail flux signal."""
    t_list = CONFIG.pulse.make_time(0, 200)
    # Square wave from t=10 to t=50, then exponential tail
    signal = np.zeros(len(t_list))
    signal[(t_list >= 10) & (t_list <= 50)] = 0.01  # square pulse
    # exponential tail after t=50
    tail_mask = t_list > 50
    signal[tail_mask] = 0.01 * np.exp(-(t_list[tail_mask] - 50) / 30)
    return FluxSignal(type=8, t_list=t_list, signal=signal)


@dataclass
class DelayRamseyExperiment(Experiment):
    """Delay Ramsey for measuring flux pulse tails.

    Slides a short Ramsey sequence across the flux signal's tail,
    measuring the accumulated phase at each delay via IQ readout.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit at flux-sensitive bias point.
    flux_signal : FluxSignal or None
        Flux signal with tail to measure. If None, creates a default
        square-wave + exponential-tail signal.
    t_d_list : np.ndarray or None
        Delay times after falling edge (ns). Default arange(0, 200, dt).
    tau_R : float
        Ramsey free evolution time (ns). Default from CONFIG.
    t_rabi : np.ndarray
        Pi/2 pulse time axis (ns). Default from CONFIG.
    omega_d : float or None
        Drive frequency. If None, uses qubit.frequency.
    run_baseline : bool
        If True, subtract baseline phase (no-flux reference). Default True.
    """

    qubit: object
    flux_signal: FluxSignal | None = None
    t_d_list: np.ndarray | None = None
    tau_R: float = field(
        default_factory=lambda: CONFIG.reconstruction.delay_ramsey_tau
    )
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    omega_d: float | None = None

    # -- P9.B --
    control_line: object | None = None
    run_baseline: bool = False  # Deprecated.  Model-guided unwrap
                                # (via sqc.reconstruction.dispersion)
                                # now sets the absolute-phase reference
                                # analytically — no baseline measurement
                                # is needed.  Kept for API compatibility.
    t_fall: float = 0.0
    """Falling-edge time in the flux signal (ns). t_d is relative to this."""

    def __post_init__(self):
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency
        if self.flux_signal is None:
            self.flux_signal = _default_tail_signal()
        if self.t_d_list is None:
            self.t_d_list = CONFIG.pulse.make_time(0, 200)

    def build_sequence(self):
        """Sequence is constructed per-delay in run()."""
        return None

    def run(self) -> ExperimentResult:
        """Execute delay Ramsey experiment.

        Returns
        -------
        ExperimentResult
            data["varphi"]      : baseline-subtracted phase at each t_d
            data["varphi_raw"]  : raw phase before baseline subtraction
            data["varphi_base"] : baseline phase (None if run_baseline=False)
            data["p_e_I"]       : I-channel probabilities
            data["p_e_Q"]       : Q-channel probabilities
            axes["t_d"]         : delay times
        """
        t_total = 2 * self.t_rabi[-1] + self.tau_R
        t_sig = CONFIG.pulse.make_time(0, t_total)

        readout = IQReadoutModel(
            tau=self.tau_R,
            t_rabi=self.t_rabi,
            omega_d=self.omega_d,
        )

        p_e_I_list: list[float] = []
        p_e_Q_list: list[float] = []

        for t_d in self.t_d_list:
            # Window the tail segment centering at t_d, but ZERO outside(i.e., no flux during π/2 pulses) to avoid detuning and corrupting the Ramsey sequence. The tail segment is defined as the flux signal during
            # the free-evolution window [t_rabi[-1], t_rabi[-1]+tau_R].
            # Flux during π/2 pulses would detune them and corrupt φ.
            signal = np.zeros(len(t_sig), dtype=float)
            free_mask = (
                (t_sig >= self.t_rabi[-1])
                & (t_sig <= self.t_rabi[-1] + self.tau_R)
            )
            signal[free_mask] = np.array(
                [self.flux_signal.value_at(self.t_fall + t_d + float(t) - self.t_rabi[-1] - self.tau_R/2)
                 for t in t_sig[free_mask]],
                dtype=float,
            )
            phi_windowed = FluxSignal(type=8, t_list=t_sig, signal=signal)
            phi_windowed = self._route_flux(phi_windowed)

            self.qubit.qubit_in_mag(
                phi_windowed, frame=1, omega_d=self.omega_d,
            )

            measured = readout.measure(self.qubit)
            p_e_I_list.append(measured["p_e_I"])
            p_e_Q_list.append(measured["p_e_Q"])

        p_e_I = np.asarray(p_e_I_list, dtype=float)
        p_e_Q = np.asarray(p_e_Q_list, dtype=float)
        t_d_arr = np.asarray(self.t_d_list, dtype=float)

        # Raw wrapped phase
        varphi_raw = np.arctan2(0.5 - p_e_I, p_e_Q - 0.5)

        # Model-guided unwrap shared with DelayRamseyCalibration so the
        # absolute-phase convention is consistent across the two paths.
        # Theory anchor: ∫_{t_d}^{t_d+τ_R} (ω_q(Φ_bias + h_tail(s)) − ω_d) ds
        # where h_tail(s) = flux_signal.value_at(t_fall + s).
        varphi_theory = cumulative_phase_theory(
            self.qubit,
            np.asarray(self.flux_signal.t_list, dtype=float) - self.t_fall,
            np.asarray(self.flux_signal.signal, dtype=float),
            t_d_arr, omega_d=self.omega_d,
            window=(0.0, float(self.tau_R)),
        )
        varphi = unwrap_phase_with_model(varphi_raw, varphi_theory)

        # Keep legacy keys for downstream compatibility.
        varphi_base: float | None = None

        return ExperimentResult(
            data={
                "varphi": varphi,
                "varphi_raw": varphi_raw,
                "varphi_base": varphi_base,
                "p_e_I": p_e_I,
                "p_e_Q": p_e_Q,
            },
            axes={
                "t_d": np.asarray(self.t_d_list, dtype=float),
            },
            metadata={
                "experiment": "DelayRamseyExperiment",
                "tau_R": self.tau_R,
            },
            config={
                "tau_R": self.tau_R,
                "t_rabi": self.t_rabi.copy(),
                "omega_d": self.omega_d,
                "run_baseline": self.run_baseline,
            },
        )
