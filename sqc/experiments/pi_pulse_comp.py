"""sqc.experiments.pi_pulse_comp — PiPulseCompensationExperiment.

Pi-pulse compensation: measure flux pulse tails by applying compensation
flux pulses of varying height z at varying delays tau, combined with a
pi-pulse.  When z exactly compensates the tail, the qubit is on resonance
and flips to |1> (P_e = 1).

Protocol:
  1. Apply a square wave flux pulse at flux-sensitive bias point.
  2. At delay tau after falling edge, apply compensation flux pulse
     (width T_pi, height z) simultaneously with a pi-pulse.
  3. 2D scan over (tau, z) → P_e(tau, z).
  4. Extract z*(tau) = argmax_z P_e(tau, z) — the optimal compensation.
  5. Tail waveform ≈ -z*(tau) (with bias offset).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_pulse
from sqc.config import CONFIG
from sqc.experiments.base import Experiment
from sqc.simulation.result import ExperimentResult


def _default_tail_signal():
    """Create default square-wave + exponential tail flux signal."""
    t_list = CONFIG.pulse.make_time(0, 200)
    signal = np.zeros(len(t_list))
    signal[(t_list >= 10) & (t_list <= 50)] = 0.01  # square pulse
    tail_mask = t_list > 50
    signal[tail_mask] = 0.01 * np.exp(-(t_list[tail_mask] - 50) / 30)
    return FluxSignal(type=8, t_list=t_list, signal=signal)


@dataclass
class PiPulseCompensationExperiment(Experiment):
    """Pi-pulse compensation for measuring flux pulse tails.

    2D scan over (tau, z): for each delay tau after the flux pulse
    falling edge, apply a compensation flux of height z together with
    a pi-pulse.  The optimal z*(tau) that maximizes P_e gives the
    tail waveform.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit at flux-sensitive bias point.
    flux_signal : FluxSignal or None
        Flux signal with tail to measure. If None, creates default.
    tau_list : np.ndarray or None
        Delay times (ns). Default arange(0, 100, dt).
    z_list : np.ndarray or None
        Compensation heights (Phi_0). Default linspace(-0.01, 0.01, 21).
    T_pi : float
        Pi-pulse width (ns). Default from CONFIG.
    omega_bias : float or None
        Qubit frequency at bias point (rad*GHz). If None, uses
        qubit.frequency (caller must ensure correct bias).
    t_rabi : np.ndarray
        Pi-pulse time axis (ns). Default from CONFIG.
    """

    qubit: object
    flux_signal: FluxSignal | None = None
    tau_list: np.ndarray | None = None
    z_list: np.ndarray | None = None
    T_pi: float = field(
        default_factory=lambda: CONFIG.reconstruction.pi_pulse_T_pi
    )
    omega_bias: float | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )

    def __post_init__(self):
        if self.omega_bias is None:
            self.omega_bias = self.qubit.frequency
        if self.flux_signal is None:
            self.flux_signal = _default_tail_signal()
        if self.tau_list is None:
            self.tau_list = CONFIG.pulse.make_time(0, 100)
        if self.z_list is None:
            self.z_list = np.linspace(-0.01, 0.01, 21)

    def build_sequence(self):
        """Sequence is constructed per-point in run()."""
        return None

    def run(self) -> ExperimentResult:
        """Execute pi-pulse compensation experiment.

        Returns
        -------
        ExperimentResult
            data["p_e"]    : 2D array (n_tau, n_z) of P_e values
            data["z_star"] : optimal z*(tau) — the tail waveform
            axes["tau"]    : delay times
            axes["z"]      : compensation heights
        """
        n_tau = len(self.tau_list)
        n_z = len(self.z_list)

        # Total simulation time: pi-pulse + buffer
        t_total = self.t_rabi[-1] + self.T_pi + self.t_rabi[-1]
        t_sig = CONFIG.pulse.make_time(0, t_total)
        t_global = CONFIG.pulse.t_global.copy()

        psi_e = basis(self.qubit.n_levels, 1)

        p_e_2d = np.zeros((n_tau, n_z), dtype=float)

        for i, tau in enumerate(self.tau_list):
            for j, z in enumerate(self.z_list):
                # Build composite flux: tail + compensation, but ZERO
                # outside the pi-pulse window to keep qubit on resonance
                # during pulse edges. Flux there would detune the π-pulse.
                signal = np.zeros(len(t_sig), dtype=float)
                # Pi-pulse window
                pulse_mask = (
                    (t_sig >= self.t_rabi[-1])
                    & (t_sig <= self.t_rabi[-1] + self.T_pi)
                )
                # Tail contribution during window (only where pi-pulse acts)
                tail = np.array(
                    [self.flux_signal.value_at(tau + float(t))
                     for t in t_sig[pulse_mask]],
                    dtype=float,
                )
                # Compensation + tail during window; zero elsewhere
                signal[pulse_mask] = tail + float(z)

                phi_composite = FluxSignal(
                    type=8, t_list=t_sig, signal=signal,
                )

                self.qubit.qubit_in_mag(
                    phi_composite, frame=1, omega_d=self.omega_bias,
                )

                # Pi-pulse Hamiltonian
                H_pi = create_pulse(
                    self.qubit,
                    frame=1,
                    type=1,
                    t_list=np.asarray(self.t_rabi, dtype=float),
                    omega_d=self.omega_bias,
                    phase=0.0,
                    angle=np.pi,
                )

                H_total = (
                    QobjEvo(
                        self.qubit.H_list,
                        tlist=self.qubit.mag_signal.t_list,
                        order=1,
                    )
                    + H_pi
                )

                result = mesolve(
                    H_total,
                    self.qubit.state,
                    t_global,
                    [],
                    e_ops=[psi_e * psi_e.dag()],
                )
                p_e_2d[i, j] = float(result.expect[0][-1])

        # Extract optimal z for each tau
        z_star = np.asarray(
            [self.z_list[int(np.argmax(p_e_2d[i]))] for i in range(n_tau)],
            dtype=float,
        )

        return ExperimentResult(
            data={
                "p_e": p_e_2d,
                "z_star": z_star,
            },
            axes={
                "tau": np.asarray(self.tau_list, dtype=float),
                "z": np.asarray(self.z_list, dtype=float),
            },
            metadata={
                "experiment": "PiPulseCompensationExperiment",
                "T_pi": self.T_pi,
                "omega_bias": self.omega_bias,
            },
            config={
                "T_pi": self.T_pi,
                "t_rabi": self.t_rabi.copy(),
                "omega_bias": self.omega_bias,
            },
        )
