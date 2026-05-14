"""sqc.reconstruction.delay_ramsey_calib — DelayRamseyCalibration.

Calibrate phase-vs-flux-height φ_cal(z) for delay Ramsey tail measurement.
Uses IQ readout (two Ramsey sequences with π/2 phase offset)
to extract signed phase via arctan2, avoiding the arccos
sign-ambiguity problem.

The slope k = dφ/dz = tau_R * kappa is used by TailReconstruction
to convert tail phase measurements to flux units.

Moved from sqc.calibration.delay_ramsey (P5 refactor).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.calibration.base import Calibration, CalibrationTable
from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.hardware.readout import IQReadoutModel


@dataclass
class DelayRamseyCalibration(Calibration):
    """Calibrate φ_cal(z) = tau_R * kappa * z for delay Ramsey.

    Scans known flux heights z, measures Ramsey phase at each,
    and fits the linear slope k used to convert tail phase to flux.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit at flux-sensitive bias point.
    z_list : np.ndarray or None
        Flux heights to scan (Phi_0). Default linspace(-0.02, 0.02, 41).
    tau_R : float
        Ramsey free evolution time (ns). Default from CONFIG.
    t_rabi : np.ndarray
        Pi/2 pulse time axis (ns). Default from CONFIG.
    omega_d : float or None
        Drive frequency. If None, uses qubit.frequency.
    """

    qubit: object
    z_list: np.ndarray | None = None
    tau_R: float = field(
        default_factory=lambda: CONFIG.reconstruction.delay_ramsey_tau
    )
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    omega_d: float | None = None

    def __post_init__(self):
        if self.omega_d is None:
            self.omega_d = self.qubit.frequency
        if self.z_list is None:
            self.z_list = np.linspace(-0.02, 0.02, 41)

    def calibrate(self) -> CalibrationTable:
        """Scan z, measure Ramsey phase φ via IQ readout, fit slope k.

        Returns
        -------
        CalibrationTable
            kind="phi_z", inputs=z_list, outputs=varphi_list,
            fit_params={"k": slope, "tau_R": tau_R}
        """
        t_total = 2 * self.t_rabi[-1] + self.tau_R
        t_sig = CONFIG.pulse.make_time(0, t_total)

        readout = IQReadoutModel(
            tau=self.tau_R,
            t_rabi=self.t_rabi,
            omega_d=self.omega_d,
        )

        varphi_list: list[float] = []

        for z in self.z_list:
            signal = np.zeros(len(t_sig))
            evo_start = self.t_rabi[-1]
            evo_end = self.t_rabi[-1] + self.tau_R
            signal[(t_sig >= evo_start) & (t_sig <= evo_end)] = float(z)

            Phi = FluxSignal(type=8, t_list=t_sig, signal=signal)
            self.qubit.qubit_in_mag(Phi, frame=1, omega_d=self.omega_d)

            measured = readout.measure(self.qubit)
            phi = float(np.arctan2(
                0.5 - measured["p_e_I"], measured["p_e_Q"] - 0.5,
            ))
            varphi_list.append(phi)

        varphi_arr = np.asarray(varphi_list, dtype=float)
        varphi_arr = np.unwrap(varphi_arr)

        k, intercept = np.polyfit(self.z_list, varphi_arr, 1)

        return CalibrationTable(
            name="delay_ramsey_phi_z",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="phi_z",
            inputs=np.asarray(self.z_list, dtype=float),
            outputs=varphi_arr,
            fit_params={
                "k": float(k),
                "intercept": float(intercept),
                "tau_R": self.tau_R,
            },
            metadata={"method": "delay_ramsey"},
        )
