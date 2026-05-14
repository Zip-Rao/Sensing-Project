"""sqc.reconstruction.cryoscope_calib — CryoscopeCalibration.

Build φ(h) lookup table via square-pulse scan + IQ readout.
Extracted from FluxResponseCalibration._calibrate_cryoscope (P3a).

This is a reconstruction pre-calibration step: the resulting
CalibrationTable(kind="phi_h") feeds into CryoscopeReconstruction.

See Gao 2021 §V for the physical model.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from sqc.calibration.base import Calibration, CalibrationTable
from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.hardware.readout import IQReadoutModel


@dataclass
class CryoscopeCalibration(Calibration):
    """Build φ(h) lookup table via cryoscope method.

    Scans square-pulse heights h, extracts φ(h) using IQ readout
    with model-guided phase unwrapping.

    Parameters
    ----------
    qubit : TransmonQubit
    h_list : np.ndarray or None
        Flux heights to scan (Φ₀). Default linspace(-0.03, 0.03, 51).
    tau : float
        Square-pulse duration (ns). Default 100.0.
    t_rabi : np.ndarray
        π/2 pulse time axis (ns).
    """

    qubit: object
    h_list: np.ndarray | None = None
    tau: float = 100.0
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )

    def __post_init__(self):
        if self.h_list is None:
            self.h_list = np.linspace(-0.03, 0.03, 51)

    def calibrate(self) -> CalibrationTable:
        """Scan square-pulse heights, extract φ(h) via IQ readout.

        Returns
        -------
        CalibrationTable
            kind="phi_h", inputs=h_list, outputs=varphi.
        """
        omega_d = self.qubit.frequency
        t_pi2_end = self.t_rabi[-1]
        t_total = 2 * t_pi2_end + self.tau

        p_e_I_list: list[float] = []
        p_e_Q_list: list[float] = []

        readout = IQReadoutModel(
            tau=self.tau, t_rabi=self.t_rabi, omega_d=omega_d,
        )

        for h in self.h_list:
            t_sig = CONFIG.pulse.make_time(0, t_total)
            signal = np.zeros_like(t_sig, dtype=float)
            mask = (t_sig >= t_pi2_end) & (t_sig <= t_pi2_end + self.tau)
            signal[mask] = float(h)
            Phi = FluxSignal(type=8, t_list=t_sig, signal=signal)

            self.qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

            result = readout.measure(self.qubit)
            p_e_I_list.append(result["p_e_I"])
            p_e_Q_list.append(result["p_e_Q"])

        p_e_I = np.asarray(p_e_I_list, dtype=float)
        p_e_Q = np.asarray(p_e_Q_list, dtype=float)
        varphi_raw = np.arctan2(0.5 - p_e_I, p_e_Q - 0.5)

        varphi = self._unwrap_with_model(varphi_raw, omega_d)

        return CalibrationTable(
            name="cryoscope_phi_h",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="phi_h",
            inputs=np.asarray(self.h_list, dtype=float),
            outputs=np.asarray(varphi, dtype=float),
            fit_params={
                "method": "cryoscope",
                "tau": self.tau,
                "omega_d": float(omega_d),
            },
            metadata={},
        )

    def _unwrap_with_model(
        self, varphi_raw: np.ndarray, omega_d: float,
    ) -> np.ndarray:
        """Resolve 2π ambiguities using analytical Transmon model.

        φ_th(h) = (ω_q(Φ_bias + h) - ω_d) · τ
        """
        EC = self.qubit.EC
        EJ0 = getattr(self.qubit, "EJ_0", self.qubit.EJ)
        flux_bias = (
            getattr(self.qubit, "flux_bias", None)
            or getattr(self.qubit, "flux", None)
            or 0.0
        )

        total_flux = flux_bias + np.asarray(self.h_list, dtype=float)
        omega_q = np.sqrt(8.0 * EJ0 * np.abs(np.cos(np.pi * total_flux)) * EC) - EC
        delta_omega = omega_q - omega_d
        varphi_theory = delta_omega * self.tau

        n_wraps = np.round((varphi_theory - varphi_raw) / (2.0 * np.pi))
        return varphi_raw + 2.0 * np.pi * n_wraps
