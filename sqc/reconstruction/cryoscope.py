"""sqc.reconstruction.cryoscope — Cryoscope reconstruction + calibration.

  - CryoscopeReconstruction: φ(t_d) → h(t) waveform reconstruction
  - CryoscopeCalibration:   φ(h) lookup table via square-pulse scan + IQ readout

See Gao 2021 §V for physical model.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from sqc.calibration.base import Calibration, CalibrationTable
from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.hardware.readout import IQReadoutModel
from sqc.reconstruction.base import Reconstruction
from sqc.reconstruction.dispersion import qubit_inverse_frequency


# ===================================================================
# CryoscopeReconstruction
# ===================================================================

@dataclass
class CryoscopeReconstruction(Reconstruction):
    """Cryoscope waveform reconstruction from phase-vs-truncation data.

    Core physics:  dφ/dt = 2π·Δf(h(t))

    Two inversion strategies:
    - ``"calibration"``: use φ(h) lookup table via CalibrationTable.inverse().
    - ``"response"``: use analytical Transmon frequency-flux dispersion.

    Parameters
    ----------
    tau : float
        Calibration square-pulse length (ns).
    inversion : str
        ``"calibration"`` or ``"response"``.
    calibration : CalibrationTable or None
        Required for inversion="calibration".
    qubit : TransmonQubit or None
        Required for inversion="response".
    use_sg_filter : bool
        Savitzky-Golay pre-smoothing on phase derivative.
    sg_window : int
        SG window length (odd). Only used when use_sg_filter=True.
    sg_poly : int
        SG polynomial order.
    """

    tau: float = field(
        default_factory=lambda: CONFIG.reconstruction.cryoscope_tau
    )
    inversion: Literal["calibration", "response"] = "calibration"
    calibration: CalibrationTable | None = None
    qubit: object | None = None
    use_sg_filter: bool = False
    sg_window: int = 7
    sg_poly: int = 2

    def reconstruct(self, measurement, kernel=None, calibration=None,
                    dt: float | None = None) -> FluxSignal:
        varphi = np.asarray(measurement.data["varphi"], dtype=float)
        trunc = np.asarray(measurement.axes["trunc"], dtype=float)

        if len(trunc) > 1 and trunc[0] > trunc[-1]:
            trunc = trunc[::-1]

        if dt is None:
            dt = float(trunc[1] - trunc[0])

        if self.use_sg_filter:
            from scipy.signal import savgol_filter
            window = self.sg_window
            if window >= len(varphi):
                window = max(3, len(varphi) // 2 * 2 - 1)
            dphi_dt = savgol_filter(
                varphi, window_length=window, polyorder=self.sg_poly,
                deriv=1, delta=dt,
            )
        else:
            dphi_dt = np.gradient(varphi, trunc)

        if self.inversion == "calibration":
            cal = calibration if calibration is not None else self.calibration
            if cal is None:
                raise ValueError(
                    "inversion='calibration' requires a CalibrationTable."
                )
            phi_equiv = dphi_dt * self.tau
            h_recon = cal.inverse(phi_equiv)
        elif self.inversion == "response":
            if self.qubit is None:
                raise ValueError(
                    "inversion='response' requires qubit=... (TransmonQubit)."
                )
            h_recon = qubit_inverse_frequency(dphi_dt, self.qubit)
        else:
            raise ValueError(f"Unknown inversion: {self.inversion}")

        return FluxSignal(type=8, t_list=trunc, signal=h_recon)


# ===================================================================
# CryoscopeCalibration
# ===================================================================

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
                "method": "cryoscope", "tau": self.tau,
                "omega_d": float(omega_d),
            },
            metadata={},
        )

    def _unwrap_with_model(self, varphi_raw: np.ndarray, omega_d: float) -> np.ndarray:
        EC = self.qubit.EC
        EJ0 = getattr(self.qubit, "EJ_0", self.qubit.EJ)
        flux_bias = (
            getattr(self.qubit, "flux_bias", None)
            or getattr(self.qubit, "flux", None)
            or 0.0
        )
        total_flux = flux_bias + np.asarray(self.h_list, dtype=float)
        omega_q = np.sqrt(8.0 * EJ0 * np.abs(np.cos(np.pi * total_flux)) * EC) - EC
        varphi_theory = (omega_q - omega_d) * self.tau
        n_wraps = np.round((varphi_theory - varphi_raw) / (2.0 * np.pi))
        return varphi_raw + 2.0 * np.pi * n_wraps
