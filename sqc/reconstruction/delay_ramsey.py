"""sqc.reconstruction.delay_ramsey — delay Ramsey reconstruction + calibration.

  - DelayRamseyReconstruction: φ(t_d) → Φ_tail(t_d) waveform reconstruction
  - DelayRamseyCalibration:   φ_cal(z) slope calibration for tail flux conversion
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
from sqc.reconstruction.dispersion import (
    cryoscope_phase_theory,
    qubit_inverse_frequency,
    unwrap_phase_with_model,
)


# ===================================================================
# DelayRamseyReconstruction
# ===================================================================

@dataclass
class DelayRamseyReconstruction(Reconstruction):
    """Delay Ramsey tail-flux waveform reconstruction.

    Parameters
    ----------
    inversion : str
        ``"response"`` — analytical Transmon dispersion (needs qubit).
        ``"calibration"`` — φ(z) lookup table via CalibrationTable.inverse().
    qubit : TransmonQubit or None
        Required for inversion="response".
    calibration : CalibrationTable or None
        Required for inversion="calibration" (kind="phi_z").
    tau_R : float or None
        Ramsey free evolution time (ns). Read from measurement config
        if None.  Default from CONFIG.
    """

    inversion: Literal["response", "calibration"] = "response"
    qubit: object | None = None
    calibration: CalibrationTable | None = None
    tau_R: float | None = field(
        default_factory=lambda: CONFIG.reconstruction.delay_ramsey_tau
    )

    def reconstruct(self, measurement, kernel=None, calibration=None,
                    dt: float | None = None) -> FluxSignal:
        cal = calibration if calibration is not None else self.calibration
        varphi = np.asarray(measurement.data["varphi"], dtype=float)
        t_axis = np.asarray(measurement.axes["t_d"], dtype=float)
        tau = float(measurement.config.get("tau_R", self.tau_R))

        match self.inversion:
            case "response":
                return self._via_response(varphi, t_axis, tau)
            case "calibration":
                return self._via_calibration(varphi, t_axis, cal)
            case _:
                raise ValueError(f"Unknown inversion: {self.inversion}")

    def _via_response(self, varphi, t_axis, tau) -> FluxSignal:
        if self.qubit is None:
            raise ValueError("inversion='response' requires qubit=...")
        h = qubit_inverse_frequency(
            np.asarray(varphi, dtype=float) / tau, self.qubit,
        )
        return FluxSignal(type=8, t_list=t_axis, signal=h)

    def _via_calibration(self, varphi, t_axis, cal) -> FluxSignal:
        if cal is None:
            raise ValueError("inversion='calibration' requires a CalibrationTable.")
        # Upstream now uses model-guided unwrap so phi and cal.outputs
        # share the same absolute-phase convention — no 2π adjustment.
        flux = cal.inverse(np.asarray(varphi, dtype=float))
        return FluxSignal(type=8, t_list=t_axis, signal=flux)


# ===================================================================
# DelayRamseyCalibration
# ===================================================================

@dataclass
class DelayRamseyCalibration(Calibration):
    """Calibrate φ_cal(z) = tau_R * kappa * z for delay Ramsey.

    Scans known flux heights z, measures Ramsey phase at each via IQ
    readout, and fits the linear slope k used to convert tail phase to flux.

    Parameters
    ----------
    qubit : TransmonQubit
    z_list : np.ndarray or None
        Flux heights to scan (Φ₀). Default linspace(-0.02, 0.02, 41).
    tau_R : float
        Ramsey free evolution time (ns). Default from CONFIG.
    t_rabi : np.ndarray
        π/2 pulse time axis (ns). Default from CONFIG.
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
        t_total = 2 * self.t_rabi[-1] + self.tau_R
        t_sig = CONFIG.pulse.make_time(0, t_total)

        readout = IQReadoutModel(
            tau=self.tau_R, t_rabi=self.t_rabi, omega_d=self.omega_d,
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
                0.5 - measured["p_e_I"], measured["p_e_Q"] - 0.5
            ))
            varphi_list.append(phi)

        # Model-guided unwrap shared with DelayRamseyExperiment so the
        # absolute-phase convention is consistent across the two paths.
        varphi_raw = np.asarray(varphi_list, dtype=float)
        varphi_theory = cryoscope_phase_theory(
            self.qubit, np.asarray(self.z_list, dtype=float),
            tau=self.tau_R, omega_d=self.omega_d,
        )
        varphi_arr = unwrap_phase_with_model(varphi_raw, varphi_theory)
        k, intercept = np.polyfit(self.z_list, varphi_arr, 1)

        return CalibrationTable(
            name="delay_ramsey_phi_z",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="phi_z",
            inputs=np.asarray(self.z_list, dtype=float),
            outputs=varphi_arr,
            fit_params={
                "k": float(k), "intercept": float(intercept),
                "tau_R": self.tau_R,
            },
            metadata={"method": "delay_ramsey"},
        )
