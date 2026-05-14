"""sqc.reconstruction.tail — TailReconstruction.

Reconstruct tail flux waveform from delay Ramsey or pi-pulse compensation
measurement data.

Two inversion strategies for delay Ramsey:
- "response": use analytical Transmon dispersion (needs qubit). No calibration.
- "calibration": use phi_z calibration table (needs calibration).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from sqc.calibration.base import CalibrationTable
from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.reconstruction.base import Reconstruction


def _qubit_inverse_frequency(dphi_dt: np.ndarray, qubit) -> np.ndarray:
    """Map frequency shift Δω → flux h via analytical Transmon dispersion.

    f_Q(Φ) = sqrt(8·EJ·|cos(π·Φ)|·EC) - EC   (ħ=1)
    Inverting: cos(π·Φ) = (f_Q + EC)² / (8·EJ·EC)
               Φ = arccos(clip(ratio, 0, 1)) / π
               h = Φ - flux_bias

    Parameters
    ----------
    dphi_dt : np.ndarray
        Phase derivative / frequency shift (rad/ns), equals Δω in natural units.
    qubit
        TransmonQubit providing frequency, EC, EJ, and flux.

    Returns
    -------
    np.ndarray
        Flux offset h (Φ₀) corresponding to dphi_dt.
    """
    f_q = qubit.frequency
    f_target = f_q + dphi_dt
    EC = qubit.EC
    EJ0 = getattr(qubit, "EJ_0", qubit.EJ)
    ratio = (f_target + EC) ** 2 / (8.0 * EC * EJ0)
    ratio = np.clip(ratio, 0.0, 1.0)
    total_flux = np.arccos(ratio) / np.pi
    bias = getattr(qubit, "flux_bias", getattr(qubit, "flux", 0.0))
    return total_flux - bias


@dataclass
class TailReconstruction(Reconstruction):
    """Reconstruct flux tail waveform from delay Ramsey or pi-pulse data.

    Parameters
    ----------
    method : str
        "delay_ramsey" — reconstruct from delay Ramsey phase data.
        "pi_pulse_comp" — negate the optimal z* values.
    inversion : str
        For delay_ramsey:
        - "response" — use analytical qubit dispersion. Requires qubit.
        - "calibration" — use phi_z calibration table. Requires calibration.
    qubit : TransmonQubit or None
        Required for inversion="response".
    calibration : CalibrationTable or None
        Required for inversion="calibration" (kind="phi_z").
    tau_R : float or None
        Ramsey free evolution time (ns). Read from measurement config
        if None. Default from CONFIG.
    """

    method: Literal["delay_ramsey", "pi_pulse_comp"] = "delay_ramsey"
    inversion: Literal["response", "calibration"] = "response"
    qubit: object | None = None
    calibration: CalibrationTable | None = None
    tau_R: float | None = field(
        default_factory=lambda: CONFIG.reconstruction.delay_ramsey_tau
    )

    def reconstruct(
        self,
        measurement,
        kernel=None,
        calibration=None,
        dt: float | None = None,
    ) -> FluxSignal:
        """Reconstruct tail flux waveform.

        Parameters
        ----------
        measurement : ExperimentResult
            delay_ramsey: data["varphi"], axes["t_d"]
            pi_pulse_comp: data["z_star"], axes["tau"]
        calibration : CalibrationTable, optional
            Override instance calibration table.
        kernel : optional
            Not used.
        dt : float, optional
            Not used.

        Returns
        -------
        FluxSignal
            Reconstructed tail flux waveform (type=8).
        """
        cal = calibration if calibration is not None else self.calibration

        match self.method:
            case "delay_ramsey":
                return self._reconstruct_delay_ramsey(measurement, cal)
            case "pi_pulse_comp":
                return self._reconstruct_pi_pulse_comp(measurement)
            case _:
                raise ValueError(f"Unknown method: {self.method}")

    # ------------------------------------------------------------------
    # delay Ramsey
    # ------------------------------------------------------------------

    def _reconstruct_delay_ramsey(
        self, measurement, cal: CalibrationTable | None
    ) -> FluxSignal:
        """Reconstruct tail flux from delay Ramsey phase data."""
        varphi = np.asarray(measurement.data["varphi"], dtype=float)
        t_axis = np.asarray(measurement.axes["t_d"], dtype=float)

        # Use tau_R from measurement config if available, else instance default
        tau = float(measurement.config.get("tau_R", self.tau_R))

        match self.inversion:
            case "response":
                return self._via_response(varphi, t_axis, tau)
            case "calibration":
                return self._via_calibration(varphi, t_axis, cal)
            case _:
                raise ValueError(f"Unknown inversion: {self.inversion}")

    def _via_response(
        self, varphi: np.ndarray, t_axis: np.ndarray, tau: float
    ) -> FluxSignal:
        """Φ ≈ f_Q⁻¹(f_q + φ/τ) - Φ_bias (analytical inversion).

        Shifts varphi by -2π when the unwrap placed it on a positive
        branch (should be negative for κ < 0 and positive Φ).
        """
        if self.qubit is None:
            raise ValueError(
                "inversion='response' requires qubit=... (TransmonQubit)."
            )
        phi = np.asarray(varphi, dtype=float).copy()
        # κ < 0 near typical bias (f ↓ as Φ ↑), so positive h → negative φ.
        # If unwrap placed varphi on the wrong (positive) branch, fix it.
        if np.median(phi) > 0:
            phi = phi - 2.0 * np.pi
        dphi_dt = phi / tau
        h = _qubit_inverse_frequency(dphi_dt, self.qubit)
        return FluxSignal(type=8, t_list=t_axis, signal=h)

    def _via_calibration(
        self, varphi: np.ndarray, t_axis: np.ndarray,
        cal: CalibrationTable | None,
    ) -> FluxSignal:
        """Φ via calibration interpolation.

        Aligns varphi to calibration range by 2π shifts, then uses
        cal.inverse() which interpolates φ → z.
        """
        if cal is None:
            raise ValueError(
                "inversion='calibration' requires a CalibrationTable."
            )
        phi = np.asarray(varphi, dtype=float).copy()
        cal_range = cal.outputs.max() - cal.outputs.min()
        # Align experiment φ median to calibration φ median via 2π shifts
        cal_center = 0.5 * (cal.outputs.min() + cal.outputs.max())
        shift = cal_center - np.median(phi)
        n2pi = np.round(shift / (2.0 * np.pi))
        phi = phi + n2pi * (2.0 * np.pi)
        flux = cal.inverse(phi)
        return FluxSignal(type=8, t_list=t_axis, signal=flux)

    # ------------------------------------------------------------------
    # pi-pulse compensation
    # ------------------------------------------------------------------

    def _reconstruct_pi_pulse_comp(self, measurement) -> FluxSignal:
        """Φ_tail(tau) = -z*(tau)."""
        z_star = np.asarray(measurement.data["z_star"], dtype=float)
        t_axis = np.asarray(measurement.axes["tau"], dtype=float)
        flux = -z_star
        return FluxSignal(type=8, t_list=t_axis, signal=flux)
