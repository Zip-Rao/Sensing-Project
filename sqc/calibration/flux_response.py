"""sqc.calibration.flux_response — FluxResponseCalibration.

Calibrate the flux-frequency or flux-phase response of a Transmon qubit.

Three methods:
  - "ramsey": scan DC flux steps, measure frequency via Ramsey (case 1)
  - "cryoscope": scan square-pulse heights, build φ(h) lookup (case 3)
  - "transient": unknown transient signal → fit Δω(Φ) polynomial (case 2)

Currently "ramsey" and "cryoscope" are implemented. "transient" requires
Track B 1.2 completion (see _TODO_master.md 1.2).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.calibration.base import Calibration, CalibrationTable
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.hardware.readout import IQReadoutModel


@dataclass
class FluxResponseCalibration(Calibration):
    """Calibrate f(Phi) or φ(h) via Ramsey, Cryoscope, or transient sensing.

    Replaces src/protocal.py:Calibration.calibrate cases 1/2/3.

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit to calibrate.
    method : str
        Calibration method:
        - "ramsey": scan DC flux steps, measure frequency via Ramsey.
        - "cryoscope": scan square-pulse heights, build φ(h).
        - "transient": unknown signal polynomial fitting. (Track B 1.2)
    h_list : np.ndarray or None
        Flux/height values to scan. Default linspace(-0.03, 0.03, 21).
    tau : float
        Free precession / square-pulse time (ns). Default 100.0.
    t_rabi : np.ndarray or None
        Rabi pulse time axis (ns). Default linspace(0, 10, 20).
    """

    qubit: object  # TransmonQubit
    method: Literal["ramsey", "cryoscope", "transient"] = "ramsey"
    h_list: np.ndarray | None = None
    tau: float = 100.0
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )

    def __post_init__(self):
        if self.h_list is None:
            self.h_list = np.linspace(-0.03, 0.03, 21)

    def calibrate(self) -> CalibrationTable:
        """Run the calibration and return a CalibrationTable.

        Returns
        -------
        CalibrationTable
            kind="f_phi" for ramsey, kind="phi_h" for cryoscope.

        Raises
        ------
        NotImplementedError
            For method="transient" (requires Track B 1.2).
        """
        match self.method:
            case "ramsey":
                return self._calibrate_ramsey()
            case "cryoscope":
                return self._calibrate_cryoscope()
            case "transient":
                return self._calibrate_transient()

    # ------------------------------------------------------------------
    # Ramsey-based flux response calibration
    # ------------------------------------------------------------------

    def _calibrate_ramsey(self) -> CalibrationTable:
        """Scan flux Φ in static steps, fit Ramsey detuning at each step.

        Builds f(Φ) lookup table.

        For each flux value in h_list:
        1. Create a constant flux signal with amplitude h
        2. Couple to qubit and run Ramsey
        3. Measure p_e(tau) and fit the oscillation frequency
        4. Record f(Φ) = omega_d + fitted_detuning
        """
        omega_d = self.qubit.frequency
        tau_list = CONFIG.pulse.make_time(0, 200)
        frequency_list: list[float] = []

        psi_e = basis(self.qubit.n_levels, 1)
        t_global = CONFIG.pulse.t_global.copy()

        for h in self.h_list:
            # Create constant flux signal
            t_sig = CONFIG.pulse.make_time(0, 300)
            Phi = FluxSignal(
                type=1,
                t_list=t_sig,
                amplitude=float(h),
                offset=0.0,
            )

            # Couple flux to qubit
            self.qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

            # Run Ramsey: measure p_e at each tau
            p_e_vals = np.zeros(len(tau_list))
            for i, tau in enumerate(tau_list):
                ctrl = create_ramsey_pulse(
                    self.t_rabi, tau,
                    omega_d=omega_d,
                    phase1=0.0, phase2=0.0,
                    qubit=self.qubit,
                )
                ctrl.t_list = ctrl.t_list - self.t_rabi[-1]

                H = (
                    QobjEvo(
                        self.qubit.H_list,
                        tlist=self.qubit.mag_signal.t_list,
                        order=1,
                    )
                    + QobjEvo(
                        ctrl.hamiltonian,
                        tlist=ctrl.t_list,
                        order=1,
                    )
                )
                result = mesolve(
                    H, self.qubit.state, t_global, [],
                    e_ops=[psi_e * psi_e.dag()],
                )
                p_e_vals[i] = result.expect[0][-1]

            # Fit Ramsey fringes: p_e(tau) = 0.5*(1 - cos(Delta*tau))
            # → frequency = omega_d + Delta
            p_e_arr = np.asarray(p_e_vals, dtype=float)
            # Simple estimate: peak-to-peak gives contrast
            C = float(np.ptp(p_e_arr))
            if C < 0.01:
                # No oscillation — qubit is at omega_d
                fitted_freq = float(omega_d)
            else:
                # Fit via FFT peak
                p_centered = p_e_arr - np.mean(p_e_arr)
                from numpy.fft import rfft, rfftfreq
                n_fft = 2048
                dt_val = float(tau_list[1] - tau_list[0])
                fft_vals = rfft(p_centered, n=n_fft)
                freqs = rfftfreq(n_fft, d=dt_val)
                peak_idx = int(np.argmax(np.abs(fft_vals)))
                detuning = float(2 * np.pi * freqs[peak_idx])
                fitted_freq = float(omega_d) + detuning

            frequency_list.append(fitted_freq)

        return CalibrationTable(
            name=f"flux_response_{self.method}",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f_phi",
            inputs=np.asarray(self.h_list, dtype=float),
            outputs=np.asarray(frequency_list, dtype=float),
            fit_params={
                "method": "ramsey",
                "omega_d": float(omega_d),
            },
            metadata={},
        )

    # ------------------------------------------------------------------
    # Cryoscope-based flux response calibration
    # ------------------------------------------------------------------

    def _calibrate_cryoscope(self) -> CalibrationTable:
        """Scan square-pulse height h, build φ(h) lookup table.

        Direct port of src/protocal.py:Calibration.calibrate case 3.

        Algorithm (per _cryoscope_implementation.md §3.3):
        1. For each h in h_list, apply a square-wave flux pulse:
           - h during free evolution (t_rabi[-1] to t_rabi[-1] + tau)
           - 0 during the π/2 pulses
        2. IQ readout via two Ramsey sequences (π/2 on X and Y axes)
        3. Extract φ = arctan2(p_e_Q - 0.5, p_e_I - 0.5)
        4. Unwrap phases → CalibrationTable(kind="phi_h")

        Physics: φ(h) = 2π · Δf_Q(h) · tau  (Gao 2021 §V)
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
            # Build square-wave flux: h during free evolution, 0 elsewhere
            t_sig = CONFIG.pulse.make_time(0, t_total)
            signal = np.zeros_like(t_sig, dtype=float)
            mask = (t_sig >= t_pi2_end) & (t_sig <= t_pi2_end + self.tau)
            signal[mask] = float(h)
            Phi = FluxSignal(type=8, t_list=t_sig, signal=signal)

            self.qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

            result = readout.measure(self.qubit)
            p_e_I_list.append(result["p_e_I"])
            p_e_Q_list.append(result["p_e_Q"])

        # Extract phase and unwrap
        p_e_I = np.asarray(p_e_I_list, dtype=float)
        p_e_Q = np.asarray(p_e_Q_list, dtype=float)
        varphi = np.arctan2(p_e_Q - 0.5, p_e_I - 0.5)
        varphi = np.unwrap(varphi, period=np.pi)

        return CalibrationTable(
            name=f"flux_response_{self.method}",
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

    # ------------------------------------------------------------------
    # Transient-based flux response calibration (Track B 1.2)
    # ------------------------------------------------------------------

    def _calibrate_transient(self) -> CalibrationTable:
        """Scan with known transient signal, fit Δω(Φ) polynomial.

        Replaces Track B's case 8 implementation.

        **Requires Track B 1.2** (transient frequency calibration).

        Raises
        ------
        NotImplementedError
            Until Track B 1.2 is complete.
        """
        raise NotImplementedError(
            "FluxResponseCalibration._calibrate_transient: "
            "requires Track B 1.2 (transient frequency calibration "
            "case 8, see _TODO_master.md 1.2)."
        )
