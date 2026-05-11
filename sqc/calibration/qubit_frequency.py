"""sqc.calibration.qubit_frequency — Qubit frequency calibration.

Calibrate qubit f_01 via Ramsey free precession (case 0) and
transient-based Δω(Φ) polynomial fitting (case 8).

Currently only QubitFrequencyCalibration (case 0) is implemented.
TransientFrequencyCalibration requires Track B 1.2.

See _TODO_master.md 1.2.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qutip import QobjEvo, basis, mesolve

from sqc.config import CONFIG
from sqc.calibration.base import Calibration, CalibrationTable
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse


@dataclass
class QubitFrequencyCalibration(Calibration):
    """Calibrate f_01 via Ramsey free precession.

    Replaces src/protocal.py:Calibration.calibrate case 0.

    The algorithm:
    1. Prepare |0⟩, apply π/2 pulse
    2. Wait variable time tau, apply second π/2 pulse
    3. Measure p_e(tau), fit Ramsey fringes
    4. Extract f_01 from fitted oscillation frequency

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit to calibrate.
    tau_list : np.ndarray or None
        Free evolution times (ns). Default linspace(0, 200, 100).
    t_rabi : np.ndarray or None
        Rabi pulse time axis (ns). Default linspace(0, 10, 20).
    t_global : np.ndarray or None
        Global evolution time (ns). Default linspace(-50, 400, 900).
    """

    qubit: object  # TransmonQubit
    tau_list: np.ndarray | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    t_global: np.ndarray | None = None

    def __post_init__(self):
        if self.tau_list is None:
            self.tau_list = CONFIG.pulse.make_time(0, 200)
        if self.t_global is None:
            self.t_global = CONFIG.pulse.t_global.copy()

    def calibrate(self) -> CalibrationTable:
        """Run Ramsey frequency calibration.

        Returns
        -------
        CalibrationTable
            With kind="f01", outputs = [fitted_frequency].
        """
        omega_d = self.qubit.frequency
        n_levels = self.qubit.n_levels

        # Create a zero-flux signal (qubit in idle)
        t_sig = CONFIG.pulse.make_time(0, 300)
        Phi = FluxSignal(type=0, t_list=t_sig)
        self.qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

        psi_e = basis(n_levels, 1)
        p_e_vals = np.zeros(len(self.tau_list), dtype=float)

        for i, tau in enumerate(self.tau_list):
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
                H, self.qubit.state, self.t_global, [],
                e_ops=[psi_e * psi_e.dag()],
            )
            p_e_vals[i] = float(result.expect[0][-1])

        # Fit Ramsey fringes: p_e(tau) = 0.5*(1 - C*cos(Delta*tau + phi0))
        p_centered = p_e_vals - np.mean(p_e_vals)
        C = float(np.ptp(p_e_vals))

        if C < 0.01:
            # No clear oscillation — use nominal frequency
            fitted_freq = float(omega_d)
        else:
            # Estimate detuning via FFT
            from numpy.fft import rfft, rfftfreq
            n_fft = 2048
            dt_val = float(self.tau_list[1] - self.tau_list[0])
            fft_vals = rfft(p_centered, n=n_fft)
            freqs = rfftfreq(n_fft, d=dt_val)
            peak_idx = int(np.argmax(np.abs(fft_vals)))
            detuning = float(2 * np.pi * freqs[peak_idx])
            fitted_freq = float(omega_d) + detuning

        return CalibrationTable(
            name="qubit_frequency",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f01",
            inputs=np.array([0.0]),  # flux = 0
            outputs=np.array([fitted_freq]),
            fit_params={
                "method": "ramsey",
                "omega_d": float(omega_d),
                "contrast": float(C),
            },
            metadata={},
        )


@dataclass
class TransientFrequencyCalibration(Calibration):
    """Calibrate Δω(Φ) polynomial coefficients via transient sensing.

    Replaces Track B's case 8 implementation.

    **Requires Track B 1.2** (transient frequency calibration case 8).

    Parameters
    ----------
    qubit : TransmonQubit
        Qubit to calibrate.
    polynomial_order : int
        Order of polynomial fit for Δω(Φ). Default 3.
    test_signal_kind : str
        Kind of test signal: "ramp", "sine", or "gaussian". Default "ramp".
    """

    qubit: object  # TransmonQubit
    polynomial_order: int = 3
    test_signal_kind: str = "ramp"

    def calibrate(self) -> CalibrationTable:
        """Run transient frequency calibration.

        Raises
        ------
        NotImplementedError
            Until Track B 1.2 is complete.
        """
        raise NotImplementedError(
            "TransientFrequencyCalibration.calibrate: "
            "requires Track B 1.2 (transient frequency calibration "
            "case 8, see _TODO_master.md 1.2)."
        )
