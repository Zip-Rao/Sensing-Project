"""sqc.calibration.frequency — Qubit frequency calibration.

Two major categories:
  - FluxResponseCalibration: f(Φ) curve via flux scan
  - SinglePointFrequencyCalibration: f01 at a single flux working point

See _sensing theory.md §闭环反馈控制 for closed-loop algorithm reference
(Vepsalainen 2022).
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


# ---------------------------------------------------------------------------
# Shared helper: Ramsey FFT frequency fitting
# ---------------------------------------------------------------------------

def _fit_ramsey_frequency(
    qubit: object,
    omega_d: float,
    tau_list: np.ndarray,
    t_rabi: np.ndarray,
    t_global: np.ndarray,
    flux: float = 0.0,
) -> float:
    """Run Ramsey at a given DC flux, return fitted f01 via FFT peak.

    Parameters
    ----------
    qubit : TransmonQubit
    omega_d : float
        Drive frequency (angular, GHz·2π).
    tau_list : np.ndarray
        Free evolution times (ns).
    t_rabi : np.ndarray
        Rabi pulse time axis (ns).
    t_global : np.ndarray
        Global evolution time axis (ns).
    flux : float
        DC flux offset (Φ₀). Default 0.

    Returns
    -------
    float
        Fitted qubit frequency (angular, GHz·2π).
    """
    n_levels = qubit.n_levels

    t_sig = CONFIG.pulse.make_time(0, 300)
    Phi = FluxSignal(
        type=1 if flux != 0.0 else 0,
        t_list=t_sig,
        amplitude=float(flux),
        offset=0.0,
    )
    qubit.qubit_in_mag(Phi, frame=1, omega_d=omega_d)

    psi_e = basis(n_levels, 1)
    p_e_vals = np.zeros(len(tau_list), dtype=float)

    for i, tau in enumerate(tau_list):
        ctrl = create_ramsey_pulse(
            t_rabi, tau,
            omega_d=omega_d,
            phase1=0.0, phase2=0.0,
            qubit=qubit,
        )
        ctrl.t_list = ctrl.t_list - t_rabi[-1]

        H = (
            QobjEvo(qubit.H_list, tlist=qubit.mag_signal.t_list, order=1)
            + QobjEvo(ctrl.hamiltonian, tlist=ctrl.t_list, order=1)
        )
        result = mesolve(
            H, qubit.state, t_global, [],
            e_ops=[psi_e * psi_e.dag()],
        )
        p_e_vals[i] = float(result.expect[0][-1])

    p_centered = p_e_vals - np.mean(p_e_vals)
    C = float(np.ptp(p_e_vals))

    if C < 0.01:
        return float(omega_d)

    from numpy.fft import rfft, rfftfreq
    n_fft = 2048
    dt_val = float(tau_list[1] - tau_list[0])
    fft_vals = rfft(p_centered, n=n_fft)
    freqs = rfftfreq(n_fft, d=dt_val)
    peak_idx = int(np.argmax(np.abs(fft_vals)))
    detuning = float(2 * np.pi * freqs[peak_idx])
    return float(omega_d) + detuning


# ===================================================================
# FluxResponseCalibration — f(Φ) curve
# ===================================================================

@dataclass
class FluxResponseCalibration(Calibration):
    """Magnetic flux response calibration — measure f(Φ) lookup table.

    Replaces src/protocal.py:Calibration.calibrate case 1 (ramsey scan).

    Methods
    -------
    - ``"ramsey"``: scan DC flux steps, measure frequency via Ramsey at each
      point → f(Φ) CalibrationTable (kind="f_phi").
    - ``"transient"``: unknown transient signal → polynomial fit of Δω(Φ).
      **Requires Track B 1.2.**

    Parameters
    ----------
    qubit : TransmonQubit
    method : str
        "ramsey" (complete) or "transient" (Track B 1.2 stub).
    h_list : np.ndarray or None
        Flux values to scan (Φ₀). Default linspace(-0.03, 0.03, 51).
    tau : float
        Free precession time per flux point (ns). Default 100.0.
    t_rabi : np.ndarray
        Rabi pulse time axis (ns).
    """

    qubit: object
    method: Literal["ramsey", "transient"] = "ramsey"
    h_list: np.ndarray | None = None
    tau: float = 100.0
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )

    def __post_init__(self):
        if self.h_list is None:
            self.h_list = np.linspace(-0.03, 0.03, 51)

    def calibrate(self) -> CalibrationTable:
        match self.method:
            case "ramsey":
                return self._calibrate_ramsey()
            case "transient":
                return self._calibrate_transient()

    # ------------------------------------------------------------------
    def _calibrate_ramsey(self) -> CalibrationTable:
        """Scan flux Φ, fit Ramsey detuning at each step → f(Φ) table."""
        omega_d = self.qubit.frequency
        tau_list = CONFIG.pulse.make_time(0, 200)
        t_global = CONFIG.pulse.t_global.copy()
        frequency_list: list[float] = []

        for h in self.h_list:
            f01 = _fit_ramsey_frequency(
                self.qubit, omega_d, tau_list,
                self.t_rabi, t_global, flux=float(h),
            )
            frequency_list.append(f01)

        return CalibrationTable(
            name="flux_response_ramsey",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f_phi",
            inputs=np.asarray(self.h_list, dtype=float),
            outputs=np.asarray(frequency_list, dtype=float),
            fit_params={"method": "ramsey", "omega_d": float(omega_d)},
            metadata={},
        )

    # ------------------------------------------------------------------
    def _calibrate_transient(self) -> CalibrationTable:
        raise NotImplementedError(
            "FluxResponseCalibration(method='transient'): "
            "requires Track B 1.2 (transient frequency calibration)."
        )


# ===================================================================
# SinglePointFrequencyCalibration — single-point f01
# ===================================================================

@dataclass
class SinglePointFrequencyCalibration(Calibration):
    """Single-point qubit frequency calibration.

    Measures f01 at a single flux working point. Replaces
    src/protocal.py:Calibration.calibrate case 0 (ramsey).

    Methods
    -------
    - ``"ramsey"``: single Ramsey FFT → f01.
    - ``"closed_loop"``: secant-method iterative convergence to a target
      frequency (Vepsalainen 2022 closed-loop feedback). The per-iteration
      frequency measurement is controlled by ``measure_method``.
    - ``"transient"``: single-shot transient measurement → f01.
      **Requires Track B 1.2.**

    Parameters
    ----------
    qubit : TransmonQubit
    method : str
        "ramsey", "closed_loop", or "transient".
    tau_list : np.ndarray or None
        Free evolution times for Ramsey (ns). Default make_time(0, 200).
    t_rabi : np.ndarray
        Rabi pulse time axis (ns).
    t_global : np.ndarray or None
        Global evolution time axis (ns).

    Closed-loop parameters
    ----------------------
    f_target : float or None
        Target qubit frequency (angular, GHz·2π). Required for closed_loop.
    epsilon_f : float
        Convergence tolerance (GHz·2π). Default 1e-4.
    V_a : float or None
        Lower voltage bound bracketing target. Required for closed_loop.
    V_b : float or None
        Upper voltage bound bracketing target. Required for closed_loop.
    max_iter : int
        Maximum secant iterations. Default 20.
    measure_method : str
        Per-iteration frequency measurement technique:
        - ``"ramsey"``: Ramsey FFT (default).
        - ``"transient"``: transient measurement (Track B 1.2).
    """

    qubit: object
    method: Literal["ramsey", "closed_loop", "transient"] = "ramsey"

    # Ramsey params
    tau_list: np.ndarray | None = None
    t_rabi: np.ndarray = field(
        default_factory=lambda: CONFIG.pulse.t_rabi.copy()
    )
    t_global: np.ndarray | None = None

    # Closed-loop params
    f_target: float | None = None
    epsilon_f: float = 1e-4
    V_a: float | None = None
    V_b: float | None = None
    max_iter: int = 20
    measure_method: Literal["ramsey", "transient"] = "ramsey"

    def __post_init__(self):
        if self.tau_list is None:
            self.tau_list = CONFIG.pulse.make_time(0, 200)
        if self.t_global is None:
            self.t_global = CONFIG.pulse.t_global.copy()

    # ------------------------------------------------------------------
    def calibrate(self) -> CalibrationTable:
        match self.method:
            case "ramsey":
                return self._calibrate_ramsey()
            case "closed_loop":
                return self._calibrate_closed_loop()
            case "transient":
                return self._calibrate_transient()

    # ------------------------------------------------------------------
    # Single-shot Ramsey
    # ------------------------------------------------------------------

    def _calibrate_ramsey(self) -> CalibrationTable:
        """Single Ramsey FFT → f01 at zero flux."""
        omega_d = self.qubit.frequency
        fitted_freq = _fit_ramsey_frequency(
            self.qubit, omega_d, self.tau_list,
            self.t_rabi, self.t_global, flux=0.0,
        )
        return CalibrationTable(
            name="frequency_ramsey",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f01",
            inputs=np.array([0.0]),
            outputs=np.array([fitted_freq]),
            fit_params={
                "method": "ramsey",
                "omega_d": float(omega_d),
            },
            metadata={},
        )

    # ------------------------------------------------------------------
    # Closed-loop feedback (Vepsalainen 2022)
    # ------------------------------------------------------------------

    def _calibrate_closed_loop(self) -> CalibrationTable:
        """Secant-method closed-loop frequency calibration.

        Iteratively adjusts flux bias voltage to drive qubit frequency to
        f_target.  Uses the secant method for root-finding on
        r(V) = f_Q(V) - f_target = 0.

        Requires f_target, V_a, V_b to be set.  V_a and V_b must bracket
        the target: r(V_a) · r(V_b) < 0.
        """
        if self.f_target is None:
            raise ValueError("f_target is required for closed_loop method.")
        if self.V_a is None or self.V_b is None:
            raise ValueError(
                "V_a and V_b are required for closed_loop method. "
                "Run FluxResponseCalibration first to bracket the target."
            )

        f_target = self.f_target
        V_a, V_b = self.V_a, self.V_b
        epsilon = self.epsilon_f

        # Initial midpoint
        V_n = (V_a + V_b) / 2.0
        r_n = self._measure_frequency(V_n) - f_target
        n_iter = 0

        # Store previous step for secant
        V_prev = V_a
        r_prev = self._measure_frequency(V_a) - f_target

        history: list[dict] = []

        while abs(r_n) > epsilon and n_iter < self.max_iter:
            # Secant update: V_next = V_n - r_n * (V_n - V_prev) / (r_n - r_prev)
            denom = r_n - r_prev
            if abs(denom) < 1e-15:
                # Secant would divide by zero — fall back to bisection
                V_next = (V_a + V_b) / 2.0
            else:
                V_next = V_n - r_n * (V_n - V_prev) / denom

            # Bound check: fall back to midpoint if out of bounds
            if V_next < min(V_a, V_b) or V_next > max(V_a, V_b):
                V_n = (V_a + V_b) / 2.0
            else:
                V_prev, r_prev = V_n, r_n
                V_n = V_next

            f_n = self._measure_frequency(V_n)
            r_n = f_n - f_target
            n_iter += 1

            history.append({
                "iter": n_iter, "V": float(V_n),
                "f": float(f_n), "residual": float(r_n),
            })

            # Tighten bracket
            if r_n > 0:
                V_b = V_n
            else:
                V_a = V_n

        return CalibrationTable(
            name="frequency_closed_loop",
            qubit_name=getattr(self.qubit, "name", "qubit"),
            kind="f01",
            inputs=np.array([V_n]),
            outputs=np.array([f_target + r_n]),
            fit_params={
                "method": "closed_loop",
                "measure_method": self.measure_method,
                "f_target": float(f_target),
                "V_opt": float(V_n),
                "n_iter": n_iter,
                "residual": float(r_n),
                "converged": abs(r_n) <= epsilon,
                "history": history,
            },
            metadata={},
        )

    def _measure_frequency(self, flux: float) -> float:
        """Measure f01 at a given flux offset.

        Dispatches to Ramsey or transient measurement based on
        self.measure_method.
        """
        if self.measure_method == "ramsey":
            omega_d = self.qubit.frequency
            return _fit_ramsey_frequency(
                self.qubit, omega_d, self.tau_list,
                self.t_rabi, self.t_global, flux=flux,
            )
        elif self.measure_method == "transient":
            return self._measure_frequency_transient(flux)
        else:
            raise ValueError(f"Unknown measure_method: {self.measure_method}")

    def _measure_frequency_transient(self, flux: float) -> float:
        """Transient-based single-point frequency measurement.

        Requires Track B 1.2.
        """
        raise NotImplementedError(
            "Transient frequency measurement requires Track B 1.2."
        )

    # ------------------------------------------------------------------
    # Transient (stub)
    # ------------------------------------------------------------------

    def _calibrate_transient(self) -> CalibrationTable:
        raise NotImplementedError(
            "SinglePointFrequencyCalibration(method='transient'): "
            "requires Track B 1.2."
        )
