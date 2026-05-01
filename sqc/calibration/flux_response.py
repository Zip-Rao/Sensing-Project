"""Flux-response calibration routines."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from sqc.calibration.base import Calibration, CalibrationTable
from sqc.control.flux_signal import FluxSignal
from sqc.hardware.readout import IQReadoutModel


@dataclass
class FluxResponseCalibration(Calibration):
    """Calibrate flux-to-phase/frequency response."""

    qubit: object
    method: Literal["ramsey", "cryoscope", "transient"] = "cryoscope"
    h_list: np.ndarray | None = None
    tau: float = 100.0

    def __post_init__(self) -> None:
        if self.h_list is None:
            self.h_list = np.linspace(-0.03, 0.03, 21)

    def calibrate(self) -> CalibrationTable:
        """Run the selected calibration method."""
        if self.method == "cryoscope":
            return self._calibrate_cryoscope()
        if self.method == "ramsey":
            return self._calibrate_ramsey()
        if self.method == "transient":
            raise NotImplementedError(
                "Transient flux-response calibration awaits Track B case 8"
            )
        raise ValueError(f"unknown calibration method: {self.method}")

    def _calibrate_cryoscope(self) -> CalibrationTable:
        """Scan square-pulse height and measure phi(h)."""
        readout = IQReadoutModel(tau=self.tau)
        varphi_list = []

        def make_signal(h):
            t_list = np.linspace(0, self.tau + 20, 240)
            samples = np.zeros_like(t_list)
            samples[(t_list >= 10) & (t_list <= self.tau + 10)] = h
            return FluxSignal(type=8, t_list=t_list, signal=samples)

        for h in self.h_list:
            phi = make_signal(h)
            self.qubit.qubit_in_mag(phi, frame=1, omega_d=self.qubit.frequency)
            measured = readout.measure(self.qubit)
            varphi_list.append(
                np.arctan2(measured["p_e_Q"] - 0.5, measured["p_e_I"] - 0.5)
            )

        varphi = np.unwrap(np.asarray(varphi_list), period=np.pi)
        return CalibrationTable(
            qubit_name=self.qubit.spec().name,
            kind="phi_h",
            inputs=np.asarray(self.h_list),
            outputs=varphi,
            fit_params={"tau": self.tau, "method": "cryoscope"},
            metadata={},
        )

    def _calibrate_ramsey(self) -> CalibrationTable:
        """Use the analytic transmon model as a static Ramsey calibration table."""
        frequencies = np.asarray([self.qubit.spec().frequency(h) for h in self.h_list])
        return CalibrationTable(
            qubit_name=self.qubit.spec().name,
            kind="frequency_flux",
            inputs=np.asarray(self.h_list),
            outputs=frequencies,
            fit_params={"method": "ramsey_static_model"},
            metadata={},
        )
