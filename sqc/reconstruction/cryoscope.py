"""sqc.reconstruction.cryoscope — CryoscopeReconstruction.

Cryoscope waveform reconstruction from phase-vs-truncation data.

Full implementation requires Track B 1.1 (Cryoscope case 6/7).
Currently only the class skeleton is available.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from sqc.calibration.base import CalibrationTable
from sqc.reconstruction.base import Reconstruction


@dataclass
class CryoscopeReconstruction(Reconstruction):
    """Cryoscope waveform reconstruction.

    Replaces Analysis.get_signal_from_cryoscope + get_h_from_phi.

    Full implementation requires Track B 1.1 (Cryoscope case 6/7).
    See _TODO_master.md 1.1.

    Parameters
    ----------
    calibration : CalibrationTable
        Phase-to-flux calibration table (phi(h) with inverse).
    tau : float
        Calibration square-pulse length (ns).
    method : str
        Reconstruction method:
        - "calib_inverse": use calibration inverse directly
        - "SG_diff": Savitzky-Golay derivative
        - "diff": simple gradient
    """

    calibration: CalibrationTable
    tau: float
    method: Literal["calib_inverse", "SG_diff", "diff"] = "calib_inverse"

    def reconstruct(
        self,
        measurement,
        kernel=None,
        calibration=None,
        dt: float | None = None,
    ) -> np.ndarray:
        """Reconstruct flux signal h(t) from cryoscope phase data.

        **Full implementation requires Track B 1.1** (Cryoscope case 6/7
        calibration). The case 5 data path measures φ(t_d) at truncation
        delays; case 6/7 adds the full reconstruction pipeline including
        calibration-curve inversion and optional SG/diff differentiation.

        Currently raises NotImplementedError.

        Parameters
        ----------
        measurement : ExperimentResult
            Must contain data["varphi"] and axes["trunc"].
        kernel : optional
            Not used in cryoscope reconstruction.
        calibration : CalibrationTable, optional
            Override calibration table.
        dt : float, optional
            Time step. Computed from trunc axis if not given.

        Returns
        -------
        np.ndarray
            Reconstructed flux signal h(t).

        Raises
        ------
        NotImplementedError
            Always until Track B 1.1 is complete.
        """
        raise NotImplementedError(
            "CryoscopeReconstruction.reconstruct: requires Track B 1.1 "
            "(Cryoscope case 6/7 calibration, see _TODO_master.md 1.1). "
            "The measurement-only CryoscopeExperiment (case 5) is "
            "available in sqc.experiments.cryoscope."
        )
