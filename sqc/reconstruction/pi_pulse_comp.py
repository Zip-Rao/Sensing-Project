"""sqc.reconstruction.pi_pulse_comp — π-pulse compensation reconstruction.

    Φ_tail(τ) = -z*(τ)

Note: due to the finite pi-pulse width T, z*(τ) measures the
moving average of the tail over [τ, τ+T], not the instantaneous
value.  For exponential tails this gives a constant attenuation
factor; the decay time-constant is preserved.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sqc.control.flux_signal import FluxSignal
from sqc.reconstruction.base import Reconstruction


@dataclass
class PiPulseCompReconstruction(Reconstruction):
    """π-pulse compensation waveform reconstruction.

    From the 2D scan P_e(τ, z), the optimal compensation height z*(τ)
    directly negates the tail:  Φ_tail(τ) = -z*(τ).
    """

    def reconstruct(self, measurement, **kwargs) -> FluxSignal:
        """Reconstruct tail flux from pi-pulse compensation data.

        Parameters
        ----------
        measurement : ExperimentResult
            Must contain data["z_star"] and axes["tau"].

        Returns
        -------
        FluxSignal
            Reconstructed tail flux waveform (type=8).
        """
        z_star = np.asarray(measurement.data["z_star"], dtype=float)
        t_axis = np.asarray(measurement.axes["tau"], dtype=float)
        return FluxSignal(type=8, t_list=t_axis, signal=-z_star)
