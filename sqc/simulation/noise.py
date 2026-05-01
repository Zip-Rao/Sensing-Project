"""sqc.simulation.noise — 1/f noise generator.

Verbatim port from TransmonQubit.generate_1f_noise().
"""
from __future__ import annotations

import numpy as np


def generate_1f_noise(
    t_list: np.ndarray,
    amplitude: float,
    f_min: float,
    f_max: float,
    seed: int | None = None,
) -> np.ndarray:
    """Generate 1/f noise time series via FFT.

    Parameters
    ----------
    t_list : np.ndarray
        Time array (ns).
    amplitude : float
        Noise amplitude.
    f_min : float
        Minimum frequency (Hz).
    f_max : float
        Maximum frequency (Hz).
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    np.ndarray
        Real-valued noise time series.
    """
    if seed is not None:
        np.random.seed(seed)
    dt = t_list[1] - t_list[0]
    n = len(t_list)
    freqs = np.fft.fftfreq(n, dt)
    spectrum = np.zeros(n, dtype=complex)
    for i in range(1, n // 2):
        f = abs(freqs[i])
        if f_min <= f <= f_max:
            spectrum[i] = amplitude / np.sqrt(f) * (
                np.random.normal() + 1j * np.random.normal()
            )
    # Enforce Hermitian symmetry for real-valued IFFT
    spectrum[n // 2 + 1:] = np.conj(spectrum[1 : n // 2][::-1])
    return np.fft.ifft(spectrum).real
