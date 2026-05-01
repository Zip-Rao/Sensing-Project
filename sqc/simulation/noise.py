"""Noise helpers."""
from __future__ import annotations

import numpy as np


def generate_1f_noise(t_list, amplitude, f_min, f_max, seed: int | None = None):
    """Generate a real-valued 1/f noise trace."""
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
    spectrum[n // 2 + 1:] = np.conj(spectrum[1:n // 2][::-1])
    return np.fft.ifft(spectrum).real
