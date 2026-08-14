"""Shared numerical helpers for the P9 short-time distortion benchmark."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.signal import lfilter, savgol_filter


DT_NS = 0.5
N_SAMPLES = 512
N_SHOTS = 4096
N_SEEDS = 20
SEEDS = np.arange(202608050, 202608050 + N_SEEDS, dtype=int)
FIR_TAPS = 64
FIR_DELAY = 24
FIR_RIDGE = 3e-3
TRANSIENT_LAMBDA = 10.0
CRYO_SG_WINDOW = 9
CRYO_PHASE_GAIN = 6.0
READOUT_CONTRAST = 0.90
TRAIN_AMPLITUDES = (0.003, 0.006, 0.009)


@dataclass(frozen=True)
class DistortionCase:
    name: str
    family: str
    parameter: float
    impulse: np.ndarray


def _normalise_dc(h: np.ndarray) -> np.ndarray:
    h = np.asarray(h, dtype=float)
    total = float(np.sum(h))
    if abs(total) < 1e-12:
        raise ValueError("distortion impulse has zero DC gain")
    return h / total


def _cascade(*stages: np.ndarray, n: int = N_SAMPLES) -> np.ndarray:
    h = np.array([1.0])
    for stage in stages:
        h = np.convolve(h, stage)
    return _normalise_dc(h[:n])


def identity_impulse(n: int = N_SAMPLES) -> np.ndarray:
    h = np.zeros(n)
    h[0] = 1.0
    return h


def lowpass_impulse(tau_ns: float, n: int = N_SAMPLES) -> np.ndarray:
    alpha = np.exp(-DT_NS / tau_ns)
    return (1.0 - alpha) * alpha ** np.arange(n)


def echo_impulse(delay_ns: float, amplitude: float = 0.04,
                 n: int = N_SAMPLES) -> np.ndarray:
    h = np.zeros(n)
    delay = int(round(delay_ns / DT_NS))
    h[0] = 1.0 - amplitude
    h[delay] = amplitude
    return h


def slow_tail_impulse(amplitude: float = 0.01, tau_ns: float = 200.0,
                      n: int = N_SAMPLES) -> np.ndarray:
    t = np.arange(n) * DT_NS
    step = 1.0 - amplitude * np.exp(-t / tau_ns)
    h = np.empty(n)
    h[0] = step[0]
    h[1:] = np.diff(step)
    return _normalise_dc(h)


def ring_impulse(freq_mhz: float, amplitude: float = 0.025,
                 decay_ns: float = 15.0, n: int = N_SAMPLES) -> np.ndarray:
    t = np.arange(n) * DT_NS
    perturb = amplitude * np.exp(-t / decay_ns) * np.sin(
        2.0 * np.pi * freq_mhz * 1e-3 * t
    )
    perturb[0] -= np.sum(perturb)
    h = identity_impulse(n) + perturb
    return _normalise_dc(h)


def make_cases() -> list[DistortionCase]:
    cases = [DistortionCase("D0_identity", "D0", 0.0, identity_impulse())]
    for tau in (2.0, 5.0, 10.0, 20.0):
        cases.append(DistortionCase(
            f"D1_lowpass_{tau:g}ns", "D1", tau, lowpass_impulse(tau)
        ))
    for delay in (3.0, 6.0, 12.0, 20.0):
        cases.append(DistortionCase(
            f"D2_echo_{delay:g}ns", "D2", delay, echo_impulse(delay)
        ))
    for freq in (50.0, 100.0, 200.0, 300.0):
        cases.append(DistortionCase(
            f"D3_ring_{freq:g}MHz", "D3", freq, ring_impulse(freq)
        ))
    main = _cascade(
        lowpass_impulse(10.0), echo_impulse(6.0), slow_tail_impulse()
    )
    cases.append(DistortionCase("D4_main_mixed", "D4", 6.0, main))
    return cases


def apply_lti(x: np.ndarray, h: np.ndarray) -> np.ndarray:
    return lfilter(h, [1.0], np.asarray(x, dtype=float))


def multisine(seed: int, amplitude: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = np.arange(N_SAMPLES)
    # Exact DFT-bin tones from 7.8 to 351.6 MHz avoid spectral leakage.
    bins = np.arange(2, 91)
    phases = rng.uniform(0.0, 2.0 * np.pi, len(bins))
    x = np.sum(np.cos(2.0 * np.pi * bins[:, None] * n / N_SAMPLES
                      + phases[:, None]), axis=0)
    x -= np.mean(x)
    x *= amplitude / np.max(np.abs(x))
    return x


def gain_probe() -> np.ndarray:
    t = np.arange(N_SAMPLES) * DT_NS
    x = np.exp(-0.5 * ((t - 85.0) / 10.0) ** 2)
    x -= np.exp(-0.5 * ((t - 125.0) / 10.0) ** 2)
    return 0.004 * x / np.max(np.abs(x))


def heldout_waveforms() -> dict[str, np.ndarray]:
    t = np.arange(N_SAMPLES) * DT_NS
    dual = 0.12 * (
        np.exp(-0.5 * ((t - 65.0) / 3.0) ** 2)
        - 0.75 * np.exp(-0.5 * ((t - 88.0) / 5.0) ** 2)
    )
    flat = np.zeros_like(t)
    rise = (t >= 55.0) & (t < 65.0)
    plateau = (t >= 65.0) & (t < 105.0)
    fall = (t >= 105.0) & (t < 115.0)
    flat[rise] = 0.12 * np.sin(0.5 * np.pi * (t[rise] - 55.0) / 10.0) ** 2
    flat[plateau] = 0.12
    flat[fall] = 0.12 * np.cos(0.5 * np.pi * (t[fall] - 105.0) / 10.0) ** 2
    return {"dual_pulse": dual, "short_flat_top": flat}


def load_transient_kernel(repo: Path) -> np.ndarray:
    cache = repo / (
        "result_sqc/predistortion/P8_transient_report_figures/"
        "fig_E1_transient_characterization_data.npz"
    )
    with np.load(cache) as data:
        return np.asarray(data["kernel"], dtype=float)


def _binomial_probability(rng: np.random.Generator, p: np.ndarray,
                          shots: int) -> np.ndarray:
    return rng.binomial(shots, np.clip(p, 1e-6, 1.0 - 1e-6)) / shots


def transient_observe(y: np.ndarray, kernel: np.ndarray, seed: int,
                      shots: int = N_SHOTS) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    delta = np.convolve(kernel, y, mode="full") * DT_NS
    p_signal = 0.5 + 0.5 * delta
    p_zero = 0.5 - 0.5 * delta
    observed = _binomial_probability(rng, p_signal, shots) - _binomial_probability(
        rng, p_zero, shots
    )
    return observed, delta


def transient_reconstruct(observed: np.ndarray, kernel: np.ndarray,
                          gain: float = 1.0) -> np.ndarray:
    n = len(observed)
    padded = np.zeros(n)
    padded[:len(kernel)] = kernel
    yf = np.fft.fft(observed)
    kf = np.fft.fft(padded)
    xf = np.conj(kf) * yf / (np.abs(kf) ** 2 + TRANSIENT_LAMBDA ** 2) / DT_NS
    n_out = n - len(kernel) + 1
    return gain * np.fft.ifft(xf).real[:n_out]


def cryoscope_observe(y: np.ndarray, seed: int,
                      shots: int = N_SHOTS) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    phase = CRYO_PHASE_GAIN * np.cumsum(y) * DT_NS
    p_i = 0.5 - 0.5 * READOUT_CONTRAST * np.sin(phase)
    p_q = 0.5 + 0.5 * READOUT_CONTRAST * np.cos(phase)
    pi_obs = _binomial_probability(rng, p_i, shots)
    pq_obs = _binomial_probability(rng, p_q, shots)
    return pi_obs, pq_obs, phase


def cryoscope_reconstruct(p_i: np.ndarray, p_q: np.ndarray,
                          gain: float = 1.0) -> np.ndarray:
    phase = np.unwrap(np.arctan2(0.5 - p_i, p_q - 0.5))
    phase = savgol_filter(phase, CRYO_SG_WINDOW, 2)
    return gain * np.gradient(phase, DT_NS) / CRYO_PHASE_GAIN


def calibrate_protocol_gain(protocol: str, kernel: np.ndarray, seed: int,
                            shots: int = N_SHOTS) -> float:
    x = gain_probe()
    if protocol == "transient":
        obs, _ = transient_observe(x, kernel, seed, shots=shots)
        raw = transient_reconstruct(obs, kernel)
    elif protocol == "cryoscope":
        pi, pq, _ = cryoscope_observe(x, seed, shots=shots)
        raw = cryoscope_reconstruct(pi, pq)
    else:
        raise ValueError(protocol)
    denom = float(np.dot(raw, raw))
    return float(np.dot(raw, x) / denom) if denom > 1e-15 else 1.0


def estimate_h(xs: list[np.ndarray], ys: list[np.ndarray], ridge: float = 1e-3):
    xf = np.stack([np.fft.rfft(x) for x in xs])
    yf = np.stack([np.fft.rfft(y) for y in ys])
    power = np.sum(np.abs(xf) ** 2, axis=0)
    reg = ridge * np.max(power)
    h = np.sum(yf * np.conj(xf), axis=0) / (power + reg)
    h[0] = 1.0 + 0.0j
    support = power >= 0.05 * np.max(power)
    support[0] = True
    return h, support, power


def true_frequency_response(h: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    freq_ghz = np.fft.rfftfreq(N_SAMPLES, d=DT_NS)
    return freq_ghz, np.fft.rfft(np.pad(h[:N_SAMPLES], (0, 0)))


def qualified_mask(h_est: np.ndarray, h_true: np.ndarray,
                   support: np.ndarray) -> np.ndarray:
    amp_error = np.abs(np.abs(h_est) - np.abs(h_true)) / np.maximum(
        np.abs(h_true), 1e-9
    )
    phase_error = np.abs(np.angle(h_est * np.conj(h_true)))
    mask = support & (amp_error <= 0.05) & (phase_error <= 0.10)
    mask[0] = True
    return mask


def effective_bandwidth_mhz(mask: np.ndarray, support: np.ndarray) -> float:
    freq = np.fft.rfftfreq(N_SAMPLES, d=DT_NS) * 1000.0
    supported_idx = np.flatnonzero(support & (freq <= 500.0))
    if not len(supported_idx):
        return 0.0
    good = []
    for idx in supported_idx:
        if mask[idx]:
            good.append(idx)
        else:
            break
    return float(freq[good[-1]]) if good else 0.0


def design_inverse_fir(h_est: np.ndarray, qualified: np.ndarray,
                       support: np.ndarray) -> np.ndarray:
    freq = np.fft.rfftfreq(N_SAMPLES, d=DT_NS)
    omega = 2.0 * np.pi * freq * DT_NS
    use = qualified & support
    use[0] = True
    # Smoothly keep the contiguous qualified band; every protocol uses this rule.
    idx = np.flatnonzero(use)
    if len(idx) < 4:
        idx = np.arange(min(8, len(h_est)))
    k = np.arange(FIR_TAPS)
    a = h_est[idx, None] * np.exp(-1j * omega[idx, None] * k[None, :])
    target = np.exp(-1j * omega[idx] * FIR_DELAY)
    a_real = np.vstack([a.real, a.imag])
    b_real = np.concatenate([target.real, target.imag])
    lhs = a_real.T @ a_real + FIR_RIDGE * np.eye(FIR_TAPS)
    rhs = a_real.T @ b_real
    taps = np.linalg.solve(lhs, rhs)
    dc = float(np.sum(taps))
    if abs(dc) > 1e-12:
        taps /= dc
    return taps


def compensate(target: np.ndarray, h_true: np.ndarray,
               inverse_taps: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    padded = np.pad(target, (0, FIR_DELAY + FIR_TAPS))
    awg = lfilter(inverse_taps, [1.0], padded)
    chip = lfilter(h_true, [1.0], awg)
    aligned = chip[FIR_DELAY:FIR_DELAY + len(target)]
    return awg, aligned


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def edge_rmse(a: np.ndarray, b: np.ndarray) -> float:
    derivative = np.abs(np.gradient(b))
    threshold = 0.15 * max(float(np.max(derivative)), 1e-12)
    mask = derivative >= threshold
    mask = np.convolve(mask.astype(int), np.ones(17, dtype=int), mode="same") > 0
    return rmse(np.asarray(a)[mask], np.asarray(b)[mask])


def hardware_metrics(awg: np.ndarray, taps: np.ndarray) -> dict[str, float | bool]:
    peak = float(np.max(np.abs(awg)))
    slew = float(np.max(np.abs(np.diff(awg))) / DT_NS)
    noise_gain = float(np.max(np.abs(np.fft.rfft(taps, 2048))))
    return {
        "awg_peak": peak,
        "awg_slew_phi0_per_ns": slew,
        "filter_noise_gain": noise_gain,
        "constraints_pass": bool(peak <= 0.25 and slew <= 0.08 and noise_gain <= 3.0),
    }
