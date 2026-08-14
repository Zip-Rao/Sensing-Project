"""Core models for the frequency-deembedded pre-CZ benchmark.

The module keeps protocol response R(f), control-line response H(f), and the
inverse filter C(f) as separate objects. No test-line truth enters protocol
calibration or filter selection.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.signal import lfilter, savgol_filter


DT_NS = 0.5
N_SAMPLES = 512
# Lowest D0-only resource-grid point at which both protocols validate at least
# MIN_QUALIFIED_BINS continuous non-zero frequency bins.
MAIN_SHOTS = 262144
N_SEEDS = 20
SEEDS = np.arange(202608150, 202608150 + N_SEEDS, dtype=int)
TUNING_SEEDS = np.arange(202608120, 202608128, dtype=int)
TRAIN_AMPLITUDES = (0.003, 0.006, 0.009)
N_PROBES = len(TRAIN_AMPLITUDES)
READOUT_CONTRAST = 0.90
CRYO_PHASE_GAIN = 6.0
R_FLOOR_FRACTION = 0.10
MIN_QUALIFIED_BINS = 8
COMMON_WAVEFORM_BAND_MHZ = 50.0
FILTER_LENGTHS = (16, 32, 64, 72)
FILTER_RIDGES = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1)
FILTER_DELAY = 24


@dataclass(frozen=True)
class DistortionCase:
    name: str
    family: str
    parameter: float
    impulse: np.ndarray


@dataclass(frozen=True)
class ProtocolChoice:
    protocol: str
    parameter: float
    validation_bandwidth_mhz: float
    validation_error: float


@dataclass(frozen=True)
class FilterChoice:
    taps: np.ndarray
    n_taps: int
    ridge: float
    validation_rmse: float
    awg_peak: float
    awg_slew: float
    noise_gain: float


def _normalise_dc(h: np.ndarray) -> np.ndarray:
    h = np.asarray(h, dtype=float)
    total = float(np.sum(h))
    if abs(total) < 1e-12:
        raise ValueError("zero-DC distortion impulse")
    return h / total


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
    h[0] = 1.0 - amplitude
    h[int(round(delay_ns / DT_NS))] = amplitude
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
    p = amplitude * np.exp(-t / decay_ns) * np.sin(
        2.0 * np.pi * freq_mhz * 1e-3 * t
    )
    p[0] -= np.sum(p)
    return _normalise_dc(identity_impulse(n) + p)


def cascade(*stages: np.ndarray, n: int = N_SAMPLES) -> np.ndarray:
    h = np.array([1.0])
    for stage in stages:
        h = np.convolve(h, stage)
    return _normalise_dc(h[:n])


def make_cases() -> list[DistortionCase]:
    cases = [DistortionCase("D0_identity", "D0", 0.0, identity_impulse())]
    cases += [DistortionCase(f"D1_lowpass_{v:g}ns", "D1", v,
                             lowpass_impulse(v)) for v in (2., 5., 10., 20.)]
    cases += [DistortionCase(f"D2_echo_{v:g}ns", "D2", v,
                             echo_impulse(v)) for v in (3., 6., 12., 20.)]
    cases += [DistortionCase(f"D3_ring_{v:g}MHz", "D3", v,
                             ring_impulse(v)) for v in (50., 100., 200., 300.)]
    main = cascade(lowpass_impulse(10.0), echo_impulse(6.0),
                   slow_tail_impulse())
    cases.append(DistortionCase("D4_main_mixed", "D4", 6.0, main))
    return cases


def apply_lti(x: np.ndarray, h: np.ndarray) -> np.ndarray:
    return lfilter(np.asarray(h, float), [1.0], np.asarray(x, float))


def apply_lti_periodic(x: np.ndarray, h: np.ndarray, periods: int = 8) -> np.ndarray:
    """Return one steady-state period for a periodic system-ID probe."""
    tiled = np.tile(np.asarray(x, float), int(periods))
    out = lfilter(np.asarray(h, float), [1.0], tiled)
    return out[-len(x):]


def multisine(seed: int, amplitude: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = np.arange(N_SAMPLES)
    bins = np.arange(1, 91)  # 3.90625--351.5625 MHz, every DFT bin excited
    phases = rng.uniform(0.0, 2.0 * np.pi, len(bins))
    x = np.sum(np.cos(2.0 * np.pi * bins[:, None] * n / N_SAMPLES
                      + phases[:, None]), axis=0)
    x -= np.mean(x)
    return amplitude * x / np.max(np.abs(x))


def probe_set(seed: int, split: int) -> list[np.ndarray]:
    return [multisine(seed + 10000 * split + 997 * j, amp)
            for j, amp in enumerate(TRAIN_AMPLITUDES)]


def heldout_waveforms() -> dict[str, np.ndarray]:
    t = np.arange(N_SAMPLES) * DT_NS
    dual = 0.12 * (np.exp(-0.5 * ((t - 65.0) / 3.0) ** 2)
                   - 0.75 * np.exp(-0.5 * ((t - 88.0) / 5.0) ** 2))
    flat = np.zeros_like(t)
    rise = (t >= 55.0) & (t < 65.0)
    top = (t >= 65.0) & (t < 105.0)
    fall = (t >= 105.0) & (t < 115.0)
    flat[rise] = 0.12 * np.sin(0.5 * np.pi * (t[rise] - 55.0) / 10.0) ** 2
    flat[top] = 0.12
    flat[fall] = 0.12 * np.cos(0.5 * np.pi * (t[fall] - 105.0) / 10.0) ** 2
    return {"validation_dual_pulse": dual, "test_short_flat_top": flat}


def load_transient_kernel(repo: Path) -> np.ndarray:
    path = repo / ("result_sqc/predistortion/P8_transient_report_figures/"
                   "fig_E1_transient_characterization_data.npz")
    with np.load(path) as z:
        return np.asarray(z["kernel"], dtype=float)


def _sample(rng: np.random.Generator, p: np.ndarray, shots: int) -> np.ndarray:
    return rng.binomial(shots, np.clip(p, 1e-7, 1 - 1e-7)) / shots


def transient_observe(y: np.ndarray, kernel: np.ndarray, seed: int,
                      shots: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    delta = transient_forward_periodic(y, kernel)
    return (_sample(rng, 0.5 + 0.5 * delta, shots)
            - _sample(rng, 0.5 - 0.5 * delta, shots))


def transient_forward_periodic(y: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """One steady-state measurement period for a periodically repeated probe."""
    kp = np.zeros(len(y))
    kp[:len(kernel)] = kernel
    return np.fft.ifft(np.fft.fft(kp) * np.fft.fft(y)).real * DT_NS


def transient_reconstruct(obs: np.ndarray, kernel: np.ndarray,
                          lambda_reg: float) -> np.ndarray:
    n = len(obs)
    kp = np.zeros(n)
    kp[:len(kernel)] = kernel
    kf = np.fft.fft(kp)
    xf = np.conj(kf) * np.fft.fft(obs) / (
        np.abs(kf) ** 2 + lambda_reg ** 2
    ) / DT_NS
    return np.fft.ifft(xf).real


def cryoscope_observe(y: np.ndarray, seed: int, shots: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    phase = CRYO_PHASE_GAIN * np.cumsum(y) * DT_NS
    p_i = 0.5 - 0.5 * READOUT_CONTRAST * np.sin(phase)
    p_q = 0.5 + 0.5 * READOUT_CONTRAST * np.cos(phase)
    return _sample(rng, p_i, shots), _sample(rng, p_q, shots)


def cryoscope_reconstruct(obs: tuple[np.ndarray, np.ndarray], window: int) -> np.ndarray:
    p_i, p_q = obs
    phase = np.unwrap(np.arctan2(0.5 - p_i, p_q - 0.5))
    phase = savgol_filter(phase, int(window), 2)
    return np.gradient(phase, DT_NS) / CRYO_PHASE_GAIN


def observe_and_reconstruct(protocol: str, y: np.ndarray, kernel: np.ndarray,
                            seed: int, shots: int, parameter: float):
    if protocol == "transient":
        obs = transient_observe(y, kernel, seed, shots)
        return obs, transient_reconstruct(obs, kernel, float(parameter))
    if protocol == "cryoscope":
        obs = cryoscope_observe(y, seed, shots)
        return np.stack(obs), cryoscope_reconstruct(obs, int(parameter))
    raise ValueError(protocol)


def estimate_transfer(xs: list[np.ndarray], ys: list[np.ndarray],
                      ridge: float = 1e-8):
    xf = np.stack([np.fft.rfft(x) for x in xs])
    yf = np.stack([np.fft.rfft(y) for y in ys])
    power = np.sum(np.abs(xf) ** 2, axis=0)
    h = np.sum(yf * np.conj(xf), axis=0) / (
        power + ridge * np.max(power)
    )
    # Multisine probes are zero mean. DC is supplied by the independent static
    # gain calibration and normalised to unity for this benchmark.
    h[0] = 1.0 + 0.0j
    support = power >= 0.05 * np.max(power)
    support[0] = True
    return h, support, power


def calibrate_protocol_response(protocol: str, parameter: float,
                                kernel: np.ndarray, seed: int, shots: int,
                                split: int = 0):
    xs = probe_set(seed, split)
    reconstructed, raw = [], []
    for j, x in enumerate(xs):
        obs, rec = observe_and_reconstruct(
            protocol, x, kernel, seed + 100 * split + j, shots, parameter
        )
        raw.append(obs)
        reconstructed.append(rec)
    response, support, power = estimate_transfer(xs, reconstructed)
    reliable = support & (np.abs(response) >= R_FLOOR_FRACTION
                          * np.max(np.abs(response[support])))
    reliable[0] = True
    return response, reliable, xs, reconstructed, raw, power


def deembed_transfer(measured: np.ndarray, response: np.ndarray,
                     reliable: np.ndarray) -> np.ndarray:
    out = np.full_like(measured, np.nan + 1j * np.nan)
    out[reliable] = measured[reliable] / response[reliable]
    out[0] = 1.0 + 0.0j
    return out


def deembed_waveform(reconstructed: np.ndarray, response: np.ndarray,
                     reliable: np.ndarray) -> np.ndarray:
    yf = np.fft.rfft(reconstructed)
    xf = np.zeros_like(yf)
    xf[reliable] = yf[reliable] / response[reliable]
    # Cosine taper over the final four retained bins limits ringing.
    idx = np.flatnonzero(reliable)
    if len(idx) >= 5:
        tail = idx[-4:]
        xf[tail] *= np.array([0.95, 0.75, 0.40, 0.10])
    return np.fft.irfft(xf, n=len(reconstructed))


def truth_frequency_response(h: np.ndarray) -> np.ndarray:
    return np.fft.rfft(np.asarray(h, float), n=N_SAMPLES)


def qualified_mask(h_est: np.ndarray, h_true: np.ndarray,
                   reliable: np.ndarray) -> np.ndarray:
    valid = reliable & np.isfinite(h_est)
    amp = np.full(len(h_est), np.inf)
    phase = np.full(len(h_est), np.inf)
    amp[valid] = np.abs(np.abs(h_est[valid]) - np.abs(h_true[valid])) / np.maximum(
        np.abs(h_true[valid]), 1e-12
    )
    phase[valid] = np.abs(np.angle(h_est[valid] * np.conj(h_true[valid])))
    out = valid & (amp <= 0.05) & (phase <= 0.10)
    out[0] = True
    return out


def contiguous_mask(mask: np.ndarray, reliable: np.ndarray) -> np.ndarray:
    out = np.zeros_like(mask, dtype=bool)
    out[0] = True
    for idx in range(1, min(len(mask), 129)):  # up to 500 MHz
        if mask[idx] and reliable[idx]:
            out[idx] = True
        else:
            break
    return out


def bandwidth_mhz(contiguous: np.ndarray) -> float:
    idx = np.flatnonzero(contiguous)
    return float(idx[-1] / (N_SAMPLES * DT_NS) * 1000.0) if len(idx) else 0.0


def transfer_error(h_est: np.ndarray, h_true: np.ndarray,
                   reliable: np.ndarray) -> float:
    valid = reliable & np.isfinite(h_est)
    valid[0] = False
    if not np.any(valid):
        return np.inf
    amp = np.median(np.abs(np.abs(h_est[valid]) - np.abs(h_true[valid]))
                    / np.maximum(np.abs(h_true[valid]), 1e-12))
    phase = np.median(np.abs(np.angle(h_est[valid] * np.conj(h_true[valid]))))
    return float(amp + phase)


def choose_protocol_parameter(protocol: str, kernel: np.ndarray,
                              shots: int = MAIN_SHOTS) -> ProtocolChoice:
    candidates = ((0.25, 0.5, 1., 2., 5., 10., 20.) if protocol == "transient"
                  else (5., 7., 9., 11., 15.))
    summaries = []
    truth = np.ones(N_SAMPLES // 2 + 1, dtype=complex)
    for parameter in candidates:
        bws, errors = [], []
        for seed in TUNING_SEEDS:
            response, reliable, *_ = calibrate_protocol_response(
                protocol, parameter, kernel, int(seed), shots, split=0
            )
            xs = probe_set(int(seed), split=1)
            rec = []
            for j, x in enumerate(xs):
                _, y = observe_and_reconstruct(
                    protocol, x, kernel, int(seed + 1000 + j), shots, parameter
                )
                rec.append(y)
            measured, support, _ = estimate_transfer(xs, rec)
            rel = reliable & support
            h_id = deembed_transfer(measured, response, rel)
            q = qualified_mask(h_id, truth, rel)
            bws.append(bandwidth_mhz(contiguous_mask(q, rel)))
            errors.append(transfer_error(h_id, truth, rel))
        summaries.append((float(np.median(bws)), float(np.median(errors)), parameter))
    # Maximise validated continuous bandwidth, then minimise common-band error.
    best = sorted(summaries, key=lambda v: (-v[0], v[1], v[2]))[0]
    return ProtocolChoice(protocol, float(best[2]), best[0], best[1])


def design_inverse_fir(h_est: np.ndarray, band: np.ndarray, n_taps: int,
                       ridge: float) -> np.ndarray:
    idx = np.flatnonzero(band)
    idx = idx[idx > 0]
    if len(idx) < MIN_QUALIFIED_BINS:
        raise ValueError(
            f"insufficient contiguous qualified band: {len(idx)} bins; "
            f"need {MIN_QUALIFIED_BINS}"
        )
    omega = 2.0 * np.pi * np.fft.rfftfreq(N_SAMPLES, d=DT_NS) * DT_NS
    k = np.arange(n_taps)
    a = h_est[idx, None] * np.exp(-1j * omega[idx, None] * k[None, :])
    target = np.exp(-1j * omega[idx] * FILTER_DELAY)
    ar = np.vstack([a.real, a.imag])
    br = np.concatenate([target.real, target.imag])
    taps = np.linalg.solve(ar.T @ ar + ridge * np.eye(n_taps), ar.T @ br)
    dc = float(np.sum(taps))
    if abs(dc) < 1e-12:
        raise ValueError("inverse FIR has zero DC gain")
    return taps / dc


def compensate(target: np.ndarray, h_true: np.ndarray,
               taps: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    padded = np.pad(target, (0, FILTER_DELAY + len(taps)))
    awg = lfilter(taps, [1.0], padded)
    chip = lfilter(h_true, [1.0], awg)
    return awg, chip[FILTER_DELAY:FILTER_DELAY + len(target)]


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def edge_rmse(a: np.ndarray, b: np.ndarray) -> float:
    d = np.abs(np.gradient(b))
    mask = d >= 0.15 * max(float(np.max(d)), 1e-12)
    mask = np.convolve(mask.astype(int), np.ones(17, int), mode="same") > 0
    return rmse(np.asarray(a)[mask], np.asarray(b)[mask])


def hardware_metrics(awg: np.ndarray, taps: np.ndarray) -> tuple[float, float, float, bool]:
    peak = float(np.max(np.abs(awg)))
    slew = float(np.max(np.abs(np.diff(awg))) / DT_NS)
    gain = float(np.max(np.abs(np.fft.rfft(taps, 2048))))
    return peak, slew, gain, bool(peak <= 0.25 and slew <= 0.08 and gain <= 3.0)


def choose_filter(h_est: np.ndarray, band: np.ndarray, h_true: np.ndarray,
                  validation: np.ndarray) -> FilterChoice | None:
    choices = []
    for n_taps in FILTER_LENGTHS:
        for ridge in FILTER_RIDGES:
            try:
                taps = design_inverse_fir(h_est, band, n_taps, ridge)
            except ValueError:
                continue
            awg, chip = compensate(validation, h_true, taps)
            peak, slew, gain, passed = hardware_metrics(awg, taps)
            if passed:
                choices.append(FilterChoice(taps, n_taps, ridge,
                                             rmse(chip, validation),
                                             peak, slew, gain))
    return min(choices, key=lambda c: (c.validation_rmse, c.n_taps, c.ridge)) \
        if choices else None
