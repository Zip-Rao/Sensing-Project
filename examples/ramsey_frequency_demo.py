"""Gaussian Ramsey sequence, fringe measurement, and frequency validation."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from sqc.calibration.frequency import FrequencyMeasurement
from sqc.config import CONFIG
from sqc.control.flux_signal import FluxSignal
from sqc.control.sequence import create_ramsey_pulse
from sqc.devices import TransmonQubit
from sqc.experiments import RamseyExperiment


OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def make_qubit() -> TransmonQubit:
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15.0,
        T1=10_000,
        T2=8_000,
        n_levels=2,
    )


def fit_fringe_frequency(tau_ns: np.ndarray, p_e: np.ndarray) -> float:
    """Return the unsigned Ramsey detuning in GHz from p_e(tau)."""
    dt_ns = float(tau_ns[1] - tau_ns[0])
    centered = p_e - np.mean(p_e)
    frequencies = np.fft.rfftfreq(len(tau_ns), d=dt_ns)
    spectrum = np.abs(np.fft.rfft(centered))
    peak = int(np.argmax(spectrum[1:]) + 1)
    initial_frequency = float(frequencies[peak])

    def fringe(t, offset, amplitude, frequency, phase):
        return offset + amplitude * np.cos(2 * np.pi * frequency * t + phase)

    params, _ = curve_fit(
        fringe,
        tau_ns,
        p_e,
        p0=(float(np.mean(p_e)), 0.5, initial_frequency, 0.0),
        bounds=([0.0, -1.0, 0.0, -2 * np.pi], [1.0, 1.0, 0.5, 2 * np.pi]),
        maxfev=20_000,
    )
    return float(params[2])


def plot_sequence(qubit: TransmonQubit, t_rabi: np.ndarray) -> Path:
    sequence = create_ramsey_pulse(
        t_rabi=t_rabi,
        tau=40.0,
        omega_d=qubit.frequency,
        qubit=qubit,
        rotation_angle=np.pi / 2,
        envelope="gaussian",
        envelope_sigma=2.0,
    )
    figure, axis = plt.subplots(figsize=(8, 3.5))
    for pulse in sequence.pulses:
        time = pulse.trigger + pulse.Omega.t_list
        axis.plot(time, pulse.Omega.signal, color="tab:blue", linewidth=2)
    axis.set(xlabel="Time (ns)", ylabel="Rabi rate (rad/ns)", title="Gaussian Ramsey sequence")
    axis.grid(alpha=0.25)
    figure.tight_layout()
    path = OUTPUT_DIR / "ramsey_gaussian_sequence.png"
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return path


def run_ramsey(qubit: TransmonQubit):
    dt = CONFIG.awg.dt
    t_rabi = CONFIG.pulse.make_time(0, 10)
    tau_list = np.arange(0.0, 200.0, 4 * dt)
    t_global = CONFIG.pulse.make_time(-10, 230)
    zero_flux = FluxSignal(type=0, t_list=CONFIG.pulse.make_time(0, 230))
    known_detuning_ghz = 0.015
    omega_d = qubit.frequency - 2 * np.pi * known_detuning_ghz

    result = RamseyExperiment(
        qubit=qubit,
        flux_signal=zero_flux,
        omega_d=omega_d,
        t_rabi=t_rabi,
        tau_list=tau_list,
        t_global=t_global,
        envelope="gaussian",
        envelope_sigma=2.0,
    ).run()

    figure, axis = plt.subplots(figsize=(8, 3.5))
    axis.plot(result.axes["tau"], result.data["p_e"], "o-", markersize=3)
    axis.set(xlabel="Free-evolution time tau (ns)", ylabel="Excited-state population", title="Ramsey fringes")
    axis.grid(alpha=0.25)
    figure.tight_layout()
    path = OUTPUT_DIR / "ramsey_fringes.png"
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return result, known_detuning_ghz, omega_d, t_rabi, t_global, path


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    qubit = make_qubit()
    sequence_path = plot_sequence(qubit, CONFIG.pulse.make_time(0, 10))
    result, expected_detuning, omega_d, t_rabi, t_global, fringe_path = run_ramsey(qubit)

    fitted_detuning = fit_fringe_frequency(result.axes["tau"], result.data["p_e"])
    measurement = FrequencyMeasurement(
        qubit=qubit,
        method="ramsey",
        omega_d=omega_d,
        tau_list=result.axes["tau"],
        t_rabi=t_rabi,
        t_global=t_global,
        f_artificial=0.1,
        envelope="gaussian",
        envelope_sigma=2.0,
    )
    measured_omega = measurement.measure()
    expected_omega = qubit.frequency
    fringe_error_mhz = abs(fitted_detuning - expected_detuning) * 1_000
    frequency_error_mhz = abs(measured_omega - expected_omega) / (2 * np.pi) * 1_000

    print(f"Expected detuning: {expected_detuning * 1_000:.6f} MHz")
    print(f"Fringe-fit detuning: {fitted_detuning * 1_000:.6f} MHz")
    print(f"Expected f01: {expected_omega / (2 * np.pi):.9f} GHz")
    print(f"Measured f01: {measured_omega / (2 * np.pi):.9f} GHz")
    print(f"Fringe-fit error: {fringe_error_mhz:.6f} MHz")
    print(f"FrequencyMeasurement error: {frequency_error_mhz:.6f} MHz")
    print(f"Sequence figure: {sequence_path}")
    print(f"Fringe figure: {fringe_path}")

    assert fringe_error_mhz < 0.5
    assert frequency_error_mhz < 1.0


if __name__ == "__main__":
    main()
