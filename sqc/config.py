"""sqc.config — Global unified configuration for the sensing simulation platform.

Layered configuration matching the real cQED control stack. Parameters flow
from hardware upward: AWG sample rate determines the global time quantum `dt`,
from which all pulse and signal time axes are derived.

Usage::

    from sqc.config import CONFIG

    dt = CONFIG.awg.dt                # global time step (ns)
    t_rabi = CONFIG.pulse.t_rabi       # standard pi-pulse time axis
    qubit = TransmonQubit(**CONFIG.transmon.to_dict())
    t_global = CONFIG.pulse.t_global   # simulation time window

All time axes use ``np.arange()`` with ``dt`` derived from the AWG sample rate
(not ``np.linspace``), ensuring integer sample counts and exact dt matching real
hardware behaviour.

See ``docs/architecture.md`` §"Global Configuration" for the design rationale.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# ══════════════════════════════════════════════════════════════════════════
# 1 — Hardware: AWG & ControlLine
# ══════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class AWGConfig:
    """Arbitrary waveform generator parameters.

    The **sample rate** is the root of all time discretisation. Every
    time axis in the simulation is an integer multiple of `dt`.

    Typical values for common AWGs:
      - Keysight M8190A: 12 GSa/s → dt = 0.083 ns
      - Tektronix AWG70001A: 50 GSa/s → dt = 0.020 ns
      - Budget / slower: 2 GSa/s → dt = 0.5 ns (project default)
    """

    sample_rate: float = 2.0          #: GSa/s
    voltage_range: float = 2.0        #: Vpp
    resolution: int = 14              #: bits

    @property
    def dt(self) -> float:
        """Global time step (ns). Derived from sample_rate in GSa/s.

        dt = 1 / sample_rate.  For 2 GSa/s, dt = 0.5 ns.
        """
        return 1.0 / self.sample_rate


@dataclass(frozen=True)
class ControlLineDefaults:
    """Default physical parameters for a flux (Z) control line."""
    impedance: float = 50.0           #: Ohm
    attenuation_db: float = 20.0      #: dB
    delay: float = 0.0                #: ns propagation delay


# ══════════════════════════════════════════════════════════════════════════
# 2 — Devices: Transmon qubit defaults
# ══════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class TransmonDefaults:
    """Default physical parameters for a flux-tunable Transmon qubit.

    All values are chosen to lie within the recommended ranges from
    Gao 2021 §III.B (see docstring of ``QubitSpec`` for details).
    """

    EC: float = 0.2                   #: GHz·2π (→ EC/h = 200 MHz)
    EJ: float = 15.0                  #: GHz·2π (→ EJ/h = 15 GHz)
    T1: float = 10_000.0              #: ns
    T2: float = 8_000.0               #: ns
    flux_bias: float = 0.0            #: Φ₀
    n_levels: int = 3                 #: Fock truncation

    # Recommended ranges (Gao 2021 §III.B)
    F01_RANGE: tuple = (4.0, 8.0)     #: GHz
    EJ_EC_RANGE: tuple = (40, 80)
    ALPHA_RANGE: tuple = (200, 300)   #: MHz

    def to_dict(self) -> dict:
        """Return kwargs dict for ``TransmonQubit(...)`` constructor."""
        return dict(EC=2 * np.pi * self.EC,
                    EJ=2 * np.pi * self.EJ,
                    T1=self.T1, T2=self.T2,
                    flux=self.flux_bias,
                    n_levels=self.n_levels)


# ══════════════════════════════════════════════════════════════════════════
# 3 — Control: pulse & time-axis parameters
# ══════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class PulseConfig:
    """Pulse timing and time-axis configuration.

    All time axes are computed from ``dt`` (from ``AWGConfig``) via
    ``np.arange()`` so that every point lands on an AWG sample boundary.
    """

    dt: float                          #: from AWGConfig.dt

    # ── pulse durations ────────────────────────────────────────────────
    t_rabi_duration: float = 10.0      #: ns — standard π-pulse window
    t_pi2_duration: float = 10.0       #: ns — standard π/2-pulse window
    # (π/2 uses the same time window as π; only the amplitude is halved.)

    # ── global simulation window ───────────────────────────────────────
    t_global_start: float = -50.0      #: ns — time before first pulse
    t_global_end: float = 400.0        #: ns — end of simulation window

    # ── flux-signal window ─────────────────────────────────────────────
    t_signal_duration: float = 250.0   #: ns — default flux signal duration
    t_signal_start: float = 0.0        #: ns

    # ── Ramsey / echo scan stride ──────────────────────────────────────
    tau_stride: int = 1                #: down-sampling factor for tau_list

    @property
    def t_rabi(self) -> np.ndarray:
        """Standard π-pulse time axis (N = duration/dt points)."""
        return np.arange(0, self.t_rabi_duration, self.dt)

    @property
    def t_pi2(self) -> np.ndarray:
        """Standard π/2-pulse time axis."""
        return np.arange(0, self.t_pi2_duration, self.dt)

    @property
    def t_global(self) -> np.ndarray:
        """Full simulation time window [t_start, t_end) at dt resolution."""
        return np.arange(self.t_global_start, self.t_global_end, self.dt)

    @property
    def t_signal(self) -> np.ndarray:
        """Default flux-signal time axis [0, duration) at dt resolution."""
        return np.arange(self.t_signal_start, self.t_signal_duration, self.dt)

    @property
    def tau_list(self) -> np.ndarray:
        """Default Ramsey free-precession scan list."""
        return self.t_signal[::self.tau_stride]

    # ── convenience: linspace-style aliases for backward compat ────────
    def make_time(self, t_start: float, t_end: float) -> np.ndarray:
        """Return ``arange(t_start, t_end, dt)``."""
        return np.arange(t_start, t_end, self.dt)

    def make_time_with_count(self, t_start: float, n_points: int) -> np.ndarray:
        """Return ``arange(t_start, t_start + n*dt, dt)`` with *n_points*."""
        return np.arange(t_start, t_start + n_points * self.dt, self.dt)


# ══════════════════════════════════════════════════════════════════════════
# 4 — Simulation: QuTiP solver defaults
# ══════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class SimulationConfig:
    """Default QuTiP mesolve options."""
    store_states: bool = False
    atol: float = 1e-8
    rtol: float = 1e-6


# ══════════════════════════════════════════════════════════════════════════
# 5 — Reconstruction: algorithm default hyper-parameters
# ══════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ReconstructionConfig:
    """Default hyper-parameters for waveform reconstruction algorithms."""

    # Wiener deconvolution. 5.0 matches the frozen src reference calibration
    # (src/.../generate_single_signal_data.py uses lambdas=5.0); higher values
    # over-regularise and suppress the recovered amplitude.
    lambda_reg: float = 5.0

    # Kernel estimation (stimulus perturbation)
    stim_amplitude: float = 0.0215     #: Φ₀ — hardcoded legacy default
    stim_width: float = 3.0            #: ns

    # Levenberg-Marquardt
    lm_n_basis: int = 100
    lm_max_iter: int = 10
    lm_tol: float = 1e-6
    lm_mu_init: float = 1e-3
    lm_basis_type: str = "fourier"
    lm_lambda: float = 100.0

    # Cryoscope
    cryoscope_tau: float = 100.0       #: ns — calibration square-pulse length

    # delay Ramsey
    delay_ramsey_tau: float = 20.0     #: ns — default Ramsey free evolution for delay Ramsey

    # pi-pulse compensation
    pi_pulse_T_pi: float = 10.0        #: ns — default pi-pulse width for compensation


# ══════════════════════════════════════════════════════════════════════════
# 6 — Top-level Config (singleton)
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class Config:
    """Global configuration aggregating all sub-system defaults.

    All time axes in the project should ultimately be derived from
    ``config.awg.dt`` via ``config.pulse.*`` properties.  Per-experiment
    overrides are encouraged—pass a custom ``PulseConfig(t_rabi_duration=20)``
    rather than calling ``np.linspace`` directly.
    """

    awg: AWGConfig = field(default_factory=AWGConfig)
    control_line: ControlLineDefaults = field(default_factory=ControlLineDefaults)
    transmon: TransmonDefaults = field(default_factory=TransmonDefaults)
    pulse: PulseConfig = field(default_factory=lambda: PulseConfig(dt=AWGConfig().dt))
    simulation: SimulationConfig = field(default_factory=SimulationConfig)
    reconstruction: ReconstructionConfig = field(default_factory=ReconstructionConfig)


# singleton — import this everywhere
CONFIG = Config()


# ══════════════════════════════════════════════════════════════════════════
# Convenience: configure AWG sample_rate and have all pulse times follow
# ══════════════════════════════════════════════════════════════════════════

def reconfigure(
    # ── AWG ────────────────────────────────────────────────────────────
    sample_rate: float | None = None,
    # ── Pulse ──────────────────────────────────────────────────────────
    t_rabi_duration: float | None = None,
    t_global_start: float | None = None,
    t_global_end: float | None = None,
    # ── Reconstruction ─────────────────────────────────────────────────
    lambda_reg: float | None = None,
    stim_amplitude: float | None = None,
    stim_width: float | None = None,
    lm_n_basis: int | None = None,
    lm_max_iter: int | None = None,
    lm_tol: float | None = None,
    lm_mu_init: float | None = None,
    lm_basis_type: str | None = None,
    lm_lambda: float | None = None,
    cryoscope_tau: float | None = None,
    delay_ramsey_tau: float | None = None,
    pi_pulse_T_pi: float | None = None,
    # ── Simulation ─────────────────────────────────────────────────────
    atol: float | None = None,
    rtol: float | None = None,
    # ── Transmon ───────────────────────────────────────────────────────
    EC: float | None = None,
    EJ: float | None = None,
    T1: float | None = None,
    T2: float | None = None,
    flux_bias: float | None = None,
    n_levels: int | None = None,
    # ── ControlLine ────────────────────────────────────────────────────
    impedance: float | None = None,
    attenuation_db: float | None = None,
    delay: float | None = None,
    # ── Generic fall-through ───────────────────────────────────────────
    **kwargs,
) -> Config:
    """Return a **new** ``Config`` with updated parameters.

    The global ``CONFIG`` singleton is *not* modified.  Pass only the
    parameters you wish to change; all others inherit from the current
    ``CONFIG``.

    Examples
    --------
    Single-point tuning::

        cfg = reconfigure(sample_rate=4.0, lambda_reg=5.0, n_levels=3)

    Then use ``cfg.pulse.t_rabi``, ``cfg.reconstruction.*``,
    ``cfg.transmon.to_dict()`` etc. downstream::

        qubit = TransmonQubit(**cfg.transmon.to_dict())
        exp = RamseyExperiment(qubit=qubit, t_rabi=cfg.pulse.t_rabi)

    Verification
    ------------
    - ``reconfigure(lambda_reg=5.0).reconstruction.lambda_reg == 5.0``
    - ``reconfigure().reconstruction == CONFIG.reconstruction``
    - ``CONFIG`` itself is unchanged.
    """
    # ── AWG ────────────────────────────────────────────────────────────
    awg = AWGConfig(
        sample_rate=sample_rate if sample_rate is not None else CONFIG.awg.sample_rate,
        voltage_range=kwargs.pop("voltage_range", CONFIG.awg.voltage_range),
        resolution=kwargs.pop("resolution", CONFIG.awg.resolution),
    )

    # ── Pulse ──────────────────────────────────────────────────────────
    pulse = PulseConfig(
        dt=awg.dt,
        t_rabi_duration=t_rabi_duration if t_rabi_duration is not None else CONFIG.pulse.t_rabi_duration,
        t_global_start=t_global_start if t_global_start is not None else CONFIG.pulse.t_global_start,
        t_global_end=t_global_end if t_global_end is not None else CONFIG.pulse.t_global_end,
    )

    # ── Reconstruction ─────────────────────────────────────────────────
    reconstruction = CONFIG.reconstruction
    rec_fields = {}
    if lambda_reg is not None:
        rec_fields["lambda_reg"] = lambda_reg
    if stim_amplitude is not None:
        rec_fields["stim_amplitude"] = stim_amplitude
    if stim_width is not None:
        rec_fields["stim_width"] = stim_width
    if lm_n_basis is not None:
        rec_fields["lm_n_basis"] = lm_n_basis
    if lm_max_iter is not None:
        rec_fields["lm_max_iter"] = lm_max_iter
    if lm_tol is not None:
        rec_fields["lm_tol"] = lm_tol
    if lm_mu_init is not None:
        rec_fields["lm_mu_init"] = lm_mu_init
    if lm_basis_type is not None:
        rec_fields["lm_basis_type"] = lm_basis_type
    if lm_lambda is not None:
        rec_fields["lm_lambda"] = lm_lambda
    if cryoscope_tau is not None:
        rec_fields["cryoscope_tau"] = cryoscope_tau
    if delay_ramsey_tau is not None:
        rec_fields["delay_ramsey_tau"] = delay_ramsey_tau
    if pi_pulse_T_pi is not None:
        rec_fields["pi_pulse_T_pi"] = pi_pulse_T_pi
    if rec_fields:
        reconstruction = CONFIG.reconstruction.__class__(
            **{f.name: rec_fields.get(f.name, getattr(CONFIG.reconstruction, f.name))
               for f in CONFIG.reconstruction.__dataclass_fields__.values()}
        )

    # ── Simulation ─────────────────────────────────────────────────────
    simulation = CONFIG.simulation
    sim_fields = {}
    if atol is not None:
        sim_fields["atol"] = atol
    if rtol is not None:
        sim_fields["rtol"] = rtol
    if sim_fields:
        simulation = CONFIG.simulation.__class__(
            **{f.name: sim_fields.get(f.name, getattr(CONFIG.simulation, f.name))
               for f in CONFIG.simulation.__dataclass_fields__.values()}
        )

    # ── Transmon ───────────────────────────────────────────────────────
    transmon = CONFIG.transmon
    tmon_fields = {}
    if EC is not None:
        tmon_fields["EC"] = EC
    if EJ is not None:
        tmon_fields["EJ"] = EJ
    if T1 is not None:
        tmon_fields["T1"] = T1
    if T2 is not None:
        tmon_fields["T2"] = T2
    if flux_bias is not None:
        tmon_fields["flux_bias"] = flux_bias
    if n_levels is not None:
        tmon_fields["n_levels"] = n_levels
    if tmon_fields:
        transmon = CONFIG.transmon.__class__(
            **{f.name: tmon_fields.get(f.name, getattr(CONFIG.transmon, f.name))
               for f in CONFIG.transmon.__dataclass_fields__.values()
               if f.name not in ("F01_RANGE", "EJ_EC_RANGE", "ALPHA_RANGE")}
        )

    # ── ControlLine ────────────────────────────────────────────────────
    control_line = CONFIG.control_line
    cl_fields = {}
    if impedance is not None:
        cl_fields["impedance"] = impedance
    if attenuation_db is not None:
        cl_fields["attenuation_db"] = attenuation_db
    if delay is not None:
        cl_fields["delay"] = delay
    if cl_fields:
        control_line = CONFIG.control_line.__class__(
            **{f.name: cl_fields.get(f.name, getattr(CONFIG.control_line, f.name))
               for f in CONFIG.control_line.__dataclass_fields__.values()}
        )

    return Config(
        awg=awg,
        pulse=pulse,
        reconstruction=reconstruction,
        simulation=simulation,
        transmon=transmon,
        control_line=control_line,
    )
