"""sqc.control.sequence — Pulse sequence factories.

Verbatim port of src/pulse.py factory functions:
    create_pulse, create_ramsey_pulse, create_diff_echo_pulse,
    create_echo_pulse, create_cpmg_pulse, create_cryoscope_pulse.

Also defines PulseSequence dataclass (placeholder for Phase 2).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from qutip import QobjEvo

from sqc.config import CONFIG
from .flux_signal import FluxSignal as Signal
from .pulse import Pulse, CompositePulse


# Use global dt from AWG config for all free-evolution gap time axes.
_GT = CONFIG.awg.dt


def _gap(duration: float, min_points: int = 2) -> np.ndarray:
    """Return a time axis covering [0, duration) at the global dt.

    Ensures at least ``min_points`` points when duration > 0.
    Falls back to a single-point zero-duration stub if duration ≤ 0.
    """
    N = max(int(duration / _GT), min_points if duration > 0 else 1)
    return np.arange(0, N * _GT, _GT)


# ---------------------------------------------------------------------------
# PulseSequence dataclass
# ---------------------------------------------------------------------------

@dataclass
class PulseSequence:
    """Ordered collection of pulses forming an experimental sequence.

    Phase 1: dataclass only. Full implementation in Phase 2.
    """
    pulses: list[Pulse] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Factory functions
# ---------------------------------------------------------------------------

def create_pulse(
    qubit, frame, type, t_list, omega_d, phase, angle=None, trigger=0.0,
    Omega_Q=None,
    **kwargs,
):
    """Create a single Pulse with auto-calibrated amplitude.

    Parameters
    ----------
    qubit : TransmonQubit
    frame : int
        0=lab, 1=rotating.
    type : int
        Signal type (0-8).
    t_list : array-like
        Time axis (ns).
    omega_d : float
        Drive frequency (rad*GHz).
    phase : float
        Rotation axis phase (rad).
    angle : float or None
        Target rotation angle (rad). If given, amplitude is auto-scaled.
    trigger : float
        Global start time (ns).
    Omega_Q : FluxSignal or None
        Q quadrature envelope for DRAG pulses.
    **kwargs
        Passed to Signal constructor.

    Returns
    -------
    QobjEvo
        Time-dependent Hamiltonian operator.
    """
    Omega_signal = Signal(type=type, t_list=t_list, **kwargs)
    Omega_pulse = Pulse(
        frame, omega_d, phase, Omega=Omega_signal, Omega_Q=Omega_Q,
        is_rwa=True, qubit=qubit, trigger=trigger,
    )
    current_angle = Omega_pulse.get_angle_simple()
    # Adjust amplitude for target angle
    if angle is not None:
        if abs(current_angle) < 1e-12:
            raise ValueError(
                "create_pulse: base envelope integrates to ~0 "
                f"(current_angle={current_angle:.3e}); cannot auto-scale "
                "amplitude to the requested angle. Check Signal type/amplitude."
            )
        scale = angle / current_angle
        Omega_signal.params["amplitude"] *= scale
        # Preserve the DRAG I/Q ratio: the Q-quadrature envelope must be
        # scaled by the same factor as the I envelope, otherwise the DRAG
        # correction is no longer matched to the drive amplitude.
        if (
            Omega_Q is not None
            and hasattr(Omega_Q, "update_signal")
            and "amplitude" in getattr(Omega_Q, "params", {})
        ):
            Omega_Q.update_signal(amplitude=Omega_Q.params["amplitude"] * scale)
    # Re-generate
    kwargs["amplitude"] = Omega_signal.params["amplitude"]
    signal = Signal(type=type, t_list=t_list, **kwargs)
    pulse = Pulse(
        frame, omega_d, phase, Omega=signal, Omega_Q=Omega_Q,
        is_rwa=True, qubit=qubit, trigger=trigger,
    )
    # NOTE: the returned QobjEvo is defined on the *local* t_list, so the
    # ``trigger`` offset is not baked in here — it only takes effect when the
    # pulse is embedded on a global time axis via Pulse.hamiltonian_on() /
    # CompositePulse. ``trigger`` is forwarded to the Pulse for that use.
    H_t = pulse.hamiltonian
    return QobjEvo(H_t, tlist=t_list)


def create_ramsey_pulse(
    t_rabi,
    tau,
    omega_d=0.0,
    phase1=np.pi / 2,
    phase2=0.0,
    qubit=None,
    trigger=0.0,
):
    """Build Ramsey sequence: pi/2 - tau - pi/2.

    Each sub-pulse is assigned an absolute ``trigger`` (relative to
    global t=0).  The first pi/2 starts at ``trigger``; the second
    pi/2 starts at ``trigger + (t_rabi[-1] - t_rabi[0]) + tau``.

    Parameters
    ----------
    t_rabi : array-like
        Time axis for pi/2 pulse (ns).
    tau : float
        Free evolution time (ns).
    omega_d : float
        Drive frequency (rad*GHz).
    phase1 : float
        First pi/2 phase (rad).
    phase2 : float
        Second pi/2 phase (rad).
    qubit : TransmonQubit or None
        Qubit for n_levels resolution. If None, defaults to 2-level.
    trigger : float
        Global start time of the sequence (ns).

    Returns
    -------
    CompositePulse
    """
    if tau != 0.0:
        Omega_0 = Signal(type=0, t_list=_gap(tau))
    else:
        Omega_0 = None

    Omega_1 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=(np.pi / 2.0) / (t_rabi[-1] - t_rabi[0]),
    )
    cur = float(trigger)
    pulses = []
    # First pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase1,
            Omega=Omega_1,
            is_rwa=True,
            qubit=qubit,
            trigger=cur,
        )
    )
    # Advance by the pulse *duration* (t_rabi[-1] - t_rabi[0]), consistent
    # with create_diff_echo_pulse / create_echo_pulse.  Numerically identical
    # to t_rabi[-1] when t_rabi starts at 0 (the default), but correct when
    # t_rabi[0] != 0.
    cur += float(t_rabi[-1] - t_rabi[0])
    # Gap (free evolution)
    if tau != 0.0 and Omega_0 is not None:
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=0.0,
                Omega=Omega_0,
                is_rwa=True,
                qubit=qubit,
                trigger=cur,
            )
        )
        cur += tau
    # Second pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase2,
            Omega=Omega_1,
            is_rwa=True,
            qubit=qubit,
            trigger=cur,
        )
    )
    return CompositePulse(pulses)


def create_diff_echo_pulse(
    t_rabi,
    tau,
    t_int,
    t_rep,
    k,
    omega_d,
    phase1=0.0,
    phase2=np.pi / 2,
    phase3=np.pi / 2,
    qubit=None,
    trigger=0.0,
):
    """Build differential echo sequence.

    pi/2 - [tau - pi - tau' - (tau+t_int) - pi - tau'']^k - pi/2

    Parameters
    ----------
    t_rabi : array-like
        Time axis for pulses.
    tau : float
        Free evolution (ns).
    t_int : float
        Interaction time (ns).
    t_rep : float
        Repetition period (ns).
    k : int
        Number of repetitions.
    omega_d : float
        Drive frequency.
    phase1, phase2, phase3 : float
        Pulse phases (rad).
    qubit : TransmonQubit or None
    trigger : float
        Global start time of the sequence (ns).

    Returns
    -------
    CompositePulse
    """
    Omega_0 = Signal(type=0, t_list=_gap(tau))
    Omega_01 = Signal(
        type=0,
        t_list=_gap(t_rep + t_int - t_rabi[-1] + t_rabi[0]),
    )
    Omega_02 = Signal(
        type=0,
        t_list=_gap(t_rep - tau - t_int - (t_rabi[-1] - t_rabi[0])),
    )
    Omega_1 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=(np.pi / 2.0) / (t_rabi[-1] - t_rabi[0]),
    )
    Omega_2 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=np.pi / (t_rabi[-1] - t_rabi[0]),
    )
    cur = float(trigger)
    pulses = []
    # First pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase1,
            Omega=Omega_1,
            is_rwa=True,
            qubit=qubit,
            trigger=cur,
        )
    )
    cur += float(t_rabi[-1] - t_rabi[0])
    for _ in range(k):
        # tau gap
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=0.0,
                Omega=Omega_0,
                is_rwa=True,
                qubit=qubit,
                trigger=cur,
            )
        )
        cur += tau
        # pi
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=phase2,
                Omega=Omega_2,
                is_rwa=True,
                qubit=qubit,
                trigger=cur,
            )
        )
        cur += float(t_rabi[-1] - t_rabi[0])
        # tau'
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=0.0,
                Omega=Omega_01,
                is_rwa=True,
                qubit=qubit,
                trigger=cur,
            )
        )
        cur += t_rep + t_int - t_rabi[-1] + t_rabi[0]
        # pi
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=phase2,
                Omega=Omega_2,
                is_rwa=True,
                qubit=qubit,
                trigger=cur,
            )
        )
        cur += float(t_rabi[-1] - t_rabi[0])
        # tau''
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=0.0,
                Omega=Omega_02,
                is_rwa=True,
                qubit=qubit,
                trigger=cur,
            )
        )
        cur += t_rep - tau - t_int - (t_rabi[-1] - t_rabi[0])
    # Final pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase3,
            Omega=Omega_1,
            is_rwa=True,
            qubit=qubit,
            trigger=cur,
        )
    )
    return CompositePulse(pulses)


def create_echo_pulse(
    t_rabi, tau, omega_d=0.0, phase1=0.0, phase2=0.0, phase3=0.0,
    trigger=0.0,
):
    """Build spin-echo sequence: pi/2 - tau - pi - tau - pi/2.

    Parameters
    ----------
    t_rabi : array-like
    tau : float
    omega_d : float
    phase1, phase2, phase3 : float
    trigger : float
        Global start time of the sequence (ns).

    Returns
    -------
    CompositePulse
    """
    Omega_0 = Signal(type=0, t_list=_gap(tau))
    Omega_1 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=(np.pi / 2.0) / (t_rabi[-1] - t_rabi[0]),
    )
    Omega_2 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=np.pi / (t_rabi[-1] - t_rabi[0]),
    )
    cur = float(trigger)
    pulses = []
    # pi/2
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=phase1,
            Omega=Omega_1, is_rwa=True, trigger=cur,
        )
    )
    cur += float(t_rabi[-1] - t_rabi[0])
    # tau
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=0.0,
            Omega=Omega_0, is_rwa=True, trigger=cur,
        )
    )
    cur += tau
    # pi
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=phase2,
            Omega=Omega_2, is_rwa=True, trigger=cur,
        )
    )
    cur += float(t_rabi[-1] - t_rabi[0])
    # tau
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=0.0,
            Omega=Omega_0, is_rwa=True, trigger=cur,
        )
    )
    cur += tau
    # pi/2
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=phase3,
            Omega=Omega_1, is_rwa=True, trigger=cur,
        )
    )
    return CompositePulse(pulses)


def create_cpmg_pulse(
    t_rabi, tau, n, omega_d, phase1=0.0, phase2=np.pi / 2, phase3=np.pi / 2,
    trigger=0.0,
):
    """Build CPMG sequence: pi/2 - [tau/2 - pi - tau - pi - tau/2] - pi/2.

    Parameters
    ----------
    t_rabi : array-like
    tau : float
    n : int
        Number of pi pulses.
    omega_d : float
    phase1, phase2, phase3 : float
    trigger : float
        Global start time of the sequence (ns).

    Returns
    -------
    CompositePulse
    """
    Omega_0 = Signal(type=0, t_list=_gap(tau / 2))
    Omega_1 = Signal(type=0, t_list=_gap(tau))
    Omega_2 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=(np.pi / 2.0) / (t_rabi[-1] - t_rabi[0]),
    )
    Omega_3 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=np.pi / (t_rabi[-1] - t_rabi[0]),
    )
    cur = float(trigger)
    pulses = []
    # pi/2
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=phase1,
            Omega=Omega_2, is_rwa=True, trigger=cur,
        )
    )
    cur += float(t_rabi[-1] - t_rabi[0])
    # tau/2
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=0.0,
            Omega=Omega_0, is_rwa=True, trigger=cur,
        )
    )
    cur += tau / 2
    for _ in range(n):
        # pi
        pulses.append(
            Pulse(
                frame=1, omega_d=omega_d, phase=phase2,
                Omega=Omega_3, is_rwa=True, trigger=cur,
            )
        )
        cur += float(t_rabi[-1] - t_rabi[0])
        # tau
        pulses.append(
            Pulse(
                frame=1, omega_d=omega_d, phase=0.0,
                Omega=Omega_1, is_rwa=True, trigger=cur,
            )
        )
        cur += tau
    # tau/2
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=0.0,
            Omega=Omega_0, is_rwa=True, trigger=cur,
        )
    )
    cur += tau / 2
    # pi/2
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=phase3,
            Omega=Omega_2, is_rwa=True, trigger=cur,
        )
    )
    return CompositePulse(pulses)


def create_cryoscope_pulse(
    t_rabi, tau, omega_d, phase1=np.pi / 2, phase2=0.0,
    trigger=0.0,
):
    """Build Cryoscope sequence: pi/2 - tau - pi/2.

    Parameters
    ----------
    t_rabi : array-like
    tau : float
    omega_d : float
    phase1 : float
        First pi/2 phase (rad).
    phase2 : float
        Second pi/2 phase (rad).
    trigger : float
        Global start time of the sequence (ns).

    Returns
    -------
    CompositePulse
    """
    Omega_0 = Signal(
        type=1,
        t_list=_gap(tau),
        amplitude=0.0,
    )
    Omega_1 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=(np.pi / 2.0) / (t_rabi[-1] - t_rabi[0]),
    )
    cur = float(trigger)
    pulses = []
    # pi/2-Y
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=phase1,
            Omega=Omega_1, is_rwa=True, trigger=cur,
        )
    )
    cur += float(t_rabi[-1] - t_rabi[0])
    # tau
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=0.0,
            Omega=Omega_0, is_rwa=True, trigger=cur,
        )
    )
    cur += tau
    # pi/2
    pulses.append(
        Pulse(
            frame=1, omega_d=omega_d, phase=phase2,
            Omega=Omega_1, is_rwa=True, trigger=cur,
        )
    )
    return CompositePulse(pulses)


# ---------------------------------------------------------------------------
# Pi-pulse compensation factory
# ---------------------------------------------------------------------------


def create_pi_pulse_compensation_pulse(
    t_rabi, T_pi, omega_d, phase=0.0, qubit=None, trigger=0.0,
):
    """Build a pi-pulse for the pi-pulse compensation protocol.

    This creates only the microwave drive Hamiltonian (pi-pulse).
    The compensation flux is applied separately via the qubit's flux
    signal (qubit_in_mag) — it is NOT part of this control sequence.

    Delegates to create_pulse(angle=pi).
    """
    return create_pulse(
        qubit,
        frame=1,
        type=1,
        t_list=np.asarray(t_rabi, dtype=float),
        omega_d=omega_d,
        phase=phase,
        angle=np.pi,
        trigger=trigger,
    )
