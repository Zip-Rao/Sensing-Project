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

from .flux_signal import FluxSignal as Signal
from .pulse import Pulse, CompositePulse


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
    qubit, frame, type, t_list, omega_d, phase, angle=None, **kwargs
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
    **kwargs
        Passed to Signal constructor.

    Returns
    -------
    QobjEvo
        Time-dependent Hamiltonian operator.
    """
    Omega_signal = Signal(type=type, t_list=t_list, **kwargs)
    Omega_pulse = Pulse(
        frame, omega_d, phase, Omega=Omega_signal, is_rwa=True, qubit=qubit
    )
    current_angle = Omega_pulse.get_angle_simple()
    # Adjust amplitude for target angle
    if angle is not None:
        Omega_signal.params["amplitude"] *= angle / current_angle
    # Re-generate
    kwargs["amplitude"] = Omega_signal.params["amplitude"]
    signal = Signal(type=type, t_list=t_list, **kwargs)
    pulse = Pulse(frame, omega_d, phase, Omega=signal, is_rwa=True, qubit=qubit)
    H_t = pulse.hamiltonian
    return QobjEvo(H_t, tlist=t_list)


def create_ramsey_pulse(
    t_rabi, tau, omega_d=0.0, phase1=np.pi / 2, phase2=0.0
):
    """Build Ramsey sequence: pi/2 - tau - pi/2.

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

    Returns
    -------
    CompositePulse
    """
    if tau != 0.0:
        Omega_0 = Signal(type=0, t_list=np.linspace(0, tau, 100))
    else:
        Omega_0 = None

    Omega_1 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=(np.pi / 2.0) / (t_rabi[-1] - t_rabi[0]),
    )
    pulses = []
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase1,
            Omega=Omega_1,
            is_rwa=True,
        )
    )
    if tau != 0.0 and Omega_0 is not None:
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=0.0,
                Omega=Omega_0,
                is_rwa=True,
            )
        )
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase2,
            Omega=Omega_1,
            is_rwa=True,
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

    Returns
    -------
    CompositePulse
    """
    Omega_0 = Signal(type=0, t_list=np.linspace(0, tau, 100))
    Omega_01 = Signal(
        type=0,
        t_list=np.linspace(
            0, t_rep + t_int - t_rabi[-1] + t_rabi[0], 100
        ),
    )
    Omega_02 = Signal(
        type=0,
        t_list=np.linspace(
            0, t_rep - tau - t_int - (t_rabi[-1] - t_rabi[0]), 100
        ),
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
    pulses = []
    # First pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase1,
            Omega=Omega_1,
            is_rwa=True,
        )
    )
    for _ in range(k):
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=0.0,
                Omega=Omega_0,
                is_rwa=True,
            )
        )
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=phase2,
                Omega=Omega_2,
                is_rwa=True,
            )
        )
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=0.0,
                Omega=Omega_01,
                is_rwa=True,
            )
        )
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=phase2,
                Omega=Omega_2,
                is_rwa=True,
            )
        )
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=0.0,
                Omega=Omega_02,
                is_rwa=True,
            )
        )
    # Final pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase3,
            Omega=Omega_1,
            is_rwa=True,
        )
    )
    return CompositePulse(pulses)


def create_echo_pulse(
    t_rabi, tau, omega_d=0.0, phase1=0.0, phase2=0.0, phase3=0.0
):
    """Build spin-echo sequence: pi/2 - tau - pi - tau - pi/2.

    Parameters
    ----------
    t_rabi : array-like
    tau : float
    omega_d : float
    phase1, phase2, phase3 : float

    Returns
    -------
    CompositePulse
    """
    Omega_0 = Signal(type=0, t_list=np.linspace(0, tau, 100))
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
    pulses = []
    # pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase1,
            Omega=Omega_1,
            is_rwa=True,
        )
    )
    # tau
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=0.0,
            Omega=Omega_0,
            is_rwa=True,
        )
    )
    # pi
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase2,
            Omega=Omega_2,
            is_rwa=True,
        )
    )
    # tau
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=0.0,
            Omega=Omega_0,
            is_rwa=True,
        )
    )
    # pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase3,
            Omega=Omega_1,
            is_rwa=True,
        )
    )
    return CompositePulse(pulses)


def create_cpmg_pulse(
    t_rabi, tau, n, omega_d, phase1=0.0, phase2=np.pi / 2, phase3=np.pi / 2
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

    Returns
    -------
    CompositePulse
    """
    Omega_0 = Signal(type=0, t_list=np.linspace(0, tau / 2, 100))
    Omega_1 = Signal(type=0, t_list=np.linspace(0, tau, 100))
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
    pulses = []
    # pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase1,
            Omega=Omega_2,
            is_rwa=True,
        )
    )
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=0.0,
            Omega=Omega_0,
            is_rwa=True,
        )
    )
    for _ in range(n):
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=phase2,
                Omega=Omega_3,
                is_rwa=True,
            )
        )
        pulses.append(
            Pulse(
                frame=1,
                omega_d=omega_d,
                phase=0.0,
                Omega=Omega_1,
                is_rwa=True,
            )
        )
    # tau/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=0.0,
            Omega=Omega_0,
            is_rwa=True,
        )
    )
    # pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase3,
            Omega=Omega_2,
            is_rwa=True,
        )
    )
    return CompositePulse(pulses)


def create_cryoscope_pulse(
    t_rabi, tau, omega_d, phase1=np.pi / 2, phase2=0.0
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

    Returns
    -------
    CompositePulse
    """
    Omega_0 = Signal(
        type=1,
        t_list=np.linspace(0, tau, 100),
        amplitude=0.0,
    )
    Omega_1 = Signal(
        type=1,
        t_list=t_rabi,
        amplitude=(np.pi / 2.0) / (t_rabi[-1] - t_rabi[0]),
    )
    pulses = []
    # pi/2-Y
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase1,
            Omega=Omega_1,
            is_rwa=True,
        )
    )
    # tau
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=0.0,
            Omega=Omega_0,
            is_rwa=True,
        )
    )
    # pi/2
    pulses.append(
        Pulse(
            frame=1,
            omega_d=omega_d,
            phase=phase2,
            Omega=Omega_1,
            is_rwa=True,
        )
    )
    return CompositePulse(pulses)
