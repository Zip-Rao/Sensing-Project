"""Generate physics regression baselines.

Run before refactor work begins:

    python -m tests.regression.generate_baselines

Regenerate only when intentionally changing physics.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from src.protocal import Protocal
from src.qubit import TransmonQubit
from tests.conftest import save_baseline


def _disable_interactive_plots() -> None:
    """Keep legacy plot calls from blocking baseline generation."""
    plt.show = lambda *args, **kwargs: None


def _make_default_qubit() -> TransmonQubit:
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=2,
    )


def _baseline_qubit_static() -> dict:
    q = _make_default_qubit()
    return {
        "frequency": float(q.frequency),
        "anharmonicity": float(q.anharmonicity),
        "sensitivity_at_zero": float(q.frequency_sensitivity(0.0)),
        "sensitivity_at_optimal": float(
            q.frequency_sensitivity(q.optimal_work_point() / np.pi)
        ),
    }


def _baseline_ramsey() -> dict:
    q = _make_default_qubit()
    proto = Protocal(type=1)
    proto.initialize(q, state=0)
    phi, tau_list, p_e_list = proto.evolve(q)
    plt.close("all")
    return {
        "Phi_signal": np.asarray(phi.signal),
        "Phi_t_list": np.asarray(phi.t_list),
        "tau_list": np.asarray(tau_list),
        "p_e_list": np.asarray(p_e_list),
    }


def _baseline_diff_echo() -> dict:
    q = _make_default_qubit()
    proto = Protocal(type=2)
    proto.initialize(q, state=0)
    phi, tau_list, p_e_list, k, t_int = proto.evolve(q)
    plt.close("all")
    return {
        "Phi_signal": np.asarray(phi.signal),
        "Phi_t_list": np.asarray(phi.t_list),
        "tau_list": np.asarray(tau_list),
        "p_e_list": np.asarray(p_e_list),
        "k": int(k),
        "t_int": float(t_int),
    }


def _baseline_transient() -> dict:
    q = _make_default_qubit()
    proto = Protocal(type=4)
    proto.initialize(q, state=0)
    t_samples, kernel, scan_list, delta_p, p_e, phi, _ctrl = proto.evolve(q)
    plt.close("all")
    return {
        "t_samples": np.asarray(t_samples),
        "kernel": np.asarray(kernel),
        "scan_list": np.asarray(scan_list),
        "delta_p": np.asarray(delta_p),
        "p_e": np.asarray(p_e),
        "Phi_signal": np.asarray(phi.signal),
        "Phi_t_list": np.asarray(phi.t_list),
    }


def _baseline_cryoscope() -> dict:
    q = _make_default_qubit()
    proto = Protocal(type=5)
    proto.initialize(q, state=0)
    trunc_list, varphi, phi, p_e_list = proto.evolve(q)
    plt.close("all")
    return {
        "trunc_list": np.asarray(trunc_list),
        "varphi": np.asarray(varphi),
        "Phi_signal": np.asarray(phi.signal),
        "Phi_t_list": np.asarray(phi.t_list),
        "p_e_I": np.asarray(p_e_list[0]),
        "p_e_Q": np.asarray(p_e_list[1]),
    }


def _baseline_predistortion() -> dict:
    from sqc.calibration.predistortion import PredistortionDesigner
    from sqc.control.waveform import Waveform
    from sqc.hardware.distortion import SingleExponentialDistortion
    from sqc.workflows.predistortion_validation import PredistortionValidationWorkflow

    q = _make_default_qubit()
    t = np.linspace(0.0, 120.0, 1024)
    target = Waveform(
        t_list=t,
        samples=0.5
        * (np.tanh((t - 20.0) / 2.0) - np.tanh((t - 80.0) / 2.0)),
    )
    result = PredistortionValidationWorkflow(
        qubit=q,
        target_waveform=target,
        true_distortion=SingleExponentialDistortion(amplitude=0.15, tau=12.0),
        designer=PredistortionDesigner(method="iir_inverse"),
        fit_type="iir",
    ).run()
    return {
        "t": t,
        "target": result["target"].samples,
        "on_chip_uncorrected": result["on_chip_uncorrected"].samples,
        "awg_predistorted": result["awg_predistorted"].samples,
        "on_chip_corrected": result["on_chip_corrected"].samples,
        "metrics": result["metrics"],
    }


def _baseline_z_crosstalk() -> dict:
    from sqc.control.waveform import Waveform
    from sqc.devices.chip import ChipTopology
    from sqc.devices.transmon import TransmonQubit
    from sqc.hardware.transfer_matrix import TransferMatrix
    from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow

    t = np.linspace(0.0, 120.0, 1024)
    pulse = Waveform(
        t_list=t,
        samples=0.5
        * (np.tanh((t - 20.0) / 2.0) - np.tanh((t - 80.0) / 2.0)),
    )
    transfer = TransferMatrix.from_dc_matrix(
        dc_matrix=np.array([[1.0, 0.0], [0.04, 1.0]]),
        source_names=["QA", "QB"],
        target_names=["QA", "QB"],
    )
    chip = ChipTopology(
        qubits=[
            TransmonQubit(
                EC=2 * np.pi * 0.2,
                EJ=2 * np.pi * 15,
                T1=10000.0,
                T2=8000.0,
                n_levels=2,
                name="QA",
            ),
            TransmonQubit(
                EC=2 * np.pi * 0.2,
                EJ=2 * np.pi * 15,
                T1=10000.0,
                T2=8000.0,
                n_levels=2,
                name="QB",
            ),
        ],
        transfer_matrix=transfer,
    )
    result = ZCrosstalkWorkflow(
        chip=chip,
        flux_pulse_on_A=pulse,
        true_transfer_matrix=transfer,
    ).run()
    return {
        "t": t,
        "pulse_A": pulse.samples,
        "phi_B_true": result["phi_B_true"].samples,
        "phi_B_reconstructed": result["phi_B_reconstructed"].samples,
        "compensation_pulse": result["compensation_pulse"].samples,
        "phi_B_after_compensation": result["phi_B_after_compensation"].samples,
        "H_BA_estimated_real": result["H_BA_estimated"].real,
        "H_BA_estimated_imag": result["H_BA_estimated"].imag,
        "H_BA_true_real": result["H_BA_true"].real,
        "H_BA_true_imag": result["H_BA_true"].imag,
        "metrics": result["metrics"],
    }


BASELINES = {
    "qubit_static": _baseline_qubit_static,
    "ramsey_default": _baseline_ramsey,
    "diff_echo_default": _baseline_diff_echo,
    "transient_default": _baseline_transient,
    "cryoscope_default": _baseline_cryoscope,
    "predistortion_default": _baseline_predistortion,
    "z_crosstalk_default": _baseline_z_crosstalk,
}


def main() -> None:
    _disable_interactive_plots()
    print("Generating physics regression baselines ...")
    for name, fn in BASELINES.items():
        print(f"  - {name} ... ", end="", flush=True)
        path = save_baseline(name, fn())
        print(f"saved to {path.relative_to(Path.cwd())}")
    print("Done.")


if __name__ == "__main__":
    main()
