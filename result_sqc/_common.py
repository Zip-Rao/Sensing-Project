"""Shared helpers for result_sqc/ reproduction scripts.

Provides the canonical qubit factory, work-point constants, and small
save / metric / overlay utilities so every generator script stays a thin
skeleton around the sqc API.

All scripts import from here and are expected to run from the repo root.
"""
from __future__ import annotations

import os
import sys

import numpy as np

# --- make repo root importable (sqc + src) regardless of CWD ------------
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# --- canonical physics parameters (match result/ src scripts) -----------
TWO_PI = 2 * np.pi
EC = 0.2 * TWO_PI
EJ = 10.0 * TWO_PI
T1 = 100.0e3      # ns
T2 = 50.0e3       # ns
N_LEVELS = 2

# optimal_work_point() of the frozen src qubit (max-sensitivity, linear)
OPTIMAL_FLUX = 0.9553166181245093

# time grid shared by all transient scripts (src: linspace(0,200,400))
T_LIST = np.linspace(0, 200, 400)

# reconstruction calibration (src harness values)
WIENER_LAMBDA = 5.0
LM_LAMBDA = 100.0
LM_N_BASIS = 100
LM_MAX_ITER = 20

# src signal_type → human name (matches generate_single_signal_data.py)
SIGNAL_NAMES = {2: "sine", 4: "step", 5: "double_peak", 7: "complex"}


# --- qubit factory ------------------------------------------------------
def make_sqc_qubit(flux=OPTIMAL_FLUX):
    """sqc TransmonQubit at the given work point (default src optimal)."""
    from sqc.devices.transmon import TransmonQubit
    return TransmonQubit(EC=EC, EJ=EJ, T1=T1, T2=T2, flux=flux,
                         n_levels=N_LEVELS)


def make_src_qubit(flux=OPTIMAL_FLUX):
    """Frozen src TransmonQubit at the given work point (for baseline)."""
    from src.qubit import TransmonQubit as SrcQubit
    q = SrcQubit(EC=EC, EJ=EJ, T1=T1, T2=T2, flux=0.0, state=0,
                 n_levels=N_LEVELS)
    q.change_flux(flux)
    return q


# --- metrics ------------------------------------------------------------
def rmse(recon, truth):
    """Root-mean-square error over the overlapping length."""
    n = min(len(recon), len(truth))
    return float(np.sqrt(np.mean((np.asarray(recon)[:n]
                                  - np.asarray(truth)[:n]) ** 2)))


def peak_ratio(recon, truth):
    """|recon|max / |truth|max — 1.0 = perfect amplitude recovery."""
    return float(np.max(np.abs(recon)) / np.max(np.abs(truth)))


# --- io -----------------------------------------------------------------
def out_dir(subdir):
    """Absolute path of a result_sqc/<subdir>, created if missing."""
    d = os.path.join(os.path.dirname(__file__), subdir)
    os.makedirs(d, exist_ok=True)
    return d


def save_npz(subdir, name, **arrays):
    """Save arrays to result_sqc/<subdir>/<name>.npz, return the path."""
    path = os.path.join(out_dir(subdir), f"{name}.npz")
    np.savez(path, **arrays)
    return path


# --- LM adapter (IMPORTANT) --------------------------------------------
def adapt_for_lm(res):
    """Adapt a TransientSensingExperiment result for the LM reconstructor.

    WHY THIS EXISTS: TransientSensingExperiment.run() emits
    ``data["p_e"]`` / ``axes["t_flux"]``, but TransientReconstruction's
    ``_reconstruct_lm`` reads ``data["p_meas"]`` / ``axes["t_signal"]``
    (a pre-existing key-name mismatch in sqc, commit 059bbe68 — NOT related
    to the wiener /dt work). Passing the raw result to method="lm" raises
    ``KeyError: 'p_meas'``. This helper renames the keys so LM runs.

    Returns a new ExperimentResult usable with
    ``TransientReconstruction(method="lm", qubit=..., control_pulse=exp.control_pulse)
      .reconstruct(adapted, kernel=res.data["kernel"])`` which returns a
    ``(FluxSignal, history)`` tuple.
    """
    from sqc.simulation.result import ExperimentResult
    return ExperimentResult(
        data={"p_meas": res.data["p_e"], "kernel": res.data["kernel"]},
        axes={"t_signal": res.axes["t_flux"]},
        metadata=dict(getattr(res, "metadata", {}) or {}),
    )
