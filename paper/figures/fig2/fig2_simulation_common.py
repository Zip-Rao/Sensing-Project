"""Shared paths and frozen protocol metadata for the simulation-only Fig. 2."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


FIGURE_DIR = Path(__file__).resolve().parent
PAPER_DIR = FIGURE_DIR.parents[1]
REPO_ROOT = PAPER_DIR.parent
DATA_DIR = PAPER_DIR / "data" / "fig2-simulation"

RESPONSE_DATA = DATA_DIR / "response.npz"
ESTIMATOR_DATA = DATA_DIR / "estimator-errors.npz"
SUMMARY_DATA = DATA_DIR / "validation-summary.npz"
COST_DATA = DATA_DIR / "solver-cost.npz"
COST_CSV = DATA_DIR / "solver-cost.csv"

TWO_PI = 2.0 * np.pi
CANDIDATE_INTERVAL_MHZ = 14.0
MAX_MODEL_ERROR = 0.01
MAX_FREQUENCY_ERROR_MHZ = 0.5
REFERENCE_FLUX_PHI0 = 0.9553166181245093


def local_flux_from_detuning_mhz(
    delta_mhz: np.ndarray, ec_rad_ghz: float, ej_rad_ghz: float,
    reference_omega_rad_ghz: float,
) -> np.ndarray:
    """Invert the analytic transmon dispersion on the local [0.5, 1] branch."""
    omega = reference_omega_rad_ghz + TWO_PI * np.asarray(delta_mhz, dtype=float) * 1e-3
    ratio = (omega + ec_rad_ghz) ** 2 / (8.0 * ej_rad_ghz * ec_rad_ghz)
    valid = (ratio >= 0.0) & (ratio <= 1.0)
    flux = np.full_like(ratio, np.nan, dtype=float)
    flux[valid] = np.arccos(-ratio[valid]) / np.pi
    return flux


def load_npz(path: Path) -> dict[str, np.ndarray]:
    """Load an NPZ without allowing object deserialization."""
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def save_npz(path: Path, data: dict[str, np.ndarray]) -> None:
    """Write one evidence product; callers may only write their own stage."""
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **data)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
