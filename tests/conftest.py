"""Global pytest fixtures and helpers."""
from __future__ import annotations

import os
import pickle
from pathlib import Path

import numpy as np
import pytest
from qutip import basis

# Path constants
PROJECT_ROOT = Path(__file__).parent.parent
BASELINE_DIR = Path(__file__).parent / "baselines"

# Numerical tolerance for physics regression
RTOL_PHYSICS = 1e-6
ATOL_PHYSICS = 1e-9


# ----------------------------- Fixtures -------------------------------------

@pytest.fixture(scope="session")
def baseline_dir() -> Path:
    """Directory holding pickled physics baselines."""
    return BASELINE_DIR


@pytest.fixture(scope="session")
def project_root() -> Path:
    return PROJECT_ROOT


@pytest.fixture(scope="session")
def qubit_default():
    """Reference TransmonQubit used for all baselines.

    Parameters chosen to match Simulation.ipynb defaults.
    DO NOT change these values; baselines are pickled against them.
    """
    from src.qubit import TransmonQubit
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=2,
    )


@pytest.fixture(scope="session")
def t_list_default():
    return np.linspace(0, 250, 500)


@pytest.fixture(scope="session")
def rng():
    """Reproducible random generator."""
    return np.random.default_rng(seed=42)


# ----------------------------- Helpers --------------------------------------

def assert_array_close(actual, expected, rtol=RTOL_PHYSICS, atol=ATOL_PHYSICS,
                        name="<unnamed>"):
    """Strict array comparison with helpful error message."""
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    if actual.shape != expected.shape:
        raise AssertionError(
            f"[{name}] shape mismatch: actual {actual.shape} vs "
            f"expected {expected.shape}"
        )
    diff = np.abs(actual - expected)
    tol = atol + rtol * np.abs(expected)
    bad = diff > tol
    if np.any(bad):
        idx = np.argmax(diff / (tol + 1e-30))
        raise AssertionError(
            f"[{name}] max relative diff {diff[bad].max():.3e} "
            f"exceeds rtol={rtol:.0e}, atol={atol:.0e}. "
            f"Worst at index {idx}: actual={actual.flat[idx]:.6g}, "
            f"expected={expected.flat[idx]:.6g}"
        )


def load_baseline(name: str) -> dict:
    """Load a pickled baseline by name (without .pkl extension)."""
    path = BASELINE_DIR / f"{name}.pkl"
    if not path.exists():
        raise FileNotFoundError(
            f"Baseline {name} not found at {path}. "
            "Run `python -m tests.regression.generate_baselines` to create it."
        )
    with open(path, "rb") as f:
        return pickle.load(f)


def save_baseline(name: str, data: dict) -> Path:
    """Pickle a baseline dict to baselines/<name>.pkl."""
    BASELINE_DIR.mkdir(parents=True, exist_ok=True)
    path = BASELINE_DIR / f"{name}.pkl"
    with open(path, "wb") as f:
        pickle.dump(data, f, protocol=5)
    return path
