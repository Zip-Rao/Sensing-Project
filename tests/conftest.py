"""Global pytest fixtures and helpers for physics regression tests."""
from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np

try:
    import pytest
except ModuleNotFoundError:  # pragma: no cover - only used before P0 deps install.
    class _PytestShim:
        @staticmethod
        def fixture(*args, **kwargs):
            def decorator(func):
                return func

            return decorator

    pytest = _PytestShim()


PROJECT_ROOT = Path(__file__).parent.parent
BASELINE_DIR = Path(__file__).parent / "baselines"

RTOL_PHYSICS = 1e-6
ATOL_PHYSICS = 1e-9


@pytest.fixture(scope="session")
def baseline_dir() -> Path:
    """Directory holding pickled physics baselines."""
    return BASELINE_DIR


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Repository root path."""
    return PROJECT_ROOT


@pytest.fixture(scope="session")
def qubit_default():
    """Reference TransmonQubit used for all baselines."""
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
def t_list_default() -> np.ndarray:
    """Reference time axis for Ramsey-style baselines."""
    return np.linspace(0, 250, 500)


@pytest.fixture(scope="session")
def rng() -> np.random.Generator:
    """Reproducible random generator."""
    return np.random.default_rng(seed=42)


def assert_array_close(
    actual,
    expected,
    rtol: float = RTOL_PHYSICS,
    atol: float = ATOL_PHYSICS,
    name: str = "<unnamed>",
) -> None:
    """Compare arrays with physics-regression tolerances."""
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
        idx = int(np.argmax(diff / (tol + 1e-30)))
        raise AssertionError(
            f"[{name}] max diff {diff[bad].max():.3e} exceeds "
            f"rtol={rtol:.0e}, atol={atol:.0e}. "
            f"Worst at index {idx}: actual={actual.flat[idx]:.6g}, "
            f"expected={expected.flat[idx]:.6g}"
        )


def load_baseline(name: str) -> dict:
    """Load a pickled baseline by name, without the .pkl suffix."""
    path = BASELINE_DIR / f"{name}.pkl"
    if not path.exists():
        raise FileNotFoundError(
            f"Baseline {name} not found at {path}. "
            "Run `python -m tests.regression.generate_baselines` first."
        )
    with open(path, "rb") as f:
        return pickle.load(f)


def save_baseline(name: str, data: dict) -> Path:
    """Pickle a baseline dict to tests/baselines/<name>.pkl."""
    BASELINE_DIR.mkdir(parents=True, exist_ok=True)
    path = BASELINE_DIR / f"{name}.pkl"
    with open(path, "wb") as f:
        pickle.dump(data, f, protocol=5)
    return path
