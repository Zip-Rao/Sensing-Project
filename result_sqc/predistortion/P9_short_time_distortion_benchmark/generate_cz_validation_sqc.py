#!/usr/bin/env python3
"""Apply the pre-registered gate before running an expensive CZ scan.

The validation plan forbids a gate-level advantage figure when G1/G2/G3 fail.
This script records that decision in machine-readable form and creates a small
NPZ status artifact; it deliberately does not call simulate_cz_from_flux.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def main() -> None:
    metrics_path = HERE / "metrics.json"
    if not metrics_path.exists():
        raise FileNotFoundError("Run generate_short_time_benchmark_sqc.py first")
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    gates = {name: bool(metrics[name]["passed"]) for name in ("G1", "G2", "G3")}
    entered = all(gates.values())
    status = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "entered_gate_validation": entered,
        "preregistered_gates": gates,
        "reason": (
            "G1-G3 passed; gate simulation may proceed."
            if entered else
            "Stopped before G4 as pre-registered: at least one of G1-G3 failed."
        ),
        "existing_reference": (
            "P7 remains a Cryoscope-only full-Hamiltonian reference and is not "
            "evidence for transient-derived compensation."
        ),
    }
    (HERE / "gate_validation_status.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    np.savez_compressed(
        HERE / "cz_validation_data.npz",
        entered_gate_validation=np.array([entered]),
        gate_names=np.array(list(gates)),
        gate_passed=np.array(list(gates.values()), dtype=bool),
        reason=np.array([status["reason"]]),
    )
    print(json.dumps(status, indent=2))
    if entered:
        raise RuntimeError(
            "All upstream gates passed, but this guarded script intentionally "
            "requires the full CZ scan implementation before claiming G4."
        )


if __name__ == "__main__":
    main()
