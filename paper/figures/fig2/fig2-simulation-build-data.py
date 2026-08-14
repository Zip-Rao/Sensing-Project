#!/usr/bin/env python3
"""Build the simulation evidence chain and render all current Fig. 2 panels."""
from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from fig2_simulation_common import DATA_DIR, FIGURE_DIR, write_json


STAGES = (
    ("response.npz", "fig2a-response-generate.py"),
    ("estimator-errors.npz", "fig2b-error-generate.py"),
    ("validation-summary.npz", "fig2c-summary-generate.py"),
    ("solver-cost.npz", "fig2d-cost-generate.py"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest() -> None:
    config_path = DATA_DIR / "config.json"
    evidence = {
        config_path.name: {
            "role": "frozen protocol and evidence conventions",
            "source_script": None,
            "sha256": sha256(config_path),
        }
    }
    for name, script in STAGES:
        path = DATA_DIR / name
        evidence[name] = {
            "role": {
                "response.npz": "frozen raw response and fitted response model",
                "estimator-errors.npz": "derived estimator outputs and errors",
                "validation-summary.npz": "derived held-out summary statistics",
                "solver-cost.npz": "derived solver-call benchmark",
            }[name],
            "source_script": script,
            "sha256": sha256(path),
        }
    csv_path = DATA_DIR / "solver-cost.csv"
    evidence[csv_path.name] = {
        "role": "human-readable projection of solver-cost.npz",
        "source_script": "fig2d-cost-generate.py",
        "sha256": sha256(csv_path),
    }
    write_json(DATA_DIR / "manifest.json", {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generation_command": f"{sys.executable} {Path(__file__).as_posix()}",
        "dependency_order": [name for name, _ in STAGES],
        "mutation_rule": "each stage writes only its own evidence product",
        "random_seed": None,
        "evidence": evidence,
    })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    for _, script in STAGES:
        command = [sys.executable, str(FIGURE_DIR / script)]
        if args.plot_only and script != "fig2c-summary-generate.py":
            command.append("--plot-only")
        subprocess.run(command, check=True)
    subprocess.run([
        sys.executable, str(FIGURE_DIR / "figS2-response-simulation-generate.py")
    ], check=True)
    subprocess.run([
        sys.executable, str(FIGURE_DIR / "fig2-simulation-evidence-assemble.py")
    ], check=True)
    write_manifest()


if __name__ == "__main__":
    main()
