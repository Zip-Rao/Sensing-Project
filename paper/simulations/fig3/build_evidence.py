#!/usr/bin/env python3
"""Build independent simulation evidence for Fig. 3 and Figs. S5/S6/S8."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "paper" / "data" / "fig3-simulation"
UPSTREAM = REPO / "paper" / "data" / "fig2-simulation"
sys.path.insert(0, str(REPO))

from result_sqc import _common as C  # noqa: E402
from sqc.config import CONFIG  # noqa: E402
from sqc.workflows.frequency_backends import (  # noqa: E402
    FiniteShotSQCExecutor,
    ProcessSimulationExecutor,
    SQCExecutor,
)
from sqc.workflows.frequency_runtime import FrequencyCalibrationRuntime  # noqa: E402
from sqc.workflows.frequency_state_machine import FrequencyCalibrationConfig  # noqa: E402
from sqc.workflows.frequency_state_machine import TrackFrequency  # noqa: E402


class FixedDriveExecutor:
    """Ablation wrapper that freezes only the Track measurement drive."""

    def __init__(self, inner, fixed_drive: float):
        self.inner = inner
        self.fixed_drive = fixed_drive

    def estimate_cost(self, cmd):
        return self.inner.estimate_cost(cmd)

    def set_machine_state(self, state):
        setter = getattr(self.inner, "set_machine_state", None)
        if callable(setter):
            setter(state)

    def execute(self, cmd):
        if isinstance(cmd, TrackFrequency):
            cmd = replace(cmd, predicted_drive=self.fixed_drive)
        return self.inner.execute(cmd)

    def rng_state(self):
        getter = getattr(self.inner, "rng_state", None)
        return getter() if callable(getter) else {}

    def restore_rng_state(self, state):
        restore = getattr(self.inner, "restore_rng_state", None)
        if callable(restore):
            restore(state)

    def provenance(self):
        return {
            "backend": type(self).__name__,
            "inner_backend": type(self.inner).__name__,
            "fixed_drive": self.fixed_drive,
        }


class SimulatedCrashExecutor:
    """Raise at a command boundary after a fixed number of completed calls."""

    def __init__(self, inner, completed_calls: int):
        self.inner = inner
        self.completed_calls = completed_calls
        self.calls = 0

    def estimate_cost(self, cmd):
        return self.inner.estimate_cost(cmd)

    def set_machine_state(self, state):
        setter = getattr(self.inner, "set_machine_state", None)
        if callable(setter):
            setter(state)

    def execute(self, cmd):
        if self.calls >= self.completed_calls:
            raise RuntimeError("intentional checkpoint boundary")
        event = self.inner.execute(cmd)
        self.calls += 1
        return event

    def provenance(self):
        return {
            "backend": type(self).__name__,
            "inner_backend": type(self.inner).__name__,
            "completed_calls_before_crash": self.completed_calls,
        }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, default=str), encoding="utf-8")


def pulse_protocol() -> dict:
    t_rabi = CONFIG.pulse.t_rabi.copy()
    u = (t_rabi - t_rabi[0]) / (t_rabi[-1] - t_rabi[0])
    return {"t_rabi": t_rabi, "rotation_angle": np.pi / 2,
            "envelope": np.sin(np.pi * u) ** 2}


def make_config() -> FrequencyCalibrationConfig:
    return FrequencyCalibrationConfig(
        epsilon_hold=2 * np.pi * 0.05e-3,
        epsilon_final=2 * np.pi * 0.5e-3,
        epsilon_enter=2 * np.pi * 2e-3,
        Delta_val=2 * np.pi * 14e-3,
        N_verify=2, max_commands=20, max_solver_calls=20_000,
        stop_after_lock_cycles=1, first_bias_step=0.005,
        max_bias_step=0.01, S_min=0.01, S_max=1e4,
    )


def finite_shots_by_role(strategy: str, verify_shots: int) -> dict[str, int]:
    return {
        "acquire": 1024,
        "track": 1024 if strategy == "short_pulse" else verify_shots,
        "verify": verify_shots,
        "monitor": verify_shots,
    }


def make_executor(strategy: str, qubit, config, target, *, seed: int | None = None,
                  shots: int = 1024, assignment_error: float = 0.0,
                  drift_mhz: float = 0.0, jump_mhz: float = 0.0,
                  drive_mode: str = "track"):
    """Construct a Track executor without exposing analytic truth to it."""
    track = ({"order": 3, "g3_source": "kernel_full",
              "estimated_solver_calls": 200}
             if strategy == "short_pulse"
             else {"order": 1, "estimated_solver_calls": 10})
    kwargs = dict(
        qubit=qubit, config=config, f_target=target,
        pulse_protocol=pulse_protocol(), track_measurement=track,
        track_method="transient" if strategy == "short_pulse" else "ramsey",
    )
    if seed is None:
        inner = SQCExecutor(**kwargs)
    else:
        shots_by_role = finite_shots_by_role(strategy, shots)
        inner = FiniteShotSQCExecutor(
            **kwargs, master_seed=seed, shots_per_circuit=shots,
            shots_by_role=shots_by_role,
            assignment_p01=assignment_error, assignment_p10=assignment_error,
        )
    if drive_mode == "fixed":
        inner = FixedDriveExecutor(inner, target)
    elif drive_mode != "track":
        raise ValueError(f"unknown drive mode: {drive_mode}")
    if drift_mhz or jump_mhz:
        inner = ProcessSimulationExecutor(
            inner=inner, drift_per_command=2 * np.pi * drift_mhz * 1e-3,
            jump_schedule={5: 2 * np.pi * jump_mhz * 1e-3} if jump_mhz else {},
        )
    return inner


def run_case(strategy: str, target_mhz: float, run_id: str, *,
             seed: int | None = None, shots: int = 1024,
             assignment_error: float = 0.0, drift_mhz: float = 0.0,
             jump_mhz: float = 0.0, max_solver_calls: int = 20_000,
             drive_mode: str = "track") -> dict:
    qubit = C.make_sqc_qubit()
    target = qubit.frequency + 2 * np.pi * target_mhz * 1e-3
    config = make_config()
    config.max_solver_calls = max_solver_calls
    executor = make_executor(
        strategy, qubit, config, target, seed=seed, shots=shots,
        assignment_error=assignment_error, drift_mhz=drift_mhz,
        jump_mhz=jump_mhz, drive_mode=drive_mode,
    )
    result = FrequencyCalibrationRuntime(
        qubit, target, config, executor, run_id=run_id,
    ).run()
    result.update(
        strategy=strategy, target_mhz=target_mhz, seed=seed,
        shots_by_role=(finite_shots_by_role(strategy, shots)
                       if seed is not None else None),
        assignment_error=assignment_error, drift_mhz_per_command=drift_mhz,
        jump_mhz=jump_mhz, drive_mode=drive_mode, protocol_config=asdict(config),
    )
    return result


def save_run(result: dict) -> Path:
    directory = DATA / "runs" / result["run_id"]
    directory.mkdir(parents=True, exist_ok=True)
    write_json(directory / "result.json", {
        key: value for key, value in result.items()
        if key not in {"journal", "transition_log", "cost_ledger"}
    })
    for name, rows in (("state-history.jsonl", result["transition_log"]),
                       ("measurements.jsonl", result["journal"]),
                       ("cost-ledger.jsonl", result["cost_ledger"])):
        with (directory / name).open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(row, default=str) + "\n")
    return directory


def run_checkpoint_resume(target_mhz: float = -10.0) -> dict:
    """Compare a command-boundary restart with an uninterrupted reference."""
    strategy = "short_pulse"
    reference = run_case(
        strategy, target_mhz, "checkpoint__reference__m10.0mhz__det",
    )
    qubit = C.make_sqc_qubit()
    target = qubit.frequency + 2 * np.pi * target_mhz * 1e-3
    config = make_config()
    inner = make_executor(strategy, qubit, config, target)
    interrupted = FrequencyCalibrationRuntime(
        qubit, target, config, SimulatedCrashExecutor(inner, completed_calls=4),
        run_id="checkpoint__resumed__m10.0mhz__det",
    )
    try:
        interrupted.run()
    except RuntimeError as exc:
        if str(exc) != "intentional checkpoint boundary":
            raise
    checkpoint_dir = DATA / "checkpoints" / "short_pulse_command_4"
    interrupted.save_run(str(checkpoint_dir))
    resume_qubit = C.make_sqc_qubit()
    resume_executor = make_executor(strategy, resume_qubit, config, target)
    resumed_runtime = FrequencyCalibrationRuntime.load_run(
        str(checkpoint_dir), qubit=resume_qubit, executor=resume_executor,
    )
    resumed = resumed_runtime.run()
    comparison = {
        "schema_version": 1,
        "crash_boundary_completed_commands": 4,
        "reference_completed_verified_lock": reference["completed_verified_lock"],
        "resumed_completed_verified_lock": resumed["completed_verified_lock"],
        "frequency_abs_difference": abs(reference["f_final"] - resumed["f_final"]),
        "bias_abs_difference": abs(reference["candidate_bias"] - resumed["candidate_bias"]),
        "reference_solver_calls": reference["cost_summary"]["solver_calls"],
        "resumed_solver_calls": resumed["cost_summary"]["solver_calls"],
    }
    comparison["equivalent"] = bool(
        comparison["reference_completed_verified_lock"]
        and comparison["resumed_completed_verified_lock"]
        and comparison["frequency_abs_difference"] <= 1e-9
        and comparison["bias_abs_difference"] <= 1e-9
        and comparison["reference_solver_calls"] == comparison["resumed_solver_calls"]
    )
    write_json(DATA / "checkpoint-resume.json", comparison)
    return comparison


def git_metadata() -> dict:
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                              text=True, capture_output=True, check=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=REPO,
                                text=True, capture_output=True, check=True).stdout.strip())
    return {"revision": revision, "dirty": dirty}


def case_records(smoke: bool) -> list[dict]:
    """Return a fully explicit design; smoke is a structural subset only."""
    records = []
    targets = (-10.0,) if smoke else (-12.0, -10.0, -8.0)
    for target in targets:
        for strategy in ("short_pulse", "ramsey_baseline"):
            records.append(dict(kind="deterministic", strategy=strategy,
                                target_mhz=target, seed=None))
    if smoke:
        return records
    for seed in range(12):
        for strategy in ("short_pulse", "ramsey_baseline"):
            records.append(dict(
                kind="finite_shot", strategy=strategy, target_mhz=-10.0,
                seed=20260813 + seed, shots=6, assignment_error=0.01,
            ))
    for seed in range(8):
        for strategy in ("short_pulse", "ramsey_baseline"):
            records.append(dict(
                kind="process", strategy=strategy, target_mhz=-10.0,
                seed=20260900 + seed, shots=6, assignment_error=0.01,
                drift_mhz=0.02, jump_mhz=0.8,
            ))
    for strategy in ("short_pulse", "ramsey_baseline"):
        records.append(dict(
            kind="budget", strategy=strategy, target_mhz=-10.0, seed=None,
            max_solver_calls=500,
        ))
    for drive_mode in ("track", "fixed"):
        records.append(dict(
            kind="drive_ablation", strategy="short_pulse", target_mhz=-10.0,
            seed=None, drive_mode=drive_mode,
        ))
    return records


def run_id_for(case: dict) -> str:
    target = f"{case['target_mhz']:+.1f}".replace("+", "p").replace("-", "m")
    seed = "det" if case.get("seed") is None else f"s{case['seed']}"
    drive_mode = case.get("drive_mode", "track")
    return f"{case['kind']}__{case['strategy']}__{drive_mode}__{target}mhz__{seed}"


def aggregate(rows: list[dict]) -> list[dict]:
    groups = {}
    for row in rows:
        groups.setdefault(
            (row["kind"], row["strategy"], row.get("drive_mode", "track")), []
        ).append(row)
    output = []
    for (kind, strategy, drive_mode), values in sorted(groups.items()):
        costs = np.asarray([v["solver_calls"] for v in values], dtype=float)
        locked = np.asarray([v["completed_verified_lock"] for v in values], dtype=float)
        output.append({
            "kind": kind, "strategy": strategy, "drive_mode": drive_mode,
            "n": len(values),
            "lock_rate": float(locked.mean()),
            "solver_calls_mean": float(costs.mean()),
            "solver_calls_std": float(costs.std(ddof=1)) if len(costs) > 1 else 0.0,
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--keep-existing", action="store_true")
    args = parser.parse_args()
    if DATA.exists() and not args.keep_existing:
        shutil.rmtree(DATA)
    DATA.mkdir(parents=True, exist_ok=True)
    upstream = json.loads((UPSTREAM / "config.json").read_text(encoding="utf-8"))
    write_json(DATA / "protocol.json", {
        "schema_version": 1, "scope": "simulation-only; not experimental evidence",
        "upstream_config_sha256": sha256(UPSTREAM / "config.json"),
        "upstream_manifest_sha256": sha256(UPSTREAM / "manifest.json"),
        "upstream": upstream, "completion_event": "verify_to_lock",
        "analytic_truth_role": "withheld scoring only", "git": git_metadata(),
        "finite_shot_design": {
            "shots_per_circuit_by_role": {
                "acquire": 1024,
                "short_pulse_track": 1024,
                "ramsey_track": 6,
                "verify": 6,
                "monitor": 6,
            },
            "verify_derivation": "floor(10000 shots / (800 circuits * 2 verifies))",
            "assignment_p01": 0.01,
            "assignment_p10": 0.01,
        },
    })
    records = []
    for case in case_records(args.smoke):
        run_id = run_id_for(case)
        result = run_case(
            case["strategy"], case["target_mhz"], run_id,
            **{key: value for key, value in case.items()
               if key not in {"kind", "strategy", "target_mhz"}},
        )
        directory = save_run(result)
        records.append({
            "run_id": run_id, "kind": case["kind"],
            "strategy": case["strategy"], "target_mhz": case["target_mhz"],
            "seed": case.get("seed"), "drive_mode": case.get("drive_mode", "track"),
            "completed_verified_lock": result["completed_verified_lock"],
            "completion_event": result["completion_event"],
            "solver_calls": result["cost_summary"]["solver_calls"],
            "physical_shots": result["cost_summary"]["shots"],
            "directory": directory.relative_to(REPO).as_posix(),
        })
    with (DATA / "run-index.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    summary = aggregate(records)
    with (DATA / "aggregate.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)
    write_json(DATA / "aggregate.json", summary)
    np.savez_compressed(
        DATA / "aggregate.npz",
        solver_calls=np.asarray([row["solver_calls"] for row in records]),
        completed_verified_lock=np.asarray(
            [row["completed_verified_lock"] for row in records], dtype=bool),
    )
    if not args.smoke:
        checkpoint = run_checkpoint_resume()
        if not checkpoint["equivalent"]:
            raise RuntimeError("checkpoint/resume differs from continuous reference")
    # A checksum manifest cannot contain its own final checksum.  Exclude it
    # explicitly so reruns never retain the previous manifest as an input.
    files = [
        path for path in DATA.rglob("*")
        if path.is_file() and path != DATA / "manifest.json"
    ]
    write_json(DATA / "manifest.json", {
        "schema_version": 1,
        "manifest_excludes_self": True,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generation_command": " ".join(sys.argv),
        "files": {path.relative_to(DATA).as_posix(): sha256(path) for path in files},
    })


if __name__ == "__main__":
    main()
