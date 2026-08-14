#!/usr/bin/env python3
"""Generate simulation-only Fig. 3 and Supplemental Figs. S5/S6/S8."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "paper" / "data" / "fig3-simulation"
RUNS = DATA / "runs"
DERIVED = DATA / "derived"
HERE = Path(__file__).resolve().parent
STYLE = ROOT.parent / "scholaraio" / ".claude" / "skills" / "draw" / "nature_pub.mplstyle"
TWO_PI = 2.0 * np.pi
COLORS = {"short_pulse": "#2E8B57", "ramsey_baseline": "#3572B0"}
LABELS = {"short_pulse": "Short pulse", "ramsey_baseline": "Ramsey"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def load_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]


def read_index() -> list[dict]:
    with (DATA / "run-index.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["target_mhz"] = float(row["target_mhz"])
        row["solver_calls"] = int(row["solver_calls"])
        row["physical_shots"] = int(row["physical_shots"])
        row["locked"] = row["completed_verified_lock"] == "True"
    return rows


def verify_source_manifest() -> None:
    manifest = load_json(DATA / "manifest.json")
    mismatches = [
        relative for relative, expected in manifest["files"].items()
        if sha256(DATA / relative) != expected
    ]
    if mismatches:
        raise RuntimeError(f"source manifest mismatch: {mismatches}")


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_npz(path: Path, rows: list[dict], numeric: tuple[str, ...]) -> None:
    payload = {key: np.asarray([row[key] for row in rows]) for key in numeric}
    payload["run_id"] = np.asarray([row["run_id"] for row in rows])
    np.savez_compressed(path, **payload)


def measurements(run_id: str) -> list[dict]:
    return load_jsonl(RUNS / run_id / "measurements.jsonl")


def target_frequency(run_id: str) -> float:
    rows = measurements(run_id)
    verify = next(row for row in rows if row["command"] == "VerifyFrequency")
    return float(verify["commanded_drive"])


def residual_mhz(row: dict) -> float:
    return abs(float(row["residual"])) * 1e3 / TWO_PI


def configure_style() -> None:
    if STYLE.exists():
        plt.style.use(STYLE)
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 7.0,
        "axes.linewidth": 0.65, "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })


def export(fig, stem: str) -> None:
    for suffix, kwargs in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
        fig.savefig(HERE / f"{stem}.{suffix}", **kwargs)
    plt.close(fig)


def fig3(rows: list[dict]) -> None:
    deterministic = [r for r in rows if r["kind"] == "deterministic" and r["target_mhz"] == -10]
    trajectories = []
    for item in deterministic:
        science = [x for x in measurements(item["run_id"]) if x.get("residual") is not None]
        for index, event in enumerate(science):
            trajectories.append({
                "run_id": item["run_id"], "strategy": item["strategy"],
                "step": index, "command": event["command"],
                "residual_mhz": residual_mhz(event),
                "state_after": event["state_after"],
            })
    write_csv(DERIVED / "fig3-trajectory.csv", trajectories)
    save_npz(DERIVED / "fig3-trajectory.npz", trajectories, ("step", "residual_mhz"))

    finite = [r for r in rows if r["kind"] == "finite_shot"]
    process = [r for r in rows if r["kind"] == "process"]
    summary = []
    for kind, values in (("finite_shot", finite), ("process", process)):
        for strategy in ("short_pulse", "ramsey_baseline"):
            group = [r for r in values if r["strategy"] == strategy]
            summary.append({
                "run_id": f"{kind}__{strategy}", "kind": kind,
                "strategy": strategy, "n": len(group),
                "lock_rate": np.mean([r["locked"] for r in group]),
                "solver_calls_mean": np.mean([r["solver_calls"] for r in group]),
                "solver_calls_std": np.std([r["solver_calls"] for r in group], ddof=1),
            })
    write_csv(DERIVED / "fig3-summary.csv", summary)
    save_npz(
        DERIVED / "fig3-summary.npz", summary,
        ("n", "lock_rate", "solver_calls_mean", "solver_calls_std"),
    )

    fig, axes = plt.subplots(1, 3, figsize=(7.10, 2.35))
    ax = axes[0]
    markers = {"AcquireFrequency": "s", "TrackFrequency": "o",
               "VerifyFrequency": "D", "MonitorFrequency": "x"}
    for strategy in ("short_pulse", "ramsey_baseline"):
        points = [r for r in trajectories if r["strategy"] == strategy]
        ax.plot([r["step"] for r in points], [r["residual_mhz"] for r in points],
                color=COLORS[strategy], lw=1.0, label=LABELS[strategy])
        for command, marker in markers.items():
            subset = [r for r in points if r["command"] == command]
            ax.scatter([r["step"] for r in subset], [r["residual_mhz"] for r in subset],
                       s=18, marker=marker, color=COLORS[strategy], zorder=3)
    ax.axhline(0.5, color="#777777", lw=0.7, ls="--", label="Verify threshold")
    ax.set_yscale("log")
    ax.set_xlabel("Science-command index")
    ax.set_ylabel("Absolute residual (MHz)")
    ax.legend(frameon=False, fontsize=5.8, loc="upper right")
    ax.set_title("Deterministic convergence", fontsize=7.2)

    ax = axes[1]
    finite_summary = [r for r in summary if r["kind"] == "finite_shot"]
    for index, item in enumerate(finite_summary):
        ax.bar(index, item["solver_calls_mean"], yerr=item["solver_calls_std"],
               width=0.58, color=COLORS[item["strategy"]], capsize=2)
        ax.text(index, item["solver_calls_mean"] * 1.03,
                f"{int(item['lock_rate'] * item['n'])}/{item['n']}",
                ha="center", va="bottom", fontsize=6)
    ax.set_xticks([0, 1], [LABELS[r["strategy"]] for r in finite_summary])
    ax.set_ylabel("Solver calls (mean +/- SD)")
    ax.set_title("Finite-shot realizations", fontsize=7.2)

    ax = axes[2]
    process_summary = [r for r in summary if r["kind"] == "process"]
    for index, item in enumerate(process_summary):
        successes = int(item["lock_rate"] * item["n"])
        ax.bar(index, item["lock_rate"], width=0.58, color=COLORS[item["strategy"]])
        ax.text(index, item["lock_rate"] + 0.025, f"{successes}/{item['n']}",
                ha="center", va="bottom", fontsize=6)
    ax.set_xticks([0, 1], [LABELS[r["strategy"]] for r in process_summary])
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Verified-lock fraction")
    ax.set_title("Drift + controlled jump", fontsize=7.2)
    for label, axis in zip("abc", axes):
        axis.text(-0.17, 1.04, label, transform=axis.transAxes,
                  fontsize=8.5, fontweight="bold")
        axis.grid(axis="y", color="#DDDDDD", lw=0.4)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.20, top=0.86, wspace=0.42)
    export(fig, "fig3-simulation-evidence")


def fig_s5(rows: list[dict]) -> None:
    ablation = [r for r in rows if r["kind"] == "drive_ablation"]
    points = []
    for item in ablation:
        track = [x for x in measurements(item["run_id"]) if x["command"] == "TrackFrequency"]
        result = load_json(RUNS / item["run_id"] / "result.json")
        target = (
            target_frequency(item["run_id"])
            if result["completed_verified_lock"]
            else float(result["f_final"]) - float(track[-1]["residual"])
        )
        for index, event in enumerate(track):
            points.append({
                "run_id": item["run_id"], "drive_mode": item["drive_mode"],
                "step": index, "commanded_drive_mhz":
                    (float(event["commanded_drive"]) - target) * 1e3 / TWO_PI,
                "applied_drive_mhz":
                    (float(event["applied_drive"]) - target) * 1e3 / TWO_PI,
                "residual_mhz": residual_mhz(event),
                "state_after": event["state_after"],
            })
    write_csv(DERIVED / "figS5-drive-tracking.csv", points)
    save_npz(DERIVED / "figS5-drive-tracking.npz", points,
             ("step", "commanded_drive_mhz", "applied_drive_mhz", "residual_mhz"))
    fig, axes = plt.subplots(1, 2, figsize=(5.0, 2.35))
    for mode, color, marker in (("track", "#2E8B57", "o"), ("fixed", "#B44B4B", "s")):
        subset = [r for r in points if r["drive_mode"] == mode]
        label = "Drive tracked" if mode == "track" else "Fixed drive"
        axes[0].plot([r["step"] for r in subset], [r["applied_drive_mhz"] for r in subset],
                     marker=marker, color=color, lw=1.0, ms=4, label=label)
        axes[1].plot([r["step"] for r in subset], [r["residual_mhz"] for r in subset],
                     marker=marker, color=color, lw=1.0, ms=4, label=label)
    axes[0].set_xlabel("Track index")
    axes[0].set_ylabel("Applied drive offset (MHz)")
    axes[1].set_xlabel("Track index")
    axes[1].set_ylabel("Absolute residual (MHz)")
    axes[1].set_yscale("log")
    axes[1].axhline(0.5, color="#777777", lw=0.7, ls="--")
    axes[0].legend(frameon=False, fontsize=6)
    for label, axis in zip("ab", axes):
        axis.text(-0.17, 1.04, label, transform=axis.transAxes,
                  fontsize=8.5, fontweight="bold")
        axis.grid(axis="y", color="#DDDDDD", lw=0.4)
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.20, top=0.90, wspace=0.38)
    export(fig, "figS5-drive-tracking-simulation")


def fig_s6(rows: list[dict]) -> None:
    selected = [r for r in rows if r["kind"] in {"deterministic", "finite_shot", "process"}]
    write_csv(DERIVED / "figS6-solver-cost.csv", selected)
    save_npz(DERIVED / "figS6-solver-cost.npz", selected,
             ("solver_calls", "physical_shots", "locked"))
    fig, ax = plt.subplots(figsize=(3.5, 2.55))
    kinds = ("deterministic", "finite_shot", "process")
    x = np.arange(len(kinds))
    width = 0.34
    for offset, strategy in ((-width / 2, "short_pulse"), (width / 2, "ramsey_baseline")):
        means, errors = [], []
        for kind in kinds:
            values = [r["solver_calls"] for r in selected
                      if r["kind"] == kind and r["strategy"] == strategy]
            means.append(np.mean(values))
            errors.append(np.std(values, ddof=1))
        ax.bar(x + offset, means, width, yerr=errors, capsize=2,
               color=COLORS[strategy], label=LABELS[strategy])
    ax.set_xticks(x, ["Deterministic", "Finite shot", "Drift + jump"])
    ax.set_ylabel("Instrumented solver calls")
    ax.legend(frameon=False, fontsize=6)
    ax.grid(axis="y", color="#DDDDDD", lw=0.4)
    fig.text(0.98, 0.035, "mesolve + sesolve; not acquisition time",
             ha="right", va="bottom", fontsize=5.3, color="#555555")
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.23, top=0.94)
    export(fig, "figS6-solver-cost-simulation")


def fig_s8(rows: list[dict]) -> None:
    process = [r for r in rows if r["kind"] == "process"]
    states = ("acquire", "track", "verify", "lock", "safe_stop")
    state_y = {state: index for index, state in enumerate(states)}
    state_rows = []
    for item in process:
        for index, event in enumerate(measurements(item["run_id"])):
            state_rows.append({
                "run_id": item["run_id"], "strategy": item["strategy"],
                "seed": item["seed"], "step": index,
                "state_before": event["state_before"], "state_after": event["state_after"],
                "locked": item["locked"],
            })
    write_csv(DERIVED / "figS8-fsm-realizations.csv", state_rows)
    save_npz(DERIVED / "figS8-fsm-realizations.npz", state_rows, ("step", "locked"))
    fig, axes = plt.subplots(2, 1, figsize=(5.0, 3.65), sharex=False)
    for axis, strategy in zip(axes, ("short_pulse", "ramsey_baseline")):
        group = [r for r in process if r["strategy"] == strategy]
        for lane, item in enumerate(group):
            seq = [r for r in state_rows if r["run_id"] == item["run_id"]]
            xs = [r["step"] for r in seq]
            ys = [state_y[r["state_after"]] + lane * 0.035 for r in seq]
            axis.plot(xs, ys, color=COLORS[strategy], lw=0.7,
                      alpha=0.75 if item["locked"] else 0.35)
            axis.scatter(xs[-1], ys[-1], s=16,
                         marker="o" if item["locked"] else "x",
                         color=COLORS[strategy])
            lock = next((r for r in seq if r["state_after"] == "lock"), None)
            if lock is not None:
                axis.scatter(lock["step"], state_y["lock"] + lane * 0.035,
                             s=28, marker="*", facecolor="white",
                             edgecolor=COLORS[strategy], linewidth=0.7, zorder=4)
        axis.set_yticks(range(len(states)), [s.replace("_", " ").title() for s in states])
        axis.set_ylabel(LABELS[strategy])
        axis.grid(color="#E2E2E2", lw=0.35)
    axes[-1].set_xlabel("Command index")
    axes[0].text(0.99, 0.03, "star: first verified lock; circle/cross: final success/failure",
                 transform=axes[0].transAxes, ha="right", va="bottom",
                 fontsize=5.2, color="#555555")
    axes[0].text(-0.12, 1.04, "a", transform=axes[0].transAxes,
                 fontsize=8.5, fontweight="bold")
    axes[1].text(-0.12, 1.04, "b", transform=axes[1].transAxes,
                 fontsize=8.5, fontweight="bold")
    fig.subplots_adjust(left=0.20, right=0.98, bottom=0.13, top=0.96, hspace=0.35)
    export(fig, "figS8-fsm-realizations-simulation")


def update_manifest() -> None:
    manifest = load_json(DATA / "manifest.json")
    manifest["files"] = {
        path.relative_to(DATA).as_posix(): sha256(path)
        for path in DATA.rglob("*")
        if path.is_file() and path != DATA / "manifest.json"
    }
    manifest["derived_figure_files"] = {
        path.relative_to(ROOT / "paper").as_posix(): sha256(path)
        for path in HERE.iterdir()
        if path.is_file() and path.suffix in {".pdf", ".png", ".svg"}
    }
    (DATA / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8",
    )


def main() -> None:
    verify_source_manifest()
    DERIVED.mkdir(parents=True, exist_ok=True)
    configure_style()
    rows = read_index()
    fig3(rows)
    fig_s5(rows)
    fig_s6(rows)
    fig_s8(rows)
    update_manifest()
    print("generated Fig. 3 and Figs. S5/S6/S8", flush=True)


if __name__ == "__main__":
    main()
