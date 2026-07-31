"""sqc.workflows.frequency_calibration — FrequencyCalibrationWorkflow.

Multi-stage closed-loop frequency tuning workflow.

Calibration is described as an ordered list of :class:`CalibrationStage`.
Each stage runs one :class:`SinglePointFrequencyCalibration` pass; the workflow
seeds every stage from the previous stage's optimum, forwards each stage's own
measurement axes / stepper config, and merges the per-stage iteration histories
into a single trace with cumulative measurement cost.

The classic **transient→Ramsey hybrid** — a cheap transient+gradient coarse
phase that drives the frequency close, then a robust Ramsey fine phase that
refines to tight tolerance — is the default 2-stage preset, built automatically
from the flat constructor parameters when ``stages`` is not given. So existing
callers keep working unchanged while power users can supply an arbitrary
pipeline via ``stages=[...]``.

A stage advances to the next once its own ``epsilon_f`` residual threshold is
met (or ``max_iter`` is hit). For gradient stages, an optional ``advance_when``
predicate enables adaptive switching on runtime signals a plain threshold
cannot express (gradient flattened, stall, divergence, cost budget); it maps to
the inner solver's ``stop_predicate`` and is only valid for gradient stages.

The single-strategy closed-loop class cannot switch measurement mode mid-loop
(``measure_method`` is fixed at construction), so the pipeline lives here at the
orchestration (workflow) layer.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Literal, Optional

import numpy as np

from sqc.config import CONFIG
from sqc.calibration.frequency import SinglePointFrequencyCalibration
from sqc.workflows.base import Workflow

# Approximate mesolve calls per single frequency measurement, used only to
# annotate the merged history with a *relative* cost trace (not a hard count).
# transient: 2 orthogonal-readout solves (+ cached G3 fit on first order-3 call).
# ramsey double-sweep: 2 sweeps x len(tau_list) solves.
_COST_TRANSIENT = 2


@dataclass
class CalibrationStage:
    """One stage of a frequency-calibration pipeline.

    A stage is a single :class:`SinglePointFrequencyCalibration` pass with its
    own measurement method, stepper, tolerance, and (optionally) measurement
    axes. Fields left as ``None`` inherit the workflow-level default.

    Parameters
    ----------
    name : str
        Label for this stage; becomes the ``phase`` tag on its history rows.
    measure_method : {"ramsey", "transient"}
        Per-iteration frequency measurement strategy.
    step_method : {"secant", "bisection", "gradient"}
        Root-finding step. Only "gradient" supports ``advance_when``.
    epsilon_f : float
        Residual tolerance (angular, rad·GHz). Reaching it advances the
        pipeline to the next stage; the final stage's value is the overall
        convergence tolerance.
    max_iter : int
        Iteration cap for this stage.
    order, g3_source
        Transient kernel order / G₃ source (used when measure_method="transient").
    bracket_tightening, damping, first_bias_step, max_bias_step
        Stepper tuning forwarded to the inner calibration.
    V_a, V_b : float or None
        Optional flux bounds for this stage (clamp for gradient, bracket for
        secant/bisection). Fall back to the workflow-level bounds if None.
    tau_list, t_rabi, t_global : np.ndarray or None
        Per-stage measurement axes. None inherits the workflow default —
        letting a coarse transient stage use a short t_global while a fine
        Ramsey stage uses a long one.
    cost_per_meas : int or None
        Override for the relative cost annotation. None auto-derives it
        (transient→2, ramsey→2·len(tau_list)).
    advance_when : callable or None
        Optional early-advance predicate (gradient stages only). Called each
        iteration with a state dict; return True to advance early. Raises if
        set on a non-gradient stage.
    """

    name: str = "stage"
    measure_method: Literal["ramsey", "transient"] = "ramsey"
    step_method: Literal["secant", "bisection", "gradient"] = "gradient"
    epsilon_f: float = 1e-4
    max_iter: int = 15

    order: int = 1
    g3_source: Literal["fit", "kernel_full"] = "fit"

    bracket_tightening: bool = True
    damping: float = 0.8
    first_bias_step: float = 0.01
    max_bias_step: float = 0.02   # keep coarse search local; see SinglePoint...

    V_a: Optional[float] = None
    V_b: Optional[float] = None

    tau_list: Optional[np.ndarray] = None
    t_rabi: Optional[np.ndarray] = None
    t_global: Optional[np.ndarray] = None

    cost_per_meas: Optional[int] = None
    advance_when: Optional[Callable[[dict], bool]] = None

    # Drive-frequency policy for this stage (gradient step only). See
    # SinglePointFrequencyCalibration.drive_policy. "sweet" (default) locks the
    # drive at the sweet spot; "target" at f_target (LOCKED); "track" predicts
    # f_q to keep the measurement linear (TRACKING); or a callable(state)->ω_d.
    drive_policy: "str | Callable[[dict], float]" = "sweet"
    sensitivity_source: Literal["secant", "model"] = "secant"
    linear_range: Optional[float] = None
    rho: float = 0.6
    converge_streak: int = 1
    # If True, seed this stage's "track" drive from the previous stage's final
    # frequency estimate (f_{d,0} = f̂_q). Essential for a track stage handed a
    # far target — otherwise iter 0 falls back to the sweet spot. Ignored unless
    # drive_policy="track". Default True (harmless for other policies).
    seed_drive_from_prev: bool = True

    def __post_init__(self):
        if self.advance_when is not None and self.step_method != "gradient":
            raise ValueError(
                f"advance_when is only supported for gradient stages; "
                f"stage {self.name!r} uses step_method={self.step_method!r}."
            )
        non_sweet = not (isinstance(self.drive_policy, str)
                         and self.drive_policy == "sweet")
        if non_sweet and self.step_method != "gradient":
            raise ValueError(
                f"drive_policy={self.drive_policy!r} requires "
                f"step_method='gradient'; stage {self.name!r} uses "
                f"step_method={self.step_method!r}."
            )


@dataclass
class FrequencyCalibrationWorkflow(Workflow):
    """Multi-stage closed-loop frequency calibration.

    Runs an ordered list of :class:`CalibrationStage`, seeding each stage from
    the previous stage's optimum. If ``stages`` is None, a 2-stage
    transient→Ramsey hybrid is built from the flat parameters below (backward
    compatible with the original hybrid workflow).

    Parameters
    ----------
    qubit : TransmonQubit
        The qubit whose frequency is tuned.
    f_target : float
        Target qubit frequency (angular, rad·GHz). Required.
    stages : list[CalibrationStage] or None
        Explicit pipeline. When None (default), the flat parameters below build
        the classic 2-stage hybrid.
    V_seed : float or None
        Starting flux bias (Φ₀) for the first stage. If None, falls back to the
        midpoint of ``V_a``/``V_b`` when both are given.
    V_a, V_b : float or None
        Optional flux bounds used to clamp/bracket the search (workflow-level
        default; a stage may override with its own V_a/V_b).

    Hybrid-preset parameters (used only when ``stages`` is None)
    -----------------------------------------------------------
    switch_residual : float
        Residual threshold at which the coarse (transient) phase hands off to
        the fine (Ramsey) phase. Default 2π·5e-3 (5 MHz).
    epsilon_f : float
        Final convergence tolerance for the fine phase. Default 1e-4.
    coarse_max_iter, fine_max_iter : int
        Iteration caps. Defaults 15 and 10.
    fine_step_method : str
        Root-finding step for the fine phase. Default "gradient".
    transient_order : int
        Kernel order for the coarse transient measurement. Default 3.
    damping, first_bias_step, max_bias_step : float
        Gradient stepper tuning for the coarse phase.

    tau_list, t_rabi, t_global
        Workflow-level default measurement axes; a stage inherits these unless
        it sets its own.
    """

    qubit: object
    f_target: Optional[float] = None

    stages: Optional[list] = None

    V_seed: Optional[float] = None
    V_a: Optional[float] = None
    V_b: Optional[float] = None

    # Hybrid-preset parameters (used only when stages is None)
    switch_residual: float = 2.0 * np.pi * 5e-3  # 5 MHz
    epsilon_f: float = 1e-4
    coarse_max_iter: int = 15
    fine_max_iter: int = 10
    fine_step_method: Literal["secant", "bisection", "gradient"] = "gradient"
    transient_order: int = 3

    damping: float = 0.8
    first_bias_step: float = 0.01
    max_bias_step: float = 0.02   # keep coarse search local; see SinglePoint...

    tau_list: Optional[np.ndarray] = None
    t_rabi: np.ndarray = field(default_factory=lambda: CONFIG.pulse.t_rabi.copy())
    t_global: Optional[np.ndarray] = None

    def _seed(self) -> float:
        if self.V_seed is not None:
            return float(self.V_seed)
        if self.V_a is not None and self.V_b is not None:
            return (float(self.V_a) + float(self.V_b)) / 2.0
        raise ValueError("V_seed is required (or provide both V_a and V_b).")

    def _default_stages(self) -> list:
        """Build the classic 2-stage transient→Ramsey hybrid from flat params."""
        return [
            CalibrationStage(
                name="coarse",
                measure_method="transient", step_method="gradient",
                epsilon_f=self.switch_residual, max_iter=self.coarse_max_iter,
                order=self.transient_order,
                damping=self.damping, first_bias_step=self.first_bias_step,
                max_bias_step=self.max_bias_step,
            ),
            CalibrationStage(
                name="fine",
                measure_method="ramsey", step_method=self.fine_step_method,
                epsilon_f=self.epsilon_f, max_iter=self.fine_max_iter,
            ),
        ]

    def run(self) -> dict:
        """Execute the multi-stage calibration.

        Returns
        -------
        dict
            ``V_final``, ``residual``, ``converged``, ``n_iter``,
            ``history`` (merged rows tagged with ``phase``/``global_iter``/
            ``cost``), ``stages`` (per-stage CalibrationTables),
            ``stage_boundaries`` (global iter at each stage's end), and
            ``metrics``. For the default 2-stage hybrid, the legacy keys
            ``switch_iter``, ``coarse``, ``fine`` and the matching ``metrics``
            entries are also provided.
        """
        if self.f_target is None:
            raise ValueError("f_target is required.")

        stages = self.stages if self.stages is not None else self._default_stages()
        if not stages:
            raise ValueError("stages must be a non-empty list.")

        tables = []
        seed = self._seed()
        f_hat = None                       # previous stage's frequency estimate
        for st in stages:
            # Seed a "track" stage's initial drive from the previous stage's
            # estimate (f_{d,0} = f̂_q), so iter 0 doesn't fall back to the
            # sweet spot far from a distant target.
            omega_d_seed = (
                f_hat if (st.seed_drive_from_prev and st.drive_policy == "track")
                else None
            )
            cal = SinglePointFrequencyCalibration(
                qubit=self.qubit, f_target=self.f_target,
                V_seed=seed,
                V_a=st.V_a if st.V_a is not None else self.V_a,
                V_b=st.V_b if st.V_b is not None else self.V_b,
                epsilon_f=st.epsilon_f, max_iter=st.max_iter,
                measure_method=st.measure_method, step_method=st.step_method,
                bracket_tightening=st.bracket_tightening,
                order=st.order, g3_source=st.g3_source,
                damping=st.damping, first_bias_step=st.first_bias_step,
                max_bias_step=st.max_bias_step,
                stop_predicate=st.advance_when,
                drive_policy=st.drive_policy,
                sensitivity_source=st.sensitivity_source,
                linear_range=st.linear_range, rho=st.rho,
                converge_streak=st.converge_streak,
                omega_d_seed=omega_d_seed,
                tau_list=st.tau_list if st.tau_list is not None else self.tau_list,
                t_rabi=st.t_rabi if st.t_rabi is not None else self.t_rabi,
                t_global=st.t_global if st.t_global is not None else self.t_global,
            ).calibrate()
            tables.append(cal)
            seed = float(cal.fit_params["V_opt"])  # hand off to next stage
            f_hat = float(cal.fit_params["residual"]) + float(self.f_target)

        return self._merge(stages, tables)

    def _stage_cost(self, st) -> int:
        """Relative mesolve cost per measurement for one stage."""
        if st.cost_per_meas is not None:
            return int(st.cost_per_meas)
        if st.measure_method == "transient":
            return _COST_TRANSIENT
        tau = st.tau_list if st.tau_list is not None else self.tau_list
        tau_n = len(tau) if tau is not None else len(CONFIG.pulse.make_time(0, 200))
        return 2 * tau_n  # ramsey double-sweep

    def _merge(self, stages, tables) -> dict:
        """Stitch per-stage histories into one trace with cumulative cost."""
        merged: list[dict] = []
        boundaries: list[int] = []
        stage_iters: list[int] = []
        stage_costs: list[int] = []
        cost = 0
        gi = 0
        for st, table in zip(stages, tables):
            per = self._stage_cost(st)
            hist = table.fit_params["history"]
            for h in hist:
                cost += per
                gi += 1
                merged.append({**h, "phase": st.name, "global_iter": gi,
                               "cost": cost})
            boundaries.append(gi)
            stage_iters.append(len(hist))
            stage_costs.append(per)

        last = tables[-1].fit_params
        V_final = float(last["V_opt"])
        residual = float(last["residual"])
        converged = abs(residual) <= stages[-1].epsilon_f

        result = {
            "V_final": V_final,
            "residual": residual,
            "converged": converged,
            "n_iter": gi,
            "history": merged,
            "stages": tables,
            "stage_boundaries": boundaries,
            "metrics": {
                "total_cost": cost,
                "stage_names": [st.name for st in stages],
                "stage_iters": stage_iters,
                "stage_costs": stage_costs,
            },
        }

        # Legacy 2-stage hybrid aliases (kept for backward compatibility).
        if self.stages is None and len(stages) == 2:
            result["switch_iter"] = boundaries[0]
            result["coarse"] = tables[0]
            result["fine"] = tables[1]
            result["metrics"].update({
                "coarse_iter": stage_iters[0],
                "fine_iter": stage_iters[1],
                "cost_transient_per_meas": stage_costs[0],
                "cost_ramsey_per_meas": stage_costs[1],
            })
        return result
