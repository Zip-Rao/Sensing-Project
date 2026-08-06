"""sqc.workflows.sensing — SensingWorkflow: unified research entry point.

The ``SensingWorkflow`` is the **single user-facing API** for the quantum
sensing simulation platform.  It combines device configuration, protocol
selection, measurement execution, waveform reconstruction, parameter
sweeping, and A/B comparison in one class.

Usage::

    from sqc.workflows import SensingWorkflow

    wf = SensingWorkflow()
    wf.configure(protocol="transient", signal_amplitude=0.01,
                 reconstruction="wiener", t_rabi_duration=20)
    result = wf.run()
    wf.plot()

    # parameter sweep
    sweep = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02])

    # method comparison
    cmp_res = wf.compare(methods=["wiener", "hammerstein", "lm"])

Protocol → class mapping (handbook 6d.5)
----------------------------------------
===============  ============================  ==========================
``protocol=``     Experiment class              Reconstruction class
===============  ============================  ==========================
``"transient"``   TransientSensingExperiment    TransientReconstruction
``"ramsey"``      RamseyExperiment              RamseyReconstruction
``"echo"``        DiffEchoExperiment            EchoReconstruction
``"cryoscope"``   CryoscopeExperiment           CryoscopeReconstruction
``"delay_ramsey"``DelayRamseyExperiment         DelayRamseyReconstruction
===============  ============================  ==========================

See ``idea/refactor/phase_6_handbook.md`` §6d for the full design.
"""
from __future__ import annotations

import time
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from sqc.config import CONFIG, reconfigure as _reconfigure
from sqc.control.flux_signal import FluxSignal
from sqc.simulation.result import ExperimentResult
from sqc.workflows.base import Workflow


# ═══════════════════════════════════════════════════════════════════════
# Result data classes
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class WorkflowResult:
    """Result of a single ``SensingWorkflow.run()`` call.

    Attributes
    ----------
    config_snapshot : dict
        Deep copy of ``_params`` at the moment ``run()`` was called.
    measurement : ExperimentResult or None
        Raw measurement output (p_e, delta_p, kernel, etc.).
    reconstructed_signal : FluxSignal or None
        Reconstructed flux waveform B(t).  None if ``reconstruct=False``.
    reconstruction_details : dict
        Algorithm-specific details (LM iteration history, etc.).
    calibration : CalibrationTable or None
        Calibration table when ``calibrate=True``.
    """
    config_snapshot: dict = field(default_factory=dict)
    measurement: ExperimentResult | None = None
    reconstructed_signal: FluxSignal | None = None
    reconstruction_details: dict = field(default_factory=dict)
    calibration: object | None = None


@dataclass
class SweepResult:
    """Result of a parameter sweep from ``sweep()``.

    Attributes
    ----------
    param_path : str
        Dotted path of the swept parameter, e.g. ``"signal.amplitude"``.
    values : list
        The values that were swept.
    results : list[WorkflowResult]
        One result per swept value.
    metrics : dict[str, list]
        Key → list of per-value metrics (``"snr"``, ``"rmse"``,
        ``"peak"``, …).
    """
    param_path: str = ""
    values: list = field(default_factory=list)
    results: list[WorkflowResult] = field(default_factory=list)
    metrics: dict[str, list] = field(default_factory=dict)


@dataclass
class CompareResult:
    """Result of reconstruction-method comparison from ``compare()``.

    Attributes
    ----------
    methods : list[str]
        Method names in order.
    signals : dict[str, FluxSignal]
        Reconstructed signal per method.
    metrics : dict[str, dict]
        Per-method metric dict, e.g. ``{"wiener": {"rmse": 0.002, ...}}``.
    best : str
        Method with the lowest RMSE (or empty string if not computable).
    """
    methods: list[str] = field(default_factory=list)
    signals: dict[str, FluxSignal] = field(default_factory=dict)
    metrics: dict[str, dict] = field(default_factory=dict)
    best: str = ""


# ═══════════════════════════════════════════════════════════════════════
# Stub return types (for future phases)
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class DiffReport:
    """(stub) Difference report between two WorkflowResults."""
    signal_rmse: float = 0.0
    snr_diff: float = 0.0
    param_diffs: dict = field(default_factory=dict)


@dataclass
class BenchmarkResult:
    """(stub) Benchmark result table."""
    table: list[dict] = field(default_factory=list)
    best: str = ""


@dataclass
class NoiseReport:
    """(stub) Noise characterization report."""
    mean: float = 0.0
    std: float = 0.0
    psd: np.ndarray | None = None


@dataclass
class CVResult:
    """(stub) Cross-validation result."""
    mean_rmse: float = 0.0
    std_rmse: float = 0.0
    per_fold: list[float] = field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════
# Protocol → (Experiment, Reconstruction) mapping
# ═══════════════════════════════════════════════════════════════════════

_PROTOCOL_MAP: dict[str, tuple[str, str, str]] = {
    "transient": ("TransientSensingExperiment",
                  "TransientReconstruction",
                  "transient"),
    "ramsey": ("RamseyExperiment",
               "RamseyReconstruction",
               "ramsey"),
    "echo": ("DiffEchoExperiment",
             "EchoReconstruction",
             "echo"),
    "cryoscope": ("CryoscopeExperiment",
                  "CryoscopeReconstruction",
                  "cryoscope"),
    "delay_ramsey": ("DelayRamseyExperiment",
                     "DelayRamseyReconstruction",
                     "delay_ramsey"),
}

_PROTOCOLS_REQUIRING_CALIBRATION = {"cryoscope", "delay_ramsey"}

_VALID_PROTOCOLS = set(_PROTOCOL_MAP.keys())


# ═══════════════════════════════════════════════════════════════════════
# SensingWorkflow
# ═══════════════════════════════════════════════════════════════════════

class SensingWorkflow(Workflow):
    """Unified research entry point for the quantum sensing platform.

    One class handles **configuration, protocol execution, parameter
    sweeping, A/B comparison, and visualisation**.  Under the hood it
    delegates to the Experiment, Reconstruction, and Calibration layers.

    Parameters
    ----------
    All parameters are set via ``configure()``.  See its docstring for
    the full list organised by group (Qubit, Protocol & Signal,
    Reconstruction, Pulse, Hardware).

    Examples
    --------
    Minimal transient sensing + Wiener reconstruction::

        wf = SensingWorkflow()
        wf.configure(protocol="transient", signal_amplitude=0.01,
                     reconstruction="wiener", t_rabi_duration=20)
        result = wf.run()
        wf.plot()

    Sweep signal amplitude::

        sweep = wf.sweep("signal.amplitude", [0.005, 0.01, 0.02, 0.05])
        print(sweep.metrics["snr"])

    Compare three reconstruction methods on the same data::

        wf.run(measure=True, reconstruct=False)
        cmp_res = wf.compare(methods=["wiener", "hammerstein", "lm"])
        print(f"Best method: {cmp_res.best}")

    Not yet available (planned — post-v1)
    -------------------------------------
    The following methods are placeholders and raise ``NotImplementedError``
    in v1; they are not part of the public v1 API (decision D4):
    ``pipeline``, ``multi_qubit``, ``crosstalk``, ``save``, ``load``,
    ``diff``, ``benchmark``, ``find_optimal_work_point``,
    ``detectability_limit``, ``noise_characterize``, ``cross_validate``.
    """

    def __init__(self) -> None:
        super().__init__()
        # ── parameter store ──────────────────────────────────────────
        self._params: dict[str, dict] = {
            "qubit": {
                "EC": 0.2, "EJ": 15.0, "T1": 10_000, "T2": 8_000,
                "flux_bias": 0.0, "n_levels": 3,
            },
            "protocol": "transient",
            "signal": {
                "type": 4, "amplitude": 0.01,
                "width": 10, "center": 100,
                "rise": 10, "fall": 10,
                "frequency": 0.01,
                "custom": None,
            },
            "reconstruction": {
                "method": "wiener", "lambda_reg": 5.0,
                "lm_n_basis": 50, "lm_max_iter": 10,
                "lm_basis_type": "fourier", "use_adjoint": True,
            },
            "pulse": {
                "t_rabi_duration": 20.0, "t_global_start": -50.0,
                "t_global_end": 400.0,
                "rotation_angle": np.pi / 2, "rabi_rate": None,
                "envelope": "square", "envelope_sigma": None,
                "phase1": np.pi / 2, "phase2": 0.0,
            },
            "hardware": {
                "sample_rate": 2.0,
            },
        }
        # ── runtime state ────────────────────────────────────────────
        self._cfg: Any = CONFIG
        self._last_result: WorkflowResult | None = None
        self._last_sweep: SweepResult | None = None
        self._last_compare: CompareResult | None = None
        # cache for measure-then-reconstruct
        self._cached_measurement: ExperimentResult | None = None

    # ==================================================================
    # configure()
    # ==================================================================

    def configure(
        self,
        # ── Qubit ────────────────────────────────────────────────────
        EC: float | None = None,
        EJ: float | None = None,
        T1: float | None = None,
        T2: float | None = None,
        flux_bias: float | None = None,
        n_levels: int | None = None,
        # ── Protocol & Signal ────────────────────────────────────────
        protocol: str | None = None,
        signal_type: int | None = None,
        signal_amplitude: float | None = None,
        signal_width: float | None = None,
        signal_center: float | None = None,
        signal_rise: float | None = None,
        signal_fall: float | None = None,
        signal_frequency: float | None = None,
        signal_custom: np.ndarray | None = None,
        # ── Reconstruction ───────────────────────────────────────────
        reconstruction: str | None = None,
        lambda_reg: float | None = None,
        lm_n_basis: int | None = None,
        lm_max_iter: int | None = None,
        lm_basis_type: str | None = None,
        use_adjoint: bool | None = None,
        # ── Pulse ────────────────────────────────────────────────────
        t_rabi_duration: float | None = None,
        t_global_start: float | None = None,
        t_global_end: float | None = None,
        rotation_angle: float | None = None,
        rabi_rate: float | None = None,
        pulse_envelope: object | None = None,
        envelope_sigma: float | None = None,
        phase1: float | None = None,
        phase2: float | None = None,
        # ── Hardware ─────────────────────────────────────────────────
        sample_rate: float | None = None,
        # ── catch-all (ignored) ──────────────────────────────────────
        **__,
    ) -> "SensingWorkflow":
        """Set or update parameters.  Only explicitly passed values are changed;
        all others keep their current value.  Returns ``self`` for chaining.

        Parameters are grouped by functional area:

        **Qubit** — ``EC``, ``EJ``, ``T1``, ``T2``, ``flux_bias``,
        ``n_levels``.  ``EC`` / ``EJ`` are in GHz (not rad*GHz).

        **Protocol & Signal** — ``protocol`` (one of ``"transient"``,
        ``"ramsey"``, ``"echo"``, ``"cryoscope"``, ``"delay_ramsey"``),
        ``signal_type`` (FluxSignal type 0-8), ``signal_amplitude``,
        ``signal_width``, ``signal_center``, ``signal_rise``,
        ``signal_fall``, ``signal_frequency``, ``signal_custom``.

        **Reconstruction** — ``reconstruction`` (method name, e.g.
        ``"wiener"``), ``lambda_reg``, ``lm_n_basis``, ``lm_max_iter``,
        ``lm_basis_type``, ``use_adjoint``.

        **Pulse** — ``t_rabi_duration``, ``t_global_start``,
        ``t_global_end`` (all in ns), ``rotation_angle``, ``rabi_rate``,
        ``pulse_envelope``, ``envelope_sigma``, ``phase1``, ``phase2``.

        **Hardware** — ``sample_rate`` (GSa/s).
        """
        qp = self._params["qubit"]
        if EC is not None:            qp["EC"] = EC
        if EJ is not None:            qp["EJ"] = EJ
        if T1 is not None:            qp["T1"] = T1
        if T2 is not None:            qp["T2"] = T2
        if flux_bias is not None:     qp["flux_bias"] = flux_bias
        if n_levels is not None:      qp["n_levels"] = n_levels

        if protocol is not None:
            if protocol not in _VALID_PROTOCOLS:
                raise ValueError(
                    f"Unknown protocol '{protocol}'. "
                    f"Valid choices: {sorted(_VALID_PROTOCOLS)}"
                )
            self._params["protocol"] = protocol

        sp = self._params["signal"]
        if signal_type is not None:        sp["type"] = signal_type
        if signal_amplitude is not None:   sp["amplitude"] = signal_amplitude
        if signal_width is not None:       sp["width"] = signal_width
        if signal_center is not None:      sp["center"] = signal_center
        if signal_rise is not None:        sp["rise"] = signal_rise
        if signal_fall is not None:        sp["fall"] = signal_fall
        if signal_frequency is not None:   sp["frequency"] = signal_frequency
        if signal_custom is not None:      sp["custom"] = signal_custom

        rp = self._params["reconstruction"]
        if reconstruction is not None:  rp["method"] = reconstruction
        if lambda_reg is not None:      rp["lambda_reg"] = lambda_reg
        if lm_n_basis is not None:      rp["lm_n_basis"] = lm_n_basis
        if lm_max_iter is not None:     rp["lm_max_iter"] = lm_max_iter
        if lm_basis_type is not None:   rp["lm_basis_type"] = lm_basis_type
        if use_adjoint is not None:     rp["use_adjoint"] = use_adjoint

        pp = self._params["pulse"]
        if t_rabi_duration is not None:  pp["t_rabi_duration"] = t_rabi_duration
        if t_global_start is not None:   pp["t_global_start"] = t_global_start
        if t_global_end is not None:     pp["t_global_end"] = t_global_end
        if rotation_angle is not None:
            pp["rotation_angle"] = rotation_angle
            pp["rabi_rate"] = None
        if rabi_rate is not None:
            pp["rabi_rate"] = rabi_rate
            pp["rotation_angle"] = None
        if pulse_envelope is not None:   pp["envelope"] = pulse_envelope
        if envelope_sigma is not None:   pp["envelope_sigma"] = envelope_sigma
        if phase1 is not None:           pp["phase1"] = phase1
        if phase2 is not None:           pp["phase2"] = phase2

        hp = self._params["hardware"]
        if sample_rate is not None:   hp["sample_rate"] = sample_rate

        return self

    # ==================================================================
    # run() — execute a single sensing pipeline
    # ==================================================================

    def run(
        self,
        measure: bool = True,
        reconstruct: bool = True,
        calibrate: bool = False,
    ) -> WorkflowResult:
        """Execute the configured sensing pipeline.

        Parameters
        ----------
        measure : bool
            If True, construct an Experiment, run ``mesolve``, and
            store the raw measurement (p_e, delta_p, kernel, etc.).
        reconstruct : bool
            If True, construct a Reconstruction and run
            ``reconstruct()`` on the measurement data.  Requires
            ``measure=True`` (or a previous measurement cached).
        calibrate : bool
            If True, run calibration **before** measurement.  Only
            supported for ``cryoscope`` and ``delay_ramsey`` protocols.

        Returns
        -------
        WorkflowResult
            Container holding measurement, reconstructed signal,
            calibration table, and a config snapshot.

        Raises
        ------
        ValueError
            If ``reconstruct=True`` but no measurement data is available.
        """
        protocol = self._params["protocol"]

        # ── 0. build config from pulse/hardware params ────────────────
        self._sync_config()

        # ── 1. build qubit ───────────────────────────────────────────
        qubit = self._build_qubit()

        # ── 2. snapshot params ───────────────────────────────────────
        result = WorkflowResult(
            config_snapshot=deepcopy(self._params),
        )

        # ── 3. calibrate (if requested + protocol supports it) ───────
        cal_table = None
        if calibrate and protocol in _PROTOCOLS_REQUIRING_CALIBRATION:
            cal_table = self._run_calibration(qubit, protocol)
            result.calibration = cal_table
        elif calibrate:
            # calibration requested but not supported — store None
            result.calibration = None

        # ── 4. measure ───────────────────────────────────────────────
        measurement: ExperimentResult | None = None
        if measure:
            measurement = self._run_measurement(qubit)
            result.measurement = measurement
            self._cached_measurement = measurement
        elif reconstruct:
            # use cached measurement if available
            measurement = self._cached_measurement

        # ── 5. reconstruct ───────────────────────────────────────────
        if reconstruct:
            if measurement is None:
                raise ValueError(
                    "Cannot reconstruct without measurement data. "
                    "Call run(measure=True) first, or pass a previous "
                    "measurement result."
                )
            rec_signal, rec_details = self._run_reconstruction(
                measurement, qubit, cal_table, protocol,
            )
            result.reconstructed_signal = rec_signal
            result.reconstruction_details = rec_details

        self._last_result = result
        return result

    # ==================================================================
    # sweep() — parameter scanning
    # ==================================================================

    def sweep(
        self,
        param_path: str,
        values: list[float],
    ) -> SweepResult:
        """Sweep a single parameter over a list of values.

        Parameters
        ----------
        param_path : str
            Dotted path ``"group.field"``.  Supported groups:
            ``signal`` (``amplitude``, ``width``, ``center``, ``rise``,
            ``fall``, ``type``), ``reconstruction`` (``lambda_reg``,
            ``lm_n_basis``, ``lm_max_iter``, ``lm_basis_type``),
            ``qubit`` (``EC``, ``EJ``, ``flux_bias``, ``n_levels``,
            ``T1``, ``T2``), ``pulse`` (``t_rabi_duration``,
            ``t_global_start``, ``t_global_end``), ``hardware``
            (``sample_rate``).
        values : list[float]
            Values to sweep over.

        Returns
        -------
        SweepResult
            ``.results`` holds one ``WorkflowResult`` per value;
            ``.metrics`` holds per-value aggregate metrics.

        Raises
        ------
        ValueError
            If ``param_path`` cannot be parsed or the group/field is
            unknown.
        """
        group, field = self._parse_param_path(param_path)
        results: list[WorkflowResult] = []
        metrics: dict[str, list] = {
            "snr": [], "rmse": [], "peak": [],
        }

        for val in values:
            # update parameter
            self._params[group][field] = val
            if group == "pulse" and field == "rotation_angle":
                self._params["pulse"]["rabi_rate"] = None
            elif group == "pulse" and field == "rabi_rate":
                self._params["pulse"]["rotation_angle"] = None
            # re-sync config (needed if pulse/hardware group)
            if group in ("pulse", "hardware"):
                self._sync_config()
            # run
            wr = self.run(measure=True, reconstruct=True)
            results.append(wr)
            # compute per-result metrics
            m = _compute_result_metrics(wr)
            for k in metrics:
                metrics[k].append(m.get(k, np.nan))

        sweep_result = SweepResult(
            param_path=param_path,
            values=list(values),
            results=results,
            metrics=metrics,
        )
        self._last_sweep = sweep_result
        return sweep_result

    # ==================================================================
    # compare() — A/B method comparison
    # ==================================================================

    def compare(
        self,
        methods: list[str],
        measurement: ExperimentResult | None = None,
        method_kwargs: dict[str, dict] | None = None,
    ) -> CompareResult:
        """Compare multiple reconstruction algorithms on the same measurement.

        Parameters
        ----------
        methods : list[str]
            Reconstruction method names (e.g., ``["wiener", "lm"]``).
        measurement : ExperimentResult or None
            Measurement data to reconstruct from.  If None, uses the
            cached measurement from the most recent ``run()`` call.
        method_kwargs : dict or None
            Per-method keyword overrides, e.g.
            ``{"wiener": {"lambda_reg": 10.0}}``.

        Returns
        -------
        CompareResult
            ``.signals`` maps method → FluxSignal;
            ``.metrics`` maps method → dict of metrics;
            ``.best`` names the lowest-RMSE method.

        Raises
        ------
        ValueError
            If no measurement data is available.
        """
        if measurement is None:
            measurement = self._cached_measurement
        if measurement is None:
            raise ValueError(
                "No measurement data available for comparison. "
                "Call run(measure=True) first or pass measurement=..."
            )

        protocol = self._params["protocol"]
        qubit = self._build_qubit()
        method_kwargs = method_kwargs or {}

        signals: dict[str, FluxSignal] = {}
        all_metrics: dict[str, dict] = {}

        # ground truth for RMSE (if available)
        ground_truth: np.ndarray | None = None
        if "flux_samples" in measurement.data:
            ground_truth = np.asarray(measurement.data["flux_samples"])

        t_start = time.perf_counter()

        for method in methods:
            self._params["reconstruction"]["method"] = method
            extra = method_kwargs.get(method, {})
            rec_signal, _rec_details = self._run_reconstruction(
                measurement, qubit, None, protocol,
                extra_kwargs=extra,
            )
            signals[method] = rec_signal
            all_metrics[method] = _compute_signal_metrics(
                rec_signal, ground_truth,
            )
            all_metrics[method]["runtime"] = time.perf_counter() - t_start
            t_start = time.perf_counter()

        # pick best by RMSE
        best_method = ""
        best_rmse = float("inf")
        for mth in methods:
            rmse = all_metrics[mth].get("rmse", float("inf"))
            if rmse < best_rmse:
                best_rmse = rmse
                best_method = mth

        compare_result = CompareResult(
            methods=list(methods),
            signals=signals,
            metrics=all_metrics,
            best=best_method,
        )
        self._last_compare = compare_result
        return compare_result

    # ==================================================================
    # plot() — auto-dispatch visualisation
    # ==================================================================

    def plot(self):
        """Auto-dispatch visualisation based on what was most recently run.

        - ``run(measure=True)``: p_e + delta_p vs scan axis
        - ``run(reconstruct=True)``: delta_p + reconstructed B(t)
        - ``sweep()``: metric vs swept-values plot
        - ``compare()``: multi-method B(t) overlay + residual subplot
        """
        import matplotlib.pyplot as plt

        # ── compare() was last ───────────────────────────────────────
        if self._last_compare is not None and self._last_compare.signals:
            self._plot_compare(self._last_compare)
            return

        # ── sweep() was last ─────────────────────────────────────────
        if self._last_sweep is not None and self._last_sweep.results:
            self._plot_sweep(self._last_sweep)
            return

        # ── run() was last ───────────────────────────────────────────
        wr = self._last_result
        if wr is None:
            print("Nothing to plot — run(), sweep(), or compare() first.")
            return

        has_measurement = wr.measurement is not None
        has_reconstruction = wr.reconstructed_signal is not None

        if has_reconstruction:
            self._plot_reconstruction(wr)
        elif has_measurement:
            self._plot_measurement(wr.measurement)
        else:
            print("Nothing to plot — result is empty.")

    # ── plot helpers ─────────────────────────────────────────────────

    def _plot_measurement(self, measurement: ExperimentResult):
        import matplotlib.pyplot as plt
        data = measurement.data
        axes = measurement.axes
        fig, axs = plt.subplots(1, 2, figsize=(10, 4))

        # p_e
        if "p_e" in data and "scan" in axes:
            axs[0].plot(axes["scan"], data["p_e"], "b.-")
            axs[0].set_xlabel("scan (ns)")
            axs[0].set_ylabel("p_e")
            axs[0].set_title("Excited-state population")

        # delta_p
        if "delta_p" in data and "scan" in axes:
            axs[1].plot(axes["scan"], data["delta_p"], "r.-")
            axs[1].set_xlabel("scan (ns)")
            axs[1].set_ylabel("Δp")
            axs[1].set_title("Delta population")
        elif "p_e" in data and "tau" in axes:
            axs[1].plot(axes["tau"], data["p_e"], "r.-")
            axs[1].set_xlabel("tau (ns)")
            axs[1].set_ylabel("p_e")
            axs[1].set_title("p_e vs free evolution")

        fig.tight_layout()
        plt.show()

    def _plot_reconstruction(self, wr: WorkflowResult):
        import matplotlib.pyplot as plt
        meas = wr.measurement
        rec = wr.reconstructed_signal
        fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=False)

        # top: delta_p
        if meas is not None and "delta_p" in meas.data:
            scan = meas.axes.get("scan",
                                 meas.axes.get("tau",
                                               meas.axes.get("trunc",
                                                             meas.axes.get("t_d"))))
            axs[0].plot(scan, meas.data["delta_p"], "b.-")
            axs[0].set_ylabel("Δp")
            axs[0].set_title("Measured Δp")
        elif meas is not None and "p_e" in meas.data:
            scan = meas.axes.get("scan",
                                 meas.axes.get("tau",
                                               meas.axes.get("trunc",
                                                             meas.axes.get("t_d"))))
            axs[0].plot(scan, meas.data["p_e"], "b.-")
            axs[0].set_ylabel("p_e")
            axs[0].set_title("Measured p_e")

        # bottom: reconstructed B(t)
        if rec is not None:
            axs[1].plot(rec.t_list, rec.signal, "r-", label="reconstructed")
            # overlay ground truth if available
            if meas is not None and "flux_samples" in meas.data:
                gt = np.asarray(meas.data["flux_samples"])
                t_flux = meas.axes.get("t_flux",
                                       meas.axes.get("t_samples"))
                if t_flux is not None:
                    axs[1].plot(t_flux, gt, "k--", alpha=0.5, label="ground truth")
            axs[1].set_xlabel("time (ns)")
            axs[1].set_ylabel("B (Φ₀)")
            axs[1].set_title("Reconstructed flux waveform")
            axs[1].legend()

        fig.tight_layout()
        plt.show()

    def _plot_sweep(self, sweep: SweepResult):
        import matplotlib.pyplot as plt
        n_metrics = len(sweep.metrics)
        if n_metrics == 0:
            print("No metrics to plot.")
            return
        fig, axs = plt.subplots(1, n_metrics, figsize=(4 * n_metrics, 3))
        if n_metrics == 1:
            axs = [axs]
        for ax, (key, vals) in zip(axs, sweep.metrics.items()):
            ax.plot(sweep.values, vals, "o-")
            ax.set_xlabel(sweep.param_path)
            ax.set_ylabel(key)
            ax.set_title(f"{key} vs {sweep.param_path}")
        fig.tight_layout()
        plt.show()

        # overlay plot of all B(t)
        fig2, ax2 = plt.subplots(figsize=(8, 4))
        for i, wr in enumerate(sweep.results):
            if wr.reconstructed_signal is not None:
                ax2.plot(wr.reconstructed_signal.t_list,
                         wr.reconstructed_signal.signal,
                         alpha=0.7, label=f"{sweep.values[i]:.3f}")
        ax2.set_xlabel("time (ns)")
        ax2.set_ylabel("B (Φ₀)")
        ax2.set_title(f"Overlay: B(t) vs {sweep.param_path}")
        ax2.legend(fontsize=8)
        fig2.tight_layout()
        plt.show()

    def _plot_compare(self, cmp_res: CompareResult):
        import matplotlib.pyplot as plt
        n = len(cmp_res.methods)
        if n == 0:
            return
        fig, (ax_sig, ax_res) = plt.subplots(2, 1, figsize=(8, 6))

        colors = plt.cm.tab10(np.linspace(0, 1, n))
        for i, method in enumerate(cmp_res.methods):
            sig = cmp_res.signals[method]
            ax_sig.plot(sig.t_list, sig.signal, color=colors[i],
                        label=method, alpha=0.8)

        ax_sig.set_xlabel("time (ns)")
        ax_sig.set_ylabel("B (Φ₀)")
        ax_sig.set_title("Reconstruction method comparison")
        ax_sig.legend()

        # residual subplot (difference from mean)
        if n >= 2:
            # compute ensemble mean
            all_signals = []
            common_t = None
            for method in cmp_res.methods:
                s = np.asarray(cmp_res.signals[method].signal)
                all_signals.append(s)
                if common_t is None:
                    common_t = cmp_res.signals[method].t_list
            all_signals = np.array(all_signals)
            mean_sig = np.mean(all_signals, axis=0)
            for i, method in enumerate(cmp_res.methods):
                residual = all_signals[i] - mean_sig
                ax_res.plot(common_t, residual, color=colors[i],
                            label=f"{method} - mean")
            ax_res.set_xlabel("time (ns)")
            ax_res.set_ylabel("residual (Φ₀)")
            ax_res.set_title("Residual vs ensemble mean")
            ax_res.legend(fontsize=8)

        fig.tight_layout()
        plt.show()

    # ==================================================================
    # Internal helpers
    # ==================================================================

    def _build_qubit(self):
        """Construct a TransmonQubit from _params['qubit']."""
        from sqc.devices.transmon import TransmonQubit
        qp = self._params["qubit"]
        return TransmonQubit(
            EC=2 * np.pi * qp["EC"],
            EJ=2 * np.pi * qp["EJ"],
            T1=qp["T1"],
            T2=qp["T2"],
            flux=qp["flux_bias"],
            n_levels=qp["n_levels"],
        )

    def _sync_config(self):
        """Rebuild _cfg from pulse/hardware params."""
        pp = self._params["pulse"]
        hp = self._params["hardware"]
        self._cfg = _reconfigure(
            sample_rate=hp["sample_rate"],
            t_rabi_duration=pp["t_rabi_duration"],
            t_global_start=pp["t_global_start"],
            t_global_end=pp["t_global_end"],
        )

    def _build_flux_signal(self):
        """Create a FluxSignal from _params['signal']."""
        sp = self._params["signal"]
        pp = self._params["pulse"]
        t_signal = np.arange(0, pp["t_global_end"] - pp["t_global_start"],
                             self._cfg.awg.dt)

        if sp["custom"] is not None:
            return FluxSignal(type=8, t_list=t_signal,
                            signal=np.asarray(sp["custom"]))

        return FluxSignal(
            type=sp["type"],
            t_list=t_signal,
            amplitude=sp["amplitude"],
            width=sp.get("width", 10),
            center=sp.get("center", 100),
            rise=sp.get("rise", 10),
            fall=sp.get("fall", 10),
            frequency=sp.get("frequency", 0.01),
        )

    def _make_experiment(self, qubit, flux_signal, protocol: str):
        """Instantiate the appropriate Experiment subclass."""
        _ = flux_signal  # some experiments ignore it and use their own default
        exp_name = _PROTOCOL_MAP[protocol][0]

        if protocol == "transient":
            from sqc.experiments.transient import TransientSensingExperiment
            return TransientSensingExperiment(
                qubit=qubit,
                flux_signal=flux_signal,
                t_rabi=self._cfg.pulse.t_rabi.copy(),
                rotation_angle=self._params["pulse"]["rotation_angle"],
                rabi_rate=self._params["pulse"]["rabi_rate"],
                envelope=self._params["pulse"]["envelope"],
                envelope_sigma=self._params["pulse"]["envelope_sigma"],
                phase1=self._params["pulse"]["phase1"],
                phase2=self._params["pulse"]["phase2"],
            )
        elif protocol == "ramsey":
            from sqc.experiments.ramsey import RamseyExperiment
            return RamseyExperiment(
                qubit=qubit,
                flux_signal=flux_signal,
                t_rabi=self._cfg.pulse.t_rabi.copy(),
                t_global=self._cfg.pulse.t_global.copy(),
                tau_list=self._cfg.pulse.tau_list.copy(),
            )
        elif protocol == "echo":
            from sqc.experiments.echo import DiffEchoExperiment
            return DiffEchoExperiment(
                qubit=qubit,
                flux_signal=flux_signal,
                t_rabi=self._cfg.pulse.t_rabi.copy(),
            )
        elif protocol == "cryoscope":
            from sqc.experiments.cryoscope import CryoscopeExperiment
            return CryoscopeExperiment(
                qubit=qubit,
                flux_signal=flux_signal,
                t_rabi=self._cfg.pulse.t_rabi.copy(),
            )
        elif protocol == "delay_ramsey":
            from sqc.experiments.delay_ramsey import DelayRamseyExperiment
            return DelayRamseyExperiment(
                qubit=qubit,
                flux_signal=flux_signal,
                t_rabi=self._cfg.pulse.t_rabi.copy(),
            )
        else:
            raise ValueError(f"Unknown protocol: {protocol}")

    def _run_measurement(self, qubit) -> ExperimentResult:
        """Build experiment + run measurement."""
        protocol = self._params["protocol"]
        flux_signal = self._build_flux_signal()
        experiment = self._make_experiment(qubit, flux_signal, protocol)
        return experiment.run()

    def _run_calibration(self, qubit, protocol: str):
        """Run calibration for protocols that require it."""
        if protocol == "cryoscope":
            from sqc.reconstruction.cryoscope import CryoscopeCalibration
            cal = CryoscopeCalibration(
                qubit=qubit,
                t_rabi=self._cfg.pulse.t_rabi.copy(),
            )
            return cal.calibrate()
        elif protocol == "delay_ramsey":
            from sqc.reconstruction.delay_ramsey import DelayRamseyCalibration
            cal = DelayRamseyCalibration(
                qubit=qubit,
                t_rabi=self._cfg.pulse.t_rabi.copy(),
            )
            return cal.calibrate()
        else:
            return None

    def _run_reconstruction(
        self, measurement, qubit, cal_table, protocol: str,
        extra_kwargs: dict | None = None,
    ) -> tuple[FluxSignal | None, dict]:
        """Run waveform reconstruction on measurement data."""
        rp = self._params["reconstruction"]
        method = rp["method"]
        extra_kwargs = extra_kwargs or {}
        details: dict = {}

        if protocol == "transient":
            from sqc.reconstruction.transient import TransientReconstruction
            kernel = measurement.data.get("kernel", None)
            rec = TransientReconstruction(
                method=method,
                lambda_reg=extra_kwargs.get("lambda_reg", rp["lambda_reg"]),
                qubit=qubit,
                n_basis=extra_kwargs.get("lm_n_basis", rp["lm_n_basis"]),
                max_iter=extra_kwargs.get("lm_max_iter", rp["lm_max_iter"]),
                basis_type=extra_kwargs.get("lm_basis_type", rp["lm_basis_type"]),
                use_adjoint=extra_kwargs.get("use_adjoint", rp["use_adjoint"]),
            )
            rec_signal = rec.reconstruct(measurement, kernel=kernel)
            return rec_signal, details

        elif protocol == "ramsey":
            from sqc.reconstruction.ramsey import RamseyReconstruction
            rec = RamseyReconstruction(
                qubit=qubit,
                method=method,
            )
            # Use the offset-corrected absolute-time axis so the wrapped
            # FluxSignal lines up with the true flux (free evolution starts
            # after the first pi/2 pulse, not at global t=0).
            if "tau" in measurement.axes:
                t_axis, B = rec.reconstruct_with_time(measurement)
                fs = FluxSignal(type=8, t_list=t_axis, signal=B)
            else:
                B = rec.reconstruct(measurement)
                tau = measurement.axes.get("scan")
                if tau is not None:
                    fs = FluxSignal(type=8, t_list=tau, signal=B)
                else:
                    fs = FluxSignal(type=8, t_list=np.arange(len(B)), signal=B)
            return fs, details

        elif protocol == "echo":
            from sqc.reconstruction.echo import EchoReconstruction
            t_int = measurement.config.get("t_int", self._cfg.pulse.t_rabi_duration * 0.5)
            k = measurement.metadata.get("k", 5)
            rec = EchoReconstruction(
                qubit=qubit,
                t_int=t_int,
                k=k,
            )
            B = rec.reconstruct(measurement)
            tau = measurement.axes.get("tau")
            if tau is not None:
                fs = FluxSignal(type=8, t_list=tau, signal=B)
            else:
                fs = FluxSignal(type=8, t_list=np.arange(len(B)), signal=B)
            return fs, details

        elif protocol == "cryoscope":
            from sqc.reconstruction.cryoscope import CryoscopeReconstruction
            rec = CryoscopeReconstruction(
                inversion=method,
                calibration=cal_table,
                qubit=qubit,
            )
            rec_signal = rec.reconstruct(measurement, calibration=cal_table)
            return rec_signal, details

        elif protocol == "delay_ramsey":
            from sqc.reconstruction.delay_ramsey import DelayRamseyReconstruction
            rec = DelayRamseyReconstruction(
                inversion=method,
                calibration=cal_table,
                qubit=qubit,
            )
            rec_signal = rec.reconstruct(measurement)
            return rec_signal, details

        else:
            raise ValueError(f"Unknown protocol: {protocol}")

    def _parse_param_path(self, param_path: str) -> tuple[str, str]:
        """Parse 'group.field' into (group, field)."""
        parts = param_path.split(".", 1)
        if len(parts) != 2:
            raise ValueError(
                f"Invalid param_path '{param_path}'. "
                "Expected format: 'group.field' (e.g., 'signal.amplitude')"
            )
        group, field = parts
        if group not in self._params:
            raise ValueError(
                f"Unknown param group '{group}'. "
                f"Valid groups: {sorted(self._params.keys())}"
            )
        if field not in self._params[group]:
            raise ValueError(
                f"Unknown field '{field}' in group '{group}'. "
                f"Valid fields: {sorted(self._params[group].keys())}"
            )
        return group, field

    # ==================================================================
    # Stub methods — planned for P6.1+
    # ==================================================================

    def pipeline(self, stages: list[tuple]) -> WorkflowResult:
        """⚠️ (stub) Chain multiple stages: calibrate → measure → reconstruct → predistort → verify.

        Parameters
        ----------
        stages : list[tuple[str, dict]]
            Ordered list of ``(stage_name, stage_kwargs)`` pairs.
            Supported stage names: ``"calibrate"``, ``"measure"``,
            ``"reconstruct"``, ``"predistort"``, ``"verify"``.

        Returns
        -------
        WorkflowResult
        """
        raise NotImplementedError(
            "pipeline() is planned for P6.1+. "
            "Use run(measure=True, reconstruct=True, calibrate=True) "
            "for the basic pipeline."
        )

    def multi_qubit(self, chip_spec: dict) -> "SensingWorkflow":
        """⚠️ (stub) Register a multi-qubit chip topology.

        Parameters
        ----------
        chip_spec : dict[str, dict]
            Mapping from qubit name to its parameters, e.g.::

                {"QA": {"EC": 0.20, "EJ": 15.0, "flux_bias": 0.0},
                 "QB": {"EC": 0.22, "EJ": 14.5, "flux_bias": 0.0}}

        Returns
        -------
        self
        """
        raise NotImplementedError("multi_qubit() is planned for P6.1+.")

    def crosstalk(self, drive: str, sense: str, **params) -> WorkflowResult:
        """⚠️ (stub) Measure Z-crosstalk transfer function between two qubits.

        Parameters
        ----------
        drive : str
            Name of the driven qubit.
        sense : str
            Name of the sensing qubit.
        **params
            Passed to ``ZCrosstalkWorkflow`` (flux_amplitude, etc.).

        Returns
        -------
        WorkflowResult
        """
        raise NotImplementedError("crosstalk() is planned for P6.1+.")

    def save(self, path: str) -> None:
        """💡 (stub) Persist the last WorkflowResult to a pickle file.

        Parameters
        ----------
        path : str
            Output file path (``.pkl``).
        """
        raise NotImplementedError("save() is planned for a future phase.")

    def load(self, path: str) -> "SensingWorkflow":
        """💡 (stub) Restore a SensingWorkflow (including its last result)
        from a pickle file.

        Parameters
        ----------
        path : str
            Input file path (``.pkl``).

        Returns
        -------
        self
        """
        raise NotImplementedError("load() is planned for a future phase.")

    def diff(self, other: WorkflowResult | str) -> DiffReport:
        """💡 (stub) Compare two WorkflowResults and produce a structured
        difference report.

        Parameters
        ----------
        other : WorkflowResult or str
            Another result object or a path to a saved result.

        Returns
        -------
        DiffReport
        """
        raise NotImplementedError("diff() is planned for a future phase.")

    def benchmark(self, configs: list[dict]) -> BenchmarkResult:
        """💡 (stub) Batch benchmark: run a set of protocol/algorithm
        combinations and produce a ranking table.

        Parameters
        ----------
        configs : list[dict]
            Each dict is a full ``configure()``-style parameter set.

        Returns
        -------
        BenchmarkResult
        """
        raise NotImplementedError("benchmark() is planned for a future phase.")

    def find_optimal_work_point(
        self, metric: str = "sensitivity",
        flux_range: tuple[float, float] = (-0.5, 0.5),
    ) -> float:
        """💡 (stub) Search for the flux bias that maximises
        ``dω/dΦ`` (sensitivity) over a given range.

        Parameters
        ----------
        metric : str
            ``"sensitivity"`` — maximise |df/dΦ|.
        flux_range : tuple[float, float]
            Search range in Φ₀.

        Returns
        -------
        float
            Optimal flux_bias value.
        """
        raise NotImplementedError(
            "find_optimal_work_point() is planned for a future phase."
        )

    def detectability_limit(
        self,
        protocol: str | None = None,
        signal_type: int | None = None,
        amplitude_range: tuple[float, float] = (1e-6, 0.1),
        confidence: float = 0.95,
    ) -> float:
        """💡 (stub) Binary search for the minimum detectable signal
        amplitude for a given protocol.

        Parameters
        ----------
        protocol : str or None
            Protocol name. Uses current if None.
        signal_type : int or None
            FluxSignal type. Uses current if None.
        amplitude_range : tuple[float, float]
            Search range in Φ₀.
        confidence : float
            Desired confidence level (0–1).

        Returns
        -------
        float
            Minimum detectable amplitude (Φ₀).
        """
        raise NotImplementedError(
            "detectability_limit() is planned for a future phase."
        )

    def noise_characterize(
        self, n_repeats: int = 100, protocol: str | None = None,
    ) -> NoiseReport:
        """💡 (stub) Characterise measurement noise by running N
        zero-signal measurements and analysing the statistics.

        Parameters
        ----------
        n_repeats : int
            Number of repeat measurements.
        protocol : str or None
            Protocol to use. Uses current if None.

        Returns
        -------
        NoiseReport
            Mean, std, PSD, and normality test results.
        """
        raise NotImplementedError(
            "noise_characterize() is planned for a future phase."
        )

    def cross_validate(
        self, k_folds: int = 5, reconstruction: str | None = None,
    ) -> CVResult:
        """💡 (stub) K-fold cross-validation of reconstruction quality.

        Parameters
        ----------
        k_folds : int
            Number of folds.
        reconstruction : str or None
            Reconstruction method. Uses current if None.

        Returns
        -------
        CVResult
        """
        raise NotImplementedError(
            "cross_validate() is planned for a future phase."
        )


# ═══════════════════════════════════════════════════════════════════════
# Metric helpers (module-level)
# ═══════════════════════════════════════════════════════════════════════

def _compute_signal_metrics(
    signal: FluxSignal,
    ground_truth: np.ndarray | None = None,
) -> dict[str, float]:
    """Compute SNR, RMSE, and peak amplitude for a FluxSignal."""
    sig = np.asarray(signal.signal, dtype=float)
    metrics: dict[str, float] = {}
    # peak
    metrics["peak"] = float(np.max(np.abs(sig)))

    # SNR: max |signal| / std of tail (or entire signal if no clear tail)
    tail = sig[int(0.8 * len(sig)):]
    noise_std = float(np.std(tail)) if len(tail) > 0 else float(np.std(sig))
    metrics["snr"] = metrics["peak"] / (noise_std + 1e-12)

    # RMSE against ground truth
    if ground_truth is not None:
        gt = np.asarray(ground_truth, dtype=float)
        # interpolate to match lengths if needed
        if len(sig) != len(gt):
            if len(sig) < len(gt):
                gt = gt[:len(sig)]
            else:
                sig_aligned = sig[:len(gt)]
                metrics["rmse"] = float(np.sqrt(np.mean((sig_aligned - gt) ** 2)))
                return metrics
        metrics["rmse"] = float(np.sqrt(np.mean((sig - gt) ** 2)))
    else:
        metrics["rmse"] = float("nan")

    return metrics


def _compute_result_metrics(wr: WorkflowResult) -> dict[str, float]:
    """Compute aggregate metrics for a WorkflowResult."""
    metrics: dict[str, float] = {"snr": float("nan"),
                                  "rmse": float("nan"),
                                  "peak": float("nan")}
    if wr.reconstructed_signal is not None:
        ground_truth = None
        if wr.measurement is not None and "flux_samples" in wr.measurement.data:
            ground_truth = np.asarray(wr.measurement.data["flux_samples"])
        return _compute_signal_metrics(wr.reconstructed_signal, ground_truth)
    return metrics
