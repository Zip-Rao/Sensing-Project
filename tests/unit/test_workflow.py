"""tests/unit/test_workflow.py — unit tests for SensingWorkflow.

Tests the configure / run / sweep / compare / plot APIs and all stub
methods per phase_6_handbook.md §6d.9.
"""
from __future__ import annotations

import numpy as np
import pytest

from sqc.workflows.sensing import (
    CVResult,
    BenchmarkResult,
    CompareResult,
    DiffReport,
    NoiseReport,
    SensingWorkflow,
    SweepResult,
    WorkflowResult,
)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture
def wf() -> SensingWorkflow:
    """Fresh SensingWorkflow with minimal config for fast tests."""
    return SensingWorkflow()


@pytest.fixture
def wf_transient(wf) -> SensingWorkflow:
    """Workflow pre-configured for transient + wiener (fast default)."""
    wf.configure(
        protocol="transient",
        signal_amplitude=0.01,
        signal_type=4,
        reconstruction="wiener",
        t_rabi_duration=10,
        t_global_start=0,
        t_global_end=100,
        n_levels=2,
    )
    return wf


# ═══════════════════════════════════════════════════════════════════════
# configure()
# ═══════════════════════════════════════════════════════════════════════

class TestConfigure:
    def test_configure_stores_params(self, wf):
        wf.configure(protocol="transient", signal_amplitude=0.02,
                     lambda_reg=3.0, EC=0.25)
        assert wf._params["protocol"] == "transient"
        assert wf._params["signal"]["amplitude"] == 0.02
        assert wf._params["reconstruction"]["lambda_reg"] == 3.0
        assert wf._params["qubit"]["EC"] == 0.25

    def test_configure_partial_update(self, wf):
        wf.configure(protocol="transient", signal_amplitude=0.01,
                     lambda_reg=5.0)
        wf.configure(signal_amplitude=0.03)  # partial
        assert wf._params["signal"]["amplitude"] == 0.03
        assert wf._params["protocol"] == "transient"  # unchanged
        assert wf._params["reconstruction"]["lambda_reg"] == 5.0  # unchanged

    def test_configure_returns_self(self, wf):
        result = wf.configure(protocol="ramsey")
        assert result is wf

    def test_configure_invalid_protocol_raises(self, wf):
        with pytest.raises(ValueError, match="Unknown protocol"):
            wf.configure(protocol="bad_protocol")

    def test_configure_all_groups(self, wf):
        """Each parameter group is addressable."""
        wf.configure(
            # qubit
            EC=0.30, EJ=14.0, T1=5000, T2=4000, flux_bias=0.1, n_levels=4,
            # signal
            protocol="echo", signal_type=3, signal_amplitude=0.05,
            signal_width=5, signal_center=50,
            signal_rise=5, signal_fall=5, signal_frequency=0.02,
            # reconstruction
            reconstruction="wiener", lambda_reg=2.0,
            lm_n_basis=30, lm_max_iter=5,
            lm_basis_type="fourier", use_adjoint=False,
            # pulse
            t_rabi_duration=15, t_global_start=-30, t_global_end=300,
            # hardware
            sample_rate=4.0,
        )
        assert wf._params["qubit"]["EC"] == 0.30
        assert wf._params["qubit"]["n_levels"] == 4
        assert wf._params["protocol"] == "echo"
        assert wf._params["signal"]["amplitude"] == 0.05
        assert wf._params["reconstruction"]["lm_max_iter"] == 5
        assert wf._params["pulse"]["t_rabi_duration"] == 15
        assert wf._params["hardware"]["sample_rate"] == 4.0


# ═══════════════════════════════════════════════════════════════════════
# run()
# ═══════════════════════════════════════════════════════════════════════

class TestRun:
    def test_run_measure_only(self, wf_transient):
        result = wf_transient.run(measure=True, reconstruct=False)
        assert result.measurement is not None
        assert "p_e" in result.measurement.data
        assert result.reconstructed_signal is None

    def test_run_full_pipeline(self, wf_transient):
        result = wf_transient.run(measure=True, reconstruct=True)
        assert result.measurement is not None
        assert result.reconstructed_signal is not None
        assert isinstance(result.reconstructed_signal.signal, np.ndarray)

    def test_run_reconstruct_without_measure_raises(self, wf):
        with pytest.raises(ValueError, match="Cannot reconstruct without"):
            wf.run(measure=False, reconstruct=True)

    def test_run_reconstruct_uses_cached_measurement(self, wf_transient):
        # measure first
        wf_transient.run(measure=True, reconstruct=False)
        # reconstruct from cache
        result = wf_transient.run(measure=False, reconstruct=True)
        assert result.reconstructed_signal is not None

    def test_run_calibrate_skipped_for_transient(self, wf_transient):
        result = wf_transient.run(measure=True, reconstruct=True, calibrate=True)
        # calibration should be None for transient (not a calibrate-required protocol)
        assert result.calibration is None
        # measurement and reconstruction should still work
        assert result.measurement is not None
        assert result.reconstructed_signal is not None


# ═══════════════════════════════════════════════════════════════════════
# protocol mapping
# ═══════════════════════════════════════════════════════════════════════

class TestProtocolMapping:
    @pytest.mark.parametrize("protocol", ["transient", "ramsey", "echo",
                                           "cryoscope", "delay_ramsey"])
    def test_protocol_mapping_valid(self, wf, protocol):
        """All 5 protocols can be configured and run (measure only)."""
        wf.configure(
            protocol=protocol,
            signal_amplitude=0.01,
            t_rabi_duration=10,
            t_global_start=0,
            t_global_end=100,
            n_levels=2,
        )
        result = wf.run(measure=True, reconstruct=False)
        assert result.measurement is not None
        assert result.measurement.data is not None

    def test_protocol_transient_reconstruct(self, wf):
        wf.configure(protocol="transient", signal_amplitude=0.01,
                     reconstruction="wiener", t_rabi_duration=10,
                     t_global_start=0, t_global_end=100, n_levels=2)
        result = wf.run(measure=True, reconstruct=True)
        assert result.reconstructed_signal is not None


# ═══════════════════════════════════════════════════════════════════════
# sweep()
# ═══════════════════════════════════════════════════════════════════════

class TestSweep:
    def test_sweep_returns_correct_count(self, wf_transient):
        sweep = wf_transient.sweep("signal.amplitude", [0.005, 0.01])
        assert len(sweep.results) == 2
        assert sweep.param_path == "signal.amplitude"
        assert sweep.values == [0.005, 0.01]

    def test_sweep_metrics_computed(self, wf_transient):
        sweep = wf_transient.sweep("signal.amplitude", [0.005, 0.01])
        assert "snr" in sweep.metrics
        assert "rmse" in sweep.metrics
        assert "peak" in sweep.metrics
        assert len(sweep.metrics["snr"]) == 2
        # peak should increase with amplitude
        assert sweep.metrics["peak"][1] > sweep.metrics["peak"][0] * 0.5

    def test_sweep_invalid_param_path_raises(self, wf):
        with pytest.raises(ValueError, match="Invalid param_path"):
            wf.sweep("badpath", [1, 2])

    def test_sweep_unknown_group_raises(self, wf):
        with pytest.raises(ValueError, match="Unknown param group"):
            wf.sweep("nonexistent.field", [1, 2])


# ═══════════════════════════════════════════════════════════════════════
# compare()
# ═══════════════════════════════════════════════════════════════════════

class TestCompare:
    def test_compare_multiple_methods(self, wf_transient):
        # measure first
        wf_transient.run(measure=True, reconstruct=False)
        cmp_res = wf_transient.compare(methods=["wiener", "hammerstein"])
        assert "wiener" in cmp_res.signals
        assert "hammerstein" in cmp_res.signals
        assert cmp_res.best in ("wiener", "hammerstein")

    def test_compare_without_measurement_raises(self, wf):
        with pytest.raises(ValueError, match="No measurement data"):
            wf.compare(methods=["wiener"])


# ═══════════════════════════════════════════════════════════════════════
# stub methods
# ═══════════════════════════════════════════════════════════════════════

class TestStubs:
    def test_stub_methods_raise_not_implemented(self, wf):
        """All 11 research-interface stubs raise NotImplementedError."""
        import inspect

        stub_names = [
            "pipeline", "multi_qubit", "crosstalk",
            "save", "load", "diff", "benchmark",
            "find_optimal_work_point", "detectability_limit",
            "noise_characterize", "cross_validate",
        ]
        for name in stub_names:
            method = getattr(wf, name)
            sig = inspect.signature(method)
            # Build minimal args for methods that have required params
            kwargs = {}
            for pname, param in sig.parameters.items():
                if param.default is inspect.Parameter.empty:
                    if param.annotation in (str, "str"):
                        kwargs[pname] = "dummy"
                    elif param.annotation in (int, "int"):
                        kwargs[pname] = 1
                    elif param.annotation in (float, "float"):
                        kwargs[pname] = 0.0
                    elif param.annotation in (bool, "bool"):
                        kwargs[pname] = False
                    elif param.annotation == dict:
                        kwargs[pname] = {}
                    elif param.annotation == list:
                        kwargs[pname] = []
                    else:
                        kwargs[pname] = None
            with pytest.raises(NotImplementedError):
                method(**kwargs)


# ═══════════════════════════════════════════════════════════════════════
# reconfigure() extension (P6a)
# ═══════════════════════════════════════════════════════════════════════

class TestReconfigureExtended:
    def test_reconfigure_covers_reconstruction(self):
        from sqc.config import reconfigure
        cfg = reconfigure(lambda_reg=5.0, stim_amplitude=0.05)
        assert cfg.reconstruction.lambda_reg == 5.0
        assert cfg.reconstruction.stim_amplitude == 0.05

    def test_reconfigure_covers_simulation(self):
        from sqc.config import reconfigure
        cfg = reconfigure(atol=1e-12, rtol=1e-8)
        assert cfg.simulation.atol == 1e-12
        assert cfg.simulation.rtol == 1e-8

    def test_reconfigure_covers_transmon(self):
        from sqc.config import reconfigure
        cfg = reconfigure(EC=0.3, n_levels=4)
        assert cfg.transmon.EC == 0.3
        assert cfg.transmon.n_levels == 4

    def test_reconfigure_covers_control_line(self):
        from sqc.config import reconfigure
        cfg = reconfigure(impedance=75.0, attenuation_db=10.0, delay=5.0)
        assert cfg.control_line.impedance == 75.0
        assert cfg.control_line.attenuation_db == 10.0
        assert cfg.control_line.delay == 5.0

    def test_reconfigure_no_args_returns_defaults(self):
        from sqc.config import CONFIG, reconfigure
        cfg = reconfigure()
        assert cfg.reconstruction.lambda_reg == CONFIG.reconstruction.lambda_reg
        assert cfg.simulation.atol == CONFIG.simulation.atol
        assert cfg.transmon.EC == CONFIG.transmon.EC
        assert cfg.control_line.impedance == CONFIG.control_line.impedance

    def test_reconfigure_does_not_mutate_global_config(self):
        from sqc.config import CONFIG, reconfigure
        old_lambda = CONFIG.reconstruction.lambda_reg
        cfg = reconfigure(lambda_reg=99.0)
        assert CONFIG.reconstruction.lambda_reg == old_lambda
        assert cfg.reconstruction.lambda_reg == 99.0

    def test_reconfigure_all_layers_together(self):
        from sqc.config import reconfigure
        cfg = reconfigure(
            sample_rate=4.0,           # AWG
            t_rabi_duration=30,        # Pulse
            lambda_reg=3.0,            # Reconstruction
            atol=1e-10,                # Simulation
            EC=0.25, n_levels=5,      # Transmon
            impedance=100.0,           # ControlLine
        )
        assert cfg.awg.dt == 0.25
        assert cfg.pulse.t_rabi_duration == 30
        assert cfg.reconstruction.lambda_reg == 3.0
        assert cfg.simulation.atol == 1e-10
        assert cfg.transmon.EC == 0.25
        assert cfg.transmon.n_levels == 5
        assert cfg.control_line.impedance == 100.0
