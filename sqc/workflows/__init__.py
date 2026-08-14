"""sqc.workflows — High-level experimental workflows."""

from .base import Workflow
from .frequency_calibration import CalibrationStage, FrequencyCalibrationWorkflow
from .predistortion_validation import PredistortionValidationWorkflow

# NOTE (D3): ZCrosstalkWorkflow is hidden from the v1 public API. The
# implementation remains in sqc/workflows/z_crosstalk.py and is importable via
# its deep path (`from sqc.workflows.z_crosstalk import ZCrosstalkWorkflow`);
# it is intentionally not re-exported here. Restore by re-adding the import and
# the "ZCrosstalkWorkflow" entry to __all__.
from .sensing import (
    SensingWorkflow,
    WorkflowResult,
    SweepResult,
    CompareResult,
    DiffReport,
    BenchmarkResult,
    NoiseReport,
    CVResult,
)

# V2 event-driven frequency calibration state machine (v2.19)
from .frequency_state_machine import (
    FrequencyCalibrationConfig,
    FrequencyStateMachine,
    FrequencyState,
    RunStatus,
    ReasonCode,
)
from .frequency_runtime import CancellationToken, FrequencyCalibrationRuntime
from .frequency_backends import (
    SQCExecutor,
    FiniteShotSQCExecutor,
    ProcessSimulationExecutor,
)

__all__ = [
    "Workflow",
    "FrequencyCalibrationWorkflow",
    "CalibrationStage",
    "PredistortionValidationWorkflow",
    # "ZCrosstalkWorkflow" hidden from v1 public API (D3); see note above.
    "SensingWorkflow",
    "WorkflowResult",
    "SweepResult",
    "CompareResult",
    "DiffReport",
    "BenchmarkResult",
    "NoiseReport",
    "CVResult",
    # V2 event-driven frequency calibration (v2.19)
    "FrequencyCalibrationConfig",
    "FrequencyStateMachine",
    "FrequencyState",
    "RunStatus",
    "ReasonCode",
    "FrequencyCalibrationRuntime",
    "CancellationToken",
    "SQCExecutor",
    "FiniteShotSQCExecutor",
    "ProcessSimulationExecutor",
]
