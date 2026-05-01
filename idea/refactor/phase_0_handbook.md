# Phase 0 Handbook — 测试基线 + 工具基础设施

> 前置阅读:[_refactor_plan.md](_refactor_plan.md) §1–§13、§0.1(本地物理参考 PDF 说明)  
> 物理参考(可选):打开本地 PDF [`./Gao 等 - 2021 - Practical Guide for Building Superconducting Quantum Devices.pdf`](./Gao%20%E7%AD%89%20-%202021%20-%20Practical%20Guide%20for%20Building%20Superconducting%20Quantum%20Devices.pdf) §I 浏览全栈架构,§II.C Transmon 物理参数推荐范围(确认 P0 baseline 选取合理性)  
> 估计工时:1–2 天  
> 触发条件:无,立即开始  
> 完成标志:`pytest tests/regression -x` 全部通过,baseline pickle 已生成并提交

---

## 1. 目标

在动任何代码之前,建立**物理回归测试基线**,使后续所有重构 phase 都可以通过运行 `pytest` 自动验证"物理结果不变"。

本 phase 不动 src/ 一行业务代码,只新增:
- `tests/` 目录及配置
- baseline 生成脚本
- pickle 化的 baseline 数据
- pytest 配置文件
- requirements 钉死

---

## 2. 前置条件

| # | 条件 | 验证方式 |
|---|---|---|
| 2.1 | git 工作区干净 | `git status` 输出 `nothing to commit` |
| 2.2 | 当前 src/ 可运行 | `python -c "from src.qubit import TransmonQubit; from src.protocal import Protocal; print('ok')"` |
| 2.3 | QuTiP 已安装 | `python -c "import qutip; print(qutip.__version__)"` |
| 2.4 | numpy/scipy/matplotlib 已安装 | `python -c "import numpy, scipy, matplotlib; print('ok')"` |

如任一条件不满足,**停止本 phase**,先解决。

---

## 3. 任务清单

### 3.1 钉死依赖版本

**文件**:[requirements.txt](../../requirements.txt)

**任务**:在文件中追加版本约束(若已有则保留)。

```text
qutip>=5.0,<6.0
numpy>=1.24,<3.0
scipy>=1.11,<2.0
matplotlib>=3.7,<4.0
gradio>=4.0,<5.0
pytest>=7.4,<9.0
pytest-xdist>=3.5,<4.0
pytest-cov>=4.1,<6.0
```

**验证**:`pip install -r requirements.txt` 成功。

### 3.2 创建 pytest 配置

**新建文件**:`pytest.ini`(项目根目录)

```ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts = 
    -ra
    --strict-markers
    --tb=short
markers =
    regression: physics regression tests (slow)
    integration: cross-module integration tests
    unit: pure unit tests (fast)
    slow: tests that take > 30s
filterwarnings =
    ignore::DeprecationWarning:qutip.*
    ignore::DeprecationWarning:scipy.*
```

### 3.3 创建 tests/ 目录骨架

**任务**:创建以下空文件(目录用 `__init__.py` 占位)。

```
tests/
├── __init__.py                         (空文件)
├── conftest.py                         (见 3.4)
├── baselines/                          (见 3.6)
│   └── README.md                       (说明 baseline 来源与刷新流程)
├── unit/
│   └── __init__.py                     (空文件)
├── integration/
│   └── __init__.py                     (空文件)
└── regression/
    ├── __init__.py                     (空文件)
    ├── generate_baselines.py           (见 3.5)
    └── test_physics_baseline.py        (见 3.7)
```

### 3.4 全局 fixture 文件

**新建文件**:`tests/conftest.py`

```python
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
        n_levels=3,
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
```

### 3.5 baseline 生成脚本

**新建文件**:`tests/regression/generate_baselines.py`

```python
"""Generate physics regression baselines.

Run once, before any refactor work begins:

    python -m tests.regression.generate_baselines

Re-run only when intentionally changing physics (require justification
in commit message).
"""
from __future__ import annotations

import sys
from pathlib import Path

# allow standalone execution
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
from qutip import basis

from src.qubit import TransmonQubit
from src.signal import Signal, CompositeSignal
from src.pulse import (
    create_ramsey_pulse,
    create_diff_echo_pulse,
)
from src.protocal import Protocal
from tests.conftest import save_baseline


def _make_default_qubit() -> TransmonQubit:
    return TransmonQubit(
        EC=2 * np.pi * 0.2,
        EJ=2 * np.pi * 15,
        T1=10000.0,
        T2=8000.0,
        flux=0.0,
        state=0,
        n_levels=3,
    )


def _baseline_qubit_static() -> dict:
    """Frequency, anharmonicity, sensitivity at default flux."""
    q = _make_default_qubit()
    return {
        "frequency": float(q.frequency),
        "anharmonicity": float(q.anharmonicity),
        "sensitivity_at_zero": float(q.frequency_sensitivity(0.0)),
        "sensitivity_at_optimal": float(
            q.frequency_sensitivity(q.optimal_work_point() / np.pi)
        ),
    }


def _baseline_ramsey() -> dict:
    """Run Protocal.evolve case 1 and snapshot outputs."""
    q = _make_default_qubit()
    proto = Protocal(type=1)
    proto.initialize(q, state=0)
    
    # NOTE: Protocal.evolve case 1 returns Phi, tau_list, p_e_list
    Phi, tau_list, p_e_list = proto.evolve(q)
    return {
        "Phi_signal": np.asarray(Phi.signal),
        "Phi_t_list": np.asarray(Phi.t_list),
        "tau_list": np.asarray(tau_list),
        "p_e_list": np.asarray(p_e_list),
    }


def _baseline_diff_echo() -> dict:
    """Run Protocal.evolve case 2."""
    q = _make_default_qubit()
    proto = Protocal(type=2)
    proto.initialize(q, state=0)
    Phi, tau_list, p_e_list, k, t_int = proto.evolve(q)
    return {
        "Phi_signal": np.asarray(Phi.signal),
        "Phi_t_list": np.asarray(Phi.t_list),
        "tau_list": np.asarray(tau_list),
        "p_e_list": np.asarray(p_e_list),
        "k": int(k),
        "t_int": float(t_int),
    }


def _baseline_transient() -> dict:
    """Run Protocal.evolve case 4 (sliding measurement)."""
    q = _make_default_qubit()
    proto = Protocal(type=4)
    proto.initialize(q, state=0)
    t_samples, kernel, scan_list, delta_p, p_e, Phi, ctrl = proto.evolve(q)
    return {
        "t_samples": np.asarray(t_samples),
        "kernel": np.asarray(kernel),
        "scan_list": np.asarray(scan_list),
        "delta_p": np.asarray(delta_p),
        "p_e": np.asarray(p_e),
        "Phi_signal": np.asarray(Phi.signal),
        "Phi_t_list": np.asarray(Phi.t_list),
    }


BASELINES = {
    "qubit_static": _baseline_qubit_static,
    "ramsey_default": _baseline_ramsey,
    "diff_echo_default": _baseline_diff_echo,
    "transient_default": _baseline_transient,
}


def main() -> None:
    print("Generating physics regression baselines ...")
    for name, fn in BASELINES.items():
        print(f"  • {name} ... ", end="", flush=True)
        data = fn()
        path = save_baseline(name, data)
        print(f"saved to {path.relative_to(path.parent.parent.parent)}")
    print("Done.")


if __name__ == "__main__":
    main()
```

**说明**:
- 这个脚本**只在 P0 跑一次**(以及 Track B 物理算法变化时)。
- 输出的 pickle 文件包含 Protocal case 1/2/4 的全部结果数组。
- 如果 case 5 (Cryoscope) 当前可运行,**也加进 BASELINES**,具体由执行者根据 case 5 当前是否能跑通决定;若 evolve case 5 报错,跳过它(后续 P3 内化时再生成)。

### 3.6 baseline 目录说明文件

**新建文件**:`tests/baselines/README.md`

```markdown
# Physics Regression Baselines

This directory holds pickled outputs from `src/protocal.Protocal.evolve()`
at fixed parameters. They define the **frozen physics behavior** the refactor
must preserve.

## Files

| File | Source | Generated by |
|---|---|---|
| `qubit_static.pkl` | TransmonQubit basic properties | `_baseline_qubit_static` |
| `ramsey_default.pkl` | Protocal(type=1).evolve() | `_baseline_ramsey` |
| `diff_echo_default.pkl` | Protocal(type=2).evolve() | `_baseline_diff_echo` |
| `transient_default.pkl` | Protocal(type=4).evolve() | `_baseline_transient` |

## Regenerating

Only regenerate when **intentionally** changing physics. Always justify
in commit message and bump the version field in baseline.

```bash
python -m tests.regression.generate_baselines
```

## Tolerance

Regression tests use `rtol=1e-6, atol=1e-9` (set in `tests/conftest.py`).
```

### 3.7 物理回归测试

**新建文件**:`tests/regression/test_physics_baseline.py`

```python
"""Physics regression: refactored code must reproduce frozen baselines.

These tests are SLOW (each runs the full Protocal.evolve). Mark with
@pytest.mark.regression so they can be skipped during fast iterations:

    pytest -m "not regression"   # skip
    pytest -m regression         # only regression
"""
from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import (
    assert_array_close,
    load_baseline,
)


pytestmark = pytest.mark.regression


def test_qubit_static_properties(qubit_default):
    bl = load_baseline("qubit_static")
    assert qubit_default.frequency == pytest.approx(
        bl["frequency"], rel=1e-12
    )
    assert qubit_default.anharmonicity == pytest.approx(
        bl["anharmonicity"], rel=1e-12
    )
    assert qubit_default.frequency_sensitivity(0.0) == pytest.approx(
        bl["sensitivity_at_zero"], rel=1e-9
    )


def test_ramsey_default_baseline():
    """Protocal(type=1).evolve must match frozen Ramsey output."""
    from src.qubit import TransmonQubit
    from src.protocal import Protocal
    
    q = TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, state=0, n_levels=3,
    )
    proto = Protocal(type=1)
    proto.initialize(q, state=0)
    Phi, tau_list, p_e_list = proto.evolve(q)
    
    bl = load_baseline("ramsey_default")
    assert_array_close(Phi.signal, bl["Phi_signal"], name="Phi.signal")
    assert_array_close(np.asarray(tau_list), bl["tau_list"], name="tau_list")
    assert_array_close(np.asarray(p_e_list), bl["p_e_list"], name="p_e_list")


def test_diff_echo_default_baseline():
    from src.qubit import TransmonQubit
    from src.protocal import Protocal
    
    q = TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, state=0, n_levels=3,
    )
    proto = Protocal(type=2)
    proto.initialize(q, state=0)
    Phi, tau_list, p_e_list, k, t_int = proto.evolve(q)
    
    bl = load_baseline("diff_echo_default")
    assert k == bl["k"]
    assert t_int == pytest.approx(bl["t_int"], rel=1e-12)
    assert_array_close(np.asarray(p_e_list), bl["p_e_list"], name="p_e_list")


def test_transient_default_baseline():
    from src.qubit import TransmonQubit
    from src.protocal import Protocal
    
    q = TransmonQubit(
        EC=2*np.pi*0.2, EJ=2*np.pi*15, T1=10000, T2=8000,
        flux=0.0, state=0, n_levels=3,
    )
    proto = Protocal(type=4)
    proto.initialize(q, state=0)
    t_samples, kernel, scan_list, delta_p, p_e, Phi, ctrl = proto.evolve(q)
    
    bl = load_baseline("transient_default")
    assert_array_close(kernel, bl["kernel"], name="kernel")
    assert_array_close(np.asarray(delta_p), bl["delta_p"], name="delta_p")
    assert_array_close(np.asarray(p_e), bl["p_e"], name="p_e")
```

### 3.8 一个最小 unit 测试(可作为模板)

**新建文件**:`tests/unit/test_smoke.py`

```python
"""Smoke tests: verify all src/ modules import cleanly."""
from __future__ import annotations

import pytest


pytestmark = pytest.mark.unit


def test_import_qubit():
    from src.qubit import TransmonQubit, Cavity, Coupled_System
    from src.qubit import ideal_iSWAP, ideal_CZ
    assert TransmonQubit is not None


def test_import_signal():
    from src.signal import Signal, CompositeSignal
    assert Signal is not None


def test_import_pulse():
    from src.pulse import (
        Pulse, CompositePulse,
        create_pulse, create_ramsey_pulse, create_diff_echo_pulse,
        create_echo_pulse, create_cpmg_pulse, create_cryoscope_pulse,
    )
    assert Pulse is not None


def test_import_protocal():
    from src.protocal import Protocal, Calibration, IQ_readout
    assert Protocal is not None


def test_import_analysis():
    from src.analysis import (
        Analysis,
        generate_basis_functions, basis_function_decomposition, R,
        forward_simulation, compute_jacobian,
        compute_jacobian_finite_difference, levenberg_marquardt,
    )
    assert Analysis is not None


def test_qubit_construction(qubit_default):
    """qubit_default fixture must construct without error."""
    assert qubit_default.frequency > 0
    assert qubit_default.anharmonicity < 0
    assert qubit_default.n_levels == 3
```

### 3.9 (可选) git pre-commit 钩子

**新建文件**:`.git/hooks/pre-commit`(本地,不入 git;在 README.md 提示如何安装)

```bash
#!/usr/bin/env bash
# Run fast unit tests before each commit.
set -e
echo "Running fast unit tests ..."
python -m pytest tests/unit -m "not slow" -x -q || {
    echo "Unit tests failed. Commit aborted."
    exit 1
}
echo "OK."
```

不强制要求,但建议在团队约定中提倡。

---

## 4. 验收标准

每条都必须满足。**任何一条失败,本 phase 不算完成**。

| # | 验收项 | 命令 / 检查 | 期望结果 |
|---|---|---|---|
| 4.1 | 依赖安装 | `pip install -r requirements.txt` | 无错误 |
| 4.2 | pytest 可发现测试 | `pytest --collect-only` | 列出至少 6 个测试用例 |
| 4.3 | smoke 测试通过 | `pytest tests/unit -v` | 全部 pass |
| 4.4 | baseline 已生成 | `ls tests/baselines/*.pkl` | 至少 4 个 pkl 文件(qubit_static, ramsey, diff_echo, transient) |
| 4.5 | 回归测试通过 | `pytest tests/regression -m regression -v` | 全部 pass |
| 4.6 | 跑 baseline 重复一致 | 重新运行 `python -m tests.regression.generate_baselines` 后再跑回归测试 | 全部 pass |
| 4.7 | Notebook 可运行 | 打开 [Simulation.ipynb](../../Simulation.ipynb),Restart & Run All 至少能跑前 5 个 cell | 无 import 错误 |
| 4.8 | Web demo 启动 | `python web_demo.py` | Gradio 界面在浏览器加载,Ramsey 协议跑通 |

---

## 5. 测试要求

本 phase 是建测试基础设施,所以"测试要求"= 上节的验收标准。无需新写额外测试。

---

## 6. 风险与回滚

| # | 风险 | 缓解 |
|---|---|---|
| 6.1 | Protocal.evolve case 1/2/4 当前有 bug,baseline 锚定的是 bug | P0 跑通的结果即定义为"现状",后续若发现 bug,在专门 phase 修正并重生成 baseline,commit message 必须显式标注 |
| 6.2 | pickle 跨 Python 版本不兼容 | requirements.txt 钉死 Python 版本下限(>= 3.10);protocol=5 是 Python 3.8+ |
| 6.3 | qutip 版本升级使 mesolve 数值结果偏移 | requirements.txt 钉死 qutip 版本范围 `>=5.0,<6.0`(写死后绝不破例升级) |
| 6.4 | baseline pickle 体积过大 | 检查 `du -sh tests/baselines/`,若 > 1 MB,改用 `npz` 格式或考虑 git-lfs |

**回滚策略**:本 phase 全部新增,无修改,直接 `git rm -r tests/ pytest.ini` 即可回滚。

---

## 7. 输出物

### 7.1 文件清单

```
tests/
├── __init__.py
├── conftest.py
├── baselines/
│   ├── README.md
│   ├── qubit_static.pkl
│   ├── ramsey_default.pkl
│   ├── diff_echo_default.pkl
│   └── transient_default.pkl
├── unit/
│   ├── __init__.py
│   └── test_smoke.py
├── integration/
│   └── __init__.py
└── regression/
    ├── __init__.py
    ├── generate_baselines.py
    └── test_physics_baseline.py

pytest.ini                              (新建)
requirements.txt                        (修改:钉死版本)
```

### 7.2 接口快照(供后续 phase 引用)

| 符号 | 来源 |
|---|---|
| `assert_array_close(actual, expected, rtol, atol, name)` | tests/conftest.py |
| `load_baseline(name) -> dict` | tests/conftest.py |
| `save_baseline(name, data) -> Path` | tests/conftest.py |
| `qubit_default` (pytest fixture) | tests/conftest.py |
| `t_list_default` (pytest fixture) | tests/conftest.py |
| `rng` (pytest fixture) | tests/conftest.py |
| `RTOL_PHYSICS = 1e-6` | tests/conftest.py |
| `ATOL_PHYSICS = 1e-9` | tests/conftest.py |

### 7.3 后续 phase 的钩子

- P1+ 每个新 module 必须新增对应 unit test。
- P2+ 任何 case 内化都必须新增对应 integration test,验证内化后等价于 baseline。
- P3+ LM/Cryoscope 内化前,先在 src/ 完成,**重生成相应 baseline**(需 commit 显式说明)。
- P4+ 失真模型加入后,**新增** baseline(predistortion_default.pkl 等),不影响已有 baseline。
- P5+ Z-crosstalk 加入后,**新增** baseline。

---

## 8. 与 Track B 的协调点

| Track B 任务 | 影响 |
|---|---|
| _TODO_master.md 0.1 (case 1 死代码清理) | 不需要重生成 baseline,因为死代码不影响输出 |
| _TODO_master.md 0.2 (kernel 自动校准) | **需要重生成 transient_default.pkl** (kernel 数值会变);commit message 显式标注 |
| _TODO_master.md 0.3 (LM 修复) | LM 还未进入 baseline,无影响;P3 引入时再加 |
| _TODO_master.md 1.1 (Cryoscope) | 完成后,在生成脚本里加 case 5/6/7 baseline;P3 内化前必须存在 |
| _TODO_master.md 1.2 (瞬态标定) | 完成后,在生成脚本里加 case 8 baseline |
| _TODO_master.md 1.3 (失真模型) | 完成后,在生成脚本里加 distortion baseline;P4 内化前必须存在 |

---

## 9. 完成确认清单

执行者完成本 phase 时,在 commit message 引用本清单逐项打勾:

```
- [ ] requirements.txt 已钉死版本
- [ ] pytest.ini 已创建
- [ ] tests/ 目录骨架完整(8 个 .py + 1 个 README.md + 4 个 .pkl)
- [ ] tests/conftest.py 含 qubit_default/load_baseline/assert_array_close
- [ ] tests/regression/generate_baselines.py 跑通,产出 4 个 pkl
- [ ] pytest tests/unit -v 全部 pass
- [ ] pytest tests/regression -m regression -v 全部 pass
- [ ] Simulation.ipynb 前 5 个 cell 可运行
- [ ] python web_demo.py 启动并能跑 Ramsey
```

---

**下一步**:[phase_1_handbook.md](phase_1_handbook.md) — 建立 sqc/ 骨架与镜像。
