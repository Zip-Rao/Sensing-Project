# Installation

## Requirements

- Python ≥ 3.10
- [QuTiP](https://qutip.org/) ≥ 5.0 (open-source quantum dynamics library)
- NumPy ≥ 1.24, SciPy ≥ 1.11, Matplotlib ≥ 3.7

## Set up an isolated environment (recommended)

Install into a dedicated virtual environment rather than the system or base
Python. This keeps `sqc` and its dependencies from conflicting with other
projects, and everything can be removed by deleting a single folder. Either
tool works equally well:

```bash
# Option A — venv (Python standard library, zero extra tools)
python -m venv .venv

# activate — pick the line for your shell:
source .venv/bin/activate          # Linux / macOS (bash/zsh)
.venv\Scripts\Activate.ps1         # Windows PowerShell
.venv\Scripts\activate.bat         # Windows cmd
```

`source` only exists on Linux/macOS shells; on Windows use one of the two
`.venv\Scripts\` lines above. Once active, the prompt shows a `(.venv)` prefix.

```bash
# Option B — conda
conda create -n sqc python=3.11  # name it whatever you like
conda activate sqc
```

## Install from source (recommended)

With the environment activated:

```bash
git clone https://github.com/Zip-Rao/Sensing-Project.git
cd Sensing-Project
pip install -e .
```

This installs the core `sqc` package in editable mode. The `src/` legacy
implementation and test suite are not distributed and remain local only.

## Optional extras

```bash
# Interactive web demo (Gradio)
pip install -e ".[demo]"

# Documentation build tools (Sphinx + furo + nbsphinx …)
pip install -e ".[docs]"

# Test suite
pip install -e ".[test]"
```

## Verify

```python
import sqc
print(sqc.__version__)          # 1.0.0

from sqc.devices import TransmonQubit
from sqc.workflows import SensingWorkflow
print("ok")
```
