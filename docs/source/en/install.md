# Installation

## Requirements

- Python ≥ 3.10
- [QuTiP](https://qutip.org/) ≥ 5.0 (open-source quantum dynamics library)
- NumPy ≥ 1.24, SciPy ≥ 1.11, Matplotlib ≥ 3.7

## Install from source (recommended)

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

## Conda environment (optional)

Any Python ≥ 3.10 environment works — venv, conda, or system Python. If you
prefer conda:

```bash
conda create -n sqc python=3.11   # name it whatever you like
conda activate sqc
pip install -e ".[demo,test]"
```

## Verify

```python
import sqc
print(sqc.__version__)          # 1.0.0

from sqc.devices import TransmonQubit
from sqc.workflows import SensingWorkflow
print("ok")
```
