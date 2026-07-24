# Installation

## Requirements

- Python ≥ 3.10
- [QuTiP](https://qutip.org/) ≥ 5.0 (open-source quantum dynamics library)
- NumPy ≥ 1.24, SciPy ≥ 1.11, Matplotlib ≥ 3.7

## Install from source (recommended)

```bash
git clone https://github.com/Serendipity-Zip/Sensing-Project.git
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

## Conda environment (project development)

The project was developed inside a conda environment named `qutip-env`. To
reproduce it:

```bash
conda create -n qutip-env python=3.11
conda activate qutip-env
pip install -e ".[demo,test]"
```

On Windows the full Python path is:
`C:\Users\<user>\anaconda3\envs\qutip-env\python.exe`

## Verify

```python
import sqc
print(sqc.__version__)          # 1.0.0

from sqc.devices import TransmonQubit
from sqc.workflows import SensingWorkflow
print("ok")
```
