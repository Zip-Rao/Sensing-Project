# 安装

## 环境要求

- Python ≥ 3.10
- [QuTiP](https://qutip.org/) ≥ 5.0(开源量子动力学库)
- NumPy ≥ 1.24、SciPy ≥ 1.11、Matplotlib ≥ 3.7

## 从源码安装(推荐)

```bash
git clone https://github.com/Zip-Rao/Sensing-Project.git
cd Sensing-Project
pip install -e .
```

以可编辑模式安装核心 `sqc` 包。遗留实现 `src/` 与测试套件不随发行分发,仅保留在本地。

## 可选 extras

```bash
# 交互式 web demo(Gradio)
pip install -e ".[demo]"

# 文档构建工具(Sphinx + furo + nbsphinx …)
pip install -e ".[docs]"

# 测试套件
pip install -e ".[test]"
```

## Conda 环境(可选)

任何 Python ≥ 3.10 的环境都可以——venv、conda 或系统 Python 均可。若偏好 conda:

```bash
conda create -n sqc python=3.11   # 名字随意
conda activate sqc
pip install -e ".[demo,test]"
```

## 验证

```python
import sqc
print(sqc.__version__)          # 1.0.0

from sqc.devices import TransmonQubit
from sqc.workflows import SensingWorkflow
print("ok")
```
