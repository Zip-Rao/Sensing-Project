# 安装

## 环境要求

- Python ≥ 3.10
- [QuTiP](https://qutip.org/) ≥ 5.0(开源量子动力学库)
- NumPy ≥ 1.24、SciPy ≥ 1.11、Matplotlib ≥ 3.7

## 创建隔离环境(推荐)

请安装到**专用的虚拟环境**,而非系统或 base Python。这样可以避免 `sqc` 及其依赖与
其他项目冲突,日后也只需删掉一个文件夹就能彻底清理。任选一种工具即可,两者效果相同:

```bash
# 方式 A —— venv(Python 标准库,零额外工具)
python -m venv .venv

# 激活——按你的终端选对应的一行:
source .venv/bin/activate          # Linux / macOS(bash/zsh)
.venv\Scripts\Activate.ps1         # Windows PowerShell
.venv\Scripts\activate.bat         # Windows cmd
```

`source` 只存在于 Linux/macOS 的 shell;Windows 上请用上面两条 `.venv\Scripts\` 中的一条。
激活成功后,命令行提示符前会出现 `(.venv)` 前缀。

```bash
# 方式 B —— conda
conda create -n sqc python=3.11  # 名字随意
conda activate sqc
```

## 从源码安装(推荐)

激活环境后:

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

## 验证

```python
import sqc
print(sqc.__version__)          # 1.0.0

from sqc.devices import TransmonQubit
from sqc.workflows import SensingWorkflow
print("ok")
```
