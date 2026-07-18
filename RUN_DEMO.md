


# 量子传感仿真平台Web演示 - 运行说明

> **v1(sqc)当前演示 = [`web_demo_v2.py`](web_demo_v2.py)。** 下方"运行演示"及之后的章节描述的是**旧版 [`web_demo.py`](web_demo.py)(基于 `src/`)**,保留作历史参考。新用户请使用下述 sqc 版。

## 🚀 当前版本:web_demo_v2.py(sqc,推荐)

基于 v1 的 `sqc` 包,涵盖器件配置、传感协议、波形重建、预失真等标签页。

```bash
# 安装(二选一)
pip install -e ".[demo]"            # 推荐:随包安装 gradio
pip install -r requirements_demo.txt

# 运行(项目根目录下)
python web_demo_v2.py
```

浏览器打开终端提示的地址(默认 http://localhost:7860)。

- 支持的协议:Rabi / Ramsey / 差分回波 / 瞬态波形重建 / Cryoscope,以及预失真设计与验证。
- **实验性标签页(如 Z-Crosstalk)在 v1 中默认隐藏**:由 [`web_demo_v2.py`](web_demo_v2.py) 顶部的 `SHOW_EXPERIMENTAL = False` 控制,翻为 `True` 即可恢复。
- 教程见 [`Simulation_sqc.ipynb`](Simulation_sqc.ipynb);模块参考见 [`docs/architecture.md`](docs/architecture.md)。

---

## (以下为旧版 web_demo.py / `src` 版,历史参考)

## 概述

本Web演示基于现有的量子传感仿真平台代码，提供了一个交互式界面，允许用户调整Qubit参数、选择传感协议并查看仿真结果。

**重要**: 本演示不修改任何原代码，仅使用现有的功能。某些功能可能受限，详细问题请查看TODO.md。

## 系统要求

- Python 3.8 或更高版本
- 至少 4GB 内存
- 网络浏览器（Chrome、Firefox、Edge等）

## 安装步骤

### 方法1: 使用现有环境（推荐）

如果您已经安装了原项目的依赖：

```bash
# 进入项目目录
cd "Sensing project"

# 安装Web演示额外需要的Gradio
pip install gradio>=4.0.0
```

### 方法2: 全新安装

```bash
# 1. 进入项目目录
cd "Sensing project"

# 2. 安装所有依赖（包括原项目依赖和Web演示依赖）
pip install -r requirements.txt  # 原项目依赖
pip install gradio>=4.0.0        # Web界面

# 或者使用演示专用的requirements文件
pip install -r requirements_demo.txt
```

### 方法3: 使用conda环境

```bash
# 创建新环境
conda create -n sensing-demo python=3.9
conda activate sensing-demo

# 安装依赖
pip install -r requirements_demo.txt
```

## 运行演示

### 基本运行

```bash
# 确保在项目根目录下
cd "Sensing project"

# 运行Web演示
python web_demo.py
```

### 可选参数

目前演示脚本没有命令行参数，但您可以在代码中修改：
- 服务器端口：默认为7860，如需修改请编辑web_demo.py第480行
- 服务器地址：默认为0.0.0.0（所有接口）

### 验证安装

运行前可以验证环境：

```bash
# 测试核心依赖
python -c "import qutip, numpy, matplotlib; print('依赖检查通过')"

# 测试项目模块导入
python -c "import sys; sys.path.insert(0, '.'); from src.qubit import TransmonQubit; print('模块导入成功')"
```

## 使用演示

1. **启动服务器**：运行 `python web_demo.py`
2. **打开浏览器**：访问 http://localhost:7860
3. **调整参数**：
   - Qubit参数：EC、EJ、T1、T2、磁通等
   - 协议参数：选择协议类型及相关参数
4. **运行仿真**：点击"运行仿真"按钮
5. **查看结果**：右侧显示图表和详细信息

## 支持的协议

| 协议 | 名称 | 状态 | 说明 |
|------|------|------|------|
| 0 | 拉比振荡测量 | ✅ 已实现 | 测量共振驱动下的Rabi振荡 |
| 1 | Ramsey干涉测量 | ✅ 已实现 | 测量两个π/2脉冲间的干涉条纹 |
| 4 | 瞬态磁场测量 | ✅ 已实现 | 使用滑动测量和维纳反卷积重建磁场 |

## 已知问题

1. **协议2和3未实现**：自旋回波和CPMG协议在原代码中只有框架
2. **计算时间**：复杂计算可能需要几秒到几十秒
3. **参数限制**：某些参数组合可能导致数值问题
4. **浏览器兼容性**：建议使用Chrome或Firefox

## 故障排除

### 1. 导入错误

```
ModuleNotFoundError: No module named 'qutip'
```

**解决方案**：
```bash
pip install qutip
```

### 2. 端口被占用

```
OSError: [Errno 98] Address already in use
```

**解决方案**：
- 关闭占用7860端口的程序
- 修改web_demo.py中的server_port参数
- 使用命令 `lsof -i :7860` (Linux/Mac) 或 `netstat -ano | findstr :7860` (Windows) 查看占用情况

### 3. 模块导入失败

```
ImportError: cannot import name 'Protocal' from 'src.protocal'
```

**解决方案**：
- 确认在项目根目录下运行
- 检查src目录下是否有protocal.py文件
- 注意拼写：原代码使用"protocal"而不是"protocol"

### 4. 计算错误

```
RuntimeError: 数值计算失败
```

**解决方案**：
- 检查参数是否在合理范围内
- 参考Simulation.ipynb中的示例参数
- 减小仿真时间或采样点数

## 开发说明

### 文件结构

```
Sensing project/
├── src/                    # 原项目源代码
│   ├── qubit.py           # Qubit类定义
│   ├── signal.py          # 信号类定义
│   ├── pulse.py           # 脉冲类定义
│   ├── protocal.py        # 协议类定义（注意拼写）
│   └── analysis.py        # 分析工具类
├── web_demo.py            # Web演示主程序
├── TODO.md               # 代码改进建议
├── requirements_demo.txt  # 演示依赖
├── requirements.txt       # 原项目依赖
└── Simulation.ipynb      # 原项目示例笔记本
```

### 扩展演示

要添加新功能：

1. **添加新协议**：先在原代码中实现，然后在web_demo.py的run_simulation函数中添加处理逻辑
2. **添加新参数**：在Gradio界面中添加控件，更新update_ui和run_simulation函数
3. **改进可视化**：修改结果图表的生成代码

### 代码修改原则

- **不修改原代码**：所有修改应集中在web_demo.py中
- **向后兼容**：确保原项目的其他部分仍能正常工作
- **错误处理**：对可能失败的操作添加适当的异常处理

## 联系与支持

如有问题：

1. 查看控制台输出获取详细错误信息
2. 参考Simulation.ipynb中的示例
3. 检查TODO.md了解代码限制
4. 确保使用正确的Python环境和依赖版本

## 许可证

本演示基于原项目代码，遵循原项目的许可证。

---
*最后更新: 2026-03-07*
*量子传感仿真平台 Web演示*