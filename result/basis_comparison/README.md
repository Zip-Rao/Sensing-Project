# 基函数对比分析和绘图

本目录包含用于生成和绘制基函数对比分析数据的脚本，符合科研制图规范。

## 脚本说明

### 1. 数据生成脚本: `generate_basis_comparison_data.py`

**功能**：
- 运行B-spline、Fourier和Legendre三种基函数的数值反演
- 收集收敛历史数据（残差 vs 迭代次数）
- 计算最终RMSE（重建信号与原始信号的均方根误差）
- 保存结果为JSON和NPZ格式

**使用方法**：
```bash
# 基本用法（使用默认参数）
python generate_basis_comparison_data.py

# 指定输出目录
python generate_basis_comparison_data.py --output-dir result/basis_comparison

# 快速测试模式（减少参数以加快运行）
python generate_basis_comparison_data.py --quick-test

# 使用自定义配置文件
python generate_basis_comparison_data.py --config my_config.json
```

**输出文件**：
- `basis_comparison_data.json`: 包含完整元数据和结果的JSON文件
- `basis_comparison_data.npz`: 包含数值数据的压缩NPZ文件

### 2. 绘图脚本: `plot_basis_comparison.py`

**功能**：
- 绘制基函数对比综合图（包含收敛曲线、RMSE对比、信号重建）
- 绘制收敛历史细节图（每种基函数的LM算法收敛行为）
- 符合科研制图规范（高DPI、矢量图输出、合适字体大小）

**使用方法**：
```bash
# 基本用法（使用默认数据文件）
python plot_basis_comparison.py

# 指定数据文件和输出路径
python plot_basis_comparison.py data.npz --output figure.png

# 同时生成收敛历史图
python plot_basis_comparison.py --convergence-output convergence.png

# 使用LaTeX渲染文本（需要系统安装LaTeX）
python plot_basis_comparison.py --latex

# 只显示图形而不保存
python plot_basis_comparison.py --show
```

**输出图形**：
1. **综合对比图** (`basis_comparison_figure.png`):
   - (a) 收敛曲线对比：三种基函数的残差范数 vs 迭代次数
   - (b) RMSE对比：柱状图显示三种基函数的最终重建误差
   - (c) 信号重建对比：原始信号与三种基函数重建信号的对比

2. **收敛历史图** (`convergence_history.png`):
   - 三种基函数各自的LM算法收敛行为细节

## 参数配置

### 默认参数
默认参数在`generate_basis_comparison_data.py`的`set_default_parameters()`函数中设置，包括：

1. **量子比特参数**:
   - `EC=0.2 GHz`, `EJ=10.0 GHz`
   - `T1=100 μs`, `T2=50 μs`
   - 两能级系统，初始态|0⟩

2. **磁场信号参数**:
   - 瞬态信号（类型4）：上升10ns，下降10ns，中心100ns
   - 幅度0.1，时间范围0-200ns，400个点

3. **控制脉冲参数**:
   - Ramsey脉冲序列，时长20ns，40个点

4. **数值反演参数**:
   - 基函数类型：`['bspline', 'fourier', 'legendre']`
   - 基函数数量：20
   - 正则化参数：100.0
   - 最大迭代次数：30
   - 收敛容忍度：1e-6

### 自定义配置
可以通过JSON配置文件自定义参数：
```json
{
  "qubit": {
    "EC": 0.2,
    "EJ": 10.0,
    "T1": 100000.0,
    "T2": 50000.0,
    "flux": 0.0,
    "state": 0,
    "n_levels": 2
  },
  "magnetic_signal": {
    "type": 4,
    "t_min": 0,
    "t_max": 200,
    "n_points": 400,
    "amplitude": 0.1,
    "rise": 10,
    "fall": 10,
    "center": 100,
    "noise_level": 0.0001
  },
  "inversion": {
    "basis_types": ["bspline", "fourier", "legendre"],
    "n_basis": 20,
    "lambdas": 100.0,
    "max_iter": 30,
    "tol": 1e-6
  }
}
```

## 科研制图规范

绘图脚本遵循以下科研制图规范：

1. **字体和字号**:
   - 轴标签: 11pt
   - 标题: 12pt
   - 刻度标签: 10pt
   - 图例: 10pt
   - 使用serif字体（Times New Roman/STIX）

2. **图形质量**:
   - DPI: 300
   - 矢量图输出（PDF格式）
   - 紧凑布局，适当边距

3. **颜色方案**:
   - B-spline: 红色 (#e41a1c)
   - Fourier: 蓝色 (#377eb8)
   - Legendre: 绿色 (#4daf4a)
   - 原始信号: 黑色

4. **标记和线型**:
   - 清晰可区分的标记（圆形、方形、三角形）
   - 适当的线宽（2.0）和标记大小（8）

5. **多子图布局**:
   - 清晰标注子图标签（(a), (b), (c)）
   - 一致的坐标轴范围和网格

## 注意事项

1. **计算时间**:
   - 完整数据生成可能需要较长时间（取决于参数设置）
   - 使用`--quick-test`参数进行快速测试
   - 可减少`n_points`、`n_basis`和`max_iter`以加快运行

2. **依赖项**:
   - 需要安装`numpy`, `scipy`, `matplotlib`
   - 可选：`qutip`（量子仿真）
   - LaTeX渲染需要系统安装LaTeX

3. **文件路径**:
   - 确保输出目录有写入权限
   - 使用绝对路径或相对于脚本位置的路径

## 示例工作流程

```bash
# 1. 生成数据（快速测试模式）
python generate_basis_comparison_data.py --quick-test --output-dir result/basis_comparison_test

# 2. 查看生成的数据摘要
python plot_basis_comparison.py result/basis_comparison_test/basis_comparison_data.npz --summary

# 3. 绘制图形
python plot_basis_comparison.py result/basis_comparison_test/basis_comparison_data.npz \
  --output result/basis_comparison_test/comparison.png \
  --convergence-output result/basis_comparison_test/convergence.png

# 4. 使用LaTeX渲染（如果需要出版物质量）
python plot_basis_comparison.py result/basis_comparison_test/basis_comparison_data.npz --latex
```

## 故障排除

1. **导入错误**:
   - 确保`src`目录在Python路径中
   - 检查所有依赖包是否已安装

2. **内存不足**:
   - 减少`n_points`和`n_basis`参数
   - 使用较小的信号范围

3. **收敛问题**:
   - 调整正则化参数`lambdas`
   - 增加`max_iter`或调整`tol`
   - 检查初始猜测信号是否合理

4. **绘图问题**:
   - 确保`matplotlib`版本兼容
   - 检查LaTeX安装（如果使用`--latex`选项）

## 输出示例

成功运行后，将生成以下文件：
```
result/basis_comparison/
├── basis_comparison_data.json     # 完整数据（JSON格式）
├── basis_comparison_data.npz      # 数值数据（NPZ格式）
├── basis_comparison_figure.png    # 综合对比图
├── basis_comparison_figure.pdf    # 综合对比图（PDF矢量图）
├── convergence_history.png        # 收敛历史图
└── convergence_history.pdf        # 收敛历史图（PDF矢量图）
```