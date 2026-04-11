# 并行优化方案：GIL 分析与最优实现路径

## 1. 问题诊断：为什么现有方案比串行还慢，还会报错

### 1.1 当前代码的核心矛盾

`sliding_measurement` 当前使用 `backend='multiprocessing'`，并将 `H_total`（一个嵌套闭包）作为参数传入每个子进程：

```python
def H_total(t, t_delay, args):   # ← 捕获了 qubit_t, qubit, Phi_signal, control_pulse
    ...
p_e = Parallel(n_jobs=6, backend='multiprocessing')(
    delayed(self.single_measurement)(..., H_total) for t_delay in scan_list
)
```

**这在根本上就无法工作**，原因有两层：

| 问题 | 原因 | 现象 |
|------|------|------|
| 闭包不可序列化 | Python 的 `pickle` 无法序列化嵌套闭包（closure），闭包捕获的局部变量也不可 pickle | 直接抛出 `PicklingError` 或 `AttributeError` |
| QuTiP 对象跨进程 | 即使用 `cloudpickle`，QuTiP 的 `Qobj` 包含 C 层内存结构，跨进程 pickle/unpickle 会损坏或产生极大开销 | 结果错误，或序列化时间超过计算时间本身 |

**多进程的开销模型**：
```
总时间 = 序列化(qubit_t + H_total + ...) + 子进程启动 + 计算 + 反序列化结果
       ≈ 数秒序列化开销 + 原始计算时间
```
因此即使计算本身加速了，序列化开销也会让总时间反而更长。

---

### 1.2 GIL 的真实影响路径

要理解为什么 `threading` 方案也慢，需要追踪 GIL 在调用栈中的位置：

```
mesolve(H_func, ...)
  └─ VODE/RK45 (Fortran/C，释放 GIL) ← 这里释放了，但...
       └─ 每个 ODE 步调用 H_func(t, args)  ← 重新获取 GIL
            └─ get_hamiltonian_rwa()         ← 纯 Python + Qobj 构造，持有 GIL
            └─ control_pulse.get_hamiltonian(t) ← Python for 循环，持有 GIL
            └─ H_0 + H_pulse                ← Qobj 加法，持有 GIL
```

**关键数据**（以 `sliding_measurement` 为例，N≈420 个扫描点）：

| 阶段 | 调用次数估算 |
|------|-------------|
| 每次 mesolve 的 ODE 步数 | ~100~300 步（自适应 RK45） |
| 每步 H_func 评估次数 | 6 次（RK45 的 6 阶段） |
| 每次 H_func 创建的 Qobj 数 | 3 个（H_0 + H_pulse + 加法结果） |
| **总 Qobj 创建次数（N=420）** | **420 × 200 × 6 × 3 ≈ 150 万次** |

每次 Qobj 创建都需要在 Python 层申请内存、初始化矩阵、管理引用计数——全程持有 GIL。

**多线程时的实际情况**：6 个线程全都在抢同一把 GIL 来创建 Qobj，效果等同于串行但多了线程切换开销，故**比串行更慢**。

---

### 1.3 混合 CPU 架构的额外影响

Core Ultra 7 255H 的 6 P核 + 8 E核 对 `n_jobs=-1` 产生负面效应：

```
P核完成一次 mesolve:  ~0.5s
E核完成一次 mesolve:  ~1.2s  (频率低 ~40%，IPC 低)
```

由于 joblib 批量分发任务，整批的完成时间取决于**最慢的那个 E 核线程**。即使 P 核已经完成，也要等 E 核收尾，实际加速比远低于理论值。

---

## 2. 根本解决路径：消除 Python 回调

所有并行化方案的共同前提是：**H_total 内部不应该有 Python 对象创建**。否则，无论用多少线程/进程，GIL 都会成为瓶颈。

### 2.1 QuTiP 列表格式（最优方案）

QuTiP 的 `mesolve` 支持将时变 Hamiltonian 表示为**算符 + 系数数组**的列表格式：

```python
H = [
    H_static,               # 静态部分（Qobj）
    [H_op_1, coeff_array_1], # 时变部分：算符 × 标量系数数组
    [H_op_2, coeff_array_2],
    ...
]
```

当 `coeff_array` 是 numpy 数组时，QuTiP 在 C 层做线性插值，**完全绕过 Python，GIL 在整个积分过程中持续释放**。

#### 适用性分析

对于 2 能级 Transmon qubit 在 RWA 框架下，哈密顿量可以精确分解：

| 部分 | 表达式 | 算符 | 时变系数 |
|------|--------|------|---------|
| 失谐项 | `Δ(t) · n` | `n`（数算符） | `Δ(t) = ω(B(t)) - ω_d`，随 B 场变化 |
| 非谐项 | `(α/2)(n²-n)` | `n²-n` | 常数（静态） |
| 脉冲 X | `Ω(t)·cos(φ)/2 · σ_x` | `σ_x/2` | 脉冲包络 × cos(相位) |
| 脉冲 Y | `Ω(t)·sin(φ)/2 · σ_y` | `σ_y/2` | 脉冲包络 × sin(相位) |

这个分解**完全适用于当前 2 能级系统**，不需要修改物理模型。

#### 实现要点

```python
# 1. 预计算静态算符（只算一次）
n_op    = num(qubit.n_levels)          # 数算符
anh_op  = n_op**2 - n_op              # 非谐项算符
sx_half = sigmax() / 2                # or (destroy + create) / 2
sy_half = sigmay() / 2

# 2. 预计算 Δ(t) 系数数组（在 t_evolve 网格上）
#    ω(B(t)) 从预计算的 qubit_t 中读取，不需要每步调用 get_hamiltonian
delta_array = np.array([qubit_t[i].frequency - omega_d
                        for i in range(len(qubit_t))])   # shape: (n_B,)

# 3. 预计算脉冲包络系数数组
pulse_x_array = np.array([pulse_envelope(t)*np.cos(phase(t))
                           for t in t_evolve])
pulse_y_array = np.array([pulse_envelope(t)*np.sin(phase(t))
                           for t in t_evolve])

# 4. 构建 QuTiP 列表格式 H（传给 mesolve）
# QuTiP 格式: [t_list, coeff_array] 打包为 interpolating_array
H = [
    (EC_alpha / 2) * anh_op,                        # 静态
    [n_op,    [t_evolve, delta_array]],              # 失谐（随 B 变化）
    [sx_half, [t_evolve, pulse_x_array]],            # 脉冲 X 分量
    [sy_half, [t_evolve, pulse_y_array]],            # 脉冲 Y 分量
]
mesolve(H, psi0, t_evolve, ...)
```

> **注**：QuTiP 的插值系数格式为 `[t_array, coeff_array]`（QuTiP 5.x），或直接传入系数字符串。具体 API 见 [QuTiP 文档 - Time-Dependent Operators](https://qutip.readthedocs.io/en/latest/guide/dynamics/dynamics-time.html)。

---

### 2.2 为什么这是根本解决方案

```
消除 Python 回调后的 H 求值路径：

mesolve([H_static, [n_op, delta_arr], ...], ...)
  └─ VODE/RK45 (C/Fortran，释放 GIL)
       └─ 插值 delta_arr[i] (C，GIL 已释放) ← 不再回调 Python！
       └─ 矩阵运算 H = H_static + delta * n_op + ... (C)
```

**150 万次 Python Qobj 创建 → 0 次**。这是 5~10 倍加速的来源，与是否并行无关。

---

## 3. 在消除 Python 回调之后的并行化方案

完成第 2 章的改造后，每次 `mesolve` 调用不再有 Python 回调，GIL 在积分过程中持续释放。此时各 `t_i` 间的 `mesolve` 调用**完全独立**，可以真正并行。

### 3.1 方案对比

| 方案 | Backend | 可行性 | 优点 | 缺点 |
|------|---------|--------|------|------|
| **A: threading（消除回调后）** | `threading` | ✅ **推荐** | 零序列化开销；QuTiP C 层积分释放 GIL；共享内存，无需拷贝预计算数组 | E核负载不均衡 |
| **B: loky 多进程（顶层函数）** | `loky` | ✅ 可行 | 完全绕过 GIL；P核独立运行 | 需要将 worker 写为模块级函数；numpy 数组跨进程需 pickle |
| **C: multiprocessing（当前方案）** | `multiprocessing` | ❌ 无效 | — | 闭包不可 pickle；Qobj 跨进程序列化开销巨大 |
| **D: 原生 QuTiP parallel_map** | qutip | ✅ 可行 | 官方支持，处理 Qobj 更稳 | API 较低层，灵活性差 |
| **E: dask** | dask | ⚠️ 过度设计 | 适合集群 | 单机引入不必要的复杂性 |

---

### 3.2 推荐方案：threading（消除回调后）

消除 Python 回调后，`mesolve` 期间 GIL 持续释放，threading 后端可以获得**真实的并行加速**：

```
线程 1 (P核): mesolve for t_i=0   → VODE (C, GIL released) → 完成
线程 2 (P核): mesolve for t_i=1   → VODE (C, GIL released) → 完成
...
线程 6 (P核): mesolve for t_i=5   → VODE (C, GIL released) → 完成
（无 Python 回调，无 GIL 竞争）
```

**关键参数**：

```python
# n_jobs=6：只用 P核，避免 E核拖累
Parallel(n_jobs=6, backend='threading')(...)
```

优势：
- **零序列化开销**：预计算的 numpy 数组（`H_0_mats`, `H_pulse_mats`）直接被所有线程共享内存读取，无需拷贝
- **无 pickle 问题**：函数不需要可序列化
- **最低调度延迟**：线程创建比进程轻量

---

### 3.3 备选方案：loky 多进程（需重构 worker 为顶层函数）

若 threading 仍有 GIL 瓶颈（例如 Qobj 构造还未完全消除），可改用多进程。但**必须将 worker 函数写为模块级函数**（不能是方法或闭包）：

```python
# analysis.py 模块级（不在任何类或函数内部）
def _run_single_mesolve(args):
    '''必须是顶层函数，才能被 pickle'''
    t_i, H_0_mats, H_pulse_mats, t_list, qubit_state, dims, psi_e_proj = args
    # ... 构建 H 并调用 mesolve
    return p_e

# 调用端
from concurrent.futures import ProcessPoolExecutor
with ProcessPoolExecutor(max_workers=6) as ex:
    results = list(ex.map(_run_single_mesolve, args_list))
```

代价：每次调用需要 pickle `H_0_mats`（numpy 数组，约 400×2×2 = 几 KB，可接受）。

---

## 4. 完整优化路线图

```
第一步：消除 H_total 内的 Python 回调（最重要）
  ├─ 将 H 改为 QuTiP 列表格式（数组系数）
  ├─ 预计算 Δ(t) 数组、脉冲 X/Y 系数数组
  └─ 预期收益：5~10x 加速（即使不并行）

第二步：并行化（threading，n_jobs=6）
  ├─ 此时 mesolve 内 GIL 持续释放
  ├─ 6 个 P核 各自独立积分
  └─ 预期收益：额外 4~5x 加速

第三步（可选）：验证与调优
  ├─ 用 n_jobs=[1,2,4,6] 测量，确认超线性加速已消失
  ├─ 检查 E核是否参与（可用 psutil 查看 CPU 亲和性）
  └─ 如仍有问题，切换为 ProcessPoolExecutor
```

**预期总加速**：第一步 × 第二步 ≈ **20~50 倍**（从 6s → 0.1~0.3s per forward_simulation）

---

## 5. 当前代码问题清单（按优先级）

| 优先级 | 位置 | 问题 | 影响 |
|--------|------|------|------|
| 🔴 最高 | `protocal.py:sliding_measurement` | `H_total` 闭包传给多进程，pickle 失败 | 直接报错 |
| 🔴 最高 | `protocal.py:single_measurement` | H_total 每 ODE 步创建 3 个 Qobj | 150 万次创建，主要性能瓶颈 |
| 🔴 最高 | `analysis.py:forward_simulation` | 同上，且 `Options(store_states=True)` 存储所有中间态 | 内存和时间双重开销 |
| 🟠 高 | `analysis.py:forward_simulation` | `backend='multiprocessing'` 下 H_0_mats 预计算（最新版）仍返回 Qobj | 改进有限 |
| 🟡 中 | `analysis.py:compute_jacobian` | `solve_ivp` 的 `adjoint` 也有相同 Python 回调问题 | 影响 LM 迭代速度 |
| 🟡 中 | `protocal.py:sliding_measurement` | `n_jobs=-1` 包含 E核，造成负载不均衡 | 约 20~30% 性能损失 |
| 🟢 低 | `analysis.py:compute_jacobian` | `G_mat = qubit.n.full()` 在内层循环重复调用 | 小优化 |

---

## 6. 参考资料

- [QuTiP 时变 Hamiltonian 文档](https://qutip.readthedocs.io/en/latest/guide/dynamics/dynamics-time.html) — 数组格式的系数插值 API
- [QuTiP parallel_map](https://qutip.readthedocs.io/en/latest/apidoc/qutip.parallel.html) — 官方并行工具
- [Python GIL 与 NumPy/SciPy 释放机制](https://numpy.org/doc/stable/reference/c-api/array.html#threading) — NumPy C API 的 GIL 释放行为
- Intel Thread Director 白皮书 — 混合核心调度机制
