# 问题清单解决方案

对应 `parallel_optimization_plan.md` 第 5 节的 7 个问题，按优先级逐一给出具体改法。

---

## 🔴 问题 1：`protocal.py:sliding_measurement` — H_total 闭包传给多进程，pickle 失败

### 根因

```python
# sliding_measrement 内部
def H_total(t, t_delay, args):      # ← 嵌套闭包，捕获了 qubit, qubit_t 等
    ...
Parallel(n_jobs=6, backend='multiprocessing')(
    delayed(self.single_measurement)(..., H_total)   # ← H_total 无法 pickle
    for t_delay in scan_list
)
```

`H_total` 是嵌套闭包，`self.single_measurement` 是类方法，两者都无法被 `pickle` 序列化，多进程直接报 `PicklingError`。

### 解决方案

**不传 H_total，改为传 numpy 数组，让 worker 自己重建 H。** 详见 `multiprocessing_refactor.md`。

核心思路：

```python
# 在 sliding_measrement 中：
# 1. 预计算 QuTiP 对象 → numpy（主进程完成）
H_0_mats    = np.array([qt.get_hamiltonian_rwa(...).full() for qt in qubit_t])
H_pulse_mats = np.array([control_pulse.get_hamiltonian(t).full() for t in pulse_t])

# 2. 将 numpy 数组打包成 tuple 传入 worker（tuple 可 pickle）
args = (t_delay, H_0_mats, H_pulse_mats, Phi_t, pulse_t, qubit_state_np, n_levels)

# 3. Worker 是模块顶层函数（不是方法，不是闭包）
def _worker(args):           # 模块顶层，可 pickle
    t_delay, H_0_mats, ... = args
    def H_func(t, _args): ...   # worker 内部定义，不跨进程传递
    return mesolve(H_func, ...)
```

---

## 🔴 问题 2：`protocal.py:single_measurement` — H_total 每 ODE 步创建 3 个 Qobj

### 根因

```python
def H_total(t, args):
    # 以下三行每次 ODE 步被调用 6 次（RK45 的 6 阶段）
    H_0 = qubit_current.get_hamiltonian_rwa(control_pulse.omega_d)  # 创建 Qobj
    H_pulse = control_pulse.get_hamiltonian(t)                       # 创建 Qobj
    return H_0 + H_pulse                                             # 创建 Qobj（加法）
```

N≈420 个扫描点，每次 ~200 ODE 步，每步 6 次调用，每次 3 个 Qobj：**420 × 200 × 6 × 3 ≈ 150 万次 Qobj 创建**，全程持有 GIL。

### 解决方案

**将 H_total 改为 QuTiP 列表格式（数组系数），消除 Python 回调。**

QuTiP 的 `mesolve` 支持如下格式，系数为 numpy 数组时在 C 层插值，不回调 Python：

```python
H = [H_static_qobj, [H_op_1, coeff_array_1], [H_op_2, coeff_array_2], ...]
```

针对当前 RWA 框架下的 2 能级 Transmon，哈密顿量分解如下：

```python
# --- 在 sliding_measrement 中预计算（所有 t_delay 共用）---

# 1. 静态算符（只建一次）
from qutip import num, destroy, create
n_op   = num(qubit.n_levels)                   # 数算符
anh_op = n_op * n_op - n_op                   # 非谐项
a      = destroy(qubit.n_levels)
H_anh  = (qubit.EC / 2) * anh_op             # 静态部分（Qobj）

# 2. 失谐系数数组：Δ(t) = ω(B(t)) - ω_d
#    从已预计算的 qubit_t 中直接读取频率，不需要调用 get_hamiltonian
omega_d = control_pulse.omega_d
delta_arr = np.array([qt.frequency - omega_d for qt in qubit_t])   # shape: (n_phi,)

# 3. 脉冲系数数组（在 t_evolve 网格上预计算）
#    这里以 Ramsey 脉冲为例，实际根据 CompositePulse 具体结构展开
#    pulse_x_arr[j] = Ω(t_j) * cos(φ(t_j)) / 2
#    pulse_y_arr[j] = Ω(t_j) * sin(φ(t_j)) / 2

# 4. 对每个 t_delay 构建列表格式 H（t_delay 改变的只是时间轴偏移）
# 关键：delta_arr 对应 Phi_signal 的时间轴；
#       对于给定 t_delay，需将 Phi 时间轴映射到 t_evolve 网格上

# 5. 调用 mesolve（无 Python 回调）
H_list = [
    H_anh,
    [n_op, [t_evolve, delta_on_t_evolve]],    # 失谐：随磁场变化
    [(a + a.dag()) / 2, [t_evolve, pulse_x]], # 脉冲 X 分量
    [1j * (a - a.dag()) / 2, [t_evolve, pulse_y]], # 脉冲 Y 分量
]
result = mesolve(H_list, psi0, t_evolve, [], e_ops=[psi_e_proj])
```

**效果**：150 万次 Qobj 创建 → 0 次，即使不并行也有 5~10x 加速。

> **注**：`CompositePulse.get_hamiltonian` 目前在 Python 层逐脉冲查找，需要先将脉冲包络展开为 numpy 数组，才能用列表格式。这是改造的主要工作量。

---

## 🔴 问题 3：`analysis.py:forward_simulation` — `Options(store_states=True)` 存储所有中间态

### 根因

```python
options = Options(store_states=True)   # 存储所有时间步的量子态
for t_i in t_list:
    result.append(mesolve(H, qubit.state, t_evolve, [], e_ops=[...], options=options))
```

`t_evolve` 长度 ≈ `n_B + n_pulse ≈ 420`，每个时间步存一个 `Qobj`（dim×1 复数向量）。对于 N=420 个扫描点：

```
内存：420 × 420 个 Qobj ≈ 17 万个对象常驻内存
时间：每个 Qobj 的写入和引用计数管理都在 Python 层
```

`forward_simulation` 在 `levenberg_marquardt` 每轮迭代调用**两次**（正式模拟 + trial 模拟），其中 trial 模拟根本不需要中间态。

### 解决方案

**只在 `compute_jacobian` 需要时才存储状态，其余调用关闭。** 给 `forward_simulation` 增加一个参数：

```python
def forward_simulation(qubit, control_pulse, B_curr, t_list, store_states=False):
    options = Options(store_states=store_states)   # 默认关闭
    ...
```

调用处区分：

```python
# compute_jacobian 内部（需要中间态）
result = forward_simulation(qubit, control_pulse, B_curr, t_lists, store_states=True)

# levenberg_marquardt 内部（只需最终测量值）
result       = forward_simulation(qubit, control_pulse, B_curr, t_meas)          # store_states=False
result_trial = forward_simulation(qubit, control_pulse, B_trial, t_meas)         # store_states=False
```

内存节省：关闭后每次 `mesolve` 只存最终态，内存占用降低 420 倍；速度提升因 Python 对象分配减少也有明显改善。

---

## 🟠 问题 4：`analysis.py:forward_simulation` — H_func 内仍有 `Qobj()` 包装

### 根因

最新版的 `forward_simulation`（经过前次改造）中：

```python
def H_func(t, args):
    idx  = np.clip(np.searchsorted(t_list, t), 0, n_B - 1)
    H_np = H_0_mats[idx]                               # numpy，快
    t_local = t - t_i + half_pulse
    return Qobj(H_np + get_H_pulse_np(t_local), dims=dims)  # ← 每步仍创建 Qobj
```

`Qobj()` 构造本身虽比 `get_hamiltonian_rwa()` 快（省去了物理计算），但仍在 Python 层申请内存、触发 GIL。对于 N×ODE步×6次/步 的调用规模，仍然可观。

### 解决方案

与问题 2 相同：改为 QuTiP 列表格式，`H_func` 完全消失，`Qobj()` 创建次数降为 0。

这是问题 2 和问题 4 的统一解法，不需要分别处理。

---

## 🟡 问题 5：`analysis.py:compute_jacobian` — `adjoint` 内的 Python 回调

### 根因

```python
# 外层预计算（已有）
H_list = [lambda t, args, t_i=t_lists[i]: H_total(t, t_i, args).full() for i in range(N)]

# adjoint 内部（被 solve_ivp 每步调用）
def adjoint(s, mu, *args):
    t_curr = t_M - s
    H = H_list[i](t_curr, args)   # ← 调用 H_total → get_hamiltonian_rwa() → Qobj
    rhs = 1j * (H @ mu - mu @ H)
    ...
```

`H_list[i]` 调用 `H_total`，`H_total` 调用 `get_hamiltonian_rwa()`，每次都创建 Qobj，然后 `.full()` 转 numpy。`solve_ivp` 的 RK45 在 SciPy 的 Fortran 核心中运行，**但每步都要回调 Python 的 `adjoint`**，GIL 被反复获取。

### 解决方案

**在 `for i in range(N)` 循环外预计算整个 `t_evolve` 网格上的 H 矩阵为 3D numpy 数组，`adjoint` 内只做数组索引。**

```python
# --- 在 for i 循环之前，一次性预计算所有时间点的 H numpy 矩阵 ---
# 构建统一的 t_evolve（各 i 对应不同的 t_evolve，但磁场部分共用 H_0_mats）
if control_pulse.frame == 0:
    H_0_np_list = np.array([qt.get_hamiltonian(qubit.frequency).full() for qt in qubit_t])
else:
    H_0_np_list = np.array([qt.get_hamiltonian_rwa(qubit.frequency).full() for qt in qubit_t])
# 脉冲部分
H_pulse_np_list = np.array([control_pulse.get_hamiltonian(t).full() for t in control_pulse.t_list])

# --- 在 for i 循环内，预计算当前 i 对应 t_evolve 网格上的 H 矩阵 ---
n_evolve = len(t_evolve)
H_evolve_mats = np.zeros((n_evolve, dim, dim), dtype=complex)
for j, t in enumerate(t_evolve):
    idx_B = np.clip(np.searchsorted(tB_list, t), 0, len(tB_list) - 1)
    H_evolve_mats[j] = H_0_np_list[idx_B]
    t_local = t - t_i + 0.5 * control_pulse.t_list[-1]
    if 0 <= t_local <= control_pulse.t_list[-1]:
        idx_p = np.clip(np.searchsorted(control_pulse.t_list, t_local), 0, len(control_pulse.t_list) - 1)
        H_evolve_mats[j] += H_pulse_np_list[idx_p]

# --- adjoint 内只做 numpy 索引，无 Python 对象创建 ---
def adjoint(s, mu, *args):
    mu_mat = mu.reshape((dim, dim))
    t_curr = t_M - s
    j = np.clip(np.searchsorted(t_evolve, t_curr), 0, n_evolve - 1)
    H = H_evolve_mats[j]                          # 纯 numpy 数组索引
    rhs = 1j * (H @ mu_mat - mu_mat @ H)          # numpy 矩阵运算
    for k in range(len(c_ops_list)):
        rhs += (c_ops_dag_list[k] @ mu_mat @ c_ops_list[k]
                - 0.5 * (c_ops_dag_list[k] @ c_ops_list[k] @ mu_mat
                         + mu_mat @ c_ops_dag_list[k] @ c_ops_list[k]))
    return rhs.flatten()
```

`adjoint` 内**零 Qobj 创建**，`solve_ivp` 的 Fortran 核在调用 `adjoint` 时只做 numpy 运算，GIL 争抢大幅减少。

---

## 🟡 问题 6：`protocal.py:sliding_measurement` — `n_jobs=-1` 包含 E 核

### 根因

```python
Parallel(n_jobs=6, backend='multiprocessing')(...)  # 当前 n_jobs=6
```

虽然已经写 6，但如果其他地方用了 `-1`，Core Ultra 7 255H 的任务会被分发到 **6P + 8E + 2LP = 16 核**，E 核完成一次 mesolve 的时间是 P 核的 2~3 倍。joblib 批量等待最慢的 worker，整体速度被 E 核拖慢。

### 解决方案

**固定为 6，只使用 P 核：**

```python
Parallel(n_jobs=6, backend='threading')(...)   # 不用 -1
```

如果需要动态获取 P 核数量，可用：

```python
import psutil

def get_p_core_count():
    """获取 P 核数量（Intel 混合架构）"""
    # Intel 混合架构下，P 核的 CPU 亲和性通常是固定的
    # 保守做法：取物理核数的一半（避开 E 核）
    return max(1, psutil.cpu_count(logical=False) // 2)

n_p_cores = get_p_core_count()  # Core Ultra 7 255H 上返回 6
Parallel(n_jobs=n_p_cores, backend='threading')(...)
```

更精确的做法（需要 `py-cpuinfo`）：

```python
# 直接指定 CPU 亲和性，将进程绑定到 P 核
import os
os.sched_setaffinity(0, {0, 1, 2, 3, 4, 5})  # Linux 有效，Windows 需用 psutil
```

---

## 🟢 问题 7：`analysis.py:compute_jacobian` — `G_mat = qubit.n.full()` 在内层循环重复调用

### 根因

```python
for n in range(N_s):                         # 内层循环，N_s ≈ 420
    rho_n = states[i][n].full()
    lam_n = lambda_t[n]
    G_mat = qubit.n.full()                   # ← 每次循环都重新计算，结果完全相同
    trace_values[n] = np.trace(lam_n @ (G_mat @ rho_n - rho_n @ G_mat))
```

`qubit.n` 是 QuTiP 的数算符，`.full()` 将其转为 numpy 矩阵。这个结果对整个 `compute_jacobian` 的执行过程完全不变，却在 `for i in range(N)` × `for n in range(N_s)` 的双层循环内被调用 `N × N_s ≈ 420 × 420 ≈ 17 万次`。

### 解决方案

**提到两层循环之外，只计算一次：**

```python
# --- 在 for i 循环之前 ---
G_mat = qubit.n.full()   # 只算一次

for i in range(N):
    ...
    for n in range(N_s):
        rho_n = states[i][n].full()
        lam_n = lambda_t[n]
        # G_mat 直接引用外层变量，无重复创建
        trace_values[n] = np.trace(lam_n @ (G_mat @ rho_n - rho_n @ G_mat))
```

同时 `trace_values` 的计算可以向量化，消除内层循环：

```python
# 向量化版本（同时消除内层 for n 循环）
G_mat = qubit.n.full()

# rho_stack[n] = states[i][n].full()，构建批量矩阵
rho_stack = np.array([states[i][n].full() for n in range(N_s)])  # (N_s, dim, dim)
lam_stack = lambda_t[:N_s]                                         # (N_s, dim, dim)

# 批量矩阵乘法（numpy 广播）
commutator = G_mat @ rho_stack - rho_stack @ G_mat                # (N_s, dim, dim)
trace_values = np.trace(lam_stack @ commutator, axis1=-2, axis2=-1)  # (N_s,)
```

---

## 优先级与改动量汇总

| 优先级 | 问题 | 改动量 | 预期收益 |
|--------|------|--------|---------|
| 🔴 1 | sliding_measurement 闭包 pickle | 中（参考 multiprocessing_refactor.md） | 从报错→正常运行 |
| 🔴 2 | H_total Python 回调（protocal） | 大（需分解脉冲包络） | 5~10x 加速 |
| 🔴 3 | `store_states=True` 默认开 | 小（加一个参数） | 内存降低 400x |
| 🟠 4 | H_func 仍有 Qobj 包装 | 与问题 2 合并解决 | 与问题 2 合并 |
| 🟡 5 | adjoint Python 回调 | 中（预计算 H_evolve_mats） | 2~5x adjoint 加速 |
| 🟡 6 | n_jobs=-1 含 E 核 | 极小（改一个数字） | 20~30% 提升 |
| 🟢 7 | G_mat 重复调用 | 极小（移出循环） | 微小 |

**建议改动顺序**：3 → 6 → 7（改动极小，立即见效）→ 5 → 1 → 2（改动较大，效果最显著）
