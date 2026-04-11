# `numerical_inverse` 运行时间分析与优化方案

## 1. 调用链与时间构成

```
numerical_inverse (line 172)
  └─ levenberg_marquardt (line 417)  ×max_iter 次（默认50）
       ├─ forward_simulation (line 447)         ← N 次 mesolve
       ├─ compute_jacobian (line 453)
       │    ├─ forward_simulation (line 336)     ← N 次 mesolve（重复！）
       │    └─ for i in range(N):               ← N 次 solve_ivp
       │         └─ adjoint ODE
       └─ forward_simulation (line 461)         ← N 次 mesolve（trial）
```

设 N = 信号长度 + 脉冲长度 - 1（即 `t_meas` 长度），典型值 N ≈ 400~800。

### 每轮迭代的时间估算

| 步骤 | mesolve 次数 | solve_ivp 次数 | 预估时间 |
|------|-------------|---------------|---------|
| forward_simulation (正式) | N | 0 | ~6s |
| compute_jacobian 内部的 forward_simulation | N | 0 | ~6s（重复计算！） |
| compute_jacobian 的 adjoint 循环 | 0 | N | ~6~12s |
| forward_simulation (trial) | N | 0 | ~6s |
| **单轮合计** | **3N** | **N** | **~24s** |

50 轮迭代总时间：**~20 分钟**（理想情况），实际因 Python 回调开销可达 **30~60 分钟**。

---

## 2. 六大瓶颈（按影响排序）

### 瓶颈 1：forward_simulation 被重复调用（最大浪费）

`levenberg_marquardt` 第 447 行调用 `forward_simulation` 得到 `result`，第 453 行 `compute_jacobian` 内部又调用一次完全相同的 `forward_simulation`（第 336 行）。

**每轮浪费 N 次 mesolve，约 6s。**

```python
# levenberg_marquardt line 447-453
result = forward_simulation(qubit, control_pulse, B_curr, t_meas)       # 第 1 次
p_sim = np.array([res.expect[0][-1] for res in result])
J = compute_jacobian(qubit, control_pulse, B_curr, t_meas)              # 内部又调用第 2 次
```

### 瓶颈 2：forward_simulation 使用 Python 回调 H_total

```python
# forward_simulation line 281-290
def H_total(t, t_i, args):
    index = np.searchsorted(t_list, t)
    qubit_current = qubit_t[index]
    H_0_current = qubit_current.get_hamiltonian_rwa(...)  # ← 每 ODE 步创建 Qobj
    return H_0_current + control_pulse.get_hamiltonian_at(...)  # ← 又创建 Qobj
```

N=500 时，每次 forward_simulation 约 500 × 200步 × 6阶段 × 3个Qobj = **180 万次 Qobj 创建**。

### 瓶颈 3：compute_jacobian 的 adjoint 使用 Python 回调

```python
# compute_jacobian line 357-364
def adjoint(s, mu, *args):
    H = H_list[i](t_curr, args)   # ← H_list[i] 调用 H_total → get_hamiltonian_rwa → Qobj
    rhs = 1j * (H @ mu - mu @ H)
    ...
```

H_list 的构建方式（line 345）：
```python
H_list = [lambda t, args, t_i=t_lists[i]: H_total(t, t_i, args).full() for i in range(N)]
```

每次 `H_list[i]` 都要调用 `H_total`，走完 `get_hamiltonian_rwa()` → Qobj 创建 → `.full()` 的完整路径。

### 瓶颈 4：store_states=True 默认开启

```python
# forward_simulation line 279
options = Options(store_states=True)
```

`compute_jacobian` 需要中间态（用于计算 Tr[λ·[G,ρ]]），但 `levenberg_marquardt` 的两次 forward_simulation（正式 + trial）不需要。

每次存储 N × n_evolve 个 Qobj，典型值 500 × 500 = 25万个对象常驻内存。

### 瓶颈 5：compute_jacobian 内 g_values 索引错误（Bug）

```python
# compute_jacobian line 393-399
for s in range(N_s):
    if t_evolve[s] < tB_list[0]:
        g_values[s] = sensitivity[0]
    elif t_evolve[s] >= tB_list[-1]:
        g_values[s] = sensitivity[-1]
    else:
        g_values[i] *= sensitivity[i]   # ← 应该是 g_values[s] *= sensitivity[s]
```

`g_values[i]` 中的 `i` 是外层循环变量（第 i 个测量点），不是内层的 `s`。这是一个 **bug**，会导致 Jacobian 计算错误，进而让 LM 优化方向偏离，需要更多迭代才能收敛（甚至不收敛）。

### 瓶颈 6：G_mat 在双层循环内重复创建

```python
# compute_jacobian line 391
G_mat = qubit.n.full()   # 每次内层循环都重建，但结果完全相同
```

N × N_s ≈ 500 × 500 = 25 万次重复创建。

---

## 3. 优化方案

### 优化 1：消除 forward_simulation 重复调用

将 `compute_jacobian` 的接口改为接受已有的 forward_simulation 结果：

```python
def compute_jacobian(qubit, control_pulse, B_curr, t_lists, fwd_result=None):
    if fwd_result is None:
        fwd_result = forward_simulation(qubit, control_pulse, B_curr, t_lists, store_states=True)
    ...
```

`levenberg_marquardt` 中：

```python
# 只调用一次 forward_simulation（带中间态），同时供 p_sim 提取和 Jacobian 计算
result = forward_simulation(qubit, control_pulse, B_curr, t_meas, store_states=True)
p_sim = np.array([res.expect[0][-1] for res in result])
J = compute_jacobian(qubit, control_pulse, B_curr, t_meas, fwd_result=result)
```

**收益：每轮省掉 N 次 mesolve，约 6s → 总运行时间减少 ~25%。**

### 优化 2：forward_simulation 转 QobjEvo 格式

参照 `protocal.py:sliding_measurement` 已完成的改造模式，将 H_total 从 Python 回调改为 QobjEvo 列表格式：

```python
def forward_simulation(qubit, control_pulse, B_curr, t_list, store_states=False):
    B = B_curr
    qubit_t = qubit.qubit_under_mag(B)
    n_levels = qubit.n_levels
    n_op = qubit.n
    psi_e = basis(n_levels, 1)

    # 预计算静态项
    H_anh = qubit.anharmonicity * 0.5 * (n_op * n_op - n_op)

    result = []
    for t_i in t_list:
        # 构建演化时间轴
        t_m = min(B.t_list[0], t_i - 0.5 * control_pulse.t_list[-1])
        t_M = max(B.t_list[-1], t_i + 0.5 * control_pulse.t_list[-1])
        N_evolve = len(B.t_list) + len(control_pulse.t_list)
        t_evolve = np.linspace(t_m, t_M, N_evolve)

        # 构建频率系数数组（随 B 场变化）
        freq_coeffs = np.zeros(N_evolve)
        for j, t in enumerate(t_evolve):
            if B.t_list[0] <= t <= B.t_list[-1]:
                idx = np.clip(np.searchsorted(B.t_list, t), 0, len(qubit_t) - 1)
                qt = qubit_t[idx]
            else:
                qt = qubit
            if control_pulse.frame == 0:
                freq_coeffs[j] = qt.frequency
            else:
                freq_coeffs[j] = qt.frequency - control_pulse.omega_d

        # 构建 H_list
        H_list = [H_anh]
        if control_pulse.frame == 0:
            H_list.append([n_op + 0.5 * qeye(n_levels), freq_coeffs])
        else:
            H_list.append([n_op, freq_coeffs])

        # 脉冲项：将 control_pulse 的哈密顿量映射到 t_evolve 上
        t_pulse = control_pulse.t_list
        half_dur = 0.5 * t_pulse[-1]
        for op, coeffs in control_pulse.hamiltonian:
            coeff_global = np.zeros(N_evolve, dtype=complex)
            for j, t in enumerate(t_evolve):
                t_loc = t - t_i + half_dur
                if 0 <= t_loc <= t_pulse[-1]:
                    idx = np.clip(np.searchsorted(t_pulse, t_loc), 0, len(coeffs) - 1)
                    coeff_global[j] = coeffs[idx]
            H_list.append([op, coeff_global])

        H_total = QobjEvo(H_list, tlist=t_evolve, order=1)
        options = {"store_states": store_states}
        res = mesolve(H_total, qubit.state, t_evolve, [],
                      e_ops=[psi_e * psi_e.dag()], options=options)
        result.append(res)

    return result
```

**收益：消除 ~180 万次/调用 的 Qobj 创建，单次 forward_simulation 从 ~6s 降至 ~1s（5~6x 加速）。**

### 优化 3：compute_jacobian 的 H 矩阵预计算

将 `adjoint` 内的 Python 回调替换为预计算的 3D numpy 数组索引：

```python
# 预计算所有 qubit_t 对应的 H_0 numpy 矩阵（循环外，一次性完成）
if control_pulse.frame == 0:
    H_0_np = np.array([qt.get_hamiltonian(qubit.frequency).full() for qt in qubit_t])
else:
    H_0_np = np.array([qt.get_hamiltonian_rwa(qubit.frequency).full() for qt in qubit_t])

# 预计算脉冲 numpy 矩阵
t_pulse = control_pulse.t_list
H_pulse_np = np.array([control_pulse.get_hamiltonian_at(t).full() for t in t_pulse])

G_mat = qubit.n.full()  # 只算一次

for i in range(N):
    # 预计算当前 i 对应的 H_evolve numpy 矩阵
    H_evolve = np.zeros((n_evolve, dim, dim), dtype=complex)
    for j, t in enumerate(t_evolve):
        idx_B = np.clip(np.searchsorted(tB_list, t), 0, len(H_0_np) - 1)
        H_evolve[j] = H_0_np[idx_B]
        t_loc = t - t_i + 0.5 * t_pulse[-1]
        if 0 <= t_loc <= t_pulse[-1]:
            idx_p = np.clip(np.searchsorted(t_pulse, t_loc), 0, len(H_pulse_np) - 1)
            H_evolve[j] += H_pulse_np[idx_p]

    # adjoint 内只做 numpy 索引，零 Qobj 创建
    def adjoint(s, mu, *args):
        mu_mat = mu.reshape((dim, dim))
        t_curr = t_M - s
        j = np.clip(np.searchsorted(t_evolve, t_curr), 0, n_evolve - 1)
        H = H_evolve[j]                          # 纯数组索引
        rhs = 1j * (H @ mu_mat - mu_mat @ H)
        for k in range(len(c_ops_list)):
            rhs += (c_ops_dag_list[k] @ mu_mat @ c_ops_list[k]
                    - 0.5 * (c_ops_dag_list[k] @ c_ops_list[k] @ mu_mat
                             + mu_mat @ c_ops_dag_list[k] @ c_ops_list[k]))
        return rhs.flatten()
```

**收益：adjoint 内零 Qobj 创建，solve_ivp 加速 2~5x。**

### 优化 4：修复 g_values 索引 bug

```python
# 修复前（line 399）
g_values[i] *= sensitivity[i]   # bug: i 是外层循环变量

# 修复后
g_values[s] *= sensitivity[s]   # s 是内层循环变量
```

**收益：Jacobian 计算正确 → LM 收敛更快 → 迭代次数可能减半。**

### 优化 5：store_states 按需开启

```python
def forward_simulation(qubit, control_pulse, B_curr, t_list, store_states=False):
    options = {"store_states": store_states}
    ...
```

调用处区分：
```python
# compute_jacobian 内（需要中间态）
result = forward_simulation(..., store_states=True)

# levenberg_marquardt 内 trial 模拟（不需要中间态）
result_trial = forward_simulation(..., store_states=False)
```

**收益：trial 模拟内存降低 ~500x，速度因减少 Python 对象分配也有提升。**

### 优化 6：trace_values 向量化 + G_mat 外提

```python
G_mat = qubit.n.full()  # 外提到所有循环之前

# 向量化计算 trace_values
rho_stack = np.array([states[i][n].full() for n in range(N_s)])  # (N_s, dim, dim)
lam_stack = lambda_t[:N_s]
commutator = G_mat @ rho_stack - rho_stack @ G_mat
trace_values = np.trace(lam_stack @ commutator, axis1=-2, axis2=-1)
```

**收益：消除内层 for 循环，25 万次 G_mat 创建降为 0 次。**

---

## 4. 综合收益预估

| 优化 | 改动量 | 单项加速 | 适用范围 |
|------|--------|---------|---------|
| 1. 消除重复 forward_simulation | 小 | 25% 减少 | levenberg_marquardt |
| 2. forward_simulation 转 QobjEvo | 大 | 5~6x | 所有 forward_simulation 调用 |
| 3. adjoint H 预计算 | 中 | 2~5x | compute_jacobian |
| 4. 修复 g_values bug | 极小 | 迭代次数减半 | compute_jacobian |
| 5. store_states 按需 | 极小 | 内存降 500x | forward_simulation |
| 6. trace_values 向量化 | 小 | 小 | compute_jacobian |

### 优化前后对比

| 指标 | 优化前 | 优化后 |
|------|--------|--------|
| 单轮 forward_simulation | ~6s | ~1s |
| 单轮 compute_jacobian | ~12s | ~3s |
| 单轮总计 | ~24s | ~5s |
| 50 轮迭代 | ~20 min | ~4 min |
| 考虑 bug 修复后收敛加快 | 50 轮 | ~20 轮 |
| **总运行时间** | **~20 min** | **~1.5 min** |

---

## 5. 建议改动顺序

```
第 1 步：修复 g_values[i] → g_values[s] bug（1 行改动，可能让结果正确并减少迭代）
第 2 步：G_mat 外提 + store_states 参数（极小改动，立即见效）
第 3 步：消除重复 forward_simulation（接口改动，省 25% 时间）
第 4 步：forward_simulation 转 QobjEvo（主要工作量，5~6x 加速）
第 5 步：adjoint H 预计算（中等工作量，2~5x 加速）
第 6 步：trace_values 向量化（小改动，收尾优化）
```

第 1~3 步可在 30 分钟内完成，第 4~5 步需要更多时间但收益最大。
