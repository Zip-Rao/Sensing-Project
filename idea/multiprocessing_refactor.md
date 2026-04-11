# 多进程改造方案：代码结构调整说明

## 核心原则：多进程的三条铁律

多进程（multiprocessing）与多线程最本质的区别是：子进程是**全新的 Python 解释器**，它与主进程之间通过 pickle 序列化传递数据。因此所有传递给子进程的东西必须满足：

| 铁律 | 说明 | 违反后果 |
|------|------|---------|
| **Worker 必须是模块级函数** | 不能是类方法、嵌套函数、lambda | `PicklingError: can't pickle local object` |
| **参数必须可 pickle** | numpy 数组✅、基础类型✅、QuTiP Qobj ⚠️、闭包❌ | `PicklingError` 或进程卡死 |
| **QuTiP 对象在 worker 内部重建** | 主进程将 Qobj 转为 numpy 传入，worker 内部再包装成 Qobj | 避免跨进程 Qobj 序列化异常 |

> **Windows 额外注意**：Windows 使用 `spawn`（而非 `fork`）启动子进程，子进程会重新 `import` 整个模块。因此 worker 函数必须在 `.py` 文件顶层定义（Notebook 中无法直接用 multiprocessing，需将 worker 函数写在 `.py` 文件中再 import）。

---

## 数据流转变

```
当前（错误的多进程方式）：
  主进程 → pickle(H_total 闭包, self.single_measurement 方法) → 子进程
                    ↑ 不可 pickle，直接报错

改造后（正确方式）：
  主进程：QuTiP 对象 → .full() → numpy 数组
  主进程 → pickle(numpy 数组, 标量) → 子进程
  子进程：numpy 数组 → Qobj() 包装 → mesolve
```

---

## 第一处改动：`protocal.py`

### 改动位置：文件顶部（类定义之前）新增模块级 worker 函数

```python
# ========== 模块级 worker，必须在类外部 ==========
def _single_measurement_worker(args):
    """
    多进程 worker：完成单次测量的 mesolve 计算。
    接收纯 numpy/标量参数，在内部重建 QuTiP 对象。
    必须是模块级函数才能被 pickle。
    """
    (t_delay,
     H_0_mats,      # numpy (n_phi, dim, dim)：各时刻的 H_0 矩阵
     H_0_base,      # numpy (dim, dim)：无磁场时的 H_0
     H_pulse_mats,  # numpy (n_pulse, dim, dim)：脉冲 H 矩阵
     Phi_t,         # numpy (n_phi,)：磁场信号时间轴
     pulse_t,       # numpy (n_pulse,)：脉冲时间轴
     qubit_state_np,# numpy (dim, 1)：初始量子态
     n_levels,      # int
     ) = args

    # worker 内部导入 QuTiP，避免跨进程共享状态
    from qutip import Qobj, mesolve, basis, expect
    import numpy as np

    dims = [[n_levels], [n_levels]]
    n_phi = len(Phi_t)
    n_pulse = len(pulse_t)
    half_pulse = 0.5 * pulse_t[-1]
    delta = t_delay - half_pulse
    H_zero = np.zeros_like(H_0_base)

    def get_H_pulse_np(t_local):
        if t_local < 0 or t_local > pulse_t[-1]:
            return H_zero
        idx = np.clip(np.searchsorted(pulse_t, t_local), 0, n_pulse - 1)
        return H_pulse_mats[idx]

    def H_func(t, _args):
        t_sig = t + delta
        if t_sig < 0 or t_sig >= Phi_t[-1]:
            H_np = H_0_base
        else:
            idx = min(np.searchsorted(Phi_t, t_sig), n_phi - 1)
            H_np = H_0_mats[idx]
        return Qobj(H_np + get_H_pulse_np(t), dims=dims)

    t_start = min(delta, 0)
    t_end = max(Phi_t[-1], t_delay + half_pulse)
    t_evolve = np.linspace(t_start, t_end, n_phi + n_pulse - 1)

    # 在 worker 内部重建 QuTiP 对象
    psi0 = Qobj(qubit_state_np, dims=[[n_levels], [1]])
    psi_e_proj = basis(n_levels, 1) * basis(n_levels, 1).dag()

    result = mesolve(H_func, psi0, t_evolve, [], e_ops=[])
    final_state = result.states[-1]
    return float(expect(psi_e_proj, final_state))
# ========== worker 定义结束 ==========


class Protocal:
    ...
```

### 改动位置：`sliding_measrement` 方法

```python
def sliding_measrement(self, qubit: TransmonQubit, Phi_signal: Signal, control_pulse: CompositePulse):
    delay_start = Phi_signal.t_list[0] - 0.5 * control_pulse.t_list[-1]
    delay_end   = Phi_signal.t_list[-1] + 0.5 * control_pulse.t_list[-1]
    n_samples   = len(Phi_signal.t_list) + len(control_pulse.t_list) - 1
    scan_list   = np.linspace(delay_start, delay_end, n_samples)

    # ---- 预计算：QuTiP 对象 → numpy，子进程只接收 numpy ----
    qubit_t = qubit.qubit_under_mag(Phi_signal)
    if control_pulse.frame == 0:
        H_0_mats = np.array([qt.get_hamiltonian().full() for qt in qubit_t])
        H_0_base = qubit.get_hamiltonian().full()
    else:
        H_0_mats = np.array([qt.get_hamiltonian_rwa(control_pulse.omega_d).full() for qt in qubit_t])
        H_0_base = qubit.get_hamiltonian_rwa(control_pulse.omega_d).full()

    H_pulse_mats  = np.array([control_pulse.get_hamiltonian(t).full() for t in control_pulse.t_list])
    Phi_t         = np.asarray(Phi_signal.t_list)
    pulse_t       = np.asarray(control_pulse.t_list)
    qubit_state_np = qubit.state.full()

    # ---- 构建参数列表（每个元素是一个 tuple，全部可 pickle）----
    args_list = [
        (t_delay, H_0_mats, H_0_base, H_pulse_mats,
         Phi_t, pulse_t, qubit_state_np, qubit.n_levels)
        for t_delay in scan_list
    ]

    # ---- 多进程调用，backend='loky'（Windows 下比 'multiprocessing' 更稳定）----
    from joblib import Parallel, delayed
    p_e = Parallel(n_jobs=6, backend='loky')(
        delayed(_single_measurement_worker)(args) for args in args_list
    )
    return scan_list, p_e
```

---

## 第二处改动：`analysis.py`

### 改动位置：文件顶部（所有类和函数定义之前）新增模块级 worker

```python
# ========== 模块级 worker ==========
def _forward_sim_worker(args):
    """
    多进程 worker：完成单个 t_i 的 mesolve 计算（正向模拟一步）。
    """
    (t_i,
     H_0_mats,      # numpy (n_B, dim, dim)
     H_pulse_mats,  # numpy (n_pulse, dim, dim)
     t_list,        # numpy：滑动测量时间轴（用于 H_0 索引）
     B_t,           # numpy：磁场信号时间轴
     pulse_t,       # numpy
     qubit_state_np,# numpy (dim, 1)
     n_levels,      # int
     ) = args

    from qutip import Qobj, mesolve, basis, Options
    import numpy as np

    dims = [[n_levels], [n_levels]]
    n_B     = len(B_t)
    n_pulse = len(pulse_t)
    half_pulse = 0.5 * pulse_t[-1]
    H_zero  = np.zeros_like(H_0_mats[0])

    def get_H_pulse_np(t_local):
        if t_local < 0 or t_local > pulse_t[-1]:
            return H_zero
        idx = np.clip(np.searchsorted(pulse_t, t_local), 0, n_pulse - 1)
        return H_pulse_mats[idx]

    def H_func(t, _args):
        idx   = np.clip(np.searchsorted(t_list, t), 0, n_B - 1)
        H_np  = H_0_mats[idx]
        t_loc = t - t_i + half_pulse
        return Qobj(H_np + get_H_pulse_np(t_loc), dims=dims)

    t_m = min(B_t[0], t_i - half_pulse)
    t_M = max(B_t[-1], t_i + half_pulse)
    t_evolve = np.linspace(t_m, t_M, n_B + n_pulse)

    psi0        = Qobj(qubit_state_np, dims=[[n_levels], [1]])
    psi_e_proj  = basis(n_levels, 1) * basis(n_levels, 1).dag()
    options     = Options(store_states=True)

    return mesolve(H_func, psi0, t_evolve, [], e_ops=[psi_e_proj], options=options)
# ========== worker 定义结束 ==========
```

### 改动位置：`forward_simulation` 函数体

```python
def forward_simulation(qubit, control_pulse, B_curr, t_list):
    import time
    start = time.time()
    B = B_curr
    qubit_t = qubit.qubit_under_mag(B)

    # ---- 预计算 numpy 矩阵 ----
    if control_pulse.frame == 0:
        H_0_mats = np.array([qt.get_hamiltonian(qubit.frequency).full() for qt in qubit_t])
    else:
        H_0_mats = np.array([qt.get_hamiltonian_rwa(qubit.frequency).full() for qt in qubit_t])

    pulse_t      = np.asarray(control_pulse.t_list)
    H_pulse_mats = np.array([control_pulse.get_hamiltonian(t).full() for t in pulse_t])
    B_t          = np.asarray(B.t_list)
    qubit_state_np = qubit.state.full()

    # ---- 构建参数列表 ----
    args_list = [
        (t_i, H_0_mats, H_pulse_mats, t_list, B_t, pulse_t, qubit_state_np, qubit.n_levels)
        for t_i in t_list
    ]

    from joblib import Parallel, delayed
    result = Parallel(n_jobs=6, backend='loky')(
        delayed(_forward_sim_worker)(args) for args in args_list
    )

    end = time.time()
    print(f"Forward simulation time: {end - start:.2f} s")
    return result
```

---

## 为什么用 `backend='loky'` 而不是 `backend='multiprocessing'`

| 对比项 | `loky`（推荐） | `multiprocessing` |
|--------|---------------|-------------------|
| Windows 支持 | ✅ 专门针对 Windows spawn 优化 | ⚠️ 在 Jupyter 下经常卡死 |
| 序列化方式 | cloudpickle（更强） | pickle（受限） |
| worker 进程复用 | ✅ 自动复用，避免重复启动开销 | ⚠️ 需手动管理 Pool |
| 异常传播 | ✅ 清晰传回主进程 | ⚠️ 有时无法捕获 |

---

## 改造前后对比总结

| 位置 | 改造前 | 改造后 |
|------|--------|--------|
| worker 位置 | 类方法 `self.single_measurement` | 模块级函数 `_single_measurement_worker` |
| H_total 形式 | 嵌套闭包（不可 pickle） | worker 内部定义的普通函数（不跨进程传递） |
| 传给子进程的数据 | `H_total` 闭包 + QuTiP 对象 | 纯 numpy 数组 + 标量 |
| QuTiP 对象处理 | 直接传递（报错） | 主进程 `.full()` → numpy，worker 内 `Qobj()` 重建 |
| backend | `'multiprocessing'`（失败） | `'loky'`（稳定） |

---

## Jupyter Notebook 的特殊处理

Windows 的 `spawn` 模式要求 worker 函数**必须可从模块 import**。在 Notebook 中直接定义的函数无法被子进程 import，因此：

```python
# ❌ 不能这样（Notebook Cell 中定义 worker）
def my_worker(args):
    ...
Parallel(n_jobs=4, backend='loky')(delayed(my_worker)(a) for a in args)

# ✅ 应该这样（worker 在 .py 文件中定义，Notebook 中 import）
# protocal.py 中已定义 _single_measurement_worker
from src.protocal import _single_measurement_worker  # 可以被子进程 import
```

由于 `_single_measurement_worker` 定义在 `src/protocal.py` 的模块级别，子进程可以正常 import，无需额外处理。

---

## 数据量估算（pickle 开销评估）

每次任务传递的数据量：

| 数据 | 大小（n_phi=400, dim=2）|
|------|------------------------|
| `H_0_mats` | 400 × 2 × 2 × 16 bytes ≈ **25 KB** |
| `H_pulse_mats` | 20 × 2 × 2 × 16 bytes ≈ **1.3 KB** |
| `Phi_t`, `pulse_t` | < 5 KB |
| `qubit_state_np` | 2 × 1 × 16 bytes ≈ 0.03 KB |
| **单次任务总计** | **≈ 32 KB** |

loky 的进程池在整个 `Parallel` 调用期间保持，32 KB × 420 任务 = **~13 MB** 总序列化量，在 DDR5 内存下可以忽略不计（< 0.1s）。
