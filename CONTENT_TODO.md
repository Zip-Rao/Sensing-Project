# 文档站 · 你的待办清单

> 这份文件是你的私人工作清单,放在仓库根、**不在 `docs/source/` 下**,所以
> **不会被 Sphinx 发布到网站**。写完可删,或加进 `.gitignore` 不追踪。
> 更新时间基准:2026-07(以当前会话为准)。

本文件三部分:
- **Part 1 — 操作清单**(GitHub / git,我改不了的)
- **Part 2 — 写作清单**(11 页文稿,你写中文,我翻英文+接入)
- **Part 3 — 交接与预览流程**(怎么把一页交给我、怎么本地预览)

---

## Part 1 · 操作清单(手动,我改不了)

### 已完成(供你确认,无需再做)
- [x] GitHub 用户名 → `Zip-Rao`(本地 remote 已指向新名)
- [x] 默认分支 → `main`
- [x] `docs.yml` CI 触发器 → `branches: [main]`
- [x] CI 构建修复:notebook 补 cell id + `ipython3` lexer 别名(已提交)

### 仍需你手动做

**OP-1 · 提交并推送 install.md 的用户名改动**
```bash
cd "c:/Users/21034/Desktop/Workspace/Sensing project/Sensing-Project"
git add docs/source/en/install.md docs/source/zh/install.md
git commit -m "docs: update clone URL to Zip-Rao"
git push
```

**OP-2 · 放行 github-pages 环境部署 main(关键,否则 deploy 被拒)**
GitHub 仓库 → Settings → Environments → `github-pages` →
"Deployment branches and tags" → 把允许的分支改成(或新增)`main`,删掉旧分支规则。

**OP-3 · 确认 Pages 源 = GitHub Actions**
Settings → Pages → Build and deployment → Source 应为 **GitHub Actions**(不是某个 branch)。

**OP-4 · (可选)清理遗留远程分支**
现在 origin 上还有 `master`、`项目重建`、`项目重建-v2` 三个旧分支。
- `项目重建-v2` 已被重命名为 `main`(GitHub rename),若列表里还在是残留,可删。
- 旧 `master` 内容(若要保留)先推到 `branch` 再删,见下:
```bash
git push origin master:refs/heads/branch    # 保留旧 master 内容到 branch
git push origin --delete master              # 删除操作,确认 branch 已在再删
git push origin --delete 项目重建 项目重建-v2  # 删中文残留分支
```
> `--delete` 不可逆;删前去 GitHub 确认要保留的内容已在别处。

**OP-5 · (上线后验证)访问站点**
- 英文:https://zip-rao.github.io/Sensing-Project/
- 中文:https://zip-rao.github.io/Sensing-Project/zh/
- 检查:根是英文、`/zh/` 是中文、右上角语言切换器能互跳、tutorial 页 notebook 图正常。
  (Pages URL 里 owner 段是**小写** `zip-rao`。)

### 待你拍板的决策(不表态就按默认)
- **B1** tutorial notebook 出英文走查版? **默认:不做**,保留中文 + EN 页已注明"代码语言中立"。
- **B2** Overview 页(我从 README 派生的双语初稿)自己重写还是留用? **默认:留用**,你随时改。
- **B3** 清理 `sqc/` docstring 里触发 docutils 告警的 `|x|` 记号? **默认:不做**,已用 `suppress_warnings=["docutils"]` 兜住,`-W` 仍守结构问题。

---

## Part 2 · 写作清单(11 页,你写中文 → 我翻英文+接入)

**分工:你写中文,我负责翻译成英文、接入两树、修 MyST 语法、跑构建验证。**
每页的中文文件已存在为占位符(`内容即将补充`),路径见下表。你只需**覆盖写中文正文**,
交给我即可。

### 写作通用规则(所有页适用)
1. **只写正文**,不用管英文、不用管 `{toctree}`、不用管交叉引用格式——我来接。
2. 引用某个类/函数时,直接写类名即可(如 `TransmonQubit`),我会转成可点击的
   `{py:class}` 交叉引用。你也可以自己写 ``{py:class}`~sqc.devices.TransmonQubit` `` 但非必须。
3. 代码块用普通 ` ```python `,我保证能跑再放进去(拿不准的代码标一下,我帮你验证)。
4. 数学公式用 `$...$`(行内)或 `$$...$$`(独立),已启用 dollarmath。
5. 物理判断、参数取舍、协议选择这些**只有你懂的内容是重点**;结构性排版交给我。

### C1 · Building Blocks × 8

**每页统一 4 段骨架:**
```
## 这层提供什么      (1-2 句:这层在栈里负责什么)
## 关键类            (列出主要类,各一句职责)
## 最小用例          (5-15 行能跑的代码)
## 物理角色 / 扩展    (对应真实物理组件;想加自己的东西继承哪个 ABC)
```

| # | 中文文件(你写) | 关键类 | 扩展点 ABC |
|---|---|---|---|
| 1 | `docs/source/zh/building_blocks/devices.md` | `TransmonQubit`, `QubitSpec`, `Resonator`, `ChipTopology`, `CoupledSystem` | `Device` |
| 2 | `docs/source/zh/building_blocks/hardware.md` | `ControlLine`, `DistortionModel`(及子类), `TransferMatrix` | `DistortionModel`, `ReadoutModel` |
| 3 | `docs/source/zh/building_blocks/control.md` | `Waveform`, `FluxSignal`, `Pulse`, `PulseSequence`, 门工厂 | `PulseBase` |
| 4 | `docs/source/zh/building_blocks/simulation.md` | `HamiltonianBuilder`, `MesolveRunner`, `SlidingMeasurementRunner`, `ExperimentResult` | `RunnerBase` |
| 5 | `docs/source/zh/building_blocks/experiments.md` | `RabiExperiment`, `RamseyExperiment`, `DiffEchoExperiment`, `TransientSensingExperiment`, `CryoscopeExperiment` | `Experiment` |
| 6 | `docs/source/zh/building_blocks/reconstruction.md` | `KernelEstimator`, `RamseyReconstruction`, `EchoReconstruction`, `TransientReconstruction`, `CryoscopeReconstruction` | `Reconstruction` |
| 7 | `docs/source/zh/building_blocks/calibration.md` | `Calibration`, `CalibrationTable`, `WaveformCalibration`, `PredistortionDesigner` | `Calibration` |
| 8 | `docs/source/zh/building_blocks/workflows.md` | `SensingWorkflow`, `PredistortionValidationWorkflow` | `Workflow` |

### C2 · Examples × 3

**每页统一 4 段骨架:**
```
## 目标            (1 句话:这条管道解决什么问题)
## 物理原理        (背后的物理;可引用 theory 页)
## 端到端代码      (完整可跑的例子)
## 结果解读        (输出说明,配图可选)
```

| # | 中文文件(你写) | 对应主线 |
|---|---|---|
| 1 | `docs/source/zh/examples/waveform_reconstruction.md` | 主线1:从 Ramsey/Cryoscope/瞬态恢复 Φ(t) |
| 2 | `docs/source/zh/examples/frequency_calibration.md` | 主线2:Ramsey `f(Φ)`/`f₀₁` 标定 |
| 3 | `docs/source/zh/examples/predistortion.md` | 主线3:AWG→芯片传函 + 补偿滤波器 |

### 填好的样板(以 devices 层为例,照着套)

复制下面这段到 `docs/source/zh/building_blocks/devices.md`,把方括号内容换成你的:

````markdown
# 器件(devices)

## 这层提供什么

devices 层是全栈最底层,建模物理量子器件本身——[一句话:比如"超导 Transmon
量子比特及其读出谐振腔、芯片拓扑",不涉及任何控制或仿真逻辑]。

## 关键类

- `TransmonQubit` — [一句职责:比如"单个可调频 Transmon,由 QubitSpec 参数化"]
- `QubitSpec` — [不可变参数容器(frozen dataclass):E_J, E_C, ...]
- `Resonator` — [读出谐振腔]
- `ChipTopology` / `CoupledSystem` — [多比特芯片布局 / 耦合系统]

## 最小用例

```python
from sqc.devices import TransmonQubit, QubitSpec

spec = QubitSpec(E_J=15.0, E_C=0.2)      # GHz
qubit = TransmonQubit(spec)
print(qubit.frequency(flux=0.0))          # [说明输出:比如甜点频率]
```

## 物理角色 / 扩展

`TransmonQubit` 对应真实芯片上的一个物理比特;`QubitSpec` 的参数直接来自
[器件表征/流片参数]。要加新器件类型(如 fluxonium),继承 `Device` 基类并实现
[列出必须实现的方法]。
````

> 我拿到你的中文后,会:翻成英文对应文件 → 把类名转成可点击交叉引用 →
> 验证代码能跑 → 两树 `-W` 构建通过。你不用管这些。

---

## Part 3 · 交接与预览流程

### 你交一页给我的方式(任选)
- **A**(推荐):直接把中文写进对应的 `zh/.../*.md` 文件,然后跟我说"devices 写好了"。
- **B**:把中文正文贴在对话里,告诉我是哪一页,我来落文件。

我收到后当次就会:翻译英文 + 接入两树 + 修语法 + 跑构建,然后告诉你结果。

### 建议顺序
1. **先写 1 页**(推荐 devices,最基础),我跑通完整流程,确认模板合你意。
2. 顺则**批量写 building_blocks 剩余 7 页**。
3. 最后写 **examples 3 页**(依赖对各层的理解,放最后最顺)。

### 本地预览(想自己看效果时)
```bash
cd "c:/Users/21034/Desktop/Workspace/Sensing project/Sensing-Project/docs"
make html          # 构建 en + zh(Windows 用 make.bat 或直接 sphinx-build)
# 打开 build/html/index.html(英文)、build/html/zh/index.html(中文)
```
> notebook 渲染在 Windows 上较慢(几分钟),只改文字页可先只构建单页目录。

### 进度勾选(写完一页打个勾)
- [ ] building_blocks/devices
- [ ] building_blocks/hardware
- [ ] building_blocks/control
- [ ] building_blocks/simulation
- [ ] building_blocks/experiments
- [ ] building_blocks/reconstruction
- [ ] building_blocks/calibration
- [ ] building_blocks/workflows
- [x] examples/waveform_reconstruction
- [x] examples/frequency_calibration
- [x] examples/predistortion
