# RELEASE_TODO — v1 发布就绪清单

> **用途**:跟踪 Sensing-Project 首个正式版本(v1)发行前的待办与已完成项。可持续维护:完成一项就把 `- [ ]` 改成 `- [x]` 并在末尾维护记录表追加一行;新增项按优先级插入对应小节。
> **目标版本**:v1(版本号待定,见下方"待定决策")
> **最后更新**:2026-07-17
> **配套文档**:[docs/architecture.md](docs/architecture.md)(开发者)、[idea/refactor/_handoff_state.md](idea/refactor/_handoff_state.md)(重构交接)、[idea/_TODO_master.md](idea/_TODO_master.md)(Track B 功能进度)

## 图例

- 优先级:**P0** = 发布阻断(不做不能发) · **P1** = 应有(影响质量/可用性) · **P2** = 可后置 / post-v1
- 类别标签:`[发布]` 打包与发行规格 · `[质量]` 正确性/精度 · `[文档]` · `[功能]` stub/未实现 · `[工程]` 工具链
- 状态:`- [ ]` 未完成 · `- [~]` 进行中 · `- [x]` 已完成

---

## 0. 待定决策(需用户拍板,决定下方优先级)

- [x] **D1 版本号**(2026-07-17 定):**1.0.0**。首个正式版,承诺公共 API(三主线+前端)稳定,后续 semver;统一代码 0.1.0 与文档 v0.3.0 → 1.0.0。
- [x] **D2 许可证**(2026-07-17 定):**MIT**,版权行 `Copyright (c) 2026 Zip`(署名已确认)。
- [x] **D3 v1 范围**(2026-07-17 定):三主线(波形重建/Ramsey 标定/预失真)+ 前端为**正式功能**;**Z-crosstalk 与 transient 频率标定完全从 v1 公共接口隐藏**(前端 tab、公共 API `__all__`、README 均不出现),不放 experimental 标记,post-v1 再引入。⚠️ 注意区分:transient **波形重建**(pipeline A,`TransientSensingExperiment`)是核心功能,**保留**;隐藏的是 transient **频率标定**(`FluxResponseCalibration(transient)` 等)。
- [x] **D4 SensingWorkflow 11 个 stub**(2026-07-17 定):**全部标 "planned in vX" 并从公共 API/文档隐藏,不实现、不委托**(与 D3 隐藏 crosstalk 一致 —— `crosstalk()` 也保持隐藏,不委托给 ZCrosstalkWorkflow)。

---

## 1. P0 — 发布阻断项

- [ ] **`[发布]` 新增 README.md** — 用户/开发者入口,链接快速开始、教程 notebook（`Simulation_sqc.ipynb`）、前端（`web_demo_v2.py`）、架构文档。
- [ ] **`[发布]` 新增 LICENSE(MIT)** — 标准 MIT 文本 + 版权行(待署名)。
- [ ] **`[发布]` 统一版本号 → 1.0.0** — [sqc/__init__.py:5](sqc/__init__.py#L5) `0.1.0` → `1.0.0`;[docs/architecture.md](docs/architecture.md) 头部 "sqc v0.3.0" → v1.0.0(文档正文版本 v2.12 是文档自身版本,另行处理)。
- [ ] **`[发布]` 新增 pyproject.toml** — 支持 `pip install -e .`,使 `import sqc` 不依赖 cwd / PYTHONPATH。
- [ ] **`[发布]` 清理 requirements.txt** — `pandas` / `scikit-learn` 在 `sqc/` 中未使用(违反 R10);删除或说明;core 与 demo 依赖分开(已有 `requirements_demo.txt`)。
- [ ] **`[发布]` 从 v1 公共接口隐藏 crosstalk + transient 频率标定**(D3 定案) — 使 v1 对外只呈现三主线+前端:
  - 前端 [web_demo_v2.py](web_demo_v2.py):移除/停用 **Crosstalk tab**(及涉及 transient 频率标定的入口);保留 Transient **波形重建** tab。
  - 公共 API:从 `sqc/workflows/__init__.py` `__all__` 移除 `ZCrosstalkWorkflow`;从 `sqc/calibration/__init__.py` 移除 transient 频率标定的公共导出(类仍可深路径导入,只是不再是公共面)。
  - README / 文档:不列 crosstalk 与 transient 标定;它们进 post-v1 roadmap。
  - Z-crosstalk 数值 bug(`compensation_factor≈0.276`<1、phi_B 偏大 5×,测试仅断言 `>0` 掩盖)→ 记入 post-v1(见 §3)。

---

## 2. P1 — 应有项

- [ ] **`[文档]` 前端文档化** — `RUN_DEMO.md` 只讲旧 src 版 `web_demo.py`,未提 sqc 版 `web_demo_v2.py`;更新或新增说明,并加前端启动 smoke 测试(构建 `build_app()`,不 `.launch()`)。
- [ ] **`[质量]` 统一瞬态重建路由** — workflow 走 legacy `pulse.get_kernel()`,前端/测试走 `KernelEstimator`;统一路径,并加真实网格精度 sanity 测试(当前短网格幅度恢复仅 ~7%)。
- [ ] **`[功能]` SensingWorkflow 11 个 stub 全部隐藏**(D4 定案) — **不实现、不委托**;全部改成清晰的 "planned in vX" 提示,并从公共 API/文档隐藏(不出现在 README、`__all__`、前端)。位置 [sqc/workflows/sensing.py](sqc/workflows/sensing.py)(:1033/1054/1072/1082/1097/1112/1127/1148/1178/1200/1220)。当前 `NotImplementedError` 文案已接近"planned",主要工作是确保它们不在公共入口/文档暴露。
- [ ] **`[文档]` 新增"已知限制 / post-v1 路线图"章节** — 列出 v1 **不包含**的能力(D3 定案:隐藏而非标 experimental):Z-crosstalk、transient 频率标定、CPMG、`TunableCoupler`([coupler.py:17](sqc/devices/coupler.py#L17))、electronics AWG/ADC([electronics.py:22](sqc/hardware/electronics.py#L22))、`schedules.py` 空壳;以及已知限制 `collapse_operators()` 在 chip/resonator 返回 `[]` 静默丢弃耗散。定位为"计划中",不是"v1 里的实验特性"。
- [x] **`[文档]` 整顿 TODO 体系** — 仓库现有 7+ 个 TODO/追踪文件,分工混乱、状态过时、互相重叠。确立单一权威:**发布工程看本文件、重构进度看 [`_handoff_state.md`](idea/refactor/_handoff_state.md)**。(2026-07-17 完成)
  - [x] 给 `idea/_TODO_master.md` / `idea/TODO.md` / `idea/TODO_noise.md` / `idea/_output.md` 头部加时效/定位声明(2026-07-17)。
  - [x] 刷新 `_TODO_master.md` 状态:新增"进度对照"权威映射表(对照 `_handoff_state.md` P0–P12;因多数任务为 △/部分或针对冻结的 src/,采用映射表而非逐格翻 ○→✓,更准确)(2026-07-17)。
  - [x] 回填 `_output.md` 的差分回波 / CPMG 段:差分回波补方法+实现指针(结果图仍欠账),CPMG 标注未实现 + 指向 post-v1(2026-07-17)。
- [ ] **`[发布]` 新增 CITATION.cff** — 作者/版本/引用信息,为被引用复用做准备;后续可选 Zenodo DOI / 投 JOSS。
- [ ] **`[文档]` 新增 CHANGELOG.md** — 记录 v1。

---

## 3. P2 — 可后置 / post-v1

- [ ] **`[工程]` 新增 argparse CLI** — 标准库(不违反 R10),用于参数扫描 / 批量 / HPC 批处理运行。
- [ ] **`[工程]` 清理被跟踪的 `src/__pycache__`** — `.gitignore` 已覆盖 `__pycache__/` 但历史文件仍被跟踪;`git rm -r --cached src/__pycache__` 清一次,解决 handoff known issue #3。
- [ ] **`[质量]` collapse_operators() 静默丢弃耗散** — chip/resonator 返回 `[]` 时改为 `warnings.warn`,避免使用者无感知地丢掉耗散。
- [ ] **`[工程]` 新增 CI** — push 时跑 `pytest`(文档站的 build/发布 CI 归入下一条)。
- [ ] **`[文档]` 完整托管文档站(Sphinx + GitHub Pages)** — 对标 [scq-cloud.github.io](https://scq-cloud.github.io/)(PyQuafu)/ scqubits 那种可导航、可搜索的在线文档。**完整层次**,拆为:
  - [ ] Sphinx 骨架 + 主题(furo 或 sphinx-rtd-theme);`conf.py`、`docs/source/`。
  - [ ] `autodoc` + `autosummary`:从现有 docstring 自动生成 **API 参考**(逐类逐方法;素材现成、单位标注齐全)。
  - [ ] `myst-parser` 纳入现有 `README.md` / `docs/architecture.md`;`nbsphinx` 渲染 `Simulation_sqc.ipynb` 为教程页。
  - [ ] 站点结构:Overview → Install → Quickstart → User Guide → API Reference → Theory(选摘 idea/)→ Roadmap。
  - [ ] GitHub Actions:build → 发布到 **GitHub Pages**(`<user>.github.io` 风格)。
  - [ ] **依赖(R10 需批准)**:sphinx 及插件(sphinx / myst-parser / nbsphinx / furo / sphinx-autodoc-typehints)——**仅文档构建工具,非运行时**;放入 `pyproject` 的 `[docs]` extra,不污染 core。引入前需用户确认(已初步同意"完整层次",实施前再确认依赖清单)。
- [ ] **`[质量]` (post-v1)修复 Z-crosstalk 重建标度 bug** — `compensation_factor≈0.276`(<1)、重建 phi_B 偏大 ~5×;根因是小串扰下 Wiener 瞬态重建核标度误差。修复后收紧集成测试(当前仅断言 `>0`),再考虑纳入正式功能。见 [sqc/workflows/z_crosstalk.py](sqc/workflows/z_crosstalk.py)。
- [ ] **`[功能]` (post-v1)实现被隐藏的能力** — transient 频率标定、CPMG、TunableCoupler、electronics 层、SensingWorkflow 的 planned 方法,按需逐步引入。

---

## 4. 已完成(本轮代码审查与整理)

- [x] **代码审查(sqc/ 全量,9 子包 ~15k 行)** — 六轴审查(完成度/死代码/可读性/文档/逻辑/可行性),findings 见会话记录。
- [x] **第一节:5 处真实缺陷修复** — commit `89e08aa`
  - echo 重建 arcsin 前加 clip 防 NaN;`create_pulse` 修 DRAG `Omega_Q` 缩放 + 除零守卫;distortion 双线性逆 `amp→1` identity 兜底;`create_ramsey_pulse` trigger 推进与 echo 系工厂统一。
  - chip `Qobj→complex` 经 QuTiP 5.2.0 实测正常,无需改(仅记录)。
  - 新增测试:echo NaN 防护 3 + distortion amp→1 稳定性 1。
- [x] **第二节:包导出与清理** — commit `a6153631`
  - `sqc/control/__init__.py`、`sqc/simulation/__init__.py` 补全公共 API 导出;`sqc/hardware/__init__.py` 补 `CascadeDistortion`;`sqc/calibration/waveform.py` 删重复 `@staticmethod`。
- [x] **发布就绪审查** — 三主线 + 前端端到端验证(见会话);本清单生成。

安全点:改前 `190072d` · 第一节 `89e08aa` · 第二节 `a6153631`。`src/` 全程未改(R1)。

---

## 维护记录(append-only,最新在上)

| 日期 | 变更 | 关联 commit |
|---|---|---|
| 2026-07-17 | 新增 P2 项「完整托管文档站(Sphinx + GitHub Pages)」,对标 PyQuafu/scqubits;标注 R10 依赖需批准 | (pending) |
| 2026-07-17 | 定案 D1(1.0.0)/ D2(MIT);LICENSE、版本统一条目具体化 | (pending) |
| 2026-07-17 | 定案 D3(隐藏 crosstalk/transient 标定)/ D4(11 stub 全隐藏);据此细化 crosstalk/stub/限制条目,数值 bug 移 post-v1 | (pending) |
| 2026-07-17 | 新增"整顿 TODO 体系"项;给 4 个过时/易混 idea/ 文件头部加时效声明 | (pending) |
| 2026-07-17 | 初始化本清单;记入本轮代码审查修复(§4) | 89e08aa, a6153631 |
