# AI 交接文档 · sqc 文档站

> **读者:AI 助手(而非人类)。** 这份文件让任何后续会话都能一致地处理 sqc 文档站的
> 三类工作:**润色中文、翻译英文、写作之外的操作(构建/CI/git/边界)**。
> 人类的待办清单在 [CONTENT_TODO.md](CONTENT_TODO.md);本文件是给你(AI)的约定。
>
> **黄金规则:任何改动后必须 `sphinx-build -W` 两树全绿再交付。** 见 §5。

---

## 0. 你的三类职责

1. **润色中文** — 用户写中文草稿,你改语言(见 §3),**不改物理判断/参数取舍**。
2. **翻译英文** — 从中文生成 `en/` 对应页(见 §4)。英文与中文**内容对等、结构一致**。
3. **写作外操作** — 构建验证、修 MyST/交叉引用、CI 排障、git 安全、守住 v1 边界(§5–§8)。

用户负责:物理内容与取舍、GitHub 网页操作、git push(见 CONTENT_TODO Part 1)。

---

## 1. 站点结构(动手前必读)

```
docs/source/
  _shared/            共享配置 —— 两树的唯一真源
    conf_base.py        所有公共 Sphinx 设置(改这里,不要各自改)
    _templates/, _static/
  en/                 英文树(站点根 /)      conf.py 只 override language/title
  zh/                 中文树(站点 /zh/)     同上
```

- **两树是各自独立的 sourcedir**,不是 gettext 翻译。en 和 zh 各有一套完整 `.md`。
- **改公共设置只动 `_shared/conf_base.py`**;各 `conf.py` 只覆盖 `language` 和标题。
- **每页几乎都成对**:`en/X.md` 与 `zh/X.md` 内容对等。改一个通常要改另一个。
- 已双语写完的 8 页:overview / install / quickstart / architecture / theory /
  extending / roadmap / tutorial。**待写 11 页**:building_blocks×8 + examples×3
  (用户写中文,你翻英文)。

---

## 2. 环境与构建命令

- Python 环境:conda `qutip-env`。Windows 全路径:
  `C:\Users\21034\anaconda3\envs\qutip-env\python.exe`
- 文档依赖装法:`pip install -e ".[docs]"`(sphinx, furo, myst-parser, nbsphinx)。
  nbsphinx 另需系统 `pandoc`(CI 用 apt 装;本地 qutip-env 已有)。
- **构建单树(严格):**
  ```bash
  cd docs
  "<qutip-env python>" -m sphinx -W --keep-going -b html source/en build/html
  "<qutip-env python>" -m sphinx -W --keep-going -b html source/zh build/html/zh
  ```
- **构建两树(便捷):** `cd docs && make html`(Windows:`make.bat`)。
- **notebook 副本:** tutorial 的 `Simulation_sqc.ipynb` 单一真源在**仓库根**;
  Makefile/make.bat/CI 在构建前把它 copy 进 `en/` 和 `zh/`。副本已 gitignore。
  手动单树构建前需先 `cp ../Simulation_sqc.ipynb source/en/`(zh 同理)。
- ⚠️ **Windows 上 notebook 渲染很慢(约 3–4 min/树,pandoc 逐 cell)**。只改文字页时,
  后台跑构建(`run_in_background`)再回收结果,别干等。CI(Linux)则很快。

---

## 3. 中文润色约定

- **只改语言,不改物理。** 参数值、协议选择、物理论断是用户的领域,不要"顺手纠正"。
  觉得物理有误,**提出来问用户**,不要擅改。
- 语气:技术、简洁、专业。术语首次出现给中英对照(如"甜点(sweet spot)")。
- **保留原文的技术名词与代码标识符**:类名、函数名、`Φ_0`、`f_01` 等不翻译、不改写。
- 数学:行内 `$...$`,独立 `$$...$$`(已启用 `dollarmath`+`amsmath`)。
- 代码块里的注释可译中文,但**代码本身一字不动**。

---

## 4. 英译约定

- **内容对等、结构一致**:en 页与 zh 页同样的小节、同样的代码、同样的表格行数。
- **代码块完全相同**(含变量名);只有代码注释和散文译成英文。
- **数学公式完全相同**。
- 交叉引用、`{toctree}`、admonition 等指令**两树逐一对应**。
- 英文风格:与已写的 8 页保持一致(技术、平实、不夸张)。可参照 `en/architecture.md`、
  `en/theory.md` 的既有笔调。
- 术语一致性:同一概念在全站用同一英文词(如始终 "sweet spot"、"flux signal"、
  "reconstruction")。

---

## 5. 交付前的强制验证(黄金规则)

**任何 `.md`/`.ipynb`/`conf` 改动后,必须两树 `-W` 全绿才算完成。** 顺序:

1. 需要时 copy notebook 进两树(见 §2)。
2. 跑 §2 的严格构建命令(en + zh)。
3. `grep -ci warning` 两树日志都应为 0,exit code 0。
4. 有告警 → 定位并修(§6),不要用 `suppress_warnings` 掩盖结构性问题。

> `-W` 把 warning 当 error。这是**特性**:它守住 toctree、交叉引用、缺文件等结构问题。
> 已知的、非结构性的 docstring 告警已在 conf 里定向抑制(§7),不要扩大抑制范围。

---

## 6. 已解决的坑(别重蹈,也是排障字典)

| 症状 | 根因 | 正确处理 |
|---|---|---|
| CI 挂但本地过,49s 早崩 | notebook 缺 cell `id`(nbformat 4.5),本地旧 nbformat 自动兜底、CI 新版报硬错 | notebook 已补 id;若再遇新 notebook,先 `nbformat` 补 id |
| `WARNING: lexer 'ipython3' not known` ×N | notebook 代码 cell 标 `ipython3`,该 lexer 需 ipython 包,CI `[docs]` 没装 | conf 已把 `ipython3` 别名到 `PythonLexer`(§7),**不要**为此加 ipython 依赖 |
| `toc.no_title` ×N(每页侧栏一次) | 曾把 `.ipynb` 塞进 `source_suffix`,破坏 nbsphinx 标题提取 | `source_suffix` **不要**含 `.ipynb`;交给 nbsphinx 自己注册 |
| `docutils` 替换引用告警 | `sqc/` docstring 用 `|x|` 纯文本记号 | 已 `suppress_warnings=["docutils"]`;别扩大;真正修法是清 docstring(B3,待定) |
| autosummary 重复对象告警 | `:recursive:` 下潜子模块 + `__all__` 重导出同名 | 保持现有 autosummary 配置(不加 recursive);`autosummary_ignore_module_all=False` |
| `sphinx-autodoc-typehints` 崩(Sphinx 8.2) | 与新版不兼容(`short_literals`) | **不要**加回该扩展;用原生 `autodoc_typehints="description"` |

---

## 7. conf_base.py 里的关键决策(动它之前先懂 why)

改 `docs/source/_shared/conf_base.py` 前必读——这些设置都是踩坑后的定论,注释里有 why:

- **不启用 `sphinx.ext.viewcode`** — 它会渲染全部源码,泄露未公开的 stub(如
  `ZCrosstalkWorkflow`)。源码在 GitHub 仓库可见即可,文档站不铺。
- **`autosummary_ignore_module_all = False`** — 尊重每个子包的 `__all__`,故意排除的
  名字(D3 的 post-v1 功能)不进 API 页。
- **`nbsphinx_execute = "never"`** — 不重跑 notebook(QuTiP 重、要确定性)。用保存的输出。
- **`source_suffix` 不含 `.ipynb`** — 见 §6。
- **`ipython3` → `PythonLexer` 别名** — 见 §6,零依赖消告警。
- **`_skip_notimplemented_members` 钩子** — 自动跳过方法体含 `raise NotImplementedError`
  的成员(`__all__` 到不了方法级)。**名字无关**,自动覆盖当前和未来的 stub。守 v1 边界。
- **`suppress_warnings = ["docutils"]`** — 只压这一类,别加别的。

---

## 8. v1 公开边界(守住,别泄露)

文档站**只呈现 v1 公开面**。以下内容存在于代码但**不得出现在**公开文档/API 页/前端:

- `ZCrosstalkWorkflow`、Z 串扰重建
- 瞬态**频率标定**(`method="transient"`)——注意:瞬态**波形重建**是 v1 核心特性,已支持,**别混淆**
- CPMG 协议、`TunableCoupler`、显式电子学层
- `SensingWorkflow` 上抛 `NotImplementedError` 的方法(save/load/diff/benchmark/…)

这些的完整列表和"为什么隐藏"见 `en/roadmap.md` / `zh/roadmap.md`。
`__all__` 过滤 + `_skip_notimplemented_members` 钩子(§7)自动挡住大部分;
**新写内容时不要主动提及或链接这些**。拿不准是否属于 v1,查 roadmap 页或问用户。

---

## 9. git / 部署安全

- 默认分支 **`main`**;远程 `origin` = `git@github.com:Zip-Rao/Sensing-Project.git`。
- **不要擅自 commit/push**——除非用户明确要。用户自己 push(见 CONTENT_TODO Part 1)。
- CI:`.github/workflows/docs.yml`,push 到 `main` 触发 → 构建 en(根)+ zh(`/zh/`)
  → 部署 GitHub Pages。用 `${{ github.* }}` 变量,**不硬编码用户名**。
- 站点地址:`https://zip-rao.github.io/Sensing-Project/`(owner 段小写)。
- 破坏性 git 操作(force push、reset --hard、删远程分支)一律先经用户确认。

---

## 10. 交接给下一个 AI 的最短清单

接手时按序确认:
1. 读本文件 + [CONTENT_TODO.md](CONTENT_TODO.md)。
2. 用户要写作(润色/翻译)还是操作?对号入座 §3/§4 或 §5–§9。
3. 动手前:`git status` 看有无未提交改动;读要改的那页的 en+zh 两份。
4. 改完:两树 `-W` 构建全绿(§5)才交付。
5. 不确定物理/边界:问用户,别猜。