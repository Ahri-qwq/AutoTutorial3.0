# ABACUS 旧教程更新 — 分批生成执行计划

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 按照设计文档 `docs/plans/2026-03-10-tutorial-topics-design.md` 的优先级，利用 CLAUDE.md 的 7 步流程，逐批生成 ABACUS 教程，共 62 个主题（4 组）。

**Architecture:** 每篇教程独立执行一次 CLAUDE.md 完整流程（Step 0–7），工作目录按 `_workspace/YYYYMMDD_HHMMSS_主题/` 命名，最终产物存入 `data/knowledge_node/Tutorials/` 并更新知识树状态。

**Tech Stack:** CLAUDE.md 7步流程、`tools/retriever.py`（RAG检索）、`tools/case_parser.py`（案例解析）、`tools/orbital_validator.py`（轨道验证）、Markdown

---

## 使用前必读

### 任务类型对照

| 设计文档分组 | CLAUDE.md 任务类型 | 核心操作 |
|------------|-----------------|---------|
| Tier 1 ⭐⭐⭐ | 类型C（案例驱动） | 提供源文件路径，围绕案例展开 |
| Tier 2 ⭐⭐ | 类型B（基于案例） | 提供源文件路径，系统补全缺失文件 |
| 安装类 / Tier 3 ⭐ | 类型A（无案例） | 仅依赖 RAG 检索 |

### 每篇教程的标准启动方式

开启新对话，提供以下信息：
```
主题：[主题名称]
任务类型：C（案例驱动）/ B（基于案例）/ A（无案例）
案例文件：data/knowledge_source/[文件名]（可多个）
特殊要求：[如有]
```

### 最终产物存放规则

1. 在 `_workspace/<时间戳_主题>/` 内完成全部 7 步
2. 将 `07_Final_Tutorial_<标题>.md` 复制到 `data/knowledge_node/Tutorials/`
3. 更新 `data/knowledge_node/knowledge_tree.md` 和 `knowledge_tree_with_links.md`：
   - 将对应知识点的 ⬜/🔄 改为 ✅
   - 追加 `| 新教程标题 | 知识点位置 | 文件路径 |` 到新教程对应关系表
4. 更新 `knowledge_tree_with_links.md` 中的 🔗 链接行，补充新教程路径

---

## 阶段一：基础计算实践系列（Tier 1，主题 1–8）

**目标：** 覆盖基础篇最常用的计算类型，每篇约 600–1000 行。

---

### Task 1：平面波收敛性测试（k点 + 截断能）

**来源文件：**
- `data/knowledge_source/ABACUS 的平面波计算与收敛性测试.md`
- `data/knowledge_source/Al元素晶体的自洽迭代计算与平面波收敛测试及k点收敛性测试.md`

**知识点：** 基础篇 > 参数设置 > k点设置 + 平面波截断能

**Step 1：启动新对话，提交任务**
```
主题：ABACUS平面波基组收敛性测试（k点与截断能）
任务类型：C（案例驱动）
案例文件：
  - data/knowledge_source/ABACUS 的平面波计算与收敛性测试.md
  - data/knowledge_source/Al元素晶体的自洽迭代计算与平面波收敛测试及k点收敛性测试.md
特殊要求：两个案例均要展示（Si用于理论说明，Al用于数值对比）
```

**Step 2：跟随 CLAUDE.md Step 0–2，确认大纲**
- 预期大纲包含：物理背景、k点收敛理论、截断能收敛理论、案例一（Si）、案例二（Al，含数值表格）、常见问题

**Step 3：审查 Step 3 初稿中的收敛曲线数值**
- 确认 ecutwfc 扫描数据（20–80 Ry范围）已完整呈现
- 确认 k点网格与能量的数值表格来自 Al 案例文件

**Step 4：完成 Step 4–7，获取最终文件**

**Step 5：更新知识树**
```bash
# 将最终教程复制到 Tutorials 目录
cp "_workspace/<时间戳>/07_Final_Tutorial_*.md" data/knowledge_node/Tutorials/
# 更新 knowledge_tree.md：将 k点设置 ⬜ 和 平面波截断能 ⬜ 改为 ✅
```

**预期产物：** `07_Final_Tutorial_ABACUS平面波基组收敛性测试.md`

---

### Task 2：结构优化（固定晶胞 + 变胞）

**来源文件：**
- `data/knowledge_source/ABACUS 使用教程｜结构优化.md`

**知识点：** 基础篇 > 结构优化 > ionic relax + cell-relax

**Step 1：启动新对话**
```
主题：ABACUS结构优化教程（固定晶胞优化与变胞优化）
任务类型：C（案例驱动）
案例文件：data/knowledge_source/ABACUS 使用教程｜结构优化.md
特殊要求：包含两个案例（MgO ionic relax 和 H₂ cell-relax），对比两种优化方式的参数差异
```

**Step 2：确认大纲包含：** 离子弛豫原理、变胞优化原理、案例一 MgO（ionic）、案例二 H₂（cell-relax）、参数速查表

**Step 3：审查初稿**
- 确认 `relax_method`、`force_thr`、`stress_thr` 等参数值与案例一致
- 确认两个案例的 INPUT 文件均完整呈现

**Step 4：完成审查步骤并获取最终文件**

**预期产物：** `07_Final_Tutorial_ABACUS结构优化教程.md`

---

### Task 3：电子自洽迭代（SCF）

**来源文件：**
- `data/knowledge_source/ABACUS 使用教程｜电子自洽迭代.md`
- `data/knowledge_source/快速开始 ABACUS｜自洽 能带 态密度 结构优化.md`

**知识点：** 基础篇 > 电子结构 > 自洽计算（SCF）

**Step 1：启动新对话**
```
主题：ABACUS电子自洽迭代（SCF）计算教程
任务类型：C（案例驱动）
案例文件：
  - data/knowledge_source/ABACUS 使用教程｜电子自洽迭代.md
  - data/knowledge_source/快速开始 ABACUS｜自洽 能带 态密度 结构优化.md
特殊要求：体现 PW 与 LCAO 两种基组的参数差异（MgO用PW，Ti用LCAO）
```

**Step 2：确认大纲包含：** SCF算法原理、PW基组案例（MgO）、LCAO基组案例（Ti）、收敛判据与混合参数、常见不收敛处理

**Step 3：审查初稿**
- 确认 `smearing_method`、`mixing_type`、`mixing_beta` 等参数值一致
- 确认两种基组的 STRU 文件均完整（含轨道文件名）

**预期产物：** `07_Final_Tutorial_ABACUS电子自洽迭代SCF教程.md`

---

### Task 4：能带结构计算

**来源文件：**
- `data/knowledge_source/快速开始 ABACUS｜自洽 能带 态密度 结构优化.md`
- `data/knowledge_source/2024秋计算材料学-上机练习：ABACUS能带和态密度计算.md`
- `data/knowledge_source/ABACUS计算模拟实例 VIII. 基于HSE06的态密度与能带计算.md`
- `data/knowledge_source/如何正确画能带，NSCF读电荷密度.docx`（如可读）

**知识点：** 基础篇 > 电子结构 > 能带结构计算

**Step 1：启动新对话**
```
主题：ABACUS能带结构计算教程（PBE+HSE06）
任务类型：C（案例驱动）
案例文件：
  - data/knowledge_source/快速开始 ABACUS｜自洽 能带 态密度 结构优化.md
  - data/knowledge_source/2024秋计算材料学-上机练习：ABACUS能带和态密度计算.md
  - data/knowledge_source/ABACUS计算模拟实例 VIII. 基于HSE06的态密度与能带计算.md
特殊要求：包含两步工作流（SCF→NSCF/band），以及 HSE06 案例中的完整 INPUT 三件套
```

**Step 2：确认大纲包含：** 能带理论背景、k路径选择、SCF计算、能带计算（band_num/kspacing）、HSE06进阶案例、后处理可视化

**Step 3：审查初稿**
- 确认 HSE06 案例中 `E_bandgap = 1.19808 eV` 数值出现
- 确认 KPT 高对称点路径与案例一致

**预期产物：** `07_Final_Tutorial_ABACUS能带结构计算教程.md`

---

### Task 5：态密度与PDOS

**来源文件：**
- `data/knowledge_source/ABACUS+Atomkit 计算态密度和能带.md`
- `data/knowledge_source/ABACUS 计算 PDOS.md`

**知识点：** 基础篇 > 电子结构 > 态密度（DOS/PDOS）

**Step 1：启动新对话**
```
主题：ABACUS态密度与分波态密度（DOS/PDOS）计算教程
任务类型：C（案例驱动）
案例文件：
  - data/knowledge_source/ABACUS+Atomkit 计算态密度和能带.md
  - data/knowledge_source/ABACUS 计算 PDOS.md
特殊要求：包含 Al（总DOS）和 Fe（PDOS+投影能带）两个案例，介绍 Atomkit 后处理
```

**Step 2：确认大纲包含：** DOS理论、参数设置（dos_edelta/dos_sigma）、Al案例（总DOS）、Fe案例（PDOS分l角量子数）、Atomkit可视化

**预期产物：** `07_Final_Tutorial_ABACUS态密度与PDOS计算教程.md`

---

### Task 6：分子动力学基础（AIMD + NVT）

**来源文件：**
- `data/knowledge_source/ABACUS 使用教程｜分子动力学.md`
- `data/knowledge_source/ABACUS 分子动力学使用教程.md`

**知识点：** 基础篇 > 分子动力学 > AIMD基本概念 + NVT系综模拟

**Step 1：启动新对话**
```
主题：ABACUS分子动力学基础教程（AIMD与NVT系综）
任务类型：C（案例驱动）
案例文件：
  - data/knowledge_source/ABACUS 使用教程｜分子动力学.md
  - data/knowledge_source/ABACUS 分子动力学使用教程.md
特殊要求：包含 Si（固态）和 He（气态）的 AIMD 对比，完整展示 8 组 INPUT 文件中的代表性参数
```

**Step 2：确认大纲包含：** MD基本原理（NVE/NVT/NPT）、AIMD vs 经典MD、Si晶体案例、He气体案例、输出文件解读（MD_dump）

**Step 3：审查初稿**
- 确认 `md_type`、`md_tfirst`、`md_tlast`、`md_dt` 等参数值与案例文件一致
- 确认 Si 和 He 两个案例的 STRU 文件均完整

**预期产物：** `07_Final_Tutorial_ABACUS分子动力学基础教程.md`

---

### Task 7：分子动力学进阶（NVT/NPT 复杂体系）

**来源文件：**
- `data/knowledge_source/ABACUS 使用教程｜分子动力学.md`（进阶参数部分）
- `data/knowledge_source/ABACUS 分子动力学使用教程.md`

**知识点：** 进阶篇 > 分子动力学进阶 > NVT/NPT系综，复杂体系模拟

**说明：** 与 Task 6 共用来源，但聚焦 NPT 系综、热压参数、大体系注意事项，避免重复基础内容。

**Step 1：启动新对话**
```
主题：ABACUS分子动力学进阶教程（NVT/NPT多系综与参数调优）
任务类型：C（案例驱动）
案例文件：
  - data/knowledge_source/ABACUS 使用教程｜分子动力学.md
  - data/knowledge_source/ABACUS 分子动力学使用教程.md
特殊要求：重点讲 NPT 系综（barostat）参数、多系综对比、大体系 AIMD 注意事项；基础概念引用"分子动力学基础教程"不重复
```

**预期产物：** `07_Final_Tutorial_ABACUS分子动力学进阶教程.md`

---

### Task 8：静电势计算（Electrostatic Potential）

**来源文件：**
- `data/knowledge_source/en__latest__advanced__elec_properties__potential.html.md`
- `data/knowledge_source/abacus_user_guide__abacus_chg.html.md`

**知识点：** 基础篇 > 电荷密度与波函数 > 静电势

**Step 1：启动新对话**
```
主题：ABACUS静电势计算教程
任务类型：C（案例驱动）
案例文件：
  - data/knowledge_source/en__latest__advanced__elec_properties__potential.html.md
  - data/knowledge_source/abacus_user_guide__abacus_chg.html.md
特殊要求：以 Si(111) slab（40原子）为案例，展示完整 STRU；说明 out_pot=2 输出格式；介绍与功函数计算的联系
```

**Step 2：确认大纲包含：** 静电势物理意义、out_pot参数、Si(111)案例（完整STRU+INPUT）、输出文件解读（ElecStaticPot.cube）、可视化

**预期产物：** `07_Final_Tutorial_ABACUS静电势计算教程.md`

---

## 阶段二：表面计算系列（Tier 1，主题 9–14）

每篇均以对应的 `采用 ABACUS 进行表面计算（X）` 文件为案例源，结构相似，可快速连续生成。

| Task | 主题 | 案例源文件 | 核心案例 |
|------|------|----------|---------|
| 9 | 功函数（静电势+真空能级） | `表面计算（一）：静电势和功函数.md` | Al(111) slab |
| 10 | 偶极修正（Dipole correction） | `表面计算（二）：偶极修正.md` | H₂O/graphene slab |
| 11 | 表面能计算 | `表面计算（三）：表面能计算.md` | Al(100) + Si(100) |
| 12 | 表面缺陷能与吸附能 | `表面计算（四）：表面缺陷能和吸附能计算.md` | Mo表面 + Li吸附 |
| 13 | 外加电场 | `表面计算（五）：外加电场.md` | 纳米带 + 电场 |
| 14 | 补偿电荷 | `表面计算（六）：补偿电荷.md` | H₂O slab |

**通用启动模板：**
```
主题：[主题名]
任务类型：C（案例驱动）
案例文件：data/knowledge_source/采用 ABACUS 进行表面计算（X）：[主题].md
特殊要求：完整展示 slab 模型的 STRU 文件；说明与其他表面教程的前置关系
```

---

## 阶段三：大体系与工具类（Tier 1，主题 15–22）

| Task | 主题 | 案例源文件 | 任务类型 | 核心案例 |
|------|------|----------|---------|---------|
| 15 | 晶胞朝向与并行效率优化 | `eff1/2/3.html.md`（3个） | C | CNT + BN + Cu-CO |
| 16 | 数值原子轨道生成（NAO） | `数值原子轨道（二）.md` + `（三）.md` | C | H原子（模守恒赝势+优化）|
| 17 | 能带反折叠（Band Unfolding） | `ABACUS+pyatb 能带反折叠计算.md` | C | GeC超胞（63+1原子）|
| 18 | 介电函数与线性光学（PYATB） | `ABACUS+pyatb：介电函数与线性光学的性质的计算.md` | C | 晶态SiO₂（8+16原子）|
| 19 | DPGEN 机器学习势生成 | `ABACUS+DPGEN 使用教程.md` | C | SiC（3C+2H构型）|
| 20 | ABACUS+LibRI 杂化泛函 | `ABACUS+LibRI 做杂化泛函计算教程.md` | C | Si SCF + H₂O relax |
| 21 | 杂化泛函（PBE0/HSE06） | `ABACUS计算模拟实例 VIII（HSE06）.md` + `平面波杂化泛函.md` | C | Si HSE06能带+DOS |
| 22 | ABACUS+Atomkit DOS后处理 | `ABACUS+Atomkit 计算态密度和能带.md` | C | Al + Fe |

**Task 15 特殊说明：** 3个子案例可合并为一篇"晶胞朝向与并行效率优化"，分三节分别讨论一维/二维/表面体系，或拆成3篇独立教程——建议合并。

**Task 21 特殊说明：** 杂化泛函(🔄)旧教程有错误，Task 20（LibRI）和 Task 21（杂化泛函）内容高度相关，建议先生成 Task 20，Task 21 引用其结论并聚焦带隙修正的物理意义。

---

## 阶段四：Tier 2 框架案例系列（主题 23–44）

按优先级排序，每篇均为类型B（基于案例），生成时需明确告知"此来源文件缺少XX，请补全"。

### 高优先级（建议优先处理）

| Task | 主题 | 案例源文件 | 主要补全内容 |
|------|------|----------|------------|
| 23 | 主要输入文件（INPUT/STRU/KPT） | `ABACUS入门教程 - 结构文件STRU.md` + `如何转换 STRU 文件.md` | 补充完整计算流程 |
| 24 | 非自洽计算（NSCF+能带） | `如何正确画能带，NSCF读电荷密度.docx`（RAG检索）| 补充完整INPUT |
| 30 | 磁性材料（自旋极化+非共线+SOC） | `ABACUS 使用教程｜磁性材料计算.md` | 补充完整STRU/KPT |
| 32 | 铁电极化（Berry相位） | `en__latest__advanced__elec_properties__Berry_phase.html.md` | 补充PbTiO₃的STRU |
| 33 | ABACUS+Wannier90（🔄） | `ABACUS+Wannier90 使用教程.md` | 补充数值结果 |

### 接口工具类

| Task | 主题 | 案例源文件 | 主要补全内容 |
|------|------|----------|------------|
| 34 | DeePMD-kit 势函数训练 | `ABACUS+DeePMD-kit 做机器学习分子动力学模拟.md` | 补充MD结果 |
| 35 | ABACUS+DP 模型 MD | `从 DFT 到 MD｜超详细「深度势能」材料计算上手指南.md` | 补充输出分析 |
| 36 | OFDFT | `ABACUS 无轨道密度泛函理论方法使用教程.md` | 补充完整案例 |
| 37+44 | Candela + MD轨迹分析（合并）| `ABACUS+Candela 使用教程.md` + `abacus_user_guide__abacus_candela.html.md` | 补充ABACUS INPUT |
| 38 | Bader 电荷分析 | `ABACUS+Bader charge 分析教程.md` | 补充完整INPUT |
| 39 | ABACUS+USPEX | `ABACUS+USPEX 接口教程.md` | 补充结构详情 |
| 40 | ABACUS+Hefei NAMD | `ABACUS+Hefei NAMD 使用教程.md` | 补充NAMD结果 |
| 41 | ABACUS+abTEM | `使用ABACUS结合abTEM进行STEM模拟.md` | 补充ABACUS参数 |
| 42 | ASE-ABACUS使用方法简介 | `ASE-ABACUS 第一章：使用方法简介.md` | 补充完整案例 |
| 43 | 单端过渡态搜索（Dimer/ART） | `ASE-ABACUS 第三章：单端过渡态搜索.md` | 补充具体体系 |

### 其他Tier2主题

| Task | 主题 | 来源 | 说明 |
|------|------|------|------|
| 25 | 波函数可视化 | `abacus_user_guide__abacus_chg.html.md` | 与静电势教程(Task 8)配合 |
| 26 | 电荷密度分布 | 同上 | 与Task 25合并可减少重复 |
| 27 | 赝势选择 | `ABACUS赝势-轨道捆绑包APNS-PPORB-v1发布...md` | 高通量数据，以概念+表格为主 |
| 28 | SCF收敛算法参数 | `ABACUS收敛性问题解决手册.docx` + RAG | docx不可直接解析，依赖RAG |
| 29 | Bohrium平台使用 | `Bohrium 帮助文档｜ABACUS.md` | 平台操作类，时效性强 |
| 31 | 磁交换参数 | `ABACUS计算磁性相互作用参数.md` + `ABACUS磁性交换作用的计算以及与微磁学的结合.md` | 多工具耦合 |

---

## 阶段五：安装/环境类（A1–A9，按需处理）

| Task | 主题 | 关键来源 | 优先级 |
|------|------|---------|------|
| A5 | GPU 编译方法 | `ABACUS LCAO基组 GPU版本使用介绍.md` | 相对最高（有使用说明）|
| A7 | Toolchain 安装（5种）| Bohrium脚本 + 微信5篇 | 次优（资料丰富）|
| A4 | CPU 编译（GCC/Intel）| `abacus_user_guide__abacus_gcc.html.md` 等 | 中 |
| A2 | 安装与编译（官方）| 官方文档 | 中 |
| A9 | 常见报错 | 官方 FAQ | 中 |
| A1 | 架构介绍 | ABACUS论文 + 微信文章 | 低 |
| A3 | 依赖库 | `abacus_user_guide__abacus_hpc.html.md` | 低 |
| A6 | DCU编译 | `abacus_user_guide__abacus_dcu.html.md` | 低 |
| A8 | 无MPI编译 | 飞书文档 | 低 |

---

## 阶段六：Tier 3 文档/概念类（45–53，暂缓）

待 Tier 1 和 Tier 2 高优先级完成后，根据需要处理。以类型A（无案例）生成，建议先确认 RAG 检索质量再启动。

---

## 知识树更新操作（每篇完成后执行）

**Step 1：复制最终教程**
```bash
cp "_workspace/<时间戳_主题>/07_Final_Tutorial_<标题>.md" \
   "data/knowledge_node/Tutorials/"
```

**Step 2：更新 knowledge_tree.md**
- 找到对应知识点行，将 `⬜` 或 `🔄` 改为 `⬜ ✅`（保留旧状态标记）
- 在"新教程对应关系"表追加一行：
```
| <教程标题> | <知识点位置> |
```
- 更新覆盖统计数字

**Step 3：更新 knowledge_tree_with_links.md**
- 在对应知识点下追加 🔗 链接行：
```
  - 🔗 [原有链接] | [新教程](data/knowledge_node/Tutorials/<文件名>.md)
```
- 在新教程对应关系表追加文件路径列

**Step 4：提交**
```bash
git add data/knowledge_node/Tutorials/<文件名>.md \
        data/knowledge_node/knowledge_tree.md \
        data/knowledge_node/knowledge_tree_with_links.md
git commit -m "feat(tutorials): 新增 <主题> 教程"
```

---

## 进度跟踪

建议在此文件末尾维护生成进度：

| # | 主题 | 状态 | 工作目录 | 完成日期 |
|---|------|------|---------|---------|
| 1 | 平面波收敛性测试 | ⬜ 待开始 | — | — |
| 2 | 结构优化 | ⬜ 待开始 | — | — |
| 3 | 电子自洽迭代（SCF）| ⬜ 待开始 | — | — |
| 4 | 能带结构计算 | ⬜ 待开始 | — | — |
| 5 | 态密度与PDOS | ⬜ 待开始 | — | — |
| 6 | 分子动力学基础 | ⬜ 待开始 | — | — |
| 7 | 分子动力学进阶 | ⬜ 待开始 | — | — |
| 8 | 静电势 | ⬜ 待开始 | — | — |
| 9 | 功函数 | ⬜ 待开始 | — | — |
| 10 | 偶极修正 | ⬜ 待开始 | — | — |
| 11 | 表面能计算 | ⬜ 待开始 | — | — |
| 12 | 表面缺陷能与吸附能 | ⬜ 待开始 | — | — |
| 13 | 外加电场 | ⬜ 待开始 | — | — |
| 14 | 补偿电荷 | ⬜ 待开始 | — | — |
| 15 | 晶胞朝向与并行效率 | ⬜ 待开始 | — | — |
| 16 | NAO生成 | ⬜ 待开始 | — | — |
| 17 | 能带反折叠 | ⬜ 待开始 | — | — |
| 18 | 介电函数与线性光学 | ⬜ 待开始 | — | — |
| 19 | DPGEN | ⬜ 待开始 | — | — |
| 20 | ABACUS+LibRI杂化泛函 | ⬜ 待开始 | — | — |
| 21 | 杂化泛函（HSE06）| ⬜ 待开始 | — | — |
| 22 | ABACUS+Atomkit | ⬜ 待开始 | — | — |
| 23–44 | Tier 2 系列 | ⬜ 待开始 | — | — |
| A1–A9 | 安装/环境类 | ⬜ 待开始 | — | — |
| 45–53 | Tier 3 概念类 | ⬜ 暂缓 | — | — |
