# ABACUS 知识树生成 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 生成 `data/knowledge_node/knowledge_tree.md`，将 knowledge_node.md 中所有知识点整理成带状态标签的树状结构。

**Architecture:** 纯内容生成任务。读取 knowledge_node.md 和 Tutorials/ 目录，按设计文档中的树状结构和图例，手动写出完整 Markdown 知识树文件。

**Tech Stack:** Markdown（无代码依赖）

**设计参考：** `docs/plans/2026-03-09-knowledge-tree-design.md`

---

## 图例（必须包含在文件头部）

```
✅ 新文章已生成   — data/knowledge_node/Tutorials/ 中已有对应新教程
⬜ 有旧教程       — knowledge_node.md 中附有外部链接的旧教程
🔄 迭代中旧教程   — 飞书"迭代中教程列表"，含错误、正在更新
❌ 缺少教程       — 无任何教程，尚待编写
（知识点可同时标注多个符号）
```

---

### Task 1: 生成知识树文件

**Files:**
- Create: `data/knowledge_node/knowledge_tree.md`

**Step 1: 写入文件头部（图例 + 说明）**

文件开头包含：
- 标题、生成日期
- 图例表格
- 知识点计数统计（✅ N个 / ⬜ N个 / 🔄 N个 / ❌ N个）

**Step 2: 写入一、基础篇**

按以下结构写入（使用缩进列表表示层级）：

```markdown
## 一、基础篇：安装与基本计算

### 1. 软件与环境

- 架构介绍（发展历史与定位）⬜
- 安装与编译 ⬜
  - 依赖库（MPI、LibXC、ELPA 等）⬜
  - CPU 编译方法 ⬜
  - GPU 编译方法 ⬜
  - DCU 编译方法 ⬜
  - Toolchain 安装教程 ⬜
  - 无 MPI 编译 ⬜
- 常见报错与解决方法 ⬜

### 2. 输入输出体系

- 主要输入文件：INPUT / STRU / KPT ⬜
- 运行输出解读（能量收敛、电子步、力与应力）⬜
- 参数设置经验
  - k 点设置（Monkhorst-Pack 网格、收敛性测试）⬜
  - 平面波截断能 ⬜
  - 赝势选择（ONCV vs HGH，元素库差异）⬜
  - 基组选择（PW vs NAO）⬜
    - SZ / DZP / TZDP 等基组差异 ⬜
    - 基组收敛测试 ❌
    - 基组与赝势匹配原则 ❌
    - 数值原子轨道生成 ⬜
  - SCF 收敛算法相关参数 ⬜

### 3. 基础计算实践

- 结构优化 ⬜
  - 固定晶胞优化 ⬜
  - 变胞优化 ⬜
- 电子结构
  - 自洽计算（SCF）⬜
  - 非自洽计算（NSCF）⬜
  - 能带结构计算 ⬜
  - 态密度（DOS / PDOS）⬜
- 分子动力学（MD）
  - AIMD 基本概念 ❌
  - NVT 系综模拟 ⬜
  - 输出文件与轨迹分析 ⬜
- 电荷密度与波函数
  - 电荷密度分布 ⬜
  - Mulliken 电荷分析 ⬜
  - 静电势 ⬜
  - 波函数可视化 ⬜
  - 电子局域函数（ELF）⬜ ✅
- Bohrium 平台运用 ⬜
```

**Step 3: 写入二、进阶篇 — 1. 晶体材料**

```markdown
## 二、进阶篇：高级功能与物性计算

### 1. 晶体材料

- 对称性分析（周期性与对称性、点阵、点群、空间群）❌
- 结合能（能量最低构型搜索）❌
- 带隙计算
  - 基础 PBE 带隙计算流程 ❌
  - DFT+U 带隙修正 🔄 ✅
  - 杂化泛函（Hybrid functional）🔄 ⬜
  - DFT+1/2 ❌
- 拓扑计算
  - 材料拓扑性质基础知识 ⬜
  - ABACUS+PYATB 计算流程 ⬜
  - ABACUS+Wannier90+WannierTools 🔄 ⬜
- 磁性计算 ⬜
  - 自旋极化 ❌
  - 非共线自旋 ❌
  - SOC 效应 ❌
  - 磁交换参数与各向异性能量 ❌
- 铁电极化（Berry 相位方法）⬜
- 光学性质分析 ⬜ ✅
  - 光吸收谱 ✅
  - ABACUS+PYATB 介电函数与线性光学 ⬜
  - GW 方法 ❌
- 实时含时密度泛函（RT-TDDFT）⬜ ✅
  - 速度规范 / 长度规范 / 混合规范 ✅
  - 外场激发（泵浦-探测）❌
  - 超快动力学 ❌
- 声子谱与热力学性质 ⬜ ✅
- 机器学习势函数
  - DPGEN ⬜
  - DeePMD-kit ⬜
  - 模型微调与蒸馏 ❌
  - ABACUS+DP 模型分子动力学 ⬜
```

**Step 4: 写入进阶篇 — 2. 表面与界面材料**

```markdown
### 2. 表面与界面材料

- 表面建模（真空层的影响）❌
- 功函数 ⬜
- 偶极修正 ⬜
- 表面能计算 ⬜
- 表面缺陷能与吸附能 ⬜
- 外加电场 ⬜
- 补偿电荷 ⬜
- NEB 过渡态搜索（ATST-Tools）⬜ ✅
- 层错能 ❌
- 隐式溶剂模型 ⬜ ✅
```

**Step 5: 写入进阶篇 — 3. 大体系材料计算**

```markdown
### 3. 大体系材料计算

- 缺陷与掺杂
  - 点缺陷建模（超胞方法）❌
  - 形成能与电荷态转变 ❌
  - 宽禁带半导体缺陷调控 ❌
- 随机波函数 DFT（SDFT）⬜ ✅
- 能带反折叠（Band Unfolding）⬜
- 无轨道密度泛函理论（OFDFT）⬜
- 晶格热导率
  - ABACUS+Phonopy ⬜ ✅
  - ABACUS+ShengBTE ⬜ ✅
  - ABACUS+Phono3py ⬜
```

**Step 6: 写入进阶篇 — 4. 预/后处理**

```markdown
### 4. 预/后处理：接口与工具

- abacustest 🔄 ✅
  - 快速准备 ABACUS 输入文件 ✅
  - 弹性常数计算 🔄 ✅
- ABACUS+Atomkit（DOS 与能带后处理）⬜
- ABACUS+LibRI（LCAO 基组杂化泛函）🔄 ⬜
- ABACUS+Candela（MD 轨迹分析）⬜
- ABACUS+Bader Charge 分析 ⬜
- ABACUS+pymatgen（弹性常数）⬜
- ABACUS+USPEX（晶体结构预测）⬜
- ABACUS+Hefei NAMD（激发态动力学）⬜
- ABACUS+abTEM（STEM 模拟）⬜
- ASE-ABACUS ⬜
  - 使用方法简介 ⬜
  - NEB 过渡态计算 ✅
  - 单端过渡态搜索 ⬜
- ABACUS+DPA-3 ⬜
- ABACUS+DPGEN ⬜
```

**Step 7: 在文件头部填写统计数量**

统计各类标签数量，填入文件头部的统计表格：
- 统计树中所有叶节点（最底层知识点）的标签分布

**Step 8: Commit**

```bash
git add data/knowledge_node/knowledge_tree.md docs/plans/2026-03-09-knowledge-tree-design.md docs/plans/2026-03-09-knowledge-tree-plan.md
git commit -m "feat(knowledge): 生成 ABACUS 知识树，标注教程覆盖状态"
```

---

## 验证方法

打开 `data/knowledge_node/knowledge_tree.md`，检查：
1. 图例在文件顶部清晰可见
2. 树状结构层级正确（使用缩进列表）
3. 每个知识点都有且仅有正确的状态标签
4. 统计数字与实际节点数一致
