# 任务简报

## 基本信息
- **任务时间**: 2026-03-30 17:03:42
- **教程主题**: ABACUS 态密度（DOS/PDOS）计算
- **任务类型**: 类型B - 基于案例的教程
- **系列位置**: 第4篇（共4篇系列教程）

## 输入材料
- **案例文件**: `data\input\dos_MgO.md`
  - 材料体系: MgO
  - 计算类型: DOS/PDOS
  - 基组类型: LCAO
  - 计算流程: SCF自洽计算
- **参考资料**: `docs\case\331力学_电子结构\` 文件夹
- **系列规划**: `docs\case\331力学_电子结构\新的要求.md`

## 系列教程背景
前3篇教程内容：
1. ABACUS基本介绍（软件背景、文件格式、结构优化、输出解读）
2. 力学性质计算（弹性模量原理、abacustest前后处理）
3. 能带结构计算（能带知识、ABACUS计算、abacustest作图）

本篇（第4篇）需要避免与前3篇重复，重点关注DOS/PDOS特有内容。

## 教程目标
1. DOS/PDOS基本知识介绍
2. 如何使用ABACUS输出DOS/PDOS（INPUT参数设置）
3. 如何使用abacustest对结果作图

## 案例关键信息
- **INPUT关键参数**:
  - `calculation = scf`
  - `basis_type = lcao`
  - `out_dos = 2` (输出总DOS和分波态密度PDOS)
  - `init_chg = file` (从文件读取初始电荷密度)
  - `smearing_method = gauss`, `smearing_sigma = 0.015`

- **STRU文件**:
  - MgO岩盐结构
  - 轨道文件: `Mg_gga_10au_100Ry_2s1p.orb`, `O_gga_6au_100Ry_2s2p1d.orb`

- **KPT文件**:
  - Gamma中心, 18×18×18 k点网格

## 特殊要求
- 详细程度：尽量详细
- 内容规划：避免与前3篇重复，可以提及但不展开
- 风格：参考 `data/reference_materials/style_references/` 中的文章

## 执行计划
按照CLAUDE.md的7步流程：
1. ✅ Step 0: 任务初始化（当前步骤）
2. Step 1: 知识调研与案例解析
3. Step 2: 大纲讨论（提供3个方案）
4. Step 3: 撰写完整初稿
5. Step 4: 内容审查
6. Step 5: 案例审查
7. Step 6: 风格审查
8. Step 7: 最终输出
