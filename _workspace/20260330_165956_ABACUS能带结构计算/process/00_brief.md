# 任务简报

## 任务信息
- **主题**：ABACUS 能带结构计算
- **任务类型**：类型B（基于案例的教程）
- **案例文件**：data\input\band_Si_diamond.md
- **参考资料**：docs\case\331力学_电子结构\
- **创建时间**：2026-03-30 16:59:56
- **工作目录**：_workspace/20260330_165956_ABACUS能带结构计算/

## 任务背景
这是4篇系列教程中的第3篇，前两篇分别为：
1. ABACUS基本介绍（INPUT/STRU/KPT文件格式、结构优化、输出文件解读）
2. 力学性质计算（弹性模量计算原理、abacustest前后处理）

本篇聚焦能带结构计算，后续第4篇为DOS/PDOS计算。

## 任务要求
1. 内容详细但不重复前两篇已讲内容（可简单提及）
2. 不提前讲第4篇DOS/PDOS的内容
3. 使用Si金刚石结构案例（band_Si_diamond.md）
4. 参考docs\case\331力学_电子结构\中的资料
5. 使用abacustest进行结果作图

## 案例概述
- **材料体系**：Si（金刚石结构，8原子）
- **计算类型**：NSCF（非自洽计算）
- **基组类型**：平面波（pw）
- **关键参数**：
  - calculation = nscf
  - symmetry = 0
  - init_chg = file（读入电荷密度）
  - out_band = 1
  - kpoint_file = KPT.nscf（高对称路径）

## 执行计划
按照CLAUDE.md的7步流程：
1. Step 1：知识调研与案例解析
2. Step 2：大纲讨论（提供3个方案）
3. Step 3：撰写完整初稿
4. Step 4：内容审查
5. Step 5：案例审查
6. Step 6：风格审查
7. Step 7：最终输出

## Think Aloud
我对任务的理解：
- 这是系列教程的第3篇，需要承上启下
- 前两篇已讲过INPUT/STRU/KPT基础，本篇重点在能带计算的特殊设置
- 核心是SCF→NSCF两步流程，以及高对称路径KPT的设置
- 案例是Si间接带隙半导体，PBE带隙约0.57 eV（不是实验值1.17 eV）
- 需要介绍abacustest作图工具
- 避免与第4篇DOS/PDOS内容重叠
