# 任务简报

## 任务信息
- **主题：** ABACUS 力学性质
- **任务类型：** B（知识主题 + 已验证案例数据）
- **案例来源：** 数据库 + 联网检索 + 知识库现有教程（无需新跑）
- **参考简报：** docs/case/力学_电子结构/主题知识点覆盖简报.md

## 必须覆盖的内容
1. 弹性常数张量（Cij，Voigt 表示）
2. 体弹/剪切/杨氏模量/泊松比
3. abacustest 自动化弹性流程（job.json + pymatgen）
4. 结构优化（ionic relax）+ 变胞优化（cell-relax）+ 收敛判据

## 核心案例（已验证，无需新跑）
- **Si 弹性常数**（ElasticPlugin 验证）：C11=165.5, C12=58.6, C44=82.5 GPa；BV=94.19, GV=70.89, EV=170.02 GPa, ν=0.199
- **BN cell-relax**（RelaxPlugin 验证）：96原子六方BN体系，h-BN层状结构

## 已有知识来源
- `data/knowledge_node/Tutorials/07_Final_Tutorial_使用abacustest计算晶体弹性性质.md` — 完整abacustest弹性流程（Si + TiO₂）
- `ABACUS+pymatgen 计算弹性常数.md` — pymatgen方式的弹性常数计算
- `abacus_user_guide__abacus_elastic.html.md` — 官方手册弹性常数章节
- `ABACUS 使用教程｜结构优化.md` — relax/cell-relax完整教程
- `abacus_user_guide__abacus_eff2.html.md` — BN cell-relax 96原子案例（完整INPUT/STRU）
- `en__latest__advanced__opt.html.md` — 官方几何优化文档

## 写作方向
本教程面向有ABACUS基础的研究者，围绕"从结构优化到弹性性质计算"的完整计算链条展开：
先讲结构优化（relax/cell-relax）的原理与参数，再介绍弹性常数的物理意义，
最后展示用abacustest/pymatgen自动化完成整个工作流。

## 预期长度控制
- 参考文章平均：345行、883词
- 本教程目标：600-900行（中等长度）
