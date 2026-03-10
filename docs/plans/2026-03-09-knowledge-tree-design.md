# ABACUS 知识树设计文档

**日期：** 2026-03-09
**任务：** 梳理 knowledge_node.md 中所有知识点，生成带状态标签的树状结构文件

---

## 输入

- `data/knowledge_node/knowledge_node.md` — 人工整理的 ABACUS 知识点文档（含旧教程链接）
- `data/knowledge_node/Tutorials/` — 当前项目新生成的 11 篇教程
- `data/knowledge_node/knowledge_node.md` 末尾"迭代中的教程列表" — 含错误、正在更新的 4 篇旧教程

## 输出

- `data/knowledge_node/knowledge_tree.md` — 带状态标签的知识树 Markdown 文件

---

## 图例系统

| 符号 | 含义 |
|------|------|
| ✅   | **新文章已生成** — `data/knowledge_node/Tutorials/` 中已有对应新教程 |
| ⬜   | **有旧教程** — knowledge_node.md 中附有外部链接的旧教程 |
| 🔄   | **迭代中旧教程** — 飞书"迭代中教程列表"，含错误、正在更新 |
| ❌   | **缺少教程** — 无任何教程，尚待编写 |

一个知识点可同时标注多个符号（如 `🔄 ✅` = 有迭代中旧教程且新文章已生成）。

---

## 知识树结构（设计稿）

### 一、基础篇

**1. 软件与环境**
- 架构介绍 ⬜
- 安装与编译 ⬜
  - 依赖库（MPI、LibXC、ELPA 等）⬜
  - CPU 编译 ⬜
  - GPU 编译 ⬜
  - DCU 编译 ⬜
  - Toolchain 安装 ⬜
- 常见报错与解决方法 ⬜

**2. 输入输出体系**
- 主要输入文件（INPUT / STRU / KPT）⬜
- 运行输出解读 ⬜
- 参数设置经验
  - k 点设置 ⬜
  - 平面波截断能 ⬜
  - 赝势选择 ⬜
  - 基组选择（PW vs NAO）⬜
    - SZ / DZP / TZDP 基组差异 ⬜
    - 基组收敛测试 ❌
    - 基组与赝势匹配原则 ❌
    - 数值原子轨道生成 ⬜
  - SCF 收敛算法参数 ⬜

**3. 基础计算实践**
- 结构优化（固定晶胞 / 变胞）⬜
- 电子结构：SCF ⬜ / NSCF ⬜ / 能带 ⬜ / DOS ⬜
- 分子动力学：AIMD 概念 ❌ / NVT 系综 ⬜ / 轨迹分析 ⬜
- 电荷密度与波函数：电荷分布 ⬜ / Mulliken ⬜ / 静电势 ⬜ / 波函数 ⬜ / ELF ⬜ ✅
- Bohrium 平台运用 ⬜

### 二、进阶篇

**1. 晶体材料**
- 对称性分析 ❌ / 结合能 ❌
- 带隙计算：PBE ❌ / DFT+U 🔄 ✅ / 杂化泛函 🔄 ⬜ / DFT+1/2 ❌
- 拓扑计算：基础 ⬜ / PYATB ⬜ / Wannier90 🔄 ⬜
- 磁性计算 ⬜：自旋极化 ❌ / 非共线自旋 ❌ / SOC ❌ / 磁交换参数 ❌
- 铁电极化（Berry 相位）⬜
- 光学性质 ⬜ ✅：光吸收谱 ✅ / PYATB ⬜ / GW 方法 ❌
- RT-TDDFT ⬜ ✅：混合规范 ✅ / 外场激发 ❌ / 超快动力学 ❌
- 声子谱与热力学性质 ⬜ ✅
- 机器学习势函数：DPGEN ⬜ / DeePMD ⬜ / 微调蒸馏 ❌ / ABACUS+DP MD ⬜

**2. 表面与界面材料**
- 表面建模 ❌ / 功函数 ⬜ / 偶极修正 ⬜ / 表面能 ⬜
- 表面缺陷能与吸附能 ⬜ / 外加电场 ⬜ / 补偿电荷 ⬜
- NEB 过渡态（ATST-Tools）⬜ ✅ / 层错能 ❌ / 隐式溶剂模型 ⬜ ✅

**3. 大体系材料计算**
- 缺陷与掺杂：点缺陷建模 ❌ / 形成能 ❌ / 宽禁带半导体 ❌
- SDFT ⬜ ✅ / Band Unfolding ⬜ / OFDFT ⬜
- 晶格热导率：Phonopy ⬜ ✅ / ShengBTE ⬜ ✅ / Phono3py ⬜

**4. 预/后处理：接口与工具**
- abacustest 🔄 ✅（快速准备输入 ✅ / 弹性常数 🔄 ✅）
- Atomkit ⬜ / LibRI（LCAO 杂化泛函）🔄 ⬜ / Candela ⬜
- Bader Charge ⬜ / USPEX ⬜ / Hefei NAMD ⬜ / abTEM ⬜
- ASE-ABACUS ⬜（简介 ⬜ / NEB ✅ / 单端过渡态 ⬜）
- DPA-3 ⬜

---

## 新教程映射关系（11 篇）

| 新教程文件名 | 对应知识点 |
|------------|-----------|
| 使用abacustest计算晶体弹性性质 | 预/后处理 > abacustest > 弹性常数 |
| ABACUS_DFT+U强关联体系计算 | 晶体材料 > 带隙计算 > DFT+U |
| ABACUS隐式溶剂模型使用 | 表面与界面 > 隐式溶剂模型 |
| 使用ABACUS+Phonopy计算声子谱 | 晶体材料 > 声子谱；大体系 > 晶格热导率 > Phonopy |
| 使用ABACUS+ShengBTE计算Si晶格热导率 | 大体系 > 晶格热导率 > ShengBTE |
| ABACUS随机密度泛函理论SDFT使用教程 | 大体系 > SDFT |
| 混合规范RTTDDFT | 晶体材料 > RT-TDDFT > 混合规范 |
| ABACUS使用教程-ELF电子局域函数计算与可视化 | 基础 > 电荷密度与波函数 > ELF |
| 使用ABACUS和ATST-Tools进行NEB-AutoNEB过渡态搜索 | 表面与界面 > NEB；预/后处理 > ASE-ABACUS > NEB |
| 使用abacustest快速准备ABACUS输入文件 | 预/后处理 > abacustest > 快速准备输入 |
| 光学性质计算 | 晶体材料 > 光学性质 |

---

## 迭代中旧教程映射（4 篇）

| 飞书标题 | 对应知识点 | 新文章？ |
|---------|-----------|---------|
| DFT+U 计算 | 晶体材料 > DFT+U | ✅ 已生成 |
| ABACUS+Wannier90 | 晶体材料 > 拓扑 > Wannier90 | ❌ 未生成 |
| LCAO 原子轨道基组杂化泛函 | 晶体材料 > 杂化泛函；预/后处理 > LibRI | ❌ 未生成 |
| 弹性常数 | 预/后处理 > abacustest > 弹性常数 | ✅ 已生成 |
