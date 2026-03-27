# Step 1: 知识调研与案例解析

## 1.1 RAG 知识检索结果

### 检索1：ABACUS 软件背景

**关键信息：**
- ABACUS（Atomic-orbital Based Ab-initio Computation at UStc）中文名"原子算筹"
- 国产开源密度泛函理论（DFT）软件
- 支持平面波（PW）和数值原子轨道（NAO/LCAO）两种基组
- GitHub 仓库：https://github.com/deepmodeling/abacus-develop
- 官方文档：https://abacus.deepmodeling.com/en/latest/
- 主要功能：电子结构自洽迭代、原子结构优化、分子动力学计算
- 支持 LDA、GGA、meta-GGA、hybrid functionals
- 高级功能：DFT+U、VdW 修正、隐式溶剂模型等

### 检索2：INPUT 文件参数

**核心参数：**
- `suffix`：系统名称，默认为 ABACUS
- `ntype`：单位晶胞中元素的种类数目
- `pseudo_dir`：赝势文件目录
- `orbital_dir`：轨道文件目录（LCAO 基组需要）
- `ecutwfc`：波函数展开的平面波能量截止（单位：Rydberg）
- `scf_thr`：电荷密度收敛阈值（单位：Rydberg）
- `basis_type`：基函数类型（pw 或 lcao）
- `calculation`：计算类型（scf、relax、cell-relax、md 等）
- `symmetry`：是否考虑对称性（1/0/-1）
- `smearing_method`：展宽方法（gauss、mp 等）
- `smearing_sigma`：展宽参数

**重要说明：**
- INPUT 参数列表以 `INPUT_PARAMETERS` 开头
- 以 `#` 或 `/` 开头的行被忽略
- 布尔参数可用 True/False、1/0、T/F 设置（大小写不敏感）
- 文件名"INPUT"不可更改
- 参数设置原则：简洁明了，默认值通常已是好选择

### 检索3：STRU 文件格式

**文件结构：**
```
ATOMIC_SPECIES
元素名 原子质量 赝势文件

NUMERICAL_ORBITAL（LCAO 基组需要）
轨道文件名

LATTICE_CONSTANT
晶格常量（缩放系数）

LATTICE_VECTORS
晶格矢量（3行，每行3个数）

ATOMIC_POSITIONS
坐标系类型（Direct/Cartesian）
元素名
磁矩初值
原子数
原子坐标 移动标志(1/0) 磁矩
```

**关键说明：**
- 赝势和轨道文件可从 http://abacus.ustc.edu.cn/pseudo/list.htm 下载
- LATTICE_CONSTANT 单位：1.8897259886 Bohr = 1.0 Angstrom
- Direct 坐标：晶格坐标系；Cartesian：笛卡尔坐标系
- 移动标志：1 表示可移动，0 表示固定

### 检索4：KPT 文件格式

**文件结构：**
```
K_POINTS
0           # k点总数，0表示自动生成
Gamma       # Monkhorst-Pack方法类型（Gamma或MP）
4 4 4 0 0 0 # 前3个数：各方向细分数；后3个数：网格偏移
```

**关键说明：**
- Gamma：以 Gamma 点为中心的 Monkhorst-Pack 方法
- MP：标准 Monkhorst-Pack 方法
- 立方晶格各方向 k 点数应相同
- 大体系可设置 `gamma_only = 1` 只用 Gamma 点计算

### 检索5：输出文件

**主要输出：**
- `OUT.suffix/`：输出文件夹（suffix 在 INPUT 中设置）
- `OUT.suffix/INPUT`：包含所有实际使用的参数（用户指定+系统默认）
- `OUT.suffix/running_scf.log`：包含几乎所有函数调用信息
- 屏幕输出：实时显示计算进度和关键信息

### 检索6：结构优化（relax/cell-relax）

**计算类型：**
- `relax`：固定晶胞，优化原子位置
- `cell-relax`：同时优化晶胞参数和原子位置

**关键参数：**
- `force_thr_ev`：原子受力收敛阈值（单位：eV/Angstrom）
- `stress_thr`：应力收敛阈值（单位：kBar）
- `relax_nmax`：最大离子迭代步数
- `out_stru`：是否输出优化后的结构（1/0）

**收敛判据：**
- relax：所有原子最大受力 < force_thr_ev
- cell-relax：最大受力 < force_thr_ev 且 总应力 < stress_thr

### 检索7：abacustest 工具

**功能：**
- 准备 ABACUS 输入文件
- 高通量任务计算
- 后处理和结果提取

**安装：**
```bash
# 通过 pip 安装
pip install abacustest

# 或从源码安装
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

**使用：**
```bash
# 准备输入文件
abacustest model inputs -h

# 查看帮助
abacustest -h
```

**检索质量评估：**
- ✅ 检索结果覆盖了所有9个核心内容点
- ✅ 参数说明详细，包含默认值和单位
- ✅ 文件格式清晰，有完整示例
- ⚠️ abacustest 具体使用细节需要补充
- ⚠️ 需要实际案例来演示完整流程

---

## 1.2 参考资料阅读

### 已有教程覆盖的知识点（避免重复）

根据 `知识点覆盖报告.md`，已有三篇教程：

**教程1：输入输出体系**
- STRU/KPT/INPUT 三文件详解（已覆盖）
- SCF 收敛算法参数（mixing 系列）
- ecutwfc 收敛测试
- 运行输出解读

**教程2：力学性质**
- 几何优化（relax/cell-relax）
- 弹性常数计算
- abacustest 弹性流程

**教程3：电子结构**
- 能带结构计算
- 态密度（DOS/PDOS）
- abacus-plot 和 abacustest 后处理

**本篇需要避免的内容：**
- ❌ 不深入讲解弹性常数计算（属于教程2）
- ❌ 不深入讲解能带/DOS 计算（属于教程3）
- ❌ 不讲解 abacus-plot 工具（属于教程3）
- ✅ 可以简单提及 relax 作为基础功能
- ✅ 可以介绍 abacustest 准备输入文件和提取基本结果

### abacustest 工具使用

**功能：**
- `abacustest model inputs`：准备输入文件
- `abacustest model dos-pdos`：DOS/PDOS 后处理（但本篇不深入）

**dos-pdos 命令参数：**
- `-j JOB`：ABACUS 输入路径
- `--range`：能量范围（默认 -10,10 eV）
- `--plot-type`：投影类型（species/shell/orbital/atom）
- `--suffix`：输出文件后缀

---

## 1.3 风格参考学习

读取了 `ABACUS 使用教程｜磁性材料计算.md` 前100行，提取风格特征：

**开头方式：**
- 直接进入主题，简洁介绍背景
- 例："材料的磁性对应的微观电子结构为电子自旋的方向并不均匀分布"

**结构特点：**
- 使用 ## 和 ### 层级标题
- 表格呈现参数对比（清晰直观）
- 代码块包含注释说明

**语言特征：**
- 简洁直接，不用"在当今"、"值得注意的是"等AI腔
- 技术术语准确（共线磁矩、泡利矩阵等）
- 用"该"、"其"等指代词保持流畅

**代码呈现：**
- 代码块标注语言类型（```INPUT）
- 行内注释说明关键参数
- 完整示例（不省略）

---

## Think Aloud 总结

**检索质量：** ✅ 优秀
- 覆盖所有9个核心内容点
- 参数说明详细（单位、默认值、推荐值）
- 文件格式清晰（完整示例）

**内容边界：** ✅ 明确
- 已有教程覆盖：ecutwfc 收敛、mixing 参数、弹性计算、能带/DOS
- 本篇重点：软件背景、三文件基础、基本参数分类、relax 入门、输出解读、abacustest 基础

**风格特点：** ✅ 已掌握
- 简洁直接，技术准确
- 表格+代码块结合
- 避免AI腔

**下一步：** 进入 Step 2，设计3个大纲方案供用户选择
