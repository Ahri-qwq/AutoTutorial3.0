# 大纲方案

## 选定方案：C（问题导向型）

**总体策略：** 以新手常见问题为引导，逐一解答，降低门槛。叙述自然，重点突出。

---

## 最终大纲

### 第一章：ABACUS 是什么（~40行）
- 软件来源：国产开源，USTC → DeepModeling
- 核心能力：DFT、PW + LCAO 双基组
- 适用场景：晶体、表面、分子
- 获取途径：GitHub、文档地址

### 第二章：运行一次计算需要什么（~30行）
- 三个必需文件概述（INPUT / KPT / STRU）
- 赝势文件和轨道文件的角色
- 文件目录结构示例

### 第三章：INPUT 文件详解（~220行）
- 格式规范（INPUT_PARAMETERS 关键词、注释、布尔值）
- 通用参数（suffix、ntype、calculation、basis_type、ecutwfc、pseudo_dir、orbital_dir）
- SCF 收敛参数（scf_thr、scf_nmax）
- Smearing 参数（smearing_method、smearing_sigma）
- Mixing 参数（mixing_type、mixing_beta、mixing_ndim）
- 结构计算参数（cal_force、cal_stress、force_thr_ev、stress_thr、relax_nmax）
- 输出控制参数（out_stru、out_chg 等）
- 完整 INPUT 示例（Si SCF）

### 第四章：STRU 文件详解（~100行）
- 格式规范（各块说明）
- ATOMIC_SPECIES / NUMERICAL_ORBITAL / LATTICE_CONSTANT / LATTICE_VECTORS / ATOMIC_POSITIONS
- Direct vs Cartesian 坐标系
- 移动标志（0/1）的含义
- 完整 STRU 示例（Si，PW 和 LCAO 两种写法）

### 第五章：KPT 文件详解（~80行）
- 格式规范
- 自动生成模式（Gamma/MP）
- kspacing 替代方案
- gamma_only 参数
- 能带计算 line mode（简要提及，不展开，指向后续教程）

### 第六章：计算结束后看什么（~100行）
- OUT.suffix/ 目录结构
- 屏幕输出解读（版本信息、SCF 迭代行格式、收敛判断）
- running_scf.log 五段结构
- 常用 grep 命令（FINAL_ETOT、DRHO、E_bandgap）
- 其他常用输出文件简介

### 第七章：用 abacustest 提高效率（~120行）
- abacustest 简介与安装
- model inputs prepare：快速准备输入文件
  - 从结构文件生成 STRU
  - 指定赝势和轨道文件
  - 指定 INPUT 模板
- model inputs post：从输出中抓取关键结果
  - 抓取 etot、band_gap 等量
- 完整使用示例

### 附录（~50行）
- 常用参数速查表
- 参考文档链接
