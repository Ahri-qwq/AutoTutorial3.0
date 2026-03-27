# Step 2: 大纲方案（修订版）

## 用户反馈

用户要求：
- 这是培训宣讲用的教案文章
- 希望按照 **INPUT → KPT → STRU** 的顺序组织
- 理由：先讲控制参数，再讲采样，最后讲结构

## Think Aloud

重新设计思路：
1. 培训宣讲场景：需要清晰的讲解顺序，便于讲师授课
2. INPUT → KPT → STRU 顺序的优势：
   - INPUT 是计算的"大脑"，先讲让学员理解计算类型和参数
   - KPT 是采样策略，承上启下
   - STRU 是具体结构，最后讲可以结合前面的参数理解
3. 保持覆盖9个核心内容点
4. 篇幅控制在 800-1000 行

---

## 方案A：培训宣讲型（推荐）

**特点：** 按 INPUT → KPT → STRU 顺序，适合培训授课

**结构：**

### 前言（约40行）
- 培训目标：帮助学员快速掌握 ABACUS 基本使用
- 适用对象：第一性原理计算初学者、科研人员
- 前置知识：基本 Linux 命令、DFT 基础概念
- 教程结构：软件介绍 → INPUT 参数 → KPT 采样 → STRU 结构 → 实战案例
- 后续课程预告：力学性质、能带结构、态密度

### 第一章：ABACUS 软件简介（约80行）
- 1.1 软件背景
  - 国产开源 DFT 软件"原子算筹"
  - 开发团队：中国科学技术大学等
  - GitHub 与官方文档
- 1.2 核心功能
  - 支持的基组：平面波（PW）、数值原子轨道（LCAO）
  - 支持的计算：SCF、结构优化、分子动力学
  - 高级功能：DFT+U、杂化泛函、隐式溶剂等
- 1.3 ABACUS 计算文件体系
  - 必需的三个文件：INPUT、KPT、STRU
  - 辅助文件：赝势（.upf）、轨道（.orb）
  - 输出文件：OUT.suffix/ 文件夹
- 1.4 获取帮助
  - 官方文档：https://abacus.deepmodeling.com
  - GitHub 仓库：https://github.com/deepmodeling/abacus-develop
  - 社区支持

### 第二章：INPUT 文件详解（约200行）
- 2.1 INPUT 文件格式与规则
  - INPUT_PARAMETERS 开头
  - 参数格式：关键字 + 空格 + 值（不用等号）
  - 注释：# 或 / 开头
  - 布尔值：1/0、True/False、T/F（大小写不敏感）
  - 文件名固定为 INPUT，不可更改
- 2.2 通用参数（必需参数）
  - suffix：输出文件夹后缀（默认 ABACUS）
  - ntype：元素种类数
  - pseudo_dir：赝势文件目录
  - orbital_dir：轨道文件目录（LCAO 需要）
  - calculation：计算类型
    - scf：自洽电子结构计算
    - relax：固定晶胞的结构优化
    - cell-relax：晶胞和原子位置同时优化
    - md：分子动力学
  - basis_type：基组类型
    - pw：平面波基组
    - lcao：数值原子轨道基组
- 2.3 精度控制参数
  - ecutwfc：平面波截断能（单位：Rydberg）
    - 物理意义：控制平面波基组大小
    - 典型值：50-100 Ry
    - 即使用 LCAO 也需要设置（局部赝势用平面波）
  - scf_thr：SCF 收敛阈值（单位：Rydberg）
    - 电荷密度收敛判据
    - 典型值：1e-6 到 1e-8
  - scf_nmax：最大 SCF 迭代步数（默认 100）
- 2.4 常用功能参数
  - smearing_method：能级展宽方法
    - gauss：高斯展宽（绝缘体/半导体推荐）
    - mp：Methfessel-Paxton（金属推荐）
  - smearing_sigma：展宽宽度（单位：Ry，典型值 0.01-0.02）
  - symmetry：对称性开关（1 打开 / 0 关闭时间反演外 / -1 全关）
  - nspin：自旋极化（1 无磁 / 2 共线磁矩）
  - cal_force：计算原子受力（0/1）
  - cal_stress：计算晶胞应力（0/1）
- 2.5 结构优化专用参数
  - force_thr_ev：力收敛阈值（单位：eV/Å，推荐 0.01-0.05）
  - stress_thr：应力收敛阈值（单位：kBar，推荐 0.5-5）
  - relax_nmax：最大离子迭代步数（推荐 50-200）
  - out_stru：是否输出优化后结构（0/1）
- 2.6 完整 INPUT 示例
  - 示例1：Si 的 SCF 计算（PW 基组）
  - 示例2：MgO 的结构优化（LCAO 基组）

### 第三章：KPT 文件详解（约100行）
- 3.1 布里渊区采样基础
  - k 点的物理意义：描述晶体周期性
  - 采样密度与计算精度的关系
  - 绝缘体与金属的不同需求
- 3.2 KPT 文件格式
  - 第1行：K_POINTS（固定关键字）
  - 第2行：k 点总数（0 = 自动生成）
  - 第3行：采样方法
    - Gamma：以 Γ 点为中心的 Monkhorst-Pack 方法
    - MP：标准 Monkhorst-Pack 方法
  - 第4行：网格参数 `n1 n2 n3 s1 s2 s3`
    - n1/n2/n3：各方向细分数
    - s1/s2/s3：网格偏移（通常取 0 0 0）
- 3.3 k 点设置建议
  - 立方晶系：各向同性，如 4×4×4
  - 六方晶系：适当调整，如 6×6×4
  - 大体系（>100原子）：可用 gamma_only 1
- 3.4 完整 KPT 示例
  - 立方体系（Gamma 方法）
  - 六方体系
  - gamma_only 模式

### 第四章：STRU 文件详解（约150行）
- 4.1 STRU 文件格式概览
  - 各块的顺序与含义
  - PW 与 LCAO 写法的区别
- 4.2 ATOMIC_SPECIES 块
  - 格式：元素名 原子质量 赝势文件名
  - 赝势文件获取：http://abacus.ustc.edu.cn/pseudo/list.htm
  - 命名规范：元素_ONCV_PBE-1.0.upf
- 4.3 NUMERICAL_ORBITAL 块（LCAO 专用）
  - 每个元素一行轨道文件
  - 轨道文件命名规范：元素_gga_Xau_YRy_轨道.orb
  - 与赝势文件的匹配关系
- 4.4 LATTICE_CONSTANT 块
  - 晶格常数（缩放系数）
  - 单位：1.8897259886 Bohr = 1 Angstrom
  - 常见设置：直接写 1.8897259886（以 Angstrom 为单位）
- 4.5 LATTICE_VECTORS 块
  - 3 行 3 列的矩阵
  - 实际长度 = 矩阵数值 × LATTICE_CONSTANT
  - 支持非正交晶格
- 4.6 ATOMIC_POSITIONS 块
  - 坐标系：Direct（晶格坐标）或 Cartesian（笛卡尔）
  - 每个元素的写法：
    - 元素名
    - 初始磁矩（通常为 0.0）
    - 原子数
    - 每个原子坐标 + 移动标志（1=可移动，0=固定）
- 4.7 完整 STRU 示例
  - 示例1：Si 晶体（PW 基组，Direct 坐标）
  - 示例2：MgO 晶体（LCAO 基组，Direct 坐标）

### 第五章：结构优化实战（约120行）
- 5.1 结构优化的物理意义
  - 势能面与稳定构型
  - relax 与 cell-relax 的区别
- 5.2 完整案例：Si 晶体 cell-relax
  - INPUT 文件（含 relax 相关参数）
  - STRU 文件
  - KPT 文件
  - 运行命令
- 5.3 判断优化是否收敛
  - 屏幕输出中的收敛信息
  - LARGEST GRAD（力）
  - TOTAL-PRESSURE（应力）
- 5.4 优化结果
  - STRU_ION_D 文件：记录每步离子构型
  - 最终结构：OUT.suffix/STRU_READIN_ADJUST.stru

### 第六章：输出文件解读（约130行）
- 6.1 输出文件结构
  - time.json：计时文件
  - OUT.suffix/ 文件夹（suffix 对应 INPUT 中的 suffix 参数）
  - OUT.suffix/INPUT：实际使用的所有参数（含默认值）
  - OUT.suffix/running_scf.log：详细运行日志
- 6.2 屏幕输出解读
  - 启动信息：版本号、并行配置
  - 参数回显：确认关键参数
  - SCF 迭代过程
    - 迭代次数（ITER）
    - 总能量变化（ETOT）
    - 电荷密度变化（DRHO）
    - 收敛判断：DRHO < scf_thr
  - 计算结果汇总
    - 总能量（Total Energy）
    - 费米能级（Fermi Energy）
    - 磁矩（磁性体系）
  - 各模块耗时统计
- 6.3 running_scf.log 文件
  - 与屏幕输出基本相同
  - 但保存在文件中便于后续查阅
  - 包含更多细节
- 6.4 提取关键结果
  - 总能量提取（grep "!FINAL_ETOT_IS"）
  - 费米能级提取
  - 力的提取（grep "TOTAL-FORCE"）

### 第七章：使用 abacustest 准备输入与提取结果（约100行）
- 7.1 abacustest 简介
  - 功能定位：ABACUS 前后处理工具
  - 安装：pip install abacustest
  - GitHub：https://github.com/pxlxingliang/abacus-test
- 7.2 准备输入文件
  - abacustest model inputs 命令
  - 常用参数说明
  - 从 CIF 文件生成 STRU
  - 批量设置 INPUT 参数
- 7.3 提取计算结果
  - 支持的物理量
  - 输出格式（JSON 等）
- 7.4 实战：从结构文件到计算结果
  - 步骤1：下载 CIF 文件
  - 步骤2：abacustest 生成输入文件
  - 步骤3：运行 ABACUS
  - 步骤4：abacustest 提取结果

### 附录（约40行）
- A.1 参考资料
  - ABACUS 官方文档
  - 教程链接
- A.2 常见问题
  - 赝势和轨道文件如何选择
  - 计算参数推荐
- A.3 进阶方向
  - ecutwfc 和 k 点收敛性测试
  - 后续教程导航

**总计：约 960 行**

---

## 方案B：精简宣讲型

**特点：** 篇幅更紧凑，适合时间有限的培训场合

**结构（章节标题）：**

### 前言（约30行）

### 第一章：ABACUS 软件简介（约60行）
- 背景、功能、文件体系

### 第二章：INPUT 文件（约150行）
- 格式规则
- 通用参数（calculation、basis_type、ecutwfc、scf_thr）
- 功能参数（smearing、relax 参数）
- 示例：Si SCF

### 第三章：KPT 文件（约70行）
- 格式与参数
- 示例（Gamma/MP）

### 第四章：STRU 文件（约100行）
- 各块说明
- 示例：Si（PW）+ MgO（LCAO）

### 第五章：一个完整的计算：结构优化（约100行）
- INPUT + KPT + STRU 三文件汇总
- 运行与收敛判断
- 优化结果查看

### 第六章：输出文件解读（约80行）
- 屏幕输出关键信息
- running_scf.log
- 关键物理量提取

### 第七章：abacustest 工具（约80行）
- 安装与使用
- 准备输入文件
- 提取结果

### 附录（约30行）

**总计：约 700 行**

**优点：** 简洁，适合1-2小时培训
**缺点：** 部分内容不够深入

---

## 方案对比

| 维度 | 方案A（完整宣讲型） | 方案B（精简宣讲型） |
|------|---------------------|---------------------|
| 篇幅 | 约960行 | 约700行 |
| 深度 | 较深，参数说明详尽 | 适中，重点突出 |
| 适合场合 | 半天培训 / 自学教材 | 1-2小时培训 |
| INPUT 顺序 | ✅ | ✅ |
| KPT 顺序 | ✅ | ✅ |
| STRU 顺序 | ✅ | ✅ |

两个方案均按照 **INPUT → KPT → STRU** 顺序组织，请选择适合的篇幅深度：
- **方案A**：内容详尽，适合教材级别
- **方案B**：重点突出，适合宣讲培训
