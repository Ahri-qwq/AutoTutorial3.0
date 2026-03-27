# 调研结果

## 1.1 RAG 检索结果总结

### ABACUS 软件背景
- **全称**：Atomic-orbital Based Ab-initio Computation at UStc（原子算筹）
- **性质**：国产开源密度泛函理论软件
- **功能**：电子结构自洽迭代、结构优化、分子动力学，支持平面波（PW）和数值原子轨道（LCAO）两种基组
- **GitHub**：https://github.com/deepmodeling/abacus-develop（DeepModeling 社区）
- **文档**：https://abacus.deepmodeling.com/en/latest/
- **赝势**：模守恒赝势 + 周期性边界条件
- **来源**：中国科学技术大学（UStc）开发，现由 DeepModeling 社区维护

### 计算文件体系
- 三个基本输入文件：INPUT（参数）、KPT（k点）、STRU（结构）
- 赝势文件（.upf）和轨道文件（.orb，LCAO时需要）
- 输出目录：OUT.suffix/

### INPUT 文件格式
- 必须以 `INPUT_PARAMETERS` 开头
- `#` 或 `/` 开头的行为注释
- 布尔值可用 True/False、1/0、T/F
- 文件名不可更改

**通用参数：**
- `suffix`：输出目录后缀，默认 ABACUS → OUT.ABACUS/
- `ntype`：元素种类数
- `calculation`：计算类型（scf/relax/cell-relax/md）
- `basis_type`：基组（pw/lcao）
- `ecutwfc`：平面波截断能（Ry）
- `pseudo_dir`：赝势目录
- `orbital_dir`：轨道目录（LCAO需要）
- `symmetry`：是否考虑对称性（0/1）

**SCF 迭代参数：**
- `scf_thr`：电荷密度收敛阈值（默认 1e-6 Ry）
- `scf_nmax`：SCF 最大步数（默认 100）
- `scf_ene_thr`：能量收敛阈值（默认 1，即不起效）

**Smearing 参数：**
- `smearing_method`：gauss（通用）、mp（金属推荐）、mp2、mv/cold、fixed（绝缘体）、fd
- `smearing_sigma`：展宽幅度（默认 0.015 Ry）

**Mixing 参数（影响SCF收敛）：**
- `mixing_type`：broyden（默认）/pulay/plain
- `mixing_beta`：混合系数（默认 0.8/0.4）
- `mixing_ndim`：历史保存数（默认 8）
- `mixing_gg0`：Kerker 预处理强度（默认 1.0）

**结构优化参数（relax/cell-relax）：**
- `cal_force`：是否计算力（0/1）
- `cal_stress`：是否计算应力（0/1）
- `force_thr_ev`：力收敛阈值（eV/Angstrom）
- `stress_thr`：应力收敛阈值（kBar）
- `relax_nmax`：最大离子步数

**输出参数：**
- `out_stru`：输出优化后结构（0/1）
- `out_chg`：输出电荷密度
- `out_band`：输出能带
- `out_dos`：输出态密度
- `out_mul`：输出 Mulliken 电荷

**k点相关：**
- `gamma_only`：只用 Gamma 点（大体系加速）
- `kspacing`：自动生成 KPT 的 k 点间距（可替代 KPT 文件）

### KPT 文件格式
```
K_POINTS
0           # 0 表示自动生成
Gamma       # Gamma 或 MP 方法
4 4 4 0 0 0 # 网格数 + 偏移
```
- `line mode`：用于能带计算的高对称路径模式

### STRU 文件格式
**必须包含的部分：**
1. `ATOMIC_SPECIES`：元素名 原子质量 赝势文件
2. `NUMERICAL_ORBITAL`（LCAO时）：轨道文件
3. `LATTICE_CONSTANT`：晶格常量（Bohr，1 Å = 1.8897259886 Bohr）
4. `LATTICE_VECTORS`：晶格向量（3x3矩阵，与 LATTICE_CONSTANT 相乘得实际长度）
5. `ATOMIC_POSITIONS`：坐标系类型（Direct 或 Cartesian），然后每种元素的原子坐标

**原子位置块格式：**
```
元素名
初始磁矩
原子数
x y z  移动标志(1/0)  [mag 磁矩]
```

### 输出文件体系
**OUT.suffix/ 下的主要文件：**
- `INPUT`：实际使用的所有参数（含默认值）
- `running_scf.log`：主日志文件，记录所有计算过程
- `KPT.info`：k 点信息
- `CHARGE-DENSITY.restart`：电荷密度重启文件
- `STRU.cif`：输出的结构文件（CIF 格式）
- `istate.info`：能级信息

**running_scf.log 结构：**
1. ABACUS 版本信息和时间戳
2. 输入参数和警告信息
3. 计算设置（网格、k点、基组初始化）
4. SCF 迭代过程：每步显示 ETOT、EDIFF、DRHO、时间
5. 最终能量：`!FINAL_ETOT_IS -XXXX.XXX eV`
6. 各模块运行时间分析

**收敛判断：** DRHO < scf_thr 即视为收敛

### abacustest 功能
- **安装**：`pip install abacustest`
- **GitHub**：https://github.com/pxlxingliang/abacus-test
- **准备输入文件**：`abacustest model inputs prepare`
  - `-f`：结构文件
  - `--pp`：赝势
  - `--orb`：轨道文件（LCAO）
  - `--input`：INPUT 文件
- **后处理结果**：`abacustest model inputs post`
- **DOS/PDOS 作图**：`abacustest model dos-pdos`

---

## 1.2 风格参考总结（基于《磁性材料计算》等3篇参考文章）

**开头方式：**
- 直接从物理概念切入，不铺垫"随着XX发展..."
- 磁性文章开头："材料的磁性对应的微观电子结构为电子自旋的方向并不均匀分布。"
- 简洁，用一句话定义核心概念

**结构偏好：**
- 使用 `##` 章节，`###` 小节，`####` 子小节
- 参数说明常用表格（避免表格内加粗）
- 代码块直接给文件内容

**语言特征：**
- 无"综上所述"、"值得注意的是"等 AI 腔
- 中英文混排，参数名保留英文
- 参数解释简洁，格式：`参数名`：含义（单位，默认值）

**代码呈现：**
- 代码块加文件名注释
- 参数按逻辑分组（General/Smearing/Mixing 等）

---

## 1.3 风险点

1. **abacustest 的 post 命令**：检索结果中关于如何抓取关键结果（如 etot、band_gap）的命令示例不充分，需要在文中介绍时以实际命令为准，不可猜测参数
2. **STRU 文件 Cartesian vs Direct**：两种坐标系格式均有说明，文中需清楚区分
3. **后续内容边界**：收敛测试（ecutwfc/k点收敛）属于后续篇章范畴，本篇只介绍文件格式，不展开收敛测试流程
