# Step 1 研究报告

## 1.1 RAG检索结果汇总

### A. INPUT 文件

**文件结构（来源：ABACUS 的平面波计算与收敛性测试.md）**
- 第一行必须是 `INPUT_PARAMETERS`，之前的内容被忽略
- `#` 或 `/` 开头的行被忽略（注释）
- 参数不区分大小写
- 布尔值支持：True/False、1/0、T/F

**五大参数分组（来源：ABACUS 的平面波计算与收敛性测试.md + abacus_user_guide__abacus_pw.html.md）**

1. **基本参数**：suffix、calculation、symmetry、pseudo_dir、orbital_dir、basis_type
2. **SCF迭代**：ecutwfc、scf_nmax、scf_thr
3. **KS方程求解**：nbands、ks_solver
4. **展宽（Smearing）**：smearing_method、smearing_sigma
5. **电荷密度混合**：mixing_type、mixing_beta、mixing_gg0

**完整Si SCF示例（来源：abacus_user_guide__abacus_pw.html.md）**
```
INPUT_PARAMETERS
#Parameters (1.General)
suffix                  Si
calculation             scf
symmetry                1
pseudo_dir              .
orbital_dir             .
basis_type              pw
ecutwfc                 100
#Parameters (2. SCF iterations)
scf_nmax                100
scf_thr                 1e-8
#Parameters (3. Solve KS equation)
nbands                  26
ks_solver               cg
#Parameters (4.Smearing)
smearing_method         gauss
smearing_sigma          0.01
#Parameters (5.Mixing)
mixing_type             broyden
mixing_beta             0.7
mixing_gg0              0
```

**calculation 类型**
- `scf`：自洽电子结构计算
- `relax`：离子弛豫（固定晶胞）
- `cell-relax`：变胞弛豫
- `md`：分子动力学
- `nscf`：非自洽计算（需 `init_chg file`）

---

### B. SCF收敛算法参数（来源：ABACUS收敛性问题解决手册.docx + 电荷密度混合算法介绍.md）

**mixing_type**
- 可选：broyden、pulay、plain
- 默认：broyden
- 速度：Broyden > Pulay > plain
- 建议：无特殊理由保持默认

**mixing_beta**
- 作用：新旧电荷密度混合时新电荷的比例（0-1）
- 默认：0.8（nspin=1），0.4（nspin=2|4）
- 越大收敛越快，但不收敛风险也越大
- 经验值：能带 >1 eV 的绝缘体设 0.7；能带 <1 eV 的金属设 0.2

**mixing_ndim**
- 默认：8
- 保存历史电荷密度步数，越大收敛越稳定，但内存消耗线性增加

**mixing_gg0（Kerker预处理）**
- 默认：1.0（开启）
- 金属体系受益显著；绝缘体/原子体系可关闭（设0.0）

**smearing_method**
- gauss/gaussian：最通用
- mp：金属体系推荐
- fd：用于SDFT
- fixed：绝缘体

**smearing_sigma**
- 默认：0.015 Ry
- 越小计算越准确，也可能更难收敛；越大更易收敛但可能收敛到错误基态

---

### C. KPT 文件（来源：en__latest__advanced__input_files__kpt.html.md + 多个教程）

**MP网格格式：**
```
K_POINTS
0           # 0表示自动生成
Gamma       # Gamma 或 MP
4 4 4 0 0 0 # 沿三个倒格矢方向的细分数 + 偏移
```

**路径KPT格式（能带计算）：**
```
K_POINTS
8           # k点段数（高对称点连线数）
Line
0.00000000 0.00000000 0.00000000 25 # GAMMA
0.50000000 0.00000000 0.50000000 9  # X
...
```

**Al k点收敛测试数据（来源：Al元素晶体的自洽迭代计算与平面波收敛测试及k点收敛性测试.md）**
- k=1: -93.8589 eV
- k=2: -93.3732 eV
- k=3: -93.3681 eV
- k=4: -93.3671 eV
- k=5: -55.0161 eV（突变！—实为换了k点网格密度）
- k=8: -54.8113 eV
- k=12: -54.8056 eV
- k=16: -54.6115 eV
- 结论：k点从4到5能量有跳跃，从6以后收敛于约-54eV至-55eV附近

---

### D. STRU 文件（来源：多个教程）

**完整结构**：
```
ATOMIC_SPECIES
Al  26.982  Al_ONCV_PBE-1.0.upf

# LCAO 专有（PW不需要）：
NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb

LATTICE_CONSTANT
1.88972612546          # 单位Bohr（1 Bohr = 0.529177 Å）

LATTICE_VECTORS
4.03459549706 0 0
0 4.03459549706 0
0 0 4.03459549706

ATOMIC_POSITIONS
Direct                 # 或 Cartesian
Al                     # 元素名
0                      # 初始磁矩
4                      # 原子数
0.0  0.0  0.0  0 0 0  # x y z  move_x move_y move_z
0.5  0.5  0.0  0 0 0
0.5  0.0  0.5  0 0 0
0.0  0.5  0.5  0 0 0
```

**PW vs LCAO 区别：**
- PW：无需 NUMERICAL_ORBITAL 块；INPUT中 basis_type pw
- LCAO：需要 NUMERICAL_ORBITAL 块列出轨道文件；INPUT中 basis_type lcao

**Al FCC结构（来源：Al元素晶体...md + ABACUS+Phonopy.md）**
- 晶格常数：a = 4.03459549706 Å（已优化结构）或 a = 4.038930 Å（原始）
- 赝势：Al_ONCV_PBE-1.0.upf 或 Al_ONCV_PBE-1.2.upf
- 轨道（LCAO）：Al_gga_7au_100Ry_4s4p1d.orb

**坐标系选择**（来源：ABACUS 使用教程）
- `Direct`：分数坐标（推荐，与晶格无关）
- `Cartesian`：直角坐标，单位为 LATTICE_CONSTANT

**位置约束（0/1）**：move_x move_y move_z（0=固定，1=自由），结构优化时用

---

### E. 截断能收敛测试（ecutwfc）

**理论背景（来源：ABACUS 的平面波计算与收敛性测试.md）**
- ecutwfc 控制平面波基矢量个数，越大越精确但成本越高
- 仅对 PW 有效，LCAO不适用
- 收敛标准：相邻ecutwfc差 < 1 meV/atom

**Al ecutwfc 收敛数据（来源：Al元素晶体...md）**
- ecutwfc=20 Ry: -57.197 eV
- ecutwfc=25 Ry: -56.2212 eV
- ecutwfc=30 Ry: -55.016 eV
- ecutwfc=40 Ry: -54.0526 eV
- ecutwfc=45 Ry: -53.9783 eV
- ecutwfc=50 Ry: -53.9617 eV
- ecutwfc=55 Ry: -53.9591 eV
- ecutwfc=70 Ry: -53.9602 eV
- ecutwfc=80 Ry: -53.9606 eV
- 结论：ecutwfc=50 Ry时曲线已基本收敛，认为损失较小

**Si ecutwfc 收敛（来源：ABACUS 的平面波计算与收敛性测试.md）**
- 在 ecut=60 Ry 时认为收敛（与50 Ry差 < 1 meV/atom）

---

### F. 输出文件（来源：en__latest__quick_start__output.html.md + 多教程）

- `OUT.suffix/`：主输出目录（suffix = INPUT中的suffix参数值）
- `OUT.suffix/running_scf.log`：主日志，含收敛过程和最终能量
- `OUT.suffix/INPUT`：实际使用的所有参数（含默认值）
- `OUT.suffix/istate.info`：能级信息
- 提取总能量：`grep FINAL_ETOT_IS OUT.suffix/running_scf.log`

---

## 1.2 风格分析（参考文章：ABACUS 使用教程｜磁性材料计算.md）

- **开头方式**：直接说明材料性质的物理背景，无过多铺垫
- **结构**：介绍 → 准备（参数）→ 案例（文件+运行）→ 结果分析
- **表格使用**：参数对比常用表格（但表格内不用加粗）
- **代码块**：有文件名标签，参数对齐但不强求等宽
- **语言特征**：直接陈述，短句，不用"值得注意的是"等AI腔
- **参考文章行数**：362行（单篇）

## 1.3 检索质量评估

| 子主题 | 质量 | 主要来源 |
|--------|------|----------|
| INPUT参数详解 | ✅ 充足 | ABACUS的平面波计算与收敛性测试.md |
| SCF收敛算法参数 | ✅ 充足 | ABACUS收敛性问题解决手册.docx、电荷密度混合算法介绍.md |
| KPT采样策略 | ✅ 充足 | en__latest__advanced__input_files__kpt.html.md，多教程 |
| Al k点收敛数据 | ✅ 有具体数据 | Al元素晶体的自洽迭代计算与平面波收敛测试.md |
| STRU定义（PW+LCAO） | ✅ 充足 | 多教程（MgO、Al等示例） |
| ecutwfc收敛测试（Al） | ✅ 有具体数据 | Al元素晶体的自洽迭代计算...md |
| 输出文件解读 | ✅ 片段充足（附录级别） | 多教程 |
