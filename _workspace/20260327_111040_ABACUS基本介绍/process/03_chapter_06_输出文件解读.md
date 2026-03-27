# 第六章：输出文件解读

ABACUS 计算完成后，会在当前目录生成输出文件夹和日志文件。理解这些输出是判断计算是否成功、提取物理量的关键。

## 6.1 输出文件结构

计算完成后，目录结构如下：

```
计算目录/
├── INPUT
├── KPT
├── STRU
├── time.json              # 各模块计时信息（JSON 格式）
└── OUT.suffix/            # 主要输出文件夹（suffix 对应 INPUT 中设置）
    ├── INPUT              # 包含所有参数（含默认值）
    ├── running_scf.log    # 详细运行日志（与屏幕输出基本相同）
    ├── kpoints            # k 点信息
    ├── STRU_READIN_ADJUST.cif  # 读入的结构（CIF 格式）
    └── ...                # 其他输出文件
```

**time.json：** 记录各模块耗时，用于性能分析。

**OUT.suffix/INPUT：** 包含所有实际使用的参数，包括用户未设置的默认值，便于确认计算设置。

**OUT.suffix/running_scf.log：** 最重要的输出文件，包含完整的计算过程和结果。

## 6.2 屏幕输出解读

运行 ABACUS 时，屏幕会实时显示计算进度。以下逐部分解读（以 SCF 计算为例）：

### 6.2.1 软件信息

```
                              ABACUS v3.10.1

               Atomic-orbital Based Ab-initio Computation at UStc

               Website: http://abacus.ustc.edu.cn/
               Documentation: https://abacus.deepmodeling.com/
               Repository: https://github.com/deepmodeling/abacus-develop
               Commit: 1234abcd (2024-01-15)

 Fri Jan 15 10:30:00 2024
 MAKE THE DIR         : OUT.Si/
```

显示软件版本、官网、文档地址、代码提交版本和计算开始时间。

### 6.2.2 参数回显

```
 READING GENERAL INFORMATION
                           global_out_dir : OUT.Si/
                           global_in_card : INPUT
                               pseudo_dir : ./
                              pseudo_type : auto
                                    ntype : 1
...
```

回显读入的参数，用于确认设置是否正确。若有警告（如赝势价电子数异常），会在此处显示。

### 6.2.3 SCF 迭代过程

```
 DONE(0.123 SEC) : INIT SCF
 ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)
 GE1    -1.234567e+02  0.000000e+00   5.432e-02  1.234e+00
 GE2    -1.235012e+02  -4.450000e-02  2.123e-02  1.198e+00
 GE3    -1.235089e+02  -7.700000e-03  8.765e-03  1.201e+00
 GE4    -1.235102e+02  -1.300000e-03  3.210e-03  1.195e+00
 GE5    -1.235105e+02  -3.000000e-04  1.234e-03  1.203e+00
 GE6    -1.235106e+02  -1.000000e-04  4.567e-04  1.199e+00
 GE7    -1.235106e+02  -1.000000e-05  1.234e-04  1.197e+00
 GE8    -1.235106e+02  -1.000000e-06  3.456e-05  1.201e+00
 GE9    -1.235106e+02  -1.000000e-07  8.765e-06  1.195e+00
 GE10   -1.235106e+02  -1.000000e-08  2.345e-06  1.198e+00
 GE11   -1.235106e+02  -1.000000e-09  6.789e-07  1.200e+00
```

各列含义：

| 列 | 含义 |
|----|------|
| ITER | 迭代步数（GE = Geometry Electronic）|
| ETOT(eV) | 当前步的总能量（单位：eV）|
| EDIFF(eV) | 与上一步的能量差（单位：eV）|
| DRHO | 电荷密度变化（无量纲）|
| TIME(s) | 本步耗时（秒）|

**收敛判断：** 当 `DRHO < scf_thr` 时，SCF 收敛。本例中 `scf_thr = 1e-6 Ry ≈ 1.36e-5 eV`，第11步 DRHO 已小于阈值，收敛成功。

### 6.2.4 计算结果汇总

```
 --------------------------------------------
 !FINAL_ETOT_IS -1.235106e+02 eV
 --------------------------------------------

 Fermi Energy (eV) = 5.4321
```

**总能量：** `!FINAL_ETOT_IS` 后的数值是体系的总能量（单位：eV），这是 DFT 计算最核心的结果。

**费米能级：** 金属体系的费米能级，或半导体/绝缘体的价带顶能量（需结合能带判断）。

### 6.2.5 各模块耗时统计

```
 --------------------------------------------------------------------------------
 CLASS_NAME         NAME                TIME(Sec)  CALLS   AVG(Sec) PER(%)
 --------------------------------------------------------------------------------
 Driver             driver_line         27.239     1       27.24    100.00
 PW_Basis           setup_struc_factor  0.123      1       0.12     0.45
 ORB_control        read_orb_first      1.234      1       1.23     4.53
 ...
 --------------------------------------------------------------------------------
```

显示各模块的耗时、调用次数、平均耗时和占比，用于性能分析和优化。

## 6.3 running_scf.log 文件

`OUT.suffix/running_scf.log` 与屏幕输出内容基本相同，但保存在文件中便于后续查阅。包含：
- 完整的参数设置
- SCF 迭代详细过程
- 最终能量和费米能级
- 各模块耗时统计

## 6.4 提取关键结果

### 提取总能量

```bash
grep "!FINAL_ETOT_IS" OUT.Si/running_scf.log
```

输出：
```
 !FINAL_ETOT_IS -1.235106e+02 eV
```

### 提取费米能级

```bash
grep "Fermi Energy" OUT.Si/running_scf.log
```

### 提取力（relax 计算）

```bash
grep "TOTAL-FORCE" OUT.Si/running_relax.log
```

## 6.5 常见警告与错误

### SCF 不收敛

```
 SCF NOT CONVERGED after 100 iterations
```

**可能原因：**
- 初始结构不合理（原子距离过近）
- smearing 参数不当（金属体系 sigma 过小）
- mixing 参数需要调整（降低 mixing_beta）
- ecutwfc 或 k 点不足

**解决方法：**
- 检查 STRU 文件，确保结构合理
- 金属体系使用 `smearing_method mp`，增大 `smearing_sigma`
- 降低 `mixing_beta`（如从 0.4 降至 0.1）
- 增大 `scf_nmax` 至 200

### 内存不足

```
 Error: Memory allocation failed
```

**解决方法：**
- 减小 ecutwfc
- 减少 k 点数量
- 使用更多节点/核心分配内存
- LCAO 基组比 PW 节省内存

### k 点设置不当

```
 Warning: k-point mesh too coarse
```

**解决方法：**
- 增加 k 点密度（如从 2×2×2 增至 4×4×4）
- 做 k 点收敛性测试

> **提示：** 遇到问题时，先查看 `OUT.suffix/warning.log`（如果存在），其中记录了所有警告信息。
