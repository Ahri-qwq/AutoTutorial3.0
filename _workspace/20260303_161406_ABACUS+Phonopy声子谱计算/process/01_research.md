# 调研笔记

## 一、案例解析

### 1.1 案例基本信息
- **题目：** ABACUS+Phonopy 计算声子谱
- **体系：** FCC Al（面心立方铝）
- **方法：** 有限位移方法（Finite Displacement Method）
- **工具链：** ABACUS（受力计算）+ Phonopy（声子处理）+ gnuplot（绘图）

### 1.2 案例完整文件清单

| 文件 | 说明 |
|------|------|
| `STRU` | FCC Al 初始结构（已优化） |
| `INPUT` | ABACUS SCF 计算输入参数 |
| `psp/Al_ONCV_PBE-1.0.upf` | Al 赝势文件 |
| `psp/Al_gga_7au_100Ry_4s4p1d.orb` | Al 数值原子轨道文件 |
| `STRU-001` | Phonopy 生成的微扰构型（命令产生） |
| `band.conf` | Phonopy 声子谱计算配置 |
| `plot_pho.gp` | gnuplot 绘图脚本 |

### 1.3 完整计算流程

```
1. 准备 FCC Al 的优化结构（STRU）
2. Phonopy 生成超胞+微扰构型：phonopy -d --dim="2 2 2" --abacus → STRU-001
3. ABACUS 对 STRU-001 做 SCF 计算，输出原子受力（cal_force=1）
4. Phonopy 生成 FORCE_SETS：phonopy -f ./disp-001/OUT*/running_scf.log
5. 配置 band.conf，计算声子谱：phonopy -p band.conf --abacus → band.yaml
6. 导出 gnuplot 格式并绘图：phonopy-bandplot --gnuplot > pho.dat && gnuplot plot_pho.gp
```

### 1.4 关键文件内容

**STRU 文件（完整）：**
```
ATOMIC_SPECIES
Al 26.982 Al_ONCV_PBE-1.0.upf upf201

NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb

LATTICE_CONSTANT
1.88972612546

LATTICE_VECTORS
4.03459549706 0 0 #latvec1
0 4.03459549706 0 #latvec2
0 0 4.03459549706 #latvec3

ATOMIC_POSITIONS
Direct

Al #label
0 #magnetism
4 #number of atoms
0  0  0  m  0  0  0
0.5  0.5  0  m  0  0  0
0.5  0  0.5  m  0  0  0
0  0.5  0.5  m  0  0  0
```

**INPUT 文件（完整）：**
```
INPUT_PARAMETERS
#Parameters (1.General)
suffix          Al-fcc
calculation     scf
esolver_type    ksdft
symmetry        1
pseudo_dir      ./psp
orbital_dir     ./psp
cal_stress      1
cal_force       1
stru_file       STRU-001

#Parameters (2.Iteration)
ecutwfc         100
scf_thr         1e-7
scf_nmax        50

#Parameters (3.Basis)
basis_type      lcao
gamma_only      0

#Parameters (4.Smearing)
smearing_method mp
smearing_sigma  0.015

#Parameters (5.Mixing)
mixing_type     pulay
mixing_beta     0.7
mixing_gg0      1.5
```

**band.conf 文件（完整）：**
```
ATOM_NAME = Al
DIM = 2 2 2
MESH = 8 8 8
PRIMITIVE_AXES = 0 1/2 1/2  1/2 0 1/2  1/2 1/2 0
BAND = 1 1 1  1/2 1/2 1  3/8 3/8 3/4  0 0 0   1/2 1/2 1/2
BAND_POINTS = 101
BAND_CONNECTION = .TRUE.
```

**plot_pho.gp 文件（完整）：**
```gnuplot
set terminal pngcairo size 1920, 1080 font 'Arial, 36'
set output "Al-FCC_plot.png"
set ylabel 'Frequency (THz)'
set ytics 2
unset key

x1 = 0.13115990
x2 = 0.17753200
x3 = 0.31664810
xmax = 0.43023590
ymin = 0
ymax = 12

set xrange [0:xmax]
set yrange [ymin:ymax]
set xtics ("{/Symbol G}" 0, "X" x1, "K" x2, "{/Symbol G}" x3, "L" xmax)
set arrow 1 nohead from x1,ymin to x1,ymax lt 2
set arrow 2 nohead from x2,ymin to x2,ymax lt 2
set arrow 3 nohead from x3,ymin to x3,ymax lt 2

plot 'pho.dat' using 1:($2) w l lw 3
```

### 1.5 案例完整性检查清单
- [x] STRU 文件（含晶格常数、原子坐标、赝势/轨道文件名）
- [x] INPUT 文件（含所有必要参数）
- [x] 赝势文件名：Al_ONCV_PBE-1.0.upf
- [x] 轨道文件名：Al_gga_7au_100Ry_4s4p1d.orb
- [x] Phonopy 命令（扩胞、生成FORCE_SETS、计算声子谱）
- [x] band.conf（含所有参数）
- [x] 绘图脚本 plot_pho.gp（含所有高对称点坐标）
- [x] FCC Al 的 k 路径：Γ → X → K → Γ → L

### 1.6 特殊设置要点
- `stru_file STRU-001`：用于指定 Phonopy 生成的微扰结构
- `symmetry = 1`：必须开启，Phonopy 扩胞后对称性信息需要
- `gamma_only = 0`：超胞计算需要多 k 点
- `cal_force = 1`：必须开启，否则无受力输出
- FCC 的 `PRIMITIVE_AXES`：0 1/2 1/2  1/2 0 1/2  1/2 1/2 0（从常规 FCC 胞到原胞的变换）
- 扩胞经验：三个方向 cell 长度均在 10-20 Å 为宜

---

## 二、RAG 检索结果

### 2.1 主要检索结果
- 找到 ABACUS 官方 Phonopy 教程的两个版本（来源和 HTML 抓取版）
- 找到 ABACUS+ShengBTE 教程（包含 Phonopy 的类似步骤，Si 体系）
- 找到 ShengBTE ABACUS 文档中的 Phonopy+Si 案例（扩胞参数：2 2 2，FORCE_CONSTANTS 写入）
- 声子谱物理：RAG 检索到弹性常数文章中的声子概念（非直接，但有声子色散与弹性模量关系）

### 2.2 补充知识（ShengBTE 案例中的 Phonopy 细节）
- 当需要将 FORCE_CONSTANTS 输出（供 ShengBTE 使用）时，需在 band.conf 中添加：
  - `FORCE_CONSTANTS = WRITE`
  - `FULL_FORCE_CONSTANTS = .TRUE.`
- 本教程不涉及 ShengBTE，上述参数不需要

---

## 三、风格参考总结

参考文章：
- 磁性材料计算：直接说明，表格呈现参数，少废话
- DeePMD-kit 教程：一、二、三 编号，子节用 2.1/2.2，代码块紧跟说明

**风格特征提取：**
- 开头直接说明"本教程介绍..."，不用排比式引言
- 使用一、二、三 中文数字节序
- 代码块用反引号，文件名加粗或反引号
- 参数说明简洁，直接列出作用+推荐值
- 不用"值得注意的是"、"综上所述"等表达
