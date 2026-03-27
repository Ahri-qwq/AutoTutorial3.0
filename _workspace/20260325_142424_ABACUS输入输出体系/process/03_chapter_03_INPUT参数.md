## 第三章 INPUT：计算控制参数

`INPUT` 文件控制计算的全部行为：计算类型、精度、算法选择。文件必须命名为 `INPUT`（不可更改）。

### 3.1 文件语法

```
INPUT_PARAMETERS        # 第一行必须是这个关键字，之前的内容被忽略

# 以 # 或 / 开头的行是注释，会被忽略

suffix   Si             # 参数名与值之间用空格分隔
ecutwfc  60             # 参数名不区分大小写
```

- 第一行必须是 `INPUT_PARAMETERS`，否则文件不被识别
- 参数不区分大小写（`ecutwfc` 和 `ECUTWFC` 等价）
- 布尔值支持多种写法：`1`/`0`、`true`/`false`、`T`/`F`
- 运行结束后，`OUT.suffix/INPUT` 文件会列出所有实际使用的参数值（含默认值），可用于核查

**设置原则**：参数不是越多越好，未设置的参数使用默认值。默认值通常已能满足常见计算，只在有明确需要时才修改。

### 3.2 参数分组与常用参数

INPUT 参数可分为五组：

#### 组一：基本参数（General）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| suffix | ABACUS | 输出目录后缀，结果保存在 OUT.suffix/ |
| calculation | scf | 计算类型，见下表 |
| basis_type | pw | 基组类型：pw 或 lcao |
| pseudo_dir | ./ | 赝势文件目录 |
| orbital_dir | ./ | 轨道文件目录（lcao 专用） |
| symmetry | 1 | 是否利用晶体对称性（推荐开启） |
| ntype | — | 元素种类数（通常可省略，ABACUS 自动从 STRU 读取） |

**calculation 可选值**：

| 值 | 含义 |
|----|------|
| scf | 自洽电子结构计算（Self-Consistent Field） |
| relax | 离子弛豫（固定晶胞，优化原子位置） |
| cell-relax | 变胞弛豫（同时优化晶胞参数和原子位置） |
| md | 分子动力学 |
| nscf | 非自洽计算（读入已有电荷密度，计算能带/DOS） |

#### 组二：SCF 迭代

| 参数 | 默认值 | 说明 |
|------|--------|------|
| ecutwfc | — | 平面波截断能（Ry），PW 计算必填，见 3.5 节 |
| scf_nmax | 100 | SCF 最大迭代步数 |
| scf_thr | 1e-6 (pw) / 1e-7 (lcao) | 电荷密度收敛阈值（Ry），两步之间密度差 Δρ < scf_thr 时判断收敛 |

**scf_thr 建议值**：
- 一般计算：`1e-6`（PW）或 `1e-7`（LCAO）
- 需要高精度（如声子、弹性常数）：`1e-8`

#### 组三：KS 方程求解

| 参数 | 默认值 | 说明 |
|------|--------|------|
| nbands | 自动 | 计算的能带数，默认根据价电子数自动设置，需要空带（能带、DOS）时手动调大 |
| ks_solver | cg (pw) / genelpa (lcao) | KS 方程的对角化算法 |

#### 组四：展宽（Smearing）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| smearing_method | gauss | 展宽方法，见下方说明 |
| smearing_sigma | 0.015 Ry | 展宽宽度（Ry），约 0.2 eV |

**smearing_method 选择**：

| 方法 | 适用体系 |
|------|---------|
| gauss / gaussian | 通用，适合绝缘体和金属 |
| mp | 金属体系推荐，精度更高 |
| fixed | 绝缘体（不展宽，需要带隙严格大于 smearing_sigma） |
| fd | Fermi-Dirac，用于高温计算（SDFT） |

**smearing_sigma 建议**：
- 绝缘体：0.01 Ry 以下
- 金属：0.02 Ry，若收敛困难可适当增大到 0.05 Ry

#### 组五：电荷密度混合（Mixing）

混合参数控制 SCF 如何从旧电荷密度更新到新电荷密度，是影响 SCF 收敛速度和稳定性的关键。详见 3.4 节。

---

### 3.3 一个完整的 Si SCF 示例

以 8 原子金刚石结构 Si 为例，下面是参数设置和注释：

```
INPUT_PARAMETERS
# === 基本参数 ===
suffix                  Si          # 输出到 OUT.Si/
calculation             scf
symmetry                1
pseudo_dir              ./          # 赝势目录
basis_type              pw

# === SCF 迭代 ===
ecutwfc                 60          # 截断能，60 Ry 对 Si 已收敛
scf_nmax                100
scf_thr                 1e-8        # 严格收敛阈值

# === KS 方程求解 ===
nbands                  26          # Si 8 原子：4 价电子/原子×8=32 个电子，需 16 条占据带+空带

# === 展宽 ===
smearing_method         gauss
smearing_sigma          0.01        # Si 是半导体，展宽可以小些

# === 混合 ===
mixing_type             broyden
mixing_beta             0.7         # Si 带隙 ~1 eV，用 0.7
mixing_gg0              0           # 绝缘体关闭 Kerker 预处理
```

---

### 3.4 电荷密度混合参数

SCF 迭代通过"混合"旧电荷密度 ρ_in 和本步算出的 ρ_out，生成下一步的输入电荷密度。混合策略直接影响收敛速度和稳定性。

#### mixing_type

选择混合算法。

| 值 | 算法 | 说明 |
|----|------|------|
| broyden | Broyden 准牛顿法 | 默认，收敛最快 |
| pulay | Pulay DIIS | 略慢于 broyden |
| plain | 线性混合 | 最慢，仅用于调试 |

无特殊理由保持默认 `broyden`。

#### mixing_beta

新电荷密度占混合结果的比例，取值范围 0~1。

```
ρ_new = mixing_beta × ρ_out + (1 - mixing_beta) × ρ_in
```

- 默认值：0.8（nspin=1），0.4（nspin=2/4）
- **越大，收敛越快，但不收敛的风险越高**
- **越小，收敛越慢，但更稳定**

经验选取：

| 体系 | 推荐值 |
|------|--------|
| 绝缘体 / 半导体（带隙 > 1 eV） | 0.7 |
| 金属 / 过渡金属 | 0.2~0.4 |
| 收敛困难的体系 | 逐步减小，试 0.3、0.2、0.1 |

#### mixing_ndim

保存历史电荷密度的步数，默认 8。Broyden 和 Pulay 算法利用历史步信息构建更好的更新方向。增大 mixing_ndim 可提升收敛稳定性，代价是内存消耗线性增加。

#### mixing_gg0（Kerker 预处理）

默认 1.0（开启）。Kerker 方法在倒空间中衰减长波长（低频）的电荷密度分量，对金属体系收敛帮助很大。

```
mixing_gg0   1.0    # 金属体系（默认，推荐保持）
mixing_gg0   0.0    # 绝缘体/原子/分子体系（可关闭）
```

#### 收敛困难时的调参策略

| 情形 | 建议操作 |
|------|---------|
| drho 振荡不降 | 减小 mixing_beta（试 0.3、0.2、0.1） |
| 收敛太慢（>100 步）| 增大 mixing_ndim（试 12、20） |
| 金属收敛困难 | 确认 mixing_gg0=1.0，减小 smearing_sigma |
| 绝缘体/原子体系 | 关闭 Kerker：mixing_gg0=0.0 |
| DFT+U 体系 | 尝试开启 mixing_dmr，配合 mixing_restart 1e-3 |

---

### 3.5 平面波截断能（ecutwfc）与收敛测试

`ecutwfc` 控制展开波函数的平面波个数上限，单位 Ry。动能等于 ecutwfc 的平面波是"最高频"的基函数——ecut 越大，基组越完整，计算结果越精确，但计算量与 ecutwfc^(3/2) 成正比。

**只有 PW 计算才需要测试 ecutwfc**。LCAO 计算中 ecutwfc 仅控制辅助格点密度，不决定基组大小。

#### Al FCC 截断能收敛数据

以 Al FCC（4 原子原胞，k 点 6×6×6）为例，逐步增大 ecutwfc，记录总能量：

| ecutwfc (Ry) | 总能量 (eV) |
|-------------|------------|
| 20 | -57.197 |
| 25 | -56.2212 |
| 30 | -55.016 |
| 40 | -54.0526 |
| 45 | -53.9783 |
| 50 | -53.9617 |
| 55 | -53.9591 |
| 70 | -53.9602 |
| 80 | -53.9606 |

> 数据来源：Al 元素晶体的自洽迭代计算与平面波收敛测试及 k 点收敛性测试

从 20 Ry 到 40 Ry，能量迅速下降（~3 eV），说明基组严重不完整。40 Ry 之后趋于平坦，50 Ry 与 55 Ry 的差值仅 0.0026 eV（2.6 meV），与 70 Ry 的差值为 0.001 eV（1 meV），满足收敛标准。**Al 使用 50 Ry 即可，精度要求高时取 60 Ry**。

**收敛判断标准**：相邻 ecutwfc 的能量差 < 1 meV/atom。

#### Si 的截断能收敛参考

Si 8 原子体系在 ecut=60 Ry 时收敛（与 50 Ry 差 < 1 meV/atom）。

#### 各元素 ecutwfc 推荐范围

ecutwfc 的合理范围主要由赝势的"硬度"决定：

| 元素类型 | 典型 ecutwfc (Ry) |
|----------|-----------------|
| 简单金属（Al、Si、Na）| 40~60 |
| 过渡金属（Fe、Ni）| 60~100 |
| 含 O、N 的氧化物 | 60~80 |
| 第一行元素（C、N、O）| 60~80 |

具体值应以收敛测试为准，上表仅供参考。

#### 批量提交收敛测试的脚本思路

```bash
# 依次运行不同 ecutwfc 的 SCF 计算
for ecut in 20 30 40 50 60 70 80; do
    mkdir -p ecut_${ecut}
    cp STRU KPT ecut_${ecut}/
    sed "s/ecutwfc.*/ecutwfc  ${ecut}/" INPUT > ecut_${ecut}/INPUT
    cd ecut_${ecut}
    mpirun -n 4 abacus > log 2>&1
    cd ..
done

# 提取各 ecutwfc 下的总能量
for ecut in 20 30 40 50 60 70 80; do
    echo -n "ecutwfc=${ecut}: "
    grep FINAL_ETOT_IS ecut_${ecut}/OUT.Al/running_scf.log
done
```

---

### 3.6 Al FCC 完整 INPUT 示例

综合以上讨论，Al FCC SCF 计算的 INPUT 文件：

```
INPUT_PARAMETERS
# === 基本参数 ===
suffix                  Al
calculation             scf
symmetry                1
pseudo_dir              ./

# === 基组 ===
basis_type              pw
ecutwfc                 50          # Al 收敛值

# === SCF 迭代 ===
scf_nmax                100
scf_thr                 1e-8

# === 展宽（Al 是金属，用 mp 展宽）===
smearing_method         mp
smearing_sigma          0.02

# === 混合（金属，mixing_beta 取较小值）===
mixing_type             broyden
mixing_beta             0.4
mixing_gg0              1.0         # 金属开启 Kerker 预处理
```
