# 第三章：SCF自洽计算

## 3.1 SCF的作用

自洽计算（SCF）是能带计算的第一步，其核心目标是获得收敛的电荷密度。

### 为什么需要SCF

在DFT中，电子的哈密顿量依赖于电荷密度，而电荷密度又由电子波函数决定。这形成了一个自洽问题：
1. 给定初始电荷密度
2. 构建哈密顿量
3. 求解本征值和本征态
4. 计算新的电荷密度
5. 重复2-4直到收敛

收敛的电荷密度是体系的基态电荷分布，是后续NSCF计算的基础。

### SCF与NSCF的区别

| 项目 | SCF | NSCF |
|------|-----|------|
| 目的 | 获得收敛电荷密度 | 计算能带 |
| 电荷密度 | 迭代更新 | 固定不变 |
| k点采样 | 均匀MP网格 | 高对称路径 |
| 计算量 | 多次迭代 | 单次对角化 |
| 输出 | SPIN1_CHG.cube | BANDS_1.dat |


## 3.2 INPUT_scf文件

以下是Si的SCF计算INPUT文件：

```
INPUT_PARAMETERS

# 系统设置
suffix              Si
calculation         scf
pseudo_dir          ./
ntype               1
nbands              26

# 平面波参数
basis_type          pw
ecutwfc             40.0
pw_diag_nmax        20
pw_diag_ndim        2

# SCF收敛参数
scf_thr             1e-08
scf_nmax            100
mixing_type         broyden
mixing_beta         0.8

# K点求解器
ks_solver           dav_subspace
smearing_method     gauss
smearing_sigma      0.015

# 输出设置
out_chg             1
```


**关键参数说明**：

**系统设置**：
- `suffix = Si`：输出目录名为 OUT.Si
- `calculation = scf`：自洽计算
- `pseudo_dir = ./`：赝势文件所在目录
- `ntype = 1`：元素种类数（只有Si）
- `nbands = 26`：能带数量（Si有4个价电子，8个原子共32个价电子，占据16条能带，26条足够）

**平面波参数**：
- `basis_type = pw`：使用平面波基组
- `ecutwfc = 40.0`：截断能（单位：Ry）
- `pw_diag_nmax = 20`：对角化最大迭代次数
- `pw_diag_ndim = 2`：Davidson算法的子空间维度

**SCF收敛参数**：
- `scf_thr = 1e-08`：电荷密度收敛阈值（单位：无量纲）
- `scf_nmax = 100`：最大SCF迭代次数
- `mixing_type = broyden`：电荷密度混合方法（Broyden方法，收敛快）
- `mixing_beta = 0.8`：混合系数（0-1之间，越大收敛越快但可能不稳定）

**K点求解器**：
- `ks_solver = dav_subspace`：Davidson子空间对角化方法
- `smearing_method = gauss`：高斯展宽（处理费米面附近的占据）
- `smearing_sigma = 0.015`：展宽宽度（单位：Ry）

**输出设置**：
- `out_chg = 1`：输出电荷密度文件（二进制+cube格式）


## 3.3 KPT_scf文件

SCF计算需要密集的k点网格以准确积分布里渊区：

```
K_POINTS
0
Gamma
12 12 12 0 0 0
```

**参数说明**：
- 第1行：关键词
- 第2行：0表示自动生成k点
- 第3行：Gamma模式（k点网格包含Γ点）
- 第4行：12×12×12网格，无平移

**k点网格选择**：
- 12×12×12对Si是合理的选择（共约864个约化k点）
- 更密集的网格提高精度但增加计算量
- k点收敛性测试在第1篇教程中已讨论

## 3.4 运行SCF计算

**运行命令**：
```bash
cp INPUT_scf INPUT
cp KPT_scf KPT
mpirun -np 8 abacus | tee scf.log
```

**检查收敛性**：
```bash
grep "ETOT DIFF" OUT.Si/running_scf.log | tail -10
```

收敛的标志是ETOT DIFF（总能量差）和DRHO（电荷密度差）都小于阈值。

**输出文件**：
- `OUT.Si/SPIN1_CHG.cube`：电荷密度文件（NSCF需要）
- `OUT.Si/running_scf.log`：运行日志
- `OUT.Si/SPIN1_CHG`：二进制电荷密度文件

**提取费米能级**：
```bash
grep -i "EFERMI" OUT.Si/running_scf.log
```

输出示例：`EFERMI = 6.2345 eV`

记录这个值，绘制能带图时需要。

