# 第四章：NSCF能带计算

## 4.1 高对称路径的选择

### 布里渊区与高对称点

能带结构需要沿布里渊区的高对称路径计算。高对称点是布里渊区中具有特殊对称性的k点，能带在这些点附近通常有极值或特殊特征。

### Si金刚石结构的布里渊区

Si采用面心立方（FCC）晶格，其倒空间布里渊区是截角八面体。

**主要高对称点**：
- **Γ (Gamma)**：布里渊区中心，坐标 (0, 0, 0)
- **X**：立方边中点，坐标 (0.5, 0, 0.5)
- **L**：立方体顶点方向，坐标 (0.5, 0.5, 0.5)
- **K**：六边形面中心，坐标 (0.375, 0.375, 0.75)
- **U**：六边形边中点，坐标 (0.625, 0.25, 0.625)
- **W**：立方边与六边形面交点，坐标 (0.5, 0.25, 0.75)

### Si的标准高对称路径

根据文献（Setyawan & Curtarolo, Comp. Mater. Sci. 2010），Si的标准能带路径为：

**L → Γ → X → U,K → Γ**

这条路径覆盖了Si能带的关键特征：
- **Γ点**：价带顶
- **X点附近**：导带底（约0.85倍X点位置）
- **L点**：次高价带


## 4.2 INPUT_nscf文件

以下是NSCF计算的INPUT文件：

```
INPUT_PARAMETERS

# 系统设置
suffix              Si
calculation         nscf
pseudo_dir          ./
ntype               1
nbands              26

# 平面波参数
basis_type          pw
ecutwfc             40.0
pw_diag_nmax        20
pw_diag_ndim        2

# NSCF特殊设置
symmetry            0
init_chg            file

# SCF收敛参数（NSCF只迭代一次）
scf_thr             1e-08
scf_nmax            100
mixing_type         broyden
mixing_beta         0.8

# K点求解器
ks_solver           dav_subspace
smearing_method     gauss
smearing_sigma      0.015

# 输出设置
out_band            1
```


**与SCF的关键区别**：

**NSCF特殊设置**：
- `calculation = nscf`：非自洽计算
- `symmetry = 0`：**必须关闭对称性**（Line模式KPT要求）
- `init_chg = file`：读取SCF的电荷密度文件

**为什么symmetry必须设为0**：
- Line模式的k点路径是手动指定的，不遵循对称性
- 如果开启对称性，程序会尝试对k点进行对称操作，导致错误
- SCF可以使用`symmetry = 1`利用对称性加速，但NSCF不行

**输出设置**：
- `out_band = 1`：输出能带文件BANDS_1.dat
- 不需要`out_chg`（不更新电荷密度）

**电荷密度读取**：
- NSCF会自动在`OUT.Si/`目录查找`SPIN1_CHG.cube`
- 如果电荷密度文件在其他位置，需要设置`read_file_dir`参数

## 4.3 KPT_nscf文件

NSCF使用Line模式，沿高对称路径采样：

```
K_POINTS
6
Line
0.500  0.500  0.500  50  # L
0.000  0.000  0.000  50  # Γ
0.500  0.000  0.500  50  # X
0.625  0.250  0.625  20  # U
0.375  0.375  0.750  50  # K
0.000  0.000  0.000  1   # Γ
```

**参数说明**：
- 第1行：关键词
- 第2行：高对称点数量（6个）
- 第3行：Line模式
- 第4行起：每行为一个高对称点
  - 前3个数：k点的分数坐标
  - 第4个数：到下一个点的插值数
  - 注释：高对称点标记

**插值数选择**：
- L→Γ：50个点（较长路径）
- Γ→X：50个点
- X→U：20个点（较短路径）
- U→K：50个点
- K→Γ：1个点（终点，不插值）

插值数越多，能带曲线越光滑，但计算量越大。50个点是常用选择。

## 4.4 运行NSCF计算

**运行命令**：
```bash
cp INPUT_nscf INPUT
cp KPT_nscf KPT
mpirun -np 8 abacus | tee nscf.log
```

**注意事项**：
- 确保SCF已完成，`OUT.Si/SPIN1_CHG.cube`存在
- NSCF只需一次迭代，速度很快（通常几分钟）

**输出文件**：
- `OUT.Si/BANDS_1.dat`：能带数据文件
- `OUT.Si/running_nscf.log`：运行日志

## 4.5 BANDS_1.dat文件格式

能带数据保存在`BANDS_1.dat`中，格式如下：

```
1    0.0000    -5.234  -5.234  -3.123  -3.123  ...
2    0.0234    -5.230  -5.230  -3.120  -3.120  ...
3    0.0468    -5.218  -5.218  -3.110  -3.110  ...
...
```

**列说明**：
- 第1列：k点序号
- 第2列：k点距离（从起点累积，单位Å⁻¹）
- 第3列起：各条能带的能量（单位eV）

**提取带隙信息**：
```bash
grep -i "bandgap" OUT.Si/running_nscf.log
```

输出示例：`E_bandgap 0.5734 eV`

这就是Si的PBE带隙，约0.57 eV，远小于实验值1.17 eV。

