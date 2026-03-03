# 四、PYATB 介电函数计算

本章使用 PYATB 读取 ABACUS 输出的紧束缚矩阵，在密 k 网格上计算介电函数张量。

## 4.1 准备矩阵文件

将 ABACUS 输出的三个矩阵文件复制到 PYATB 工作目录：

```bash
cp OUT.silica/data* ./pyatb_OpticConductivity/
cd pyatb_OpticConductivity
```

目录结构：

```
pyatb_OpticConductivity/
├── Input                        # PYATB 控制文件（需手动填写）
├── data-HR-sparse_SPIN0.csr     # 从 ABACUS 复制
├── data-SR-sparse_SPIN0.csr     # 从 ABACUS 复制
└── data-rR-sparse.csr           # 从 ABACUS 复制
```

## 4.2 Input 文件详解

PYATB 的控制文件固定命名为 `Input`，由三个区块组成：

```
# Input
INPUT_PARAMETERS
{
nspin               1
package             ABACUS
fermi_energy        5.5385382545    # eV，从 SCF 日志中读取（E_Fermi 第二列）
fermi_energy_unit   eV
HR_route            data-HR-sparse_SPIN0.csr
SR_route            data-SR-sparse_SPIN0.csr
rR_route            data-rR-sparse.csr
HR_unit             Ry              # ABACUS 输出的哈密顿量单位为 Rydberg
rR_unit             Bohr            # 偶极矩阵单位为 Bohr
}

LATTICE
{
lattice_constant        1.8897261246257702
lattice_constant_unit   Bohr
lattice_vector
7.1199998856 0.0 0.0
0.0 7.1199998856 0.0
0.0 0.0 7.1199998856
}

OPTICAL_CONDUCTIVITY
{
occ_band    64      # 占据能带数，从 SCF 日志 "occupied bands" 读取
omega       0 30    # 计算能量范围：0 ~ 30 eV
domega      0.01    # 能量步长：0.01 eV，共 3000 个点
eta         0.1     # 展宽参数（eV），模拟有限寿命效应
grid        20 20 20  # 布里渊区积分网格
}
```

### 参数说明

**INPUT_PARAMETERS 区块**

| 参数 | 本案例值 | 说明 |
|------|----------|------|
| `nspin` | `1` | 与 ABACUS 的 nspin 保持一致 |
| `package` | `ABACUS` | 指定矩阵文件来源为 ABACUS |
| `fermi_energy` | `5.5385382545` | 费米能（eV），从 SCF 日志最后一行 E_Fermi 读取 |
| `HR_unit` | `Ry` | ABACUS 输出的 HR 单位固定为 Rydberg |
| `rR_unit` | `Bohr` | 偶极矩阵单位固定为 Bohr |

**LATTICE 区块**

晶格参数必须与 ABACUS STRU 文件完全一致：
- `lattice_constant`：取 STRU 中的值 `1.8897261246257702`（对应 1 Å 的 Bohr 换算因子）
- `lattice_vector`：与 STRU 中的 LATTICE_VECTORS 相同

**OPTICAL_CONDUCTIVITY 区块**

| 参数 | 本案例值 | 说明 |
|------|----------|------|
| `occ_band` | `64` | 占据能带数，需与 SCF 结果一致 |
| `omega` | `0 30` | 能量范围（eV）。上限覆盖感兴趣的光学跃迁 |
| `domega` | `0.01` | 能量步长（eV）。越小谱线越细，计算量正比增加 |
| `eta` | `0.1` | 展宽参数（eV），对应 Lorentz 展宽。越小峰形越尖锐，但需要更密 k 网格；太小会出现噪声 |
| `grid` | `20 20 20` | k 积分网格。光学性质通常需要比 SCF 更密的网格；对本案例 20³ 已足够，复杂材料可能需要 40³ 以上 |

> **关于 `eta` 的选择**：`eta` 同时影响谱线宽度和 k 网格需求。减小 `eta` 需同步加密 `grid`，否则积分不收敛，谱线出现虚假震荡。

## 4.3 运行与输出

```bash
export OMP_NUM_THREADS=1 && mpirun -np 16 pyatb
```

计算过程中不会在屏幕上输出进度信息，可以追踪日志文件：

```bash
tail -f Out/running.log
```

计算完成后，进入输出目录：

```bash
cd Out/Optical_Conductivity
ls
```

```
dielectric_function_imag_part.dat    # 介电函数虚部 ε₂
dielectric_function_real_part.dat    # 介电函数实部 ε₁
optical_conductivity.dat             # 光电导率
```

### 输出文件格式

以 `dielectric_function_imag_part.dat` 为例，文件格式如下：

```
# omega(eV)  xx    xy    xz    yx    yy    yz    zx    zy    zz
  0.00000   0.000000e+00  -6.324391e-15  ...
  0.01000   9.564624e-06  -6.324391e-15  ...
  0.02000   1.912933e-05  ...
```

- 第一列：光子能量 $\omega$（eV）
- 后 9 列：介电函数张量 $\epsilon_2^{\alpha\beta}$（$\alpha,\beta=x,y,z$）的各分量，顺序为 xx、xy、xz、yx、yy、yz、zx、zy、zz

对于立方 SiO₂，三个对角分量 xx = yy = zz，非对角分量接近 0（~10⁻¹⁴ 量级），体现材料的各向同性。

低能量端（< 光学带隙）的 $\epsilon_2$ 应接近 0；$\epsilon_1(0)$（零频极限）应接近材料的静态介电常数，本案例约为 1.90（参考下一章实部数据）。

> **验证**：SCF 收敛精度不足时，$\epsilon_2$ 在低频段可能出现非零"尾巴"，需检查 `scf_thr` 设置。
