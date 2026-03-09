# 二、从基态 SCF 到 TDDFT 演化

RT-TDDFT 计算的起点是基态 SCF。先跑 SCF，再在 SCF 基态上开始时间演化。本章以 Si 原胞为例，介绍完整的输入文件准备，并讨论规范选择的背景。

## 2.1 必要输入文件

ABACUS 的 RT-TDDFT 计算需要四类输入文件：

| 文件 | 作用 |
|------|------|
| `INPUT` | 控制计算类型、基组、时间步、电场参数等 |
| `STRU` | 晶体结构（晶格矢量 + 原子位置） |
| `KPT` | 布里渊区 k 点采样 |
| 赝势 (`.upf`) 和轨道 (`.orb`) 文件 | 元素的赝势和数值原子轨道 |

### Si 原胞结构（STRU）

Si 金刚石结构原胞含 2 个原子，采用面心立方原胞表示：

```
ATOMIC_SPECIES
Si   28.085   Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb

LATTICE_CONSTANT
10.2

LATTICE_VECTORS
   0.5000000000     0.5000000000     0.0000000000
   0.0000000000     0.5000000000     0.5000000000
   0.5000000000     0.0000000000     0.5000000000

ATOMIC_POSITIONS
Direct
Si
0
2
   0.0000000000     0.0000000000     0.0000000000  m  1  1  1
   0.2500000000     0.2500000000     0.2500000000  m  1  1  1
```

晶格常数 10.2 Bohr ≈ 5.40 Å，接近 Si 的实验值（5.43 Å）。LCAO 基组选用 `Si_gga_8au_60Ry_2s2p1d.orb`（截断半径 8 Bohr，60 Ry 截断能，2s2p1d 构型）。

### k 点设置（KPT）

```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

RT-TDDFT 对 k 点采样的要求与静态计算类似，周期体系中需合理采样布里渊区。4×4×4 的 Gamma 中心网格适用于 Si 原胞的基础计算。

### 基态 SCF 的 INPUT

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix           Si_scf
calculation      scf
pseudo_dir       ./
orbital_dir      ./

#Parameters (2.Accuracy)
ecutwfc          60
scf_thr          1e-8
scf_nmax         100

#Parameters (3.Basis)
basis_type       lcao
gamma_only       0

#Parameters (4.Smearing)
smearing_method  gauss
smearing_sigma   0.01
```

> **注意：** `gamma_only` 必须设为 `0`。RT-TDDFT 使用传播子迭代计算复数波函数，不兼容 `gamma_only = 1`（后者只存储实数）。这个设置在 SCF 阶段就必须确定，以保证 TDDFT 阶段能读取正确格式的电荷密度。

## 2.2 为什么周期体系必须用速度规范——以及它的局限性

在 RT-TDDFT 中，加入外加电场的方式取决于**规范选择**。从电动力学的基本关系出发：

$$\mathbf{E} = -\nabla V - \frac{\partial \mathbf{A}}{\partial t}$$

有两种常用选择：

**长度规范**：将电场归为纯标势，加入 $V \propto \mathbf{E} \cdot \mathbf{r}$。哈密顿量形式简单，计算快，但 $\mathbf{r}$ 在周期性体系边界处不连续——只能用锯齿状势来衔接，对周期性体系不适用。

**速度规范**：将电场归为矢势 $\mathbf{A}(t) \propto \int_0^t \mathbf{E}(t') dt'$，保持哈密顿量的周期性。这是处理周期性体系的标准方法。

**速度规范在 NAO 下的问题**

理论上，长度规范和速度规范通过如下规范变换完全等价：

$$\psi^{\text{velocity}}(\mathbf{r}, t) = e^{i\mathbf{A}(t)\cdot\mathbf{r}} \psi^{\text{length}}(\mathbf{r}, t)$$

这个相位因子 $e^{i\mathbf{A}(t)\cdot\mathbf{r}}$ 描述了矢势引起的空间相位变化。在完备基组下，两种规范的结果完全一致。

但在 NAO 基组下，波函数被展开为原子轨道的线性组合：

$$\psi(\mathbf{r}, t) = \sum_{\mu} c_\mu(t) \phi_\mu(\mathbf{r})$$

这里的原子轨道 $\phi_\mu(\mathbf{r})$ 在规范变换时不会自动获得 $e^{i\mathbf{A}(t)\cdot\mathbf{r}}$ 相位。结果是：直接使用速度规范的 NAO 展开，**忽略了矢势在每个原子轨道内部引起的相位变化**。

这一遗漏导致计算的电流和能量出现系统性误差。在响应电流的计算中，低频区域会出现明显的虚假发散；在介电函数等关键物理量上，精度随基组减小而快速降低。

增大基组可以缓解误差，但计算成本急剧上升——这正是 NAO 基组的速度优势被抵消的地方。**混合规范直接从根源上修正这一误差。**
