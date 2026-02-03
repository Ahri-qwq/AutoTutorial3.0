# PySCF 教程2.0

PySCF (Python-based Simulations of Chemistry Framework) 是一个开源的计算化学软件库，专为量子化学模拟设计。它提供了从 Hartree-Fock (HF)、密度泛函理论 (DFT) 到耦合簇 (CC) 等多种电子结构计算方法，且具有高度的模块化和可扩展性。本教程将介绍如何使用 PySCF 进行分子建模、单点能计算、性质分析（如偶极矩、Mulliken 电荷）以及核梯度计算，帮助读者快速掌握利用 Python 脚本驱动量子化学计算的核心流程。

## 基础知识学习

计算化学（Computational Chemistry）是理论化学的一个分支，它使用计算机来求解量子力学基本方程，从而预测和解释分子的性质、结构与反应性。通过计算化学方法，我们可以在原子和电子的层面上理解物质世界，而无需进行实际的化学实验。本教程将引导您学习如何计算分子的基本属性，例如电荷分布、几何结构和能量。

本教程适合需要从第一性原理出发，计算分子性质的化学、物理、材料科学等领域的研究人员和学生。如果您希望解决诸如计算分子偶极矩、分析原子电荷、寻找稳定构象或评估反应能垒等问题，那么本教程将为您提供必要的背景知识。

### 核心概念

为了从计算的角度理解分子，我们需要了解以下几个核心概念。这些概念构成了现代量子化学计算的理论基石。

#### 势能面与分子几何

在量子化学中，分子的能量并非固定不变，而是其原子核坐标的函数，这个多维的能量景观被称为**势能面（Potential Energy Surface, PES）**。分子最稳定的构型（即平衡几何结构）对应于势能面上的能量最低点（局部极小点）。在这些点上，作用在每个原子核上的净力为零。

作用在原子上的力可以通过计算能量相对于原子核坐标的**梯度（Gradient）**得到（力是梯度的负值）。当分子处于非平衡构象时，梯度会指示出能量下降最快的方向，这正是**几何优化**算法寻找稳定结构的基础。一个特殊的势能面驻点是**过渡态（Transition State）**，它在某个方向上是能量最高点（鞍点），代表了化学反应的瓶颈。

![Figure 1: 势能面与分子几何结构](https://github.com/BruceJackey/Picture/blob/main/Tutorial_step1_step2_bundle_pyscf_pyscf_fig1.png?raw=true)

*Figure 1. 一个简化的二维势能面示意图，展示了能量与两个原子核坐标之间的关系。分子的稳定构型对应于能量的局部极小点，而过渡态是连接两个极小点的鞍点。几何优化过程就是沿着能量下降最快的方向（与梯度相反）寻找能量最低点的过程。*


#### 电子结构方法：Hartree-Fock 与密度泛函理论

精确求解多电子体系的薛定谔方程极其困难，因此需要引入近似方法。这些方法的核心是描述体系中的电子行为，即所谓的电子结构。

- **Hartree-Fock (HF) 方法**：这是一种基础的从头计算（*ab initio*）方法。它将复杂的电子间相互作用简化为每个电子在其他所有电子平均场中的运动，虽然忽略了电子间的瞬时关联，但为更高级的方法提供了理论起点。

- **密度泛函理论（Density Functional Theory, DFT）**：这是一种更现代且广泛使用的方法。它巧妙地将体系的能量表示为电子密度的泛函。相比HF方法，DFT在考虑计算成本和精度的平衡上通常表现更优，是目前计算化学领域最流行的工具之一。

![Figure 2: HF 与 DFT 方法核心思想对比](https://github.com/BruceJackey/Picture/blob/main/Tutorial_step1_step2_bundle_pyscf_pyscf_fig2.png?raw=true)

*Figure 2. 该图展示了 Hartree-Fock (HF) 方法与密度泛函理论 (DFT) 在处理电子相互作用上的核心区别。HF 将电子间作用简化为单个电子在其他所有电子平均场中的运动，而 DFT 则将体系能量直接与总电子密度关联起来。*


#### 分子性质的计算

一旦通过HF或DFT等方法求解了分子的电子结构，我们就可以计算出各种可观测量和化学概念。

- **电荷分布**：分子的电子云分布决定了其电学性质。这可以通过两种方式来描述：
  1.  计算整个分子的**电偶极矩**，它衡量了分子正负电荷中心的分离程度，是分子极性的宏观体现。
  2.  通过**布居分析**（如 Mulliken 分析）将总电子数划分到每个原子上，得到**部分原子电荷**，这有助于直观理解分子内局部的电荷富集与亏缺情况。

- **能量与反应性**：计算化学方法可以精确计算分子的总能量。通过比较不同状态（如反应物、产物、过渡态）的能量，我们可以预测化学过程的能量变化。例如，**反应能垒**（或活化能）就是过渡态与反应物之间的能量差，它决定了反应速率的快慢。

![Figure 3: 反应能垒与电荷分布示意图](https://github.com/BruceJackey/Picture/blob/main/Tutorial_step1_step2_bundle_pyscf_pyscf_fig3.png?raw=true)

*Figure 3. 左侧为反应能量剖面图，展示了反应物经过过渡态生成产物的能量变化，其中反应能垒决定了反应速率。右侧以水分子为例，展示了两种描述电荷分布的方式：宏观的电偶极矩矢量和微观的原子部分电荷。*


#### 基组

在实际计算中，分子的轨道是通过一组数学函数——**基函数（Basis Functions）**——线性组合来近似描述的。这些通常以原子为中心的函数集合被称为**基组（Basis Set）**（例如 `sto-3g`, `6-31g`, `cc-pVDZ`）。基组的选择直接影响计算的精度和成本：更大的基组通常能更准确地描述电子分布，但计算也更耗时。因此，选择合适的基组是在计算精度和效率之间进行权衡的关键一步。

![Figure 4: 基组对分子轨道的近似示意图](https://github.com/BruceJackey/Picture/blob/main/Tutorial_step1_step2_bundle_pyscf_pyscf_fig4.png?raw=true)

*Figure 4. 分子轨道（左）由以原子为中心的基函数线性组合来近似。使用较少的基函数（小基组，右上）计算速度快但精度较低。使用更多的基函数（大基组，右下）可以更精确地描述电子分布，但计算成本也更高。*


### 参考链接

- [分子电学性质：偶极矩与极化率](https://www.bohrium.com/sciencepedia/article/823590)
- [部分原子电荷](https://www.bohrium.com/sciencepedia/article/775161)
- [几何优化与梯度算法](https://www.bohrium.com/sciencepedia/article/820579)
- [势垒](https://www.bohrium.com/sciencepedia/article/803805)

## 导入工具集

PySCF (Python-based Simulations of Chemistry Framework) 是一个基于 Python 的开源计算化学框架，广泛用于电子结构计算。它支持 Hartree-Fock (HF)、密度泛函理论 (DFT)、耦合簇 (CC) 等多种方法。

以下是本教程所需的模块导入：

*   `numpy`: 用于处理矩阵运算和数组操作，PySCF 的底层数据结构基于 NumPy。
*   `pyscf`: 主包。
*   `pyscf.gto`: **G**aussian **T**ype **O**rbitals，用于定义分子几何结构、基组和物理属性（如电荷、自旋）。
*   `pyscf.scf`: **S**elf-**C**onsistent **F**ield，包含 Hartree-Fock (RHF, UHF, ROHF) 方法的实现。
*   `pyscf.dft`: **D**ensity **F**unctional **T**heory，包含 Kohn-Sham DFT (RKS, UKS) 方法的实现。
*   `pyscf.grad`: 用于计算解析核梯度（Nuclear Gradients），即能量对原子坐标的导数，用于几何优化或受力分析。


```python
import numpy as np
import pyscf
from pyscf import gto, scf, dft, grad

# 打印版本信息以确认安装成功
print(f"PySCF Version: {pyscf.__version__}")
print(f"NumPy Version: {np.__version__}")
```

    PySCF Version: 2.11.0
    NumPy Version: 2.2.6
    

## 工具介绍

PySCF 的核心工作流通常遵循以下步骤：
1.  **定义分子 (Mole)**：指定原子坐标、元素、基组、电荷和自旋。
2.  **选择方法 (Mean-Field)**：选择近似方法（如 RHF 或 DFT）并构建均值场对象。
3.  **运行计算 (Kernel)**：执行自洽场迭代以求解电子结构，获得总能量和波函数。
4.  **后处理 (Analysis)**：基于收敛的波函数计算性质，如偶极矩、原子电荷、核梯度等。

### 1. 分子构建

`pyscf.gto.Mole` 是 PySCF 中最基础的类，用于存储分子的几何和物理信息。构建分子时，通常使用 `build()` 方法来初始化积分所需的辅助数据。

下面的封装函数 `create_molecule` 简化了分子的初始化过程。


```python
def create_molecule(atom_str, basis='sto-3g', charge=0, spin=0, verbose=0):
    """
    使用 PySCF 构建分子对象。

    参数:
    - atom_str (str): 原子坐标字符串。格式例如 "O 0 0 0; H 0 1 0; H 0 0 1" 或 "O 0 0 0\nH 0 1 0"。
                  默认单位为 Angstrom (埃)。
    - basis (str): 基组名称，例如 'sto-3g', '6-31g', 'cc-pvdz'。
    - charge (int): 分子总电荷数 (默认 0)。
    - spin (int): 分子自旋多重度 2S (即未配对电子数)。
                  例如，单重态 spin=0，二重态 spin=1。
    - verbose (int): 输出详细程度，0 为静默，4 为详细。

    返回:
    - mol (pyscf.gto.Mole): 已构建的 PySCF 分子对象。
    """
    mol = gto.Mole()
    mol.atom = atom_str
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.verbose = verbose
    mol.build()
    return mol
```

### 2. 能量计算 (SCF & DFT)

PySCF 通过 `scf` 和 `dft` 模块提供能量计算功能。对于闭壳层体系，通常使用 `RHF` (Restricted Hartree-Fock) 或 `RKS` (Restricted Kohn-Sham DFT)。

*   **RHF**: 适用于不包含电子相关能的平均场计算。
*   **DFT (RKS)**: 通过交换-关联泛函 (XC functional) 引入电子相关效应。常用的泛函包括 `b3lyp`, `pbe0` 等。

下面的封装函数 `run_calculation` 统一了 HF 和 DFT 的调用接口。


```python
def run_calculation(mol, method='rhf', xc='b3lyp'):
    """
    运行单点能计算 (SCF 或 DFT)。

    参数:
    - mol (pyscf.gto.Mole): 已构建的分子对象。
    - method (str): 计算方法，可选 'rhf' (Hartree-Fock) 或 'dft' (密度泛函理论)。
    - xc (str): 交换关联泛函 (仅当 method='dft' 时有效)，例如 'b3lyp', 'pbe'。

    返回:
    - mf (pyscf.scf.hf.SCF): 收敛后的均值场对象，包含波函数、轨道系数等信息。
    - energy (float): 计算得到的总能量 (单位: Hartree)。
    """
    method = method.lower()
    
    if method == 'rhf':
        mf = scf.RHF(mol)
    elif method == 'dft':
        mf = dft.RKS(mol)
        mf.xc = xc
    else:
        raise ValueError(f"Unsupported method: {method}. Use 'rhf' or 'dft'.")
    
    # 运行计算
    energy = mf.kernel()
    return mf, energy
```

### 3. 性质分析：梯度、偶极矩与电荷

计算完成后，我们通常需要从均值场对象 (`mf`) 中提取物理性质：

*   **核梯度 (Nuclear Gradients)**: 能量对原子坐标的一阶导数，用于受力分析。在 PySCF 中，可以通过 `mf.nuc_grad_method()` 获取梯度计算器。
*   **Mulliken 电荷分析**: 将电子密度分配到原子上的一种方法。PySCF 提供了 `analyze` 方法或 `mulliken_meta` 函数。
*   **偶极矩 (Dipole Moment)**: 描述分子正负电荷中心的偏离程度。`mf.dip_moment()` 可以直接计算。

以下封装函数提供了这些常用分析功能的简便接口。


```python
def compute_gradients(mf):
    """
    计算分子的核梯度 (受力)。

    参数:
    - mf (pyscf.scf.hf.SCF): 已收敛的均值场对象。

    返回:
    - gradients (numpy.ndarray): 形状为 (N_atoms, 3) 的数组，表示每个原子的梯度向量 (单位: Hartree/Bohr)。
    """
    # 获取对应的梯度计算方法对象 (如 RHF 对应 grad.RHF, RKS 对应 grad.RKS)
    grad_calc = mf.nuc_grad_method()
    gradients = grad_calc.kernel()
    return gradients

def perform_mulliken_analysis(mf):
    """
    执行 Mulliken 布居分析并返回原子电荷。

    参数:
    - mf (pyscf.scf.hf.SCF): 已收敛的均值场对象。

    返回:
    - charges (numpy.ndarray): 每个原子的 Mulliken 部分电荷数组。
    """
    # analyze() 方法会打印详细信息，并返回 (布居, 电荷)
    # 我们只需要返回电荷部分
    pop, charges = mf.analyze(verbose=0)
    return charges

def get_dipole_moment(mf):
    """
    计算分子的偶极矩。

    参数:
    - mf (pyscf.scf.hf.SCF): 已收敛的均值场对象。

    返回:
    - dipole_vec (numpy.ndarray): 偶极矩向量 [x, y, z] (单位: Debye)。
    """
    # mf.dip_moment() 默认返回 Debye 单位的向量
    dipole_vec = mf.dip_moment(verbose=0)
    return dipole_vec
```

## 综合应用

本部分将结合前文定义的工具函数，解决几个实际的计算化学问题。我们将涵盖分子偶极矩计算、原子电荷分析、核梯度计算以及反应能垒的求解。每个案例都对应一个具体的化学场景，旨在展示 PySCF 在不同任务中的应用流程。

### 1. 水分子偶极矩计算

**理论背景说明**

分子偶极矩（Dipole Moment, $\vec{\mu}$）衡量了分子内部正负电荷中心的分离程度，是描述分子极性的关键物理量。总偶极矩由核电荷贡献（$\vec{\mu}_{nuc}$）和电子贡献（$\vec{\mu}_{elec}$）两部分矢量和组成。在 Hartree-Fock 近似下，电子部分通过密度矩阵与单电子偶极积分的迹（Trace）求得。本案例将计算水分子在不同基组下的偶极矩大小，以观察基组对电子分布描述的影响。


![Figure 5: 水分子偶极矩示意图](https://github.com/BruceJackey/Picture/blob/main/Tutorial_step1_step2_bundle_pyscf_pyscf_fig5.png?raw=true)

*Figure 5. 水分子中，由于氧的电负性更强，电子云偏向氧原子，导致正负电荷中心分离。这个电荷分离的程度由偶极矩向量 μ⃗ 来描述，其方向从正电荷中心指向负电荷中心。*

**目标**：使用 RHF 方法计算水分子在 `sto-3g`, `6-31g`, `cc-pVDZ` 三种基组下的总偶极矩模长。

**代码实现**



```python
# 定义水分子坐标 (Angstrom)
water_coords = """
O 0.000000 0.000000 0.117300
H 0.000000 0.757200 -0.469200
H 0.000000 -0.757200 -0.469200
"""

# 待测试的基组列表
basis_sets = ['sto-3g', '6-31g', 'cc-pvdz']
dipole_magnitudes = []

for basis in basis_sets:
    # 1. 构建分子对象
    mol = create_molecule(water_coords, basis=basis, verbose=0)
    
    # 2. 运行 RHF 计算
    # 偶极矩计算依赖于收敛的密度矩阵
    mf, _ = run_calculation(mol, method='rhf')
    
    # 3. 获取偶极矩向量 (单位: Debye)
    # get_dipole_moment 是我们在第一部分定义的封装函数
    dip_vec = get_dipole_moment(mf)
    
    # 4. 计算向量模长
    magnitude = np.linalg.norm(dip_vec)
    dipole_magnitudes.append(round(magnitude, 4))

# 输出结果列表
print(dipole_magnitudes)
```

    [np.float64(1.7253), np.float64(2.631), np.float64(2.0574)]
    

### 2. 甲醛分子 Mulliken 电荷分析

**理论背景说明**

虽然原子电荷不是量子力学中的可观测量，但它为理解分子内的电子分布提供了一个直观的化学图像。Mulliken 布居分析（Population Analysis）是一种经典方法，它基于基函数将电子密度划分归属到各个原子上。本案例将使用密度泛函理论（DFT）中的 B3LYP 泛函，计算甲醛分子（$CH_2O$）中各原子的部分电荷，这有助于理解该分子的亲电/亲核反应位点。


![Figure 6: 甲醛分子的Mulliken电荷分布](https://github.com/BruceJackey/Picture/blob/main/Tutorial_step1_step2_bundle_pyscf_pyscf_fig6.png?raw=true)

*Figure 6. 该图展示了甲醛分子（CH₂O）的结构及其各原子的Mulliken部分电荷。氧原子因其高电负性而带有显著的负电荷（δ-），而碳原子和氢原子则带有正电荷（δ+），这揭示了分子中潜在的亲电和亲核反应位点。*

**目标**：使用 B3LYP 泛函计算甲醛分子在 `sto-3g` 和 `6-31g` 基组下的 Mulliken 原子电荷。

**代码实现**



```python
# 定义甲醛分子坐标 (Angstrom)
# 原子顺序: C, H, H, O
formaldehyde_coords = """
C 0.0000 0.0000 -0.5316
H 0.0000 0.9386 -1.1186
H 0.0000 -0.9386 -1.1186
O 0.0000 0.0000 0.6724
"""

test_cases = ['sto-3g', '6-31g']
all_charges = []

for basis in test_cases:
    # 1. 构建分子
    mol = create_molecule(formaldehyde_coords, basis=basis, verbose=0)
    
    # 2. 运行 DFT (B3LYP) 计算
    mf, _ = run_calculation(mol, method='dft', xc='b3lyp')
    
    # 3. 执行 Mulliken 分析获取电荷
    charges = perform_mulliken_analysis(mf)
    
    # 格式化保留4位小数
    formatted_charges = [round(c, 4) for c in charges]
    all_charges.append(formatted_charges)

# 输出结果: [[C, H, H, O] (sto-3g), [C, H, H, O] (6-31g)]
print(all_charges)
```

    [[np.float64(0.0), np.float64(-0.0), np.float64(-1.2931)], [np.float64(0.0), np.float64(-0.0), np.float64(-2.3612)]]
    

### 3. 氨分子非平衡构型的核梯度计算

**理论背景说明**

核梯度（Nuclear Gradient）是总能量对原子核坐标的一阶导数（$\nabla_i E$），其负值即为作用在原子上的力（$\vec{F}_i = -\nabla_i E$）。在平衡几何结构处，核梯度为零；而在非平衡结构处，梯度的方向指示了能量下降最快的方向，这正是几何优化算法寻找稳定结构的基础。本案例将计算氨分子（$NH_3$）在三种非平衡构型下的 RMS（均方根）梯度，量化分子受力的大小。


![Figure 7: 势能面、核梯度与力的关系](https://github.com/BruceJackey/Picture/blob/main/Tutorial_step1_step2_bundle_pyscf_pyscf_fig7.png?raw=true)

*Figure 7. 该图展示了原子在势能面上的运动。在非平衡构型处，总能量对原子坐标的梯度（∇E）不为零，其方向指向能量上升最快的方向。原子受到的力（F）与梯度方向相反，指向能量最低点（平衡构型），驱动分子进行几何优化。*

**目标**：计算氨分子在平面、键拉伸、角压缩三种构型下的 RMS 梯度（单位：Hartree/Bohr）。

**代码实现**



```python
# 定义三种测试构型 (Angstrom)
geometries = [
    # Case 1: Planar Geometry (平面构型)
    """
    N 0.00000000 0.00000000 0.00000000
    H 0.00000000 0.98188106 0.00000000
    H 0.85037119 -0.49094053 0.00000000
    H -0.85037119 -0.49094053 0.00000000
    """,
    # Case 2: N-H bond stretched (键拉伸)
    """
    N 0.00000000 0.00000000 0.00000000
    H 0.00000000 0.00000000 1.21000000
    H 1.01890173 0.00000000 -0.40333333
    H -0.50945086 0.88239035 -0.40333333
    """,
    # Case 3: H-N-H angle compressed (角压缩)
    """
    N 0.00000000 0.00000000 0.00000000
    H 0.00000000 0.00000000 1.01000000
    H 0.90000000 0.00000000 -0.45000000
    H -0.45000000 0.77942286 -0.45000000
    """
]

rms_gradients = []

for geom in geometries:
    # 1. 构建分子 (统一使用 sto-3g 基组)
    mol = create_molecule(geom, basis='sto-3g', verbose=0)
    
    # 2. 运行 RHF 计算
    mf, _ = run_calculation(mol, method='rhf')
    
    # 3. 计算解析梯度
    grads = compute_gradients(mf)
    
    # 4. 计算 RMS Gradient
    # 公式: g_rms = sqrt( sum(g^2) / 3N )
    N_atoms = mol.natm
    sum_sq = np.sum(grads**2)
    g_rms = np.sqrt(sum_sq / (3 * N_atoms))
    
    rms_gradients.append(round(g_rms, 6))

print(rms_gradients)
```

    [np.float64(0.014709), np.float64(0.048526), np.float64(0.019964)]
    

### 4. 氨分子翻转能垒计算

**理论背景说明**

化学反应或构象转变通常需要越过一个能量势垒，即过渡态（Transition State）与基态（Ground State）之间的能量差。氨分子的氮原子翻转（Nitrogen Inversion）是一个经典的量子力学现象，涉及从角锥形基态经过平面形过渡态的过程。本案例将利用 DFT 方法计算这一翻转过程的能垒，即 $E_{barrier} = E_{TS} - E_{GS}$。


![Figure 8: 氨分子翻转能垒示意图](https://github.com/BruceJackey/Picture/blob/main/Tutorial_step1_step2_bundle_pyscf_pyscf_fig8.png?raw=true)

*Figure 8. 该图展示了氨分子（NH3）翻转的能量变化过程。分子从角锥形的基态（GS）出发，经过一个能量较高的平面形过渡态（TS），转变为另一个等价的角锥形基态。基态与过渡态之间的能量差即为翻转能垒（E_barrier）。*

**目标**：分别使用 `B3LYP/6-31g*` 和 `PBE0/cc-pVDZ` 计算氨分子的翻转能垒（单位：kcal/mol）。

**代码实现**



```python
# 定义几何结构 (Angstrom)
# 1. 角锥形基态 (Ground State)
geom_gs = """
N   0.000000    0.000000    0.382000
H   0.000000    0.937000   -0.127333
H   0.811475   -0.468500   -0.127333
H  -0.811475   -0.468500   -0.127333
"""

# 2. 平面形过渡态 (Transition State)
geom_ts = """
N   0.000000    0.000000    0.000000
H   0.000000    1.008000    0.000000
H   0.872954   -0.504000    0.000000
H  -0.872954   -0.504000    0.000000
"""

# 测试用例: (functional, basis)
test_cases = [
    ('b3lyp', '6-31g*'),
    ('pbe0', 'cc-pvdz')
]

barriers = []
hartree_to_kcal = 627.5095

for xc, basis in test_cases:
    # 1. 计算基态能量
    mol_gs = create_molecule(geom_gs, basis=basis, verbose=0)
    _, e_gs = run_calculation(mol_gs, method='dft', xc=xc)
    
    # 2. 计算过渡态能量
    mol_ts = create_molecule(geom_ts, basis=basis, verbose=0)
    _, e_ts = run_calculation(mol_ts, method='dft', xc=xc)
    
    # 3. 计算能垒并转换为 kcal/mol
    barrier = (e_ts - e_gs) * hartree_to_kcal
    barriers.append(round(barrier, 4))

print(barriers)
```

    [np.float64(2.8047), np.float64(4.4522)]
    

## 反馈（帮助我们持续改进教程）

欢迎你在完成本 Tutorial 后进行反馈，告诉我们哪些地方好用/不好用、还希望增加哪些内容。
