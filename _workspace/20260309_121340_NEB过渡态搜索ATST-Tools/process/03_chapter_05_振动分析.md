## 第五章：振动分析与过渡态验证

### 5.1 为什么需要振动分析

NEB 计算给出能量最高的映像作为过渡态候选，但这只是能量上的判断。严格意义上的过渡态必须满足：**在反应坐标方向存在且仅存在一个虚频**。

振动分析（频率计算）的作用：
1. **验证过渡态正确性**：唯一虚频对应反应坐标方向的不稳定模式
2. **排查错误**：若虚频不止一个，说明过渡态结构不对
3. **零点能（ZPE）校正**：将电子能量转换为 0 K 下的振动包含能量
4. **有限温度自由能校正**：通过谐振近似计算振动熵和自由能

### 5.2 vib_analysis.py 配置

ATST-Tools 的 `vibration/vib_analysis.py` 使用 ASE 的有限差分方法，对指定原子施加微小位移，通过力的变化计算 Hessian 矩阵，进而得到振动模式。

**从 NEB 结果中提取过渡态和活跃原子：**

```python
from ase.io import read
from neb2vib import neb2vib

# 从 NEB 轨迹中自动识别过渡态结构和活跃原子
neb_traj = read('neb_latest.traj', index=':')
atoms, vib_indices = neb2vib(neb_traj)
```

`neb2vib` 函数自动选取 NEB 链中位移最大的原子作为振动计算对象，无需手动指定原子序号。也可以手动指定：

```python
# 手动指定：只计算 H 原子（序号 0, 1）和 C 原子（序号 37）
vib_indices = [0, 1, 37]
```

完整脚本（针对案例二配置）：

```python
# vib_analysis.py — Cy-Pt@graphene 过渡态振动分析

from ase.io import read
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.parallel import world
from neb2vib import neb2vib
import os

# ===== 读取 NEB 结果 =====
neb_traj = read('neb_latest.traj', index=':')
atoms, vib_indices = neb2vib(neb_traj)

T = 523.15  # K（振动自由能校正温度）

# ===== ABACUS 计算设置（与 NEB 保持一致）=====
abacus = "abacus"
mpi = 16
omp = 4
# ... pp, basis, parameters 与 autoneb_run.py 相同 ...

# ===== 振动计算参数 =====
vib_name = 'vib'
delta    = 0.01   # 有限位移步长（Å）
nfree    = 2      # 每个方向位移次数（正+负）

def set_calculator(abacus, parameters, mpi=1, omp=1):
    os.environ['OMP_NUM_THREADS'] = f'{omp}'
    profile = AbacusProfile(f"mpirun -np {mpi} {abacus}")
    out_directory = f"SCF-rank{world.rank}"
    return Abacus(profile=profile, directory=out_directory, **parameters)

if __name__ == "__main__":
    atoms.calc = set_calculator(abacus, parameters, mpi=mpi, omp=omp)
    vib = Vibrations(atoms, indices=vib_indices,
                     name=vib_name, delta=delta, nfree=nfree)
    vib.run()
    vib.summary()
    vib.write_mode()

    # 热力学校正
    vib_energies = vib.get_energies()
    thermo = HarmonicThermo(vib_energies, ignore_imag_modes=True)
    entropy     = thermo.get_entropy(T)
    free_energy = thermo.get_helmholtz_energy(T)
    print(f"Entropy:     {entropy:.6e} eV/K")
    print(f"Free Energy: {free_energy:.6f} eV")
```

运行：

```bash
python vib_analysis.py
```

### 5.3 计算结果解读

案例二过渡态的振动分析结果：

```
  #    meV     cm⁻¹
---------------------
  0   87.4i   705.0i   ← 唯一虚频：反应坐标方向（C-H 断裂）
  1    2.4     19.3
  2    3.5     28.3
  ...
 53  377.5   3044.9
---------------------
Zero-point energy: 4.416 eV
```

| 结果 | 数值 | 含义 |
|------|------|------|
| 虚频（imaginary） | 87.4i meV（705.0i cm⁻¹） | 对应 C-H 键断裂的反应坐标 |
| 虚频个数 | 1 | ✓ 过渡态正确（唯一虚频） |
| 零点能（ZPE） | 4.416 eV | 包含所有实振动模式的零点贡献 |
| 振动熵（523 K） | 见输出 | 用于 Helmholtz 自由能校正 |

**结果判读：**
- 虚频 705.0i cm⁻¹，对应标准 C-H 伸缩频率区间，物理上合理
- 只有 1 个虚频，确认该映像为真正的过渡态
- 查看 `vib.0.traj` 可在 ASE-GUI 中动画展示该虚频振动模式

```bash
ase -T gui vib.0.traj
```

计算完成后，`vib/` 目录中存储了所有位移结构的 DFT 结果（JSON 格式），下次修改温度重新分析时无需重新调用 ABACUS。
