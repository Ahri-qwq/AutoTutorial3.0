## 四、运行与结果分析

### 4.1 运行计算

进入工作目录后，执行：

```bash
cd ABACUS_DFT+U
OMP_NUM_THREADS=1 mpirun -np 8 abacus
```

计算完成后，输出文件位于 `OUT.NiO/` 目录。

---

### 4.2 关键结果提取

#### 总能量

在 `OUT.NiO/running_scf.log` 中搜索 `FINAL_ETOT`：

```
--------------------------------------------
!FINAL_ETOT_IS -9255.7279034240546025 eV
--------------------------------------------
```

#### 能隙

由于设置了 `out_bandgap = 1`，搜索 `E_bandgap`，取最后一次出现的值：

```
E_bandgap    +0.205369322748    +2.794192983776
```

两列分别对应 spin up 和 spin down 通道的能隙（单位 eV）。NiO 是反铁磁绝缘体，spin down 通道的能隙约 2.79 eV，与实验相比已有明显改善（纯 GGA 通常给出远小的值甚至金属态）。

#### 磁矩

搜索 `absolute magnetism`，取最后一次出现：

```
      total magnetism (Bohr mag/cell) = 0.00000000
   absolute magnetism (Bohr mag/cell) = 3.35321634
```

- **total magnetism = 0**：符合反铁磁体系整体无净磁矩的预期
- **absolute magnetism ≈ 3.35**：体系内存在大小相等、方向相反的局域磁矩

**原子分辨磁矩（Mulliken 分析）：**

由于设置了 `out_mul = 1`，查看 `OUT.NiO/mulliken.txt`，搜索 `Magnetism`：

```
Total Magnetism on atom  Ni1     1.8268646
Total Magnetism on atom  Ni2    -1.8268646
Total Magnetism on atom  O      -3.6718263e-13
Total Magnetism on atom  O       1.7330755e-13
```

两个 Ni 原子的磁矩大小相等（约 ±1.83 μB）、方向相反，O 原子磁矩接近于零，确认得到了反铁磁基态。

---

### 4.3 Occupation Matrix 的读取与解析

DFT+U 计算的一个特色输出是 **occupation matrix**，它记录了每个 +U 原子上 d 轨道的占据情况，是理解 +U 修正效果的重要信息。

在 `running_scf.log` 中搜索以 `L(S)DA+U` 开头的块，格式如下：

**头部信息**（说明哪些 species 施加了 +U 及其参数）：

```
atom_type=0  L=2  chi=0    U=5ev
atom_type=1  L=2  chi=0    U=5ev
```

与 INPUT 文件中的设定一致：atom_type 0（Ni1）和 1（Ni2）均在 l=2（d 轨道）上施加 U=5 eV。

**各原子的 occupation matrix**：

```
atoms  0          // 原子编号：第一个 Ni 原子（Ni1）
L  2              // +U 的 l 量子数（d 轨道）
zeta  0           // 基组 zeta 指标；Ni 只有一个 d 基组，故 zeta=0
spin  0           // spin up 的 5×5 矩阵
// [5×5 matrix]
spin  1           // spin down 的 5×5 矩阵
// [5×5 matrix]

atoms  1          // 第二个 Ni 原子（Ni2）
L  2
zeta  0
spin  0
// [5×5 matrix]
spin  1
// [5×5 matrix]
```

d 轨道共有 5 个（m = -2, -1, 0, +1, +2），因此每个自旋通道的 occupation matrix 是 5×5 的实对称矩阵。对角元代表各 d 轨道的占据数，0 到 1 之间。

**onsite.dm 输出文件：**

由于设置了 `out_chg 1`，ABACUS 会将最后一步 SCF 收敛后的 occupation matrix 写入 `OUT.NiO/onsite.dm`，文件格式与日志中的输出完全相同。这个文件在使用 occupation matrix control 功能时会用到（见下一章）。
