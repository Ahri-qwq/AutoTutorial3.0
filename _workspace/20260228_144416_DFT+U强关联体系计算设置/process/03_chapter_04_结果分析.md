## 4. 结果分析

计算完成后，所有输出文件位于 `OUT.NiO/` 目录下。

### 4.1 总能量

在 `running_scf.log` 中搜索 `FINAL_ETOT`：

```bash
grep "FINAL_ETOT" OUT.NiO/running_scf.log
```

预期输出：

```
 --------------------------------------------
 !FINAL_ETOT_IS -9255.7279034240546025 eV
 --------------------------------------------
```

### 4.2 能隙

INPUT 中设置了 `out_bandgap 1`，在 `running_scf.log` 中搜索 `E_bandgap`，取最后一次出现的值：

```bash
grep "E_bandgap" OUT.NiO/running_scf.log | tail -1
```

预期输出：

```
E_bandgap               +0.205369322748               +2.794192983776
```

两列分别是 spin-up 和 spin-down 的能隙（单位：eV）。NiO 是绝缘体，DFT+U 给出了非零能隙。若使用纯 GGA，NiO 通常给出 0 能隙（金属态），这正是 DFT+U 修正的意义所在。

> spin-up 能隙（0.21 eV）和 spin-down 能隙（2.79 eV）不对称，反映了两种自旋通道中能带结构的差异。

### 4.3 磁矩

**总磁矩与绝对磁矩**

搜索 `absolute magnetism`，取最后一次出现的值：

```bash
grep -A1 "total magnetism" OUT.NiO/running_scf.log | tail -4
```

预期输出：

```
      total magnetism (Bohr mag/cell) = 0.00000000
   absolute magnetism (Bohr mag/cell) = 3.35321634
```

- 总磁矩为 0：验证了体系是反铁磁态
- 绝对磁矩 3.35：两个 Ni 原子磁矩的绝对值之和

**原子磁矩（Mulliken 分析）**

INPUT 中设置了 `out_mul 1`，生成 `OUT.NiO/Mulliken.txt`。搜索 `Magnetism`：

```bash
grep "Total Magnetism" OUT.NiO/Mulliken.txt
```

预期输出：

```
Total Magnetism on atom  Ni1           1.8268646
Total Magnetism on atom  Ni2          -1.8268646
Total Magnetism on atom  O      -3.6718263e-13
Total Magnetism on atom  O       1.7330755e-13
```

Ni1 和 Ni2 的磁矩大小相等（约 1.83 μ_B）、方向相反，O 原子的磁矩几乎为零。这是反铁磁 NiO 的预期结果。

### 4.4 Occupation Matrix 的读取

DFT+U 计算会在 `running_scf.log` 中的每一步输出 occupation matrix。搜索以 `L(S)DA+U` 开头的块：

```bash
grep -A 50 "L(S)DA+U" OUT.NiO/running_scf.log | head -60
```

该块的结构如下：

```
L(S)DA+U:
atom_type=0  L=2  chi=0    U=5ev    # Ni1：d轨道，U=5 eV
atom_type=1  L=2  chi=0    U=5ev    # Ni2：d轨道，U=5 eV

atoms  0          # 第一个Ni原子（Ni1）
L  2              # d 轨道（l=2）
zeta  0           # 第一个（也是唯一一个）d基组
spin  0           # spin-up 的 occupation matrix（5×5）
  0.9xx  ...
  ...             # 5行5列
spin  1           # spin-down 的 occupation matrix
  0.0xx  ...
  ...

atoms  1          # 第二个Ni原子（Ni2）
L  2
zeta  0
spin  0
  ...
spin  1
  ...
```

occupation matrix 的物理含义：每个元素 n_{mm'} 表示 d 轨道磁量子数 m 和 m' 之间的占据，对角元代表各 d 分量的占据数（0~1）。Ni1 的 spin-up 分量占据较满，Ni2 的 spin-down 分量占据较满，体现了两者的磁矩相反。

**`onsite.dm` 文件**

计算结束后（因设置了 `out_chg 1`），在 `OUT.NiO/` 中生成 `onsite.dm` 文件，保存最后一步的 occupation matrix，格式与上述 log 中的输出相同。该文件在下一节 omc 演示中使用。
