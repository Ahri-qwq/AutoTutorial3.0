# Step 1 调研结果

## 1.1 RAG 检索结果总结

### 核心流程

**SCF → NSCF 两步框架（来源：abacus_user_guide__abacus_dos.html.md、2024秋计算材料学）**
- Step 1：SCF 自洽计算，设置 `out_chg 1`，输出 SPIN1_CHG.cube
- Step 2：NSCF 非自洽，设置 `init_chg file`，读入电荷密度，固定不更新
- NSCF 时必须设置 `symmetry 0`，否则 k 点会被对称性折叠
- 能带和 DOS 各用不同的 KPT（能带=Line模式，DOS=密网格）

---

### Si 能带案例（来源：2024秋计算材料学）

**SCF INPUT：**
```
INPUT_PARAMETERS
symmetry                1
calculation             scf
ecutwfc                 50
scf_thr                 1e-7
scf_nmax                300
basis_type              lcao
out_chg                 1
```

**NSCF INPUT（能带）：**
```
INPUT_PARAMETERS
symmetry                0
init_chg                file
calculation             nscf
ecutwfc                 50
scf_thr                 1e-7
scf_nmax                300
basis_type              lcao
out_band                1
out_bandgap             1
```

**KPT-scf（9x9x9 Gamma）：**
```
K_POINTS
0
Gamma
9 9 9 0 0 0
```

**KPT-nscf（能带路径，Si FCC 原胞布里渊区）：**
```
K_POINTS
8
Line
0.000 0.000 0.000  30  # G
0.500 0.000 0.500  20  # X
0.625 0.250 0.625  20  # U
0.375 0.375 0.750  50  # K
0.000 0.000 0.000  50  # G
0.500 0.500 0.500  20  # L
0.500 0.250 0.750  20  # W
0.500 0.000 0.500   1  # X
```

**关键数值：** Si PBE 带隙 ~0.57 eV（间接带隙，价带顶在 Γ，导带底在 Δ 方向即 G-X 之间）

---

### Al DOS/PDOS 案例（来源：abacus_user_guide__abacus_dos.html.md）

**STRU（FCC 4原子单胞）：**
```
ATOMIC_SPECIES
Al 26.982 Al_ONCV_PBE-1.0.upf upf201
LATTICE_CONSTANT
1.88972612546
LATTICE_VECTORS
4.0450551637 0.0000000000 0.0000000000 #latvec1
0.0000000000 4.0450551637 0.0000000000 #latvec2
0.0000000000 0.0000000000 4.0450551637 #latvec3
ATOMIC_POSITIONS
Direct
Al #label
0 #magnetism
4 #number of atoms
0.0000000000 0.0000000000 0.0000000000 m 0 0 0
0.5000000000 0.5000000000 0.0000000000 m 0 0 0
0.5000000000 0.0000000000 0.5000000000 m 0 0 0
0.0000000000 0.5000000000 0.5000000000 m 0 0 0
```

注：实际计算需转换为原胞（2原子→1原子），用 atomkit 命令转换。

**Al SCF INPUT（nspin=1，金属 smearing）：**
参数参考：ecutwfc 60，smearing_method gauss，smearing_sigma 0.02（金属需较大展宽），
mixing_type broyden，mixing_beta 0.4，out_chg 1

**Al NSCF INPUT（DOS）：**
增加：init_chg file，out_dos 1（或 2 同时输出 PDOS），dos_sigma 0.07（默认值），
将 symmetry 设为 0

**Al KPT 能带路径（FCC 原胞，来源：官方文档 KLINES）：**
```
K_POINTS
8
Line
0.00000000 0.00000000 0.00000000 25 # GAMMA
0.50000000 0.00000000 0.50000000  9 # X
0.62500000 0.25000000 0.62500000  1 # U
0.37500000 0.37500000 0.75000000 27 # K
0.00000000 0.00000000 0.00000000 22 # GAMMA
0.50000000 0.50000000 0.50000000 18 # L
0.50000000 0.25000000 0.75000000 12 # W
0.50000000 0.00000000 0.50000000  1 # X
```

**费米能级提取：**
```bash
grep EFERMI OUT.*/running_nscf.log
# 示例输出：EFERMI = 10.963171515 eV
```
Al 金属 EFERMI ≈ 10.963 eV（参考值）

---

### Atomkit 后处理命令（来源：多个教程）

**能带绘图（abacus-plot -b）：**
```json
{
    "bandfile": "BANDS_1.dat",
    "efermi": 10.963171515,
    "energy_range": [-6, 6],
    "bandfig": "band.png",
    "kptfile": "KPT",
    "dpi": 300
}
```
```bash
abacus-plot -b
```

**TDOS 绘图（abacus-plot -d）：**
```json
{
    "tdosfile": "TDOS",
    "efermi": 10.963171515,
    "energy_range": [-6, 6],
    "dos_range": [0, 5],
    "tdosfig": "tdos.png",
    "dpi": 300
}
```
```bash
abacus-plot -d
```

**PDOS 绘图（abacus-plot -d -p -o）：**
```json
{
    "pdosfile": "PDOS",
    "efermi": 10.963171515,
    "energy_range": [-6, 6],
    "dos_range": [0, 5],
    "figsize": [14, 10],
    "species": {"Al": [0, 1, 2]},
    "pdosfig": "pdos.png"
}
```
```bash
abacus-plot -d -p -o
```
注：species 中 [0,1,2] 分别对应 s、p、d 轨道

**PDOS 详细说明（来源：abacus_user_guide__abacus_pdos.html.md）：**
- `out_dos 2`：同时输出 DOS 和 PDOS（xml 格式）
- `abacus-plot -d -o`：只输出分元素 PDOS
- `abacus-plot -d -p -o`：输出分轨道 PDOS
- efermi 从 running_scf.log 中 EFERMI 关键字获取

---

## 1.2 风格总结（参考：磁性材料教程、DeePMD-kit教程）

- **文章结构**：一、介绍→二、准备→三、SCF 计算→四、能带计算（Si）→五、DOS/PDOS 计算（Al）→六、总结
- **开头方式**：直接介绍功能（"电子的态密度... 在凝聚态物理... 有重要用途"），无"在当今..."铺垫
- **代码块**：文件名注释、参数注释简洁，必要说明用注释行
- **参数说明**：通常在代码块后用编号列表说明关键参数
- **语言特征**：技术直接，"注意这里"提醒关键点，句式简洁
- **长度**：参考文章约 250-350 行（中等长度），本教程目标 500-700 行（含2个案例，内容更多）
