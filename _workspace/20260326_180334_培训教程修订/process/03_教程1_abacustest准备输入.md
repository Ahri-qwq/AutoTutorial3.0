# 使用 abacustest 准备输入文件

abacustest 是 ABACUS 的前后处理工具，可以快速从结构文件（CIF、POSCAR 等）生成完整的 ABACUS 输入文件夹。

## 安装 abacustest

```bash
# 方法1：通过 pip 安装
pip install abacustest

# 方法2：从源码安装
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

安装后验证：
```bash
abacustest -h
```

---

## 从单个结构文件准备输入

### 基本用法

假设你有一个 CIF 文件 `Si.cif`，想生成 ABACUS 输入文件：

```bash
abacustest model inputs prepare \
  -f Si.cif \
  --ftype cif \
  --pp /path/to/pseudopotentials \
  --orb /path/to/orbitals
```

**参数说明：**
- `-f`：结构文件路径（支持 CIF、POSCAR、XYZ 等）
- `--ftype`：文件类型（cif、poscar、xyz 等）
- `--pp`：赝势库目录
- `--orb`：轨道库目录（LCAO 计算需要，PW 计算可省略）

**赝势和轨道库要求：**
- 文件名必须以元素名开头，例如 `Si_ONCV_PBE-1.0.upf`、`Si_gga_8au_60Ry_2s2p1d.orb`
- 或者在目录下提供 `element.json` 文件，格式：`{"Si": "Si_ONCV_PBE-1.0.upf"}`

### 生成结果

执行后会在当前目录生成 `Si/` 文件夹，包含：
- `STRU`：结构文件
- `INPUT`：默认参数文件
- `KPT`：默认 k 点文件
- 赝势和轨道文件的软链接或副本

---

## 批量准备多个任务

对于多个结构文件，可以使用 `param.json` 配置文件批量准备。

### 创建 param.json

```json
{
  "prepare": {
    "strus": ["Si.cif", "Ge.cif", "GaAs.cif"],
    "stru_format": "cif",
    "input_template": "INPUT_template",
    "kpt_template": "KPT_template",
    "pp_dict": {
      "Si": "Si_ONCV_PBE-1.0.upf",
      "Ge": "Ge_ONCV_PBE-1.0.upf",
      "Ga": "Ga_ONCV_PBE-1.0.upf",
      "As": "As_ONCV_PBE-1.0.upf"
    },
    "orb_dict": {
      "Si": "Si_gga_8au_60Ry_2s2p1d.orb",
      "Ge": "Ge_gga_8au_100Ry_2s2p1d.orb",
      "Ga": "Ga_gga_8au_100Ry_2s2p2d.orb",
      "As": "As_gga_7au_100Ry_2s2p1d.orb"
    },
    "pp_path": "/path/to/pseudopotentials",
    "orb_path": "/path/to/orbitals"
  }
}
```

### 准备 INPUT 和 KPT 模板

**INPUT_template：**
```
INPUT_PARAMETERS
suffix              ABACUS
calculation         scf
basis_type          lcao
ecutwfc             100
scf_thr             1e-7
scf_nmax            100
smearing_method     gauss
smearing_sigma      0.01
mixing_type         broyden
mixing_beta         0.4
```

**KPT_template：**
```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

### 执行批量准备

```bash
abacustest prepare -p param.json -s abacustest_jobs
```

生成结果：
```
abacustest_jobs/
├── Si/
│   ├── INPUT
│   ├── STRU
│   ├── KPT
│   └── ...
├── Ge/
│   ├── INPUT
│   ├── STRU
│   ├── KPT
│   └── ...
└── GaAs/
    ├── INPUT
    ├── STRU
    ├── KPT
    └── ...
```

---

## 自动设置 ecutwfc

如果赝势库目录下有 `ecutwfc.json` 文件，abacustest 会自动为每个体系设置合适的 ecutwfc：

**ecutwfc.json 示例：**
```json
{
  "Si": 60,
  "Ge": 80,
  "Ga": 100,
  "As": 100
}
```

abacustest 会自动选择体系中所有元素的最大值作为 ecutwfc。

---

## 常见问题

**Q1：找不到赝势或轨道文件**
- 检查文件名是否以元素名开头
- 或创建 `element.json` 文件指定映射关系

**Q2：生成的 INPUT 参数不符合需求**
- 修改 `INPUT_template` 模板文件
- 或生成后手动修改各任务的 INPUT 文件

**Q3：k 点密度不合适**
- 修改 `KPT_template` 模板文件
- 或使用 `kspacing` 参数自动生成 k 点

---
