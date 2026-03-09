# Step 1 研究笔记

## 1.1 RAG检索结果总结

### 检索质量评估
- 主查询命中率高，案例文件本身已在知识库中（多个版本：官网页面+案例文件）
- 额外获取了 CONTROL 文件完整内容、Phonopy 接口文档
- CONTROL 文件 &param 段完整版已检索到

### 核心知识点

**物理背景**
- 晶格热导率来源：声子散射（声子-声子三声子过程、同位素散射）
- 方法：玻尔兹曼输运方程（BTE） + 力常数（2nd + 3rd）
- 二阶力常数：谐振势，决定声子色散关系
- 三阶力常数：非谐势，决定声子寿命/散射

**计算流程**
```
ABACUS relax → STRU（优化后结构）
     ↓
Phonopy (setting.conf --abacus -d)
     ↓
ABACUS SCF × N个微扰构型（产生受力）
     ↓
Phonopy -f → FORCE_SETS → FORCE_CONSTANTS
     ↓
au2si.py 单位转换 → FORCE_CONSTANTS_2ND (eV/Å²)
     ↓
thirdorder_vasp.py sow → 40个 3RD.POSCAR.*
     ↓
pos2stru.py (ASE) → STRU 文件
     ↓
ABACUS SCF × 40 → 原子受力
     ↓
aba2vasp.py → vasprun.xml
     ↓
thirdorder_vasp.py reap → FORCE_CONSTANTS_3RD
     ↓
ShengBTE (CONTROL + 2ND + 3RD) → κ(T)
```

**案例完整性检查清单**
- [x] STRU文件（结构优化后）- 含所有参数
- [x] setting.conf 内容：DIM=2 2 2，ATOM_NAME=Si
- [x] band.conf 完整内容（7行参数）
- [x] CONTROL 文件完整内容（4个section）
- [x] vasprun.xml 格式示例
- [x] 所有Python脚本：au2si.py, pos2stru.py, aba2vasp.py
- [x] 关键命令：10个命令全部记录
- [x] 关键参数：scf_thr 1e-8（LCAO）/1e-12（PW）
- [x] 结果数值：300K下~100 W/(m·K)，实验值~150 W/(m·K)
- [x] 特殊注意事项：不能用dpdata，必须FULL_FORCE_CONSTANTS=.TRUE.，单位转换

**CONTROL 文件完整内容（来自知识库）**
```
&allocations
    nelements=1
    natoms=2
    ngrid(:)=10 10 10
&end
&crystal
    lfactor=0.100000
    lattvec(:,1)=0 2.81594778072 2.81594778072
    lattvec(:,2)=2.81594778072 0 2.81594778072
    lattvec(:,3)=2.81594778072 2.81594778072 0
    elements="Si"
    types=1 1
    positions(:,1)=0.8750000000000000  0.8750000000000000  0.8750000000000000
    positions(:,2)=0.1250000000000000  0.1250000000000000  0.1250000000000000
    scell(:)=2 2 2
&end
&parameters
    !T=300,
    T_min=200
    T_max=500
    T_step=50
    scalebroad=1.0
&end
&flags
    !espresso=.true.
    nonanalytic=.true.,
    isotopes=.true.
&end
```

## 1.2 案例解析

**材料**：Si 金刚石结构（Fd-3m 空间群），2原子原胞

**文件结构**（案例目录）：
```
LCAO/
├── 2nd/           # 二阶力常数
│   ├── STRU       # relax后的结构
│   ├── setting.conf
│   ├── band.conf
│   └── au2si.py   # 单位转换脚本
├── 3rd/           # 三阶力常数
│   ├── POSCAR     # 由STRU转换
│   ├── pos2stru.py
│   ├── run_stru.sh
│   └── aba2vasp.py
└── shengbte/      # 最终计算
    ├── CONTROL
    ├── FORCE_CONSTANTS_2ND
    ├── FORCE_CONSTANTS_3RD
    └── Ref/       # 参考结果
```

**关键注意事项**（必须在教程中出现）：
1. FULL_FORCE_CONSTANTS = .TRUE. 必须设置（否则ShengBTE报错）
2. 单位转换：FORCE_CONSTANTS (eV/Å·au) → FORCE_CONSTANTS_2ND (eV/Å²)，用 au2si.py
3. thirdorder 只支持 VASP/QE 格式，需要 aba2vasp.py 封装
4. 不能用 dpdata 转换（会旋转晶格导致受力方向错误）
5. scf_thr：LCAO 三阶需要 1e-8，PW 三阶需要 1e-12

## 1.3 风格总结

**参考文章风格特征**：
- 开头：直接说明教程目的（"本教程旨在介绍..."或更直接的方式）
- 结构：一级标题章节（一、二、三...）
- 层级：### 用于子步骤，#### 很少用
- 代码块：包含完整命令，有时包含输出示例
- 参数说明：bullet list，参数名+说明，简洁
- 注意事项：**加粗**或"注意："提醒
- 外部链接：显式列出相关软件的文档链接
- 结语：简短，说明结果，欢迎反馈

**字数控制**：
- 参考文章均值 345行/883词
- 本教程流程复杂（3大步×多子步），目标 700-900行
- 各章长度：引言30-50行，准备50行，二阶力常数200行，三阶力常数200行，ShengBTE 100行，结果讨论80行，总结30行
