# Step 1: 知识调研与案例解析

## 一、RAG知识检索结果

### 1.1 DOS/PDOS基本概念

**态密度（DOS, Density of States）**
- 定义：在能量为E的能量附近可供电子占据的电子状态数目
- 表示：通常以每单位能量范围的态数目来表示
- 物理意义：
  - DOS的峰值对应能带结构中的范霍夫奇点
  - 费米能级附近的DOS行为对材料电学性质有显著影响
  - 金属的高DOS导致高电导率
  - 半导体和绝缘体在费米能级附近的低DOS导致低电导率

**投影态密度（PDOS, Partial DOS）**
- 定义：将电子波函数投影到每个原子的每个轨道得到的态密度值
- 作用：可以分析成键信息、导电性等
- 投影方式：按元素、按轨道（s/p/d/f）、按原子

**DOS与能带的关系**
- DOS通过对能带结构在能量上的积分得到
- 能带结构描述电子能量与动量（波矢K）的关系
- 能带结构中的价带和导带之间的能隙决定材料是导体、半导体还是绝缘体

### 1.2 ABACUS计算DOS/PDOS的流程

**标准两步流程：**

**第一步：SCF自洽计算**
- 目的：对结构弛豫后的稳定晶体结构做自洽计算
- 输出：收敛的体系基态电子密度文件
- 关键参数：
  - `calculation = scf`
  - `out_chg = 1`（输出电荷密度）

**第二步：NSCF非自洽计算**
- 目的：读入电荷密度，固定电子密度，计算态密度
- 输出：DOS和PDOS数据文件
- 关键参数：
  - `calculation = nscf`
  - `init_chg = file`（从文件读取电荷密度）
  - `out_dos = 1`（输出总DOS）
  - `out_dos = 2`（同时输出DOS和PDOS）

### 1.3 ABACUS DOS/PDOS相关INPUT参数

**必需参数：**

1. **out_dos**
   - 类型：整数
   - 取值：0（不输出）、1（输出DOS）、2（输出DOS+PDOS）
   - 说明：LCAO基组下可输出PDOS

2. **init_chg**
   - 类型：字符串
   - 取值：`atomic`（SCF）、`file`（NSCF）
   - 说明：NSCF计算必须设置为`file`

3. **calculation**
   - 类型：字符串
   - 取值：`scf`（自洽）、`nscf`（非自洽）
   - 说明：DOS计算需要先SCF后NSCF

**可选参数（影响DOS质量）：**

4. **smearing_method**
   - 类型：字符串
   - 常用值：`gauss`、`mp`（Methfessel-Paxton）、`fd`（Fermi-Dirac）
   - 说明：展宽方法，影响费米面附近的态密度平滑度

5. **smearing_sigma**
   - 类型：浮点数
   - 单位：Rydberg
   - 典型值：0.002-0.015 Ry
   - 说明：展宽参数，用于SCF收敛，不影响NSCF

6. **dos_sigma**
   - 类型：浮点数
   - 单位：eV
   - 说明：仅用于绘制DOS图时的高斯展宽，不影响SCF/NSCF结果

7. **symmetry**
   - SCF时可设为1（利用对称性减少k点）
   - NSCF时建议设为0（不使用对称性，获得完整DOS）

**SDFT方法的DOS参数（高温体系）：**
- `dos_emin_ev`：DOS能量最小范围（eV）
- `dos_emax_ev`：DOS能量最大范围（eV）
- `dos_edelta_ev`：能量间隔（eV）
- `dos_nche`：切比雪夫展开阶数

### 1.4 PDOS文件格式（LCAO基组）

**文件结构：**
- 格式：XML
- 位置：`OUT.suffix/PDOS`

**内容组成：**

1. **头部信息**
   - `<nspin>`：自旋极化（1或2）
   - `<norbitals>`：轨道总数
   - `<energy_values>`：能量范围（横坐标）

2. **轨道信息**
   - `index`：轨道编号
   - `atom_index`：原子编号
   - `species`：元素符号
   - `l`：角量子数（s=0, p=1, d=2, f=3）
   - `m`：磁量子数（0到2l）
   - `z`：径向轨道编号
   - `<data>`：该轨道的DOS数据（两列对应自旋上/下）

**轨道数量示例（Li的DZP基组4s1p）：**
- 4个s轨道：l=0, m=0, z=1,2,3,4
- 3个p轨道：l=1, m=0,1,2, z=1
- 共7个轨道

### 1.5 费米能级提取

**从输出文件提取：**
```bash
grep -i 'efermi' running_scf.log
# 或
grep -i 'efermi' running_nscf.log
```

**输出示例：**
```
EFERMI = -1.4947682735 eV
```

**注意事项：**
- 金属体系：费米能级在导带内
- 半导体/绝缘体：费米能级在带隙中
- 自旋极化体系：自旋上/下共享同一费米能级

---

## 二、案例解析（MgO DOS/PDOS）

### 2.1 案例基本信息
- **材料体系**：MgO（氧化镁，岩盐结构）
- **计算类型**：DOS/PDOS
- **基组类型**：LCAO（局域原子轨道）
- **计算流程**：SCF自洽计算

### 2.2 文件结构
- INPUT：参数设置文件
- STRU：结构文件（MgO岩盐结构）
- KPT：k点文件（18×18×18 Gamma中心网格）
- 赝势文件：
  - Mg.PD04.PBE.UPF
  - O.upf
- 轨道文件：
  - Mg_gga_10au_100Ry_2s1p.orb
  - O_gga_6au_100Ry_2s2p1d.orb

### 2.3 关键INPUT参数

**基本设置：**
- `calculation = scf`：自洽计算
- `basis_type = lcao`：LCAO基组
- `symmetry = 0`：不使用对称性
- `precision = double`：双精度

**平面波截断：**
- `ecutwfc = 100`：100 Ry（与轨道文件推荐值一致）

**电子结构：**
- `ks_solver = genelpa`：Kohn-Sham方程求解器
- `smearing_method = gauss`：高斯展宽
- `smearing_sigma = 0.015`：展宽参数0.015 Ry

**SCF收敛：**
- `mixing_type = broyden`：Broyden混合方法
- `mixing_beta = 0.8`：混合参数
- `scf_nmax = 100`：最大迭代步数
- `scf_thr = 1e-07`：收敛阈值

**DOS输出：**
- `init_chg = file`：从文件读取电荷密度
- `out_dos = 2`：输出总DOS和PDOS

**k点设置：**
- `kspacing = 0.08`：k点间距（单位：1/Bohr）
- KPT文件：18×18×18 Gamma中心网格

### 2.4 STRU文件关键信息

**晶格常数：**
- `LATTICE_CONSTANT = 1.889726`（Bohr转换因子）

**晶格矢量：**
```
2.97691954880  0.00000000000  0.00000000000
1.48845977440  2.57808795428  0.00000000000
1.48845977440  0.85936265143  2.43064463329
```

**原子位置：**
- Mg：1个原子，位于(0, 0, 0)
- O：1个原子，位于(2.977, 1.719, 1.215)

**轨道信息：**
- Mg：2s1p基组（10 au截断半径）
- O：2s2p1d基组（6 au截断半径）

### 2.5 案例完整性检查清单

**必需文件：**
- [x] INPUT文件完整
- [x] STRU文件完整
- [x] KPT文件完整
- [x] 赝势文件已指定
- [x] 轨道文件已指定

**关键参数：**
- [x] `out_dos = 2`（输出PDOS）
- [x] `init_chg = file`（读取电荷密度）
- [x] `basis_type = lcao`（LCAO基组）
- [x] k点网格密度合理（18×18×18）

**计算流程：**
- [x] 需要先运行SCF计算生成电荷密度
- [x] 然后运行NSCF计算输出DOS/PDOS

---

## 三、abacustest DOS/PDOS后处理

### 3.1 abacustest model dos-pdos命令

**命令格式：**
```bash
abacustest model dos-pdos -j <job_path> [options]
```

**主要参数：**
- `-j, --job`：ABACUS计算目录路径（必需）
- `--range RANGE RANGE`：能量范围（相对费米能级），默认-10到10 eV
- `--plot-type`：绘图类型
  - `species`：按元素（如C）
  - `shell`：按轨道壳层（如C的p轨道）
  - `orbital`：按具体轨道（如C的p_x轨道）
  - `atom`：按原子编号
- `--atom-index`：原子编号（1-based），仅在plot-type=atom时有效
- `--suffix`：输出文件后缀
- `--no-save-data`：不保存数据文件
- `--no-save-plot`：不保存图片文件

**输出文件：**
- `DOS.dat`：总态密度数据
- `DOS.png`：总态密度图
- `PDOS.dat`：投影态密度数据
- `PDOS.png`：投影态密度图

### 3.2 使用示例

**基本用法（默认参数）：**
```bash
cd <计算目录>
abacustest model dos-pdos -j .
```

**指定能量范围：**
```bash
abacustest model dos-pdos -j . --range -5 7
```

**按元素绘制PDOS：**
```bash
abacustest model dos-pdos -j . --plot-type species
```

**按轨道壳层绘制：**
```bash
abacustest model dos-pdos -j . --plot-type shell
```

### 3.3 与abacus-plot的对比

**abacus-plot（传统方法）：**
- 需要手动配置config.json
- 需要手动提取费米能级
- 命令：`abacus-plot -d -p -o`

**abacustest（推荐方法）：**
- 自动识别输出文件
- 自动提取费米能级
- 一条命令完成所有操作
- 更灵活的绘图选项

---

## 四、风格参考总结

### 4.1 参考文章分析

**参考文章：** ABACUS 使用教程｜磁性材料计算.md

**风格特征：**

1. **开头方式**
   - 直接进入主题，先介绍物理概念
   - 用简洁的语言解释微观机制
   - 不使用"在当今..."、"随着...的发展"等AI腔

2. **结构特点**
   - 清晰的层次结构（介绍→准备→案例→注意事项）
   - 使用表格对比不同情况的参数设置
   - 代码块包含完整的文件内容和注释

3. **语言特征**
   - 句式简洁，单句不超过30字
   - 使用专业术语但配有解释
   - 避免冗长的背景介绍
   - 使用"该"、"其"等指代词保持简洁

4. **代码呈现**
   - 代码块标注文件类型（如```INPUT）
   - 包含必要的注释说明参数含义
   - 参数对齐，便于阅读
   - 给出完整可运行的示例

5. **参数说明方式**
   - 使用表格列出参数及其取值
   - 每个参数配有简短说明
   - 给出典型值和推荐值
   - 说明参数之间的关系

### 4.2 本教程应遵循的风格

**开头：**
- 直接说明DOS/PDOS的物理意义
- 简要说明计算流程（SCF+NSCF）
- 不展开DFT理论背景（前面教程已讲）

**正文：**
- 参数说明使用表格
- 代码块完整且带注释
- 每节控制在150-300行
- 避免与前3篇教程重复（INPUT/STRU/KPT基础格式）

**案例：**
- 完整呈现MgO案例的所有文件
- 说明每个参数的作用
- 给出预期输出结果
- 展示abacustest后处理命令

**语言：**
- 避免AI腔（"值得注意的是"、"综上所述"等）
- 句式简洁，段落适中
- 使用主动语态
- 技术准确，不模糊表述

---

## 五、检索质量评估

### 5.1 检索结果覆盖度

**已覆盖的知识点：**
- ✅ DOS/PDOS基本概念和物理意义
- ✅ ABACUS计算DOS/PDOS的两步流程
- ✅ INPUT参数详细说明（out_dos, init_chg等）
- ✅ PDOS文件格式（XML结构）
- ✅ 费米能级提取方法
- ✅ abacustest dos-pdos命令用法
- ✅ MgO案例完整信息

**检索质量：**
- 检索结果相关性高，来源权威（官方文档、教程）
- 参数说明详细，有示例代码
- 覆盖了LCAO基组的DOS/PDOS计算
- 包含了后处理工具的使用方法

**需要补充的内容：**
- 平面波基组的DOS计算（可简要提及）
- DOS/PDOS图的解读方法
- 常见问题和注意事项

### 5.2 案例完整性评估

**MgO案例：**
- ✅ 文件结构完整（INPUT/STRU/KPT/赝势/轨道）
- ✅ 参数设置合理（LCAO基组，18×18×18 k点）
- ✅ 计算流程清晰（SCF→NSCF）
- ✅ 适合作为教程案例（结构简单，物理意义明确）

**案例优势：**
- MgO是典型的离子晶体，PDOS可清晰展示O 2p和Mg 3s轨道
- 结构简单（2原子），便于理解
- 计算量适中，适合教学演示

---

## Think Aloud - Step 1总结

**检索执行情况：**
1. 执行了3次RAG检索，获取了10+8+5=23个相关文档
2. 解析了MgO案例文件，提取了完整的参数信息
3. 阅读了风格参考文章，提取了写作风格特征

**检索结果质量：**
- DOS/PDOS的理论知识充足，来源权威
- INPUT参数说明详细，有多个示例
- abacustest工具使用方法清晰
- MgO案例完整可用

**发现的关键信息：**
1. DOS计算必须是两步流程（SCF+NSCF）
2. LCAO基组才能输出PDOS（out_dos=2）
3. NSCF时必须设置init_chg=file
4. abacustest可以一键完成后处理
5. 费米能级需要从log文件提取

**下一步计划：**
- 进入Step 2，设计3个大纲方案
- 方案需要考虑：避免与前3篇重复、突出DOS/PDOS特色、详细程度适中
