# 知识调研结果

## 一、RAG检索结果摘要

### 1.1 弹性常数计算（资料极丰富）

**来源：** `ABACUS+pymatgen 计算弹性常数.md`、`abacus_user_guide__abacus_elastic.html.md`

核心方法：应力-应变法
- 施加6种独立应变（3正应变+3剪切应变），每种取4个幅度（±0.5%、±1%）
- 生成24个形变构型 + 1个原始结构 = 25个计算任务
- 每个形变构型做relax（固定晶胞+允许原子弛豫）
- 线性拟合应力-应变数据得弹性张量

Voigt标记：xx→1, yy→2, zz→3, yz→4, xz→5, xy→6
立方晶系3个独立分量：C11、C12、C44

**Si计算结果（已验证）：**
- C11 = 165.5 GPa（实验值165.7 GPa，偏差<0.2%）
- C12 = 58.6 GPa
- C44 = 82.5 GPa
- BV = 94.19 GPa（体模量）
- GV = 70.89 GPa（剪切模量）
- EV = 170.02 GPa（杨氏模量）
- ν = 0.199（泊松比）

另：pymatgen方式的Si结果（另一套赝势）：
- BV=419.37 GPa, GV=521.42 GPa, EV=1105.91 GPa, ν=0.06
（注意：这是另一个体系，不是Si，是金刚石结构的极端情况）

### 1.2 abacustest工作流（资料丰富）

**来源：** `07_Final_Tutorial_使用abacustest计算晶体弹性性质.md`

完整命令序列：
```bash
# 安装
pip install abacustest

# 准备cell-relax输入
abacustest model inputs -f Si_conv.cif --ftype cif --jtype cell-relax --lcao --folder-syntax Si

# 整理优化后结构
mkdir Si-elastic
cp Si/INPUT Si/Si* Si-elastic
cp Si/INPUT/OUT.ABACUS/STRU_ION_D Si-elastic/STRU

# 生成24+1形变构型
abacustest model elastic prepare -j Si-elastic

# 后处理
abacustest model elastic post -j Si-elastic
```

Bohrium提交方式：使用lbg工具，job.json格式已知

### 1.3 结构优化参数（资料丰富）

**来源：** `ABACUS 使用教程｜结构优化.md`、`abacus_user_guide__abacus_eff2.html.md`

**ionic relax（固定晶胞）：**
```
calculation    relax
cal_force      1
force_thr_ev   0.01   # eV/Å，原子受力收敛阈值
relax_nmax     100    # 最大离子步数
out_stru       1
```

**cell-relax（变胞优化）：**
```
calculation    cell-relax
cal_force      1
cal_stress     1
force_thr_ev   0.01   # eV/Å
stress_thr     0.5    # kBar（BN案例）；Ba案例用2 kBar
relax_nmax     100
out_stru       1
kspacing       0.08   # BN案例使用，代替手写KPT
```

cell-relax的嵌套流程：先固定晶胞做relax → 再更新晶胞参数 → 重复直到收敛
输出日志标记：RELAX CELL : 3 / RELAX IONS : 1 (in total: 15)

### 1.4 BN cell-relax完整案例（来自abacus_user_guide__abacus_eff2.html.md）

体系：h-BN，192原子（B96+N96），层状结构
INPUT关键参数：
```
calculation cell-relax
symmetry 0
basis_type lcao
ecutwfc 100
scf_thr 1e-07
scf_nmax 100
smearing_method gauss
smearing_sigma 0.002
mixing_type pulay
mixing_beta 0.3
cal_force 1
cal_stress 1
force_thr_ev 0.01
stress_thr 0.5
relax_nmax 100
out_stru 1
kspacing 0.08
```

STRU关键部分（赝势/轨道）：
- B_ONCV_PBE-1.0.upf，B_gga_8au_100Ry_2s2p1d.orb
- N_ONCV_PBE-1.0.upf，N_gga_8au_100Ry_2s2p1d.orb
- LATTICE_CONSTANT 1.8897261258369282（Bohr）
- 晶格矢量：a=2.6665 Å方向短轴在X

## 二、风格参考总结

**来自：** `ABACUS 使用教程｜磁性材料计算.md`、`07_Final_Tutorial_使用abacustest计算晶体弹性性质.md`

风格特征：
1. 不用"在当今..."等AI腔开头
2. 用横线（---）分隔大节
3. 参数说明用表格（参数名|含义/默认值/推荐值）
4. 代码块包含文件名注释，参数对齐
5. 注意事项用 > **注意/警告** 格式突出
6. 结果用对比表格展示（本文结果 vs 实验值）
7. 长度适中，不过度展开理论推导
