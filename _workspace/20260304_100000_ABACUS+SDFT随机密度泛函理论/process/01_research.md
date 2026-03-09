# Step 1 调研结果

## 1.1 RAG检索摘要

### 核心原理（来自检索文档）

**SDFT vs KSDFT 的优势：**
- KSDFT在高温（数十~上千eV）下需要海量占据态波函数→对角化成本O(N³)爆炸
- SDFT用随机波函数替代对角化哈密顿量，计算时间与温度成反比（温度越高越快）
- MDFT（混合）：加入少量低能KS轨道，大幅加速随机误差收敛
- 参考文献：Qianrui Liu & Mohan Chen, Phys. Rev. B 106, 125132 (2022)

**SDFT并行效率（检索文档）：**
- 随机轨道操作独立，并行效率接近线性
- MDFT比SDFT稍慢（有部分KS轨道对角化）
- 已成功用于512碳原子体系，16~256 CPU线性扩展

### 核心参数表

| 参数 | 作用 | 默认值 | 备注 |
|------|------|--------|------|
| esolver_type | 求解器类型 | ksdft | 设为 sdft 启用SDFT/MDFT |
| nbands | KS确定性轨道数 | - | 0：纯SDFT；>0：MDFT |
| nbands_sto | 随机轨道数 | - | 越大误差越小，效率越低 |
| nche_sto | 切比雪夫展开阶数 | - | 与温度反比，与ecut正比；目标：Chebyshev Precision<1e-8 |
| method_sto | 计算方法 | 2 | 1：低内存慢；2：快高内存 |
| seed_sto | 随机种子 | 0（随时间） | 设大于1的整数可复现 |
| bndpar | 进程分组数 | 1 | 优先kpar并行，再测bndpar |
| smearing_method | 展宽类型 | - | SDFT只支持 fd |
| smearing_sigma | 电子温度 | - | 单位 Ry（1 Ry≈13.6 eV） |

### nbands_sto 收敛测试方法
- 用10个不同随机种子分别跑SDFT，取平均能量和误差
- 增加nbands_sto直到能量误差 < 万分之一

### ecut 测试方法
- 固定nbands_sto，对不同ecut（每次加10 Ry）跑10个随机种子
- 相邻两个ecut的平均能量差 < 万分之一

### DOS 专属参数
| 参数 | 作用 | 案例值 |
|------|------|--------|
| out_dos | 输出DOS开关 | 1 |
| dos_emin_ev | 能量下限(eV) | -20 |
| dos_emax_ev | 能量上限(eV) | 100 |
| dos_edelta_ev | 能量间隔(eV) | 0.1 |
| dos_sigma | 高斯展宽(eV) | 4 |
| dos_nche | DOS切比雪夫阶数 | 240 |
| npart_sto | 内存控制 | 2（减半内存） |
| emax_sto/emin_sto | DOS能量范围（Ry，0=自动） | 0 |

输出文件：`OUT.autotest/DOS1_smearing.dat`

### MD 专属参数
| 参数 | 作用 | 案例值 |
|------|------|--------|
| md_tfirst | 初始温度(K) | 1160400（≈100 eV） |
| md_dt | 时间步长(fs) | 0.2 |
| md_nstep | MD步数 | 10 |

---

## 1.2 案例解析

### 三个算例总览

| 案例 | 体系 | 电子温度 | 计算类型 | nbands/nbands_sto |
|------|------|----------|----------|-------------------|
| pw_Si2 | 2原子Si（金刚石）| 0.6 Ry（≈8.16 eV）| SCF（MDFT） | 4/64 |
| pw_md_Al | 16原子Al | 7.35 Ry（≈100 eV）| MD（纯SDFT） | 0/64 |
| 186_PW_SDOS_10D10S | 1原子Si | 0.6 Ry（≈8.16 eV）| SCF+DOS（MDFT） | 10/10 |

### 案例1：pw_Si2 INPUT
```
INPUT_PARAMETERS
#Parameters (General)
calculation    scf
esolver_type   sdft
pseudo_dir     ../../PP_ORB
nbands         4
nbands_sto     64
nche_sto       100
method_sto     1
#Parameters (Accuracy)
ecutwfc        50
scf_nmax       20
symmetry       1
#Parameters (Smearing)
smearing_method fd
smearing_sigma  0.6
```

### 案例2：pw_md_Al INPUT
```
INPUT_PARAMETERS
#Parameters (General)
calculation    md
esolver_type   sdft
pseudo_dir     ../../PP_ORB
nbands         0
nbands_sto     64
nche_sto       20
method_sto     2
#Parameters (Accuracy)
ecutwfc        50
scf_nmax       20
scf_thr        1e-6
symmetry       1
#Parameters (Smearing)
smearing_method fd
smearing_sigma  7.34986072
#Parameters (MD)
md_tfirst      1160400
md_dt          0.2
md_nstep       10
```

注意：
- nche_sto=20（温度100 eV，阶数可更小）
- nbands=0（纯SDFT，无KS轨道）
- md_tfirst=1160400 K = 100 eV

### 案例3：186_PW_SDOS_10D10S INPUT
```
INPUT_PARAMETERS
#Parameters (1.General)
suffix          autotest
calculation     scf
esolver_type    sdft
method_sto      2
nbands          10
nbands_sto      10
nche_sto        120
emax_sto        0
emin_sto        0
seed_sto        20000
pseudo_dir      ../../PP_ORB
symmetry        1
kpar            1
bndpar          2
#Parameters (2.Iteration)
ecutwfc         20
scf_thr         1e-6
scf_nmax        20
#Parameters (3.Basis)
basis_type      pw
#Parameters (4.Smearing)
smearing_method fd
smearing_sigma  0.6
#Parameters (5.Mixing)
mixing_type     broyden
mixing_beta     0.4
out_dos         1
dos_emin_ev     -20
dos_emax_ev     100
dos_edelta_ev   0.1
dos_sigma       4
dos_nche        240
npart_sto       2
```

---

## 1.3 风格参考总结

- 参考文章（磁性材料）：先介绍概念，再参数表格，再流程代码块
- 开头不用"在当今..."，直接进入主题
- 参数说明清晰：参数名+说明+默认值
- 代码块对齐，有注释
- 结构：介绍→准备→案例（流程）→结语
- 平均长度：约350行
