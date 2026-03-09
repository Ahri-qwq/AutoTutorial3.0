# Step 2 大纲方案

## 背景信息
- 参考文章平均：345行、883 words
- 目标教程长度：600-800行（约1.7-2.3倍参考）
- 三个案例：SCF（Si 2原子）、MD（Al 16原子）、DOS（Si 1原子）

---

## 方案A：参数手册型（约500行）

**核心思路：** 把参数讲解放在核心位置，三个案例作为参数说明的具体实例。适合已有SDFT背景、需要快速查阅参数的用户。

```
# ABACUS 随机密度泛函理论（SDFT）使用指南

## 一、SDFT 与 MDFT 简介（约30行）
- KSDFT的高温瓶颈
- SDFT的核心思路
- MDFT的进一步加速

## 二、算例下载与文件结构（约30行）
- 下载命令
- 三个算例目录说明

## 三、关键参数详解（约130行）
### 3.1 基本控制参数（esolver_type/nbands/nbands_sto）
### 3.2 精度参数（nche_sto/method_sto）
### 3.3 并行参数（kpar/bndpar）
### 3.4 DOS相关参数（out_dos/dos_**/dos_nche/npart_sto）
### 3.5 MD相关参数（md_tfirst/md_dt/md_nstep）

## 四、案例一：SCF 计算（Si，0.6 Ry，MDFT）（约80行）
- INPUT完整文件
- 运行方法
- 注意事项

## 五、案例二：分子动力学模拟（Al，100 eV，纯SDFT）（约80行）
- INPUT完整文件
- 运行方法
- 纯SDFT与MDFT的区别

## 六、案例三：态密度计算（Si，0.6 Ry，MDFT+DOS）（约90行）
- INPUT完整文件
- 输出文件说明
- npart_sto内存控制

## 七、精度收敛测试（约60行）
- nbands_sto收敛测试方法
- nche_sto选取方法
- ecut测试注意事项

## 附录（约20行）
```

**优点：** 参数组织清晰，便于查阅。
**缺点：** 参数与案例割裂，读者需来回参考；案例部分较薄。

---

## 方案B：三案例并行型（约650行）[推荐]

**核心思路：** 以三个具体算例为主轴，每个案例内部讲清背景→INPUT→参数→运行。参数在首次出现时详细介绍，后续案例中复用的参数只简要提及。结构完整，适合实际操作。

```
# ABACUS 随机密度泛函理论（SDFT）使用教程

## 一、引言（约40行）
- KSDFT在高温的局限
- SDFT的思路与优势（计算时间随温度反比）
- MDFT的混合策略
- 本教程三个案例的适用场景预览

## 二、准备工作（约30行）
- 软件版本要求（ABACUS 3.2.0+）
- 算例下载
- 赝势说明（高温时赝势可移植性注意）

## 三、案例一：SCF 计算（Si，MDFT，0.6 Ry）（约150行）
### 3.1 场景说明（2原子金刚石Si，T=8.16 eV）
### 3.2 INPUT 文件与参数详解
  - esolver_type / nbands(4) / nbands_sto(64) / nche_sto(100) / method_sto(1)
  - smearing_method=fd 的必要性
### 3.3 运行方法
### 3.4 注意事项（赝势高温可移植性、K点收敛）

## 四、案例二：分子动力学模拟（Al，纯SDFT，100 eV）（约130行）
### 4.1 场景说明（16原子Al，T=100 eV，纯SDFT）
### 4.2 INPUT 文件与参数详解
  - calculation=md / nbands=0（纯SDFT）/ nche_sto=20（高温阶数小）
  - md_tfirst=1160400 / md_dt=0.2 / md_nstep=10
### 4.3 纯SDFT vs MDFT的选择

## 五、案例三：态密度计算（Si，MDFT+DOS，0.6 Ry）（约160行）
### 5.1 场景说明（1原子Si，MDFT，T=8.16 eV）
### 5.2 INPUT 文件与参数详解
  - nbands=10 / nbands_sto=10 / bndpar=2 / seed_sto=20000
  - out_dos=1 / dos_emin_ev / dos_emax_ev / dos_edelta_ev
  - dos_sigma / dos_nche / npart_sto=2
### 5.3 输出文件（DOS1_smearing.dat）

## 六、精度控制与参数选取（约90行）
### 6.1 nbands_sto 收敛测试
### 6.2 nche_sto 选取（目标：Chebyshev Precision < 1e-8）
### 6.3 ecut 测试方法
### 6.4 并行策略（kpar优先，bndpar补充）

## 附录（约30行）
- 参考文献
- 常见问题
```

**优点：** 案例驱动，参数在上下文中学习；三个案例覆盖主要使用场景；结构清晰。
**缺点：** 参数说明稍有分散，但首次详述后续引用可解决。

---

## 方案C：渐进进阶型（约700行）

**核心思路：** 从SDFT到MDFT，从SCF到MD，从总能到DOS，展示一条学习曲线。强调方法对比（SDFT vs MDFT vs KSDFT），理论介绍比方案B略多。

```
# ABACUS 随机密度泛函理论：从入门到实战

## 一、为什么需要 SDFT（约60行）
- KSDFT高温困境（能带数爆炸、O(N³)对角化）
- SDFT的随机轨道思路（无需对角化）
- MDFT的加速原理
- 适用场景速查

## 二、快速开始：SDFT-SCF 计算（Si，0.6 Ry）（约100行）
- 最小案例（pw_Si2）
- INPUT讲解：esolver_type / nbands / nbands_sto / nche_sto
- smearing_method=fd 的说明
- 运行与输出

## 三、加入 KS 轨道：MDFT 的优势（约80行）
- nbands>0 的效果（误差随KS轨道数快速下降）
- nbands=4的含义（低能KS轨道）
- 与方案二（纯SDFT）的对比

## 四、高温大体系：MD 模拟（Al，100 eV）（约130行）
- 纯SDFT（nbands=0）的适用时机
- nche_sto与温度的关系（100 eV只需阶数20）
- MD参数（md_tfirst/md_dt/md_nstep）
- method_sto=2的内存考量

## 五、态密度计算（Si，MDFT+DOS）（约160行）
- emax_sto / emin_sto / dos_nche / npart_sto
- 完整INPUT讲解
- 输出文件说明

## 六、参数选取总结（约100行）
- 参数互相关系汇总表
- 收敛测试流程图（nbands_sto→nche_sto→ecut）
- 并行策略

## 附录（约30行）
```

**优点：** 学习路线清晰，对初学者友好；方法对比帮助理解本质。
**缺点：** 篇幅最长，理论内容稍多；对熟悉SDFT的用户冗余。
