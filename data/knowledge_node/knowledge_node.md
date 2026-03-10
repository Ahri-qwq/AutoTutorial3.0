# ABACUS 培训大纲（基础 + 进阶）

http://abacus.ustc.edu.cn/\_upload/tpl/0c/d8/3288/template3288/download/ABACUS\_Developer\_Guide.pdf（《ABACUS 软件使用与开发手册》V0.1 版本下载地址）

https://mcresearch.github.io/abacus-user-guide/ （Github上的中文手册地址）

https://abacus.deepmodeling.com/en/latest/index.html（Deepmodeling.ABACUS Documentation）

https://abacus.ustc.edu.cn/main.html（ABACUS官方网站）

https://kirk0830.github.io/ABACUS-Pseudopot-Nao-Square/about.html (APNS）

## **一、基础篇：安装与基本计算**

目标：让学员能独立安装 ABACUS 并完成最基本的计算任务。

### **1. 软件与环境**

* 架构介绍：ABACUS 的发展历史与定位，及特点

https://arxiv.org/abs/2501.08697（ABACUS 文章）

https://mp.weixin.qq.com/s/Fx5jx99g9PVjgnLAsOyQ9Q(ABACUS ：携手DeepModeling，做源自中国、开源开放的DFT软件)

https://mp.weixin.qq.com/s/V2WVr\_O2\_LJ-VJ09MaHaqQ (ABACUS：一款开源开放的AI4S时代电子结构软件包)

* 安装与编译
  * 官方完整教程

https://github.com/deepmodeling/abacus-develop（ABACUS 在 DeepModeling 社区中的 GitHub 仓库地址）

https://gitee.com/deepmodeling/abacus-develop （ABACUS 的 Gitee 镜像仓库地址）

https://abacus.deepmodeling.com/en/latest/quick\_start/easy\_install.html （ABACUS官方编译教程 Easy Installation）

https://abacus.deepmodeling.com/en/latest/advanced/install.html（ABACUS官方编译教程 Advanced Installation Options）

* 依赖库（MPI、LibXC、ELPA 等）

https://mcresearch.github.io/abacus-user-guide/abacus-hpc.html#21-abacus-%E4%BE%9D%E8%B5%96%E7%9A%84%E8%BD%AF%E4%BB%B6%E5%BA%93（ABACUS 依赖的软件库）

* CPU/GPU/DCU 编译方法

CPU

https://mcresearch.github.io/abacus-user-guide/abacus-gcc.html（GCC 编译 ABACUS 教程）

https://mcresearch.github.io/abacus-user-guide/abacus-oneapi.html（Intel oneAPI 2024/2025 编译 ABACUS 教程）

https://mcresearch.github.io/abacus-user-guide/abacus-intel.html（Intel oneAPI 编译 ABACUS 教程）

GPU

https://mcresearch.github.io/abacus-user-guide/abacus-gpu.html（编译 Nvidia GPU 版本的 ABACUS）

https://mcresearch.github.io/abacus-user-guide/abacus-gpu-lcao.html（ABACUS LCAO 基组 GPU 版本使用说明）

DCU

https://mcresearch.github.io/abacus-user-guide/abacus-dcu.html（ABACUS 在曙光 DCU 集群与 AMD 显卡上的编译与使用）

* Toolchain安装教程

https://www.bohrium.com/notebooks/5215742477 （一键配置编译ABACUS | toolchain 脚本的使用）

https://mp.weixin.qq.com/s/ypc0RT5ePm0vMlGRYI46JA（ABACUS安装教程 - Toolchain (1-GNU））

https://mp.weixin.qq.com/s/K58BVQwSoxcgcNufODo8tg（ABACUS安装教程 - Toolchain (2-Intel））

https://mp.weixin.qq.com/s/S0EfqhNn-FLh\_QXRZE2bqw (ABACUS安装教程 - Toolchain (3-AMD）)

https://mp.weixin.qq.com/s/aym7LBcgfutDtlKOR0aevA （ABACUS安装教程 - Toolchain (4) - 离线安装）

https://mp.weixin.qq.com/s/YekI3LXy7vuTN36ut5M6Rg （ABACUS安装教程 - Toolchain (5) - GPU (1)）

* 其他安装教程

[编译无MPI的ABACUS](https://xmywuqhxb0.feishu.cn/docx/JCv0dHPP6o69JdxtG33cIcTnnke)（编译无MPI的ABACUS）

https://zhuanlan.zhihu.com/p/574031713 （ABACUS 3.0 安装笔记 ）

* 常见报错与解决方法

https://abacus.deepmodeling.com/en/latest/community/faq.html (Deepmodeling.ABACUS Documentation: Frequently Asked Questions)

### **2. 输入输出体系**

* 主要输入文件：`INPUT`, `STRU`, `KPT`

https://abacus.deepmodeling.com/en/latest/quick\_start/input.html （Deepmodeling.ABACUS Documentation: Brief Introduction of the Input Files）

https://abacus.deepmodeling.com/en/latest/advanced/input\_files/index.html（Detailed Introduction of the Input Files）

https://mp.weixin.qq.com/s/FysreUHIRB5RHKDtoj99vQ （ABACUS入门教程 - 结构文件STRU）

https://www.bohrium.com/notebooks/9814968648（ABACUS 使用教程｜如何转换 STRU 文件）

* 运行输出解读：能量收敛、电子步收敛、力与应力

https://abacus.deepmodeling.com/en/latest/quick\_start/output.html（Deepmodeling.ABACUS Documentation: Brief Introduction of the Output Files）

* 参数设置经验
  
  - ​**k 点**​：Monkhorst–Pack 网格、收敛性测试

https://mcresearch.github.io/abacus-user-guide/abacus-pw.html（ABACUS 的平面波计算与收敛性测试）

https://mp.weixin.qq.com/s/auO0I5gjhPP-lxKY6evhRQ（PYATB 能带模块介绍）

- ​**平面波截断能**​：能量收敛、不同赝势的差异

https://mcresearch.github.io/abacus-user-guide/abacus-pw.html（ABACUS 的平面波计算与收敛性测试）

- ​**赝势选择**​：ONCV vs HGH，元素库差异

https://mp.weixin.qq.com/s/vUG-7uRjbMZJP2fs8Ep\_8Q (ABACUS赝势-轨道捆绑包APNS-PPORB-v1发布：精度与效率平衡的初探)

https://www.bohrium.com/notebooks/6416644691（ABACUS 无轨道密度泛函理论方法使用教程）

https://mcresearch.github.io/abacus-user-guide/abacus-ofdft.html（ABACUS 无轨道密度泛函理论方法使用教程）

- **基组选择**  （plan wave vs. NAO）
  
  - SZ / DZP / TZDP 等 NAO 基组的区别    - 如何做基组收敛测试    - 基组与赝势匹配原则    - 下载地址

http://abacus.ustc.edu.cn/pseudo/list.htm （下载地址Pseudopotental files, and corresponding optimized atomic basis sets for ABACUS）

https://mp.weixin.qq.com/s/vUG-7uRjbMZJP2fs8Ep\_8Q (ABACUS赝势-轨道捆绑包APNS-PPORB-v1发布：精度与效率平衡的初探)

https://www.bohrium.com/notebooks/9319634192 ｜https://mcresearch.github.io/abacus-user-guide/abacus-nac1.html（数值原子轨道（一）：ABACUS 中的数值原子轨道命名和使用方法）

https://www.bohrium.com/notebooks/5215642163｜https://mcresearch.github.io/abacus-user-guide/abacus-nac2.html（数值原子轨道（二）：生成给定模守恒赝势的数值原子轨道）

https://www.bohrium.com/notebooks/8841868194 ｜https://mcresearch.github.io/abacus-user-guide/abacus-nac3.html（数值原子轨道（三）：产生高精度数值原子轨道）

https://mcresearch.github.io/abacus-user-guide/abacus-upf.html（模守恒赝势生成方法简介）

https://www.bohrium.com/notebooks/7417640496 （ABACUS 使用教程｜电子自洽迭代（LCAO 基组与 PW 基组）

#### 参数介绍与设置技巧

SCF收敛算法相关参数：[ABACUS收敛性问题解决手册](https://ucoyxk075n.feishu.cn/docx/R0sqdk6T0o2RY4x5IWgcJ3RHnug?from=from_copylink)

### **3. 基础计算实践**

- **结构优化**

* 固定晶胞优化：例子 ...
* 变胞优化：例子...

https://www.bohrium.com/notebooks/9119461238（ABACUS 使用教程｜结构优化/晶格弛豫/几何优化 ）

https://www.bohrium.com/notebooks/4641406377（快速开始 ABACUS｜自洽 能带 态密度 结构优化）

https://abacus.deepmodeling.com/en/latest/advanced/opt.html(Deepmodeling.ABACUS Documentation:  Geometry Optimization)

- **电子结构**

https://abacus.deepmodeling.com/en/latest/advanced/scf/index.html（Running SCF）

* 自洽计算 vs. 非自洽计算

https://mcresearch.github.io/abacus-user-guide/algorithm-mix.html（电荷密度混合算法介绍）

[如何正确画能带，NSCF读电荷密度](https://xmywuqhxb0.feishu.cn/docx/K8GRdTst4oXQNoxnQVbcFZTmntb)（如何正确画能带，NSCF读电荷密度）

* 能带结构计算

https://abacus.deepmodeling.com/en/latest/advanced/elec\_properties/band.html（Extracting Band Structure）

https://www.bohrium.com/notebooks/4641406377（快速开始 ABACUS｜自洽 能带 态密度 结构优化）

https://mp.weixin.qq.com/s/auO0I5gjhPP-lxKY6evhRQ（PYATB 能带模块介绍）

* 态密度 (DOS) 分析

https://abacus.deepmodeling.com/en/latest/advanced/elec\_properties/dos.html（Calculating DOS and PDOS）

[如何正确画能带，NSCF读电荷密度](https://xmywuqhxb0.feishu.cn/docx/K8GRdTst4oXQNoxnQVbcFZTmntb)（如何正确画能带，NSCF读电荷密度）

https://www.bohrium.com/notebooks/1211642609 （用ABACUS-ASE自动产生能带路径 ）

[ABACUS里怎样做DOS和PDOS计算](https://xmywuqhxb0.feishu.cn/docx/ONSldj82VoNGKSxaoDQcoKBtnGh)（ABACUS里怎样做DOS和PDOS计算）

https://mcresearch.github.io/abacus-user-guide/abacus-pdos.html（ABACUS 计算 PDOS）

https://abacus.deepmodeling.com/en/latest/advanced/scf/advanced.html（SCF in Complex Environments）

- ​**分子动力学 (MD**​**) 计算**

* AIMD 的基本概念
* NVT 系综下的简单模拟（如 Si 晶体升温过程）
* 输出文件与轨迹分析

https://www.bohrium.com/notebooks/2241262724（ABACUS 使用教程｜分子动力学）

http://mcresearch.github.io/abacus-user-guide/abacus-md.html（ABACUS 分子动力学使用教程）

https://abacus.deepmodeling.com/en/latest/advanced/md.html（Deepmodeling.ABACUS Documentation:  Molecular Dynamics）

https://mcresearch.github.io/abacus-user-guide/abacus-candela.html（ABACUS+Candela 使用教程）

- **电荷密度**​**与**​**波函数**

* 电荷密度分布

https://abacus.deepmodeling.com/en/latest/advanced/elec\_properties/Mulliken.html（Mulliken Charge Analysis）

https://abacus.deepmodeling.com/en/latest/advanced/elec\_properties/potential.html（Extracting Electrostatic Potential）

https://abacus.deepmodeling.com/en/latest/advanced/elec\_properties/charge.html（Extracting Charge Density）

https://mcresearch.github.io/abacus-user-guide/abacus-chg.html（ABACUS 输出部分的电荷密度和波函数及可视化教程）

* 波函数可视化与后处理

https://abacus.deepmodeling.com/en/latest/advanced/elec\_properties/wfc.html（Extracting Wave Functions）

https://mcresearch.github.io/abacus-user-guide/abacus-elf.html（ABACUS 计算电子局域函数 ELF 使用教程）

**- 基于Bohrium平台的运用 ​**

https://www.bohrium.com/notebooks/6643676733?utm\_source=help\_doc  （Bohrium 帮助文档｜ABACUS）

https://mcresearch.github.io/abacus-user-guide/abacus-dpa3-toturial.html（ABACUS+DPA-3：基于 Bohrium 平台的安装与使用）

## **二、进阶篇：高级功能与物性计算**

目标：掌握 ABACUS 的核心功能，能独立进行科研层面的物性研究。

### **1. 晶体材料**

* 对称性分析：周期性与对称性、点阵、惯用胞与原胞、基矢、晶向与晶面、对称操作、点群与空间群以及晶体内部的对称要素
* 结合能：基于DFT的能量最低构型搜索
* 带隙计算：基础PBE带隙计算流程，几种带隙修正算法：DFT+U，杂化泛函，DFT+1/2

https://mcresearch.github.io/abacus-user-guide/abacus-exx.html（ABACUS 平面波基组下的杂化泛函）

https://www.bohrium.com/notebooks/2112617648（ABACUS 使用教程｜DFT+U 计算）

https://mp.weixin.qq.com/s/auO0I5gjhPP-lxKY6evhRQ（PYATB 能带模块介绍）

* 拓扑计算：材料拓扑性质基础知识

https://mp.weixin.qq.com/s/WQaVzaXZ86ShjAtCcFmWjQ（ABACUS应用案例分享：薄膜技术显著增强TaAs中的非线性霍尔效应）

* ABACUS+PYATB计算流程

https://mp.weixin.qq.com/s/auO0I5gjhPP-lxKY6evhRQ（PYATB 能带模块介绍）

https://www.bohrium.com/notebooks/35791425971（ABACUS+pyatb：介电函数与线性光学的性质的计算）

* ABACUS+wannier90+wannier-tools

https://mcresearch.github.io/abacus-user-guide/abacus-wannier.html（ABACUS+Wannier90 使用教程）

https://mcresearch.github.io/abacus-user-guide/algorithm-wannier.html（最大局域化 Wannier 函数方法简介）

* 磁性计算

https://mp.weixin.qq.com/s/64XicP\_bionUZOYMeXvc-g （ABACUS计算磁性相互作用参数）

https://mp.weixin.qq.com/s/8ZO682BIJy7w0ge-qJB\_9g （ABACUS磁性交换作用的计算以及与微磁学的结合）

https://www.bohrium.com/notebooks/7141761751 (ABACUS 使用教程｜磁性材料计算)

* 电子自旋-原子磁矩-晶体宏观磁性（铁磁-反铁磁-错磁-顺磁），
* 自旋极化，
* 非共线自旋，
* SOC效应，
* 磁交换参数与各向异性能量
* 铁电极化：基础概念，Berry 相位方法

https://abacus.deepmodeling.com/en/latest/advanced/elec\_properties/Berry\_phase.html（Berry Phase Calculation）

https://mp.weixin.qq.com/s/\_56PIb94LVDXOxKvqSW-lg（ABACUS+pyatb Berry Curvature Dipole计算）

* 光学性质分析：光吸收谱，ABACUS+PYATB计算流程，ABACUS+wannier90+wannier-tools，GW方法

https://mp.weixin.qq.com/s/auO0I5gjhPP-lxKY6evhRQ（PYATB 能带模块介绍）

* 实时含时密度泛函 (rt-TDDFT)
  * 基本原理：速度规范 vs 长度规范、混合规范
  * 外场激发：泵浦–探测、电流与极化响应
  * 应用示例：光学吸收谱、超快动力学 ...

https://mp.weixin.qq.com/s/uMg2jMgDSzrP1VYqLhobog（混合规范RT-TDDFT——兼具效率和精度的算法创新）

https://mcresearch.github.io/abacus-user-guide/abacus-tddft.html（ABACUS 实时演化含时密度泛函理论使用教程）

* 振动分析：声子谱与热力学性质

https://www.bohrium.com/notebooks/8741867512｜https://mcresearch.github.io/abacus-user-guide/abacus-phonopy.html （ABACUS+Phonopy 计算声子谱）

* 分子动力学 (AIMD)：NVT、NPT 系综，复杂体系模拟
* 机器学习势函数：DPGEN，DeePMD，模型微调和蒸馏，ABACUS读入DP模型的分子动力学模拟

https://mp.weixin.qq.com/s/2YGVyNXMGL6S\_MAyNhMADw （探索 AI+DFT 最前线，「原子算筹」ABACUS 3.0 重磅发布！）

https://www.bohrium.com/notebooks/26592976824（快速上手使用ABACUS + DP-GEN | init\_bulk/run）

https://mcresearch.github.io/abacus-user-guide/abacus-dpgen.html（ABACUS+DPGEN 使用教程）

https://mcresearch.github.io/abacus-user-guide/abacus-deepks-es.html（DeePKS-ES 介绍及使用教程）

https://mcresearch.github.io/abacus-user-guide/abacus-dpmd.html（ABACUS+DeePMD-kit 做机器学习分子动力学）

### **2. 表面与界面材料**

* 建模：真空层的影响
* 吸附能与催化：NEB，ATST-tools使用方法

https://www.bohrium.com/notebooks/39369325971 （ASE-ABACUS | 第二章：NEB过渡态计算与ATST-Tools）

* 层错能
* 功函数

https://mcresearch.github.io/abacus-user-guide/abacus-surface1.html（采用 ABACUS 进行表面计算（一）：静电势和功函数）

https://mcresearch.github.io/abacus-user-guide/abacus-surface2.html（采用 ABACUS 进行表面计算（二）：偶极修正）

https://mcresearch.github.io/abacus-user-guide/abacus-surface3.html（采用 ABACUS 进行表面计算（三）：表面能计算）

https://mcresearch.github.io/abacus-user-guide/abacus-surface4.html（采用 ABACUS 进行表面计算（四）：表面缺陷能和吸附能计算）

https://mcresearch.github.io/abacus-user-guide/abacus-surface5.html（采用 ABACUS 进行表面计算（五）：外加电场）

https://mcresearch.github.io/abacus-user-guide/abacus-surface6.html（采用 ABACUS 进行表面计算（六）：补偿电荷）

https://mcresearch.github.io/abacus-user-guide/abacus-sol.html（ABACUS 隐式溶剂模型使用教程）

### **3. 大体系材料计算**

* 缺陷与掺杂：
  * 点缺陷建模：超胞方法
  * 形成能与电荷态转变
  * 宽禁带半导体缺陷调控
* ABACUS大体系计算技巧

https://mp.weixin.qq.com/s/M0R3WdrpcPtVV68rlsNZuQ（ABACUS新进展：用密度泛函理论模拟千原子以上半导体电极表面反应）

https://mcresearch.github.io/abacus-user-guide/abacus-sdft.html（ABACUS 随机波函数 DFT 方法使用教程）

https://mcresearch.github.io/abacus-user-guide/abacus-sdft\_cond.html（ABACUS 随机波函数 DFT 计算电子电导热导教程）

https://mp.weixin.qq.com/s/oOF-Zel5ufqUx8ahZ0mpcw（ABACUS新进展：用混合随机密度泛函理论方法模拟极端高温物质电子性质）

https://www.bohrium.com/notebooks/5915692245 （ABACUS 随机波函数 DFT 方法使用教程 ）

* band-unfolding

https://mp.weixin.qq.com/s/auO0I5gjhPP-lxKY6evhRQ（PYATB 能带模块介绍）

https://mp.weixin.qq.com/s/cPf6nxXihtprnISUXBmbYw（ABACUS新进展：能带反折叠（band unfolding）方法及相关应用）

https://mp.weixin.qq.com/s/llKi5KG81Txa5BebnmzLng （Band unfolding 介绍及在PYATB中的使用）

https://www.bohrium.com/notebooks/2012704420（ABACUS+pyatb 能带反折叠计算 ）

使用案例：

https://www.bohrium.com/notebooks/45535812168（ABACUS使用案例：铀U及其氧化物的ABACUS计算）

https://www.bohrium.com/notebooks/92616245231（ABACUS使用教程（案例）：锰元素 by liu-ws）

https://www.bohrium.com/notebooks/16627233825（ABACUS使用案例 | 使用ABACUS氮及其氧化物进行密度泛函计算 by 王晨阳）

https://www.bohrium.com/notebooks/97366952314（ABACUS使用案例 | 零基础开始使用ABACUS对钡及其氧化物进行密度泛函计算 by roujin）

（TBC）

### **4. 预/后处理：**

https://abacus.deepmodeling.com/en/latest/advanced/interface/index.html（Interfaces to Other Softwares）

https://mcresearch.github.io/abacus-user-guide/abacus-dos.html（ABACUS+Atomkit 计算态密度和能带）

https://www.bohrium.com/notebooks/6116471155（ABACUS+Phono3py 计算晶格热导率）

https://www.bohrium.com/notebooks/8741867512（ABACUS+Phonopy 计算声子谱）

https://mcresearch.github.io/abacus-user-guide/abacus-phonopy.html （ABACUS+Phonopy 计算声子谱 ）

https://www.bohrium.com/notebooks/2712467526（ABACUS+ShengBTE 计算晶格热导率 ）https://mcresearch.github.io/abacus-user-guide/abacus-shengbte.html（ABACUS+ShengBTE 计算晶格热导率）

https://www.bohrium.com/notebooks/6116401077（从 DFT 到 MD｜超详细「深度势能」材料计算上手指南 ）https://www.bohrium.com/notebooks/8041860882（ABACUS+LibRI 做杂化泛函计算教程 ）

https://mcresearch.github.io/abacus-user-guide/abacus-libri.html（ABACUS+LibRI 做杂化泛函计算教程）

https://www.bohrium.com/notebooks/2912697542（ABACUS+Candela 分析分子动力学轨迹教程 ）

https://mcresearch.github.io/abacus-user-guide/abacus-candela.html（ABACUS+Candela 使用教程）

http://mcresearch.github.io/abacus-user-guide/abacus-bader.html （ABACUS+Bader charge 分析教程）

https://mcresearch.github.io/abacus-user-guide/abacus-elastic.html（ABACUS+pymatgen 计算弹性常数 ）https://mcresearch.github.io/abacus-user-guide/abacus-uspex.html（ABACUS+进化类晶体结构预测算法 USPEX 接口教程）

https://mcresearch.github.io/abacus-user-guide/abacus-namd.html（ABACUS+自主可控的激发态动力学软件Hefei NAMD 使用教程 ）

https://www.bohrium.com/notebooks/15169213382（使用ABACUS结合abTEM进行STEM模拟 ）

https://www.bohrium.com/notebooks/6516485694（ASE-ABACUS | 第一章：使用方法简介）

https://www.bohrium.com/notebooks/39369325971 （ASE-ABACUS | 第二章：NEB过渡态计算与ATST-Tools）

https://www.bohrium.com/notebooks/29581597682（ASE-ABACUS | 第三章：单端过渡态搜索）


# 迭代中的教程列表（内有错误内容）

[使用ABACUS软件进行DFT+U（密度泛函理论+Hubbard U校正）计算](https://ucoyxk075n.feishu.cn/wiki/UwVxwKF7biMZ1okh1vDcRj0Pnqf)

[ABACUS结合wannier90的计算教程](https://ucoyxk075n.feishu.cn/wiki/MgS7wm3EHiHC1jkFIn3cDBHXnZj)

[ABACUS原子轨道基组杂化泛函计算教程](https://ucoyxk075n.feishu.cn/wiki/P2y2w1LsaiWumPk8vr2cD5GqnPg)

[基于ABACUS的弹性常数计算方法与实践](https://ucoyxk075n.feishu.cn/wiki/A7k3wwIARiTJiOkdJyfcZedlnlh)
