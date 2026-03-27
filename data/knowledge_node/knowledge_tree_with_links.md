# ABACUS 知识树（含教程链接版）

> 生成日期：2026-03-10
> 来源文档：`data/knowledge_node/knowledge_node.md`
> 新教程目录：`data/knowledge_node/Tutorials/`（共 11 篇）
> 说明：本文件在 `knowledge_tree.md` 基础上，为每个知识点补充了对应教程链接。

---

## 图例

| 符号 | 含义 |
|------|------|
| ✅ | 新文章已生成 — `Tutorials/` 中已有对应新教程 |
| ⬜ | 有旧教程 — knowledge_node.md 中附有外部链接的旧教程 |
| 🔄 | 迭代中旧教程 — 飞书列表中含错误、正在更新的教程 |
| ❌ | 缺少教程 — 无任何教程，尚待编写 |
| 🔗 | 教程链接行（链接格式：`[教程简称](URL)`，多个链接用 ` \| ` 分隔） |
> 新教程路径为相对路径，如 `data/knowledge_node/Tutorials/07_Final_Tutorial_XXX.md`

---

## 覆盖统计

| 状态 | 知识点数 | 说明 |
|------|---------|------|
| ✅ 新文章已生成（含组合） | 16 | 其中 4 项仅有新文章，12 项同时有旧教程 |
| ⬜ 有旧教程（含组合） | 81 | 其中 69 项仅有旧教程，12 项与 ✅/🔄 组合 |
| 🔄 迭代中旧教程（含组合） | 6 | 其中 2 项已有新文章替代，4 项尚未替代 |
| ❌ 缺少教程 | 15 | 无任何教程，是当前主要空缺 |
| 合计知识点 | ~104 | 含组合标注，故各行相加超过 104 |

---

## 知识树

### 一、基础篇：安装与基本计算

#### 1. 软件与环境

- 架构介绍（ABACUS 发展历史、定位与特点）⬜
  - 🔗 [ABACUS论文 (arXiv)](https://arxiv.org/abs/2501.08697) | [携手DeepModeling (微信)](https://mp.weixin.qq.com/s/Fx5jx99g9PVjgnLAsOyQ9Q) | [AI4S时代电子结构软件 (微信)](https://mp.weixin.qq.com/s/V2WVr_O2_LJ-VJ09MaHaqQ)
- 安装与编译 ⬜
  - 🔗 [Easy Installation (官方文档)](https://abacus.deepmodeling.com/en/latest/quick_start/easy_install.html) | [Advanced Installation (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/install.html)
  - 依赖库（MPI、LibXC、ELPA 等）⬜
    - 🔗 [ABACUS 依赖软件库](https://mcresearch.github.io/abacus-user-guide/abacus-hpc.html#21-abacus-%E4%BE%9D%E8%B5%96%E7%9A%84%E8%BD%AF%E4%BB%B6%E5%BA%93)
  - CPU 编译方法（GCC / Intel oneAPI）⬜
    - 🔗 [GCC 编译](https://mcresearch.github.io/abacus-user-guide/abacus-gcc.html) | [Intel oneAPI 2024/2025](https://mcresearch.github.io/abacus-user-guide/abacus-oneapi.html) | [Intel oneAPI (旧版)](https://mcresearch.github.io/abacus-user-guide/abacus-intel.html)
  - GPU 编译方法（Nvidia CUDA）⬜
    - 🔗 [Nvidia GPU 版本编译](https://mcresearch.github.io/abacus-user-guide/abacus-gpu.html) | [LCAO GPU 版本使用](https://mcresearch.github.io/abacus-user-guide/abacus-gpu-lcao.html)
  - DCU 编译方法（曙光 DCU / AMD）⬜
    - 🔗 [曙光 DCU/AMD 编译与使用](https://mcresearch.github.io/abacus-user-guide/abacus-dcu.html)
  - Toolchain 安装教程（GNU / Intel / AMD / 离线 / GPU）⬜
    - 🔗 [一键配置脚本 (Bohrium)](https://www.bohrium.com/notebooks/5215742477) | [GNU (微信)](https://mp.weixin.qq.com/s/ypc0RT5ePm0vMlGRYI46JA) | [Intel (微信)](https://mp.weixin.qq.com/s/K58BVQwSoxcgcNufODo8tg) | [AMD (微信)](https://mp.weixin.qq.com/s/S0EfqhNn-FLh_QXRZE2bqw) | [离线安装 (微信)](https://mp.weixin.qq.com/s/aym7LBcgfutDtlKOR0aevA) | [GPU (微信)](https://mp.weixin.qq.com/s/YekI3LXy7vuTN36ut5M6Rg)
  - 无 MPI 编译 ⬜
    - 🔗 [编译无 MPI 的 ABACUS (飞书)](https://xmywuqhxb0.feishu.cn/docx/JCv0dHPP6o69JdxtG33cIcTnnke)
- 常见报错与解决方法 ⬜
  - 🔗 [Frequently Asked Questions (官方文档)](https://abacus.deepmodeling.com/en/latest/community/faq.html)
- 晶胞朝向与并行效率优化 ⬜
  - 一维材料（以碳纳米管为例）⬜
    - 🔗 [晶胞朝向与并行效率：碳纳米管](https://mcresearch.github.io/abacus-user-guide/abacus-eff1.html)
  - 二维材料（以氮化硼为例）⬜
    - 🔗 [晶胞朝向与并行效率：二维氮化硼](https://mcresearch.github.io/abacus-user-guide/abacus-eff2.html)
  - 表面体系（以铜表面 CO 吸附为例）⬜
    - 🔗 [晶胞朝向与并行效率：铜表面 CO 吸附](https://mcresearch.github.io/abacus-user-guide/abacus-eff3.html)

#### 2. 输入输出体系

- 主要输入文件：INPUT / STRU / KPT ⬜
  - 🔗 [Input Files 简介 (官方文档)](https://abacus.deepmodeling.com/en/latest/quick_start/input.html) | [Input Files 详细说明 (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/input_files/index.html) | [STRU 文件入门 (微信)](https://mp.weixin.qq.com/s/FysreUHIRB5RHKDtoj99vQ) | [STRU 文件转换 (Bohrium)](https://www.bohrium.com/notebooks/9814968648)
- 运行输出解读（能量收敛、电子步收敛、力与应力）⬜
  - 🔗 [Output Files 简介 (官方文档)](https://abacus.deepmodeling.com/en/latest/quick_start/output.html)
- 参数设置经验
  - k 点设置（Monkhorst-Pack 网格、收敛性测试）⬜
    - 🔗 [平面波计算与收敛性测试](https://mcresearch.github.io/abacus-user-guide/abacus-pw.html)
  - 平面波截断能（能量收敛、赝势差异）⬜
    - 🔗 [平面波计算与收敛性测试](https://mcresearch.github.io/abacus-user-guide/abacus-pw.html)
  - 赝势选择（ONCV vs HGH，元素库差异）⬜
    - 🔗 [APNS-PPORB-v1 赝势轨道捆绑包 (微信)](https://mp.weixin.qq.com/s/vUG-7uRjbMZJP2fs8Ep_8Q)
  - 基组选择（平面波 PW vs 数值原子轨道 NAO）⬜
    - SZ / DZP / TZDP 等基组差异 ⬜
      - 🔗 [NAO 命名与使用方法 (Bohrium)](https://www.bohrium.com/notebooks/9319634192) | [NAO 命名与使用方法 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-nac1.html)
    - 基组收敛测试 ❌
    - 基组与赝势匹配原则 ❌
    - 数值原子轨道生成（模守恒赝势 + 优化流程）⬜
      - 🔗 [NAO 生成（二）(Bohrium)](https://www.bohrium.com/notebooks/5215642163) | [NAO 生成（二）(mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-nac2.html) | [高精度 NAO（三）(Bohrium)](https://www.bohrium.com/notebooks/8841868194) | [高精度 NAO（三）(mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-nac3.html) | [模守恒赝势生成](https://mcresearch.github.io/abacus-user-guide/abacus-upf.html)
  - SCF 收敛算法相关参数 ⬜
    - 🔗 [收敛性问题解决手册 (飞书)](https://ucoyxk075n.feishu.cn/docx/R0sqdk6T0o2RY4x5IWgcJ3RHnug)

#### 3. 基础计算实践

- 结构优化
  - 固定晶胞优化（ionic relax）⬜
    - 🔗 [结构优化教程 (Bohrium)](https://www.bohrium.com/notebooks/9119461238) | [快速开始 ABACUS (Bohrium)](https://www.bohrium.com/notebooks/4641406377) | [Geometry Optimization (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/opt.html)
  - 变胞优化（cell-relax）⬜
    - 🔗 [结构优化教程 (Bohrium)](https://www.bohrium.com/notebooks/9119461238) | [快速开始 ABACUS (Bohrium)](https://www.bohrium.com/notebooks/4641406377)
- 电子结构
  - 自洽计算（SCF）⬜
    - 🔗 [Running SCF (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/scf/index.html) | [快速开始 ABACUS (Bohrium)](https://www.bohrium.com/notebooks/4641406377)
  - 非自洽计算（NSCF，读入电荷密度）⬜
    - 🔗 [如何正确画能带，NSCF 读电荷密度 (飞书)](https://xmywuqhxb0.feishu.cn/docx/K8GRdTst4oXQNoxnQVbcFZTmntb)
  - 能带结构计算 ⬜
    - 🔗 [Extracting Band Structure (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/band.html) | [快速开始 ABACUS (Bohrium)](https://www.bohrium.com/notebooks/4641406377)
  - 态密度（DOS / PDOS）⬜
    - 🔗 [Calculating DOS and PDOS (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/dos.html) | [ABACUS 计算 PDOS (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-pdos.html) | [DOS 和 PDOS 计算 (飞书)](https://xmywuqhxb0.feishu.cn/docx/ONSldj82VoNGKSxaoDQcoKBtnGh)
- 分子动力学（基础 MD）
  - AIMD 基本概念 ⬜
    - 🔗 [分子动力学使用教程 (Bohrium)](https://www.bohrium.com/notebooks/2241262724) | [分子动力学使用教程 (mcresearch)](http://mcresearch.github.io/abacus-user-guide/abacus-md.html) | [MD 官方文档](https://abacus.deepmodeling.com/en/latest/advanced/md.html)
  - NVT 系综模拟 ⬜
    - 🔗 [分子动力学使用教程 (Bohrium)](https://www.bohrium.com/notebooks/2241262724) | [分子动力学使用教程 (mcresearch)](http://mcresearch.github.io/abacus-user-guide/abacus-md.html)
  - 输出文件与轨迹分析 ⬜
    - 🔗 [Candela 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-candela.html) | [Candela 分析 MD 轨迹 (Bohrium)](https://www.bohrium.com/notebooks/2912697542)
- 电荷密度与波函数
  - 电荷密度分布 ⬜
    - 🔗 [Extracting Charge Density (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/charge.html) | [电荷密度和波函数可视化 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-chg.html)
  - Mulliken 电荷分析 ⬜
    - 🔗 [Mulliken Charge Analysis (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/Mulliken.html)
  - 静电势（Electrostatic potential）⬜
    - 🔗 [Extracting Electrostatic Potential (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/potential.html)
  - 波函数可视化与后处理 ⬜
    - 🔗 [Extracting Wave Functions (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/wfc.html) | [电荷密度和波函数可视化 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-chg.html)
  - 电子局域函数（ELF）⬜ ✅
    - 🔗 [ELF 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-elf.html) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_ABACUS使用教程-ELF电子局域函数计算与可视化.md)
- 基于 Bohrium 平台的运用 ⬜
  - 🔗 [Bohrium 帮助文档｜ABACUS](https://www.bohrium.com/notebooks/6643676733)

---

### 二、进阶篇：高级功能与物性计算

#### 1. 晶体材料

- 对称性分析（周期性、点阵、惯用胞与原胞、点群、空间群）❌
- 结合能（基于 DFT 的能量最低构型搜索）❌
- 带隙计算
  - 基础 PBE 带隙计算流程 ❌
  - DFT+U 带隙修正 🔄 ✅
    - 🔗 [旧教程 (飞书)](https://ucoyxk075n.feishu.cn/wiki/UwVxwKF7biMZ1okh1vDcRj0Pnqf) | [DFT+U 计算 (Bohrium)](https://www.bohrium.com/notebooks/2112617648) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_ABACUS_DFT+U强关联体系计算.md)
  - 杂化泛函（Hybrid functional，PBE0 / HSE06）🔄 ⬜
    - 🔗 [旧教程 (飞书)](https://ucoyxk075n.feishu.cn/wiki/P2y2w1LsaiWumPk8vr2cD5GqnPg) | [PW 基组杂化泛函 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-exx.html) | [LibRI 杂化泛函 (Bohrium)](https://www.bohrium.com/notebooks/8041860882) | [LibRI 杂化泛函 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-libri.html)
  - DFT+1/2 方法 ❌
- 拓扑计算
  - 材料拓扑性质基础知识 ⬜
    - 🔗 [薄膜技术增强非线性霍尔效应 (微信)](https://mp.weixin.qq.com/s/WQaVzaXZ86ShjAtCcFmWjQ)
  - ABACUS+PYATB 能带拓扑计算流程 ⬜
    - 🔗 [PYATB 能带模块介绍 (微信)](https://mp.weixin.qq.com/s/auO0I5gjhPP-lxKY6evhRQ)
  - ABACUS+Wannier90+WannierTools 🔄 ⬜
    - 🔗 [旧教程 (飞书)](https://ucoyxk075n.feishu.cn/wiki/MgS7wm3EHiHC1jkFIn3cDBHXnZj) | [Wannier90 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-wannier.html) | [最大局域化 Wannier 函数简介 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/algorithm-wannier.html)
- 磁性计算 ⬜
  - 🔗 [磁性材料计算 (Bohrium)](https://www.bohrium.com/notebooks/7141761751) | [磁交换参数计算 (微信)](https://mp.weixin.qq.com/s/64XicP_bionUZOYMeXvc-g) | [磁交换作用与微磁学 (微信)](https://mp.weixin.qq.com/s/8ZO682BIJy7w0ge-qJB_9g)
  - 电子自旋-原子磁矩-晶体宏观磁性（铁磁 / 反铁磁 / 亚铁磁 / 顺磁）⬜
    - 🔗 [磁性材料计算 (Bohrium)](https://www.bohrium.com/notebooks/7141761751)
  - 自旋极化计算 ⬜
    - 🔗 [磁性材料计算 (Bohrium)](https://www.bohrium.com/notebooks/7141761751)
  - 非共线自旋 ⬜
    - 🔗 [磁性材料计算 (Bohrium)](https://www.bohrium.com/notebooks/7141761751)
  - 自旋轨道耦合（SOC）⬜
    - 🔗 [磁性材料计算 (Bohrium)](https://www.bohrium.com/notebooks/7141761751)
  - 磁交换参数与各向异性能量 ⬜
    - 🔗 [磁性交换作用与微磁学结合 (微信)](https://mp.weixin.qq.com/s/8ZO682BIJy7w0ge-qJB_9g) | [磁交换参数计算 (微信)](https://mp.weixin.qq.com/s/64XicP_bionUZOYMeXvc-g)
- 铁电极化（Berry 相位方法）⬜
  - 🔗 [Berry Phase Calculation (官方文档)](https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/Berry_phase.html) | [Berry Curvature Dipole (微信)](https://mp.weixin.qq.com/s/_56PIb94LVDXOxKvqSW-lg)
- 光学性质分析 ⬜ ✅
  - 🔗 [PYATB 能带模块介绍 (微信)](https://mp.weixin.qq.com/s/auO0I5gjhPP-lxKY6evhRQ) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_光学性质计算.md)
  - 光吸收谱 ✅
    - 🔗 [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_光学性质计算.md)
  - ABACUS+PYATB 介电函数与线性光学 ⬜
    - 🔗 [介电函数与线性光学 (Bohrium)](https://www.bohrium.com/notebooks/35791425971)
  - GW 方法 ❌
- 实时含时密度泛函（RT-TDDFT）⬜ ✅
  - 🔗 [TDDFT 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-tddft.html) | [混合规范 RT-TDDFT (微信)](https://mp.weixin.qq.com/s/uMg2jMgDSzrP1VYqLhobog) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_混合规范RTTDDFT.md)
  - 速度规范 / 长度规范 / 混合规范 ✅
    - 🔗 [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_混合规范RTTDDFT.md)
  - 外场激发（泵浦-探测、电流与极化响应）❌
  - 超快动力学应用 ❌
- 声子谱与热力学性质（ABACUS+Phonopy）⬜ ✅
  - 🔗 [Phonopy 计算声子谱 (Bohrium)](https://www.bohrium.com/notebooks/8741867512) | [Phonopy 计算声子谱 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-phonopy.html) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_使用ABACUS+Phonopy计算声子谱.md)
- 分子动力学进阶（AIMD）
  - NVT / NPT 系综，复杂体系模拟 ⬜
    - 🔗 [分子动力学使用教程 (Bohrium)](https://www.bohrium.com/notebooks/2241262724) | [分子动力学使用教程 (mcresearch)](http://mcresearch.github.io/abacus-user-guide/abacus-md.html)
- 机器学习势函数
  - DPGEN（数据生成与迭代）⬜
    - 🔗 [DPGEN 快速上手 (Bohrium)](https://www.bohrium.com/notebooks/26592976824) | [DPGEN 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-dpgen.html)
  - DeePMD-kit（势函数训练）⬜
    - 🔗 [DeePMD-kit 机器学习 MD (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-dpmd.html) | [DeePKS-ES 介绍及使用 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-deepks-es.html)
  - 模型微调与蒸馏 ❌
  - ABACUS+DP 模型分子动力学 ⬜
    - 🔗 [DeePMD-kit 机器学习 MD (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-dpmd.html) | [DFT→MD 深度势能上手指南 (Bohrium)](https://www.bohrium.com/notebooks/6116401077)

#### 2. 表面与界面材料

- 表面建模（真空层的影响、slab 模型构建）❌
- 功函数（静电势 + 真空能级）⬜
  - 🔗 [表面计算（一）：静电势和功函数](https://mcresearch.github.io/abacus-user-guide/abacus-surface1.html)
- 偶极修正（Dipole correction）⬜
  - 🔗 [表面计算（二）：偶极修正](https://mcresearch.github.io/abacus-user-guide/abacus-surface2.html)
- 表面能计算 ⬜
  - 🔗 [表面计算（三）：表面能计算](https://mcresearch.github.io/abacus-user-guide/abacus-surface3.html)
- 表面缺陷能与吸附能 ⬜
  - 🔗 [表面计算（四）：表面缺陷能和吸附能](https://mcresearch.github.io/abacus-user-guide/abacus-surface4.html)
- 外加电场 ⬜
  - 🔗 [表面计算（五）：外加电场](https://mcresearch.github.io/abacus-user-guide/abacus-surface5.html)
- 补偿电荷 ⬜
  - 🔗 [表面计算（六）：补偿电荷](https://mcresearch.github.io/abacus-user-guide/abacus-surface6.html)
- NEB 过渡态搜索（ATST-Tools / AutoNEB）⬜ ✅
  - 🔗 [ASE-ABACUS NEB 过渡态计算 (Bohrium)](https://www.bohrium.com/notebooks/39369325971) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_使用ABACUS和ATST-Tools进行NEB-AutoNEB过渡态搜索.md)
- 层错能（Stacking fault energy）❌
- 隐式溶剂模型（Implicit solvation）⬜ ✅
  - 🔗 [隐式溶剂模型使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-sol.html) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_ABACUS隐式溶剂模型使用.md)

#### 3. 大体系材料计算

- 缺陷与掺杂
  - 点缺陷建模（超胞方法）❌
  - 形成能与电荷态转变能级 ❌
  - 宽禁带半导体缺陷调控 ❌
- 随机波函数 DFT（SDFT / MDFT）⬜ ✅
  - 🔗 [SDFT 使用教程 (Bohrium)](https://www.bohrium.com/notebooks/5915692245) | [SDFT 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-sdft.html) | [SDFT 计算电导热导 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-sdft_cond.html) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_ABACUS随机密度泛函理论SDFT使用教程.md)
- 能带反折叠（Band Unfolding）⬜
  - 🔗 [Band Unfolding 介绍 (微信)](https://mp.weixin.qq.com/s/cPf6nxXihtprnISUXBmbYw) | [PYATB 能带反折叠 (Bohrium)](https://www.bohrium.com/notebooks/2012704420) | [Band Unfolding 在 PYATB 中的使用 (微信)](https://mp.weixin.qq.com/s/llKi5KG81Txa5BebnmzLng)
- 无轨道密度泛函理论（OFDFT）⬜
  - 🔗 [OFDFT 使用教程 (Bohrium)](https://www.bohrium.com/notebooks/6416644691) | [OFDFT 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-ofdft.html)
- 晶格热导率
  - ABACUS+Phonopy（声子谱 + 二阶力常数）⬜ ✅
    - 🔗 [Phonopy 声子谱 (Bohrium)](https://www.bohrium.com/notebooks/8741867512) | [Phonopy 声子谱 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-phonopy.html) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_使用ABACUS+Phonopy计算声子谱.md)
  - ABACUS+ShengBTE（三阶力常数 + 热导率）⬜ ✅
    - 🔗 [ShengBTE 晶格热导率 (Bohrium)](https://www.bohrium.com/notebooks/2712467526) | [ShengBTE 晶格热导率 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-shengbte.html) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_使用ABACUS+ShengBTE计算Si晶格热导率.md)
  - ABACUS+Phono3py ⬜
    - 🔗 [Phono3py 晶格热导率 (Bohrium)](https://www.bohrium.com/notebooks/6116471155)

#### 4. 预/后处理：接口与工具

- abacustest 工具 🔄 ✅
  - 快速准备 ABACUS 输入文件 ✅
    - 🔗 [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_使用abacustest快速准备ABACUS输入文件.md)
  - 弹性常数计算（abacustest + pymatgen）🔄 ✅
    - 🔗 [旧教程 (飞书)](https://ucoyxk075n.feishu.cn/wiki/A7k3wwIARiTJiOkdJyfcZedlnlh) | [弹性常数计算 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-elastic.html) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_使用abacustest计算晶体弹性性质.md)
- ABACUS+Atomkit（DOS 与能带后处理）⬜
  - 🔗 [Atomkit 计算 DOS 和能带 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-dos.html)
- ABACUS+LibRI（LCAO 基组杂化泛函）🔄 ⬜
  - 🔗 [旧教程 (飞书)](https://ucoyxk075n.feishu.cn/wiki/P2y2w1LsaiWumPk8vr2cD5GqnPg) | [LibRI 杂化泛函 (Bohrium)](https://www.bohrium.com/notebooks/8041860882) | [LibRI 杂化泛函 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-libri.html)
- ABACUS+Candela（MD 轨迹分析）⬜
  - 🔗 [Candela 分析 MD 轨迹 (Bohrium)](https://www.bohrium.com/notebooks/2912697542) | [Candela 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-candela.html)
- ABACUS+Bader Charge 分析 ⬜
  - 🔗 [Bader Charge 分析教程 (mcresearch)](http://mcresearch.github.io/abacus-user-guide/abacus-bader.html)
- ABACUS+USPEX（进化类晶体结构预测）⬜
  - 🔗 [USPEX 接口教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-uspex.html)
- ABACUS+Hefei NAMD（激发态动力学）⬜
  - 🔗 [Hefei NAMD 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-namd.html)
- ABACUS+abTEM（STEM 模拟）⬜
  - 🔗 [abTEM STEM 模拟 (Bohrium)](https://www.bohrium.com/notebooks/15169213382)
- ASE-ABACUS 接口 ⬜
  - 使用方法简介 ⬜
    - 🔗 [ASE-ABACUS 第一章：使用方法简介 (Bohrium)](https://www.bohrium.com/notebooks/6516485694)
  - NEB 过渡态计算 ✅
    - 🔗 [ASE-ABACUS 第二章：NEB 过渡态计算与 ATST-Tools (Bohrium)](https://www.bohrium.com/notebooks/39369325971) | [新教程](data/knowledge_node/Tutorials/07_Final_Tutorial_使用ABACUS和ATST-Tools进行NEB-AutoNEB过渡态搜索.md)
  - 单端过渡态搜索（Dimer / ART 等）⬜
    - 🔗 [ASE-ABACUS 第三章：单端过渡态搜索 (Bohrium)](https://www.bohrium.com/notebooks/29581597682)
- ABACUS+DPA-3（基于 Bohrium 的安装与使用）⬜
  - 🔗 [DPA-3 安装与使用 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-dpa3-toturial.html) | [DPA-3 教程 (Bohrium)](https://www.bohrium.com/notebooks/6643676733)
- ABACUS+DPGEN（机器学习势函数生成）⬜
  - 🔗 [DPGEN 使用教程 (mcresearch)](https://mcresearch.github.io/abacus-user-guide/abacus-dpgen.html) | [DPGEN 快速上手 (Bohrium)](https://www.bohrium.com/notebooks/26592976824)

---

## 新教程对应关系（11 篇）

| 新教程 | 对应知识点位置 | 文件路径 |
|-------|--------------|---------|
| 使用abacustest快速准备ABACUS输入文件 | 进阶篇 > 预/后处理 > abacustest | `data/knowledge_node/Tutorials/07_Final_Tutorial_使用abacustest快速准备ABACUS输入文件.md` |
| 使用abacustest计算晶体弹性性质 | 进阶篇 > 预/后处理 > abacustest > 弹性常数 | `data/knowledge_node/Tutorials/07_Final_Tutorial_使用abacustest计算晶体弹性性质.md` |
| ABACUS_DFT+U强关联体系计算 | 进阶篇 > 晶体材料 > 带隙计算 > DFT+U | `data/knowledge_node/Tutorials/07_Final_Tutorial_ABACUS_DFT+U强关联体系计算.md` |
| ABACUS隐式溶剂模型使用 | 进阶篇 > 表面与界面 > 隐式溶剂模型 | `data/knowledge_node/Tutorials/07_Final_Tutorial_ABACUS隐式溶剂模型使用.md` |
| 使用ABACUS+Phonopy计算声子谱 | 进阶篇 > 晶体材料 > 声子谱；大体系 > 晶格热导率 | `data/knowledge_node/Tutorials/07_Final_Tutorial_使用ABACUS+Phonopy计算声子谱.md` |
| 使用ABACUS+ShengBTE计算Si晶格热导率 | 进阶篇 > 大体系 > 晶格热导率 > ShengBTE | `data/knowledge_node/Tutorials/07_Final_Tutorial_使用ABACUS+ShengBTE计算Si晶格热导率.md` |
| ABACUS随机密度泛函理论SDFT使用教程 | 进阶篇 > 大体系 > SDFT | `data/knowledge_node/Tutorials/07_Final_Tutorial_ABACUS随机密度泛函理论SDFT使用教程.md` |
| 混合规范RTTDDFT | 进阶篇 > 晶体材料 > RT-TDDFT > 混合规范 | `data/knowledge_node/Tutorials/07_Final_Tutorial_混合规范RTTDDFT.md` |
| ABACUS使用教程-ELF电子局域函数计算与可视化 | 基础篇 > 电荷密度与波函数 > ELF | `data/knowledge_node/Tutorials/07_Final_Tutorial_ABACUS使用教程-ELF电子局域函数计算与可视化.md` |
| 使用ABACUS和ATST-Tools进行NEB-AutoNEB过渡态搜索 | 进阶篇 > 表面与界面 > NEB；预/后处理 > ASE-ABACUS | `data/knowledge_node/Tutorials/07_Final_Tutorial_使用ABACUS和ATST-Tools进行NEB-AutoNEB过渡态搜索.md` |
| 光学性质计算 | 进阶篇 > 晶体材料 > 光学性质 | `data/knowledge_node/Tutorials/07_Final_Tutorial_光学性质计算.md` |

---

## 迭代中旧教程（4 篇，来自飞书列表）

| 飞书教程标题 | 对应知识点 | 是否有新文章替代 | 飞书链接 |
|------------|-----------|--------------|---------|
| 使用ABACUS软件进行DFT+U计算 | 晶体材料 > DFT+U | ✅ 已生成 | [飞书链接](https://ucoyxk075n.feishu.cn/wiki/UwVxwKF7biMZ1okh1vDcRj0Pnqf) |
| ABACUS结合wannier90的计算教程 | 晶体材料 > 拓扑 > Wannier90 | ❌ 未生成 | [飞书链接](https://ucoyxk075n.feishu.cn/wiki/MgS7wm3EHiHC1jkFIn3cDBHXnZj) |
| ABACUS原子轨道基组杂化泛函计算教程 | 晶体材料 > 杂化泛函；预/后处理 > LibRI | ❌ 未生成 | [飞书链接](https://ucoyxk075n.feishu.cn/wiki/P2y2w1LsaiWumPk8vr2cD5GqnPg) |
| 基于ABACUS的弹性常数计算方法与实践 | 预/后处理 > abacustest > 弹性常数 | ✅ 已生成 | [飞书链接](https://ucoyxk075n.feishu.cn/wiki/A7k3wwIARiTJiOkdJyfcZedlnlh) |

---

## 缺少教程的知识点汇总（❌ 共 15 项）

> 供后续教程选题参考：

**基础篇：**
1. 基组收敛测试
2. 基组与赝势匹配原则

**进阶篇 — 晶体材料：**
3. 对称性分析（晶体学基础）
4. 结合能计算
5. 基础 PBE 带隙计算流程
6. DFT+1/2 方法
7. GW 方法
8. 外场激发（泵浦-探测）
9. 超快动力学应用
10. 机器学习势函数微调与蒸馏

**进阶篇 — 表面与界面：**
11. 表面建模（slab 构建、真空层）
12. 层错能计算

**进阶篇 — 大体系：**
13. 点缺陷建模（超胞方法）
14. 形成能与电荷态转变能级
15. 宽禁带半导体缺陷调控
