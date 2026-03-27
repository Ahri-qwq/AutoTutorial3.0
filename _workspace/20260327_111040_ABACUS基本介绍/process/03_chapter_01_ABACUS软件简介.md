# 第一章：ABACUS 软件简介

## 1.1 软件背景

ABACUS（Atomic-orbital Based Ab-initio Computation at UStc），中文名**原子算筹**，是国内自主研发的开源第一性原理材料计算软件。软件基于密度泛函理论（Density Functional Theory，DFT），采用模守恒赝势和周期性边界条件，支持平面波（Plane Wave，PW）和数值原子轨道（Numerical Atomic Orbital，NAO/LCAO）两种基组。

ABACUS 托管于 DeepModeling 开源社区，代码完全公开，任何人均可参与贡献：

- **GitHub 仓库：** https://github.com/deepmodeling/abacus-develop
- **官方文档：** https://abacus.deepmodeling.com/en/latest/
- **中文教程：** https://mcresearch.github.io/abacus-user-guide/

软件免费开源，可在 Linux 集群、云计算平台（如 Bohrium）上运行。

## 1.2 核心功能

ABACUS 目前支持以下主要计算类型：

| 计算类型 | 关键字 | 说明 |
|----------|--------|------|
| 电子自洽迭代 | `scf` | 求解电子基态，获得总能量、电荷密度 |
| 非自洽计算 | `nscf` | 基于已有电荷密度计算能带/DOS |
| 结构优化 | `relax` | 固定晶胞，优化原子位置 |
| 晶胞优化 | `cell-relax` | 同时优化晶胞参数和原子位置 |
| 分子动力学 | `md` | 模拟原子随时间的运动轨迹 |

除基本功能外，ABACUS 还支持：
- DFT+U（处理强关联体系）
- 范德华修正（vdW correction）
- 隐式溶剂模型
- 杂化泛函（HSE06 等）
- 自旋极化与非共线磁矩计算
- 与 DeePMD-kit 联用的机器学习分子动力学

## 1.3 两种基组：PW 与 LCAO

ABACUS 同时支持两种基组，各有适用场景：

| 对比项 | 平面波（PW） | 数值原子轨道（LCAO） |
|--------|-------------|----------------------|
| 控制参数 | `basis_type pw` | `basis_type lcao` |
| 额外文件 | 仅赝势（.upf） | 赝势（.upf）+ 轨道（.orb）|
| 精度控制 | ecutwfc 截断能 | 轨道文件精度 |
| 计算速度 | 中等 | 通常更快（尤其大体系）|
| 适用场景 | 通用，精度基准 | 大体系、高通量计算 |

初学者建议先用 PW 基组入手，熟悉后再尝试 LCAO。

## 1.4 ABACUS 计算文件体系

运行一个 ABACUS 计算，至少需要准备以下三个文件：

```
计算目录/
├── INPUT        # 计算参数控制（必需）
├── KPT          # 布里渊区 k 点采样（必需）
├── STRU         # 晶体结构信息（必需）
├── *.upf        # 赝势文件（必需）
└── *.orb        # 数值轨道文件（LCAO 基组才需要）
```

计算完成后，输出文件保存在：

```
计算目录/
├── OUT.ABACUS/           # 输出文件夹（ABACUS 为默认后缀）
│   ├── INPUT             # 包含所有参数（含默认值）的完整 INPUT
│   ├── running_scf.log   # 详细运行日志
│   ├── KPT.info          # k 点信息
│   └── ...               # 其他输出文件
└── time.json             # 各模块计时信息
```

下面三章分别介绍这三个核心输入文件。
