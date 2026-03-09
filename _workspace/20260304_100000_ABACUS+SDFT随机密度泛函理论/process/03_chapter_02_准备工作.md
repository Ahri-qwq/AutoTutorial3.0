## 二、准备工作

### 2.1 软件版本

SDFT 功能从 ABACUS 2.3.0 版本开始支持，DOS 计算和多 K 点采样等完整功能从 3.2.0 版本起可用。本教程基于 ABACUS 3.2.0 版本的算例。

### 2.2 算例下载

官方算例包可从 Gitee 获取：

```bash
git clone https://gitee.com/mcresearch/abacus-user-guide.git
```

或从 GitHub 获取：

```bash
git clone https://github.com/MCresearch/abacus-user-guide.git
```

下载后进入随机波函数算例目录：

```bash
cd abacus-user-guide/examples/stochastic
```

目录中包含三个文件夹：

```
stochastic/
├── pw_Si2/                  # 案例一：Si SCF
├── pw_md_Al/                # 案例二：Al MD
└── 186_PW_SDOS_10D10S/      # 案例三：Si DOS
```

每个文件夹均包含 `INPUT`、`STRU`、`KPT` 三个输入文件，与常规 KSDFT 计算的文件结构完全相同。赝势文件位于 `../../PP_ORB/` 目录（相对于各算例文件夹）。

### 2.3 赝势说明

案例一使用 `Si.pz-vbc.UPF`，包含硅的 4 个价电子。

当电子温度特别高（如超过 100 eV）时，内壳层电子可能被热激发，常规赝势的可移植性会下降。此时需要选用包含更多内壳层电子的赝势，甚至需要定制生成新的赝势。ABACUS 目前支持模守恒赝势（NCPP）。
