# 第一章：准备工作

从结构文件准备ABACUS输入文件夹，需要先完成工具安装和环境配置。本章将介绍如何安装abacustest工具，下载推荐的赝势轨道库，并配置环境变量。

## 1.1 安装abacustest

abacustest是ABACUS的前后处理工具，支持从结构文件快速准备完整的输入文件夹。安装方式有两种：

### 方式1：通过pip安装（推荐）

```bash
pip install abacustest
```

### 方式2：从源码安装

```bash
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

安装完成后，验证安装是否成功：

```bash
abacustest --version
```

如果显示版本号，说明安装成功。

## 1.2 下载APNS赝势轨道库

APNS（ABACUS Pseudopotential-NAO Square）是ABACUS官方推荐的赝势轨道库，经过大规模测试，兼顾精度和效率。abacustest提供了一键下载功能。

在你的工作目录下执行：

```bash
abacustest model inputs --download-pporb
```

命令执行完成后，当前目录会出现两个文件夹：

```bash
$ ls
apns-orbitals-efficiency-v1  apns-pseudopotentials-v1
```

这两个文件夹包含从H到Bi共83种元素的赝势和轨道文件：
- `apns-pseudopotentials-v1`：赝势文件（UPF格式）
- `apns-orbitals-efficiency-v1`：数值原子轨道文件（ORB格式，DZP水平）

**说明：**
- efficiency系列轨道适合结构优化、分子动力学、声子谱等计算
- 如需更高精度（如能带计算、激发态计算），可下载precision系列轨道
- 镧系元素的4f电子被视为核电子，轨道截断半径为8 au

## 1.3 配置环境变量

为了让abacustest自动找到赝势和轨道文件，需要设置环境变量。

### 临时设置（当前终端有效）

```bash
export ABACUS_PP_PATH=/your/path/to/apns-pseudopotentials-v1
export ABACUS_ORB_PATH=/your/path/to/apns-orbitals-efficiency-v1
```

将 `/your/path/to/` 替换为实际的绝对路径。

### 永久设置（推荐）

将上述命令添加到 `~/.bashrc` 或 `~/.bash_profile` 文件末尾：

```bash
echo 'export ABACUS_PP_PATH=/your/path/to/apns-pseudopotentials-v1' >> ~/.bashrc
echo 'export ABACUS_ORB_PATH=/your/path/to/apns-orbitals-efficiency-v1' >> ~/.bashrc
source ~/.bashrc
```

验证环境变量是否设置成功：

```bash
echo $ABACUS_PP_PATH
echo $ABACUS_ORB_PATH
```

如果显示正确的路径，说明配置成功。

**注意：**
- 如果只做平面波（PW）计算，只需设置 `ABACUS_PP_PATH`
- LCAO计算需要同时设置 `ABACUS_PP_PATH` 和 `ABACUS_ORB_PATH`
- 也可以在使用abacustest时通过 `--pp` 和 `--orb` 选项手动指定路径，不依赖环境变量

## 1.4 赝势轨道库的文件命名规则

abacustest能够自动识别赝势和轨道文件，前提是文件命名符合规则：

### 规则1：文件名以元素名开头

例如：
- 赝势：`Mg.PD04.PBE.UPF`、`O.upf`
- 轨道：`Mg_gga_10au_100Ry_2s1p.orb`、`O_gga_6au_100Ry_2s2p1d.orb`

### 规则2：使用element.json文件

如果文件名不以元素名开头，可在赝势库目录下创建 `element.json` 文件，内容为：

```json
{
  "Mg": "my_magnesium_pseudopotential.upf",
  "O": "my_oxygen_pseudopotential.upf"
}
```

### 规则3：自动设置ecutwfc（可选）

如果在赝势库目录下提供 `ecutwfc.json` 文件，abacustest会自动设置ecutwfc为体系所有元素的最大值：

```json
{
  "Mg": 100,
  "O": 100
}
```

APNS赝势轨道库已经包含了这些配置文件，可以直接使用。

## 1.5 准备工作检查清单

完成上述步骤后，检查以下内容：

- [ ] abacustest已安装，`abacustest --version` 可以正常运行
- [ ] 已下载APNS赝势轨道库，两个文件夹存在
- [ ] 环境变量已设置，`echo $ABACUS_PP_PATH` 显示正确路径
- [ ] 环境变量已设置，`echo $ABACUS_ORB_PATH` 显示正确路径

如果以上都完成，就可以开始准备ABACUS输入文件了。
