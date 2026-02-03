# 第四章：其他工具简介

除了abacustest，还有其他工具可以用于准备ABACUS输入文件和进行结构转换。本章简要介绍ASE-ABACUS和ATOMKIT两个工具。

## 4.1 ASE-ABACUS接口

ASE（Atomic Simulation Environment）是丹麦技术大学开发的开源原子模拟Python工具库，广泛应用于计算材料科学。ASE-ABACUS是专门为ABACUS开发的接口，独立于ASE主仓库。

### 4.1.1 安装ASE-ABACUS

ASE-ABACUS需要单独下载安装：

```bash
git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
pip install .
```

**注意：** 通过 `pip install ase` 安装的官方ASE不包含ABACUS接口，必须使用上述方式安装ASE-ABACUS分支。

### 4.1.2 设置环境变量

ASE-ABACUS需要知道赝势和轨道文件的位置。环境变量设置方法参见第一章1.3节。

也可以在Python脚本中设置：

```python
import os
os.environ['ABACUS_PP_PATH'] = '/path/to/pseudopotentials'
os.environ['ABACUS_ORBITAL_PATH'] = '/path/to/orbitals'
```

### 4.1.3 CIF转STRU

使用ASE-ABACUS将CIF文件转换为STRU文件：

```python
from ase.io import read, write
from pathlib import Path

# 读取CIF文件
cif_file = 'MgO.cif'
atoms = read(cif_file, format='cif')

# 指定赝势和轨道文件
pp = {'Mg': 'Mg.PD04.PBE.UPF', 'O': 'O.upf'}
basis = {'Mg': 'Mg_gga_10au_100Ry_2s1p.orb', 'O': 'O_gga_6au_100Ry_2s2p1d.orb'}

# 写入STRU文件
write('STRU', atoms, format='abacus', pp=pp, basis=basis)
```

**说明：**
- `pp` 和 `basis` 字典指定每个元素对应的赝势和轨道文件名
- 文件名需要与环境变量指定的目录中的文件匹配
- 这只生成STRU文件，不会自动生成INPUT和KPT文件

### 4.1.4 STRU转CIF

反向转换也很简单：

```python
from ase.io import read, write

# 读取STRU文件
atoms = read('STRU', format='abacus')

# 写入CIF文件
write('output.cif', atoms, format='cif')
```

### 4.1.5 STRU转POSCAR

转换为VASP的POSCAR格式：

```python
from ase.io import read, write

atoms = read('STRU', format='abacus')
write('POSCAR', atoms, format='vasp')
```

### 4.1.6 ASE-ABACUS的优势

- **灵活性高**：可以在Python脚本中编程控制
- **支持格式多**：ASE支持几十种结构文件格式
- **可调用ABACUS计算**：不仅能转换文件，还能直接调用ABACUS进行计算
- **适合自动化流程**：可以集成到工作流中

### 4.1.7 ASE-ABACUS的局限

- 需要手动指定赝势和轨道文件
- 不会自动生成INPUT和KPT文件
- 需要一定的Python编程基础
- 对于批量处理，需要自己编写脚本

## 4.2 ATOMKIT工具

ATOMKIT是VASPKIT开发团队开发的跨平台建模与结构转换工具，提供命令行交互界面，支持多种材料模拟软件的结构文件格式。

### 4.2.1 下载和安装

从VASPKIT官网下载最新版本：

```bash
wget https://sourceforge.net/projects/vaspkit/files/Binaries/atomkit.0.9.0.linux.x64.tar.gz
tar -zxvf atomkit.0.9.0.linux.x64.tar.gz
cd atomkit.0.9.0.linux.x64
bash setup.sh
```

安装完成后，可以直接运行 `atomkit` 命令。

### 4.2.2 使用ATOMKIT转换结构

ATOMKIT提供交互式界面，使用方式如下：

**启动ATOMKIT：**

```bash
atomkit
```

**或者使用脚本模式：**

创建一个输入脚本 `convert.txt`：

```
1        # 选择功能：结构转换
113      # 选择输入格式：CIF
101      # 选择输出格式：ABACUS STRU
MgO.cif  # 输入文件名
```

然后执行：

```bash
atomkit < convert.txt
```

### 4.2.3 ATOMKIT的功能代码

ATOMKIT的主要功能包括：

**结构转换（功能1）：**
- 支持格式：CIF、POSCAR、XYZ、ABACUS STRU等
- 输入格式代码：
  - 113：CIF
  - 101：VASP POSCAR
  - 102：ABACUS STRU
- 输出格式代码：
  - 101：ABACUS STRU
  - 102：VASP POSCAR
  - 103：CIF

**其他功能：**
- 晶体结构建模
- 结构可视化
- 晶格参数调整
- 原子坐标变换

### 4.2.4 ATOMKIT的优势

- **无需编程**：命令行交互，易于使用
- **功能丰富**：不仅转换格式，还能建模和可视化
- **跨平台**：支持Linux、macOS、Windows
- **快速迭代**：开发活跃，功能持续更新

### 4.2.5 ATOMKIT的局限

- 交互式操作，不适合大规模批量处理
- 不会自动配置赝势和轨道文件
- 不会生成INPUT和KPT文件
- 需要手动输入功能代码，学习曲线稍陡

## 4.3 三种工具对比

| 特性 | abacustest | ASE-ABACUS | ATOMKIT |
|------|-----------|------------|---------|
| **安装方式** | pip安装 | git+pip | 下载解压 |
| **使用方式** | 命令行 | Python脚本 | 交互式 |
| **结构转换** | ✓ | ✓ | ✓ |
| **自动配置赝势轨道** | ✓ | ✗ | ✗ |
| **生成INPUT文件** | ✓ | ✗ | ✗ |
| **生成KPT文件** | ✓ | ✗ | ✗ |
| **批量处理** | ✓ | ✓（需编程） | ✗ |
| **磁性材料设置** | ✓ | ✗ | ✗ |
| **DFT+U设置** | ✓ | ✗ | ✗ |
| **可视化** | ✗ | ✓（需额外工具） | ✓ |
| **编程基础要求** | 低 | 高 | 低 |

## 4.4 工具选择建议

**推荐使用abacustest，如果：**
- 需要快速准备完整的ABACUS输入文件夹
- 需要批量处理多个结构
- 需要自动配置赝势和轨道
- 需要设置磁性材料或DFT+U参数

**推荐使用ASE-ABACUS，如果：**
- 需要在Python脚本中自动化处理
- 需要集成到复杂的工作流中
- 需要调用ABACUS进行计算
- 熟悉Python编程

**推荐使用ATOMKIT，如果：**
- 只需要简单的格式转换
- 需要可视化结构
- 需要建模和调整晶格参数
- 不想编写脚本

**组合使用：**

实际工作中，可以组合使用多个工具：
1. 用ATOMKIT可视化和调整结构
2. 用abacustest快速准备输入文件
3. 用ASE-ABACUS进行后处理和分析

## 4.5 其他工具小结

本章介绍了ASE-ABACUS和ATOMKIT两个工具，它们各有特点：

- **ASE-ABACUS**：适合Python用户，灵活性高，可编程控制
- **ATOMKIT**：适合交互式使用，功能丰富，易于上手

选择合适的工具可以提高工作效率。对于大多数用户，abacustest是最便捷的选择。
