# 第七章：使用 abacustest 准备输入与提取结果

abacustest 是 ABACUS 的前后处理工具，可以快速准备输入文件、批量提交任务、提取计算结果。

## 7.1 abacustest 简介

**功能定位：**
- 从 CIF/POSCAR 等格式生成 ABACUS 输入文件
- 批量设置计算参数
- 提取和可视化计算结果
- 高通量计算任务管理

**安装：**

```bash
# 方法1：通过 pip 安装
pip install abacustest

# 方法2：从源码安装
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

**验证安装：**

```bash
abacustest -h
```

## 7.2 准备输入文件

### 7.2.1 基本用法

```bash
abacustest model inputs prepare -h
```

### 7.2.2 从 CIF 文件生成 STRU

假设有一个 Si.cif 文件：

```bash
abacustest model inputs prepare -f Si.cif --ftype cif
```

会在当前目录生成：
- `STRU`：结构文件
- `INPUT`：基本参数文件
- `KPT`：k 点文件

### 7.2.3 自定义参数

```bash
abacustest model inputs prepare \
  -f Si.cif \
  --ftype cif \
  --input "calculation=scf,basis_type=pw,ecutwfc=50" \
  --kpt "4 4 4 0 0 0"
```

参数说明：
- `--input`：设置 INPUT 参数（逗号分隔）
- `--kpt`：设置 k 点网格
- `--pp`：指定赝势文件路径
- `--orb`：指定轨道文件路径（LCAO）

## 7.3 提取计算结果

### 7.3.1 提取总能量

计算完成后，可以用 abacustest 提取关键结果：

```bash
abacustest model post -j ./
```

会自动识别 `OUT.suffix/` 文件夹，提取总能量、费米能级等信息。

### 7.3.2 提取结构信息

对于 relax/cell-relax 计算，提取优化后的结构：

```bash
abacustest model post -j ./ --extract-stru
```

会生成优化后的 STRU 文件。

## 7.4 实战示例：从 CIF 到计算结果

完整流程演示：

**步骤1：下载 CIF 文件**

从 Materials Project 下载 Si 的 CIF 文件（mp-149.cif）。

**步骤2：生成输入文件**

```bash
abacustest model inputs prepare \
  -f mp-149.cif \
  --ftype cif \
  --input "calculation=scf,basis_type=pw,ecutwfc=50,scf_thr=1e-6" \
  --kpt "4 4 4 0 0 0"
```

**步骤3：运行 ABACUS**

```bash
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

**步骤4：提取结果**

```bash
abacustest model post -j ./
```

输出示例：
```
Total Energy: -123.456 eV
Fermi Energy: 5.432 eV
```

## 7.5 进阶功能预告

abacustest 还支持更多高级功能，将在后续教程中介绍：
- DOS/PDOS 绘图（`abacustest model dos-pdos`）
- 能带绘图（`abacustest model band`）
- 弹性常数计算（`abacustest model elastic`）

> **提示：** abacustest 功能丰富，完整文档见 GitHub：https://github.com/pxlxingliang/abacus-test
