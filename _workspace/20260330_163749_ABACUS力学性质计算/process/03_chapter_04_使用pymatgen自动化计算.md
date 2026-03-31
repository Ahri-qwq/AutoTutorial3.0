# 第四章：使用 pymatgen 自动化计算

手动生成 24 个应变构型并逐个计算比较繁琐。pymatgen 提供了自动化工具，可以简化整个流程。

## 4.1 pymatgen 简介

### 什么是 pymatgen

pymatgen（Python Materials Genomics）是一个用于材料科学计算的 Python 库，由加州大学圣地亚哥分校开发。

**主要功能：**
- 晶体结构操作和分析
- 与 Materials Project 数据库交互
- 高通量计算工作流
- 弹性常数计算

**官方网站：** https://pymatgen.org/

### 安装 pymatgen

```bash
pip install pymatgen dpdata monty numpy
```

**依赖库说明：**
- `pymatgen`：核心库
- `dpdata`：数据格式转换
- `monty`：工具函数
- `numpy`：数值计算

## 4.2 生成应变构型

### 4.2.1 准备工作

在工作目录下准备：
- `relax/`：包含优化后的 STRU 文件（STRU_ION_D）
- `INPUT`：relax 计算的 INPUT 文件
- `KPT`：k 点文件
- 赝势文件：Al.upf

### 4.2.2 使用 gene_dfm.py 脚本

执行以下命令生成应变构型：

```bash
python gene_dfm.py abacus
```

**脚本功能：**
1. 读取 `relax/STRU_ION_D` 文件
2. 生成 6 种应变模式 × 4 种应变大小 = 24 个构型
3. 创建 `task.000` 到 `task.023` 文件夹
4. 每个文件夹包含：INPUT、KPT、STRU、strain.json

**应变大小设置：**

脚本中默认的应变大小（可修改）：

```python
norm_strains = [-0.010, -0.005, 0.005, 0.010]
shear_strains = [-0.010, -0.005, 0.005, 0.010]
```

### 4.2.3 检查生成的文件

进入任意 task 文件夹查看：

```bash
cd task.000
ls
```

应该看到：INPUT、KPT、STRU、strain.json

查看 strain.json 内容：

```bash
cat strain.json
```

示例输出：

```json
{"strain": [0.01, 0.0, 0.0, 0.0, 0.0, 0.0]}
```

表示在 x 方向施加 1% 的正应变。

## 4.3 批量计算应力

### 4.3.1 使用 shell 脚本批量提交

创建 `run_all.sh` 脚本：

```bash
#!/bin/bash
for i in {000..023}; do
    cd task.$i
    echo "Running task.$i"
    OMP_NUM_THREADS=1 mpirun -np 8 abacus > log 2>&1
    cd ..
done
```

运行脚本：

```bash
bash run_all.sh
```

### 4.3.2 检查计算完成

确认所有 task 都已完成：

```bash
ls task.*/OUT.ABACUS/running_relax.log | wc -l
```

应该输出 24。

## 4.4 计算弹性常数

### 4.4.1 使用 compute_dfm.py 脚本

所有计算完成后，执行：

```bash
python compute_dfm.py abacus
```

**脚本功能：**
1. 读取所有 task 的应力数据
2. 读取对应的应变数据（strain.json）
3. 对每种应变状态进行线性拟合
4. 计算弹性常数矩阵和弹性模量

### 4.4.2 输出结果

屏幕输出示例：

```
# Elastic Constants in GPa
114.23  61.85  61.85   0.00   0.00   0.00
 61.85 114.23  61.85   0.00   0.00   0.00
 61.85  61.85 114.23   0.00   0.00   0.00
  0.00   0.00   0.00  31.67   0.00   0.00
  0.00   0.00   0.00   0.00  31.67   0.00
  0.00   0.00   0.00   0.00   0.00  31.67
# Bulk Modulus BV = 79.31 GPa
# Shear Modulus GV = 28.12 GPa
# Youngs Modulus EV = 72.45 GPa
# Poission Ratio uV = 0.34
```

**结果说明：**
- 第 2-7 行：弹性常数矩阵（6×6）
- Bulk Modulus：体弹模量
- Shear Modulus：剪切模量
- Youngs Modulus：杨氏模量
- Poission Ratio：泊松比

### 4.4.3 保存详细结果

精度更高的结果保存在 `elastic.json` 文件中：

```bash
cat elastic.json
```

包含完整的弹性常数矩阵和所有弹性模量。

## 4.5 完整流程总结

使用 pymatgen 计算弹性常数的完整流程：

**步骤 1：结构优化**
```bash
# 准备 INPUT（cell-relax）、STRU、KPT
mkdir relax
cd relax
# 运行 ABACUS
OMP_NUM_THREADS=1 mpirun -np 8 abacus
cd ..
```

**步骤 2：生成应变构型**
```bash
python gene_dfm.py abacus
```

**步骤 3：批量计算应力**
```bash
bash run_all.sh
```

**步骤 4：计算弹性常数**
```bash
python compute_dfm.py abacus
```

**输出文件：**
- 屏幕输出：弹性常数矩阵和弹性模量
- `elastic.json`：详细结果（JSON 格式）

