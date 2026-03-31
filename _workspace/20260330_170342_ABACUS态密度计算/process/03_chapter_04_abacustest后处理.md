# 第四章：使用 abacustest 后处理

## 4.1 abacustest dos-pdos 命令详解

### 基本语法

```bash
abacustest model dos-pdos -j <job_path> [options]
```

### 主要参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-j, --job` | ABACUS计算目录路径 | 必需 |
| `--range` | 能量范围（相对费米能级），单位eV | -10 10 |
| `--plot-type` | 绘图类型 | species |
| `--atom-index` | 原子编号（1-based） | 所有原子 |
| `--suffix` | 输出文件后缀 | 无 |
| `--no-save-data` | 不保存数据文件 | False |
| `--no-save-plot` | 不保存图片文件 | False |

### 绘图类型说明

| plot-type | 说明 | 示例 |
|-----------|------|------|
| species | 按元素投影 | Mg、O的PDOS |
| shell | 按轨道壳层投影 | Mg的s、p轨道 |
| orbital | 按具体轨道投影 | O的p_x、p_y、p_z |
| atom | 按原子编号投影 | 第1个原子的PDOS |

### 输出文件

| 文件名 | 内容 | 格式 |
|--------|------|------|
| DOS.dat | 总DOS数据 | 文本 |
| DOS.png | 总DOS图 | 图片 |
| PDOS.dat | PDOS数据 | 文本 |
| PDOS.png | PDOS图 | 图片 |

## 4.2 绘制总 DOS 图

### 基本用法

进入计算目录，执行：

```bash
cd <计算目录>
abacustest model dos-pdos -j .
```

这将自动：
1. 读取 `OUT.ABACUS/TDOS` 和 `OUT.ABACUS/PDOS` 文件
2. 从 `running_nscf.log` 提取费米能级
3. 生成 `DOS.png` 和 `PDOS.png`

### 自定义能量范围

聚焦费米能级附近：

```bash
abacustest model dos-pdos -j . --range -5 7
```

这将绘制费米能级以下5 eV到以上7 eV的范围。

### DOS 图的特征

**MgO 的总 DOS 图应显示：**
- 价带区域（负能量）：-6 eV 到 0 eV
- 带隙：约 4-5 eV（PBE泛函低估）
- 导带区域（正能量）：4 eV 以上
- 费米能级处 DOS 为零

## 4.3 绘制 PDOS 图

### 按元素投影

默认按元素投影：

```bash
abacustest model dos-pdos -j . --plot-type species
```

**MgO 的 PDOS 特征：**
- O 的 PDOS 主导价带（-6 eV 到 0 eV）
- Mg 的 PDOS 主导导带（4 eV 以上）
- 两者在价带区域几乎不重叠（离子键特征）

### 按轨道壳层投影

查看不同轨道的贡献：

```bash
abacustest model dos-pdos -j . --plot-type shell
```

**预期结果：**
- O 2p 轨道主导价带顶
- Mg 3s 轨道主导导带底
- O 2s 轨道在价带深处（约 -20 eV）

### 按具体轨道投影

查看 p 轨道的各个分量：

```bash
abacustest model dos-pdos -j . --plot-type orbital
```

这将显示 p_x、p_y、p_z 的独立贡献。

## 4.4 PDOS 图的物理意义

### 离子键的 PDOS 特征

MgO 的 PDOS 图清晰展示了离子键的特征：

**价带区域（占据态）**
- 主要由 O 2p 轨道组成
- Mg 的贡献极小
- 说明价电子主要定域在 O 原子上

**导带区域（未占据态）**
- 主要由 Mg 3s 轨道组成
- O 的贡献较小
- 说明激发态电子倾向于定域在 Mg 原子上

**PDOS 分离明显**
- Mg 和 O 的 PDOS 在价带区域几乎不重叠
- 表明电子转移完全，形成 Mg²⁺ 和 O²⁻
- 这是典型的离子键特征

### 与共价键的对比

如果是共价键材料（如 Si）：
- 价带区域不同原子的 PDOS 会显著重叠
- 表明电子在原子间共享
- PDOS 峰值能量接近

MgO 的 PDOS 分离明显，与共价键形成鲜明对比。

### 轨道贡献分析

通过 PDOS 可以确定：
- 价带顶的主要轨道成分：O 2p
- 导带底的主要轨道成分：Mg 3s
- 带隙跃迁类型：O 2p → Mg 3s

这些信息对理解材料的光学性质非常重要。
