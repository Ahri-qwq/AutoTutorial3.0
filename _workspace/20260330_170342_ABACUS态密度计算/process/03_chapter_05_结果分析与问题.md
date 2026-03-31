# 第五章：结果分析与常见问题

## 5.1 DOS/PDOS 图的解读

### 总 DOS 的分析

**判断材料类型**

通过费米能级处的 DOS 判断：
- DOS(E_F) > 0：金属
- DOS(E_F) = 0 且带隙 < 3 eV：半导体
- DOS(E_F) = 0 且带隙 > 3 eV：绝缘体

MgO 的带隙约 4-5 eV（PBE），属于绝缘体。

**带隙的确定**

从 DOS 图确定带隙：
1. 找到价带顶：费米能级以下 DOS 最后消失的位置
2. 找到导带底：费米能级以上 DOS 开始出现的位置
3. 带隙 = 导带底能量 - 价带顶能量

**注意：** PBE 泛函通常低估带隙，MgO 的实验带隙约 7.8 eV。

### PDOS 的分析

**轨道贡献识别**

通过 PDOS 确定：
- 价带主要由哪些轨道组成
- 导带主要由哪些轨道组成
- 不同元素的相对贡献

**成键性质判断**

| PDOS 特征 | 化学键类型 |
|----------|-----------|
| 不同元素 PDOS 分离 | 离子键 |
| 不同元素 PDOS 重叠 | 共价键 |
| 费米能级处 PDOS 连续 | 金属键 |

### MgO 案例的物理解释

**价带特征**
- O 2p 轨道主导，能量范围 -6 到 0 eV
- 对应 O²⁻ 的满壳层电子结构
- 6 个电子占据 3 个 p 轨道

**导带特征**
- Mg 3s 轨道主导，能量 > 4 eV
- 对应 Mg²⁺ 失去的 2 个电子的空轨道
- 激发态电子会回到 Mg 原子

**离子键形成**
- Mg 失去 2 个 3s 电子 → Mg²⁺
- O 获得 2 个电子填充 2p 轨道 → O²⁻
- PDOS 分离证实了完全的电子转移

## 5.2 常见问题排查

### Q1：PDOS 文件未生成

**症状：**
- `OUT.ABACUS/` 目录下没有 `PDOS` 文件
- 只有 `DOS1` 和 `TDOS` 文件

**可能原因：**
1. `basis_type` 不是 `lcao`
2. `out_dos` 设置为 1 而不是 2

**解决方法：**
```
# 检查 INPUT 文件
basis_type    lcao    # 必须是 lcao
out_dos       2       # 必须是 2
```

### Q2：DOS 曲线不平滑，有尖峰

**症状：**
- DOS 曲线呈锯齿状
- 有明显的尖峰

**可能原因：**
- k 点密度不够

**解决方法：**
1. 增加 k 点密度：
```
# KPT 文件
K_POINTS
0
Gamma
24 24 24 0 0 0    # 从 18×18×18 增加到 24×24×24
```

2. 或使用更小的 kspacing：
```
# INPUT 文件
kspacing    0.06    # 从 0.08 减小到 0.06
```

### Q3：费米能级提取错误

**症状：**
- grep 命令找不到 EFERMI
- 或者提取的值明显不合理

**可能原因：**
1. 查找的文件路径错误
2. 计算未正常完成

**解决方法：**
```bash
# 确认文件存在
ls OUT.ABACUS/running_nscf.log

# 使用正确的路径
grep -i 'efermi' OUT.ABACUS/running_nscf.log

# 或从 SCF 日志提取
grep -i 'efermi' OUT.ABACUS/running_scf.log
```

### Q4：NSCF 计算不收敛

**症状：**
- NSCF 计算报错或不收敛
- 提示找不到电荷密度文件

**可能原因：**
1. SCF 未完成或未输出电荷密度
2. `init_chg` 未设置为 `file`

**解决方法：**
1. 确认 SCF 已完成：
```bash
ls OUT.ABACUS/SPIN1_CHG
```

2. 检查 INPUT 参数：
```
calculation    nscf
init_chg       file    # 必须设置
```

## 5.3 进阶技巧

### 提高计算效率

**SCF 阶段：**
- 使用对称性：`symmetry = 1`
- 使用较稀疏的 k 点（12×12×12）
- 选择高效的求解器：`ks_solver = genelpa`

**NSCF 阶段：**
- 只计算需要的能量范围
- 使用并行计算

### 输出特定能量范围的 DOS

虽然 ABACUS 会输出完整能量范围的 DOS，但可以在后处理时限制范围：

```bash
abacustest model dos-pdos -j . --range -10 15
```

### 对比不同材料的 DOS

计算多个材料后，可以将 DOS 数据导出，使用绘图软件对比：

```bash
# 材料1
abacustest model dos-pdos -j material1/ --suffix mat1

# 材料2
abacustest model dos-pdos -j material2/ --suffix mat2
```

这将生成 `DOS_mat1.dat` 和 `DOS_mat2.dat`，便于对比。
