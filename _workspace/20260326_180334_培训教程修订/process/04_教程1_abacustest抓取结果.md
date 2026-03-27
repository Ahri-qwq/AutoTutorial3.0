# 使用 abacustest 抓取计算结果

ABACUS 计算完成后，结果分散在多个输出文件中。abacustest 提供了便捷的后处理功能，可以快速提取关键结果。

## 基本用法

```bash
abacustest model inputs post -j <job_directory>
```

这个命令会自动识别计算类型（SCF、relax、band、DOS 等），提取相应的关键结果。

---

## 提取 SCF 计算结果

对于自洽计算，abacustest 可以提取：
- 最终总能量
- 收敛步数
- 计算时间
- 费米能级（金属）或带隙（半导体）

**示例：**
```bash
cd Si_scf/
abacustest model inputs post -j .
```

**输出示例：**
```json
{
  "final_energy_eV": -241.234567890,
  "converged": true,
  "scf_steps": 15,
  "total_time_s": 123.45,
  "fermi_energy_eV": 5.678
}
```

---

## 批量提取多个任务结果

对于批量计算任务，可以一次性提取所有结果：

```bash
abacustest model inputs post -j job1/ job2/ job3/
```

或使用通配符：
```bash
abacustest model inputs post -j */
```

结果会汇总到一个 JSON 文件中，便于后续分析。

---

## 提取结构优化结果

对于 relax 或 cell-relax 计算：

```bash
cd Si_relax/
abacustest model inputs post -j .
```

**提取信息：**
- 优化后的晶格参数
- 优化后的原子坐标
- 最终能量
- 最大受力
- 应力张量（cell-relax）

**输出示例：**
```json
{
  "final_energy_eV": -241.234567890,
  "converged": true,
  "lattice_constant_angstrom": 5.431,
  "max_force_eV_per_angstrom": 0.008,
  "stress_kbar": [0.12, 0.12, 0.12, 0.0, 0.0, 0.0]
}
```

优化后的结构保存在 `OUT.*/STRU_ION_D` 文件中。

---

## 提取能带和 DOS 数据

**能带计算：**
```bash
cd Si_band/
abacustest model inputs post -j .
```

提取的数据包括：
- 能带路径
- 各 k 点的能量本征值
- 费米能级位置

**DOS 计算：**
```bash
cd Si_dos/
abacustest model inputs post -j .
```

提取的数据包括：
- 能量网格
- 总态密度（TDOS）
- 投影态密度（PDOS，如果计算了）

---

## 自定义提取规则

abacustest 支持自定义提取规则。创建 `extract_rules.json`：

```json
{
  "extract": {
    "final_energy": {
      "file": "OUT.*/running_scf.log",
      "pattern": "FINAL_ETOT_IS\\s+([\\-\\d\\.]+)\\s+eV",
      "type": "float"
    },
    "band_gap": {
      "file": "OUT.*/running_scf.log",
      "pattern": "E_bandgap\\s+([\\d\\.]+)\\s+eV",
      "type": "float"
    }
  }
}
```

使用自定义规则：
```bash
abacustest model inputs post -j . --rules extract_rules.json
```

---

## 结果可视化

abacustest 还提供了简单的可视化功能（需要安装 matplotlib）：

**绘制收敛曲线：**
```bash
abacustest plot scf -j Si_scf/
```

**绘制能带结构：**
```bash
abacustest plot band -j Si_band/
```

**绘制态密度：**
```bash
abacustest plot dos -j Si_dos/
```

---

## 常见问题

**Q1：找不到输出文件**
- 确保计算已完成
- 检查 `OUT.*` 目录是否存在

**Q2：提取的数据不完整**
- 检查计算是否正常收敛
- 查看 `running_scf.log` 是否有错误信息

**Q3：批量提取时部分任务失败**
- abacustest 会跳过失败的任务，继续处理其他任务
- 检查失败任务的日志文件

---
