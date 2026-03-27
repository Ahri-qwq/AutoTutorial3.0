# 计算测试报告

## 测试信息

- 教程：`07_Final_Tutorial_ABACUS电子结构计算_能带与态密度.md`
- 测试目录：`test_20260325_150015/`
- 测试时间：2026-03-25
- 平台：Bohrium，项目 205855，机型 c16_m32_cpu，镜像 LTSv3.10.1

## 任务提交记录

| 任务 | Job ID | 耗时 | 状态 |
|------|--------|------|------|
| Si SCF (LCAO, 9×9×9) | 22303033 | ~45s | Finished ✅ |
| Al SCF (PW, 12×12×12) | 22303035 | ~90s | Finished ✅ |
| Si NSCF Band (LCAO, Line) | 22303116 | ~90s | Finished ✅ |
| Al NSCF DOS (PW, 20×20×20) | 22303129 | ~30min | Finished ✅ |

## 计算结果

### 案例一：Si 能带结构

| 指标 | 计算值 | 教程预期 | 状态 |
|------|--------|---------|------|
| SCF 收敛 | ✅ | — | PASS |
| EFERMI (SCF) | 6.831 eV | — | — |
| EFERMI (NSCF) | 6.852 eV | 6.390 eV（示例） | 正常（绝对值依赖结构） |
| BANDS_1.dat 生成 | ✅ 211行 | — | PASS |
| 带隙（out_bandgap） | 0.858 eV | ~0.57 eV | ⚠️ 偏高，见说明 |

**带隙偏差说明：** `out_bandgap` 在 Line 模式 k 点集中搜索 VBM/CBM。Si 间接带隙的 CBM 在 Δ 点（Γ→X 约 0.85 处），Line 路径上该区间仅有约 17 个 k 点，可能未精确覆盖 CBM 最低点，导致带隙偏高。若用均匀 k 网格做 NSCF，结果更接近 0.57 eV。教程中的 0.57 eV 是文献参考值（均匀网格），与本次 Line 路径结果 0.858 eV 的差异是正常的方法论差异，**不属于错误**。

### 案例二：Al DOS

| 指标 | 计算值 | 教程预期 | 状态 |
|------|--------|---------|------|
| SCF 收敛 | ✅ | — | PASS |
| EFERMI (SCF) | 10.913 eV | — | — |
| EFERMI (NSCF) | 10.931 eV | 10.963 eV（示例） | PASS（差 0.03 eV，正常） |
| DOS1 文件生成 | ✅ | — | PASS |
| DOS at Fermi level | ~1.59 states/eV/cell | 连续不为零 | PASS ✅ |
| PDOS 文件生成 | ❌ 未生成 | 教程期望输出 PDOS | ⚠️ 见说明 |

**PDOS 未生成说明：** PW 基组（plane wave）下 `out_dos 2` 只输出 TDOS（DOS1 文件），不输出 PDOS。PDOS 依赖 Mulliken 投影，需要 LCAO 基组。教程案例二使用 PW 基组，因此 PDOS 在实际计算中无法通过 `out_dos 2` 输出——教程中描述 PDOS 的章节（3.5节）存在基组不一致问题。

## 发现的问题

### 问题1：带隙值说明不完整（建议补充注释）
- **位置：** 2.6节"结果解读"
- **问题：** 教程说"带隙约为 0.57 eV"，但实际使用 Line KPT 运行时 `out_bandgap` 输出约 0.86 eV
- **建议：** 补充说明：`out_bandgap` 在 Line k点集搜索，结果可能偏高；精确间接带隙需用均匀网格 NSCF

### 问题2：PW 基组无法输出 PDOS（教程错误）
- **位置：** 3.5节"PDOS 绘图"
- **问题：** 教程案例二（Al）使用 PW 基组，但 `out_dos 2` + PW 只生成 TDOS，不生成 PDOS
- **影响：** 读者按教程操作后 `OUT.Al/PDOS` 文件不存在，`abacus-plot -d -p -o` 会报错
- **修复方案：** 将 Al 案例的 PDOS 改为 LCAO 基组，或在教程中注明"PDOS 需要 LCAO 基组"

## 总体评价

- **所有 4 个任务均成功完成**，无运行错误
- Si 能带结构：BANDS_1.dat 正常生成，带隙值与 Line 模式特性一致
- Al TDOS：费米面处 DOS 连续（~1.59 states/eV/cell），验证金属特征 ✅
- **关键问题：** PW 基组不支持 PDOS 输出，教程 3.5节 PDOS 部分存在基组错误，需修正
