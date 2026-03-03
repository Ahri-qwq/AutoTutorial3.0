# 教程测试报告

**生成时间：** 2026-03-03 17:15
**测试教程：** `_workspace/20260303_161406_ABACUS+Phonopy声子谱计算/07_Final_Tutorial_使用ABACUS+Phonopy计算声子谱.md`
**测试目录：** `_workspace/20260303_161406_ABACUS+Phonopy声子谱计算/test_20260303_165531/`

---

## 测试概要

- ✅ **测试状态：** 通过
- 📊 **案例覆盖率：** 1/1（Al_FCC ✅）
- ⏱️ **总耗时：** 约 8 分钟（含任务排队）

---

## 案例 1/1：Al_FCC（FCC 铝，SCF+原子受力计算）

> **验证范围说明：** 本插件测试 ABACUS+Phonopy 工作流中 ABACUS 负责的部分——
> 即 SCF 自洽场计算 + 原子受力输出（`cal_force = 1`）。
> Phonopy 后处理（force constant 提取、声子谱绘制）为外部软件，不在验证范围内。

### 计算步骤

| 步骤 | Job ID | 状态 | 耗时 |
|------|--------|------|------|
| Al_FCC SCF+力 | 22182083 | ✅ Finished | 51 秒 |

### 结果对比

| 验证项 | 期望 | 实际值 | 状态 |
|--------|------|--------|------|
| SCF 是否收敛 | 是 | 是（第 603 行："charge density convergence is achieved"） | ✅ PASS |
| cal_force=1 是否输出了原子受力 | 是 | 是（第 952 行：`TOTAL-FORCE (eV/Angstrom)` 块存在） | ✅ PASS |

**补充信息（仅记录，不作为 PASS/FAIL 判据）：**
- 总能量：-7531.3568 eV（FCC Al 4原子常规胞，含 ONCV 半芯 2s2p 态，绝对值大为正常）
- 原子受力：Al1~Al4 均为 0.000 eV/Å（未位移对称结构，力为零符合预期）
- 应力张量：5.1303 kbar（对角元相等，各向同性，符合 FCC 立方对称性）

**总体结果：✅ 通过（2/2 验证项）**

---

## 结论

✅ **教程测试通过**，案例覆盖率 1/1。

ABACUS SCF+力 计算正常运行并收敛，`cal_force = 1` 参数生效，原子受力正确输出。
教程中的 INPUT/STRU/KPT 文件格式正确，参数设置有效。

**轨道文件修正说明：**
教程中使用的 `Al_gga_7au_100Ry_4s4p1d.orb` 在 ABACUS GitHub 仓库中不存在。
测试时自动使用 `Al_gga_8au_100Ry_4s4p1d.orb`（截断半径 8 au，命名规律一致），计算通过。
教程原文已在 Step 7 中修正此文件名。

---

## 测试详情

### 任务信息
- Job 22182083: Al_FCC_phonopy_scf（c16_m32_cpu，51 秒）

### 环境信息
- Python: 3.12.9 / Bohrium CLI: 1.1.0
- ABACUS 镜像: registry.dp.tech/dptech/abacus:LTSv3.10.1

---

**测试框架版本：** AutoTutorial 3.0
**生成工具：** test_framework_integrated.py（PhonopyPlugin v1.0）
