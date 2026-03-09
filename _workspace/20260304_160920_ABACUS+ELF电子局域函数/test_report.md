# 教程测试报告

**生成时间：** 2026-03-09 11:45
**测试教程：** `_workspace/20260304_160920_ABACUS+ELF电子局域函数/07_Final_Tutorial_ABACUS使用教程-ELF电子局域函数计算与可视化.md`
**测试目录：** `_workspace/20260304_160920_ABACUS+ELF电子局域函数/test_20260309_110541/`

---

## 测试概要

- ✅ **测试状态：** 通过
- 📊 **案例覆盖率：** 3/3（H2O_PW ✅、H2O_LCAO ✅、Fe_BCC ✅）
- ⏱️ **总耗时：** 约 25 分钟（含插件开发 ~15 分钟 + 计算 ~128s + 下载）

---

## 案例 1/3：H2O_PW（水分子 PW 基组）

### 计算步骤

| 步骤 | Job ID | 状态 | 耗时 |
|------|--------|------|------|
| scf（PW, nspin=1） | 22216368 | ✅ Finished | 128s |

### 结果对比

| 检查项 | 预期 | 实际 | 状态 |
|--------|------|------|------|
| ELF.cube 存在 | 是 | 是（57903 KB） | ✅ PASS |
| nspin | 1 | 1 | ✅ PASS |
| FFT 格点 | 180³ | 180³ | ✅ PASS |
| 格点步长 | ~0.156 Bohr | 0.155587 Bohr | ✅ PASS |
| 原子信息 | H/H/O（序数 1/1/8） | 1/1/8 | ✅ PASS |
| SCF 收敛 | 是 | 是 | ✅ PASS |

---

## 案例 2/3：H2O_LCAO（水分子 LCAO 基组）

### 计算步骤

| 步骤 | Job ID | 状态 | 耗时 |
|------|--------|------|------|
| scf（LCAO, nspin=1） | 22216369 | ✅ Finished | 72s |

### 结果对比

| 检查项 | 预期 | 实际 | 状态 |
|--------|------|------|------|
| ELF.cube 存在 | 是 | 是（57903 KB） | ✅ PASS |
| nspin | 1 | 1 | ✅ PASS |
| FFT 格点 | 与 PW 一致 | 180³（相同） | ✅ PASS |
| SCF 收敛 | 是 | 是 | ✅ PASS |

### 轨道文件修正

`H_gga_6au_100Ry_2s1p.orb` 不存在（GitHub 无此文件），使用 `H_gga_8au_100Ry_2s1p.orb` 替代（100Ry，2s1p 基组相同，截断半径 8au 为标准版本）。计算正常收敛，教程本身已注明此替代方案。

---

## 案例 3/3：Fe_BCC（BCC 铁自旋极化）

### 计算步骤

| 步骤 | Job ID | 状态 | 耗时 |
|------|--------|------|------|
| scf（PW, nspin=2） | 22216370 | ✅ Finished | 64s |

### 结果对比

| 检查项 | 预期 | 实际 | 状态 |
|--------|------|------|------|
| ELF.cube 存在 | 是 | 是（463 KB） | ✅ PASS |
| ELF_SPIN1.cube 存在 | 是（nspin=2） | 是（463 KB） | ✅ PASS |
| ELF_SPIN2.cube 存在 | 是（nspin=2） | 是（463 KB） | ✅ PASS |
| nspin | 2 | 2 | ✅ PASS |
| SCF 收敛 | 是 | 是 | ✅ PASS |
| 磁矩/atom | ~2.2 μB | 2.30 μB（实验 ~2.2 μB） | ✅ PASS |
| 总能量 | 铁磁态 | -6441.21 eV（铁磁收敛） | ✅ PASS |

---

## 轨道文件修正记录

本次测试发现并修正了以下轨道文件名错误（计算已验证通过）：

| 错误文件名 | 正确文件名 | 来源 |
|-----------|-----------|------|
| `H_gga_6au_100Ry_2s1p.orb` | `H_gga_8au_100Ry_2s1p.orb` | GitHub 候选列表 + 教程注释 |

---

## 新增插件记录

本次测试新建了 ELF 测试插件：
- **插件文件：** `tools/test_plugins/elf_plugin.py`
- **计算类型：** `elf`
- **can_handle 关键词：** `out_elf\s+1`、`ELF.cube`、`电子局域函数`
- **支持案例：** H2O_PW、H2O_LCAO、Fe_BCC（多案例单插件设计）

---

## 结论

✅ **测试通过**，案例覆盖率 3/3。

ELF 计算流程正确：三种计算场景（PW 单自旋、LCAO 单自旋、PW 自旋极化）均能正常输出 ELF cube 文件，SCF 全部收敛。Fe BCC 磁矩 2.30 μB/atom 与实验值吻合。教程内容准确可复现。

---

## 测试详情

### 任务信息
- Job 22216368: H2O_PW-scf（c16_m32_cpu，128s）
- Job 22216369: H2O_LCAO-scf（c16_m32_cpu，72s）
- Job 22216370: Fe_BCC-scf（c16_m32_cpu，64s）

### 环境信息
- Bohrium 项目：205855（【新】ABACUS功能开发与测试）
- ABACUS 镜像：registry.dp.tech/dptech/abacus:LTSv3.10.1
- 机型：c16_m32_cpu

---

**测试框架版本：** AutoTutorial 3.0
**生成工具：** test_framework_integrated.py + elf_plugin.py
