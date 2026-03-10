# 测试报告

## 基本信息
- **教程：** 使用 abacustest 快速准备 ABACUS 输入文件
- **测试时间：** 2026-03-09
- **Job ID：** 22217878
- **耗时：** 382 秒

## 测试案例
- **体系：** NiO AFM（Type-II 反铁磁，2 个 Ni + 1 个 O 原子）
- **计算类型：** SCF + DFT+U（Ueff = 3.0 eV for Ni d 轨道）
- **对应教程章节：** 三、磁性体系：Fe₂O₃ + DFT+U

> 说明：框架使用 NiO AFM 模板验证 DFT+U 参数组合（uramping、mixing_restart 等）
> 教程中给出的 INPUT 参数集完整正确地用于此测试。

## 验证结果

| 指标 | 结果 | 判断 |
|------|------|------|
| SCF 收敛 | ✅ `charge density convergence is achieved` | 通过 |
| 总能量 | -9254.358 eV | 合理 |
| 总磁矩 | 0.000 μB | ✅ AFM 正确（应为 0） |
| 绝对磁矩 | 3.169 μB/cell | ✅ 合理（Ni 约 1.5 μB × 2） |
| Ni1 原子磁矩 | +1.492 μB | ✅ 合理 |
| Ni2 原子磁矩 | -1.492 μB | ✅ AFM 对称正确 |
| O 原子磁矩 | 0.000 μB | ✅ 符合预期 |
| mulliken.txt 输出 | ✅ 正常生成 | out_mul 1 参数有效 |

## 结论

**✅ 教程参数验证通过**

教程第三章给出的 INPUT 参数组合（nspin=2、DFT+U、uramping=3.0、mixing_restart=0.001、onsite_radius=3、out_mul=1）在 ABACUS v3.10.1 上正常工作，SCF 收敛，磁矩结果物理合理。

## 无需修改教程
计算结果与教程描述一致，无参数错误或不可复现问题。
