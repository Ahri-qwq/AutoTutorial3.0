# AutoTutorial 3.0 - 测试报告

**测试时间：** 2026-03-02 15:33

**教程文件：** `_workspace/20260302_光学性质计算/07_Final_Tutorial_光学性质计算.md`

**测试目录：** `_workspace/20260302_光学性质计算/test_20260302_151725/`

---

## 测试结果

### ABACUS+PYATB 光学性质测试 - SiO2_optic

**Job ID：** 22166217
**运行时长：** 3459 秒（约 57.7 分钟）
**镜像：** `registry.dp.tech/dptech/prod-19853/abacus-pyatb-open:v0.0.1`
**机型：** c16_m32_cpu（16核32GB）

| 参数 | 预期值 | 实际值 | 误差 | 容差 | 状态 |
|------|--------|--------|------|------|------|
| total_energy | -7835.176 eV | -7832.827 eV | 0.030%（相对）| 0.1% | ✅ PASS |
| occupied_bands | 64 | 64 | 0（绝对）| <0.5 | ✅ PASS |
| fermi_energy | 5.5385382545 eV | 5.6436583054 eV | 1.898%（相对）| 5% | ✅ PASS |
| pyatb_completed | 存在 | 存在 | 布尔 | — | ✅ PASS |

**PYATB 输出文件：**
- `pyatb/Out/Optical_Conductivity/dielectric_function_imag_part.dat` ✅
- `pyatb/Out/Optical_Conductivity/dielectric_function_real_part.dat` ✅
- `pyatb/Out/Optical_Conductivity/optical_conductivity_imag_part.dat` ✅
- `pyatb/Out/Optical_Conductivity/optical_conductivity_real_part.dat` ✅
- 9 个方向分量 PDF 图（df-xx/yy/zz 等）✅

**总体结果：** ✅ 通过

---

## 总结

- 总测试数：1
- 通过数：1
- 失败数：0
- 通过率：100.0%

✅ **所有测试通过！**

---

## 备注

- 总能量偏差 -2.35 eV（相对误差 0.030%），在容差范围内。偏差原因：轨道文件映射（7au→8au），基组略有差异导致能量微小变化。
- 费米能偏差约 0.105 eV（1.898%），在 5% 容差内。LCAO 基组大小对费米能有一定影响。
- PYATB 两步流程完整运行（ABACUS 矩阵输出 → 复制 → PYATB 光学计算），介电函数和光学电导率文件均正常生成。
- 教程参数、文件结构、计算流程完全正确，可复现。
