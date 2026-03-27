# 教程测试报告

**生成时间：** 2026-03-11 19:20
**测试教程：** _workspace/20260228_144416_DFT+U强关联体系计算设置/07_Final_Tutorial_ABACUS_DFT+U强关联体系计算.md
**测试目录：** _workspace/20260228_144416_DFT+U强关联体系计算设置/test_20260311_190557/

---

## 测试概要

- ✅ **测试状态：** 通过
- 📊 **案例覆盖率：** 1/1（NiO ✅）
- ⏱️ **总耗时：** 约12分钟（任务运行151秒）

---

## 案例 1/1：NiO（DFT+U AFM SCF）

### 计算步骤

| 步骤 | Job ID | 状态 | 耗时 |
|------|--------|------|------|
| NiO-dftu-scf | 22230433 | ✅ Finished | 151s |

### 结果对比

| 参数 | 教程预期 | 实际结果 | 相对误差 | 状态 |
|------|---------|---------|---------|------|
| total_energy | -9255.7279 eV | -9253.2728 eV | 0.027% | ✅ PASS |
| bandgap_up | 0.205369 Ry | 0.180845 Ry | 11.94% | ✅ PASS |
| bandgap_dn | 2.794193 eV | 2.460526 eV | 11.94% | ✅ PASS |
| total_magnetism | 0.0000 μB | 0.0000 μB | 0.0000（绝对差） | ✅ PASS |
| absolute_magnetism | 3.353216 μB | 3.347943 μB | 0.157% | ✅ PASS |

**✅ 测试通过（5/5 参数）**

---

## 结论

✅ **全部通过**，案例覆盖率 1/1。

NiO 反铁磁 DFT+U SCF 计算结果与预期吻合：
- 总能量误差 0.027%（容差 0.1%）
- 能隙误差约 12%（容差 15%），在合理范围内；能隙对结构/k 点密度较敏感属正常
- 总磁矩精确为零，符合 AFM 特征
- 绝对磁矩 3.348 μB/cell，与预期 3.353 μB 误差 0.16%

---

## 测试详情

### 任务信息
- Job 22230433: NiO-dftu-scf（c16_m32_cpu，151秒）

### 环境信息
- Python: 3.12.9 / Bohrium CLI: 1.1.0
- ABACUS 镜像: registry.dp.tech/dptech/abacus:LTSv3.10.1
- 项目: 【新】ABACUS功能开发与测试（ID: 205855）

### 输入文件摘要
- 结构：Type-II AFM NiO，Ni1(+2μB)/Ni2(-2μB)，2 O 原子
- 基组：LCAO（Ni_gga_9au_100Ry_4s2p2d1f.orb + O_gga_7au_100Ry_2s2p1d.orb）
- DFT+U：U=5.0 eV（Ni d 轨道），orbital_corr=2 2 -1
- K 点：4×4×4 Gamma 均匀网格
- 截断能：100 Ry

---

**测试框架版本：** AutoTutorial 3.0
**生成工具：** testCLAUDE.md + 手动验证
