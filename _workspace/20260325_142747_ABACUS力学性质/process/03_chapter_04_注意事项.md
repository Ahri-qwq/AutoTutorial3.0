# 第四章：注意事项与参数速查

## 4.1 结构优化参数速查

**ionic relax（固定晶胞）最简 INPUT：**

```
INPUT_PARAMETERS
suffix                  example
ntype                   1
pseudo_dir              ./
calculation             relax
basis_type              lcao
ecutwfc                 100
scf_thr                 1e-6
cal_force               1
force_thr_ev            0.01
relax_nmax              100
out_stru                1
kspacing                0.1
```

**cell-relax（变胞优化）最简 INPUT：**

```
INPUT_PARAMETERS
suffix                  example
ntype                   1
pseudo_dir              ./
orbital_dir             ./
calculation             cell-relax
basis_type              lcao
ecutwfc                 100
scf_thr                 1e-7
cal_force               1
cal_stress              1
force_thr_ev            0.01
stress_thr              0.5
relax_nmax              100
out_stru                1
kspacing                0.08
```

---

## 4.2 常见问题汇总

| 问题 | 原因 | 解决方法 |
|------|------|---------|
| 结构优化不收敛 | 初始结构离平衡位置太远，或 K 点/截断能不足 | 检查初始结构合理性；增大 K 点密度或 ecutwfc |
| 应力收敛但力不收敛 | 原子被约束（STRU 中 0 0 0 约束标志） | 检查 STRU 文件中原子自由度设置 |
| `elastic prepare` 后重复运行丢失数据 | 命令会直接删除已有形变文件夹 | 提交计算前只运行一次 prepare |
| 弹性张量非对角元不为零 | K 点密度或截断能不足导致应力精度差 | 增大 K 点密度（减小 kspacing）或 ecutwfc |
| 与 Materials Project 结果差异大 | 赝势、泛函或截断能不同；取向不一致 | 确认取向与 IEEE 标准一致；对比赝势选择 |

---

## 4.3 计算精度建议

- **K 点密度**：弹性常数对 K 点密度敏感，建议 `kspacing` ≤ 0.1 Å⁻¹，必要时做收敛测试
- **截断能**：建议参考所用赝势的推荐截断能，不要低于赝势推荐值
- **SCF 收敛**：弹性计算的 SCF 阈值应严格（`scf_thr 1e-7` 或更严格），否则应力误差会传导到弹性张量
- **结构优化充分度**：cell-relax 的 `stress_thr` 建议 ≤ 1 kBar；若残余应力过大，拟合结果可靠性降低
