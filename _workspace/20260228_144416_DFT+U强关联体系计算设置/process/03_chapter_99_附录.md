## 附录

### A. DFT+U 参数速查表

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `dft_plus_u` | int | 0 | DFT+U 总开关（1=开，0=关） |
| `orbital_corr` | int 数组 | — | 各 species 施加+U的l量子数（-1=不施加；长度=ntype） |
| `hubbard_u` | float 数组 | — | 各 species 的U值，单位 eV（长度=ntype） |
| `omc` | int | 0 | occupation matrix control（0=标准，1=读入后更新，2=固定） |
| `yukawa_potential` | int | 0 | 自动计算U值开关（1=开） |
| `yukawa_lambda` | float | — | Yukawa screening length（手动设置，可选） |

### B. 常用 grep 命令速查

```bash
# 总能量
grep "FINAL_ETOT" OUT.NiO/running_scf.log

# 能隙（最后一步）
grep "E_bandgap" OUT.NiO/running_scf.log | tail -1

# 总磁矩和绝对磁矩（最后一步）
grep "magnetism" OUT.NiO/running_scf.log | tail -2

# 原子磁矩（Mulliken）
grep "Total Magnetism" OUT.NiO/Mulliken.txt

# occupation matrix（每SCF步）
grep -n "L(S)DA+U" OUT.NiO/running_scf.log
```

### C. 常见问题

**Q：计算报错 "orbital_corr length does not match ntype"**
A：`orbital_corr` 和 `hubbard_u` 的数组长度必须等于 STRU 中 `ATOMIC_SPECIES` 的种类数（ntype）。本例 ntype=3（Ni1/Ni2/O），所以两个数组都需要 3 个数值。

**Q：SCF 不收敛，磁矩振荡**
A：反铁磁体系 SCF 有时不稳定。尝试：减小 `mixing_beta`（0.2~0.3）、设置 `mixing_gg0`（约 1.0）或增大 `scf_nmax`。

**Q：想算铁磁态，如何设置？**
A：在 STRU 中将所有 Ni 原子设为同一 species（或两种 species 但磁矩同号），`orbital_corr` 和 `hubbard_u` 也只需对 Ni 的 species 设置。

**Q：只有一种 Ni 的情况如何设置 orbital_corr？**
A：如果 STRU 中只有 Ni 和 O（ntype=2），则 `orbital_corr 2 -1`，`hubbard_u 5.0 0.0`。

### D. 参考文献

[1] Qu X, Xu P, Jiang H, He L, Ren X. DFT+U within the framework of linear combination of numerical atomic orbitals. *The Journal of Chemical Physics*. 2022;156(23):234104. [doi:10.1063/5.0090122](https://doi.org/10.1063/5.0090122)
