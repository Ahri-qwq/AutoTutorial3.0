# Step 5 案例审查报告

## 5.1 案例完整性检查

对照 process/01_research.md 中的案例完整性检查清单：

### STRU 文件要素
- [x] ATOMIC_SPECIES（Ni1/Ni2/O 的赝势文件名）：`Ni_ONCV_PBE-1.0.upf` / `O_ONCV_PBE-1.0.upf` ✓
- [x] 不同磁矩的Ni定义为不同 species 的原因：已在2.1节解释 ✓
- [x] ATOMIC_POSITIONS 中的磁矩设置（2.0 / -2.0）：已展示 ✓

### INPUT 文件要素
- [x] dft_plus_u = 1 ✓
- [x] orbital_corr = 2 2 -1（l量子数含义已说明）✓
- [x] hubbard_u = 5.0 5.0 0.0 ✓
- [x] nspin = 2（已在INPUT中设置）✓
- [x] out_bandgap = 1（已在INPUT中设置）✓
- [x] out_mul = 1（已在INPUT中设置）✓
- [x] out_chg = 1（已在INPUT中设置）✓
- [x] omc = 2（在第5章演示中设置）✓

### 计算结果要素
- [x] 总能量：-9255.7279034240546025 eV ✓
- [x] 能隙：+0.205369322748 / +2.794192983776 eV ✓
- [x] 总磁矩：0.0；绝对磁矩：3.35321634 ✓
- [x] 原子磁矩：Ni1 +1.8268646，Ni2 -1.8268646 ✓

### omc 演示要素
- [x] onsite.dm 文件格式说明 ✓
- [x] 将 onsite.dm 复制为 initial_onsite.dm ✓
- [x] 设置 omc = 2 后重新计算 ✓
- [x] 验证结果相同（总能量相同） ✓
- [x] occupation matrix 不变（文字说明了） ✓

## 5.2 案例数据准确性检查

所有直接引用的数值均与案例文件核对无误：
- 赝势文件名：与案例完全一致
- DFT+U参数：与案例完全一致
- 计算结果数值：与案例完全一致

## 5.3 待修正项

1. **onsite.dm 的格式**：案例描述了格式，教程第4.4节的occupation matrix格式说明与案例一致。✓
2. **omc=2验证**：案例提到"结果与上一步SCF完全相同"，教程第5.2节有体现。✓

## 结论

案例完整性和准确性通过。无需进一步修改。
