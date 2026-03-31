# 第六章：参数优化与常见问题

## 6.1 常见问题与解决

### 问题1：找不到电荷密度文件

**错误信息**：
```
Can't find the file SPIN1_CHG.cube
```

**原因**：
- SCF未完成或失败
- 电荷密度文件被移动或删除
- `out_chg`未设置为1

**解决方法**：
- 检查SCF是否成功完成
- 确认`OUT.Si/SPIN1_CHG.cube`存在
- 如果文件在其他目录，在INPUT中设置`read_file_dir`

### 问题2：symmetry设置错误

**错误信息**：
```
Error in k-point symmetry operation
```

**原因**：
- NSCF使用Line模式但`symmetry = 1`

**解决方法**：
- 在INPUT_nscf中设置`symmetry = 0`

### 问题3：费米能级单位转换

**问题**：
- 旧版本ABACUS输出费米能级单位为Rydberg

**解决方法**：
- 检查输出单位
- 如果是Ry，乘以13.6058转换为eV
- 新版本（v3.0+）直接输出eV

## 6.2 进阶方向

**杂化泛函计算**：
- 使用HSE06改善带隙精度
- 计算量大，需要更多资源

**自旋极化体系**：
- 磁性材料需要设置`nspin = 2`
- 输出BANDS_1.dat（自旋向上）和BANDS_2.dat（自旋向下）

**投影能带（fat band）**：
- 显示不同轨道对能带的贡献
- 需要使用LCAO基组

第六章完成。

