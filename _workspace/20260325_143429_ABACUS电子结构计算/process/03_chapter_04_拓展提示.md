# 第四章：拓展提示

## 四、拓展提示

### 4.1 磁性体系的 DOS 与能带（Fe，nspin=2）

对于铁磁或反铁磁体系，需要在 INPUT 中增加 `nspin 2`，并在 STRU 的原子行设置非零初始磁矩：

```
# STRU 中 Fe 原子行（设置初始磁矩 4.0 μB）
Fe  #label
4.0 #magnetism（初始磁矩）
1   #number of atoms
0.0  0.0  0.0  m  0  0  0
```

`nspin 2` 时，NSCF 计算会输出两套结果：
- 能带：`BANDS_1.dat`（自旋向上）和 `BANDS_2.dat`（自旋向下）
- DOS：`DOS1`（自旋向上）和 `DOS2`（自旋向下）

两套 DOS 共享同一费米能级，绘图时通常将自旋向下的 DOS 取负值，叠加显示。

bcc Fe 的计算参考：
- 晶格常数：2.866 Å
- 磁矩：约 2.2 μB（收敛后）
- smearing_sigma 建议：0.002 Ry（磁性体系需要较小展宽以准确描述自旋极化）

完整 Al/Fe 算例可从官方仓库下载：
```bash
git clone https://gitee.com/mcresearch/abacus-user-guide.git
cd abacus-user-guide/examples/dos_band
# Al/ 和 Fe/ 两个目录
```

### 4.2 带隙修正：HSE06 杂化泛函

PBE 系统性低估带隙（Si: 0.57 eV vs 实验 1.17 eV）。若需要接近实验值的带隙，可使用 HSE06 杂化泛函：

在 INPUT 中添加：
```
dft_functional  hse
hse_omega       0.11    # 屏蔽参数，单位 Bohr^-1，HSE06 默认值 0.11
```

注意 HSE06 需要 LCAO 基组，且计算量显著大于 PBE（通常 5–10 倍）。ABACUS 的 HSE06 依赖 LibRI 库，需确认编译环境支持。

### 4.3 Atomkit 高对称路径自动生成

手动查找布里渊区高对称路径容易出错。Atomkit 可以从结构文件自动生成标准 KPT：

```bash
# 自动生成 FCC 结构的高对称路径 KPT
echo -e "3\n301\n3\n101 PRIMCELL.STRU\n0.06" | atomkit
# 生成 KLINES 文件，即为 NSCF 能带计算的 KPT
```

`0.06` 是 kspacing 参数（Å⁻¹），控制相邻 k 点的间距，值越小，能带越平滑，计算量越大。

### 4.4 NSCF 找不到电荷密度文件的处理

常见错误：NSCF 报"找不到 SPIN1_CHG.cube"。原因和解决方法：

1. **路径问题**：ABACUS 默认在 `OUT.suffix/` 目录下找电荷密度。确认 NSCF 的 `suffix` 与 SCF 完全一致
2. **读取不同目录的电荷密度**：在 INPUT 中设置 `read_file_dir /path/to/scf/OUT.suffix/`
3. **SCF 未收敛就运行 NSCF**：先检查 `running_scf.log` 是否有 `convergence is achieved`
