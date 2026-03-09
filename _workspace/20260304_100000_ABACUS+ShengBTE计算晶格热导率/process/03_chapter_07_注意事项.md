# 第七章：关键注意事项

整个计算流程中有几处容易出错，汇总如下。

## 7.1 必须设置 FULL_FORCE_CONSTANTS = .TRUE.

`band.conf` 中必须加入这一行：

```
FULL_FORCE_CONSTANTS = .TRUE.
```

不设置时，Phonopy 默认只输出精简形式的力常数文件（仅含不等价部分），ShengBTE 读取时会因格式不匹配而直接报错。这是一个不设置时不会有任何警告、但结果完全无法运行的问题。

## 7.2 二阶力常数的单位转换不能跳过

ABACUS 结合 Phonopy 输出的 `FORCE_CONSTANTS` 单位为 **eV/(Å·au)**，
而 ShengBTE 要求 `FORCE_CONSTANTS_2ND` 的单位为 **eV/Å²**（1 au = 0.52918 Å）。
如果跳过 `au2si.py` 的单位转换步骤，直接将 `FORCE_CONSTANTS` 重命名为 `FORCE_CONSTANTS_2ND`，
ShengBTE 的计算不会报错，但热导率结果的数量级完全错误。

## 7.3 结构格式转换不能用 dpdata

将 `3RD.POSCAR.*` 转换为 STRU 时，只能使用 `pos2stru.py`（调用 ASE），**不能用 dpdata**。

dpdata 在转换时会将晶格矩阵强制变换为下三角形式，这等效于对整个晶格做了一次旋转。旋转后，晶体结构的物理是等价的，但 ABACUS 计算出的原子受力方向也会对应旋转——而 thirdorder 在 `reap` 阶段并不知道这个旋转，最终提取出的三阶力常数方向完全错误。

## 7.4 SCF 收敛阈值须足够严格

计算三阶力常数时，受力精度直接影响力常数精度，进而影响声子寿命和热导率。

| 基组 | `scf_thr` 最低要求 |
|------|------------------|
| LCAO | 1e-8 |
| PW   | 1e-12 |

PW 基组对收敛阈值要求更严格（需要 1e-12），是因为平面波基组下哈密顿矩阵元的数值精度对 scf_thr 更敏感。如果不满足此要求，受力误差会使三阶力常数出现较大噪声，导致热导率结果不收敛。

## 7.5 CONTROL 中的 scell 须与 thirdorder 保持一致

CONTROL 文件中的 `scell(:)=2 2 2` 必须与 thirdorder 中使用的超胞参数（`thirdorder_vasp.py sow 2 2 2 -2` 中的 `2 2 2`）完全一致。
如果不一致，ShengBTE 在读取 `FORCE_CONSTANTS_3RD` 时会出现维度不匹配的错误，或者无声地给出错误结果。
