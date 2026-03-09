# 四、案例验证：HCO 分子（非周期体系）

本章用 HCO 分子（甲酰基，非周期体系）验证混合规范的正确性，并直观展示速度规范的系统性误差。

## 4.1 体系设置

HCO 分子放置于含真空层的超胞中，沿所有方向留出足够真空（≥10 Å），确保周期性边界条件不影响分子响应。

**INPUT（混合规范，HCO 分子）**

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix           HCO_hybrid
calculation      md
esolver_type     tddft
pseudo_dir       ./
orbital_dir      ./

#Parameters (2.Basis)
basis_type       lcao
gamma_only       0

#Parameters (3.MD / 时间步)
md_nstep         2000          # 40 fs 演化
md_dt            0.02
md_tfirst        0

#Parameters (4.外场设置)
td_vext          1
td_stype         2             # 混合规范
td_tstart        1
td_tend          2000

td_vext_dire     1             # x 方向
td_ttype         0
td_gauss_amp     0.01
td_gauss_t0      300
td_gauss_sigma   0.2
td_gauss_freq    3.66          # eV/hbar
td_gauss_phase   0.0

#Parameters (5.输出)
out_dipole       1             # 非周期体系：输出电偶极矩
out_efield       1
```

> **注意**：对非周期体系，应输出电偶极矩（`out_dipole 1`），而非电流。分子体系中也可使用长度规范（`td_stype 0`），此时 `out_dipole` 同样适用。

非周期体系的 KPT 文件只需 Gamma 点：
```
K_POINTS
0
Gamma
1 1 1 0 0 0
```

## 4.2 结果：三种规范的对比

施加高斯脉冲激发后，从输出的 `SPIN1_DIPOLE` 文件读取各时刻的电偶极矩 $D(t)$，配合 `efield_0.dat` 中的外加电场 $E(t)$，通过傅里叶变换得到吸收谱：

$$\alpha_i(\omega) = \frac{\int D_i(t) e^{i\omega t} dt}{\int E_i(t) e^{i\omega t} dt}$$

其中 $i$ 为激发方向，$\alpha_i(\omega)$ 的虚部对应吸收强度。

后处理使用 ABACUS 自带的 `tools/plot-tools` 脚本：

```bash
python tools/plot-tools/plot_absorption.py \
    --dipole SPIN1_DIPOLE \
    --efield efield_0.dat \
    --output absorption_HCO.dat
```

**三种规范的结果对比（HCO 分子）**

| 规范 | `td_stype` | 吸收谱精度 | 低频行为 |
|------|-----------|-----------|---------|
| 长度规范 | `0` | 参考基准 | 正常 |
| 混合规范 | `2` | 与长度规范完全一致 | 正常 |
| 速度规范 | `1` | 明显偏差 | 低频发散 |

赵昊天等的计算结果显示，对 HCO 分子，长度规范与混合规范的响应电流和吸收谱结果高度吻合，验证了混合规范在非周期体系中的正确性。速度规范的响应电流与前两者存在明显偏差，对应的吸收谱在低频区出现显著的虚假发散，这正是 NAO 基组下忽略轨道内部相位变化导致的系统性误差的体现。

这个对比实验说明：**混合规范可以完全替代长度规范用于非周期体系，同时为周期体系提供准确的结果**，而速度规范在 NAO 基组下无论对哪类体系都存在不可忽视的误差。

> 参考算例：https://gitee.com/mcresearch/abacus-user-guide/tree/master/examples/tddft
