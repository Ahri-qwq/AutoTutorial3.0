# 五、核心案例：Si 单胞（周期体系）

Si 是验证周期体系 RT-TDDFT 方法的标准体系。本章包含两个计算任务：
1. 弱场激发 → 线性响应 → 介电函数
2. 强场激发 → 非线性响应 → 高次谐波谱（HHG）

## 5.1 弱场：介电函数（线性响应）

### 输入文件准备

使用第二章中的 STRU 和 KPT，仅修改 INPUT 如下：

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix           Si_linear
calculation      md
esolver_type     tddft
pseudo_dir       ./
orbital_dir      ./

#Parameters (2.Basis)
basis_type       lcao
gamma_only       0

#Parameters (3.MD / 时间步)
md_nstep         4000          # 80 fs，足够分辨率
md_dt            0.02
md_tfirst        0

#Parameters (4.外场：宽带弱场激发)
td_vext          1
td_stype         2             # 混合规范

td_tstart        1
td_tend          4000

td_vext_dire     1             # x 方向
td_ttype         0
td_gauss_amp     0.001         # 弱场振幅 (V/Å)，线性响应区间
td_gauss_t0      50
td_gauss_sigma   0.5           # 较宽的高斯，覆盖宽频域
td_gauss_freq    0.0
td_gauss_phase   0.0

#Parameters (5.输出)
out_current      1             # 周期体系：输出电流
out_efield       1
```

确保场强足够弱（`td_gauss_amp ≤ 0.001 V/Å`），使体系响应处于线性区间。

### 运行

先跑 SCF，再跑 TDDFT：
```bash
# 第一步：SCF
cp -r scf_input/ tddft_linear/
cd tddft_linear/
# 替换 INPUT 为 TDDFT 版本
mpirun -n 8 abacus
```

### 后处理：电流 → 介电函数

TDDFT 演化结束后，输出的 `SPIN1_CURRENT` 文件包含各时刻的电流密度 $J(t)$。

电流 $J(t)$ 与电偶极矩的时间导数成正比：$J(t) \propto \dot{D}(t)$，因此傅里叶变换时：

$$D(\omega) \propto \frac{J(\omega)}{i\omega}$$

介电函数的虚部（对应光吸收）由此得到：

$$\varepsilon_2(\omega) \propto \frac{\text{Im}[J(\omega)]}{\omega \cdot E(\omega)}$$

使用 ABACUS 的后处理脚本：

```bash
python tools/plot-tools/plot_current.py \
    --current SPIN1_CURRENT \
    --efield  efield_0.dat \
    --output  dielectric_Si.dat
```

> **说明**：`out_current 1` 开启后，ABACUS 会在 `OUT.*` 目录中输出电流文件（请参考实际运行后的文件名）。后处理脚本的具体用法请参考 `tools/plot-tools` 目录下的说明文件，或 ABACUS 官方算例仓库中的示例。

### 结果与分析

Si 的介电函数有两个主要特征峰，分别位于约 3.4 eV 和 4.3 eV，对应 E₁ 和 E₂ 带间跃迁。

| 规范 | 介电函数精度 | 说明 |
|------|------------|------|
| 混合规范（`td_stype = 2`） | 准确 | 峰位、峰高与参考值吻合 |
| 速度规范（`td_stype = 1`） | 偏低且展宽 | 小基组下系统性低估 |

赵昊天等的计算表明，速度规范在 Si 原胞计算中存在明显的系统性误差，通过经验修正方案可部分缓解，但混合规范无需任何修正即可稳定地给出准确且物理一致的结果。

---

## 5.2 强场：高次谐波谱（HHG，非线性响应）

当场强增大到非线性区间（0.001–0.01 V/Å 量级），体系响应不再线性，电流信号中出现基频的高次谐波成分。

### 强场 INPUT

将线性响应 INPUT 中的场强参数修改如下（其余不变）：

```
#Parameters (4.外场：强场激发)
td_vext          1
td_stype         2             # 混合规范

td_tstart        1
td_tend          4000

td_vext_dire     1
td_ttype         0
td_gauss_amp     0.01          # 强场振幅，比线性响应大 10 倍
td_gauss_t0      50
td_gauss_sigma   0.5
td_gauss_freq    0.372         # 约 1.0 eV，基频激光频率
td_gauss_phase   0.0
```

> **参数说明**：`td_gauss_freq` 设置激光的中心频率（单位 eV/ℏ，具体单位请参考 ABACUS 在线文档）。高次谐波谱通过对 Si 单胞施加一系列不同场强（0.001、0.003、0.005、0.01 V/Å）的高斯脉冲来测试。

### 后处理：电流 → 高次谐波谱

HHG 谱通过对电流信号 $J(t)$ 作傅里叶变换得到：

$$I_{\text{HHG}}(\omega) \propto |\omega^2 \cdot J(\omega)|^2$$

对 $J(\omega)$ 取绝对值平方（强度谱），在频率轴上以基频 $\omega_0$ 的倍数标注，即可得到 HHG 谱。Python 处理示例：

```python
import numpy as np

# 读取电流文件（具体格式参考 ABACUS OUT.* 目录中的输出说明）
t, J = np.loadtxt('SPIN1_CURRENT', unpack=True, usecols=(0, 1))
dt = t[1] - t[0]

# FFT
N = len(J)
freq = np.fft.rfftfreq(N, d=dt)     # 频率轴
Jw   = np.fft.rfft(J) * dt

# HHG 强度谱
I_HHG = (freq**2 * np.abs(Jw))**2

# 以基频 omega0 为单位绘图
omega0 = 0.372   # 激光基频（eV，请根据实际计算参数设置）
harmonic_order = freq / omega0
```

### 5.3 场强从弱到强：谐波阶次的演化

赵昊天等对 Si 单胞进行了系统的强场测试，将高斯脉冲场强从弱到强递增：

| 场强 (V/Å) | 响应类型 | 出现的谐波阶次 |
|-----------|---------|--------------|
| ~0.001 | 线性响应 | 仅 1 次（基频） |
| ~0.003 | 弱非线性 | 1、3 次 |
| ~0.005 | 非线性 | 1、3、5 次 |
| ~0.010 | 强非线性 | 1、3、5、7 次 |

实验规律清晰：随场强增大，HHG 谱中依次出现 1、3、5、7 次谐波信号（奇次谐波，源于 Si 的反演对称性）。这一结果与 HHG 的理论预期完全吻合，也证明混合规范在强场和弱场条件下均能给出可靠结果。

> **混合规范的计算效率优势**：相比速度规范，混合规范在处理非局域势时更高效，在相同基组下精度更高。这意味着可以使用较小的基组（如 2s2p1d）完成高精度的 HHG 计算，而速度规范为达到同等精度则需要更大的基组（计算成本显著增加）。
