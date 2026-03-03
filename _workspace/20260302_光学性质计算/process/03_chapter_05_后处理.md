# 五、后处理：从介电函数到光学性质

PYATB 输出了介电函数张量的全部 9 个分量。本章用 Python 读取数据，先可视化介电函数，再计算其余光学性质。

## 5.1 介电函数的读取与可视化

以下脚本读取虚部和实部数据，分别绘制各分量随能量的变化：

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

RealPartData = 'Out/Optical_Conductivity/dielectric_function_real_part.dat'
ImagPartData = 'Out/Optical_Conductivity/dielectric_function_imag_part.dat'

# 读取数据
with open(ImagPartData, 'r') as f:
    labels_imag = f.readline().strip().split()[1:]
data_imag = np.loadtxt(ImagPartData, skiprows=1)
omega = data_imag[:, 0]
imag_parts = data_imag[:, 1:]

with open(RealPartData, 'r') as f:
    labels_real = f.readline().strip().split()[1:]
data_real = np.loadtxt(RealPartData, skiprows=1)
real_parts = data_real[:, 1:]

# 绘图
markers = ['o', '*', '^', 'd', 's', 'p', '+', 'x', 'h']
fig = plt.figure(constrained_layout=True)
gs = GridSpec(1, 2, figure=fig)

ax_imag = fig.add_subplot(gs[0, 0])
for i, (imag_part, marker) in enumerate(zip(imag_parts.T, markers), start=1):
    ax_imag.plot(omega, imag_part, linestyle='-', marker=marker,
                 markersize=5, label=labels_imag[i])
ax_imag.legend()
ax_imag.set_title('Dielectric Function Imaginary Part')
ax_imag.set_xlabel('Energy (eV)')
ax_imag.set_ylabel(r'$\epsilon_2$')

ax_real = fig.add_subplot(gs[0, 1])
for i, (real_part, marker) in enumerate(zip(real_parts.T, markers), start=1):
    ax_real.plot(omega, real_part, linestyle='-', marker=marker,
                 markersize=5, fillstyle='none', mew=0.4, label=labels_real[i])
ax_real.legend()
ax_real.set_title('Dielectric Function Real Part')
ax_real.set_xlabel('Energy (eV)')
ax_real.set_ylabel(r'$\epsilon_1$')

plt.show()
```

由于 SiO₂ 是各向同性立方结构，三个对角分量（xx、yy、zz）完全重合，非对角分量为零。

实部 $\epsilon_1(0) \approx 1.90$，对应折射率 $n(0) \approx 1.38$。虚部 $\epsilon_2$ 在约 10 eV 以下（带隙以下）接近 0，之后随能量上升出现明显吸收峰。

## 5.2 各光学性质的计算

取各向同性情况下的平均值（(xx+yy+zz)/3）进行后续计算：

```python
def read_and_average(filename):
    """读取介电函数文件，返回能量轴和各向同性平均值"""
    data = np.loadtxt(filename, skiprows=1)
    energy = data[:, 0]
    # 取 xx(列1)、yy(列5)、zz(列9) 的平均值
    xx = data[:, 1]
    yy = data[:, 5]
    zz = data[:, 9]
    return energy, (xx + yy + zz) / 3

energy, eps2 = read_and_average(ImagPartData)    # 虚部
_,     eps1 = read_and_average(RealPartData)     # 实部

# ── 计算各光学量 ──────────────────────────────────────────
eps_abs = np.sqrt(eps1**2 + eps2**2)             # |ε|

# 折射率
n = np.sqrt((eps_abs + eps1) / 2)

# 消光系数
kappa = np.sqrt((eps_abs - eps1) / 2)

# 吸收系数（单位：cm⁻¹）
# ω = E/ħ，ħ ≈ 4.136e-15 eV·s，c = 3e10 cm/s
hbar = 4.135667696e-15   # eV·s
c    = 3.0e10            # cm/s（注意：此处需用 cm/s 而非 m/s，以得到 cm⁻¹ 单位）
omega_rad = energy / hbar                         # s⁻¹
alpha = np.sqrt(2) * omega_rad / c * np.sqrt(eps_abs - eps1)

# 能量损失函数
L = eps2 / (eps1**2 + eps2**2)

# ── 绘图 ──────────────────────────────────────────────────
fig, axs = plt.subplots(2, 2, figsize=(10, 7.5))

axs[0, 0].plot(energy, n, 'b-')
axs[0, 0].set_ylabel('Refraction Index $n$')
axs[0, 0].set_xlabel('Energy (eV)')
axs[0, 0].grid(True)

axs[0, 1].plot(energy, kappa, 'm-')
axs[0, 1].set_ylabel('Extinction Coefficient $\\kappa$')
axs[0, 1].set_xlabel('Energy (eV)')
axs[0, 1].grid(True)

axs[1, 0].plot(energy, alpha, 'y-')
axs[1, 0].set_ylabel(r'Absorption Coefficient (cm$^{-1}$)')
axs[1, 0].set_xlabel('Energy (eV)')
axs[1, 0].set_yscale('log')
axs[1, 0].grid(True)

axs[1, 1].plot(energy, L, 'r-')
axs[1, 1].set_ylabel('Energy Loss Function $L$')
axs[1, 1].set_xlabel('Energy (eV)')
axs[1, 1].grid(True)

plt.suptitle('Optical Properties of SiO₂ vs Energy')
plt.tight_layout()
plt.show()
```

> **注意**：代码中各向同性平均的列索引需与输出文件列顺序对应：
> 列 1 = xx，列 5 = yy，列 9 = zz（从 0 开始计数）。

## 5.3 结果物理解读

**介电函数**

- $\epsilon_1(0) \approx 1.90$，对应静态介电常数，与 SiO₂ 的实验值（~2.1）接近（PBE 低估带隙导致轻微偏差）
- $\epsilon_2$ 在约 9 eV 以下接近零，反映 SiO₂ 宽带隙绝缘体特征（实验带隙 ~8.9 eV，PBE 计算约低估 1~2 eV）
- $\epsilon_1$ 在 $\epsilon_2$ 峰值附近出现过零点，对应强吸收区

**折射率与消光系数**

- 低能区（可见光及近紫外）折射率约 1.38，接近实验值（~1.46，偏差来自 PBE 低估带隙）
- 高能区（>10 eV，深紫外）消光系数 $\kappa$ 明显上升，对应强烈光吸收

**吸收系数**

- 吸收系数以对数坐标显示，在带隙以下接近 0（通常 < 10² cm⁻¹）
- 超过光学带隙后急剧上升，可达 10⁶ cm⁻¹ 量级
- 吸收边（$\alpha$ 开始上升的能量）可用于估读 DFT 计算的光学带隙

**能量损失函数**

- $L(\omega)$ 的峰值对应体等离激元频率（plasmon resonance）
- SiO₂ 的等离激元峰位于 ~20 eV 左右，是其典型特征

> **提高精度的建议**：
> - 增大 PYATB 的 `grid`（如 40×40×40）可改善 k 积分收敛性
> - 减小 `eta`（如 0.05 eV）可获得更尖锐的谱线，但需配合更密的 `grid`
> - 若需准确带隙，改用 HSE06 泛函或 G₀W₀ 修正
