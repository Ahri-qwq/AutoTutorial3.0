# 一、计算原理

## 1.1 介电函数与线性光学性质

介电函数 $\epsilon(\omega)$ 是复数函数，描述材料对电磁场的响应：

$$
\epsilon(\omega) = \epsilon_1(\omega) + i\epsilon_2(\omega)
$$

**实部 $\epsilon_1$** 与材料的折射、色散特性相关；**虚部 $\epsilon_2$** 与光的吸收相关，其起始上升位置对应材料的光学带隙。

由介电函数可以导出一系列线性光学性质：

| 物理量 | 符号 | 计算公式 |
|--------|------|----------|
| 折射率 | $n$ | $\displaystyle n = \sqrt{\frac{\sqrt{\epsilon_1^2+\epsilon_2^2}+\epsilon_1}{2}}$ |
| 消光系数 | $\kappa$ | $\displaystyle \kappa = \sqrt{\frac{\sqrt{\epsilon_1^2+\epsilon_2^2}-\epsilon_1}{2}}$ |
| 吸收系数 | $\alpha$ | $\displaystyle \alpha = \frac{\sqrt{2}\,\omega}{c}\sqrt{\sqrt{\epsilon_1^2+\epsilon_2^2}-\epsilon_1}$ |
| 能量损失函数 | $L$ | $\displaystyle L = \operatorname{Im}\!\left(\frac{-1}{\epsilon}\right) = \frac{\epsilon_2}{\epsilon_1^2+\epsilon_2^2}$ |
| 反射率 | $R$ | $\displaystyle R = \frac{(n-1)^2+\kappa^2}{(n+1)^2+\kappa^2}$ |

> **单位说明**：吸收系数 $\alpha$ 的单位为 cm⁻¹。将能量 $E$（eV）换算为角频率时，
> $\omega = E/\hbar$，其中 $\hbar \approx 4.136\times10^{-15}\ \text{eV·s}$；
> 光速取 $c = 3\times10^{10}\ \text{cm/s}$。

## 1.2 ABACUS + PYATB 计算方法

### 紧束缚模型的构建

ABACUS 采用数值原子轨道（NAO）作为基函数展开 Kohn-Sham 波函数。对于给定 $\mathbf{k}$ 点，Kohn-Sham 方程在 NAO 基下变为广义本征值问题：

$$
H(\mathbf{k})\,C_n(\mathbf{k}) = E_{n\mathbf{k}}\,S(\mathbf{k})\,C_n(\mathbf{k})
$$

其中哈密顿矩阵 $H(\mathbf{k})$ 和重叠矩阵 $S(\mathbf{k})$ 通过傅里叶变换由实空间矩阵得到：

$$
H_{\nu\mu}(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}\,H_{\nu\mu}(\mathbf{R}), \quad
S_{\nu\mu}(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}\,S_{\nu\mu}(\mathbf{R})
$$

ABACUS 完成自洽计算后，通过设置 `out_mat_hs2 = 1` 和 `out_mat_r = 1`，可以将以下实空间矩阵以稀疏格式输出：
$H_{\nu\mu}(\mathbf{R})$（哈密顿量）、$S_{\nu\mu}(\mathbf{R})$（重叠矩阵）以及偶极矩阵 $r_{\nu\mu,a}(\mathbf{R})$（$a=x,y,z$）。
PYATB 读取这三组矩阵，在任意 $\mathbf{k}$ 点重建哈密顿量，用于后续物性计算。

### Kubo-Greenwood 公式与介电函数

光电导率张量由 Kubo-Greenwood 公式给出：

$$
\sigma_{\alpha\beta}(\hbar\omega) = -\frac{ie^2\hbar}{NV_{\text{cell}}} \sum_{\mathbf{k}}\sum_{n,m}
\frac{f_{n\mathbf{k}}-f_{m\mathbf{k}}}{E_{n\mathbf{k}}-E_{m\mathbf{k}}}
\frac{\langle n\mathbf{k}|v_\alpha|m\mathbf{k}\rangle\langle m\mathbf{k}|v_\beta|n\mathbf{k}\rangle}
{\hbar\omega + E_{n\mathbf{k}} - E_{m\mathbf{k}} + i\eta}
$$

其中 $f_{n\mathbf{k}}$ 为 Fermi-Dirac 占据函数，$\eta$ 为展宽参数，速度矩阵元 $\langle n\mathbf{k}|v_\alpha|m\mathbf{k}\rangle$ 由偶极矩阵 $r_{\nu\mu,a}(\mathbf{k})$ 和本征矢量计算得到。

在 PYATB 中，介电函数虚部由下式直接计算（避免 $\omega=0$ 处的奇点处理问题）：

$$
\epsilon_2^{\alpha\beta}(\omega) = \frac{e^2\pi}{\epsilon_0\hbar}
\int\frac{d\mathbf{k}}{(2\pi)^3}\sum_{n,m}f_{nm}\,r^\alpha_{nm}r^\beta_{mn}\,\delta(\omega_{nm}-\omega)
$$

介电函数实部 $\epsilon_1$ 则通过 Kramers-Kronig 变换由 $\epsilon_2$ 求得：

$$
\epsilon_1^{\alpha\beta}(\omega) = \delta_{\alpha\beta} + \frac{2}{\pi}\mathcal{P}\int_0^\infty
\frac{\omega'\,\epsilon_2^{\alpha\beta}(\omega')}{\omega'^2-\omega^2}\,d\omega'
$$

整体计算流程如下：

```
ABACUS SCF 自洽计算
  │  (basis_type = lcao, out_mat_hs2 = 1, out_mat_r = 1)
  ↓
输出三组矩阵文件
  │  data-HR-sparse_SPIN0.csr  (哈密顿量)
  │  data-SR-sparse_SPIN0.csr  (重叠矩阵)
  │  data-rR-sparse.csr        (偶极矩阵)
  ↓
PYATB 计算介电函数
  │  (Kubo-Greenwood + Kramers-Kronig)
  ↓
输出 dielectric_function_*.dat
  ↓
Python 后处理
  └→ 折射率、消光系数、吸收系数、能量损失函数
```
