# 三、混合规范：一个参数解决问题

## 3.1 物理思想

混合规范的核心想法直接：既然速度规范在 NAO 展开中遗漏了矢势诱导的相位，那就显式地把它补回来。

具体做法是在每个原子轨道前引入一个依赖于矢量势的相位因子：

$$\tilde{\phi}_\mu(\mathbf{r}, t) = e^{i\mathbf{A}(t)\cdot\mathbf{R}_\mu} \phi_\mu(\mathbf{r})$$

其中 $\mathbf{R}_\mu$ 是第 $\mu$ 个原子轨道的中心位置，$\mathbf{A}(t)$ 是时变矢势。波函数展开变为：

$$\psi(\mathbf{r}, t) = \sum_{\mu} c_\mu(t) \, e^{i\mathbf{A}(t)\cdot\mathbf{R}_\mu} \phi_\mu(\mathbf{r})$$

在这个修正下，哈密顿量矩阵元和重叠矩阵均需相应修改。由于混合规范同时包含矢量势 $\mathbf{A}$ 和标量场，它本质上是一种新的规范选择——相比长度规范，混合规范下的外场项保持周期性，不破坏布洛赫定理；相比速度规范，它补偿了 NAO 基组不完备带来的相位误差。

混合规范和长度规范之间满足严格的规范变换关系，理论上结果完全等价。在 ABACUS 已有测试中，混合规范与长度规范（非周期体系）的结果高度吻合，而速度规范则存在明显偏差。

> 深入的理论推导参见原始论文：赵昊天, 何力新, *J. Chem. Theory Comput.* **2025**, 21, 3335–3341 (DOI: 10.1021/acs.jctc.5c00111)

## 3.2 在 ABACUS 中启用混合规范

启用混合规范只需设置一个参数：

```
td_stype   2     # 0=长度规范  1=速度规范  2=混合规范
```

以下是一个完整的 TDDFT INPUT 示例（Si 周期体系，弱场线性响应）：

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix           Si_tddft
calculation      md
esolver_type     tddft
pseudo_dir       ./
orbital_dir      ./

#Parameters (2.Basis)
basis_type       lcao
gamma_only       0

#Parameters (3.MD / 时间步)
md_nstep         4000          # 总时间步数 = 4000 × 0.02 fs = 80 fs
md_dt            0.02          # 时间步长 (fs)，建议 0.02–0.05
md_tfirst        0             # 离子初始温度，0 表示固定离子

#Parameters (4.外场设置)
td_vext          1             # 开启外加电场
td_stype         2             # 混合规范（重点参数）

td_tstart        1             # 第 1 步开始施加电场
td_tend          4000          # 与 md_nstep 相同（全程施加）

td_vext_dire     1             # 沿 x 方向施加电场

td_ttype         0             # 高斯型脉冲
td_gauss_amp     0.001         # 振幅 (V/Å)，弱场取 ~0.001
td_gauss_t0      50            # 脉冲中心位于第 50 步
td_gauss_sigma   0.5           # 高斯宽度 (fs)，宽带覆盖
td_gauss_freq    0.0           # 中心频率 0 (宽带激发)
td_gauss_phase   0.0           # 初相位

#Parameters (5.输出)
out_current      1             # 输出电流（周期体系/混合规范用此项）
out_efield       1             # 输出外加电场时程
```

> **与速度规范的区别：** 将 `td_stype 1` 改为 `td_stype 2`，其余所有参数保持不变。

## 3.3 外场参数速查

| 参数 | 说明 | 典型值 |
|------|------|--------|
| `td_vext` | 外场开关 | `1`（开） |
| `td_stype` | 规范选择 | `0`/`1`/`2`（长度/速度/混合） |
| `td_tstart` | 开始步数 | `1` |
| `td_tend` | 结束步数 | 通常等于 `md_nstep` |
| `td_vext_dire` | 方向 | `1`=x, `2`=y, `3`=z |
| `td_ttype` | 电场形状 | `0`=高斯, `1`=梯形, `2`=三角, `3`=阶跃 |
| `td_gauss_amp` | 振幅 (V/Å) | 线性响应 ~0.001；强场 0.01–0.1 |
| `td_gauss_t0` | 中心步数 | 脉冲中心位置 |
| `td_gauss_sigma` | 宽度 (fs) | 0.2–1.0；越小带宽越大 |
| `td_gauss_freq` | 中心频率（单位请参考 ABACUS 在线文档） | `0` 为宽带；[需查阅文档确认] |
| `out_current` | 电流输出 | `1`（周期体系用） |
| `out_dipole` | 偶极矩输出 | `1`（非周期体系用） |
| `out_efield` | 外场输出 | `1` |

**规范选择建议**

| 体系类型 | 推荐规范 | `td_stype` |
|---------|---------|------------|
| 非周期（分子、团簇） | 长度规范或混合规范 | `0` 或 `2` |
| 周期（晶体） | **混合规范**（推荐） | `2` |
| 周期（旧方法兼容） | 速度规范 | `1` |

在 ABACUS LCAO 框架下，**混合规范适用于所有体系**，精度优于速度规范，与长度规范等价（非周期）。对周期体系，混合规范是目前最可靠的选择。
