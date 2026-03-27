## 二、案例背景：为什么 NiO 需要 DFT+U

### 2.1 NiO 的电子结构特点

NiO 晶体中，Ni 以 +2 价存在，电子构型为 [Ar] 3d⁸。8 个 d 电子填充在 5 条 d 轨道上，d 轨道未满，轨道间的局域库仑排斥（on-site Coulomb interaction）很强。这种强局域化的 d 电子是 LDA/GGA 难以准确描述的典型情形。

### 2.2 纯 GGA 的问题

在标准 GGA-PBE 计算中，NiO 的 d 轨道能带过度展宽（过离域），导致：

- 预测的能隙显著低于实验值（~4 eV），甚至给出金属态
- 磁矩被低估

这不是收敛问题，而是泛函本身的系统性误差。

### 2.3 DFT+U 如何修正

DFT+U 的能量泛函可以写成：

$$E_{\rm DFT+U} = E_{\rm DFA} + E_U - E_{\rm dc}$$

其中 $E_{\rm DFA}$ 是标准 LDA/GGA 的能量，$E_U$ 是 Hubbard 修正项，$E_{\rm dc}$ 是双计数项（避免对已在 $E_{\rm DFA}$ 中部分包含的相互作用重复计数）。

修正的核心是 occupation matrix $n_{I,mm'}^\sigma$，它描述原子 $I$ 上 $d$（或 $f$）轨道的占据情况：

$$n_{I, m m^{\prime}}^\sigma=\frac{1}{N_{\mathbf{k}}} \sum_{n \mathbf{k}} f_{n \mathbf{k}}^\sigma\left\langle\psi_{n \mathbf{k}}^\sigma\left|\hat{P}_{I, m m^{\prime}}^\sigma\right| \psi_{n \mathbf{k}}^\sigma\right\rangle$$

U 修正会对高占据和低占据的轨道产生不同方向的能量移动，将它们分开，从而打开能隙。
