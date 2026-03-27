## 一、介绍

NiO 是过渡金属氧化物中最经典的强关联体系之一。实验上它是一个绝缘体，能隙约 4 eV，并具有 II 型反铁磁序。但若用标准的 LDA 或 GGA 泛函进行计算，往往得到的能隙远小于实验值，甚至错误地预测为金属态。

问题的根源在于 d 轨道电子的局域化。LDA/GGA 对电子间相互作用的描述是均匀化的，无法准确反映局域 d 电子之间强烈的库仑排斥。DFT+U 方法通过引入 Hubbard U 参数，在 d（或 f）轨道上附加一项平均场型的 Hartree-Fock 修正，抑制自相互作用带来的非物理离域化，从而更准确地描述强关联体系。

ABACUS 在 LCAO 基组（`basis_type lcao`）下实现了 DFT+U 功能，理论细节见参考文献 [1]。本教程以反铁磁 NiO 的 SCF 计算为主线，完整演示：

- 如何在 STRU 文件中设置反铁磁初始磁矩
- 如何在 INPUT 中开启 DFT+U 并配置关键参数
- 如何读取和理解 occupation matrix 输出
- 如何使用 occupation matrix control（omc）功能固定 occupation matrix

**适用读者：** 熟悉 ABACUS 基本 SCF 计算流程，希望了解如何对过渡金属体系启用 DFT+U 的用户。
