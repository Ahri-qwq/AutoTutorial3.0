## 六、结语

过渡金属氧化物、稀土化合物等强关联体系是第一性原理计算中的常见挑战。ABACUS 在 LCAO 基组下实现了 DFT+U 功能，参数设置直观，通过 `dft_plus_u`、`orbital_corr`、`hubbard_u` 三个参数即可快速上手。

本教程以反铁磁 NiO 为例，展示了从输入文件准备到结果分析的完整流程，以及 occupation matrix control 功能的使用方式。使用过程中如有问题，欢迎在 [ABACUS GitHub](https://github.com/deepmodeling/abacus-develop) 提交 Issue。

---

**参考文献**

[1] Qu X, Xu P, Jiang H, He L, Ren X. DFT+U within the framework of linear combination of numerical atomic orbitals. *The Journal of Chemical Physics*. 2022;156(23):234104. [doi:10.1063/5.0090122](https://doi.org/10.1063/5.0090122)
