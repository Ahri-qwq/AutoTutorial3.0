# 第二章：算例获取与文件结构

## 2.1 算例下载

ABACUS 官方提供了隐式溶剂模型的完整算例（以 Pt 表面板层为例），可从以下地址获取：

**国内（Gitee）：**
```bash
git clone https://gitee.com/deepmodeling/abacus-develop.git
```
下载完成后进入：
```
abacus-develop/examples/implicit_solvation_model/Pt-slab/
```

**国际（GitHub）：**
```
https://github.com/deepmodeling/abacus-develop/tree/develop/examples/implicit_solvation_model/Pt-slab
```

本教程的理论演示案例使用的是更简单的 **H₂ 分子体系**（`suffix = H2`），该体系将一个氢分子置于超胞中，模拟其处于水溶液中的状态。这一体系计算快、结构简单，便于聚焦在溶剂化参数的理解上。

## 2.2 PW 基组计算的文件结构

使用平面波（PW）基组进行隐式溶剂 SCF 计算，需要以下文件：

```
计算目录/
├── INPUT          # 计算参数（本教程重点）
├── STRU           # 晶胞和原子坐标
├── KPT            # k 点设置
└── H_ONCV_PBE-1.0.upf   # H 原子赝势文件
```

> 注意：PW 基组无需轨道文件（`.orb`），这是与 LCAO 基组的重要区别。

## 2.3 体系说明：H₂ 分子在超胞中

本案例的物理图像：将一个孤立的 H₂ 分子放在一个足够大的超胞中，四周被隐式溶剂（水）包围。超胞需要足够大，以避免周期性镜像之间的相互作用。

INPUT 文件中的几个参数直接对应这一体系：

- `ntype = 1`：体系只有一种元素（H）
- `nbands = 2`：H₂ 有 2 个价电子，对应 2 条能带
- `symmetry = 0`：关闭晶体对称性（分子体系一般不用晶体对称性）
- `pseudo_dir = ./`：赝势文件放在当前目录下

STRU 文件需要定义超胞大小和 H₂ 的原子坐标。对于孤立分子，超胞通常取立方盒子，边长约 10-15 Å。KPT 文件只需取 Γ 点（`1 1 1 0 0 0`），分子体系不需要 k 点采样。具体文件格式可参考官方 Pt-slab 算例，按分子体系调整。
