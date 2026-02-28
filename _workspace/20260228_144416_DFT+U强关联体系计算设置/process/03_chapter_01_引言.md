# ABACUS 使用教程｜DFT+U 强关联体系计算

作者：AutoTutorial 3.0

---

## 1. 引言

过渡金属氧化物（如 NiO、MnO、FeO）是凝聚态物理和材料科学中的经典研究对象。用常规的 LDA 或 GGA 泛函计算这类体系时，往往得到错误的基态：NiO 本应是绝缘体，用 GGA 却常给出金属态，能隙严重低估，磁矩也不准确。

这个失败的根源在于 GGA 无法正确描述 Ni 的 d 轨道中强烈的电子关联效应。DFT+U 方法通过在 d/f 轨道上叠加 Hubbard 模型修正，以很低的计算代价改善这一问题。

本教程以反铁磁 NiO 为算例，展示 ABACUS LCAO 基组下 DFT+U 的完整计算流程：
- 设置反铁磁初始磁矩
- 配置 DFT+U 输入参数
- 提取并验证总能量、能隙、磁矩
- 使用 occupation matrix control (omc) 功能

> **注意**：ABACUS 的 DFT+U 功能**仅支持 LCAO 基组**（`basis_type lcao`），平面波基组下不可用。
