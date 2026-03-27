# 教程资料来源说明

**教程：** 《用 ABACUS 计算能带与态密度：从硅的带隙到铝的费米面》
**生成日期：** 2026-03-25

---

## 关于 abacus_user_guide__abacus_dos.html.md 的使用

该文档（来源：https://mcresearch.github.io/abacus-user-guide/abacus-dos.html）在研究阶段被用于两处：

- **一、SCF→NSCF 核心流程**：`out_chg 1` 输出电荷密度、`init_chg file` 读入、`symmetry 0` 关闭对称性的框架说明，与该文档一并来自 `2024秋计算材料学`
- **三、Al DOS 案例**：Al FCC 4原子单胞 STRU、SCF/NSCF INPUT 参数（smearing、out_dos、dos_sigma）、KPT 设置，**主要来源**即为该文档

---

## 关于 2.5 节 Atomkit 后处理的说明

`process/01_research.md` 中标注来源为"多个教程"，具体为：

- **`快速开始 ABACUS｜自洽 能带 态密度 结构优化.md`**：提供了 `config.json` 基本结构（bandfile/efermi/energy_range/kptfile）和 `abacus-plot -b` 命令
https://www.bohrium.com/notebooks/4641406377
- **`2024秋计算材料学-上机练习：ABACUS能带和态密度计算.md`**：补充了 `bandfig` 和 `dpi` 字段，与教程最终模板一致
https://www.bohrium.com/notebooks/21913576824

两个文档均在知识库中。

---

## 实测数据来源（testCLAUDE 写回）

以下数值来自 2026-03-25 在 Bohrium 平台的实际计算结果，非 LLM 估算：

| 数值 | 来源任务 | Job ID |
|------|---------|--------|
| Si EFERMI = 6.852 eV | Si NSCF Band | 22303116 |
| Si 带隙 = 0.858 eV（Line 模式） | Si NSCF Band | 22303116 |
| Al EFERMI = 10.931 eV | Al NSCF DOS | 22303129 |
| Al DOS at Fermi ≈ 1.59 states/eV/cell | Al NSCF DOS | 22303129 |
