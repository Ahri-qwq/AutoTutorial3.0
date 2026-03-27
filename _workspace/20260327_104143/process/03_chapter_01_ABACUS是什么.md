# 第一章：ABACUS 是什么

ABACUS（Atomic-orbital Based Ab-initio Computation at UStc，中文名**原子算筹**）是国产开源第一性原理计算软件，基于密度泛函理论（DFT）实现。由中国科学技术大学（USTC）团队主导开发，现托管于 DeepModeling 社区。

## 功能定位

ABACUS 支持两种基组：

- **平面波（PW）**：精度高，参数简单，适合中小体系
- **数值原子轨道（LCAO）**：计算效率更高，适合大体系和需要输出 Hamiltonian 矩阵的场景

主要计算类型包括：电子自洽迭代（SCF）、结构优化（relax/cell-relax）、分子动力学（MD）、态密度、能带结构等。

采用模守恒赝势和周期性边界条件，可对晶格对称性、布里渊区 k 点对称性进行分析。

## 获取途径

| 资源 | 地址 |
|------|------|
| GitHub 仓库 | https://github.com/deepmodeling/abacus-develop |
| Gitee 镜像 | https://gitee.com/deepmodeling/abacus-develop |
| 官方文档 | https://abacus.deepmodeling.com/en/latest/ |
| 中文教程 | https://mcresearch.github.io/abacus-user-guide/ |

推荐通过 Docker 或 conda 安装，官方也提供了 Bohrium 云计算平台上的预配置镜像，无需本地安装即可运行。
