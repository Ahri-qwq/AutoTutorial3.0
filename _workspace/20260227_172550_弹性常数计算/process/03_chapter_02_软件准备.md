# 二、软件准备

## 2.1 安装 abacustest

abacustest 是 ABACUS 的辅助工具，提供输入文件准备、高通量任务管理、结果后处理等功能。通过 pip 安装：

```bash
pip install abacustest
```

也可从 GitHub 获取最新版本：

```bash
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

安装完成后，执行 `abacustest -h` 验证安装是否成功。

## 2.2 准备赝势和轨道文件

abacustest 在准备 LCAO 计算的输入文件时，会自动从环境变量 `ABACUS_PP_PATH` 和 `ABACUS_ORB_PATH` 指定的路径中查找赝势和数值原子轨道文件。在运行前需设置这两个环境变量：

```bash
export ABACUS_PP_PATH=/path/to/pseudopotentials
export ABACUS_ORB_PATH=/path/to/orbitals
```

本教程的案例使用 APNS-v1 赝势库与 efficiency 数值原子轨道基组。赝势和轨道文件可从 ABACUS 官方渠道获取。

> **说明：** 若计划在 Bohrium 云平台上运行，可使用平台提供的 abacustest App，无需手动配置环境变量，赝势和轨道文件由平台统一管理。
