## 一、安装与准备

### 安装 abacustest

```bash
# 通过 pip 安装
pip install abacustest

# 或从源码安装
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

安装完成后，执行 `abacustest model inputs -h` 可查看所有可用选项。

### 下载 APNS-pp-orb-v1 赝势轨道库

```bash
abacustest model inputs --download-pporb
```

命令执行完成后，当前目录下会出现两个文件夹：

```bash
$ ls
apns-orbitals-efficiency-v1  apns-pseudopotentials-v1
```

`apns-pseudopotentials-v1` 和 `apns-orbitals-efficiency-v1` 包含从 H 到 Bi 共 83 种元素推荐使用的赝势和轨道文件，均基于 APNS-pp-orb-v1。所有轨道为 DZP 水平；镧系元素的 4f 电子被视为核电子，选用截断半径 8 au 的轨道。

### 设置环境变量

将赝势和轨道路径写入环境变量，后续命令会自动读取：

```bash
export ABACUS_PP_PATH=/your/path/to/apns-pseudopotentials-v1
export ABACUS_ORB_PATH=/your/path/to/apns-orbitals-efficiency-v1
```

如果不设置环境变量，也可以在每次命令中通过 `--pp` 和 `--orb` 选项显式指定路径。
