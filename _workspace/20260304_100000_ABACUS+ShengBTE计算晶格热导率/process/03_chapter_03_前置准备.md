# 第三章：前置准备

## 3.1 案例文件

案例文件可从 Gitee 下载：

```bash
git clone https://gitee.com/mcresearch/abacus-user-guide.git
```

进入 `examples/interface_ShengBTE/LCAO/` 目录，结构如下：

```
LCAO/
├── 2nd/              # 二阶力常数计算
│   ├── STRU          # relax 后的结构文件
│   ├── INPUT         # ABACUS 输入文件
│   ├── KPT           # k 点设置
│   ├── setting.conf  # Phonopy 超胞设置
│   ├── band.conf     # Phonopy 声子谱和力常数设置
│   └── au2si.py      # 单位转换脚本
├── 3rd/              # 三阶力常数计算
│   ├── POSCAR        # 由 STRU 转换而来（已提供）
│   ├── pos2stru.py   # POSCAR→STRU 转换脚本（调用 ASE）
│   ├── run_stru.sh   # 批量提交 SCF 的脚本
│   └── aba2vasp.py   # 将 ABACUS 受力封装为 vasprun.xml 的脚本
└── shengbte/         # ShengBTE 计算
    ├── CONTROL       # ShengBTE 参数文件
    ├── FORCE_CONSTANTS_2ND   # 参考文件（已提供）
    ├── FORCE_CONSTANTS_3RD   # 参考文件（已提供）
    └── Ref/          # 参考计算结果
```

## 3.2 软件依赖

| 软件 | 说明 |
|------|------|
| ABACUS ≥3.2.0 | 第一性原理计算，提供结构优化和 SCF 受力 |
| Phonopy ≥2.19.1 | 二阶力常数计算，已支持 ABACUS 接口 |
| ASE | 原子结构格式转换（pos2stru.py 内部调用） |
| thirdorder | 三阶力常数，依赖 VASP/QE 格式的结构和受力输入 |
| ShengBTE | BTE 求解，输出 κ(T) |

相关文档：
- Phonopy 接口：https://abacus.deepmodeling.com/en/latest/advanced/interface/phonopy.html
- ASE 接口：https://abacus.deepmodeling.com/en/latest/advanced/interface/ase.html
- ShengBTE / thirdorder：https://bitbucket.org/sousaw/shengbte/src/master/
