# 第五章：Step 2 — 计算三阶力常数

进入 `3rd` 目录。本步骤使用 ABACUS + thirdorder 计算三阶力常数，最终得到 `FORCE_CONSTANTS_3RD`。

**背景说明**：thirdorder 目前只支持读取 VASP 和 QE 格式的输入/输出文件，
因此需要将 ABACUS 的结构文件转换为 POSCAR，将输出的受力封装为 vasprun.xml，再交给 thirdorder 处理。

## 5.1 产生微扰构型

首先，将 relax 后的 STRU 文件转换为 POSCAR 格式（`3rd/` 目录下已提供转换好的 POSCAR，
若需自行转换，可使用 ASE 完成）。

运行 thirdorder 产生微扰构型：

```bash
thirdorder_vasp.py sow 2 2 2 -2
```

参数说明：
- `sow`：播种模式，生成微扰超胞构型
- `2 2 2`：超胞扩展倍数（与二阶力常数计算保持一致）
- `-2`：近邻截断参数，负号表示取到第 2 近邻范围内的原子对

此命令为 Si 体系产生 **40 个** POSCAR 文件：`3RD.POSCAR.1`、`3RD.POSCAR.2`、……、`3RD.POSCAR.40`。

将这些 POSCAR 文件转换为 ABACUS 的 STRU 格式：

```bash
python pos2stru.py
```

> **注意**：转换时只能使用 `pos2stru.py`（内部调用 ASE），**不能使用 dpdata**。
> dpdata 在转换时会强制输出下三角晶格矩阵，等效于对晶格做了旋转，
> 导致 ABACUS 计算的受力方向同步旋转，三阶力常数因此完全错误。
> 这是一个容易忽视但影响致命的问题。

## 5.2 批量 SCF 计算

需要对 40 个微扰构型分别进行单点 SCF 计算，获取原子受力。
`run_stru.sh` 脚本会自动创建 `SCF-1/`、`SCF-2/`、……、`SCF-40/` 子目录并批量提交计算。

INPUT 中的关键参数：

```
calculation   scf
cal_force     1      # 必须输出受力
scf_thr       1e-8   # LCAO 基组：收敛阈值至少需要 1e-8
```

**`scf_thr` 对三阶力常数精度至关重要**：受力误差直接传递到力常数误差，进而影响声子寿命和热导率。
LCAO 基组至少需要 `1e-8`；若使用 PW 基组，需要至少 `1e-12`（原因详见第七章注意事项）。

40 个 SCF 任务彼此独立，建议并行提交到集群以节省时间。

## 5.3 封装受力为 vasprun.xml

40 个 SCF 计算完成后，执行：

```bash
python aba2vasp.py
```

该脚本读取每个 `SCF-*/` 中 ABACUS 输出的原子受力，将其封装为 thirdorder 能识别的 vasprun.xml 格式，放置在各 `SCF-*/` 目录中。

生成的 vasprun.xml 格式如下（以某一构型的受力为例）：

```xml
<modeling>
    <calculation>
        <varray name="forces">
            <v>1.865e-05 -0.04644196 -0.00153852</v>
            <v>-1.77e-05 -0.00037715 -0.00149635</v>
            <v>1.973e-05  0.002213   -0.00149461</v>
            <v>-1.976e-05 0.00065303 -0.0014804</v>
            ...（共 16 个原子的受力分量）
        </varray>
    </calculation>
</modeling>
```

每个 `<v>` 标签对应一个原子在 x、y、z 三个方向上的受力（单位：eV/Å）。

## 5.4 提取三阶力常数

执行以下命令，收集所有构型的受力并提取三阶力常数：

```bash
find SCF-* -name vasprun.xml | sort -n | thirdorder_vasp.py reap 2 2 2 -2
```

命令说明：
- `find SCF-* -name vasprun.xml | sort -n`：按编号顺序找到所有 vasprun.xml 文件
- `thirdorder_vasp.py reap`：收割模式，从受力数据中提取三阶力常数
- `2 2 2 -2`：与 `sow` 步骤使用相同的参数（超胞倍数和截断近邻数）

运行完成后生成 `FORCE_CONSTANTS_3RD`。
`shengbte/` 目录中已提供参考文件，可用于对照。

至此，两种力常数文件均已准备好，可以进入最后一步。
