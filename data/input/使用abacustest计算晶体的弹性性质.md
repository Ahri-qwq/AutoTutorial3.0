# 弹性张量的概念及其在晶体中的计算

当材料中存在内部应力时，材料的形状会相应发生变化，称为材料的应变。材料中的内部应力

$$
\sigma_{ij}
$$

与应变

$$
\varepsilon_{ij}
$$

均为二阶对称张量，均可以用3×3的对称矩阵表示。材料内部应力与应变的关系是材料的一种本构关系。在材料处于平衡位置附近，应力与应变之间的关系可以用广义胡克定律表示：

$$
\sigma_{ij} = C_{ijkl} \varepsilon_{kl}
$$

其中

$$
C_{ijkl}
$$

是材料的弹性张量，是一个四阶张量，有81个分量。但由于\$\$\\sigma\_{ij}\$\$和\$\$\\varepsilon\_{kl}\$\$都是二阶对称张量，\$\$C\_{ijkl}\$\$也具有一定的对称性：

$$
C_{ijkl} = C_{jikl}, C_{ijkl} = C_{ijlk}, C_{ijkl} = C_{klij}
$$

经过上面的对称性约化后，弹性张量

$$
C_{ijkl}
$$

中将剩余21个独立分量。利用Voigt记号可以简化应力、应变和弹性张量的表示。Voigt记号使用如下的对应：

$$
xx \rightarrow 1, yy \rightarrow 2, zz \rightarrow 3, yz \rightarrow 4, xz \rightarrow 5, xy \rightarrow 6
$$

则应力

$$
\sigma
$$

、应变\$\$\\varepsilon\$\$可用6维向量表示，弹性张量可用6维对称矩阵表示，应力与应变间的本构关系表示为

$$
\begin{bmatrix} \sigma_1 \\ \sigma_2 \\ \sigma_3 \\ \sigma_4 \\ \sigma_5 \\ \sigma_6 \\ \end{bmatrix} = \begin{bmatrix} C_{11} & C_{12} & C_{13} & C_{14} & C_{15} & C_{16}   \\ C_{21} & C_{22} & C_{23} & C_{24} & C_{25} & C_{26}   \\  C_{31} & C_{32} & C_{33} & C_{34} & C_{35} & C_{36}   \\  C_{41} & C_{42} & C_{43} & C_{44} & C_{45} & C_{46}   \\  C_{51} & C_{52} & C_{53} & C_{54} & C_{55} & C_{56}   \\  C_{61} & C_{62} & C_{63} & C_{64} & C_{65} & C_{66}   \\ \end{bmatrix} \begin{bmatrix} \varepsilon_1 \\ \varepsilon_2 \\ \varepsilon_3 \\ \varepsilon_4 \\ \varepsilon_5 \\ \varepsilon_6 \\ \end{bmatrix}
$$

目前在计算材料弹性张量时被广泛使用的方法是应力-应变法。这种方法的理论基础为上面的广义胡克定律。通过对材料施加一系列应变，并通过第一性原理计算获取材料的内部应力，即可拟合得到弹性张量中的部分矩阵元，遍历所有正应变与剪切应变的种类，即可拟合得到弹性张量中的所有矩阵元。

由于晶体的应力通常具有各向异性，弹性张量的具体形式将与晶体的取向有关，因此比较弹性张量的计算结果时需要确定晶体的取向。在Material Project数据库中记录的弹性张量数据通常会将晶体旋转到IEEE汇报介电常数张量时使用的取向[1]。此外，使用的原胞或惯用胞对弹性张量的形式也有影响。对Si、Cu等常见的立方晶系晶体，计算弹性模量时通常使用正方体形式的惯用胞，并且使用坐标轴与惯用胞晶格矢量重合的坐标系。

晶体的对称性会进一步对弹性张量的形式施加限制。例如，通过对称性分析，Si、Cu等立方晶系的晶体中只存在3个独立弹性张量元：

$$
C_{11}(=C_{22}=C_{33})
$$

，\$\$C\_{12}(=C\_{13}=C\_{23})\$\$和\$\$C\_{44}(=C\_{55}=C\_{66})\$\$，其它弹性张量元均为0。不同晶系中的独立弹性张量元可参看文献[2]。基于上面的介绍，为了计算一个晶体的弹性张量，需要的步骤为：

1. 优化晶体的晶胞和原子位置，获取无应力的晶体结构
2. （可选）将晶体旋转至IEEE 176/1987标准规定的标准取向
3. 对所有可能正应变和剪切应变，施加一批不同大小的应变后，使用第一性原理计算，获取晶体的应力
4. 拟合应力和应变，得到完整的弹性张量矩阵

# 使用abacustest计算材料的弹性性质

abacustest是一个ABACUS的辅助软件，已经包含一系列工作流，可使用ABACUS作为第一性原理计算引擎，完成输入文件准备、提交计算、结果后处理等功能。目前abacustest中已经包含了计算晶体弹性性质的工作流，以下我们以计算Si和金红石型TiO2的弹性模量为例，介绍其使用方法。

## 计算Si的弹性模量

### 生成Si的晶体结构，并使用ABACUS优化

首先需要生成一个Si晶体的结构。下面的脚本可使用ase生成一个Si晶体惯用胞的CIF文件：

```Python
from ase.build import bulk

si_conv = bulk('Si', cubic=True)
si_conv.write("Si_conv.cif")
```

接下来需要优化Si的晶体结构。abacustest中包含了使用结构文件准备ABACUS输入文件的功能，会自动使用环境变量ABACUS\_PP\_PATH和ABACUS\_ORB\_PATH中包含的赝势和数值原子轨道，根据命令生成指定任务的输入文件。本教程中采用APNS-v1中的赝势和efficiency基组。使用下面的命令生成优化Si晶体结构所需的输入文件：

```Python
abacustest model inputs -f Si_conv.cif --ftype cif --jtype cell-relax --lcao --folder-syntax Si
```

各个命令行参数的含义如下：

* `-f`: 结构文件名称。此处可使用多个结构文件，会对每个结构生成ABACUS输入文件夹
* `--ftype`: 结构文件格式
* `--jtype`: ABACUS任务类型。这里需要优化Si的晶胞结构和原子位置，因此使用cell-relax计算类型
* `--lcao`：选择是否使用更高效的LCAO基组
* `--folder-syntax`: 指定生成的ABACUS输入文件夹名称

执行完毕后，在当前目录生成了Si/文件夹，包含用于优化Si晶体结构的输入文件，其中INPUT文件中的参数适配ABACUS LTSv3.10，在很多情况下是合理的。部分参数，例如nspin、初始磁矩、DFT+U参数可以在准备输入文件夹时设置，也可以在准备好输入文件夹后再修改。

确认INPUT参数合理后，可以提交优化计算。

优化计算完成后，使用优化后的结构，重新准备一个ABACUS输入文件夹：

```Python
mkdir Si-elastic
cp Si/INPUT Si/Si* Si-elastic
cp Si/INPUT/OUT.ABACUS/STRU_ION_D Si-elastic/STRU
```

然后修改新文件夹中的INPUT文件为做SCF计算。

由于在准备Si的晶体结构文件时，已经将晶格矢量沿着坐标轴放置，因此这里不需要对晶体结构做旋转。

### 使用abacustest计算Si的弹性模量

接下来使用下面的命令准备计算弹性模量所需的输入文件：

```Python
abacustest model elastic prepare -j Si-elastic
```

该命令将在Si-elastic目录下生成一系列ABACUS输入文件夹，用于进行弹性性质计算。执行完毕后，会在Si-elastic目录下生成一系列以`deform-`开头的文件夹和`org`文件夹，分别包含对3种正应变、3种剪切应变施加4种不同大小的应变的结构的输入文件夹和计算原始结构的输入文件夹。使用`--norm`和`--shear`参数可分别指定最大正应变和最大剪切应变，默认值0.01的含义是将为每种应变生成±0.5%和±1%的变形结构。使用`--norelax`参数可选择不对变形结构做固定晶胞的优化，直接计算应力。注意重复执行该命令会直接删除已有的准备结果，不要再准备好后重复计算！

接下来可以提交计算。等待计算完成后，使用下面的命令对计算结果做后处理，得到计算结果：

```Python
abacustest model elastic post -j Si-elastic
```

在屏幕上会输出弹性张量，以及体模量、剪切模量、杨氏模量、泊松比等弹性性质的计算结果：

```Python
Model: elastic
Postprocessing elastic calculation for job: Si-elastic/
             bulk_modulus  shear_modulus  young_modulus  poisson_ratio
Si-elastic/     90.705919      65.134208     157.664088       0.210302

Si-elastic/     elastic_tensor:
            0             1             2          3          4          5
0  155.464972  5.832639e+01  5.832639e+01   0.000000   0.000000   0.000000
1   58.326393  1.554650e+02  5.832639e+01   0.000000   0.000000   0.000000
2   58.326393  5.832639e+01  1.554650e+02   0.000000   0.000000   0.000000
3    0.000000  0.000000e+00  0.000000e+00  76.177486   0.000000   0.000000
4    0.000000 -2.000000e-10  0.000000e+00   0.000000  76.177486   0.000000
5    0.000000  0.000000e+00 -2.000000e-10   0.000000   0.000000  76.177486

The postprocess is done. The metrics are saved in 'metrics.json', and the elastic results are saved in 'metrics_elastic.json'.
```

计算得到的Si的弹性张量中

$$
C_{11}
$$

=155.5 GPa，\$\$C\_{12}\$\$=58.2 GPa，\$\$C\_{44}\$\$=76.2 GPa，与Materials Project上的结果（\$\$C\_{11}\$\$=153 GPa，\$\$C\_{12}\$\$=57 GPa，\$\$C\_{44}\$\$=74 GPa，）接近，这表明成功使用ABACUS计算了Si的弹性张量。

## 计算金红石型TiO2的弹性模量

金红石型TiO2在Materials Project上的编号为mp-2657。从Materials Project上下载得到金红石型TiO2的晶体结构后，使用abacustest准备做结构优化的ABACUS输入文件：

```Python
abacustest model inputs -f TiO2.cif --ftype cif --lcao --jtype cell-relax --folder-syntax TiO2-rutile
```

提交ABACUS计算，并等待计算完成后，使用优化后的结构，重新准备一个ABACUS输入文件夹：

```Python
mkdir TiO2-rutile-elastic
cp TiO2-rutile/INPUT TiO2-rutile/Ti* TiO2-rutile/O* TiO2-rutile-elastic
cp TiO2-rutile/INPUT/OUT.ABACUS/STRU_ION_D TiO2-rutile-elastic/
```

按照与计算Si的弹性模量相同的流程，准备计算所需的输入文件，并提交计算，等待计算完成。

```Python
abacustest model elastic prepare -j TiO2-rutile-elastic
```

计算完成后，使用相同的方法进行结果处理：

```Python
abacustest model elastic post -j TiO2-rutile-elastic
```

屏幕输出结果为：

```Python
Model: elastic
Postprocessing elastic calculation for job: TiO2-rutile-elastic/
                      bulk_modulus  shear_modulus  young_modulus  poisson_ratio
TiO2-rutile-elastic/    220.705174     128.683753     323.230608       0.255911

TiO2-rutile-elastic/    elastic_tensor:
              0             1             2           3           4           5
0  2.812277e+02  1.581445e+02  1.559258e+02    0.000000    0.000000    0.000000
1  1.581445e+02  2.812277e+02  1.559258e+02    0.000000    0.000000    0.000000
2  1.558677e+02  1.558677e+02  4.840152e+02    0.000000    0.000000    0.000000
3  4.200000e-09  5.500000e-09  1.230000e-08  117.414353    0.000000    0.000000
4 -2.250000e-08 -1.770000e-08 -4.810000e-08    0.000000  117.414354    0.000000
5 -1.300000e-09  2.200002e-09  0.000000e+00    0.000000    0.000000  216.431859

The postprocess is done. The metrics are saved in 'metrics.json', and the elastic results are saved in 'metrics_elastic.json'.
```

计算得到的弹性张量中有6个独立分量（

$$
C_{11}
$$

、\$\$C\_{12}\$\$、\$\$C\_{13}\$\$、\$\$C\_{33}\$\$、\$\$C\_{44}\$\$、\$\$C\_{66}\$\$），符合金红石型TiO2所属的四方晶系（Laue类型I）对独立弹性张量的要求。下面的表格汇总了本文中的计算、Materials Project数据库，以及实验测量的金红石型TiO2的弹性张量中的独立元素。本文的计算结果比Materials Project数据库更好地符合实验测量结果[3]。

|                     | 本文结果 | Materials Project | 实验测量结果 |
| --------------------- | ---------- | ------------------- | -------------- |
| \$\$C\_{11}\$\$/GPa | 281.2    | 426               | 268.0 (1.4)  |
| \$\$C\_{12}\$\$/GPa | 158.1    | 2                 | 174.9 (1.4)  |
| \$\$C\_{13}\$\$/GPa | 155.9    | 149               | 147.4 (1.5)  |
| \$\$C\_{33}\$\$/GPa | 484.0    | 470               | 484.2 (1.8)  |
| \$\$C\_{44}\$\$/GPa | 117.4    | 113               | 123.8 (0.2)  |
| \$\$C\_{66}\$\$/GPa | 216.4    | 43                | 190.2 (0.5)  |

参考文献：

[1] M. de Jong et al., Charting the complete elastic properties of inorganic crystalline compounds, Scientific Data ​**2**​, 150009 (2015).

[2] S. Singh, L. Lang, V. Dovale-Farelo, U. Herath, P. Tavadze, F.-X. Coudert, and A. H. Romero, MechElastic: A Python library for analysis of mechanical and elastic properties of bulk and 2D materials, Computer Physics Communications ​**267**​, 108068 (2021).

[3] D. G. Isaak, J. D. Carnes, O. L. Anderson, H. Cynn, and E. Hake, Elasticity of TiO2 rutile to 1800 K, Physics and Chemistry of Minerals ​**26**​, 31 (1998).
