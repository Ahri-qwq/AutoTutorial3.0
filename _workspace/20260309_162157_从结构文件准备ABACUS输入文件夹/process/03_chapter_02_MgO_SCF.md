## 二、基础用法：MgO SCF 计算

以 MgO 的 CIF 文件为例，生成 LCAO 基组的 SCF 计算输入文件夹：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao --folder-syntax MgO
```

各参数含义：

| 参数 | 说明 |
|------|------|
| `-f MgO.cif` | 指定输入的结构文件 |
| `--ftype cif` | 结构文件格式（支持 `cif`、`poscar` 等 dpdata 支持的格式） |
| `--lcao` | 使用 LCAO 基组；不加此选项则使用 PW 基组（不需要轨道文件） |
| `--folder-syntax MgO` | 指定生成的任务文件夹名称；不设置则从 `000000` 开始编号 |

命令执行后，目录结构如下：

```
MgO
├── INPUT
├── Mg_gga_10au_100Ry_2s1p.orb -> /home/abc/apns-orbitals-efficiency-v1/Mg_gga_10au_100Ry_2s1p.orb
├── Mg.PD04.PBE.UPF -> /home/abc/apns-pseudopotentials-v1/Mg.PD04.PBE.UPF
├── O_gga_6au_100Ry_2s2p1d.orb -> /home/abc/apns-orbitals-efficiency-v1/O_gga_6au_100Ry_2s2p1d.orb
├── O.upf -> /home/abc/apns-pseudopotentials-v1/O.upf
├── STRU
└── struinfo.txt
MgO.cif
run.sh
setting.json
struinfo.json
```

几点说明：

- 赝势（`.UPF`）和轨道（`.orb`）文件是软链接，指向 `ABACUS_PP_PATH` / `ABACUS_ORB_PATH` 中的原始文件。如需复制文件到目录中（而非软链接），加 `--copy-pp-orb` 选项。
- `run.sh` 和 `setting.json` 用于通过 dpdispatcher 提交计算，不需要可以删除。
- `struinfo.txt` 和 `struinfo.json` 记录了原始结构路径，不需要也可以删除。

生成的 `INPUT` 文件内容如下：

```
calculation     scf
symmetry     1
ecutwfc     100
scf_thr     1e-07
scf_nmax     100
smearing_method     gauss
smearing_sigma     0.015
mixing_type     broyden
mixing_beta     0.8
basis_type     lcao
ks_solver     genelpa
precision     double  # or single
#cal_force     1
#cal_stress     1
kspacing     0.14 # unit in 1/bohr
#gamma_only     0
```

默认使用 `kspacing = 0.14`（单位 1/Bohr）自动生成 K 点，这组参数对很多体系适用，可按需修改。
