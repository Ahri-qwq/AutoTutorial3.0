# AutoTutorial 3.0 - 持久记忆

## 项目结构

- **根目录：** `C:\MyCode\AutoTutorial3.0\`
- **教程生成流程：** `CLAUDE.md`（7步）
- **计算测试流程：** `testCLAUDE.md`（7步）
- **工作区：** `_workspace/YYYYMMDD_HHMMSS_主题/`（必须含时间戳，例：`20260228_144416_DFT+U强关联体系计算设置`）

## 测试框架关键文件

- `tools/test_framework_integrated.py` — 主框架，注册所有插件
- `tools/test_plugins/` — 插件目录
  - `base_plugin.py` — 抽象基类（BaseTestPlugin, TestInfo, ValidationResult）
  - `relax_plugin.py`, `elastic_plugin.py`, `band_plugin.py`, `dos_plugin.py`
  - `dftu_plugin.py` — DFT+U 插件（2026-02-28 新增，NiO AFM SCF）
  - `optic_plugin.py` — 光学性质插件（2026-03-02 新增，SiO₂）
  - `solvation_plugin.py` — 隐式溶剂插件（2026-03-02 新增，H₂@水溶液）
  - `phonopy_plugin.py` — Phonopy声子谱插件（2026-03-03 新增，FCC Al SCF+力）
  - `PLUGIN_REGISTRY.md` — 插件注册表（新增插件时必须更新）
- `tools/orbital_validator.py` — 验证并修复轨道文件名 + nbands auto
- `tools/fix_stru.py` — 修复 STRU 文件格式

## 框架自修复机制（核心设计）

当 `analysis.json` 为空（`detected_types = []`）时：
1. 查阅 `tools/test_plugins/PLUGIN_REGISTRY.md`
2. 阅读教程找 INPUT 关键词，判断计算类型
3. 若无对应插件 → Step 1.5 自主创建插件（无需等用户）
4. 创建完毕后重跑 prepare，验证非空后继续

## 已支持的计算类型

| 插件 | calc_type | 关键词 |
|------|-----------|--------|
| RelaxPlugin | relax | `calculation = relax/cell-relax`（排除 phonopy 教程） |
| ElasticPlugin | elastic | `弹性常数`、`elastic` |
| BandPlugin | band | `能带结构`、`band structure`（排除 phonopy 教程） |
| DOSPlugin | dos | `态密度`、`density of states`（排除 phonopy 教程） |
| DFTUPlugin | dftu | `dft_plus_u = 1`、`DFT+U` |
| OpticPlugin | optic | `out_mat_hs2` + `OPTICAL_CONDUCTIVITY` |
| SolvationPlugin | solvation | `imp_sol = 1`、`隐式溶剂` |
| PhonopyPlugin | phonopy | `phonopy`、`声子谱`、`FORCE_SETS`、`有限位移方法` |

## DFT+U NiO 案例关键参数（已验证通过）

- STRU：Type-II AFM，a1=(9.6226,0,0), a2=(7.9999,5.3468,0), a3=(7.9999,2.4270,4.7629) Bohr
- 轨道：`Ni_gga_9au_100Ry_4s2p2d1f.orb`（双d基组），`O_gga_7au_100Ry_2s2p1d.orb`
- 赝势：`Ni_ONCV_PBE-1.0.upf`，`O_ONCV_PBE-1.0.upf`
- Job 22150516：133秒完成，5项指标全部通过

## 验证容差原则

| 物理量 | 容差 | 理由 |
|--------|------|------|
| 总能量（大负数）| 相对 0.1% | 绝对差几 eV 物理上微小 |
| 能隙 | 相对 15% | 晶格微差导致能隙~10%变化 |
| 总磁矩（应为0）| 绝对 0.1 | 直接判断近零 |
| 绝对磁矩 | 相对 5% | 结构相近时稳定 |

## Bohrium 平台

- 项目 ID：205855（【新】ABACUS功能开发与测试）
- 默认机型：c16_m32_cpu（16核32GB）
- 默认镜像：`registry.dp.tech/dptech/abacus:LTSv3.10.1`
- 命令格式：`OMP_NUM_THREADS=1 mpirun -np 8 abacus`

## 重要约定

- testCLAUDE.md 中的框架测试命令：`python tools/test_framework_integrated.py "<path>" --test-dir "<dir>" --phase prepare`（注意是位置参数，不是 `--tutorial`）
- 新增插件后：更新 `test_framework_integrated.py`、`PLUGIN_REGISTRY.md`、testCLAUDE.md 附录D 三处
- **工作目录命名**：Step 0 必须用完整格式 `_workspace/YYYYMMDD_HHMMSS_主题/`，**不能只写日期**（历史错误案例：20260302_隐式溶剂模型、20260302_光学性质计算）
