# AutoTutorial 3.0 - 持久记忆

## 项目结构

- **根目录：** `C:\MyCode\AutoTutorial3.0\`
- **教程生成流程：** `CLAUDE.md`（7步 + Step 8 测试询问，含 compact 暂停点）
- **计算测试流程：** `testCLAUDE.md`（~1346行，核心流程 Step 0–7）
- **testCLAUDE 模块目录：** `tools/testCLAUDE/`（按需读取，减少上下文压力）
  - `bohrium_setup.md` — Bohrium 配置流程 + CLI 命令速查（Step 0 未配置时加载）
  - `plugin_dev_guide.md` — 插件自主创建流程 Step 1.5（detected_types=[] 时加载）
  - `plugins_history.md` — 插件开发历史表（新增插件只追加此文件一行）
  - `troubleshooting.md` — 8 类故障排除清单（Step 4.4/6.3 出错时加载）
  - `file_formats.md` — job.json/analysis.json/STRU 格式参考（Step 2 时加载）
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
  - `sdft_plugin.py` — SDFT/MDFT 插件（2026-03-04 新增，Si SCF+Al MD+Si DOS）
  - `elf_plugin.py` — ELF 电子局域函数插件（2026-03-09 新增，H₂O PW/LCAO + Fe BCC）
  - `tddft_plugin.py` — RT-TDDFT 混合规范插件（2026-03-09 新增，Si 原胞 50步）
  - `neb_plugin.py` — NEB/ATST-Tools 插件（2026-03-09 新增，Li-Si IS SCF 验证）
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
| ElasticPlugin | elastic | `弹性常数`、`elastic`（排除 NEB 教程，因"Nudged Elastic Band"含 elastic 词）|
| BandPlugin | band | `能带结构`、`band structure`（排除 phonopy 教程） |
| DOSPlugin | dos | `态密度`、`density of states`（排除 phonopy 教程） |
| DFTUPlugin | dftu | `dft_plus_u = 1`、`DFT+U` |
| OpticPlugin | optic | `out_mat_hs2` + `OPTICAL_CONDUCTIVITY` |
| SolvationPlugin | solvation | `imp_sol = 1`、`隐式溶剂` |
| PhonopyPlugin | phonopy | `phonopy`、`声子谱`、`FORCE_SETS`、`有限位移方法` |
| SDFTPlugin | sdft | `esolver_type = sdft`（INPUT代码块中，优先于 DOSPlugin）|
| ELFPlugin | elf | `out_elf 1`、`ELF.cube`、`电子局域函数`（多案例单插件，案例存入 expected_results['cases']）|
| TDDFTPlugin | tddft | `esolver_type = tddft`（INPUT代码块中）|
| NEBPlugin | neb | `AbacusNEB`、`AbacusAutoNEB`、`ATST-Tools`、`dyneb_run`、`autoneb_run`；测试 IS 单点 SCF（全 NEB 太贵） |

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

## 写作格式约定

- **表格中禁止使用 `**加粗**`**：表格单元格内不得使用双星号加粗，保持纯文本

## 重要约定

- testCLAUDE.md 中的框架测试命令：`python tools/test_framework_integrated.py "<path>" --test-dir "<dir>" --phase prepare`（注意是位置参数，不是 `--tutorial`）
- 新增插件后：更新 `test_framework_integrated.py`、`PLUGIN_REGISTRY.md`、`tools/testCLAUDE/plugins_history.md` 三处（不再需要修改 testCLAUDE.md 主文件）
- **工作目录命名**：Step 0 必须用完整格式 `_workspace/YYYYMMDD_HHMMSS_主题/`，**不能只写日期**（历史错误案例：20260302_隐式溶剂模型、20260302_光学性质计算）
- **时间戳必须从 Bash 获取**：系统上下文只提供日期，不提供时间。Step 0 创建目录前必须执行 `date +"%Y%m%d_%H%M%S"` 获取真实时间戳，**严禁编造时间**（历史错误案例：20260304_**100000**_ABACUS+ShengBTE、20260304_**100000**_ABACUS+SDFT，均因未执行 date 命令导致时间伪造）

## ABACUS TDDFT 参数注意事项（2026-03-09 验证）

- `calculation = md` + `esolver_type = tddft` 时，必须显式设置 `md_type nve`（字符串）避免 NHC 报错
- ABACUS v3.10.1 的 `md_type` 取字符串值：`nve`、`nhc` 等（不是整数 0/-1）
- `md_tfirst 0` 配合 `md_type nve` → 固定离子的纯电子 TDDFT
- 电流输出文件名为 `current_total.dat`（不是 `SPIN1_CURRENT`）
- `continue` 模式 bug 已修复：须先用 state 文件恢复 tutorial_path，再调 phase1_analyze() 填充 plugin_tests，再 _load_state() 恢复 job_ids

## 轨道文件备用下载源

`tools/test_framework_phase3_7_impl.py` 的 `get_file` 已扩展：
- 主源：`deepmodeling/abacus-develop/tests/PP_ORB/`（不含 Li 等轻元素特殊轨道）
- 备用源：`abacusmodeling/ABACUS-orbitals/SG15_v1.0/Orbitals_v2.0/<Element>_DZP/`
  - URL 格式：`https://raw.githubusercontent.com/abacusmodeling/ABACUS-orbitals/main/SG15_v1.0/Orbitals_v2.0/{Element}_{quality}/{filename}`
  - 已验证：`Li_gga_8au_100Ry_4s1p.orb`（4s1p=DZP）已缓存到 `tools/orbitals/`
