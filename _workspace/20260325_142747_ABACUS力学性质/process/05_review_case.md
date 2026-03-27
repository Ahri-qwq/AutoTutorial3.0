# 案例审查报告

## 5.1 案例一：h-BN cell-relax

**文件结构完整性：**
- [x] INPUT 文件完整（21个参数，全部来自知识库 abacus_user_guide__abacus_eff2.html.md）
- [x] STRU 关键头部（ATOMIC_SPECIES, NUMERICAL_ORBITAL, LATTICE_CONSTANT, LATTICE_VECTORS）已呈现
- [x] 明确说明"完整192原子坐标从官方GitHub获取"，避免编造坐标
- [ ] KPT：本案例使用 `kspacing 0.08` 代替手写KPT文件，已在INPUT中体现，无单独KPT文件，符合实际

**参数准确性核查：**
- `suffix BN` ✓ （来源文档使用 suffix BN）
- `ecutwfc 100` ✓
- `scf_thr 1e-07` ✓
- `smearing_method gauss` ✓
- `smearing_sigma 0.002` ✓
- `mixing_type pulay` ✓（来源文档）
- `mixing_beta 0.3` ✓
- `force_thr_ev 0.01` ✓
- `stress_thr 0.5` ✓
- `kspacing 0.08` ✓
- LATTICE_CONSTANT: `1.8897261258369282` ✓
- 晶格矢量 a1: `2.6665300000` ✓，a2: `16.970664000` ✓
- a3 c方向：改为 `28.0218`（知识库原文截断，无法确认完整值）→ 已标注为参考值

**轨道赝势：**
- `B_ONCV_PBE-1.0.upf` + `B_gga_8au_100Ry_2s2p1d.orb` ✓（来源文档）
- `N_ONCV_PBE-1.0.upf` + `N_gga_8au_100Ry_2s2p1d.orb` ✓

**收敛日志示例：**
- TOTAL-PRESSURE=-2.070e+00 和 -4.350e-01 两个值来自 `ABACUS 使用教程｜结构优化.md` 的参考示例 ✓
- 已加"示意性输出"说明 ✓

## 5.2 案例二：Si 弹性常数

**数据完整性：**
- [x] 弹性张量完整输出（6×6矩阵）完整引用
- [x] C11=165.5, C12=58.6, C44=82.5 GPa ✓（与原始文件一致）
- [x] BV=94.19, GV=70.89, EV=170.02, ν=0.199 ✓（与原始文件 94.191857, 70.892488, 170.022310, 0.199156 一致）
- [x] 与Materials Project对比（C11:153, C12:57, C44:74）✓（来源文件已验证）
- [x] 与实验值对比（C11实验=165.7 GPa）✓

**流程完整性：**
- [x] ase生成结构代码 ✓
- [x] abacustest model inputs 命令 ✓
- [x] 整理目录步骤（cp命令）✓
- [x] abacustest model elastic prepare 命令 ✓
- [x] job.json格式（Bohrium平台）✓
- [x] abacustest model elastic post 命令 ✓
- [x] 完整屏幕输出 ✓（直接引用原文）

**job.json参数：**
- project_id: 205855 来自 MEMORY.md（用户项目ID）✓
- machine_type: c16_m32_cpu ✓（MEMORY.md默认机型）
- image_address: LTSv3.10.1 ✓（MEMORY.md默认镜像）

## 审查结论

两个案例数据准确，来源可追溯，无编造数值。
h-BN STRU第三晶格矢量c分量已标注为参考值，诚实告知需从官方获取完整坐标。
案例审查通过，进入 Step 6 风格审查。
