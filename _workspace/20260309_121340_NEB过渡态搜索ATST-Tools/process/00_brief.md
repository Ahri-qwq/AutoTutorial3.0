# 任务简报

## 基本信息
- **任务类型：** C（案例驱动教程）
- **主题：** NEB / AutoNEB 过渡态搜索，基于 ATST-Tools + ABACUS
- **案例文件：** data/input/NEB过渡态计算与ATST-Tools.md
- **拟定标题：** 待大纲确定后定

## 案例核心内容

### 工具链
- **ATST-Tools**（Advanced ASE Transition State Tools for ABACUS）v1.5.0
- 基于 ASE + ASE-ABACUS 接口 + ABACUS 计算后端
- 核心脚本：`neb_make.py`、`neb_run.py`、`autoneb_run.py`、`neb_post.py`、`vib_analysis.py`

### 案例一：Li在Si中的扩散（串并行NEB）
- 体系：Li diffusion in Si
- 插值方法：Pymatgen IDPP，3个映像
- 串行DyNEB参数：mpi=16, omp=1, fmax=0.05, algorism=improvedtangent, climb=True, dyneb=True, parallel=False
- 并行NEB参数：mpi=5, omp=1, fmax=0.05, parallel=True
- 赝势：Li_ONCV_PBE-1.2.upf, Si_ONCV_PBE-1.2.upf
- 轨道：Li_gga_8au_100Ry_4s1p.orb, Si_gga_8au_100Ry_2s2p1d.orb
- kpts: [2,2,2]
- 计算结果：能垒 0.618 eV（DyNEB）/ 0.618 eV（并行NEB）

### 案例二：Cy-Pt@graphene AutoNEB（环己烷在Pt负载石墨烯上的解离）
- 体系：Cyclohexane C-H bond dissociation on Pt-doped graphene
- 插值方法：Pymatgen IDPP，4个初始映像
- AutoNEB参数：mpi=16, omp=4, algorism=improvedtangent, climb=True,
  fmax=[0.20, 0.05], n_simul=4, n_images=10, k=0.10
- 赝势：C_ONCV_PBE-1.0.upf, H_ONCV_PBE-1.0.upf, Pt_ONCV_PBE-1.0.upf
- 轨道：C_gga_7au_100Ry_2s2p1d.orb, H_gga_6au_100Ry_2s1p.orb, Pt_gga_7au_100Ry_4s2p2d1f.orb
- kpts: [2,1,2]
- AutoNEB迭代过程：6步（iter001~iter006），从6个映像动态增至10个，最后CI-NEB
- 计算结果：能垒 1.328 eV，反应能差 0.389 eV
- AutoNEB iter006（CI-NEB）：88步FIRE收敛，fmax=0.040

### 振动分析
- 方法：有限差分法（ASE Vibrations + ABACUS）
- 脚本：vib_analysis.py（neb2vib从NEB链提取参与原子）
- 结果：单一虚频 87.4i meV（705.0i cm⁻¹），ZPE = 4.416 eV
- 自由能校正：T=523.15 K，HarmonicThermo

## 风格目标
- 简洁、无AI腔；代码完整；以案例为中心
- 参考：ABACUS磁性材料教程（直接给参数，少废话）
