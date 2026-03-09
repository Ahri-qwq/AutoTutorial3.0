# 调研结果

## 知识检索评估
- 检索结果主要来自案例文件本身（知识库存有多个副本）
- 无额外补充文档，但案例内容完整，包含：
  - 过渡态理论（NEB/AutoNEB原理）
  - 完整代码脚本（neb_run.py, autoneb_run.py, vib_analysis.py）
  - 两个实际计算案例（Li-diffu-Si, Cy-Pt@graphene）
  - 计算输出日志（完整收敛过程）

## 案例要点提炼

### 工具链
- ATST-Tools v1.5.0 = ASE + ASE-ABACUS + 过渡态工作流脚本
- PYTHONPATH: /opt/ATST-Tools/source
- 目录结构：neb/, autoneb/, dimer/, sella/, vibration/, relax/, source/

### 完整工作流步骤
1. 准备初末态结构（结构优化完毕）
2. neb_make.py 生成初猜 NEB 链（IDPP 插值）→ init_neb_chain.traj
3. neb_run.py / autoneb_run.py 运行计算 → neb.traj / run_autoneb???.traj
4. neb_post.py 后处理 → nebplots.pdf, neb_latest.traj
5. vib_analysis.py 振动分析验证过渡态

### 案例一参数（Li-diffu-Si）
- 赝势/轨道：Li_ONCV_PBE-1.2.upf, Si_ONCV_PBE-1.2.upf
- Li_gga_8au_100Ry_4s1p.orb, Si_gga_8au_100Ry_2s2p1d.orb
- kpts: [2,2,2]; nspin=1; basis_type=lcao
- DyNEB: mpi=16, omp=1, dyneb=True, parallel=False
- 并行NEB: mpi=5, omp=1, parallel=True（映像数+2进程）
- 能垒: 0.618 eV

### 案例二参数（Cy-Pt@graphene AutoNEB）
- 赝势/轨道：C/H/Pt ONCV + gga orb
- kpts: [2,1,2]; nspin=2; vdw_method=d3_bj; efield_flag/dip_cor_flag=1
- fmax=[0.20, 0.05]; n_images=10; n_simul=4
- 运行命令：mpirun -np 4 gpaw python autoneb_run.py
- 能垒: 1.328 eV

### 振动分析
- 虚频: 87.4i meV (705.0i cm⁻¹) → 验证唯一虚频
- ZPE: 4.416 eV
- 自由能校正: T=523.15K, HarmonicThermo

## 风格特征（参考文章）
- 磁性教程：直接给参数表格，几乎无废话，代码块完整
- DeePMD教程：分章节，代码配说明，命令行操作清晰
- 共同特点：没有"综上所述"等套话，语言简洁直接
