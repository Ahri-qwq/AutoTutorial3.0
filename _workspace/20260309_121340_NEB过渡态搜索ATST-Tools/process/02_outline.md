# 大纲（方案C：案例驱动渐进型）

## 标题
使用 ABACUS 和 ATST-Tools 进行 NEB/AutoNEB 过渡态搜索

## 结构

### 0. 引言（约30行）
- 本教程干什么、需要什么、学完能做什么

### 1. 过渡态与 NEB 方法（约100行）
- 过渡态与 MEP
- IT-NEB + CI-NEB 核心思路（含公式，第一次详细说）
- AutoNEB：动态映像的改进

### 2. ATST-Tools 工作流（约50行）
- 环境配置
- 脚本目录结构
- 通用工作流：neb_make → neb_run/autoneb_run → neb_post

### 3. 案例一：Li 在 Si 中的扩散（NEB 入门）（约260行）
- 体系与目标
- 初末态准备 + neb_make.py 插值（完整命令+参数说明）
- 串行 DyNEB 运行（dyneb_run.py 完整脚本+逐参解析）
- 并行 NEB 对比运行
- 后处理 + 计算结果（0.618 eV）

### 4. 案例二：单原子催化体系的 AutoNEB（约300行）
- 为什么需要 AutoNEB
- 初末态准备 + 插值
- autoneb_run.py 关键参数（重点讲与 neb_run.py 的区别）
- AutoNEB 迭代过程（6步收敛，动态映像添加过程）
- 后处理 + 计算结果（1.328 eV）

### 5. 振动分析验证过渡态（约80行）
- 唯一虚频判据
- vib_analysis.py 配置与运行
- 结果解读（705i cm⁻¹，ZPE 4.416 eV）

### 附录（约40行）
- 辅助脚本与续算
- 参考资料
