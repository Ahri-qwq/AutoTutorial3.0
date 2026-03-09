# Step 2 大纲（已确认）

## 选定方案：方案C + 适量理论背景

### 标题
使用 ABACUS+ShengBTE 计算 Si 的晶格热导率

### 结构
1. 概述（约30行）— ASCII流程图 + 案例简介
2. 物理背景（约90行）— 声子图像、力常数作用、有限位移法、ShengBTE角色
3. 前置准备（约40行）— 案例下载、目录结构、软件依赖
4. Step 1：二阶力常数（约140行）— relax→Phonopy→SCF→FORCE_CONSTANTS→单位转换
5. Step 2：三阶力常数（约160行）— thirdorder→POSCAR→STRU→批量SCF→aba2vasp→reap
6. Step 3：运行 ShengBTE（约100行）— CONTROL完整呈现 + 运行 + 结果
7. 关键注意事项（约70行）— 5个坑 + LCAO vs PW对比表
8. 进阶与展望（约30行）— 收敛性测试建议

### 目标总行数：约660行
