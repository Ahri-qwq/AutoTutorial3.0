# 故障排除清单

> 本文件由 testCLAUDE.md 按需 read，在 Step 4.4（任务失败）或 Step 6.3（结果偏差）时加载。

---

## 附录B: 故障排除清单

### 问题1：Bohrium CLI未找到

**症状：**
```
'bohr' 不是内部或外部命令
```

**解决方案：**
1. 确认已安装Bohrium CLI
2. 重启VSCode（PATH更新）
3. 在PowerShell中运行`refreshenv`

---

### 问题2：ACCESS_KEY未设置

**症状：**
```
echo $env:ACCESS_KEY
# 输出为空
```

**解决方案：**
1. 在PowerShell中运行：`setx ACCESS_KEY "your_key"`
2. 完全关闭并重新打开VSCode（不是重新加载窗口）
3. 验证：`echo $env:ACCESS_KEY`

---

### 问题3：API连接失败

**症状：**
```
bohr project list
Error: Failed to connect to API
```

**解决方案：**
1. 检查AccessKey是否正确
2. 检查网络连接
3. 尝试访问 https://bohrium.dp.tech/
4. 检查防火墙设置

---

### 问题4：任务提交失败

**症状：**
```
Error: Job submission failed
```

**可能原因：**
1. 项目ID未设置或错误
2. 项目配额不足
3. 机型或镜像不可用
4. job.json格式错误

**解决方案：**
1. 检查项目ID：`bohr config list`
2. 检查项目配额：`bohr project list`
3. 检查机型可用性：`bohr machine list`
4. 验证job.json格式

---

### 问题5：轨道文件下载为HTML

**症状：**
下载的`.orb`文件内容是HTML（404页面）

**解决方案：**
- 自动应用文件映射表（已内置）
- 常见映射：
  - `Si_gga_7au_100Ry_2s2p1d.orb` → `Si_gga_8au_100Ry_2s2p1d.orb`
  - `Ti_gga_8au_100Ry_2s2p2d1f.orb` → `Ti_gga_8au_100Ry_4s2p2d1f.orb`
  - `O_gga_6au_100Ry_2s2p1d.orb` → `O_gga_7au_100Ry_2s2p1d.orb`

---

### 问题6：STRU文件格式错误

**症状：**
```
ERROR: autotest2006 has detected a ABACUS STRU file format error
```

**解决方案：**
- 自动修复工具已内置（`fix_stru.py`）
- 新格式要求：
  - `ATOMIC_SPECIES`只包含3个字段
  - 轨道文件单独放在`NUMERICAL_ORBITAL`块

---

### 问题7：计算结果与预期差异大

**症状：**
相对误差超过5%

**可能原因：**
1. 教程预期结果不准确
2. 计算参数不同（K点、截断能）
3. ABACUS版本差异
4. 赝势或轨道文件不同

**解决方案：**
1. 检查教程中的预期结果来源
2. 对比INPUT文件参数
3. 查看ABACUS版本（教程生成时 vs 当前测试）
4. 检查赝势和轨道文件版本

---

### 问题8：任务长时间Waiting

**症状：**
任务提交后一直处于Waiting状态

**可能原因：**
1. 平台资源紧张
2. 机型不可用
3. 项目优先级低

**解决方案：**
1. 等待（通常5-15分钟）
2. 更换机型（如c8_m16_cpu）
3. 联系平台管理员
