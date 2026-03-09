# 文件格式参考

> 本文件由 testCLAUDE.md 按需 read，在 Step 2 生成输入文件或调试格式问题时加载。

---

## 附录C: 文件格式参考

### job.json格式

```json
{
  "job_name": "任务名称",
  "command": "执行命令",
  "log_file": "日志文件路径（支持通配符）",
  "project_id": 205855,
  "platform": "ali",
  "job_type": "container",
  "machine_type": "c16_m32_cpu",
  "image_address": "registry.dp.tech/dptech/abacus:LTSv3.10.1"
}
```

**字段说明：**
- `job_name`: 任务名称（用于识别）
- `command`: 容器内执行的命令
- `log_file`: 日志文件路径（用于判断任务完成）
- `project_id`: Bohrium项目ID
- `platform`: 平台类型（ali/aws）
- `job_type`: 任务类型（container）
- `machine_type`: 机型（c8_m16_cpu, c16_m32_cpu等）
- `image_address`: Docker镜像地址

**常用机型：**
- `c8_m16_cpu`: 8核16GB（适合小计算）
- `c16_m32_cpu`: 16核32GB（推荐）
- `c32_m64_cpu`: 32核64GB（大计算）

**ABACUS镜像：**
- `registry.dp.tech/dptech/abacus:3.10.1`: 稳定版
- `registry.dp.tech/dptech/abacus:LTSv3.10.1`: 长期支持版（推荐）
- `registry.dp.tech/dptech/abacus:latest`: 最新版

---

### analysis.json格式

`cases` 数组有几个案例，就有几个条目。**单案例教程只有 1 个条目，多案例教程有 N 个**，结构相同：

```json
{
  "tutorial_title": "教程标题",
  "tutorial_path": "教程路径",
  "cases": [
    {
      "name": "<案例名>",
      "type": "elastic",
      "steps": [
        {
          "step": "relax",
          "input_files": { "INPUT": "...", "STRU": "...", "KPT": "..." },
          "pseudopotentials": ["<赝势文件>"],
          "orbitals": ["<轨道文件>"]
        }
      ],
      "expected_results": {
        "<参数名>": <数值>
      },
      "tolerance": 0.05
    }
    // 若有多个案例，在此继续添加相同结构的条目
  ]
}
```

**注意：** `cases` 数组必须包含教程中**所有**案例，缺一不可。

---

### STRU文件新格式

```
ATOMIC_SPECIES
Si 28.085 Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
./Si_gga_8au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
10.2

LATTICE_VECTORS
1 0 0
0 1 0
0 0 1

ATOMIC_POSITIONS
Cartesian
Si
0.0
2
0.0 0.0 0.0 1 1 1
0.25 0.25 0.25 1 1 1
```

**关键要点：**
- `ATOMIC_SPECIES`只包含3个字段：元素、质量、赝势
- `NUMERICAL_ORBITAL`独立块，每行一个轨道文件
- 轨道文件路径需要`./`前缀
