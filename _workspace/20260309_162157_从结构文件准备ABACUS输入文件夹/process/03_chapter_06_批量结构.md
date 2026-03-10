## 六、批量结构准备

`-f` 参数支持同时传入多个结构文件，配合 `--folder-syntax` 可以批量生成各自的输入文件夹。

以 Pd(100) 表面能随层数的收敛性测试为例，结构文件如下：

```
Pd100_1layer.vasp  Pd100_3layer.vasp  Pd100_5layer.vasp  Pd100_7layer.vasp
Pd100_2layer.vasp  Pd100_4layer.vasp  Pd100_6layer.vasp  Pd100_8layer.vasp
```

使用以下命令为所有结构准备优化计算的输入文件夹：

```bash
abacustest model inputs -f Pd100_*layer.vasp --ftype poscar --lcao --jtype relax --folder-syntax "x[:-5]"
```

`--folder-syntax "x[:-5]"` 中的 `x` 代表结构文件名，`[:-5]` 是 Python 字符串切片，含义为去掉文件名末尾 5 个字符（即 `.vasp`）。执行后，每个结构生成同名的任务文件夹：

```
Pd100_1layer/
Pd100_2layer/
...
Pd100_8layer/
```

除了字符串切片，`--folder-syntax` 支持任意合法的 Python 字符串表达式，`x` 代表原文件名（含扩展名），例如：

- `"x[:-5]"`：去掉 `.vasp` 后缀
- `"x.split('.')[0]"`：取点号前的部分
- `"MgO_" + x[:-4]"`：添加前缀
