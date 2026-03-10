## 五、自定义 INPUT 和 K 点

默认情况下，INPUT 中使用 `kspacing = 0.14` 自动生成 K 点。如需自定义参数，可通过 `--input` 指定 INPUT 模板，通过 `--kpt` 显式设置 K 点网格。

以 Fe₂O₃ 为例，将 `smearing_sigma` 改为 0.001，`mixing_beta` 降至 0.2，并使用 5×5×5 K 点。

先准备 INPUT 模板文件（`INPUT_template`），只写需要覆盖的参数：

```
INPUT_PARAMETERS
smearing_sigma     0.001
mixing_beta     0.2
```

然后执行：

```bash
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 --input INPUT_template --kpt 5 5 5 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

`INPUT_template` 中的参数会替换默认值，同时生成以 Gamma 为中心、采样密度为 5×5×5 的 KPT 文件，并去掉 INPUT 中的 `kspacing` 参数。
