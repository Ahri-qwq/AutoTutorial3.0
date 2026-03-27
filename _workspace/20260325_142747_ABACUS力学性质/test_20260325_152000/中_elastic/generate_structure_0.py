from ase.build import bulk

si_conv = bulk('Si', cubic=True)  # 8原子立方惯用胞
si_conv.write("Si_conv.cif")
