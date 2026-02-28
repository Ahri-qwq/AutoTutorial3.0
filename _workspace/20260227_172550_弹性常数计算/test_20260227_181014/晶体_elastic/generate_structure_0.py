from ase.build import bulk

si_conv = bulk('Si', cubic=True)
si_conv.write("Si_conv.cif")
