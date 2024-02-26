from matid import SymmetryAnalyzer
from ase.spacegroup import crystal
from ase.io import write
from ase.visualize import view

a = 5.64
nacl1 = crystal(
	['Na', 'Cl'],
	[(0, 0, 0), (0.5, 0.5, 0.5)],
	spacegroup=225,
    cellpar=[a, a, a, 90, 90, 90]
)
# view(nacl1)
write(f"NaCl1.eps", nacl1, rotation="60y,20x,0y", show_unit_cell=2, scale=50, maxwidth=5000)

a = 5.64
nacl2 = crystal(
	['Cl', 'Na'],
	[(0, 0, 0), (0.5, 0.5, 0.5)],
	spacegroup=225,
    cellpar=[a, a, a, 90, 90, 90]
)
# view(nacl2)
write(f"NaCl2.eps", nacl2, rotation="60y,20x,0y", show_unit_cell=2, scale=50, maxwidth=5000)

# Let's retrieve the conventional cell as ASE Atoms
analyzer = SymmetryAnalyzer(nacl1)
conventional_cell = analyzer.get_conventional_system()
# view(conventional_cell)

