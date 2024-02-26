from matid import SymmetryAnalyzer
from ase.spacegroup import crystal
from ase.visualize import view

# Create an initial system with ASE
a = 5.64
nacl = crystal(
	['Na', 'Cl'],
	[(0, 0, 0), (0.5, 0.5, 0.5)],
	spacegroup=225,
    cellpar=[a, a, a, 90, 90, 90]
) * [3, 2, 1]
view(nacl)

# Let's retrieve the conventional cell as ASE Atoms
analyzer = SymmetryAnalyzer(nacl)
conventional_cell = analyzer.get_conventional_system()
view(conventional_cell)