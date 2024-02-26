from matid import SymmetryAnalyzer
from ase.spacegroup import crystal
from ase.io import write
from ase import Atoms
from ase.visualize import view
import spglib

def get_spglib_conventional(system):
    cell = (
        system.get_cell(),
        system.get_scaled_positions(),
        system.get_atomic_numbers(),
    )
    dataset = spglib.get_symmetry_dataset(cell)

    return Atoms(
        symbols=dataset["std_types"],
        scaled_positions=dataset["std_positions"],
        cell=dataset["std_lattice"],
    )

# Lets define two variants of NaCl in rocksalt structure
a = 5.64
nacl1 = crystal(
	['Na', 'Cl'],
	[(0, 0, 0), (0.5, 0.5, 0.5)],
	spacegroup=225,
    cellpar=[a, a, a, 90, 90, 90]
)
write(f"original1.eps", nacl1, rotation="60y,20x,0y", show_unit_cell=2, scale=50, maxwidth=5000)
a = 5.64
nacl2 = crystal(
	['Cl', 'Na'],
	[(0, 0, 0), (0.5, 0.5, 0.5)],
	spacegroup=225,
    cellpar=[a, a, a, 90, 90, 90]
)
write(f"original2.eps", nacl2, rotation="60y,20x,0y", show_unit_cell=2, scale=50, maxwidth=5000)

# Lets retrieve the conventional cell using spglib
conventional1_spglib = get_spglib_conventional(nacl1)
conventional2_spglib = get_spglib_conventional(nacl2)
write(f"spglib1.eps", conventional1_spglib, rotation="60y,20x,0y", show_unit_cell=2, scale=50, maxwidth=5000)
write(f"spglib2.eps", conventional2_spglib, rotation="60y,20x,0y", show_unit_cell=2, scale=50, maxwidth=5000)

# Let's retrieve the conventional cell as ASE Atoms using MatID
analyzer = SymmetryAnalyzer(nacl1)
conventional1_matid = analyzer.get_conventional_system()
print(analyzer.get_wyckoff_sets_conventional())
analyzer = SymmetryAnalyzer(nacl2)
print(analyzer.get_wyckoff_sets_conventional())
conventional2_matid = analyzer.get_conventional_system()
write(f"matid1.eps", conventional1_matid, rotation="60y,20x,0y", show_unit_cell=2, scale=50, maxwidth=5000)
write(f"matid2.eps", conventional2_matid, rotation="60y,20x,0y", show_unit_cell=2, scale=50, maxwidth=5000)

