from matid import SymmetryAnalyzer
from ase.spacegroup import crystal
from ase.visualize import view

# Create an initial system with ASE
a = 9.04
skutterudite = crystal(
    ('Co', 'Sb'),
    basis=[(0.25, 0.25, 0.25), (0.0, 0.335, 0.158)],
    spacegroup=204,
    cellpar=[a, a, a, 90, 90, 90]
)
view(skutterudite)

# Let's retrieve the material id
analyzer = SymmetryAnalyzer(skutterudite)
material_id = analyzer.get_material_id()
print(material_id)
