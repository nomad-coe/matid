from matid import SymmetryAnalyzer
from ase.io import read

# Load clathrate structure from file
clathrate = read('clathrate.xyz')

# Print the wyckoff sets
symm = SymmetryAnalyzer(clathrate)
wyckoff_sets_conv = symm.get_wyckoff_sets_conventional()

for i_group, group in enumerate(wyckoff_sets_conv):
    print("Set {}".format(i_group))
    print("    Letter: {}".format(group.wyckoff_letter))
    print("    Element: {}".format(group.element))
    print("    Indices: {}".format(group.indices))
    print("    Multiplicity: {}".format(group.multiplicity))
    print("    Repr.: {}".format(group.representative))
    print("    x: {}".format(group.x))
    print("    y: {}".format(group.y))
    print("    z: {}".format(group.z))
